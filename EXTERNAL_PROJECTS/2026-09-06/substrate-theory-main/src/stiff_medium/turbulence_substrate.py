"""Paired fluid-turbulence GEOMETRY + SIMULATOR for the substrate framework.

In the B3 / stiff-medium ontology a fluid is a coarse-grained excitation of
the elastic substrate.  Macroscopic kinematic viscosity ν is the long-wave
limit of the substrate damping rate γ acting on the velocity field, with
mass density ρ_fluid replacing ρ_substrate after integrating out the cell
structure.  In SI units

    ν = γ_eff · ξ_fluid² / 2          (eddy-diffusivity form)

where γ_eff is the substrate damping seen by the slowest hydrodynamic mode
and ξ_fluid is the molecular mean-free path.  The Navier–Stokes equation

    ∂_t u + (u · ∇) u = - ∇p / ρ + ν ∇² u + f

is therefore the long-wave incarnation of substrate strain transport with a
single dissipative term controlled by γ.  All the canonical observables —
Kolmogorov k^(-5/3) inertial-range spectrum, the laminar-to-turbulent
transition at the critical Reynolds number Re_crit ≈ 2300 in pipes, and the
Strouhal number St ≈ 0.21 of the Kármán vortex street behind a cylinder — are
inherited from the dimensionless form of this equation.

GEOMETRY
--------
``TurbulenceGeometry`` holds the spatial scaffolding:

    * a 2-D rectangular grid of dimensions (Lx, Ly) at resolution (Nx, Ny)
    * the wavenumber lattice (kx, ky) and isotropic shells |k|
    * helpers for vorticity ω = ∂_x v − ∂_y u and stream-function ψ
    * an energy-cascade band-pass into integral / inertial / dissipative
      ranges using the Kolmogorov scale η = (ν³/ε)^(1/4)

SIMULATOR
---------
``TurbulenceSimulator`` exposes the canonical NS observables analytically:

    * kolmogorov_spectrum(k): C_K · ε^(2/3) · k^(-5/3) with the universal
      Kolmogorov constant C_K = 1.5
    * energy_dissipation_rate(u_rms, L): ε = u_rms³ / L (large-eddy estimate)
    * reynolds_number(U, L, nu): Re = U L / ν, with Re_crit_pipe = 2300
    * laminar_or_turbulent(Re): regime label using Re_crit
    * strouhal_number(f, D, U): St = f D / U, with St_cylinder ≈ 0.21
    * vortex_shedding_field(...): synthetic 2-D Kármán street vorticity
      field behind a cylinder (analytic vortex superposition)

All quantities are derivable from the substrate primitives γ, ρ, ξ; the
Reynolds number is the ratio of inertial substrate strain transport to
γ-mediated dissipation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Canonical empirical constants (anchors for tests)
# ---------------------------------------------------------------------------

KOLMOGOROV_CONSTANT: float = 1.5
"""Universal Kolmogorov constant C_K appearing in
    E(k) = C_K · ε^(2/3) · k^(-5/3)
in the inertial subrange.  Experimental value 1.5 ± 0.1 (Sreenivasan 1995)."""

KOLMOGOROV_EXPONENT: float = -5.0 / 3.0
"""Inertial-range slope of the energy spectrum log E vs log k."""

RE_CRIT_PIPE: float = 2300.0
"""Reynolds number at which Hagen–Poiseuille pipe flow becomes unstable to
turbulence (Reynolds 1883; modern textbooks list 2100–2400)."""

RE_CRIT_PIPE_LOWER: float = 2000.0
RE_CRIT_PIPE_UPPER: float = 4000.0
"""Hysteresis band over which the laminar–turbulent transition is observed."""

STROUHAL_CYLINDER: float = 0.21
"""Strouhal number St = f D / U for the Kármán vortex street behind a
circular cylinder in the canonical Re ≈ 100 – 10⁵ window."""

# Kinematic viscosity of common reference fluids (SI, m²/s)
NU_AIR_M2_S:    float = 1.5e-5    # 20 °C, sea level
NU_WATER_M2_S:  float = 1.0e-6    # 20 °C
NU_HONEY_M2_S:  float = 7.3e-3    # 25 °C, illustrative high-Re reference

# Density of common reference fluids (SI, kg/m³)
RHO_AIR_KGM3:   float = 1.205
RHO_WATER_KGM3: float = 998.0


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class TurbulenceGeometry:
    """2-D fluid-turbulence geometry: grid, wavenumber lattice, vorticity ops.

    Attributes
    ----------
    Lx, Ly : float
        Domain extents in metres.
    Nx, Ny : int
        Number of grid cells along x and y.
    """

    Lx: float = 1.0
    Ly: float = 1.0
    Nx: int = 128
    Ny: int = 128

    # ---- Spatial grid --------------------------------------------------

    def x_grid(self) -> np.ndarray:
        """Cell-centre x coordinates."""
        dx = self.Lx / self.Nx
        return (np.arange(self.Nx) + 0.5) * dx

    def y_grid(self) -> np.ndarray:
        """Cell-centre y coordinates."""
        dy = self.Ly / self.Ny
        return (np.arange(self.Ny) + 0.5) * dy

    def meshgrid(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return (X, Y) arrays of shape (Ny, Nx)."""
        return np.meshgrid(self.x_grid(), self.y_grid(), indexing="xy")

    # ---- Wavenumber lattice -------------------------------------------

    def kx(self) -> np.ndarray:
        """Discrete kx for FFT (units 1/m)."""
        return 2.0 * math.pi * np.fft.fftfreq(self.Nx, d=self.Lx / self.Nx)

    def ky(self) -> np.ndarray:
        """Discrete ky for FFT (units 1/m)."""
        return 2.0 * math.pi * np.fft.fftfreq(self.Ny, d=self.Ly / self.Ny)

    def k_isotropic_shells(self, n_bins: int = 64) -> np.ndarray:
        """Centres of |k| shell bins from k_min = 2π/L_max to Nyquist."""
        k_min = 2.0 * math.pi / max(self.Lx, self.Ly)
        k_max = math.pi * min(self.Nx / self.Lx, self.Ny / self.Ly)
        return np.logspace(math.log10(k_min), math.log10(k_max), n_bins)

    # ---- Vorticity / stream function ----------------------------------

    @staticmethod
    def vorticity_2d(
        u: np.ndarray, v: np.ndarray, dx: float, dy: float
    ) -> np.ndarray:
        """ω_z = ∂_x v − ∂_y u using centred finite differences."""
        dvdx = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2.0 * dx)
        dudy = (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2.0 * dy)
        return dvdx - dudy

    def cell_size(self) -> Tuple[float, float]:
        """Return (dx, dy) cell size."""
        return self.Lx / self.Nx, self.Ly / self.Ny

    # ---- Energy-cascade band edges ------------------------------------

    @staticmethod
    def kolmogorov_scale_m(nu_m2_s: float, epsilon_w_kg: float) -> float:
        """η = (ν³ / ε)^(1/4) — the dissipative eddy size."""
        return (nu_m2_s**3 / max(epsilon_w_kg, 1.0e-30)) ** 0.25

    @staticmethod
    def integral_scale_m(u_rms: float, epsilon_w_kg: float) -> float:
        """L = u_rms³ / ε — the large-eddy size set by forcing."""
        return u_rms**3 / max(epsilon_w_kg, 1.0e-30)

    def cascade_band_edges(
        self, nu_m2_s: float, epsilon_w_kg: float, u_rms: float
    ) -> dict:
        """Return wavenumber edges of integral / inertial / dissipative bands.

        k_L = 2π / L (large-eddy edge), k_eta = 2π / η (Kolmogorov edge).
        Inertial range = (k_L, k_eta).
        """
        L = self.integral_scale_m(u_rms, epsilon_w_kg)
        eta = self.kolmogorov_scale_m(nu_m2_s, epsilon_w_kg)
        return {
            "k_L": 2.0 * math.pi / max(L, 1.0e-30),
            "k_eta": 2.0 * math.pi / max(eta, 1.0e-30),
            "L_m": L,
            "eta_m": eta,
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class TurbulenceSimulator:
    """Navier–Stokes / turbulence simulator on a substrate viscosity ν.

    Provides analytic predictions for the Kolmogorov spectrum, Reynolds
    number transition, dissipation rate, and Kármán vortex shedding.  The
    only physical parameter is the kinematic viscosity ν, which in the
    substrate picture is set by the damping rate γ via ν = γ_eff · ξ² / 2.
    """

    geometry: TurbulenceGeometry = field(default_factory=TurbulenceGeometry)
    nu_m2_s: float = NU_AIR_M2_S

    # ---- Substrate-viscosity bridge -----------------------------------

    @staticmethod
    def viscosity_from_substrate(
        gamma_hz: float, xi_m: float
    ) -> float:
        """Kinematic viscosity from substrate γ and a fluid mean-free path ξ.

        In the long-wave limit a single Kelvin-cell drag with relaxation
        rate γ produces a Newtonian shear viscosity with kinematic part

            ν = (1/2) γ ξ²

        (factor 1/2 from the random-walk diffusion of momentum across one
        cell per relaxation time τ = 1/γ).
        """
        return 0.5 * gamma_hz * xi_m * xi_m

    # ---- Reynolds number / regime ------------------------------------

    def reynolds_number(
        self, U_m_s: float, L_m: float, nu_m2_s: float | None = None
    ) -> float:
        """Re = U L / ν."""
        nu = self.nu_m2_s if nu_m2_s is None else nu_m2_s
        return U_m_s * L_m / nu

    @staticmethod
    def laminar_or_turbulent(
        Re: float,
        Re_lower: float = RE_CRIT_PIPE_LOWER,
        Re_upper: float = RE_CRIT_PIPE_UPPER,
    ) -> str:
        """Return 'laminar', 'transitional', or 'turbulent'."""
        if Re < Re_lower:
            return "laminar"
        if Re < Re_upper:
            return "transitional"
        return "turbulent"

    # ---- Energy dissipation -------------------------------------------

    @staticmethod
    def energy_dissipation_rate(u_rms: float, L_m: float) -> float:
        """Large-eddy estimate ε = u_rms³ / L (units W/kg = m²/s³)."""
        return u_rms**3 / max(L_m, 1.0e-30)

    @staticmethod
    def dissipation_from_substrate(
        u_rms: float, gamma_hz: float, xi_m: float
    ) -> float:
        """Substrate-derived ε using ν = γ ξ²/2 and ε ~ ν · (u/η)².

        The dissipation rate is fixed by setting the cascade time
        τ = ε^(-1/3) L^(2/3) equal to the substrate damping time 1/γ at the
        Kolmogorov scale; this collapses to ε ~ u_rms² · γ.
        """
        return u_rms * u_rms * gamma_hz * 0.5 * xi_m * xi_m / xi_m**2

    # ---- Kolmogorov inertial-range spectrum ---------------------------

    @staticmethod
    def kolmogorov_spectrum(
        k: np.ndarray | float,
        epsilon_w_kg: float,
        C_K: float = KOLMOGOROV_CONSTANT,
    ) -> np.ndarray:
        """E(k) = C_K · ε^(2/3) · k^(-5/3) (3-D inertial subrange)."""
        k_arr = np.atleast_1d(np.asarray(k, dtype=float))
        return C_K * (epsilon_w_kg ** (2.0 / 3.0)) * k_arr ** (-5.0 / 3.0)

    @staticmethod
    def model_spectrum_full(
        k: np.ndarray,
        epsilon_w_kg: float,
        L_m: float,
        eta_m: float,
        C_K: float = KOLMOGOROV_CONSTANT,
        p0: float = 4.0,
        beta: float = 5.2,
    ) -> np.ndarray:
        """Pope (2000) model spectrum: low-k roll-off + dissipative cut-off.

        E(k) = C_K ε^(2/3) k^(-5/3) · f_L(kL) · f_η(kη)

        with
            f_L(x) = (x / (x² + c_L)^(1/2))^(5/3 + p0)
            f_η(y) = exp(-β · y)

        Reproduces the textbook turbulence spectrum across the energy-
        containing, inertial, and dissipation ranges.
        """
        k = np.asarray(k, dtype=float)
        c_L = 6.78  # Pope (2000) tabulated constant
        x = k * L_m
        f_L = (x / np.sqrt(x * x + c_L)) ** (5.0 / 3.0 + p0)
        y = k * eta_m
        f_eta = np.exp(-beta * y)
        return (
            C_K
            * (epsilon_w_kg ** (2.0 / 3.0))
            * np.power(k, -5.0 / 3.0)
            * f_L
            * f_eta
        )

    @staticmethod
    def fit_inertial_slope(
        k: np.ndarray, E: np.ndarray, k_min: float, k_max: float
    ) -> float:
        """Fit log E vs log k slope inside [k_min, k_max].  Should give -5/3."""
        k = np.asarray(k); E = np.asarray(E)
        mask = (k >= k_min) & (k <= k_max) & (E > 0)
        if mask.sum() < 2:
            return float("nan")
        coef = np.polyfit(np.log(k[mask]), np.log(E[mask]), 1)
        return float(coef[0])

    # ---- Kármán vortex street -----------------------------------------

    def strouhal_number(
        self, shedding_freq_hz: float, diameter_m: float, U_m_s: float
    ) -> float:
        """St = f D / U.  For a circular cylinder St ≈ 0.21 in the
        Re ≈ 100 – 10⁵ regime."""
        return shedding_freq_hz * diameter_m / U_m_s

    @staticmethod
    def shedding_frequency(
        U_m_s: float, diameter_m: float, St: float = STROUHAL_CYLINDER
    ) -> float:
        """f = St · U / D using the canonical Strouhal number."""
        return St * U_m_s / diameter_m

    def vortex_shedding_field(
        self,
        U_m_s: float = 1.0,
        diameter_m: float = 0.1,
        n_pairs: int = 6,
        wake_offset_m: float = 0.06,
        cylinder_x_m: float | None = None,
        circulation_strength: float = 1.0,
        core_radius_m: float | None = None,
    ) -> dict:
        """Synthetic Kármán-street vorticity field for visualization.

        Builds a row of alternating-sign Lamb–Oseen vortices spaced according
        to the canonical Strouhal-number wavelength

            λ = U / f = D / St

        downstream of a cylinder centred at ``(cylinder_x_m, Ly/2)``.
        Returns a dict with the (X, Y, U, V, ω) fields and metadata for
        plotting.

        Parameters
        ----------
        U_m_s : float
            Free-stream velocity (mean horizontal flow).
        diameter_m : float
            Cylinder diameter D.
        n_pairs : int
            Number of vortex pairs (top + bottom) in the wake.
        wake_offset_m : float
            Vertical separation h between the two rows (matches the canonical
            value h/λ ≈ 0.281 for Kármán streets when set close to that).
        cylinder_x_m : float or None
            x-coordinate of the cylinder centre (default Lx/4).
        circulation_strength : float
            Peak vortex circulation Γ (in m²/s, positive).
        core_radius_m : float or None
            Lamb–Oseen core radius a.  Defaults to D/4.
        """
        geom = self.geometry
        X, Y = geom.meshgrid()
        if cylinder_x_m is None:
            cylinder_x_m = 0.25 * geom.Lx
        if core_radius_m is None:
            core_radius_m = 0.25 * diameter_m

        cyl_y = 0.5 * geom.Ly
        f_shed = self.shedding_frequency(U_m_s, diameter_m)
        wavelength = U_m_s / f_shed

        # Build vortex centres — alternating top / bottom downstream
        centres = []  # (x, y, sign)
        for i in range(n_pairs):
            xc_top = cylinder_x_m + 0.6 * diameter_m + (i + 0.0) * wavelength
            xc_bot = (
                cylinder_x_m + 0.6 * diameter_m + (i + 0.5) * wavelength
            )
            centres.append((xc_top, cyl_y + 0.5 * wake_offset_m, +1))
            centres.append((xc_bot, cyl_y - 0.5 * wake_offset_m, -1))

        # Free-stream flow
        u = np.full_like(X, U_m_s, dtype=float)
        v = np.zeros_like(Y, dtype=float)
        omega = np.zeros_like(X, dtype=float)

        # Add vortices (Lamb–Oseen velocity + Gaussian-core vorticity)
        Gamma = circulation_strength
        a2 = core_radius_m * core_radius_m
        for xc, yc, sign in centres:
            if xc > geom.Lx + 0.5 * wavelength:
                continue
            dx = X - xc
            dy = Y - yc
            r2 = dx * dx + dy * dy
            decay = 1.0 - np.exp(-r2 / a2)
            # tangential velocity v_theta = Γ/(2π r) (1 - exp(-r²/a²))
            denom = 2.0 * math.pi * (r2 + 1.0e-30)
            v_theta_over_r = sign * Gamma * decay / denom
            u += -dy * v_theta_over_r
            v += dx * v_theta_over_r
            # core vorticity: ω = (Γ / (π a²)) · exp(-r²/a²)
            omega += sign * (Gamma / (math.pi * a2)) * np.exp(-r2 / a2)

        # Mask out cylinder body itself
        cyl_mask = (X - cylinder_x_m) ** 2 + (Y - cyl_y) ** 2 < (
            0.5 * diameter_m
        ) ** 2
        u[cyl_mask] = 0.0
        v[cyl_mask] = 0.0
        omega[cyl_mask] = 0.0

        return {
            "X": X,
            "Y": Y,
            "U_field": u,
            "V_field": v,
            "vorticity": omega,
            "cylinder_xy": (cylinder_x_m, cyl_y),
            "diameter_m": diameter_m,
            "U_inf_m_s": U_m_s,
            "shedding_freq_hz": f_shed,
            "wavelength_m": wavelength,
            "cylinder_mask": cyl_mask,
        }

    # ---- Synthetic spectrum on the grid -------------------------------

    def synthetic_kolmogorov_field(
        self,
        epsilon_w_kg: float = 1.0,
        L_m: float | None = None,
        seed: int = 0,
        random: np.random.Generator | None = None,
    ) -> dict:
        """Generate a 2-D random velocity field with a k^(-5/3) shell spectrum.

        The field is constructed in Fourier space:

            û(k) = A · k^(-(5/3 + 1)/2) · exp(i φ_k)
            v̂(k) = solenoidal projection of û(k)

        so that the radially binned shell energy follows E(k) ∝ k^(-5/3).
        Useful for verification of the inertial-range slope without running a
        full DNS.
        """
        geom = self.geometry
        if L_m is None:
            L_m = 0.5 * min(geom.Lx, geom.Ly)
        rng = random if random is not None else np.random.default_rng(seed)

        kx = geom.kx()[None, :]                     # (1, Nx)
        ky = geom.ky()[:, None]                     # (Ny, 1)
        k_mag = np.sqrt(kx * kx + ky * ky)
        k_safe = np.where(k_mag > 0, k_mag, 1.0)

        # Amplitude that produces |û|² ~ k^(-5/3 - 1) so the radially
        # integrated shell energy E(k) = 2π k |û|² ~ k^(-5/3).
        amp = k_safe ** (-(5.0 / 3.0 + 1.0) / 2.0)
        # Energy-containing roll-off below k_L
        k_L = 2.0 * math.pi / L_m
        amp *= 1.0 / np.sqrt(1.0 + (k_L / k_safe) ** 4)
        # Random phase
        phi = rng.uniform(0.0, 2.0 * math.pi, size=k_mag.shape)
        u_hat = amp * np.exp(1j * phi)
        # Solenoidal projection (incompressibility)
        kdotu = (kx * u_hat) / k_safe   # use the same scalar amp for both comps
        u_hat_x = u_hat - (kx * kdotu) / k_safe
        u_hat_y = - (ky * kdotu) / k_safe
        u_hat_x[k_mag == 0] = 0.0
        u_hat_y[k_mag == 0] = 0.0

        u = np.fft.ifft2(u_hat_x).real
        v = np.fft.ifft2(u_hat_y).real

        # Renormalise so that the dissipation rate matches input ε
        u_rms = float(np.sqrt(np.mean(u * u + v * v)))
        if u_rms > 0:
            target_u_rms = (epsilon_w_kg * L_m) ** (1.0 / 3.0)
            scale = target_u_rms / u_rms
            u *= scale
            v *= scale
        return {
            "u": u,
            "v": v,
            "k_mag": k_mag,
            "kx": kx,
            "ky": ky,
            "L_m": L_m,
            "epsilon_w_kg": epsilon_w_kg,
        }

    @staticmethod
    def shell_average_spectrum(
        u: np.ndarray, v: np.ndarray, kx: np.ndarray, ky: np.ndarray,
        n_bins: int = 48,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Radially bin |û(k)|² + |v̂(k)|² into shell-averaged E(k)."""
        Ny, Nx = u.shape
        u_hat = np.fft.fft2(u) / (Nx * Ny)
        v_hat = np.fft.fft2(v) / (Nx * Ny)
        power = np.abs(u_hat) ** 2 + np.abs(v_hat) ** 2
        k_mag = np.sqrt(kx * kx + ky * ky)
        k_flat = k_mag.flatten()
        p_flat = power.flatten()
        k_min = float(np.min(k_flat[k_flat > 0]))
        k_max = float(np.max(k_flat))
        bin_edges = np.logspace(
            math.log10(k_min), math.log10(k_max), n_bins + 1
        )
        bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        E = np.zeros(n_bins)
        for i in range(n_bins):
            mask = (k_flat >= bin_edges[i]) & (k_flat < bin_edges[i + 1])
            if mask.sum() > 0:
                # Shell-averaged energy density: sum power / dk (per shell)
                dk = bin_edges[i + 1] - bin_edges[i]
                E[i] = p_flat[mask].sum() / max(dk, 1.0e-30)
        return bin_centres, E


__all__ = [
    # constants
    "KOLMOGOROV_CONSTANT",
    "KOLMOGOROV_EXPONENT",
    "RE_CRIT_PIPE",
    "RE_CRIT_PIPE_LOWER",
    "RE_CRIT_PIPE_UPPER",
    "STROUHAL_CYLINDER",
    "NU_AIR_M2_S",
    "NU_WATER_M2_S",
    "NU_HONEY_M2_S",
    "RHO_AIR_KGM3",
    "RHO_WATER_KGM3",
    # classes
    "TurbulenceGeometry",
    "TurbulenceSimulator",
]
