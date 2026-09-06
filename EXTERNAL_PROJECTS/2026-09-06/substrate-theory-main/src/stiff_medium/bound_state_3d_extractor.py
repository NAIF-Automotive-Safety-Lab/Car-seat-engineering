"""bound_state_3d_extractor.py
=================================

Paired GEOMETRY + SIMULATION + VISUALIZATION module that **extracts** a
3D substrate bound state directly from the field u(x, y, z, t).

A particle, in MODEL.md §2.5, is a localized strain pattern of the
substrate field whose internal cone-bouncing oscillation sets its rest
mass via

    m c² = ℏ ω_b,         ω_b² = (c/ξ)² + (γ / 2ρ)² .

This module shows that a particle EMERGES from the substrate dynamics:
seed the field with a kink-like radial bump centred at the origin,
evolve the full §18.11 Lagrangian, and read the rest energy back from
the breathing (radial) mode.

Components
----------

* ``BoundState3DGeometry`` — analytic kink-like radial bump, plus the
  cone-bouncing frequency and rest energy that the bump should breathe
  at.

* ``BoundState3DSimulator`` — wraps ``LatticeSubstrate3D``: seeds the
  bump, evolves the full sine-Gordon Lagrangian, and *measures*
  localization radius R(t), centre-amplitude u₀(t), and the breathing
  frequency by FFT of u₀(t).  Mass is then ℏ ω_meas / c² and rest
  energy E = ℏ ω_meas.

* Three canonical seeds — ground (radial bump), first-excited (one
  radial node), breathing-mode amplitude jitter.

Conventions
-----------

* Natural units (K = ρ = ξ = 1) are used by default so that
  c = √(K/ρ) = 1, ω_b = 1, and the analytic period T = 2π.  This keeps
  the simulator small (32³) and tests fast (<30 s).

* Optional SI conversion routines plug into ``cone_bouncing_protocol``
  and reproduce the electron rest mass m_e c² = ℏ · ω_b at the
  Compton anchor ξ = ξ_e.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, List

import numpy as np
from numpy.typing import NDArray

from src.stiff_medium.lattice_substrate_3d import LatticeSubstrate3D
from src.stiff_medium.cone_bouncing_protocol import (
    omega_b_canonical,
    mass_from_omega,
    HBAR,
    C,
    MEV,
    XI_E,
)


# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------


@dataclass
class BoundState3DGeometry:
    """Analytic 3D substrate-strain bump — the bound-state seed.

    The kink-like radial bump is

        u(r) = A · sech(r / ξ)                         (ground)

    which is the spherical analogue of the 1D sine-Gordon kink derivative
    (a localized soliton sitting in the cosine well at u = 0).  The bump
    is a *small-amplitude excitation* — its breathing frequency is the
    canonical cone-bouncing frequency:

        ω_b² = (c/ξ)² + (γ / 2ρ)²,    c ≡ √(K/ρ) .

    Parameters
    ----------
    K, rho, xi, gamma : float
        Substrate primitives.  Defaults are the natural-baseline
        K = ρ = ξ = 1, γ = 0 so c = 1, ω_b = 1, T = 2π.
    sigma_max : float
        Saturation cap (default 0.5 per §2.6).  Sets the peak amplitude
        of the bump as A = amplitude_frac · sigma_max.
    L : float
        Half-width of the cubic simulation box.  The bump must vanish
        well inside ±L.  Default 6.0 ξ.
    N : int
        Grid points along each axis (cubic Nx = Ny = Nz = N).  Default
        32 (8 192 / 32 768 cells, fast enough for the unit tests).
    """

    K: float = 1.0
    rho: float = 1.0
    xi: float = 1.0
    gamma: float = 0.0
    sigma_max: float = 0.5
    L: float = 12.0
    N: int = 32

    # ---- derived ----
    @property
    def c_eff(self) -> float:
        """Substrate wave speed √(K/ρ)."""
        return float(np.sqrt(self.K / self.rho))

    @property
    def dx(self) -> float:
        """Grid spacing 2L / N."""
        return float(2.0 * self.L / self.N)

    @property
    def dt(self) -> float:
        """Time step set by the 3D CFL with safety factor 0.4."""
        return float(0.4 * self.dx / (self.c_eff * np.sqrt(3.0)))

    @property
    def omega_b(self) -> float:
        """Canonical cone-bouncing angular frequency [rad/s or rad/τ]."""
        return omega_b_canonical(self.K, self.rho, self.xi, self.gamma)

    @property
    def period(self) -> float:
        """Bounce / breathing period T = 2π / ω_b."""
        omega = self.omega_b
        return float("inf") if omega <= 0 else float(2.0 * np.pi / omega)

    @property
    def mass_kg(self) -> float:
        """Rest mass m = ℏ ω_b / c² (SI when primitives are SI)."""
        return mass_from_omega(self.omega_b)

    @property
    def rest_energy(self) -> float:
        """Rest energy E_rest = ℏ ω_b (SI when primitives are SI;
        in natural units this is just ω_b in dimensionless form)."""
        return float(HBAR * self.omega_b)

    # ---- seed profiles ----

    @property
    def envelope_width(self) -> float:
        """Width σ_env of the Gaussian envelope used for bound-state seeds.

        Chosen wide compared to ξ so that the wavevector content of the
        bump satisfies |k| ≪ 1/ξ.  In that limit the small-amplitude
        oscillation frequency reduces cleanly to the cone-bouncing
        baseline ω₀ = c/ξ (the rest-frequency of the §18.11 cosine
        well's small-amplitude expansion).  Default = 4 ξ.
        """
        return float(4.0 * self.xi)

    def ground_bump(
        self,
        X: NDArray[np.float64],
        Y: NDArray[np.float64],
        Z: NDArray[np.float64],
        amplitude_frac: float = 0.05,
    ) -> NDArray[np.float64]:
        """Ground-mode small-amplitude radial Gaussian bump.

            u(r) = A · exp(−r² / (2 σ_env²))

        Centred at the origin, spherically symmetric.  Amplitude
        A = amplitude_frac · sigma_max is intentionally small so the
        cosine potential linearises to a quadratic well; the bump then
        oscillates vertically (centre amplitude up/down) at the rest
        frequency ω₀ = c/ξ.

        At amplitude_frac → 0 this is the cleanest realisation of the
        cone-bouncing rest mode; finite amplitude introduces small
        anharmonic corrections.
        """
        if not 0.0 < amplitude_frac < 1.0:
            raise ValueError("amplitude_frac must lie in (0, 1).")
        A = amplitude_frac * self.sigma_max
        r2 = X * X + Y * Y + Z * Z
        return A * np.exp(-r2 / (2.0 * self.envelope_width ** 2))

    def first_excited_bump(
        self,
        X: NDArray[np.float64],
        Y: NDArray[np.float64],
        Z: NDArray[np.float64],
        amplitude_frac: float = 0.05,
    ) -> NDArray[np.float64]:
        """First radial-excited bump:  A · (1 − r²/σ²) · exp(−r²/(2σ²)).

        Has a spherical node at r = σ_env (the 3D analogue of the
        first excited harmonic-oscillator radial mode, ψ_1 ∝ (1 − ρ²) e^(-ρ²/2)).
        Width and amplitude inherit from the ground mode so the FFT
        peak still sits at ω₀ = c/ξ.
        """
        if not 0.0 < amplitude_frac < 1.0:
            raise ValueError("amplitude_frac must lie in (0, 1).")
        A = amplitude_frac * self.sigma_max
        r2 = X * X + Y * Y + Z * Z
        s2 = self.envelope_width ** 2
        return A * (1.0 - r2 / s2) * np.exp(-r2 / (2.0 * s2))

    def breathing_bump(
        self,
        X: NDArray[np.float64],
        Y: NDArray[np.float64],
        Z: NDArray[np.float64],
        amplitude_frac: float = 0.05,
        breath_frac: float = 0.10,
    ) -> NDArray[np.float64]:
        """Breathing-mode seed: ground bump amplified by (1 + δ).

        Re-using the matched envelope width but with the *amplitude*
        scaled up by ``breath_frac`` so the centre cell oscillates with
        the desired AC-DC ratio while keeping the fundamental frequency
        ω₀ = c/ξ unchanged.

        This is the simplest seed that exercises the time-evolved
        centre amplitude well above lattice-noise floor without
        introducing spurious wavevector content (which would shift
        the FFT peak above ω₀).
        """
        if not 0.0 < amplitude_frac < 1.0:
            raise ValueError("amplitude_frac must lie in (0, 1).")
        if abs(breath_frac) >= 1.0:
            raise ValueError("|breath_frac| must be < 1.")
        A = amplitude_frac * self.sigma_max * (1.0 + breath_frac)
        r2 = X * X + Y * Y + Z * Z
        return A * np.exp(-r2 / (2.0 * self.envelope_width ** 2))


# ---------------------------------------------------------------------------
# 2. Simulator
# ---------------------------------------------------------------------------


@dataclass
class BoundState3DSimulator:
    """Time-evolve a 3D substrate field carrying a bound-state bump.

    Wraps ``LatticeSubstrate3D`` with the full §18.11 sine-Gordon
    Lagrangian and exposes shape / energy / frequency extractors that
    treat the field AS the particle.

    Attributes
    ----------
    geometry : BoundState3DGeometry
        Source of K, ρ, ξ, γ, σ_max plus the analytic seeds.
    mode : str
        One of ``'ground'``, ``'excited'``, ``'breathing'``.  Sets
        which seed is used by ``run``.
    n_periods : float
        Number of analytic periods to evolve (default 4).  The FFT
        frequency-extractor needs ≥ 2 periods to resolve the peak.
    samples_per_period : int
        Number of stored snapshots per analytic period (default 32).
        Must obey Nyquist for the FFT.
    """

    geometry: BoundState3DGeometry = field(default_factory=BoundState3DGeometry)
    mode: str = "breathing"
    n_periods: float = 4.0
    samples_per_period: int = 32

    # populated by run()
    _result: Optional[Dict[str, NDArray[np.float64]]] = field(default=None, init=False, repr=False)

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def _make_lattice(self) -> LatticeSubstrate3D:
        g = self.geometry
        sub = LatticeSubstrate3D(
            Nx=g.N, Ny=g.N, Nz=g.N,
            dx=g.dx, dt=g.dt,
            rho=g.rho, K=g.K, xi=g.xi, gamma=g.gamma,
            bc="periodic",
        )
        return sub

    def seed(self, sub: LatticeSubstrate3D) -> None:
        """Apply the configured analytic seed in-place to ``sub.u``."""
        g = self.geometry
        if self.mode == "ground":
            sub.u = g.ground_bump(sub.X, sub.Y, sub.Z)
        elif self.mode == "excited":
            sub.u = g.first_excited_bump(sub.X, sub.Y, sub.Z)
        elif self.mode == "breathing":
            sub.u = g.breathing_bump(sub.X, sub.Y, sub.Z)
        else:
            raise ValueError(
                f"mode must be 'ground' | 'excited' | 'breathing'; got {self.mode!r}"
            )
        sub.v.fill(0.0)

    # ------------------------------------------------------------------
    # Time evolution + extraction
    # ------------------------------------------------------------------

    def run(self) -> Dict[str, NDArray[np.float64]]:
        """Evolve the field, sampling diagnostics every Δt_sample.

        Returns
        -------
        dict with keys
            ``t``         : time samples [period units when natural]
            ``u_center``  : centre-cell strain u(0,0,0,t)
            ``u_max``     : max(|u|) over the box
            ``radius``    : second-moment radius √(⟨r²⟩_|u|)
            ``energy``    : total Lagrangian energy
            ``snapshots`` : list of full u(x,y,z) at t = 0, T/2, T
        """
        g = self.geometry
        sub = self._make_lattice()
        self.seed(sub)

        T_anal = g.period
        T_total = self.n_periods * T_anal
        n_samples = int(self.n_periods * self.samples_per_period) + 1
        dt_sample = T_total / (n_samples - 1)
        steps_between = max(1, int(round(dt_sample / sub.dt)))

        t_arr     = np.empty(n_samples, dtype=float)
        u_center  = np.empty(n_samples, dtype=float)
        u_max     = np.empty(n_samples, dtype=float)
        radius    = np.empty(n_samples, dtype=float)
        energy    = np.empty(n_samples, dtype=float)
        snapshots: List[NDArray[np.float64]] = []
        snap_times = [0.0, 0.5 * T_anal, 1.0 * T_anal]
        snap_taken = [False, False, False]

        cx, cy, cz = sub.Nx // 2, sub.Ny // 2, sub.Nz // 2
        # |r|² weight for the second-moment radius
        R2 = sub.X ** 2 + sub.Y ** 2 + sub.Z ** 2

        for k in range(n_samples):
            t_arr[k]    = sub.t
            u_center[k] = float(sub.u[cx, cy, cz])
            absu        = np.abs(sub.u)
            u_max[k]    = float(absu.max())
            mass        = absu.sum()
            radius[k]   = float(np.sqrt((R2 * absu).sum() / max(mass, 1e-30)))
            energy[k]   = float(sub.total_energy())

            # Take canonical snapshots at t = 0, T/2, T
            for s_idx, target in enumerate(snap_times):
                if not snap_taken[s_idx] and sub.t >= target - 0.5 * sub.dt:
                    snapshots.append(sub.u.copy())
                    snap_taken[s_idx] = True

            if k < n_samples - 1:
                sub.evolve(steps_between)

        # Backstop: ensure 3 snapshots even if loop missed one
        while len(snapshots) < 3:
            snapshots.append(sub.u.copy())

        self._result = {
            "t":         t_arr,
            "u_center":  u_center,
            "u_max":     u_max,
            "radius":    radius,
            "energy":    energy,
            "snapshots": np.asarray(snapshots[:3]),
            "x":         sub.x,
            "y":         sub.y,
            "z":         sub.z,
            "dt":        np.array([sub.dt]),
            "dx":        np.array([sub.dx]),
            "T_period":  np.array([T_anal]),
        }
        return self._result

    # ------------------------------------------------------------------
    # Extractors
    # ------------------------------------------------------------------

    @property
    def result(self) -> Dict[str, NDArray[np.float64]]:
        if self._result is None:
            self.run()
        assert self._result is not None
        return self._result

    def measured_omega(self) -> float:
        """FFT the centre-cell strain trace and return the peak ω.

        Uses a Hann window and zero-DC scrub.  For the breathing seed
        the peak should sit within ~10% of the analytic ω_b.
        """
        res = self.result
        u = res["u_center"]
        t = res["t"]
        if len(u) < 8:
            return float("nan")
        # De-mean and apply Hann window
        u0 = u - u.mean()
        win = np.hanning(len(u0))
        spec = np.fft.rfft(u0 * win)
        freqs_hz = np.fft.rfftfreq(len(u0), d=float(t[1] - t[0]))
        omega = 2.0 * np.pi * freqs_hz
        power = np.abs(spec) ** 2
        # Skip DC bin
        if len(power) <= 1:
            return float("nan")
        idx = 1 + int(np.argmax(power[1:]))
        # Quadratic interpolation around the peak for sub-bin accuracy
        if 1 < idx < len(power) - 1:
            y0, y1, y2 = power[idx - 1], power[idx], power[idx + 1]
            denom = (y0 - 2.0 * y1 + y2)
            shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        else:
            shift = 0.0
        omega_peak = omega[idx] + shift * (omega[1] - omega[0])
        return float(omega_peak)

    def measured_mass_kg(self) -> float:
        """Rest mass m = ℏ ω_meas / c²."""
        return mass_from_omega(self.measured_omega())

    def measured_rest_energy(self) -> float:
        """E_rest = ℏ ω_meas (joules, SI when primitives are SI)."""
        return float(HBAR * self.measured_omega())

    def localization_metric(self) -> Dict[str, float]:
        """Summary of how well-localized the bound state stays.

        Returns dict with:
          ``r_init``, ``r_max``, ``r_mean``, ``r_growth_pct``.
        """
        res = self.result
        r = res["radius"]
        r_init = float(r[0])
        r_max  = float(r.max())
        r_mean = float(r.mean())
        r_growth_pct = 100.0 * (r_max - r_init) / max(r_init, 1e-30)
        return {
            "r_init":       r_init,
            "r_max":        r_max,
            "r_mean":       r_mean,
            "r_growth_pct": r_growth_pct,
        }

    def report(self) -> Dict[str, float]:
        """One-shot diagnostic dict — measured ω, predicted ω_b, error."""
        omega_meas = self.measured_omega()
        omega_pred = self.geometry.omega_b
        rel = abs(omega_meas - omega_pred) / max(omega_pred, 1e-30)
        loc = self.localization_metric()
        return {
            "omega_pred":      omega_pred,
            "omega_meas":      omega_meas,
            "omega_rel_err":   rel,
            "rest_E_pred":     float(HBAR * omega_pred),
            "rest_E_meas":     float(HBAR * omega_meas),
            "mass_pred_kg":    mass_from_omega(omega_pred),
            "mass_meas_kg":    mass_from_omega(omega_meas),
            "r_init":          loc["r_init"],
            "r_max":           loc["r_max"],
            "r_growth_pct":    loc["r_growth_pct"],
        }


# ---------------------------------------------------------------------------
# 3. Multi-mode helper for visualisation
# ---------------------------------------------------------------------------


def run_three_modes(
    geometry: Optional[BoundState3DGeometry] = None,
    n_periods: float = 3.0,
    samples_per_period: int = 24,
) -> Dict[str, Dict[str, NDArray[np.float64]]]:
    """Run all three canonical modes side-by-side.

    Returns a dict keyed by ``'ground' | 'excited' | 'breathing'``,
    each value being the corresponding ``Simulator.result`` dict.
    Used by the visualiser to draw the side-by-side modes panel.
    """
    geom = geometry or BoundState3DGeometry()
    out: Dict[str, Dict[str, NDArray[np.float64]]] = {}
    for mode in ("ground", "excited", "breathing"):
        sim = BoundState3DSimulator(
            geometry=geom, mode=mode,
            n_periods=n_periods, samples_per_period=samples_per_period,
        )
        out[mode] = sim.run()
        out[mode]["_omega_meas"] = np.array([sim.measured_omega()])
        out[mode]["_omega_pred"] = np.array([geom.omega_b])
    return out


# ---------------------------------------------------------------------------
# 4. Anchor convenience: electron primitives
# ---------------------------------------------------------------------------


def electron_geometry() -> BoundState3DGeometry:
    """Return a geometry whose bound state is the electron at the
    Compton anchor ξ = ξ_e, K = ρ c² so c_eff = c, γ = 0.

    This is for SI cross-checks; the simulator itself runs natively in
    natural units by default for speed.
    """
    rho = 1.0
    K = rho * C * C
    return BoundState3DGeometry(K=K, rho=rho, xi=XI_E, gamma=0.0,
                                 sigma_max=0.5, L=12.0 * XI_E, N=24)


def electron_rest_energy_mev(geom: Optional[BoundState3DGeometry] = None) -> float:
    """E_rest of the electron geometry expressed in MeV — should match
    the PDG m_e c² = 0.510999 MeV exactly when geom uses ξ_e and γ=0."""
    g = geom or electron_geometry()
    return float(HBAR * g.omega_b / MEV)
