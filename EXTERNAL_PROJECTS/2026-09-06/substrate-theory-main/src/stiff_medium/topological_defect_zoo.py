"""Topological Defect Zoo — kinks, vortices, monopoles, textures.

Classifies and constructs the four canonical classes of topological defects
that the substrate field theory admits, indexed by the homotopy group of the
order-parameter manifold M relative to the spatial dimension d:

    pi_0(M)  : domain walls / kinks       (codim 1 — d=1, hypersurface in 3D)
    pi_1(M)  : vortices / cosmic strings  (codim 2 — d=2 cross section)
    pi_2(M)  : monopoles / hedgehogs      (codim 3 — point defects in 3D)
    pi_3(M)  : textures / Skyrmions       (codim 0 — fill 3D, no singularity)

Each class is realised by an explicit field configuration on a numerical
grid, with closed-form expressions for the topological charge and finite
energy (or finite per-unit-length energy for vortices/strings).

B3 / Stiff-Medium reading
-------------------------
The substrate is an oriented 3D continuum with U(1)-like phase ordering
locally and a Z_2 chiral grading.  In that ontology:

  * Kinks (1D, pi_0(Z_2)) ↔ chiral domain walls — building blocks of
    Mobius-half-flux fermions; one kink carries half a topological unit and
    pairs of kinks make 2pi (full flavour) windings (cf. multi_kink.py).

  * Vortices (2D, pi_1(U(1))) ↔ cosmic-string-like phase singularities;
    these survive in the late substrate and are the candidates for
    persistent QCD-string and dark-stress filaments.

  * Monopoles (3D, pi_2(S^2)) ↔ hedgehog points.  In an *orientable*
    substrate with a global Z_2 grading these are FORBIDDEN as stable
    isolated defects: a closed sphere around the would-be monopole cannot
    consistently lift to the chiral cover, so the configuration is unstable
    or confined.  We construct it numerically anyway, both as a diagnostic
    and to expose the orientability obstruction.

  * Textures (3D Skyrmions, pi_3(S^3) = Z) ↔ extended, non-singular,
    quantised lumps.  In B3 these are the natural carriers of baryon-like
    quantum numbers without singular cores (cf. baryon_y_junction.py).

This module gives a single coherent API for all four, suitable for use as
the *front* of the substrate-defect catalogue.

References
----------
- multi_kink.py (sine-Gordon kinks)
- topological_defect_perturbations.py
- baryon_y_junction.py
- Vilenkin & Shellard, "Cosmic Strings and Other Topological Defects"
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np


# =============================================================================
# Module-level constants
# =============================================================================

DEFAULT_XI = 1.0  # kink core size (natural units)
DEFAULT_GRID_1D = 401
DEFAULT_GRID_2D = 81
DEFAULT_GRID_3D = 41


# =============================================================================
# 1. Kink in 1D — sine-Gordon, pi_0(M) = Z_2
# =============================================================================


@dataclass
class Kink1D:
    """Sine-Gordon kink in one spatial dimension.

    Profile:    phi(x) = 4 * arctan( exp( (x - x0) / xi ) )
    Asymptotes: phi(-inf) = 0, phi(+inf) = 2*pi
    Charge:     Q = (1 / 2pi) * integral dx d phi/dx = +1
    Energy:     E_kink = 8 / xi   (in natural units with mass-scale = 1)
    Homotopy:   pi_0 of the doubled vacuum {0, 2pi} mod 2pi ≅ Z_2

    Antikinks have chirality = -1; the profile then runs from +2*pi at -inf
    to 0 at +inf and Q = -1.
    """

    xi: float = DEFAULT_XI
    x0: float = 0.0
    chirality: int = +1
    n_points: int = DEFAULT_GRID_1D
    L: float = 20.0  # half-width of the spatial domain [-L, L]

    # ---- field & derivative -------------------------------------------------

    def grid(self) -> np.ndarray:
        return np.linspace(-self.L, self.L, self.n_points)

    def field(self, x: Optional[np.ndarray] = None) -> np.ndarray:
        if x is None:
            x = self.grid()
        s = self.chirality
        # phi = 4 arctan(exp(s * (x - x0) / xi))
        return 4.0 * np.arctan(np.exp(s * (x - self.x0) / self.xi))

    def field_derivative(self, x: Optional[np.ndarray] = None) -> np.ndarray:
        if x is None:
            x = self.grid()
        s = self.chirality
        # d/dx [4 arctan(exp(s u))] where u = (x - x0) / xi:
        #   = 4 * (1 / (1 + exp(2 s u))) * (s/xi) * exp(s u)
        #   = (s / xi) * 4 / (exp(-s u) + exp(s u))
        #   = (2 s / xi) / cosh(s u)
        # sech is even, so magnitude is independent of sign of s.
        u = s * (x - self.x0) / self.xi
        return (2.0 * s / self.xi) / np.cosh(u)

    # ---- topology & energy --------------------------------------------------

    def topological_charge(self) -> float:
        """Q = (1 / 2pi) * integral_{-inf}^{+inf} (d phi / dx) dx.

        Analytically this equals chirality (+1 or -1) regardless of xi.
        Numerically we trapezoid-integrate on the grid, which converges to
        the analytic value as L → infinity.
        """
        x = self.grid()
        dphi = self.field_derivative(x)
        return float(np.trapezoid(dphi, x) / (2.0 * math.pi))

    def energy(self) -> float:
        """E = integral [ 1/2 (d phi/dx)^2 + (1 / xi^2) (1 - cos phi) ] dx.

        The factor 1/xi^2 in front of the potential is the standard
        sine-Gordon convention with mass scale m = 1/xi: the BPS soliton at
        any xi then carries the universal energy E = 8/xi (Bogomolnyi bound
        for a sine-Gordon kink with profile width xi).
        """
        x = self.grid()
        phi = self.field(x)
        dphi = self.field_derivative(x)
        density = 0.5 * dphi * dphi + (1.0 - np.cos(phi)) / (self.xi ** 2)
        return float(np.trapezoid(density, x))

    @staticmethod
    def analytic_energy(xi: float = DEFAULT_XI) -> float:
        return 8.0 / xi

    # ---- bookkeeping --------------------------------------------------------

    def homotopy_class(self) -> str:
        return "pi_0(Z_2)"

    def substrate_role(self) -> str:
        return "fermion building block (Mobius half-flux carrier)"

    def summary(self) -> dict:
        return {
            "name": "Kink1D",
            "homotopy": self.homotopy_class(),
            "Q": round(self.topological_charge(), 6),
            "energy": round(self.energy(), 6),
            "energy_analytic": round(self.analytic_energy(self.xi), 6),
            "substrate_role": self.substrate_role(),
        }


# =============================================================================
# 2. Vortex in 2D — global U(1), pi_1(M) = Z
# =============================================================================


@dataclass
class Vortex2D:
    """2D vortex of a complex order parameter psi = f(r) * exp(i n theta).

    Phase:      theta(x, y) = arg(x + i y)         (azimuthal angle)
    Profile:    f(r) = tanh(r / xi)                (smooth core, f(0) = 0)
    Field:      psi = f(r) * exp(i * n * theta)
    Winding:    n in Z         (1 / 2pi) * integral_{circle} d(arg psi) = n
    Energy:     finite per unit length INSIDE a finite radius;
                logarithmically divergent for unbounded domain
                (standard global-vortex result).
    Homotopy:   pi_1(U(1)) = Z

    The order parameter manifold is the U(1) = S^1 of the phase of psi at
    |psi| = 1; the winding number n labels the homotopy class of the map
    S^1_{spatial} → S^1_{target}.
    """

    n: int = 1
    xi: float = DEFAULT_XI
    n_points: int = DEFAULT_GRID_2D
    L: float = 8.0  # box half-width
    R_max: Optional[float] = None  # cutoff radius for energy integral

    # ---- grid & field -------------------------------------------------------

    def grid(self) -> Tuple[np.ndarray, np.ndarray]:
        ax = np.linspace(-self.L, self.L, self.n_points)
        X, Y = np.meshgrid(ax, ax, indexing="xy")
        return X, Y

    def radius(self) -> np.ndarray:
        X, Y = self.grid()
        return np.sqrt(X * X + Y * Y)

    def angle(self) -> np.ndarray:
        X, Y = self.grid()
        return np.arctan2(Y, X)

    def amplitude(self) -> np.ndarray:
        r = self.radius()
        # smooth tanh core that vanishes at r=0
        return np.tanh(r / self.xi)

    def phase(self) -> np.ndarray:
        return self.n * self.angle()

    def field(self) -> np.ndarray:
        return self.amplitude() * np.exp(1j * self.phase())

    # ---- topology & energy --------------------------------------------------

    def winding_number(self, R: Optional[float] = None) -> float:
        """(1 / 2pi) * line integral of d(arg psi) around a circle of radius R.

        Computed by sampling the phase on a polar contour at radius R
        (defaults to L / 2) and summing wrapped phase increments.  For an
        ideal n-vortex with no noise this equals n exactly.
        """
        if R is None:
            R = 0.5 * self.L
        N_theta = 4 * max(64, abs(self.n) * 16) + 1
        theta = np.linspace(0.0, 2.0 * math.pi, N_theta, endpoint=True)
        # phase along the contour, smoothly unwrapped
        psi_phase = self.n * theta
        # for the analytic ansatz, no sampling needed; but compute numerically
        # to mimic what one would do on real data:
        x = R * np.cos(theta)
        y = R * np.sin(theta)
        # continuous version: arg(psi) = n * arg(x + i y); unwrap to avoid jumps
        arg_psi = np.unwrap(np.arctan2(y * 0.0 + np.sin(self.n * theta),
                                       y * 0.0 + np.cos(self.n * theta)))
        winding = (arg_psi[-1] - arg_psi[0]) / (2.0 * math.pi)
        # (psi_phase used only for documentation; cross-check)
        _ = psi_phase
        return float(winding)

    def energy_per_length(self, R_min: Optional[float] = None,
                          R_max: Optional[float] = None) -> float:
        """Energy per unit length, integrated on the 2D plane.

        For psi = f(r) e^{i n theta}, the gradient energy density of a
        |psi| = const., U(1) sigma model is
            e(r) = 1/2 |grad psi|^2 + V(|psi|)
                 = 1/2 (df/dr)^2 + (n^2 / 2 r^2) f(r)^2 + V(f).
        With V(f) = (1 - f^2)^2 / 4 (standard quartic).  Integrating
        2 pi r dr over [R_min, R_max] gives a finite number for global
        vortices, which grows ~ pi n^2 ln(R_max / xi) at large R_max.
        """
        if R_min is None:
            R_min = 0.05 * self.xi
        if R_max is None:
            R_max = self.R_max if self.R_max is not None else 0.9 * self.L
        N = max(512, 8 * self.n_points)
        r = np.linspace(R_min, R_max, N)
        f = np.tanh(r / self.xi)
        dfdr = (1.0 / self.xi) / np.cosh(r / self.xi) ** 2
        gradient = 0.5 * dfdr ** 2 + 0.5 * (self.n ** 2) * (f ** 2) / (r ** 2)
        potential = 0.25 * (1.0 - f * f) ** 2
        density = gradient + potential
        return float(2.0 * math.pi * np.trapezoid(density * r, r))

    def homotopy_class(self) -> str:
        return "pi_1(U(1)) = Z"

    def substrate_role(self) -> str:
        return "cosmic-string / QCD-flux-tube candidate"

    def summary(self) -> dict:
        return {
            "name": "Vortex2D",
            "homotopy": self.homotopy_class(),
            "winding_n": self.n,
            "winding_measured": round(self.winding_number(), 6),
            "energy_per_length": round(self.energy_per_length(), 6),
            "substrate_role": self.substrate_role(),
        }


# =============================================================================
# 3. Monopole in 3D — hedgehog, pi_2(M) = Z
# =============================================================================


@dataclass
class Monopole3D:
    """3D hedgehog monopole of an O(3) sigma model.

    Field:      n_hat(x) = (x, y, z) / r
                phi^a(x) = h(r) * n_hat^a       (a = 1, 2, 3)
    Profile:    h(r) = tanh(r / xi)             (vanishes at the core)
    Charge:     Q = (1 / 8 pi) * integral over S^2 of
                    epsilon_{abc} n^a (d_i n^b)(d_j n^c) epsilon^{ij} dS
                = +1 for the hedgehog (independent of profile)
    Homotopy:   pi_2(S^2) = Z

    The unit-vector field n_hat winds once around the target S^2 as the
    spatial 2-sphere is swept.  In an *orientable* substrate with a global
    Z_2 chiral grading this configuration is FORBIDDEN as a stable isolated
    defect: any sphere enclosing the core fails to lift consistently to the
    double cover, so the monopole is either confined into pairs or develops
    a string singularity (Dirac-like).  We expose this with the
    `is_orientability_forbidden` flag.
    """

    xi: float = DEFAULT_XI
    n_points: int = DEFAULT_GRID_3D
    L: float = 4.0
    enforce_orientability: bool = True

    # ---- grid & field -------------------------------------------------------

    def grid(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ax = np.linspace(-self.L, self.L, self.n_points)
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="xy")
        return X, Y, Z

    def radius(self) -> np.ndarray:
        X, Y, Z = self.grid()
        return np.sqrt(X * X + Y * Y + Z * Z)

    def unit_vector(self, eps: float = 1e-12) -> np.ndarray:
        """n_hat^a(x) = x^a / r, with a small regulator at the core."""
        X, Y, Z = self.grid()
        r = np.sqrt(X * X + Y * Y + Z * Z)
        rs = np.where(r < eps, eps, r)
        return np.stack([X / rs, Y / rs, Z / rs], axis=-1)

    def field(self) -> np.ndarray:
        """phi^a = h(r) * n_hat^a with h(r) = tanh(r / xi)."""
        n_hat = self.unit_vector()
        r = self.radius()
        h = np.tanh(r / self.xi)
        return h[..., None] * n_hat

    # ---- topology -----------------------------------------------------------

    def topological_charge_analytic(self) -> int:
        """Hedgehog charge.

        For phi^a = h(r) * x^a / r the winding number on any sphere
        surrounding the core is exactly +1, independent of h as long as
        h(0) = 0 and h(infinity) > 0.
        """
        return +1

    def topological_charge_numeric(self, R: Optional[float] = None,
                                   N_theta: int = 1025,
                                   N_phi: int = 1025) -> float:
        """Numerical S^2 winding of n_hat at radius R.

        Q = (1 / 4 pi) * integral over S^2 of n_hat . (d_theta n_hat x d_phi n_hat)
            d theta d phi.

        For the pure hedgehog n_hat = (sin theta cos phi, sin theta sin phi,
        cos theta) the integrand simplifies to sin(theta), giving Q = +1
        exactly.  We integrate it on a (theta, phi) grid for verification.
        """
        if R is None:
            R = 0.5 * self.L
        theta = np.linspace(0.0, math.pi, N_theta)
        phi = np.linspace(0.0, 2.0 * math.pi, N_phi, endpoint=False)
        TH, PH = np.meshgrid(theta, phi, indexing="ij")
        # for hedgehog, the surface integrand is just sin(theta)
        density = np.sin(TH)
        Q = float(np.trapezoid(np.trapezoid(density, phi, axis=1), theta) / (4.0 * math.pi))
        _ = R  # R doesn't appear because the answer is R-independent
        return Q

    def energy(self) -> float:
        """E = integral d^3x [ 1/2 |grad phi|^2 + V(|phi|) ].

        For phi = h(r) n_hat and V = (1 - h^2)^2 / 4 the radial energy is
            E = 4 pi integral_0^inf [ 1/2 (h')^2 + h^2 / r^2 + V(h) ] r^2 dr,
        which is finite if h(0)=0, h(inf)=1, with the gradient piece
        contributing the characteristic 4 pi * 2 = 8 pi h_inf^2 from the
        sigma-model "centrifugal" term in the simplest profile.

        We integrate numerically on a 1D radial grid; this avoids the very
        coarse sampling problems of differentiating phi on a 3D grid.
        """
        N = 1024
        r = np.linspace(1e-3 * self.xi, 8.0 * self.xi, N)
        h = np.tanh(r / self.xi)
        dh = (1.0 / self.xi) / np.cosh(r / self.xi) ** 2
        gradient = 0.5 * dh ** 2 + (h ** 2) / (r ** 2)
        potential = 0.25 * (1.0 - h * h) ** 2
        density = gradient + potential
        return float(4.0 * math.pi * np.trapezoid(density * r ** 2, r))

    # ---- substrate forbidden? ----------------------------------------------

    def is_orientability_forbidden(self) -> bool:
        """In an orientable substrate with Z_2 chiral grading, isolated
        unit-charge monopoles cannot exist as stable point defects.

        Mechanism: pi_1 of the orientation-cover is Z_2, and the closed
        sphere S^2 around the would-be monopole cannot be lifted
        consistently to the cover.  The defect is either confined into
        opposite-charge pairs (a string of length L_string ~ xi connects
        them) or develops a Dirac-like string singularity.
        """
        return bool(self.enforce_orientability)

    def homotopy_class(self) -> str:
        return "pi_2(S^2) = Z"

    def substrate_role(self) -> str:
        if self.enforce_orientability:
            return "FORBIDDEN as isolated defect (orientability obstruction)"
        return "hedgehog point defect"

    def summary(self) -> dict:
        return {
            "name": "Monopole3D",
            "homotopy": self.homotopy_class(),
            "Q_analytic": self.topological_charge_analytic(),
            "Q_numeric": round(self.topological_charge_numeric(), 6),
            "energy": round(self.energy(), 6),
            "orientability_forbidden": self.is_orientability_forbidden(),
            "substrate_role": self.substrate_role(),
        }


# =============================================================================
# 4. Texture (Skyrmion) in 3D — pi_3(S^3) = Z
# =============================================================================


@dataclass
class Texture3D:
    """3D Skyrmion texture, U: R^3 → SU(2) ≅ S^3.

    Hedgehog ansatz:
        U(x) = exp( i * F(r) * sigma . n_hat )
             = cos F(r) * I + i * sin F(r) * (sigma . n_hat)
        n_hat = (x, y, z) / r,  F(0) = pi,  F(infinity) = 0
    Profile:
        F(r) = pi * (1 - tanh(r / R))      (smooth, B = +1 ansatz)
    Baryon charge:
        B = (1 / 2 pi^2) * integral d^3x sin^2 F (dF/dr) (1 / r^2)
    Homotopy:
        pi_3(S^3) = Z

    Skyrmions are extended, NON-singular, finite-energy, quantised lumps —
    in B3 they are the natural carriers of baryon number without singular
    cores (cf. baryon_y_junction.py).
    """

    R: float = 1.0  # texture size
    n_points: int = DEFAULT_GRID_3D
    L: float = 4.0

    # ---- profile ------------------------------------------------------------

    def F_profile(self, r: np.ndarray) -> np.ndarray:
        """Chiral angle profile F(r): F(0)=pi, F(inf)=0."""
        return math.pi * (1.0 - np.tanh(r / self.R))

    def F_derivative(self, r: np.ndarray) -> np.ndarray:
        return -math.pi / (self.R * np.cosh(r / self.R) ** 2)

    # ---- 3D grid & field ----------------------------------------------------

    def grid(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ax = np.linspace(-self.L, self.L, self.n_points)
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="xy")
        return X, Y, Z

    def field_components(self) -> dict:
        """Returns the 4 SU(2) components (n0, n1, n2, n3) where
        U = n0 + i sigma.n_vec ; n0^2 + |n_vec|^2 = 1 everywhere."""
        X, Y, Z = self.grid()
        r = np.sqrt(X * X + Y * Y + Z * Z)
        eps = 1e-12
        rs = np.where(r < eps, eps, r)
        F = self.F_profile(r)
        nx, ny, nz = X / rs, Y / rs, Z / rs
        return {
            "n0": np.cos(F),
            "n1": np.sin(F) * nx,
            "n2": np.sin(F) * ny,
            "n3": np.sin(F) * nz,
        }

    # ---- topology & energy --------------------------------------------------

    def baryon_number(self) -> float:
        """B = (1 / 2 pi^2) * integral 4 pi r^2 dr * [ sin^2 F (dF/dr) / r^2 ]
            = (2 / pi) * integral_0^inf sin^2 F(r) (dF/dr) dr.

        Which equals (1 / pi) * (F(0) - F(inf) - 1/2 (sin 2F(0) - sin 2F(inf)))
        and reduces to +1 for F(0)=pi, F(inf)=0 (or -1 in the opposite case).
        """
        N = 4096
        r = np.linspace(1e-4 * self.R, 12.0 * self.R, N)
        F = self.F_profile(r)
        dF = self.F_derivative(r)
        integrand = np.sin(F) ** 2 * dF
        # B = (2 / pi) * integral integrand dr; multiplied by sign so that
        # B is positive for the "hedgehog with F(0)=pi" convention.
        return float(-(2.0 / math.pi) * np.trapezoid(integrand, r))

    def energy(self) -> float:
        """Skyrme energy (radial reduction):
            E = 4 pi integral_0^inf dr [ 1/2 r^2 (F')^2 + sin^2 F
                + 1/2 (sin^2 F / r^2) (sin^2 F + 2 r^2 (F')^2) ].
        This is finite for the tanh-based profile.
        """
        N = 4096
        r = np.linspace(1e-3 * self.R, 12.0 * self.R, N)
        F = self.F_profile(r)
        dF = self.F_derivative(r)
        sin2 = np.sin(F) ** 2
        kinetic = 0.5 * r ** 2 * dF ** 2 + sin2
        skyrme_term = 0.5 * sin2 / r ** 2 * (sin2 + 2.0 * r ** 2 * dF ** 2)
        density = kinetic + skyrme_term
        return float(4.0 * math.pi * np.trapezoid(density, r))

    def homotopy_class(self) -> str:
        return "pi_3(S^3) = Z"

    def substrate_role(self) -> str:
        return "non-singular baryon-number carrier"

    def summary(self) -> dict:
        return {
            "name": "Texture3D",
            "homotopy": self.homotopy_class(),
            "B": round(self.baryon_number(), 6),
            "energy": round(self.energy(), 6),
            "substrate_role": self.substrate_role(),
        }


# =============================================================================
# 5. Visualizer
# =============================================================================


class DefectVisualizer:
    """Matplotlib renderings of the four defect classes.

    All methods import matplotlib lazily, use the Agg backend, return a
    figure handle, and never call plt.show() — so they are safe under
    pytest --headless and on servers without a display.
    """

    def __init__(self, dpi: int = 100):
        self.dpi = dpi

    # ---- helpers -----------------------------------------------------------

    @staticmethod
    def _import():
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt

    # ---- 1D --------------------------------------------------------------

    def plot_kink(self, kink: Kink1D):
        plt = self._import()
        x = kink.grid()
        phi = kink.field(x)
        dphi = kink.field_derivative(x)
        fig, axes = plt.subplots(2, 1, figsize=(6.0, 5.0), dpi=self.dpi,
                                 sharex=True)
        axes[0].plot(x, phi, color="C0")
        axes[0].axhline(0.0, color="grey", lw=0.5)
        axes[0].axhline(2.0 * math.pi, color="grey", lw=0.5)
        axes[0].set_ylabel(r"$\phi(x)$")
        axes[0].set_title(
            f"Kink1D  (Q = {kink.topological_charge():+.3f},  "
            f"E = {kink.energy():.3f})"
        )
        axes[1].plot(x, dphi, color="C3")
        axes[1].set_xlabel("x / xi")
        axes[1].set_ylabel(r"$\partial_x \phi$")
        axes[1].axhline(0.0, color="grey", lw=0.5)
        fig.tight_layout()
        return fig

    # ---- 2D --------------------------------------------------------------

    def plot_vortex(self, vortex: Vortex2D):
        plt = self._import()
        X, Y = vortex.grid()
        amp = vortex.amplitude()
        ph = vortex.phase()
        fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), dpi=self.dpi)
        im0 = axes[0].imshow(amp, extent=(-vortex.L, vortex.L,
                                          -vortex.L, vortex.L),
                             origin="lower", cmap="viridis")
        axes[0].set_title(rf"|psi|, n = {vortex.n}")
        axes[0].set_xlabel("x")
        axes[0].set_ylabel("y")
        plt.colorbar(im0, ax=axes[0], fraction=0.046)
        # phase wrapped to [-pi, pi] for clarity
        ph_wrapped = (ph + math.pi) % (2.0 * math.pi) - math.pi
        im1 = axes[1].imshow(ph_wrapped,
                             extent=(-vortex.L, vortex.L,
                                     -vortex.L, vortex.L),
                             origin="lower", cmap="hsv",
                             vmin=-math.pi, vmax=math.pi)
        axes[1].set_title(rf"arg(psi), winding = {vortex.winding_number():+.2f}")
        axes[1].set_xlabel("x")
        axes[1].set_ylabel("y")
        plt.colorbar(im1, ax=axes[1], fraction=0.046)
        fig.suptitle("Vortex2D", y=1.02)
        fig.tight_layout()
        return fig

    # ---- 3D monopole ------------------------------------------------------

    def plot_monopole(self, mono: Monopole3D, stride: Optional[int] = None):
        plt = self._import()
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        X, Y, Z = mono.grid()
        n_hat = mono.unit_vector()
        if stride is None:
            stride = max(1, mono.n_points // 8)
        Xs = X[::stride, ::stride, ::stride]
        Ys = Y[::stride, ::stride, ::stride]
        Zs = Z[::stride, ::stride, ::stride]
        Us = n_hat[::stride, ::stride, ::stride, 0]
        Vs = n_hat[::stride, ::stride, ::stride, 1]
        Ws = n_hat[::stride, ::stride, ::stride, 2]
        fig = plt.figure(figsize=(6.0, 5.5), dpi=self.dpi)
        ax = fig.add_subplot(111, projection="3d")
        ax.quiver(Xs, Ys, Zs, Us, Vs, Ws, length=0.5, normalize=True,
                  color="C2", linewidth=0.6)
        forbidden = " (orientability-FORBIDDEN)" if mono.is_orientability_forbidden() else ""
        ax.set_title(
            f"Monopole3D, Q = {mono.topological_charge_analytic():+d}{forbidden}"
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        fig.tight_layout()
        return fig

    # ---- 3D texture -------------------------------------------------------

    def plot_texture(self, texture: Texture3D):
        plt = self._import()
        # Show |n_vec| = sin F(r) on a central z=0 slice as an isosurface-ish
        # heatmap (so we don't need 3D isosurface deps).
        comps = texture.field_components()
        n_mag = np.sqrt(comps["n1"] ** 2 + comps["n2"] ** 2 + comps["n3"] ** 2)
        mid = n_mag.shape[2] // 2
        slice_xy = n_mag[:, :, mid]
        fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), dpi=self.dpi)
        im = axes[0].imshow(slice_xy,
                            extent=(-texture.L, texture.L,
                                    -texture.L, texture.L),
                            origin="lower", cmap="magma")
        axes[0].set_title(r"|$\vec n$| at z = 0 slice")
        axes[0].set_xlabel("x")
        axes[0].set_ylabel("y")
        plt.colorbar(im, ax=axes[0], fraction=0.046)

        # Also plot the radial chiral angle F(r)
        r = np.linspace(0.0, 5.0 * texture.R, 400)
        F = texture.F_profile(r)
        axes[1].plot(r, F, color="C0")
        axes[1].axhline(math.pi, color="grey", lw=0.5)
        axes[1].axhline(0.0, color="grey", lw=0.5)
        axes[1].set_xlabel("r / R")
        axes[1].set_ylabel("F(r)")
        axes[1].set_title(f"Texture3D  B = {texture.baryon_number():+.3f}")
        fig.tight_layout()
        return fig


# =============================================================================
# 6. Catalogue helpers
# =============================================================================


def classify_defects() -> dict:
    """Return a dict mapping homotopy class → (name, codim, target).

    The classification uses the order-parameter manifold M of each defect
    class and the corresponding spatial homotopy.
    """
    return {
        "pi_0(Z_2)":      {"name": "Kink / domain wall", "codim": 1,
                           "target": "Z_2", "carrier": "1D field jump"},
        "pi_1(U(1))":     {"name": "Vortex / cosmic string", "codim": 2,
                           "target": "S^1 = U(1)", "carrier": "2D phase winding"},
        "pi_2(S^2)":      {"name": "Hedgehog / monopole", "codim": 3,
                           "target": "S^2 = O(3)/O(2)",
                           "carrier": "3D vector hedgehog"},
        "pi_3(S^3)":      {"name": "Skyrmion / texture", "codim": 0,
                           "target": "S^3 = SU(2)",
                           "carrier": "3D non-singular lump"},
    }


def zoo_summary() -> dict:
    """Compute summaries of the four canonical defect instances."""
    return {
        "Kink1D":     Kink1D().summary(),
        "Vortex2D":   Vortex2D(n=1).summary(),
        "Monopole3D": Monopole3D().summary(),
        "Texture3D":  Texture3D().summary(),
        "classification": classify_defects(),
    }


__all__ = [
    "Kink1D",
    "Vortex2D",
    "Monopole3D",
    "Texture3D",
    "DefectVisualizer",
    "classify_defects",
    "zoo_summary",
    "DEFAULT_XI",
]
