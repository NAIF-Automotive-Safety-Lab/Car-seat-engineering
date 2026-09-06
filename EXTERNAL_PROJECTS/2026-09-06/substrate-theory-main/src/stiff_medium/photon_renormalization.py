"""Photon wave-function renormalization Z_photon on the Möbius bundle.

This module computes the one-loop radiative correction to the photon propagator
that arises from integrating out Jackiw-Rebbi fermion zero-modes living in the
kink background on the Möbius U(1) bundle (spec §§18.10, 18.34, 18.45, 18.48).

Physical Setup
--------------
The Lagrangian (spec §18.45) contains a Dirac fermion ψ coupled to a sine-Gordon
scalar φ via a Yukawa interaction g_Y φ ψ̄ψ.  Around the kink background:

    φ_kink(x) = (2/ξ) arctan(exp(x/ξ))    [continuum form of tanh kink]
              ≃ M_kink · tanh(x/ξ)          [mass × profile, same kink]

Jackiw and Rebbi (1976) showed that this Yukawa-kink system possesses an
exact fermionic zero-mode — a normalizable solution of the Dirac equation with
eigenvalue zero:

    H_Dirac ψ_0 = 0    with   ψ_0(x) ∝ exp(−∫₀ˣ m(y) dy)
                          = exp(−M_kink ξ ln cosh(x/ξ))
                          = cosh^(−M_kink ξ)(x/ξ)

This is the fermion whose loop corrects the photon propagator.

The 1-loop correction to the photon self-energy Π_μν(q) from the zero-mode
fermion is (schematically, in 1+1D Euclidean):

    Π^(1-loop)_μν(q) = g² ∫ d²k / (2π)²
                          Tr[γ_μ S_F(k) γ_ν S_F(k+q)]

where S_F is the zero-mode Feynman propagator restricted to the zero-mode
sector.  In the zero-mode approximation (valid for energies ≪ M_kink) the
propagator is:

    S_F(k) = i ρ_0(k) / (k² + ε)    with   ρ_0(k) = |ψ̃_0(k)|²

the momentum-space weight of the zero-mode wave function.

The photon wave-function renormalization at q=0 is then:

    Z_photon = 1 + Π(0) / m_photon²

In the massless photon case (m_photon → 0), the relevant quantity is the
residue of the photon pole, extracted as the coefficient of the transverse
projector in Π_μν(q) at small q.  In 1+1D:

    Π_T(q²) = (g²/π²) ∫ d²k ρ_0(k) ρ_0(k+q) / [(k²+μ²)((k+q)²+μ²)]

and Z_photon = 1 − dΠ_T/dq² |_{q²=0}  in the on-shell / Wilsonian sense.

Implementation Strategy
-----------------------
1. Symbolic: express ψ_0(x) = N · cosh^(−Mξ)(x/ξ) using sympy; compute its
   Fourier transform ψ̃_0(k) analytically (it is expressible via the Euler beta
   function / gamma function for rational Mξ).

2. Numerical: define the momentum-space weight ρ_0(k) = |ψ̃_0(k)|² by
   numerically evaluating the Fourier transform on a grid, then use that weight
   in the loop integral.

3. Dimensional regularization analog: the loop integral is UV-finite for the
   zero-mode profile (the cosh profile is in the Schwartz class; its FT decays
   exponentially).  No explicit UV cutoff is needed beyond the kink mass M_kink
   which acts as the natural short-distance regulator.

4. IR: use a small photon mass μ² = m_photon²  (set to M_kink²×10⁻⁶ for numerics)
   to regulate the zero-momentum singularity if needed.

5. Bundle factor: on the Möbius half-flux bundle the coupling of the zero-mode
   to the U(1) gauge field carries an extra topological weight from the holonomy:

       g_bundle = g × c_1    with   c_1 = 1/2   (first Chern class, §18.10)

   This is the geometric origin of the "half-charge" of a kink zero-mode
   relative to a bulk fermion.  The loop integral is multiplied by g_bundle².

6. Scan over β²: using the Coleman bosonization relation (spec §18.9)
       g_thirring = π(4π/β² − 1)
   and identifying the Yukawa coupling g_Y with g_thirring (same Lagrangian),
   compute Z_photon(β²) for the set {3π, 4π, 4.5π, 5π}.

7. Renormalized α: define
       α_ren(β²) = α_bare(β²) × Z_photon(β²)
   and ask which β² gives α_ren = 1/137.

Caveats and Honest Assessment
------------------------------
This is a *model calculation*, not the complete multi-loop field theory.  The
specific approximations are:

- **Zero-mode dominance**: we keep only the Jackiw-Rebbi zero-mode in the
  fermion loop.  Continuum-mode contributions (energies ~ M_kink) are
  suppressed by 1/M_kink² and are genuine higher-order corrections.

- **1+1D reduction**: the zero-mode is localized along the kink direction; the
  transverse (2+1D) directions are integrated assuming free momentum.  This
  gives a finite, well-defined contribution.

- **Möbius bundle topology**: the factor c_1 = 1/2 from the half-flux holonomy
  multiplies the coupling.  This is exact in the Möbius bundle topology (§18.10).

- **No crossed diagrams**: we compute only the bubble (1PI, one fermion loop).
  Box/crossed diagrams enter at two-loop order.

References
----------
    Jackiw & Rebbi (1976) Phys. Rev. D 13, 3398.
    Coleman (1975) Phys. Rev. D 11, 2088.
    spec §18.10, §18.34, §18.45, §18.48.7 item 1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np
from scipy import integrate

# ---------------------------------------------------------------------------
# Import project modules
# ---------------------------------------------------------------------------

from stiff_medium.mobius_bundle import coleman_bosonization_g, chern_class_first
from stiff_medium.alpha_derivation import (
    ALPHA_TARGET,
    INV_ALPHA_TARGET,
    M_KINK_GEV,
    Q_SUBSTRATE_GEV,
    bosonization_alpha,
)
from stiff_medium.rg_running import RGRunning, M_E_GEV

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

PI: Final[float] = math.pi

#: First Chern class of the Möbius bundle (§18.10): c_1 = 1/2.
C1_MOBIUS: Final[float] = chern_class_first()

#: Kink half-width ξ (natural units where M_kink = 1 sets the scale).
#: The profile m(x) = M_kink tanh(x/ξ) has ξ = 1/M_kink in natural units.
XI_NATURAL: Final[float] = 1.0  # ξ = 1 in units of 1/M_kink

#: Product M_kink × ξ = 1 in natural units.  This is the dimensionless coupling
#: that controls the spatial extent of the zero-mode:
#:     ψ_0(x) ∝ cosh^(−M_kink ξ)(x/ξ) = cosh^(−1)(x)
#: In our natural-unit convention we set M_kink = 1 (everything measured in
#: units of the kink mass), so ξ = 1/M_kink = 1 and M_kink × ξ = 1 exactly.
#: The physical value M_KINK_GEV = 27 GeV sets the overall energy scale but
#: does not enter dimensionless loop integrals.
M_KINK_XI: Final[float] = 1.0  # dimensionless, = M_kink × ξ in natural units

#: Small IR regulator mass (in units of M_kink = 1).
MU_IR: Final[float] = 1e-3

#: Integration cutoff in momentum space (in units of M_kink = 1).
#: The zero-mode FT decays exponentially for |k| ≫ M_kink, so 20 M_kink
#: captures > 99.99% of the spectral weight.
K_CUTOFF: Final[float] = 20.0

#: Number of points for momentum-space quadrature (must be odd for Simpson).
N_K_POINTS: Final[int] = 801


# ---------------------------------------------------------------------------
# §1: Jackiw-Rebbi zero-mode wave function
# ---------------------------------------------------------------------------


def zero_mode_position(x: float, m_kink_xi: float = M_KINK_XI) -> float:
    """Jackiw-Rebbi zero-mode wave function in position space.

    The normalized zero-mode of a Dirac fermion in the kink background
    m(x) = M_kink tanh(x/ξ) is (Jackiw & Rebbi 1976):

        ψ_0(x) = N · cosh^(−Mξ)(x/ξ)

    where the normalization N is determined by ∫|ψ_0|² dx = 1.

    In natural units (M_kink = ξ = 1, so Mξ = 1):

        ψ_0(x) = N · cosh^(−1)(x) = N · sech(x)

    For general Mξ: ψ_0(x) ∝ cosh^(−Mξ)(x).

    The normalization integral is:

        ∫_{-∞}^{∞} cosh^(−2Mξ)(x) dx = √π Γ(Mξ) / Γ(Mξ + 1/2)
                                        = √π / (Mξ · B(Mξ, 1/2))^(−1)

    For Mξ = 1: ∫ sech²(x) dx = 2, so N = 1/√2.

    Args:
        x: Position coordinate in units of ξ.
        m_kink_xi: Dimensionless product M_kink × ξ.  Controls the
            "width" of the zero-mode (larger Mξ → sharper profile).

    Returns:
        Unnormalized wave function value ψ_0(x).  Caller is responsible
        for normalization.
    """
    # Avoid overflow for large |x|: cosh(x) ~ exp(|x|)/2 → cosh^(-Mξ) → 0
    arg = abs(x)
    if arg > 500.0 / m_kink_xi:
        return 0.0
    return math.cosh(x) ** (-m_kink_xi)


def normalization_constant(m_kink_xi: float = M_KINK_XI) -> float:
    """Normalization constant N for the Jackiw-Rebbi zero-mode.

    Computes N such that ∫_{-∞}^{∞} |N · cosh^(−Mξ)(x)|² dx = 1.

    The integral ∫ cosh^(−2Mξ)(x) dx = √π Γ(Mξ) / Γ(Mξ + 1/2)
    using the standard result ∫_{-∞}^{∞} sech^(2ν)(x) dx = √π Γ(ν)/Γ(ν+1/2).

    Args:
        m_kink_xi: Dimensionless product M_kink × ξ.

    Returns:
        Normalization constant N = 1/√(∫|ψ_unnorm|² dx).
    """
    from scipy.special import gamma

    # ∫_{-∞}^{∞} cosh^{-2 m_kink_xi}(x) dx  = √π Γ(m_kink_xi) / Γ(m_kink_xi + 0.5)
    norm_sq_inv = math.sqrt(PI) * gamma(m_kink_xi) / gamma(m_kink_xi + 0.5)
    return 1.0 / math.sqrt(norm_sq_inv)


# ---------------------------------------------------------------------------
# §2: Momentum-space zero-mode weight ρ_0(k)
# ---------------------------------------------------------------------------


def build_zero_mode_ft(
    k_array: np.ndarray,
    m_kink_xi: float = M_KINK_XI,
    x_cutoff: float = 60.0,
    n_x: int = 4001,
) -> np.ndarray:
    """Compute the Fourier transform ψ̃_0(k) of the Jackiw-Rebbi zero-mode.

    The Fourier transform is computed numerically via direct integration:

        ψ̃_0(k) = ∫_{-∞}^{∞} ψ_0(x) e^{−ikx} dx

    For Mξ = 1: ψ_0(x) = N·sech(x) and its FT is known analytically:

        ψ̃_0(k) = N · π · sech(πk/2)

    For general Mξ: the FT is expressible in terms of the digamma function
    but is easiest computed numerically for the loop integral purpose.

    We use the half-symmetry (ψ_0 is even in x) so:

        ψ̃_0(k) = 2 ∫_0^∞ ψ_0(x) cos(kx) dx    (real and even in k)

    Args:
        k_array: 1-D array of momentum values (in units of M_kink = 1).
        m_kink_xi: M_kink × ξ (default 1.0 in natural units).
        x_cutoff: Upper limit for position-space integration.
        n_x: Number of quadrature points in position space.

    Returns:
        1-D array of |ψ̃_0(k)|² (the spectral weight, non-negative).
        Shape matches k_array.
    """
    N_norm = normalization_constant(m_kink_xi)

    # Position-space grid (half-line, use symmetry)
    x_grid = np.linspace(0.0, x_cutoff, n_x)
    dx = x_grid[1] - x_grid[0]

    # Unnormalized profile on the grid
    psi_x = np.array([math.cosh(x) ** (-m_kink_xi) for x in x_grid])
    psi_x *= N_norm

    # For each k, compute ψ̃_0(k) = 2 ∫_0^∞ ψ_0(x) cos(kx) dx  (even FT)
    ft_sq = np.zeros(len(k_array), dtype=float)
    for i, k in enumerate(k_array):
        integrand = psi_x * np.cos(k * x_grid)
        ft_k = 2.0 * np.trapezoid(integrand, x_grid)  # factor 2 for the full line
        ft_sq[i] = ft_k**2

    return ft_sq


# ---------------------------------------------------------------------------
# §3: One-loop photon self-energy — position-space overlap integral
# ---------------------------------------------------------------------------
#
# Physical derivation:
# --------------------
# The photon self-energy from the zero-mode current-current correlator is, in
# 1+1D Euclidean position space (per the definition in the module docstring):
#
#   Π_00 ∝ ∫d²x d²y  J_0(x) D(x-y) J_0(y)
#
# where:
#   J_0(x) = ψ_0^†(x) ψ_0(x) = |ψ_0(x)|²    [charge density of the zero-mode]
#   D(x-y) = −(1/2π) ln(|x-y| μ_IR)          [1+1D massless scalar propagator]
#
# The zero-mode ψ_0(x) is localized at the kink; J_0 = |ψ_0|² is therefore a
# smooth, rapidly decaying function.  Both position-space integrals are
# dominated by |x|, |y| ≲ ξ = 1/M_kink.
#
# In our 1+1D (kink-axis) reduction, the transverse coordinates are trivially
# integrated: d²x → dx (the zero-mode is uniform in the transverse direction
# with a transverse normalization already absorbed into N).
#
# The integral becomes:
#
#   I_bundle = ∫_{-L}^{L} dx ∫_{-L}^{L} dy  ρ(x) D(x-y) ρ(y)
#
# where ρ(x) = |ψ_0(x)|² = N² cosh^{-2Mξ}(x) is the normalized zero-mode density.
#
# The propagator has a logarithmic divergence at x = y, regulated by the kink
# width ξ (a physical short-distance cutoff).  We use:
#
#   D(r) = −(1/2π) ln(max(|r|, ξ_reg) × μ_IR)
#
# where ξ_reg = 0.01 (in units of ξ = 1) regularizes the contact term.
#
# The photon wave-function renormalization is then:
#
#   Z_photon = 1 + (g_yukawa²/π²) × I_bundle
#
# The sign: for the U(1) gauge theory, the vacuum polarization gives a POSITIVE
# correction to the photon kinetic term (screening), so Z_photon > 1 and
# α_ren = α_bare × Z_photon > α_bare (standard screening).
# Note: whether this increases or decreases 1/α depends on the convention.
# We use: α_ren = α_bare / Z_photon  (renormalization REDUCES the coupling
# as seen at low energies — the physical α is LESS than the bare coupling).
# This is the standard QED convention where Z_photon > 1 implies screening.


def _bundle_propagator_1d(r: float, xi_reg: float = 0.01, mu_ir: float = MU_IR) -> float:
    """1+1D Euclidean massless scalar propagator, regulated.

    D(r) = −(1/2π) ln(max(|r|, ξ_reg) × μ_IR)

    The positive sign convention: D(r) > 0 for r ≪ 1/μ_IR (i.e. for separations
    smaller than the IR scale), giving a positive contribution to Z_photon.

    We flip the overall sign so that D(r) = +(1/2π) ln(1/(max(|r|,ξ_reg) μ_IR))
    which is positive for μ_IR × |r| < 1, i.e. at sub-IR separations.

    Args:
        r: Separation |x-y| in units of ξ = 1/M_kink.
        xi_reg: Short-distance regulator (contact-term cutoff) in units of ξ.
        mu_ir: IR mass regulator in units of M_kink (= μ_IR × ξ).

    Returns:
        D(r) [dimensionless].
    """
    sep = max(abs(r), xi_reg)
    # Standard 1+1D propagator: D = -(1/2π) ln(sep × mu_ir)
    # For sep × mu_ir < 1 (typical short-distance physics), ln < 0, so -ln > 0.
    return -(1.0 / (2.0 * PI)) * math.log(sep * mu_ir)


def compute_i_bundle(
    m_kink_xi: float = M_KINK_XI,
    x_cutoff: float = 30.0,
    n_x: int = 201,
    xi_reg: float = 0.01,
    mu_ir: float = MU_IR,
) -> float:
    """Compute the position-space overlap integral I_bundle.

    I_bundle = ∫dx ∫dy ρ(x) D(x-y) ρ(y)

    where ρ(x) = |ψ_0(x)|² is the zero-mode density and D is the 1+1D
    bundle propagator.

    The double integral is evaluated via 2D Simpson quadrature on a grid.

    Args:
        m_kink_xi: Dimensionless M_kink × ξ (= 1 in natural units).
        x_cutoff: Half-width of integration domain [in units of ξ].
        n_x: Number of grid points per axis (odd for Simpson's rule).
        xi_reg: Contact-term regulator for D(r) at r→0.
        mu_ir: IR mass regulator [M_kink units].

    Returns:
        I_bundle [dimensionless].
    """
    N_norm = normalization_constant(m_kink_xi)

    # Position grid
    x_grid = np.linspace(-x_cutoff, x_cutoff, n_x)
    dx = x_grid[1] - x_grid[0]

    # Zero-mode density ρ(x) = N² cosh^{-2Mξ}(x)
    rho = np.array([N_norm**2 * math.cosh(x) ** (-2.0 * m_kink_xi) for x in x_grid])

    # Double integral: I = ∫∫ ρ(x) D(x-y) ρ(y) dx dy
    # Exploit symmetry D(x-y) = D(y-x) and ρ(x) = ρ(-x) to compute upper triangle.
    # Use separability: precompute the inner integral
    #   f(x) = ∫ D(x-y) ρ(y) dy
    # then I = ∫ ρ(x) f(x) dx

    f = np.zeros(n_x)
    for i, xi in enumerate(x_grid):
        # Inner integral: f(xi) = ∫ D(xi - y) ρ(y) dy
        inner = np.array([
            _bundle_propagator_1d(xi - x_grid[j], xi_reg=xi_reg, mu_ir=mu_ir) * rho[j]
            for j in range(n_x)
        ])
        f[i] = np.trapezoid(inner, x_grid)

    # Outer integral: I = ∫ ρ(x) f(x) dx
    i_bundle = float(np.trapezoid(rho * f, x_grid))
    return i_bundle


# ---------------------------------------------------------------------------
# §4: (Momentum-space FT interpolator — kept for verification / other uses)
# ---------------------------------------------------------------------------


def build_rho_interpolator(
    m_kink_xi: float = M_KINK_XI,
    k_cutoff: float = K_CUTOFF,
    n_k: int = 401,
) -> tuple[np.ndarray, np.ndarray]:
    """Pre-compute ρ_0(k) = |ψ̃_0(k)|² on a momentum grid for verification.

    Args:
        m_kink_xi: M_kink × ξ in natural units.
        k_cutoff: Maximum |k| to evaluate [M_kink units].
        n_k: Number of momentum points.

    Returns:
        Tuple (k_grid, rho_grid) where both are 1-D arrays.
    """
    k_grid = np.linspace(-k_cutoff, k_cutoff, n_k)
    k_pos_grid = np.linspace(0.0, k_cutoff, (n_k + 1) // 2)
    rho_pos = build_zero_mode_ft(k_pos_grid, m_kink_xi=m_kink_xi)
    rho_neg = rho_pos[::-1][:-1]
    rho_full = np.concatenate([rho_neg, rho_pos])
    assert len(rho_full) == n_k, f"Length mismatch: {len(rho_full)} vs {n_k}"
    return k_grid, rho_full


# ---------------------------------------------------------------------------
# §5: Main Z_photon computation (position-space overlap integral)
# ---------------------------------------------------------------------------


@dataclass
class ZPhotonResult:
    """Result of Z_photon computation for a given (β², M_kink) pair.

    Attributes:
        beta_squared: Sine-Gordon coupling β².
        m_kink_xi: Dimensionless product M_kink × ξ (= 1 in natural units).
        g_thirring: Coleman-bosonized Thirring coupling g = π(4π/β² − 1).
        g_yukawa: Yukawa coupling used in the loop = g_thirring × c_1 (Möbius factor).
        alpha_bare: Bare fine-structure constant from bosonization (convention 1).
        i_bundle: Position-space overlap integral I_bundle = ∫∫ρ(x)D(x-y)ρ(y)dxdy.
        z_photon: Photon wave-function renormalization:
            Z_photon = 1 + (g_yukawa²/π²) × I_bundle.
            Z_photon > 1 implies screening (α_ren < α_bare).
        alpha_renormalized: Renormalized α = α_bare / Z_photon.
        inv_alpha_renormalized: 1/α_ren.
        delta_from_target: |1/α_ren − 137.036|.
        satisfies_alpha: True if delta_from_target < 5 (within ~4% of 1/137).
        note: Human-readable interpretation.
    """

    beta_squared: float
    m_kink_xi: float
    g_thirring: float
    g_yukawa: float
    alpha_bare: float
    i_bundle: float
    z_photon: float
    alpha_renormalized: float
    inv_alpha_renormalized: float
    delta_from_target: float
    satisfies_alpha: bool
    note: str


def compute_z_photon(
    beta_squared: float,
    m_kink_xi: float = M_KINK_XI,
    mu_ir: float = MU_IR,
    x_cutoff: float = 30.0,
    n_x: int = 151,
    xi_reg: float = 0.01,
) -> ZPhotonResult:
    """Compute Z_photon(β², M_kink) via the position-space overlap integral.

    Algorithm:
    1. Get g_thirring from Coleman bosonization.
    2. Apply Möbius bundle factor: g_yukawa = g_thirring × c_1 = g_thirring / 2.
    3. Compute I_bundle = ∫∫ρ(x) D(x-y) ρ(y) dx dy  (position-space double integral)
       where ρ(x) = |ψ_0(x)|² (Jackiw-Rebbi zero-mode density)
       and D(r) = −(1/2π) ln(|r| μ_IR)  (1+1D massless bundle propagator).
    4. Z_photon = 1 + (g_yukawa²/π²) × I_bundle.
    5. α_ren = α_bare / Z_photon.
       (Z > 1 screens: high-energy bare coupling > low-energy physical coupling.)

    Physical interpretation of Z_photon:
    - Z_photon > 1: loop correction screens the photon — the zero-mode fermion
      cloud around the kink reduces the effective charge at long distances.
      This is analogous to QED vacuum polarization, and α_ren < α_bare.
    - Z_photon < 1: would imply anti-screening (would require I_bundle < 0,
      which cannot happen for the positive-definite propagator at sub-IR scales).

    Args:
        beta_squared: Sine-Gordon coupling β².
        m_kink_xi: M_kink × ξ in natural units (default 1.0).
        mu_ir: IR regulator mass [M_kink units].  Controls the IR log cutoff.
        x_cutoff: Half-width of position-space integration [ξ units].
        n_x: Number of grid points per axis in the double integral.
        xi_reg: Contact-term UV regulator for D(r) at r→0.

    Returns:
        ZPhotonResult with all intermediate quantities.

    Raises:
        ValueError: If beta_squared is not in the physical range (0, 8π).
    """
    if not (0.0 < beta_squared < 8.0 * PI):
        raise ValueError(
            f"beta_squared={beta_squared:.4f} outside physical range (0, 8π={8*PI:.4f}). "
            "β² must be in the soliton phase (0 < β² < 8π)."
        )

    # --- Step 1: Coleman bosonization ---
    bos = bosonization_alpha(beta_squared)
    g_thirring = bos.g_thirring
    alpha_bare = bos.alpha_bare

    if g_thirring <= 0.0:
        # Free or attractive Thirring regime: no positive coupling, no loop.
        # For β² ≥ 4π the Coleman duality gives g ≤ 0 (free / attractive fermion).
        # The zero-mode is still present (it exists for any kink background) but
        # the coupling to the gauge field vanishes (g_yukawa = 0).
        return ZPhotonResult(
            beta_squared=beta_squared,
            m_kink_xi=m_kink_xi,
            g_thirring=g_thirring,
            g_yukawa=0.0,
            alpha_bare=0.0,
            i_bundle=0.0,
            z_photon=1.0,
            alpha_renormalized=0.0,
            inv_alpha_renormalized=float("inf"),
            delta_from_target=float("inf"),
            satisfies_alpha=False,
            note=(
                f"β²={beta_squared/PI:.4f}π: g_Thirring={g_thirring:.6f} ≤ 0 "
                "(free or attractive Thirring regime). No positive α_bare from "
                "bosonization. Physical EM domain requires β² < 4π."
            ),
        )

    # --- Step 2: Möbius bundle coupling factor ---
    # The zero-mode couples to A_μ with strength g_bundle = g_thirring × c_1
    # where c_1 = 1/2 is the first Chern class (spec §18.10).
    g_yukawa = g_thirring * C1_MOBIUS  # = g_thirring / 2

    # --- Step 3: Position-space overlap integral ---
    i_bundle = compute_i_bundle(
        m_kink_xi=m_kink_xi,
        x_cutoff=x_cutoff,
        n_x=n_x,
        xi_reg=xi_reg,
        mu_ir=mu_ir,
    )

    # --- Step 4: Z_photon ---
    # From the problem statement and 1+1D Euclidean field theory:
    #
    #   Z_photon = 1 + (g²/π²) × I_bundle
    #
    # where g = g_yukawa (includes Möbius factor) and I_bundle is the
    # position-space current-current correlator integrated with the bundle
    # propagator.
    #
    # Sign: I_bundle is positive (D(r) > 0 at short separations, ρ > 0).
    # So Z_photon > 1, which means the photon kinetic term is enhanced by the
    # loop, and the physical (pole) coupling is REDUCED:
    #   α_ren = α_bare / Z_photon < α_bare
    # This is standard QED-style charge screening.

    z_photon = 1.0 + (g_yukawa**2 / PI**2) * i_bundle

    # Protect against non-physical values
    if z_photon <= 0.0:
        z_photon = 1e-8

    # --- Step 5: Renormalized α ---
    if alpha_bare > 0.0:
        alpha_ren = alpha_bare / z_photon
        inv_alpha_ren = 1.0 / alpha_ren
        delta = abs(inv_alpha_ren - INV_ALPHA_TARGET)
        satisfies = delta < 5.0
    else:
        alpha_ren = 0.0
        inv_alpha_ren = float("inf")
        delta = float("inf")
        satisfies = False

    # --- Compose note ---
    if alpha_bare > 0.0:
        note = (
            f"β²={beta_squared/PI:.4f}π | "
            f"g_Thirring={g_thirring:.6f} | "
            f"g_Yukawa(×c₁)={g_yukawa:.6f} | "
            f"α_bare=1/{1/alpha_bare:.2f} | "
            f"I_bundle={i_bundle:.6f} | "
            f"Z_γ={z_photon:.6f} | "
            f"α_ren=1/{inv_alpha_ren:.3f} | "
            f"Δ(1/α)={delta:.4f}"
        )
    else:
        note = f"β²={beta_squared/PI:.4f}π | g_Thirring={g_thirring:.6f} | g≤0, no physical α"

    return ZPhotonResult(
        beta_squared=beta_squared,
        m_kink_xi=m_kink_xi,
        g_thirring=g_thirring,
        g_yukawa=g_yukawa,
        alpha_bare=alpha_bare,
        i_bundle=i_bundle,
        z_photon=z_photon,
        alpha_renormalized=alpha_ren,
        inv_alpha_renormalized=inv_alpha_ren,
        delta_from_target=delta,
        satisfies_alpha=satisfies,
        note=note,
    )


# ---------------------------------------------------------------------------
# §6: Scan over β² values
# ---------------------------------------------------------------------------


def scan_beta_squared(
    beta_sq_values: list[float] | None = None,
    m_kink_xi: float = M_KINK_XI,
    verbose: bool = False,
) -> list[ZPhotonResult]:
    """Compute Z_photon for a list of β² values.

    Default scan: {3π, 4π, 4.5π, 5π} as specified in the task.
    Additional values around the Higgs/W constraint β² ≈ 4.54π are included.

    Args:
        beta_sq_values: List of β² values to scan.  If None, uses the default set.
        m_kink_xi: M_kink × ξ (natural units).
        verbose: If True, print progress to stdout.

    Returns:
        List of ZPhotonResult, one per β² value.
    """
    if beta_sq_values is None:
        beta_sq_values = [
            3.0 * PI,    # below free-fermion point
            4.0 * PI,    # free-fermion point (Coleman duality point)
            4.5 * PI,    # between free-fermion and Higgs/W constraint
            4.54 * PI,   # Higgs/W constraint β² (spec §18.52.4)
            5.0 * PI,    # above Higgs/W constraint, near KT
        ]

    results: list[ZPhotonResult] = []
    for beta_sq in beta_sq_values:
        if verbose:
            print(f"  Computing Z_photon at β²={beta_sq/PI:.4f}π ...")
        try:
            result = compute_z_photon(beta_sq, m_kink_xi=m_kink_xi, n_x=201)
        except ValueError as exc:
            if verbose:
                print(f"    Skipped: {exc}")
            continue
        results.append(result)
        if verbose:
            print(f"    {result.note}")

    return results


# ---------------------------------------------------------------------------
# §7: Find the β² that best satisfies α_ren = 1/137
# ---------------------------------------------------------------------------


def find_best_z_photon_beta(
    beta_sq_min: float = 0.5 * PI,
    beta_sq_max: float = 4.0 * PI - 0.01,
    n_scan: int = 500,
    m_kink_xi: float = M_KINK_XI,
    n_x: int = 101,
) -> dict[str, object]:
    """Scan β² finely and find the value that minimises |1/α_ren − 137.036|.

    Only scans β² < 4π (repulsive Thirring / positive g regime), since
    the bosonization only gives positive α_bare in that range.

    Note: I_bundle does NOT depend on β² (it depends only on m_kink_xi and the
    zero-mode profile, which is fixed). Only g_yukawa(β²) changes.  We therefore
    compute I_bundle once and reuse it for all scan points.

    Args:
        beta_sq_min: Lower bound of β² scan.
        beta_sq_max: Upper bound (just below 4π).
        n_scan: Number of scan points.
        m_kink_xi: M_kink × ξ.
        n_x: Grid size for position-space integral (coarser = faster for scan).

    Returns:
        Dictionary with 'best_beta_sq', 'best_result', 'all_results'.
    """
    # I_bundle is independent of β²: compute it once.
    i_bundle_cached = compute_i_bundle(
        m_kink_xi=m_kink_xi, x_cutoff=30.0, n_x=n_x, xi_reg=0.01, mu_ir=MU_IR
    )

    beta_sq_arr = np.linspace(beta_sq_min, beta_sq_max, n_scan)
    all_results: list[ZPhotonResult] = []

    for beta_sq in beta_sq_arr:
        try:
            bos = bosonization_alpha(float(beta_sq))
        except ValueError:
            continue
        g_t = bos.g_thirring
        if g_t <= 0.0:
            continue
        a_bare = bos.alpha_bare
        if a_bare <= 0.0:
            continue

        g_yuk = g_t * C1_MOBIUS
        z = 1.0 + (g_yuk**2 / PI**2) * i_bundle_cached
        if z <= 0.0:
            z = 1e-8
        a_ren = a_bare / z
        inv_a_ren = 1.0 / a_ren
        delta = abs(inv_a_ren - INV_ALPHA_TARGET)

        all_results.append(ZPhotonResult(
            beta_squared=float(beta_sq),
            m_kink_xi=m_kink_xi,
            g_thirring=g_t,
            g_yukawa=g_yuk,
            alpha_bare=a_bare,
            i_bundle=i_bundle_cached,
            z_photon=z,
            alpha_renormalized=a_ren,
            inv_alpha_renormalized=inv_a_ren,
            delta_from_target=delta,
            satisfies_alpha=(delta < 5.0),
            note=f"β²={float(beta_sq)/PI:.4f}π | scan point",
        ))

    if not all_results:
        return {"best_beta_sq": None, "best_result": None, "all_results": []}

    best = min(all_results, key=lambda r: r.delta_from_target)

    return {
        "best_beta_sq": best.beta_squared,
        "best_beta_sq_in_pi": best.beta_squared / PI,
        "best_result": best,
        "all_results": all_results,
        "n_physical": len(all_results),
    }


# ---------------------------------------------------------------------------
# §8: Symbolic Jackiw-Rebbi profile (sympy, for documentation + verification)
# ---------------------------------------------------------------------------


def symbolic_zero_mode():
    """Return the symbolic Jackiw-Rebbi zero-mode using sympy.

    This documents the analytical form used in the numerical computation.
    The symbolic expression is:

        ψ_0(x) = N · cosh^(−Mξ)(x/ξ)

    with the kink profile m(x) = M tanh(x/ξ).

    Returns:
        Dictionary with sympy expressions for profile, zero-mode, and normalization.
    """
    try:
        import sympy as sp

        x, M, xi, k = sp.symbols("x M xi k", real=True)
        m_kink_xi_sym = M * xi  # dimensionless product

        # Kink profile
        kink_profile = M * sp.tanh(x / xi)

        # Jackiw-Rebbi zero-mode (unnormalized)
        psi_0 = sp.cosh(x / xi) ** (-m_kink_xi_sym)

        # Normalization integral (symbolic, for general Mξ)
        # ∫ cosh^{-2Mξ}(x/ξ) dx  = ξ √π Γ(Mξ) / Γ(Mξ + 1/2)
        norm_integral = xi * sp.sqrt(sp.pi) * sp.gamma(m_kink_xi_sym) / sp.gamma(
            m_kink_xi_sym + sp.Rational(1, 2)
        )

        # Momentum-space FT (analytic, for Mξ = 1 special case)
        # ψ̃_0(k) = N · π · sech(π k ξ / 2)  [only for Mξ=1 exactly]
        psi_0_ft_mxi1 = sp.pi * sp.sech(sp.pi * k * xi / 2)

        return {
            "kink_profile": kink_profile,
            "zero_mode_unnorm": psi_0,
            "norm_integral": norm_integral,
            "ft_mxi_equals_1": psi_0_ft_mxi1,
            "symbols": {"x": x, "M": M, "xi": xi, "k": k},
        }
    except ImportError:
        return {"error": "sympy not available; numerical path used exclusively"}
