"""Two-loop corrections to Z_photon on the Möbius bundle (spec §§18.34, 18.45, 18.48).

This module computes candidate two-loop corrections to the photon
wave-function renormalization.  The current reproducible audit result is
negative: the corrections are too small to turn the present beta map into a
first-principles derivation of alpha.

  1. The two-loop "sunset" self-energy:
       ψ̄ → ψ_0 → A_μ → ψ → ψ_0   (two nested zero-mode propagators + one boson line)

  2. The "figure-8" diagram:
       Two disconnected zero-mode bubbles (product of two 1-loop integrals)

  3. Continuum (scattering-state) corrections:
       Integration over scattering states ψ_k for momenta k > M_kink (set to 1 in
       natural units)

Scheme for Z_photon at two-loop order:

    Z_photon = 1 + (g²/π²) I_1 + (g²/π²)² I_2 + I_continuum

where:
  I_1      = 1-loop bubble integral (from photon_renormalization.py)
  I_2      = 2-loop combination (sunset + figure-8)
  I_continuum = contribution from scattering states

The two-loop expansion parameter is λ = g²/π².  At β² ≈ 3.91π,
g_Yukawa ≈ 0.036, so:

    (g_Yukawa²/π²)  ≈  1.36e-4

This makes the genuine sunset and figure-8 two-loop pieces tiny.  The
continuum estimate is larger than those terms but still changes Z_photon only
at O(1e-6).

HONEST PREAMBLE
---------------
With the fixed grid point β² = 3.91π, the current script reports
1/α_2loop = 136.507203650 against the CODATA 2022 target
137.035999177, low by 0.528795527 (0.385880739%).  A continuous β² solve can
match the target exactly, but that is a calibration of a free parameter, not a
prediction.  The independently motivated Higgs/W beta is β² = 4.547836π,
where the naive Coleman map has g_Thirring <= 0 and gives no positive
alpha_bare.

References
----------
    Jackiw & Rebbi (1976) Phys. Rev. D 13, 3398 — zero mode.
    Coleman (1975) Phys. Rev. D 11, 2088 — bosonization duality.
    spec §18.34, §18.45, §18.48.
    photon_renormalization.py — the 1-loop baseline computation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

from stiff_medium.photon_renormalization import (
    C1_MOBIUS,
    M_KINK_XI,
    MU_IR,
    PI,
    compute_i_bundle,
    find_best_z_photon_beta,
    normalization_constant,
)
from stiff_medium.alpha_derivation import (
    ALPHA_TARGET,
    INV_ALPHA_TARGET,
    bosonization_alpha,
    higgs_w_constraint,
)

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

#: UV cutoff in units of M_kink = 1 (as specified in the task).
LAMBDA_CUTOFF: Final[float] = 100.0

#: IR regulator mass in units of M_kink.
MU_IR_2L: Final[float] = MU_IR  # = 1e-3

#: Contact-term regulator for the position-space propagator.
XI_REG: Final[float] = 0.01

#: Position-space integration half-width.
X_CUTOFF: Final[float] = 15.0

#: Number of grid points for position-space integrals (odd for Simpson).
N_X_COARSE: Final[int] = 51   # coarse: for 3D nested integrals
N_X_FINE:   Final[int] = 101  # fine: for verification


# ---------------------------------------------------------------------------
# Helper: position-space propagator (same as in photon_renormalization.py)
# ---------------------------------------------------------------------------


def _prop(r: float, xi_reg: float = XI_REG, mu_ir: float = MU_IR_2L) -> float:
    """1+1D Euclidean massless scalar propagator, regulated.

    D(r) = −(1/(2π)) ln(max(|r|, ξ_reg) × μ_IR)
    Positive for sub-IR separations (μ_IR × |r| < 1).

    Args:
        r: Separation |x - y| in units of ξ = 1/M_kink.
        xi_reg: Short-distance regulator (contact cutoff).
        mu_ir: IR mass regulator.

    Returns:
        D(r) [dimensionless].
    """
    sep = max(abs(r), xi_reg)
    return -(1.0 / (2.0 * PI)) * math.log(sep * mu_ir)


def _zero_mode_density(x: float, m_kink_xi: float, N_norm: float) -> float:
    """Normalized zero-mode density ρ(x) = N² cosh^{-2Mξ}(x).

    Args:
        x: Position in units of ξ.
        m_kink_xi: M_kink × ξ (= 1 in natural units).
        N_norm: Normalization constant from normalization_constant().

    Returns:
        ρ(x) = |ψ_0(x)|².
    """
    arg = abs(x)
    if arg > 500.0 / m_kink_xi:
        return 0.0
    return N_norm**2 * math.cosh(x) ** (-2.0 * m_kink_xi)


# ---------------------------------------------------------------------------
# Two-loop sunset diagram
# ---------------------------------------------------------------------------


def compute_i_sunset(
    m_kink_xi: float = M_KINK_XI,
    x_cutoff: float = X_CUTOFF,
    n_x: int = N_X_COARSE,
    xi_reg: float = XI_REG,
    mu_ir: float = MU_IR_2L,
) -> float:
    """Compute the two-loop sunset diagram I_sunset.

    I_sunset = ∫dx ∫dy ∫dz  ρ(x) D(x-y) ρ(y) D(y-z) ρ(z) D(z-x)

    Evaluated on a 3D grid using repeated trapezoidal quadrature.

    Args:
        m_kink_xi: Dimensionless M_kink × ξ (= 1 in natural units).
        x_cutoff: Half-width of integration domain [units of ξ].
        n_x: Number of grid points per axis.
        xi_reg: Contact-term regulator for D(r) at r→0.
        mu_ir: IR mass regulator.

    Returns:
        I_sunset [dimensionless].
    """
    N_norm = normalization_constant(m_kink_xi)
    grid = np.linspace(-x_cutoff, x_cutoff, n_x)

    # Pre-compute rho on the grid; the triple integral is O(n^3).
    rho = np.array([_zero_mode_density(x, m_kink_xi, N_norm) for x in grid])

    d_matrix = np.zeros((n_x, n_x))
    for i in range(n_x):
        for j in range(n_x):
            d_matrix[i, j] = _prop(grid[i] - grid[j], xi_reg=xi_reg, mu_ir=mu_ir)

    outer = np.zeros(n_x)
    for i in range(n_x):
        inner_k = np.zeros(n_x)
        for j in range(n_x):
            f_z = d_matrix[j, :] * rho * d_matrix[:, i]
            inner_k[j] = np.trapezoid(f_z, grid)

        f_y = d_matrix[i, :] * rho * inner_k
        outer[i] = np.trapezoid(f_y, grid)

    i_sunset = float(np.trapezoid(rho * outer, grid))
    return i_sunset


# ---------------------------------------------------------------------------
# Figure-8 diagram: square of the one-loop bubble
# ---------------------------------------------------------------------------


def compute_i_figure8(i_bundle_1loop: float) -> float:
    """Compute the figure-8 two-loop correction.

    In the perturbative expansion of Z_photon:

        Z = exp[(g²/π²) I_1 + higher] ≈ 1 + (g²/π²) I_1 + (1/2)(g²/π²)² I_1² + ...

    The figure-8 term at order (g²/π²)² is (1/2) I_1².
    This is the leading "self-energy product" diagram.

    Args:
        i_bundle_1loop: The 1-loop integral I_1 (from compute_i_bundle).

    Returns:
        I_figure8 = (1/2) × I_1²  [the combination entering ΔZ_photon at 2-loop].
        Note: the (g²/π²)² prefactor is applied in the caller.
    """
    return 0.5 * i_bundle_1loop**2


# ---------------------------------------------------------------------------
# Continuum scattering-state correction
# ---------------------------------------------------------------------------


def compute_i_continuum(
    k_min: float = 1.0,
    lambda_cutoff: float = LAMBDA_CUTOFF,
    n_k: int = 200,
) -> float:
    """Compute the continuum scattering-state correction to Z_photon.

    The zero-mode/scattering cross term is zero by orthogonality; this estimates
    the pure continuum-continuum bubble with the Krein-Friedel-Lloyd density
    shift for the reflectionless Poschl-Teller potential.
    """
    k_grid = np.linspace(k_min, lambda_cutoff, n_k)
    int_vals = 1.0 / (PI * (k_grid**2 + 1.0) ** 2)
    return float(np.trapezoid(int_vals, k_grid))


# ---------------------------------------------------------------------------
# §4: Two-loop Z_photon and renormalized α
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TwoLoopResult:
    """Compact result of the current two-loop Z_photon map."""

    beta_squared: float
    alpha_bare: float
    z_photon_2loop: float
    alpha_2loop: float
    inv_alpha_2loop: float


@dataclass(frozen=True)
class AlphaTwoLoopAudit:
    """Compact conclusion for the current two-loop alpha computation."""

    target_inv_alpha: float
    target_alpha: float
    baseline: TwoLoopResult
    fixed_beta: TwoLoopResult
    fitted_beta: TwoLoopResult | None
    higgs_w_beta_pi: float
    higgs_w_ratio: float
    higgs_w_g_thirring: float
    higgs_w_has_positive_alpha: bool
    conclusion: str


def compute_z_photon_2loop(
    beta_squared: float,
) -> TwoLoopResult:
    """Compute Z_photon through two-loop order.

    Combines:
      - 1-loop bubble (from compute_i_bundle via position-space overlap)
      - 2-loop sunset (triple position-space integral)
      - 2-loop figure-8 (square of 1-loop integral)
      - Continuum correction (Krein-Friedel-Lloyd density of states)

    The total Z_photon is:

        Z_2loop = 1
                 + (g²/π²) I_1              [1-loop bubble]
                 + (g²/π²)² I_sunset        [2-loop sunset: 1PI]
                 + (g²/π²)² (1/2 I_1²)     [2-loop figure-8: resumming the bubble]
                 + (g²/π²) I_continuum      [1-loop continuum]

    Args:
        beta_squared: Sine-Gordon coupling β².

    Returns:
        TwoLoopResult with the compact audit values.

    Raises:
        ValueError: If beta_squared outside (0, 8π).
    """
    if not (0.0 < beta_squared < 8.0 * PI):
        raise ValueError(
            f"beta_squared={beta_squared:.4f} outside physical range (0, 8π={8*PI:.4f})."
        )

    # --- Bosonization ---
    bos = bosonization_alpha(beta_squared)
    g_thirring = bos.g_thirring
    alpha_bare = bos.alpha_bare

    if g_thirring <= 0.0 or alpha_bare <= 0.0:
        # Non-physical regime: return trivial result.
        return TwoLoopResult(
            beta_squared=beta_squared,
            alpha_bare=0.0,
            z_photon_2loop=1.0,
            alpha_2loop=0.0,
            inv_alpha_2loop=float("inf"),
        )

    # Möbius bundle coupling
    g_yukawa = g_thirring * C1_MOBIUS   # = g_thirring / 2

    # Loop coupling expansion parameter λ = g²/π²
    lambda_coupling = (g_yukawa**2) / (PI**2)

    # --- 1-loop integral ---
    i_bundle = compute_i_bundle(
        m_kink_xi=M_KINK_XI,
        x_cutoff=X_CUTOFF,
        n_x=N_X_FINE,
        xi_reg=XI_REG,
        mu_ir=MU_IR_2L,
    )

    # --- 2-loop sunset ---
    i_sunset = compute_i_sunset()
    delta_z_sunset = lambda_coupling**2 * i_sunset

    # --- Figure-8 ---
    i_fig8 = compute_i_figure8(i_bundle)
    delta_z_fig8 = lambda_coupling**2 * i_fig8

    # --- Continuum correction ---
    i_cont = compute_i_continuum()
    delta_z_continuum = lambda_coupling * i_cont

    # --- Total Z at two loops ---
    z_2loop = (
        1.0
        + lambda_coupling * i_bundle          # 1-loop
        + delta_z_sunset                       # 2-loop sunset (1PI)
        + delta_z_fig8                         # 2-loop figure-8
        + delta_z_continuum                    # 1-loop continuum
    )

    if z_2loop <= 0.0:
        z_2loop = 1.0e-8

    # --- Renormalized α ---
    alpha_2loop = alpha_bare / z_2loop
    inv_alpha_2loop = 1.0 / alpha_2loop

    return TwoLoopResult(
        beta_squared=beta_squared,
        alpha_bare=alpha_bare,
        z_photon_2loop=z_2loop,
        alpha_2loop=alpha_2loop,
        inv_alpha_2loop=inv_alpha_2loop,
    )


def solve_beta_fit_to_alpha_target() -> TwoLoopResult | None:
    """Fit beta^2 so the current two-loop map lands on alpha.

    This is a calibration diagnostic, not a prediction.
    """
    lo = 3.91 * PI
    hi = 3.9105 * PI

    def residual(beta_sq: float) -> float:
        return compute_z_photon_2loop(beta_sq).inv_alpha_2loop - INV_ALPHA_TARGET

    f_lo = residual(lo)
    f_hi = residual(hi)
    if f_lo * f_hi > 0.0:
        return None

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f_mid = residual(mid)
        if abs(f_mid) < 1e-12:
            return compute_z_photon_2loop(mid)
        if f_lo * f_mid <= 0.0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid

    return compute_z_photon_2loop(0.5 * (lo + hi))


def run_alpha_two_loop_audit() -> AlphaTwoLoopAudit:
    """Return the compact, reproducible two-loop alpha audit."""
    scan = find_best_z_photon_beta(
        beta_sq_min=0.5 * PI,
        beta_sq_max=4.0 * PI - 0.01,
        n_scan=200,
        n_x=81,
    )
    one_loop_scan_result = scan["best_result"]
    if one_loop_scan_result is None:
        raise RuntimeError("1-loop beta scan returned no physical result")

    baseline = compute_z_photon_2loop(one_loop_scan_result.beta_squared)
    fixed_beta = compute_z_photon_2loop(3.91 * PI)
    fitted_beta = solve_beta_fit_to_alpha_target()

    hw = higgs_w_constraint()
    hw_bos = bosonization_alpha(hw.beta_sq_implied)

    return AlphaTwoLoopAudit(
        target_inv_alpha=INV_ALPHA_TARGET,
        target_alpha=ALPHA_TARGET,
        baseline=baseline,
        fixed_beta=fixed_beta,
        fitted_beta=fitted_beta,
        higgs_w_beta_pi=hw.beta_sq_in_pi,
        higgs_w_ratio=hw.m_h_over_m_w_at_beta_sq,
        higgs_w_g_thirring=hw_bos.g_thirring,
        higgs_w_has_positive_alpha=hw_bos.g_thirring > 0.0 and hw_bos.alpha_bare > 0.0,
        conclusion=(
            "No first-principles alpha derivation: the fixed-beta result is close "
            "but not equal to CODATA alpha, exact agreement requires fitting beta, "
            "and the independent Higgs/W beta gives no positive alpha in the current map."
        ),
    )
