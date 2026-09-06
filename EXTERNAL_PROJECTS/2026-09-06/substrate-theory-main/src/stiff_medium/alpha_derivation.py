"""Audit the current attempt to derive alpha from the substrate model.

This module keeps the reproducible diagnostics motivated by spec §§18.9,
18.34, 18.45, 18.46, 18.48. None of these is a complete derivation; the spec
is explicit that a rigorous result requires multi-loop bundle field theory
(§18.48.7 item 1). What we do here is:

  1. Compute everything the one-loop Coleman-bosonization route gives us.
  2. Invert the RG equations to ask "which bare coupling flows to α(0) = 1/137?"
  3. Use the breather/W-mass ratio m_H/m_W = 2 sin(β²/16) to constrain β².
  4. Combine the constraints and look for a self-consistent β².

The module is deliberately honest about what is and is not derived:

  - α(0) = 1/137.035999177 is a target number we want to predict.
  - Coleman bosonization relates β² ↔ g_Thirring ↔ α_bare, but gives α_bare ≈ 0
    near the free-fermion point — far from 1/137.
  - RG running can tell us what α_bare must be at the substrate scale, but that
    is a CONSTRAINT on β², not a derivation of it.
  - The Möbius bundle fixes the TOPOLOGY (charge quantisation, spin-½) but does
    NOT fix the MAGNITUDE of α without additional Lagrangian dynamics.
  - The breather formula pins β² ≈ 4.54 π from the observed m_H/m_W = 1.558, but
    that β² gives α_bare far from 1/137 via naive bosonization.

CONCLUSION (as of this computation): no single naive β² satisfies both
  α(0) = 1/137 AND m_H/m_W = 1.558 simultaneously from the leading-order
  bosonization relation.  The gap between the two constraints is large — not a
  small numerical accident.  What is needed to close it is stated explicitly at
  the end of the module (multi-loop renormalisation on the Möbius bundle, not
  perturbative patches).

References
----------
    spec §18.9   : α = e²/(Kξ⁴) dimensional route
    spec §18.34  : structural correspondence to QED → RG running inherited
    spec §18.45  : encompassing Lagrangian
    spec §18.46  : derived constants from substrate primitives
    spec §18.48  : breather masses M_n = 2 M_K sin(nβ²/16)
    Coleman 1975 : sine-Gordon / Thirring duality, Phys. Rev. D 11, 2088
    Dashen-Hasslacher-Neveu 1975: breather spectrum

Modules used
------------
    stiff_medium.mobius_bundle : coleman_bosonization_g
    stiff_medium.rg_running    : RGRunning, M_E_GEV, run_alpha_to_M_Z
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

from stiff_medium.mobius_bundle import coleman_bosonization_g
from stiff_medium.physical_constants import (
    ALPHA_THOMSON_CODATA_2022,
    INV_ALPHA_THOMSON_CODATA_2022,
)
from stiff_medium.rg_running import (
    RGRunning,
    M_E_GEV,
    run_alpha_to_M_Z,
)

# ---------------------------------------------------------------------------
# Physical and model constants
# ---------------------------------------------------------------------------

PI: Final[float] = math.pi

# Target values we are trying to derive
ALPHA_TARGET: Final[float] = ALPHA_THOMSON_CODATA_2022
INV_ALPHA_TARGET: Final[float] = INV_ALPHA_THOMSON_CODATA_2022

# W-boson mass [GeV]
M_W_GEV: Final[float] = 80.379
# Higgs-boson mass [GeV] — ATLAS 2025 most-precise measurement (spec §18.52.4)
M_H_GEV: Final[float] = 125.22
# Observed Higgs-to-W mass ratio
M_H_OVER_M_W_MEASURED: Final[float] = M_H_GEV / M_W_GEV  # ≈ 1.558

# Kink mass estimate (per spec §18.22, numerical analysis consistent with observed α and m_e)
M_KINK_GEV: Final[float] = 27.0

# Substrate scale — conjectured scale at which the kink crystallises.
# Identifying with the kink mass: Q_substrate ~ M_kink (a natural UV cutoff for the
# effective theory before the kink forms).
Q_SUBSTRATE_GEV: Final[float] = M_KINK_GEV

# ---------------------------------------------------------------------------
# Route 1 — Coleman bosonization: β² → g_Thirring → α_bare
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BosonizationResult:
    """Bare-alpha map from Coleman bosonization."""

    g_thirring: float
    alpha_bare: float


def bosonization_alpha(beta_squared: float) -> BosonizationResult:
    """Compute α_bare from Coleman bosonization at the given β².

    Derivation chain (spec §18.9, Coleman 1975):

        g_Thirring = π(4π/β² − 1)          [Coleman duality]
        α_bare     = g / π²                  [coupling extraction, convention 1]

    At β² = 4π (free fermion): g = 0, both α_bare vanish — no interaction.
    At β² < 4π: g > 0, repulsive Thirring regime; α_bare > 0.
    At β² > 4π: g < 0, attractive; unphysical for α (clamped to 0).

    Args:
        beta_squared: Sine-Gordon coupling constant β².  Must be positive.

    Returns:
        BosonizationResult with g_Thirring and the convention-1 bare alpha.

    Raises:
        ValueError: If beta_squared is not positive.
    """
    if beta_squared <= 0.0:
        raise ValueError(f"beta_squared must be positive; got {beta_squared!r}")

    g = coleman_bosonization_g(beta_squared)

    if g <= 0.0:
        alpha_bare = 0.0
    else:
        alpha_bare = g / (PI**2)

    return BosonizationResult(
        g_thirring=g,
        alpha_bare=alpha_bare,
    )


def beta_sq_for_alpha_target_bosonization(
    alpha_target: float = ALPHA_TARGET,
) -> float:
    """Invert the bosonization relation to find β² that gives α_target.

    Solves α_target = g/π² for g, then inverts
    g = π(4π/β² − 1) to find β².

    This gives the β² that would be NEEDED if naive one-loop bosonization
    were the complete story.  The spec is clear it is not — but computing
    this value makes the gap explicit.

    Args:
        alpha_target: Target α value.  Default: CODATA 2022 alpha(0).

    Returns:
        β² that solves the equation.
    """
    g_needed = alpha_target * PI**2
    # g = π(4π/β² − 1) → 4π/β² = g/π + 1 → β² = 4π / (g/π + 1)
    beta_sq = 4.0 * PI / (g_needed / PI + 1.0)
    return beta_sq


# ---------------------------------------------------------------------------
# Route 2 — RG running: find α_bare at substrate scale that flows to α(0)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RGRouteResult:
    """RG-required bare alpha and round-trip check."""

    alpha_bare_at_substrate: float
    inv_alpha_bare: float
    alpha_at_me: float


def rg_route(
    q_substrate_gev: float = Q_SUBSTRATE_GEV,
    alpha_target: float = ALPHA_TARGET,
) -> RGRouteResult:
    """Find the α_bare value at the substrate scale that flows to α(0)=1/137.

    Strategy: the RGRunning object runs α from any reference (Q_ref, α_ref)
    to any target Q.  We want to run DOWNWARD: from Q_substrate to m_e.

    Key insight: since the one-loop QED beta function is known analytically,
    we can simply set the REFERENCE to m_e (where α = 1/137) and ask for
    α at Q_substrate by running UPWARD.  That α(Q_substrate) is exactly the
    α_bare that the substrate must supply.

    This is NOT a derivation of 1/137 — it's a CONSTRAINT on what α_bare
    must be.  We then convert that α_bare to a β² via the bosonization
    relation and ask whether that β² makes sense.

    Note on the round-trip residual: this should be numerically zero apart
    from floating-point roundoff.  The RG step is implemented symmetrically
    across fermion thresholds and across the hadronic VP step.

    Args:
        q_substrate_gev: Substrate energy scale [GeV].  Default: M_kink ≈ 27 GeV.
        alpha_target: Low-energy α to match.  Default: CODATA 2022 alpha(0).

    Returns:
        RGRouteResult with the full diagnostics.
    """
    rg = RGRunning(q_ref_gev=M_E_GEV, alpha_ref=alpha_target)

    # Run FROM m_e upward to Q_substrate → this gives α(Q_substrate)
    alpha_subst = rg.alpha_at_scale(q_substrate_gev)
    inv_alpha_subst = 1.0 / alpha_subst

    # Also compute the round-trip: run down from Q_substrate to m_e
    rg_down = RGRunning(q_ref_gev=q_substrate_gev, alpha_ref=alpha_subst)
    alpha_at_me_roundtrip = rg_down.alpha_at_scale(M_E_GEV)

    return RGRouteResult(
        alpha_bare_at_substrate=alpha_subst,
        inv_alpha_bare=inv_alpha_subst,
        alpha_at_me=alpha_at_me_roundtrip,
    )


# ---------------------------------------------------------------------------
# Route 3 — Higgs/W mass ratio: m_H/m_W = 2 sin(β²/16) → β² → α
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HiggsWConstraintResult:
    """Beta implied by the Higgs/W breather-ratio constraint."""

    m_h_over_m_w_at_beta_sq: float
    beta_sq_implied: float
    beta_sq_in_pi: float


def higgs_w_constraint(
    n_breather: int = 1,
    m_h_gev: float = M_H_GEV,
    m_w_gev: float = M_W_GEV,
) -> HiggsWConstraintResult:
    """Constrain β² from the observed m_H/m_W ratio via the breather formula.

    Uses DHN breather formula M_n = 2 M_K sin(n β²/16) with:
      - Breather mode n (default 1 for fundamental Higgs-like state)
      - M_1 identified with Higgs mass m_H
      - M_K identified with W-boson mass m_W (the electroweak scale kink)

    Args:
        n_breather: Breather mode number (1 = first/lightest state).
        m_h_gev: Higgs-boson mass [GeV].  Default: 125.22 (ATLAS 2025).
        m_w_gev: W-boson mass [GeV].  Default: 80.379.
    Returns:
        HiggsWConstraintResult with all diagnostics.

    Raises:
        ValueError: If m_h_gev > 2 * n_breather * m_w_gev (arcsin argument > 1).
    """
    ratio = m_h_gev / m_w_gev
    arcsin_arg = ratio / 2.0

    if arcsin_arg > 1.0:
        raise ValueError(
            f"m_H / (2 m_W) = {arcsin_arg:.4f} > 1: no real β² satisfies "
            f"the breather formula for n={n_breather}."
        )

    beta_sq = 16.0 * math.asin(arcsin_arg / n_breather)

    # Verify round-trip
    ratio_predicted = 2.0 * n_breather * math.sin(beta_sq / 16.0)

    return HiggsWConstraintResult(
        m_h_over_m_w_at_beta_sq=ratio_predicted,
        beta_sq_implied=beta_sq,
        beta_sq_in_pi=beta_sq / PI,
    )


# ---------------------------------------------------------------------------
# Route 4 — Self-consistency scan: find β² satisfying all constraints
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConsistencyScanSummary:
    """Aggregate result of the alpha/Higgs-W beta scan."""

    n_fully_consistent: int
    best_inv_alpha: float


@dataclass(frozen=True)
class AlphaDerivationAudit:
    """Compact conclusion for the current alpha-derivation attempt."""

    target_inv_alpha: float
    target_alpha: float
    inv_alpha_mz_computed: float
    inv_alpha_mz_measured: float
    rg_roundtrip_inv_alpha_me: float
    rg_roundtrip_delta_alpha: float
    inv_alpha_bare_required_27gev: float
    beta_alpha_bare_pi: float
    beta_rg_pi: float
    beta_higgs_w_pi: float
    beta_gap_pi: float
    higgs_w_ratio: float
    higgs_w_g_thirring: float
    higgs_w_has_positive_alpha: bool
    n_fully_consistent_scan_points: int
    best_scan_inv_alpha: float
    conclusion: str


def summarize_consistency_scan(
    beta_sq_min: float = 0.5 * PI,
    beta_sq_max: float = 8.0 * PI,
    n_steps: int = 2000,
    q_substrate_gev: float = Q_SUBSTRATE_GEV,
    alpha_tol_inv_alpha: float = 5.0,
    ratio_tol: float = 0.05,
) -> ConsistencyScanSummary:
    """Return the two scan facts used by the audit."""
    if n_steps < 2:
        raise ValueError(f"n_steps must be at least 2; got {n_steps!r}")

    best_delta_inv = float("inf")
    best_inv_alpha = float("inf")
    n_fully_consistent = 0

    for i in range(n_steps):
        beta_sq = beta_sq_min + i * (beta_sq_max - beta_sq_min) / (n_steps - 1)
        bos = bosonization_alpha(beta_sq)
        a_bare_c1 = bos.alpha_bare

        if a_bare_c1 > 1e-10:
            rg = RGRunning(q_ref_gev=q_substrate_gev, alpha_ref=a_bare_c1)
            a_me = rg.alpha_at_scale(M_E_GEV)
        else:
            a_me = 0.0

        inv_a_me = 1.0 / a_me if a_me > 1e-20 else float("inf")
        delta_inv = abs(inv_a_me - INV_ALPHA_TARGET) if math.isfinite(inv_a_me) else float("inf")

        if a_me > 0.0 and delta_inv < best_delta_inv:
            best_delta_inv = delta_inv
            best_inv_alpha = inv_a_me

        breather_ratio = 2.0 * math.sin(beta_sq / 16.0)
        delta_rat = abs(breather_ratio - M_H_OVER_M_W_MEASURED)
        if delta_inv < alpha_tol_inv_alpha and delta_rat < ratio_tol:
            n_fully_consistent += 1

    return ConsistencyScanSummary(
        n_fully_consistent=n_fully_consistent,
        best_inv_alpha=best_inv_alpha,
    )


def run_alpha_derivation_audit(n_steps: int = 2000) -> AlphaDerivationAudit:
    """Return the compact, reproducible alpha-derivation audit."""
    rg_check = run_alpha_to_M_Z()
    rg27 = rg_route(Q_SUBSTRATE_GEV)

    beta_alpha_bare = beta_sq_for_alpha_target_bosonization(ALPHA_TARGET)
    beta_rg = beta_sq_for_alpha_target_bosonization(rg27.alpha_bare_at_substrate)

    hw = higgs_w_constraint()
    hw_bos = bosonization_alpha(hw.beta_sq_implied)

    scan = summarize_consistency_scan(n_steps=n_steps)

    return AlphaDerivationAudit(
        target_inv_alpha=INV_ALPHA_TARGET,
        target_alpha=ALPHA_TARGET,
        inv_alpha_mz_computed=float(rg_check["inv_alpha_MZ_predicted"]),
        inv_alpha_mz_measured=float(rg_check["inv_alpha_MZ_measured"]),
        rg_roundtrip_inv_alpha_me=1.0 / rg27.alpha_at_me,
        rg_roundtrip_delta_alpha=rg27.alpha_at_me - ALPHA_TARGET,
        inv_alpha_bare_required_27gev=rg27.inv_alpha_bare,
        beta_alpha_bare_pi=beta_alpha_bare / PI,
        beta_rg_pi=beta_rg / PI,
        beta_higgs_w_pi=hw.beta_sq_in_pi,
        beta_gap_pi=abs(beta_alpha_bare / PI - hw.beta_sq_in_pi),
        higgs_w_ratio=hw.m_h_over_m_w_at_beta_sq,
        higgs_w_g_thirring=hw_bos.g_thirring,
        higgs_w_has_positive_alpha=hw_bos.g_thirring > 0.0 and hw_bos.alpha_bare > 0.0,
        n_fully_consistent_scan_points=scan.n_fully_consistent,
        best_scan_inv_alpha=scan.best_inv_alpha,
        conclusion=(
            "No first-principles alpha derivation: the alpha-matching beta values "
            "are constraints/calibrations, and the independent Higgs/W beta is "
            "outside the positive-alpha branch of the current Coleman map."
        ),
    )
