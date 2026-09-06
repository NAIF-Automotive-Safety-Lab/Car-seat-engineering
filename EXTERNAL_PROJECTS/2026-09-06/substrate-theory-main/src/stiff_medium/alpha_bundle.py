"""Direct bundle/back-reaction derivation of alpha — bypassing sine-Gordon bosonization.

The existing route in `alpha_derivation.py` (Coleman bosonization →
Thirring g → α_bare) cannot simultaneously satisfy:

    RG constraint:    α_bare(27 GeV) ≈ 1/129    (must run down to α(0)=1/137)
    Higgs/W constraint: β² = 4.55π             (from m_H/m_W = 2 sin(β²/16))

because β² = 4.55π is in the attractive Thirring regime where g ≤ 0,
giving α_bare = 0 — far from the required 1/129.

This module pursues a *different* route:

    α = q²_substrate / (4π K ℏ c)

with ℏ = K ξ⁴/c (the framework's working derivation, MODEL.md §1.4).
This reduces to a clean closed form:

    α = β² / (4π)

where β = q_substrate / (K ξ²) is the dimensionless substrate charge in
Heaviside-Lorentz natural units.

This is the *same* relation as standard QED in HL units (e_HL = √(4πα)
≈ 0.3028), expressed in substrate variables. The framework's derivation
of α therefore reduces to *deriving the dimensionless charge β = 0.3028
from substrate topology* — specifically from the Möbius bundle's
holonomy structure.

This module:
  1. Computes α from a given β (forward direction)
  2. Inverts α(observed) → β_required = 0.3028
  3. Compares β_required against three candidate structural values from
     the Möbius bundle topology
  4. Reports which is closest and what residual the framework still owes

References
----------
  spec §1.4   : ℏ = K ξ⁴/c
  spec §2.2   : back-reaction Coulomb-like potential
  spec §2.4   : Möbius half-flux bundle
  spec §18.9  : original bosonization route (which fails the dual constraint)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

PI: Final[float] = math.pi

# CODATA 2022 fine-structure constant
ALPHA_CODATA: Final[float] = 7.2973525643e-3
INV_ALPHA_CODATA: Final[float] = 1.0 / ALPHA_CODATA  # 137.035999177

# Heaviside-Lorentz natural-units charge: e_HL = sqrt(4π α)
E_HL_OBSERVED: Final[float] = math.sqrt(4.0 * PI * ALPHA_CODATA)


# ---------------------------------------------------------------------------
# Forward and inverse relations
# ---------------------------------------------------------------------------


def alpha_from_beta(beta: float) -> float:
    """Compute α = β²/(4π) from the dimensionless substrate charge β.

    Args:
        beta: Dimensionless substrate charge (q_substrate / (K ξ²)),
            equivalently the Heaviside-Lorentz electric charge.

    Returns:
        Fine-structure constant α.
    """
    return beta * beta / (4.0 * PI)


def beta_from_alpha(alpha: float = ALPHA_CODATA) -> float:
    """Invert: solve β² = 4π α for β.

    Args:
        alpha: Target fine-structure constant. Default: CODATA 2022.

    Returns:
        β = √(4π α). For α(0) = 1/137.036, β = 0.30282...
    """
    return math.sqrt(4.0 * PI * alpha)


# ---------------------------------------------------------------------------
# Candidate structural values for β from Möbius bundle topology
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BetaCandidate:
    """A candidate value of β from substrate topology.

    Attributes:
        name: Identifier of the candidate (e.g., "1/(π·X)", "Möbius half-flux × Y").
        value: Numerical value of β.
        alpha: α = β²/(4π) computed from this β.
        inv_alpha: 1/α.
        residual_pct: Percentage residual vs CODATA α(0).
        rationale: One-line description of the structural origin attempted.
    """

    name: str
    value: float
    alpha: float
    inv_alpha: float
    residual_pct: float
    rationale: str


def _candidate(name: str, value: float, rationale: str) -> BetaCandidate:
    a = alpha_from_beta(value)
    return BetaCandidate(
        name=name,
        value=value,
        alpha=a,
        inv_alpha=1.0 / a,
        residual_pct=100.0 * abs(a - ALPHA_CODATA) / ALPHA_CODATA,
        rationale=rationale,
    )


def candidate_beta_values() -> list[BetaCandidate]:
    """Return a curated list of candidate β values from substrate topology.

    Each candidate is a structurally-motivated guess for what β might be
    from the Möbius bundle + substrate dynamics. None is proven; the goal
    is to identify which forms come *closest* to β_observed = 0.3028
    and therefore guide the next round of derivation work.

    Returns:
        List of BetaCandidate, sorted by |residual_pct| ascending.
    """
    cands: list[BetaCandidate] = [
        # --- Pure π-based forms ---
        _candidate(
            "1/π",
            1.0 / PI,
            "Möbius holonomy phase (π) per traversal, charge inverse",
        ),
        _candidate(
            "1/(π·√(4/π))",
            1.0 / (PI * math.sqrt(4.0 / PI)),
            "Möbius holonomy normalized by free-photon dispersion",
        ),
        _candidate(
            "√(2/π)/√7",
            math.sqrt(2.0 / PI) / math.sqrt(7.0),
            "Substrate strain integral over half-flux bundle, 7-vertex normalisation",
        ),
        # --- Trigonometric forms ---
        _candidate(
            "sin(π/10)",
            math.sin(PI / 10.0),
            "10-fold substrate symmetry (decagonal hex-packing)",
        ),
        _candidate(
            "sin(π/(7+α^-1·0))",  # placeholder, just sin(π/7)
            math.sin(PI / 7.0),
            "7-fold substrate stress-tensor mode",
        ),
        # --- Square-root forms ---
        _candidate(
            "1/√11",
            1.0 / math.sqrt(11.0),
            "11-vertex saddle (5+6 substrate cell topology)",
        ),
        _candidate(
            "1/√(4π+1)",
            1.0 / math.sqrt(4.0 * PI + 1.0),
            "Self-consistent: β² = 1/(4π+1) → α(4π+1) = 1, near-fixed-point",
        ),
        # --- Exact match form (target) ---
        _candidate(
            "β_observed",
            E_HL_OBSERVED,
            "Inverted from α_CODATA = 1/137.036 (target, not derived)",
        ),
    ]
    cands.sort(key=lambda c: c.residual_pct)
    return cands


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AlphaBundleAudit:
    """Compact audit of the bundle-route derivation.

    Attributes:
        beta_observed: β = √(4π α_CODATA) = 0.3028...
        alpha_observed: α(0) from CODATA 2022.
        candidates: Sorted list of candidate β values.
        best_candidate: Closest candidate to β_observed.
        best_residual_pct: |α(best) - α_obs|/α_obs × 100.
        conclusion: Human-readable status.
    """

    beta_observed: float
    alpha_observed: float
    candidates: list[BetaCandidate]
    best_candidate: BetaCandidate
    best_residual_pct: float
    conclusion: str


def audit() -> AlphaBundleAudit:
    """Run the full bundle-route audit and return diagnostics.

    Returns:
        AlphaBundleAudit with sorted candidates and best-match identification.
    """
    cands = candidate_beta_values()
    # Exclude the target itself from "best derived" — it's tautological
    derived_cands = [c for c in cands if c.name != "β_observed"]
    best = min(derived_cands, key=lambda c: c.residual_pct)

    if best.residual_pct < 0.5:
        verdict = (
            f"Best structural candidate ('{best.name}') matches α to {best.residual_pct:.3f}% — "
            "strong evidence this is the right form, pending derivation of why this candidate."
        )
    elif best.residual_pct < 2.0:
        verdict = (
            f"Best structural candidate ('{best.name}') is within {best.residual_pct:.3f}% — "
            "promising, needs one geometric correction or fine-tuning of the substrate cell topology."
        )
    elif best.residual_pct < 10.0:
        verdict = (
            f"Best structural candidate ('{best.name}') off by {best.residual_pct:.3f}% — "
            "pattern suggests right family, missing structural factor."
        )
    else:
        verdict = (
            f"No structural candidate within 10% of β_observed = {E_HL_OBSERVED:.4f}. "
            "The framework owes a fresh derivation route — neither pure-π, trig, "
            "nor square-root forms with low integers reach observed α from substrate alone."
        )

    return AlphaBundleAudit(
        beta_observed=E_HL_OBSERVED,
        alpha_observed=ALPHA_CODATA,
        candidates=cands,
        best_candidate=best,
        best_residual_pct=best.residual_pct,
        conclusion=verdict,
    )


def main() -> None:  # pragma: no cover
    """Print the audit report."""
    a = audit()
    print(f"α_CODATA       = {a.alpha_observed:.10f} = 1/{1/a.alpha_observed:.6f}")
    print(f"β_observed     = √(4π α) = {a.beta_observed:.6f}")
    print()
    print("Candidate β values (sorted by |residual|):")
    print(f"  {'Name':<25s} {'β':>10s} {'1/α':>12s} {'residual %':>12s}  rationale")
    for c in a.candidates:
        marker = "  ← target" if c.name == "β_observed" else ""
        print(f"  {c.name:<25s} {c.value:>10.6f} {c.inv_alpha:>12.4f} {c.residual_pct:>12.4f}  {c.rationale}{marker}")
    print()
    print(f"Conclusion: {a.conclusion}")


if __name__ == "__main__":  # pragma: no cover
    main()
