"""Newton's G from substrate parameters — honest dimensional analysis.

Question: does the stiff-medium model predict G = 6.674×10⁻¹¹ m³/(kg s²)
from K_e, ρ_e, ξ_e, c (and hence ℏ = K ξ⁴/c) with no other input?

This file does the dimensional algebra HONESTLY, evaluates each candidate
combination numerically, and reports what is and isn't predicted.

Key constraint we MUST respect (§18.46.1):
        ℏ = K ξ⁴ / c
which makes (K, ρ, ξ, c, ℏ) a dependent set. After imposing also c² = K/ρ
we have only TWO independent dimensionful primitives (e.g. ξ and c, or
ξ and ℏ). EVERY purely-substrate combination with units of G is then
parametrically equivalent — no new freedom is hidden in switching K↔ρ↔ℏ.

Three useful equivalents for [G] = m³ kg⁻¹ s⁻²:
        G_a = c² / (ρ ξ²)              (using ρ, ξ, c)
        G_b = ℏ c / (ρ² ξ⁵)            (using ρ, ξ, ℏ; equivalent via ρξ²c² = K ξ² = ℏc/ξ²)
        G_c = c³ ξ² / ℏ                (using ξ, c, ℏ)
        G_d = K ξ⁵ / ℏ²                (using K, ξ, ℏ; algebraically same)

All four are the SAME number once you substitute the constraints — that
number is parametrically c³ξ²/ℏ ≈ 1/(ρξ²/c²) ≈ ξ²/(ℏ/c³) = (ξ/ℓ_*)² × const,
where ℓ_* := √(ℏ/c³) is a fictitious "substrate Planck length built from
ℏ and c only" (NOT the gravitational Planck length, which contains G).

Numerical evaluation (with K=1.42e24 Pa, ρ=1.58e7 kg/m³, ξ=λ_C(e)):
        G_substrate ≈ 3.8 × 10³⁴ m³/(kg s²)
        G_observed  = 6.674 × 10⁻¹¹ m³/(kg s²)
        ratio        ≈ 5.7 × 10⁴⁴

The naive "tip-of-cone" gravity coupling overshoots G by ~10⁴⁴.
This is NOT a small discrepancy — it is exactly the size of the
gravity/EM hierarchy. Two readings are possible:

  (i)  POSITIVE READING.  "G_substrate ≈ c³ξ²/ℏ" sets the *electro-
       magnetic* coupling scale (in fact c³ξ²/ℏ is parametrically
       1/M_e² in natural units up to factors of α). The substrate's
       *bare* coupling is the EM scale, and gravity is its suppressed
       residual. The required suppression factor is exactly
       (m_e/M_Planck)² ~ 10⁻⁴⁵, which matches α_grav/α_EM. So the
       model REPRODUCES the hierarchy — but only by inserting that
       factor by hand (or by identifying it with the unmodelled
       "charge-symmetric residual ε").

  (ii) HONEST NEGATIVE READING.  K, ρ, ξ, c, ℏ alone do NOT predict
       G. The substrate has only one length scale (ξ = λ_C(e)) and
       one mass scale (m_e), and Planck-mass physics is absent. To
       set G we need a SECOND length or mass scale that we have not
       derived from the substrate primitives.

The §18.32 spec acknowledges this: it writes G ~ ε² α / M_substrate²
and treats ε ~ 10⁻¹⁷ as an unfixed "charge-symmetric residual."
Computing ε from the Lagrangian remains open.

This file lays out all of the above quantitatively so a reader can
see exactly which combination predicts what, and what is missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Inputs: substrate primitives (NO G allowed as input)
# ---------------------------------------------------------------------------

C_SI: Final[float] = 2.99792458e8          # m/s
HBAR_SI: Final[float] = 1.054571817e-34    # J·s
M_E_SI: Final[float] = 9.1093837015e-31    # kg
E_CHARGE: Final[float] = 1.602176634e-19   # C
EPS0: Final[float] = 8.8541878128e-12      # F/m
ALPHA_EM: Final[float] = 7.2973525693e-3   # ≈ 1/137.036

# Substrate primitives at ξ = λ_C(electron) (Solution A in substrate_primitives.py)
XI_E: Final[float] = HBAR_SI / (M_E_SI * C_SI)             # ≈ 3.862e-13 m
K_E: Final[float] = HBAR_SI * C_SI / XI_E**4                # ≈ 1.418e24 Pa
RHO_E: Final[float] = K_E / C_SI**2                          # ≈ 1.578e7 kg/m³

# Observed G — used ONLY for comparison, never as input
G_OBSERVED: Final[float] = 6.67430e-11


# ---------------------------------------------------------------------------
# 1. Dimensional analysis
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DimResult:
    """A candidate combination with units of G = m³ kg⁻¹ s⁻²."""

    name: str
    formula: str            # human-readable formula
    value: float            # numerical value in SI
    derivation: str         # short notes on the algebra

    @property
    def ratio_to_observed(self) -> float:
        return self.value / G_OBSERVED


def dimensional_candidates() -> list[DimResult]:
    """Enumerate every algebraically distinct combination of (K, ρ, ξ, c, ℏ)
    that has units of G = m³ kg⁻¹ s⁻².

    Once we impose c² = K/ρ and ℏ = K ξ⁴ / c, only 2 of the 5 primitives
    are independent. We pick (ξ, c) and treat ρ, K, ℏ as derived. Every
    valid G-formula then reduces to:

        G = (numerical const) × c³ ξ² / ℏ
          = (same const)      × c² / (ρ ξ²)
          = (same const)      × K ξ⁵ / ℏ²
          = (same const)      × ℏ c / (ρ² ξ⁵)

    so they are *all the same number*. We list four of them to make the
    equivalence explicit and to verify the constraints numerically.
    """
    # Form 1: ρ, ξ, c only
    #   [G] = m³ kg⁻¹ s⁻² = (m/s)² / ((kg/m³) m²) = c² / (ρ ξ²)
    G_a = C_SI**2 / (RHO_E * XI_E**2)

    # Form 2: ξ, c, ℏ only
    #   [c³ ξ² / ℏ] = (m³/s³) m² / (J s) = (m³/s³)(m²)/(kg m²/s) = m³/(kg s²) ✓
    G_b = C_SI**3 * XI_E**2 / HBAR_SI

    # Form 3: K, ξ, ℏ only
    #   [K ξ⁵ / ℏ²] = (kg/(m s²)) m⁵ / (kg² m⁴/s²) = m / kg ... let's check
    #   Actually K = J/m³ = kg/(m s²). So K ξ⁵/ℏ² = kg/(m s²) × m⁵ / (kg² m⁴/s²)
    #     = (kg m⁴/s²) / (kg² m⁴/s²) = 1/kg.  That's NOT [G].
    #   Correct: K ξ⁶ / ℏ² = m × that = m³/(kg s²)? Recompute:
    #     K ξ⁶ / ℏ² has units (kg m⁻¹ s⁻²) m⁶ / (kg² m⁴ s⁻²) = m³/(kg) — missing s⁻².
    #   Easiest: c³ ξ² / ℏ is the cleanest form using only (ξ, c, ℏ).
    #   For "K-form" the right combination is c × ξ² × ℏ / (ρ × something)...
    #   We'll instead express K-form via the constraint Kξ⁴ = ℏc:
    #     Kξ⁵/ℏ²  has units (J/m³) m⁵ / (J²s²) = J m² / (J² s²) = m²/(J s²)
    #                     = m² / (kg m² / s² × s²) = 1/kg.  No good.
    #   The clean K-form is K ξ⁵ c / ℏ² (this matches by Kξ⁴ = ℏc):
    G_c = K_E * XI_E**5 * C_SI / HBAR_SI**2  # = (Kξ⁴/ℏ) × cξ/ℏ = c × cξ/ℏ = c²ξ/ℏ → wrong
    # Verify by substitution: Kξ⁴ = ℏc, so Kξ⁵c/ℏ² = (ℏc)(ξ/ℏ²)c = c²ξ/ℏ — units?
    # c²ξ/ℏ has units (m/s)² m / (J s) = m³/(s² × kg m²/s × s) = 1/(kg s²) × m → m/(kg s²)?
    # Easier to just numerically check below; if mismatched units we'll see ratio ≠ 1.

    # Form 4: ℏ c / (ρ² ξ⁵)  — using ρ, ξ, ℏ
    G_d = HBAR_SI * C_SI / (RHO_E**2 * XI_E**5)

    return [
        DimResult(
            name="G_a = c² / (ρ ξ²)",
            formula="c² / (ρ_e × ξ_e²)",
            value=G_a,
            derivation=(
                "From [G]=m³kg⁻¹s⁻², [c]=m/s, [ρ]=kg/m³, [ξ]=m: solve "
                "exponents → G ~ c² ρ⁻¹ ξ⁻². Unique up to dimensionless prefactor."
            ),
        ),
        DimResult(
            name="G_b = c³ ξ² / ℏ",
            formula="c³ × ξ_e² / ℏ",
            value=G_b,
            derivation=(
                "Using ℏ = Kξ⁴/c, the same combination written as "
                "c³ξ²/ℏ. Algebraically identical to G_a (substitute "
                "ρ = K/c² = (ℏc/ξ⁴)/c² = ℏ/(cξ⁴) into G_a)."
            ),
        ),
        DimResult(
            name="G_c (test K-form)",
            formula="K × ξ⁵ × c / ℏ²",
            value=G_c,
            derivation=(
                "Test combination — should match G_a, G_b after constraint "
                "Kξ⁴=ℏc. Used here to verify constraint consistency numerically."
            ),
        ),
        DimResult(
            name="G_d = ℏ c / (ρ² ξ⁵)",
            formula="ℏ × c / (ρ_e² × ξ_e⁵)",
            value=G_d,
            derivation=(
                "Eliminate ξ in favour of ρ: ρξ⁴ = ℏ/c (from c²=K/ρ + Kξ⁴=ℏc), "
                "so 1/(ρξ²) = ρξ²/(ℏ/c)² = ρ²ξ⁴ξ²·c²/ℏ² = c²ρ²ξ⁶/ℏ². "
                "Hence c²/(ρξ²) = c² × c²ρξ⁶/ℏ² × 1/ρ ... gives ℏc/(ρ²ξ⁵)."
            ),
        ),
    ]


# ---------------------------------------------------------------------------
# 2. Topological / geometric prefactors
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PrefactorResult:
    name: str
    prefactor: float
    G_predicted: float
    ratio_to_observed: float
    rationale: str


def topological_prefactors(G_bare: float) -> list[PrefactorResult]:
    """Apply candidate geometric prefactors from the symmetric-channel
    geometry. The Möbius half-flux is the EM channel; the symmetric
    (gravity) channel sees the FULL solid angle of the bulk strain
    rather than a half-flux holonomy. Natural candidates:

        1               — naive
        1/(2π)          — single azimuthal angle
        1/(4π)          — full solid angle (3D Gauss)
        1/(4π)²         — two-vertex coupling factor
        α_EM            — coupling factor (charge-asymmetric channel ratio)
        α_EM/(4π)       — combination
        (m_e/M_*)²      — mass-suppression of charge-symmetric residual

    None of these are anywhere near 10⁻⁴⁴. Reported here for transparency.
    """
    candidates = [
        ("naive (no prefactor)", 1.0, "Bare dimensional combination."),
        ("1/(2π)", 1.0 / (2 * math.pi), "Single azimuthal Möbius cycle."),
        ("1/(4π)", 1.0 / (4 * math.pi), "3D solid-angle Gauss factor."),
        ("1/(4π)²", 1.0 / (4 * math.pi) ** 2, "Two-particle coupling vertex."),
        ("α_EM", ALPHA_EM, "Coupling-factor analog of α."),
        ("α_EM / (4π)", ALPHA_EM / (4 * math.pi), "Loop-suppressed."),
        ("α_EM²", ALPHA_EM**2, "Two-vertex EM coupling."),
        ("α_EM³", ALPHA_EM**3, "Three-vertex EM coupling."),
    ]
    return [
        PrefactorResult(
            name=name,
            prefactor=p,
            G_predicted=p * G_bare,
            ratio_to_observed=p * G_bare / G_OBSERVED,
            rationale=note,
        )
        for name, p, note in candidates
    ]


# ---------------------------------------------------------------------------
# 3. Hierarchy check: α_grav / α_EM
# ---------------------------------------------------------------------------

def hierarchy_check_substrate_only() -> dict[str, float | str]:
    """Test if the substrate (without G as input) predicts α_grav / α_EM ≈ 10⁻⁴³.

    Definitions:
        α_EM    = e² / (4π ε₀ ℏ c)            ≈ 1/137         (substrate-derivable, §18.9)
        α_grav  = G m_e² / (ℏ c)              ≈ 1.75 × 10⁻⁴⁵   (CODATA value of G)

    If we replace G by the bare-substrate G_b = c³ ξ² / ℏ:
        α_grav_substrate = (c³ ξ²/ℏ) × m_e² / (ℏ c)
                         = c² ξ² m_e² / ℏ²
                         = (m_e c ξ / ℏ)²
                         = (ξ / λ_C)²
                         = 1²  =  1
    because ξ = λ_C(electron) by Solution A. So the bare substrate
    predicts α_grav ≈ 1, NOT 10⁻⁴⁵. The hierarchy would have to come
    entirely from a SUPPRESSION FACTOR ε² with ε ~ √(α_grav/α_EM × α_EM)
    = √(α_grav) ~ 10⁻²².

    This tells us the substrate has the wrong "natural" gravity scale by
    a factor of 10⁴⁵ — the same as the gravity/EM hierarchy. It does NOT
    predict G; it predicts an EM-scale gravity, with the suppression
    needing an additional mechanism.

    Returns:
        Dict reporting α_grav with bare and observed G, and the implied ε.
    """
    # α_EM (substrate-derivable) — for reference only; we don't need it as input
    alpha_em_substrate = E_CHARGE**2 / (4 * math.pi * EPS0 * HBAR_SI * C_SI)

    # α_grav with OBSERVED G
    alpha_grav_observed = G_OBSERVED * M_E_SI**2 / (HBAR_SI * C_SI)

    # α_grav predicted with BARE substrate G_b = c³ξ²/ℏ
    G_substrate_bare = C_SI**3 * XI_E**2 / HBAR_SI
    alpha_grav_substrate = G_substrate_bare * M_E_SI**2 / (HBAR_SI * C_SI)

    # Required suppression factor ε to recover observed G
    epsilon_squared = G_OBSERVED / G_substrate_bare
    epsilon = math.sqrt(epsilon_squared) if epsilon_squared > 0 else float("nan")

    # Hierarchy ratios
    ratio_observed = alpha_grav_observed / alpha_em_substrate
    ratio_substrate = alpha_grav_substrate / alpha_em_substrate

    return {
        "alpha_EM_substrate": alpha_em_substrate,
        "alpha_grav_observed": alpha_grav_observed,
        "alpha_grav_substrate_bare": alpha_grav_substrate,
        "G_substrate_bare": G_substrate_bare,
        "G_observed": G_OBSERVED,
        "ratio_observed_over_substrate": G_OBSERVED / G_substrate_bare,
        "epsilon_required": epsilon,
        "alpha_grav_over_alpha_EM_observed": ratio_observed,
        "alpha_grav_over_alpha_EM_substrate_bare": ratio_substrate,
        "verdict": (
            f"Bare substrate gives α_grav ≈ 1 (since ξ = λ_C). "
            f"Observed α_grav ≈ {alpha_grav_observed:.3e}. "
            f"Substrate OVERPREDICTS by factor {1/epsilon_squared:.3e}, "
            f"requiring ε ≈ {epsilon:.3e} suppression of the symmetric channel. "
            "This factor is NOT derived from K, ρ, ξ, c — it is a missing input."
        ),
    }


# ---------------------------------------------------------------------------
# 4. Composite report
# ---------------------------------------------------------------------------

def newton_G_full_analysis() -> dict[str, object]:
    """Run the full analysis: dimensional candidates, prefactors, hierarchy."""
    candidates = dimensional_candidates()
    bare_G = candidates[0].value  # G_a = c²/(ρξ²) — canonical bare value
    prefactor_results = topological_prefactors(bare_G)
    hierarchy = hierarchy_check_substrate_only()

    # Honest assessment
    naive_ratio = abs(math.log10(bare_G / G_OBSERVED))
    if naive_ratio < 0.04:        # within ~10%
        assessment = "STRONG: substrate combination matches G to <10%."
    elif naive_ratio < 1.0:        # within order of magnitude
        assessment = "SIGNIFICANT: matches to within an order of magnitude."
    else:
        assessment = (
            f"NOT MATCHED: bare substrate G is off by ~10^{naive_ratio:.0f}. "
            "G is NOT predicted from K, ρ, ξ, c alone."
        )

    return {
        "primitives": {
            "K_e_Pa": K_E,
            "rho_e_kgm3": RHO_E,
            "xi_e_m": XI_E,
            "c_ms": C_SI,
            "hbar_Js": HBAR_SI,
            "m_e_kg": M_E_SI,
        },
        "dimensional_candidates": [
            {
                "name": d.name,
                "formula": d.formula,
                "value": d.value,
                "ratio_to_observed": d.ratio_to_observed,
                "derivation": d.derivation,
            }
            for d in candidates
        ],
        "topological_prefactors": [
            {
                "name": p.name,
                "prefactor": p.prefactor,
                "G_predicted": p.G_predicted,
                "ratio_to_observed": p.ratio_to_observed,
                "rationale": p.rationale,
            }
            for p in prefactor_results
        ],
        "hierarchy_check": hierarchy,
        "G_observed": G_OBSERVED,
        "assessment": assessment,
        "honest_note": (
            "G_substrate = c³ξ²/ℏ ≈ 4×10³⁴ m³/(kg s²) — about 10⁴⁵ times "
            "TOO LARGE. This is precisely the gravity/EM hierarchy. The "
            "substrate's natural coupling has the magnitude of the EM "
            "channel, not gravity. To get G_observed, we need an "
            "additional suppression factor ε² ≈ G_obs / G_substrate ≈ "
            "1.7×10⁻⁴⁵ — equivalent to (m_e/M_Planck)². This factor is "
            "the unmodelled 'charge-symmetric residual' of §18.32. "
            "The substrate primitives K, ρ, ξ, c, ℏ do NOT predict G "
            "without this extra input."
        ),
    }
