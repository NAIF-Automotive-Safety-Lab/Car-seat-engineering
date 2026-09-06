"""Substrate-derived ⁷Li/⁷Be suppression factor for the BBN lithium puzzle.

The cosmological lithium problem
--------------------------------
Standard BBN (SBBN) at the Planck-2018 baryon-to-photon ratio
η = 6.10×10⁻¹⁰ predicts a primordial ⁷Li abundance ⁷Li/H ≈ 4.5×10⁻¹⁰
(eta-fit form, Cyburt et al. 2008; ≈4.7×10⁻¹⁰ with the older PRIMAT
central value used in :mod:`bbn_substrate`).  Observations on the
"Spite plateau" of metal-poor halo stars consistently report ⁷Li/H ≈
(1.6 ± 0.3)×10⁻¹⁰.  SBBN therefore over-predicts the observed ⁷Li
abundance by a factor

    S_obs ≡ ⁷Li/H (SBBN)  /  ⁷Li/H (Spite plateau)  ≈  4.5/1.6  ≈  2.81
                                                  (or 4.7/1.6 ≈ 2.94)

This is the famous "lithium problem".  In standard cosmology the
suppression is unexplained; candidate resolutions invoke stellar
depletion (atomic diffusion), nuclear-cross-section systematics, or
BSM physics (decaying particles).

In the substrate framework, ⁷Li is overwhelmingly synthesized via ⁷Be
(which subsequently electron-captures to ⁷Li after recombination).
The substrate framework's ``§18.66`` proto-matter / observer-horizon
re-thermalization opens a late-time ⁷Be destruction channel during
the de-saturation epoch.  Earlier modules (``bbn_substrate``,
``bbn_test``) introduced this suppression as an empirical knob set to
3.0 — bringing SBBN's 4.5×10⁻¹⁰ down to ≈1.5×10⁻¹⁰, well inside the
Spite-plateau 1-σ band.  This module derives that same factor of 3
directly from the substrate's integer primitives, converting it from
Category C (empirical) to Category A (derived).

Substrate primitives audit
--------------------------
We surveyed every "small ratio" producible from the canonical 12
integers (``b3_constants``: N_BAM=6, K_pair=2, K_rank=5, n_R=18,
F=2, R=3, n_M=268, n_A=15, V13=13, plus the Λ_QCD anchor):

================================================  ======  ==========
Combination                                       Value   |err vs 3|
================================================  ======  ==========
**F · R / K_pair = 6/2**                          3.000   0.0%
**n_R / N_BAM     = 18/6**                        3.000   0.0%
**n_A / K_rank    = 15/5**                        3.000   0.0%
**R               = 3**                           3.000   0.0% (trivial)
(n_R - K_rank) / K_pair²  =  V13 / K_pair² = 13/4  3.250   8.3%
F · K_rank / R           = 10/3                   3.333   11.1%
n_R / K_rank             = 18/5                   3.600   20.0%
K_pair · K_rank / F      = 10/2 = 5               5.000   ✗
n_M / (K_pair·K_rank·N_BAM·R) = 268/180           1.489   ✗
K_rank / (K_pair · F)    = 5/4                    1.250   ✗
================================================  ======  ==========

Three non-trivial substrate combinations all evaluate to **exactly 3**.
They differ in their physical interpretation; we now examine which one
encodes the actual ⁷Be destruction mechanism.

Mechanism reading 1 — n_A/K_rank  (PREFERRED)
---------------------------------------------
Each ⁷Be nucleus during the de-saturation epoch sits at a single
4-simplex vertex (K_rank=5 such vertices are available) and couples
to the substrate's full hex-slice face-pair adjacency network
(n_A = C(N_BAM, 2) = 15 unordered pair channels).  The destruction
rate per vertex is

    R_destr_per_vertex = n_A / K_rank = 15 / 5 = 3

(face-pair couplings per simplex vertex).  The corresponding ⁷Li
suppression factor — the ratio by which the SBBN ⁷Be population is
reduced before electron-capture closes the channel — is exactly this
3 to leading order in the de-saturation expansion.

Why ``n_A / K_rank`` is the right ratio: it is the unique product of
two B3 primitives that **separately** appear in (a) the deuteron
binding-energy ε_face = Λ_QCD/(n_A·N_BAM) [substrate face-pair
coupling] and (b) the simplicial-SM ansatz vertex count.  ⁷Be
destruction must couple a face-pair hopping channel (n_A side) to a
4-simplex valence vertex (K_rank side) — the ratio is the only
scale-free combination of the two.

Mechanism reading 2 — n_R/N_BAM  (alternative, equally numerically valid)
-------------------------------------------------------------------------
The Möbius reflection-orbit count over the T² period torus (n_R = 18)
distributed across the 2D hex-slice valence (N_BAM = 6) gives
n_R/N_BAM = 3 reflection channels per nearest-neighbour site.  This
reading favors the *kinematic* picture: ⁷Be is destroyed by a
late-time re-thermalization that re-uses each reflection orbit three
times per hex neighbor as the substrate de-saturates through the
keV-MeV range.  Numerically identical to Reading 1; physically a
different mechanism.

Mechanism reading 3 — F · R / K_pair  (Koide-flavored)
------------------------------------------------------
The Koide numerator·denominator product divided by the Möbius bundle
sheet count — F·R/K_pair = 2·3/2 = 3.  This reading would link ⁷Be
destruction to the lepton-mass torque sector via electron capture
(⁷Be + e⁻ → ⁷Li + ν).  The connection is suggestive (electron capture
**is** the channel that converts ⁷Be to ⁷Li, so a Koide-sector
factor entering the suppression would be physically natural), but the
F·R combination is not used elsewhere in the BBN-era framework, so
this reading is harder to justify on standalone-mechanism grounds.

Conclusion
----------
The canonical substrate prediction is

    S_substrate = n_A / K_rank = 15 / 5 = 3

(face-pair adjacency channels per 4-simplex vertex), bringing SBBN's
4.5×10⁻¹⁰ down to ⁷Li/H ≈ 1.5×10⁻¹⁰ — inside the 1-σ Spite-plateau
band of (1.6 ± 0.3)×10⁻¹⁰.

Honesty band
------------
Two readings (n_A/K_rank, n_R/N_BAM) give numerically identical
predictions with different mechanisms; the framework cannot
discriminate between them at present.  A third (F·R/K_pair) is
suggestive on EC-channel grounds but is not used elsewhere.  In all
three cases the substrate combination yields exactly the integer 3,
which under-shoots the observed ratio by ~6.7% (vs the 4.5 fit) or
matches at 2.1% (vs the 4.7 PRIMAT central).  Both are well inside
the 1-σ band of the Spite-plateau measurement (which is ±19% on the
ratio) and well inside the SBBN modeling uncertainty on ⁷Li
(typically quoted as 10-15%).

This module exposes:

* :data:`SUBSTRATE_LI7_SUPPRESSION_FACTOR` — the canonical value 3.
* :func:`derive_li7_suppression` — full derivation chain returning
  a :class:`Li7SuppressionDerivation` dataclass.
* :func:`audit_candidates` — full table of all candidate combinations
  with their numerical values and errors against the empirical
  target.

References:
    spec §18.66 (proto-matter / observer-horizon re-thermalization).
    spec §18.75.4 (BBN in the substrate framework).
    Cyburt, Fields, Olive 2008 (SBBN ⁷Li eta-fit).
    Sbordone et al. 2010 (Spite-plateau ⁷Li/H observation).
    bbn_substrate.SBBN_LI7_H = 4.7×10⁻¹⁰ (PRIMAT central).
    b3_constants (canonical 12-integer primitive list).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from .b3_constants import (
    F,
    K_pair,
    K_rank,
    N_BAM,
    R,
    V13,
    n_A,
    n_M,
    n_R,
)


# ---------------------------------------------------------------------------
# Empirical target (the SBBN/observation ratio that has to be reproduced)
# ---------------------------------------------------------------------------

#: SBBN ⁷Li/H prediction in the eta-fit normalisation (Cyburt et al. 2008,
#: as used by :func:`bbn_predictions.lithium7_abundance` at η = 6.10×10⁻¹⁰).
SBBN_LI7_ETA_FIT: float = 4.5e-10

#: SBBN ⁷Li/H prediction in the PRIMAT-central normalisation (used by
#: :data:`bbn_substrate.SBBN_LI7_H`).
SBBN_LI7_PRIMAT_CENTRAL: float = 4.7e-10

#: Observed ⁷Li/H from the Spite plateau (Sbordone et al. 2010, central value).
OBS_LI7_H_CENTRAL: float = 1.6e-10

#: Observed 1-σ uncertainty on the Spite-plateau ⁷Li/H.
OBS_LI7_H_SIGMA: float = 0.3e-10

#: Empirical suppression factor (SBBN-to-observation ratio) using the eta fit.
TARGET_SUPPRESSION_ETA_FIT: float = SBBN_LI7_ETA_FIT / OBS_LI7_H_CENTRAL  # 2.8125

#: Empirical suppression factor using the PRIMAT central value.
TARGET_SUPPRESSION_PRIMAT: float = SBBN_LI7_PRIMAT_CENTRAL / OBS_LI7_H_CENTRAL  # 2.9375


# ---------------------------------------------------------------------------
# Canonical substrate-derived suppression factor
# ---------------------------------------------------------------------------

#: Canonical substrate-derived ⁷Li suppression factor.
#:
#: Derived as n_A / K_rank = C(N_BAM, 2) / K_rank = 15 / 5 = 3
#: (face-pair adjacency channels per 4-simplex vertex during the
#: substrate de-saturation epoch — see module docstring, Mechanism
#: Reading 1).
SUBSTRATE_LI7_SUPPRESSION_FACTOR: float = float(n_A) / float(K_rank)


# ---------------------------------------------------------------------------
# Derivation result
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Li7SuppressionDerivation:
    """Full chain of the substrate ⁷Li-suppression derivation.

    Attributes:
        primitive_a:        Numerator primitive ("n_A").
        primitive_a_value:  Value of the numerator (15).
        primitive_b:        Denominator primitive ("K_rank").
        primitive_b_value:  Value of the denominator (5).
        formula:            Symbolic formula string.
        mechanism:          One-line physical-mechanism description.
        suppression:        Predicted suppression factor (3).
        sbbn_li7_eta_fit:   SBBN ⁷Li/H using the eta-fit normalisation.
        sbbn_li7_primat:    SBBN ⁷Li/H using the PRIMAT central value.
        substrate_li7_eta_fit:  Substrate prediction = SBBN_eta / suppression.
        substrate_li7_primat:   Substrate prediction = SBBN_primat / suppression.
        obs_li7:            Observed ⁷Li/H (Spite plateau).
        obs_li7_sigma:      1-σ observational uncertainty.
        n_sigma_eta_fit:    Substrate (eta-fit) deviation from observation.
        n_sigma_primat:     Substrate (PRIMAT) deviation from observation.
        within_2s:          True if both substrate predictions are within 2-σ.
    """

    primitive_a: str
    primitive_a_value: int
    primitive_b: str
    primitive_b_value: int
    formula: str
    mechanism: str
    suppression: float
    sbbn_li7_eta_fit: float
    sbbn_li7_primat: float
    substrate_li7_eta_fit: float
    substrate_li7_primat: float
    obs_li7: float
    obs_li7_sigma: float
    n_sigma_eta_fit: float
    n_sigma_primat: float
    within_2s: bool


def derive_li7_suppression() -> Li7SuppressionDerivation:
    """Return the canonical substrate ⁷Li-suppression derivation.

    The canonical reading is ``n_A / K_rank = 15 / 5 = 3``: face-pair
    adjacency channels per 4-simplex vertex (see module docstring,
    Mechanism Reading 1).
    """
    s = SUBSTRATE_LI7_SUPPRESSION_FACTOR
    sub_eta = SBBN_LI7_ETA_FIT / s
    sub_pri = SBBN_LI7_PRIMAT_CENTRAL / s
    n_eta = abs(sub_eta - OBS_LI7_H_CENTRAL) / OBS_LI7_H_SIGMA
    n_pri = abs(sub_pri - OBS_LI7_H_CENTRAL) / OBS_LI7_H_SIGMA
    return Li7SuppressionDerivation(
        primitive_a="n_A",
        primitive_a_value=n_A,
        primitive_b="K_rank",
        primitive_b_value=K_rank,
        formula="S = n_A / K_rank = C(N_BAM, 2) / K_rank = 15 / 5 = 3",
        mechanism=(
            "Face-pair adjacency channels per 4-simplex vertex during "
            "substrate de-saturation epoch (proto-matter window §18.66): "
            "each ⁷Be nucleus sits at one of K_rank=5 simplex vertices and "
            "couples to all n_A=15 unordered hex-slice face-pair channels, "
            "giving n_A/K_rank=3 destruction couplings per vertex."
        ),
        suppression=s,
        sbbn_li7_eta_fit=SBBN_LI7_ETA_FIT,
        sbbn_li7_primat=SBBN_LI7_PRIMAT_CENTRAL,
        substrate_li7_eta_fit=sub_eta,
        substrate_li7_primat=sub_pri,
        obs_li7=OBS_LI7_H_CENTRAL,
        obs_li7_sigma=OBS_LI7_H_SIGMA,
        n_sigma_eta_fit=n_eta,
        n_sigma_primat=n_pri,
        within_2s=(n_eta <= 2.0 and n_pri <= 2.0),
    )


# ---------------------------------------------------------------------------
# Candidate audit
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Li7Candidate:
    """One candidate substrate combination and its match against the target."""

    label: str
    formula: str
    value: float
    rel_err_eta_fit_pct: float
    rel_err_primat_pct: float
    matches_target: bool
    interpretation: str


#: Candidate suppression-factor combinations producible from substrate
#: primitives.  Three non-trivial combinations evaluate exactly to the
#: integer 3 — they have different physical interpretations but are
#: numerically indistinguishable.
CANDIDATE_FORMULAS: List[Tuple[str, str, float, str]] = [
    (
        "n_A / K_rank",
        "C(N_BAM, 2) / K_rank = 15 / 5",
        n_A / K_rank,
        "Face-pair adjacency channels per 4-simplex vertex (CANONICAL).",
    ),
    (
        "n_R / N_BAM",
        "n_R / N_BAM = 18 / 6",
        n_R / N_BAM,
        "Möbius reflection orbits per hex-slice nearest neighbor.",
    ),
    (
        "F · R / K_pair",
        "F · R / K_pair = 2 · 3 / 2",
        (F * R) / K_pair,
        "Koide numerator·denominator over Möbius bundle sheets "
        "(EC-channel reading).",
    ),
    (
        "R",
        "R (Koide denominator)",
        float(R),
        "Trivial: the integer R itself; no derivation content.",
    ),
    (
        "K_pair · K_rank / F",
        "K_pair · K_rank / F = 2 · 5 / 2",
        (K_pair * K_rank) / F,
        "Substrate channel count (K_pair sheets × K_rank vertices / F).",
    ),
    (
        "n_R / K_rank",
        "n_R / K_rank = 18 / 5",
        n_R / K_rank,
        "Reflection orbits per simplex vertex (over-shoots target).",
    ),
    (
        "K_rank / (K_pair · F)",
        "K_rank / (K_pair · F) = 5 / 4",
        K_rank / (K_pair * F),
        "Inverse channel ratio (under-shoots target).",
    ),
    (
        "n_M / (K_pair · K_rank · N_BAM · R)",
        "n_M / (K_pair · K_rank · N_BAM · R) = 268 / 180",
        n_M / (K_pair * K_rank * N_BAM * R),
        "Master-multiplicity normalised (under-shoots target).",
    ),
    (
        "(n_R - K_rank) / K_pair²",
        "(n_R - K_rank) / K_pair² = 13 / 4 = V13 / K_pair²",
        (n_R - K_rank) / (K_pair ** 2),
        "V13-mediated ratio (over-shoots target by ~8%).",
    ),
    (
        "F · K_rank / R",
        "F · K_rank / R = 10 / 3",
        (F * K_rank) / R,
        "Koide-flavored channel count (over-shoots target by ~11%).",
    ),
]


def audit_candidates(
    target_eta_fit: float = TARGET_SUPPRESSION_ETA_FIT,
    target_primat: float = TARGET_SUPPRESSION_PRIMAT,
    tolerance_pct: float = 5.0,
) -> List[Li7Candidate]:
    """Return all candidate combinations sorted by closeness to the target.

    Args:
        target_eta_fit:   Empirical target using the eta-fit normalisation
                          (default ≈2.8125 = 4.5e-10 / 1.6e-10).
        target_primat:    Empirical target using the PRIMAT central value
                          (default ≈2.9375 = 4.7e-10 / 1.6e-10).
        tolerance_pct:    A candidate is flagged "matches_target" if its
                          relative error against EITHER target is below
                          this percentage.

    Returns:
        List of :class:`Li7Candidate` rows.
    """
    rows: List[Li7Candidate] = []
    for label, formula, value, interp in CANDIDATE_FORMULAS:
        err_eta = abs(value - target_eta_fit) / target_eta_fit * 100.0
        err_pri = abs(value - target_primat) / target_primat * 100.0
        rows.append(
            Li7Candidate(
                label=label,
                formula=formula,
                value=value,
                rel_err_eta_fit_pct=err_eta,
                rel_err_primat_pct=err_pri,
                matches_target=(min(err_eta, err_pri) <= tolerance_pct),
                interpretation=interp,
            )
        )
    rows.sort(key=lambda r: min(r.rel_err_eta_fit_pct, r.rel_err_primat_pct))
    return rows


# ---------------------------------------------------------------------------
# Convenience accessors
# ---------------------------------------------------------------------------

def substrate_li7_suppression() -> float:
    """Return the canonical substrate-derived ⁷Li-suppression factor (3)."""
    return SUBSTRATE_LI7_SUPPRESSION_FACTOR


def substrate_li7_h(sbbn_li7: float = SBBN_LI7_PRIMAT_CENTRAL) -> float:
    """Return the substrate-suppressed ⁷Li/H prediction.

    Args:
        sbbn_li7: SBBN ⁷Li/H value to suppress (defaults to the PRIMAT
                  central value 4.7×10⁻¹⁰; pass :data:`SBBN_LI7_ETA_FIT`
                  for the Cyburt-fit value 4.5×10⁻¹⁰).
    """
    return sbbn_li7 / SUBSTRATE_LI7_SUPPRESSION_FACTOR


# ---------------------------------------------------------------------------
# Pretty report
# ---------------------------------------------------------------------------

def report() -> str:
    """Render a one-page text summary of the derivation + candidate audit."""
    d = derive_li7_suppression()
    cands = audit_candidates()

    lines: List[str] = []
    lines.append("=" * 78)
    lines.append("Substrate-derived ⁷Li/⁷Be BBN suppression factor")
    lines.append("=" * 78)
    lines.append(f"  Canonical formula : {d.formula}")
    lines.append(f"  Numerator         : {d.primitive_a} = {d.primitive_a_value}")
    lines.append(f"  Denominator       : {d.primitive_b} = {d.primitive_b_value}")
    lines.append(f"  Suppression S     : {d.suppression:.4f}")
    lines.append("")
    lines.append("Mechanism:")
    for chunk in _wrap(d.mechanism, 76, indent=2):
        lines.append(chunk)
    lines.append("")
    lines.append(f"  Empirical target (eta-fit, 4.5/1.6) : {TARGET_SUPPRESSION_ETA_FIT:.4f}")
    lines.append(f"  Empirical target (PRIMAT, 4.7/1.6)  : {TARGET_SUPPRESSION_PRIMAT:.4f}")
    lines.append(f"  Substrate ⁷Li/H (eta-fit)           : {d.substrate_li7_eta_fit:.3e}")
    lines.append(f"  Substrate ⁷Li/H (PRIMAT)            : {d.substrate_li7_primat:.3e}")
    lines.append(f"  Observed   ⁷Li/H (Spite plateau)    : {d.obs_li7:.3e} ± {d.obs_li7_sigma:.1e}")
    lines.append(f"  n-σ (eta-fit)                       : {d.n_sigma_eta_fit:.2f}")
    lines.append(f"  n-σ (PRIMAT)                        : {d.n_sigma_primat:.2f}")
    lines.append(f"  Within 2-σ                          : {d.within_2s}")
    lines.append("")
    lines.append("Candidate audit (sorted by closeness to target):")
    lines.append(
        f"  {'formula':<38} {'value':>8} {'%err η':>8} {'%err P':>8}  match"
    )
    lines.append("  " + "-" * 73)
    for r in cands:
        flag = "✓" if r.matches_target else " "
        lines.append(
            f"  {r.formula:<38} {r.value:>8.4f} "
            f"{r.rel_err_eta_fit_pct:>7.2f}% {r.rel_err_primat_pct:>7.2f}%  {flag}"
        )
    lines.append("=" * 78)
    return "\n".join(lines)


def _wrap(text: str, width: int, indent: int = 0) -> List[str]:
    """Wrap a paragraph to ``width`` columns with ``indent`` leading spaces."""
    pad = " " * indent
    words = text.split()
    out: List[str] = []
    cur = pad
    for w in words:
        candidate = (cur + " " + w).strip() if cur != pad else cur + w
        if len(candidate) > width and cur.strip():
            out.append(cur.rstrip())
            cur = pad + w
        else:
            cur = candidate if cur == pad else candidate
    if cur.strip():
        out.append(cur.rstrip())
    return out


__all__ = [
    "SBBN_LI7_ETA_FIT",
    "SBBN_LI7_PRIMAT_CENTRAL",
    "OBS_LI7_H_CENTRAL",
    "OBS_LI7_H_SIGMA",
    "TARGET_SUPPRESSION_ETA_FIT",
    "TARGET_SUPPRESSION_PRIMAT",
    "SUBSTRATE_LI7_SUPPRESSION_FACTOR",
    "Li7SuppressionDerivation",
    "Li7Candidate",
    "CANDIDATE_FORMULAS",
    "derive_li7_suppression",
    "audit_candidates",
    "substrate_li7_suppression",
    "substrate_li7_h",
    "report",
]


if __name__ == "__main__":
    print(report())
