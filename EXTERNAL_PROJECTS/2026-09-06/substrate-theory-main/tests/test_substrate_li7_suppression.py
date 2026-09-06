"""Tests for stiff_medium.substrate_li7_suppression."""

from __future__ import annotations

import math

import pytest

from stiff_medium.b3_constants import K_pair, K_rank, N_BAM, R, F, n_A, n_R
from stiff_medium.substrate_li7_suppression import (
    CANDIDATE_FORMULAS,
    Li7Candidate,
    Li7SuppressionDerivation,
    OBS_LI7_H_CENTRAL,
    OBS_LI7_H_SIGMA,
    SBBN_LI7_ETA_FIT,
    SBBN_LI7_PRIMAT_CENTRAL,
    SUBSTRATE_LI7_SUPPRESSION_FACTOR,
    TARGET_SUPPRESSION_ETA_FIT,
    TARGET_SUPPRESSION_PRIMAT,
    audit_candidates,
    derive_li7_suppression,
    report,
    substrate_li7_h,
    substrate_li7_suppression,
)


# ---------------------------------------------------------------------------
# Canonical constant
# ---------------------------------------------------------------------------


def test_canonical_suppression_is_3() -> None:
    """The substrate-derived suppression must equal 3 exactly."""
    assert math.isclose(SUBSTRATE_LI7_SUPPRESSION_FACTOR, 3.0, rel_tol=1e-12)


def test_canonical_suppression_equals_n_A_over_K_rank() -> None:
    """The canonical formula is n_A / K_rank = 15 / 5 = 3."""
    assert math.isclose(SUBSTRATE_LI7_SUPPRESSION_FACTOR,
                        n_A / K_rank, rel_tol=1e-12)


def test_substrate_li7_suppression_helper() -> None:
    assert substrate_li7_suppression() == SUBSTRATE_LI7_SUPPRESSION_FACTOR


# ---------------------------------------------------------------------------
# Empirical targets
# ---------------------------------------------------------------------------


def test_target_eta_fit() -> None:
    """4.5e-10 / 1.6e-10 = 2.8125."""
    assert math.isclose(TARGET_SUPPRESSION_ETA_FIT, 2.8125, rel_tol=1e-9)


def test_target_primat() -> None:
    """4.7e-10 / 1.6e-10 = 2.9375."""
    assert math.isclose(TARGET_SUPPRESSION_PRIMAT, 2.9375, rel_tol=1e-9)


def test_substrate_within_observational_one_sigma() -> None:
    """Substrate predicts ⁷Li/H = 1.5 to 1.567 ×10⁻¹⁰; observed 1.6 ± 0.3."""
    pred_eta = SBBN_LI7_ETA_FIT / SUBSTRATE_LI7_SUPPRESSION_FACTOR
    pred_pri = SBBN_LI7_PRIMAT_CENTRAL / SUBSTRATE_LI7_SUPPRESSION_FACTOR
    # Both must be within 1-σ of observation.
    assert abs(pred_eta - OBS_LI7_H_CENTRAL) <= OBS_LI7_H_SIGMA
    assert abs(pred_pri - OBS_LI7_H_CENTRAL) <= OBS_LI7_H_SIGMA


# ---------------------------------------------------------------------------
# Derivation dataclass
# ---------------------------------------------------------------------------


def test_derive_returns_dataclass() -> None:
    d = derive_li7_suppression()
    assert isinstance(d, Li7SuppressionDerivation)


def test_derivation_uses_canonical_primitives() -> None:
    d = derive_li7_suppression()
    assert d.primitive_a == "n_A"
    assert d.primitive_a_value == n_A == 15
    assert d.primitive_b == "K_rank"
    assert d.primitive_b_value == K_rank == 5


def test_derivation_suppression_value() -> None:
    d = derive_li7_suppression()
    assert math.isclose(d.suppression, 3.0, rel_tol=1e-12)


def test_derivation_predictions_match_arithmetic() -> None:
    d = derive_li7_suppression()
    assert math.isclose(d.substrate_li7_eta_fit,
                        SBBN_LI7_ETA_FIT / 3.0, rel_tol=1e-12)
    assert math.isclose(d.substrate_li7_primat,
                        SBBN_LI7_PRIMAT_CENTRAL / 3.0, rel_tol=1e-12)


def test_derivation_within_2_sigma() -> None:
    """Substrate prediction must land within 2-σ of the Spite plateau."""
    d = derive_li7_suppression()
    assert d.within_2s is True
    assert d.n_sigma_eta_fit <= 2.0
    assert d.n_sigma_primat <= 2.0


def test_derivation_mechanism_text_present() -> None:
    d = derive_li7_suppression()
    assert "face-pair" in d.mechanism.lower() or "adjacency" in d.mechanism.lower()
    assert "vertex" in d.mechanism.lower() or "simplex" in d.mechanism.lower()


def test_derivation_formula_string() -> None:
    d = derive_li7_suppression()
    assert "n_A" in d.formula
    assert "K_rank" in d.formula
    assert "15" in d.formula
    assert "5" in d.formula


# ---------------------------------------------------------------------------
# Candidate audit
# ---------------------------------------------------------------------------


def test_audit_returns_nonempty_list() -> None:
    rows = audit_candidates()
    assert len(rows) >= len(CANDIDATE_FORMULAS)
    for r in rows:
        assert isinstance(r, Li7Candidate)


def test_audit_sorted_by_closeness() -> None:
    rows = audit_candidates()
    for i in range(len(rows) - 1):
        a = min(rows[i].rel_err_eta_fit_pct, rows[i].rel_err_primat_pct)
        b = min(rows[i + 1].rel_err_eta_fit_pct, rows[i + 1].rel_err_primat_pct)
        assert a <= b


def test_audit_contains_canonical_combination() -> None:
    rows = audit_candidates()
    labels = [r.label for r in rows]
    assert "n_A / K_rank" in labels


def test_audit_canonical_has_lowest_or_tied_error() -> None:
    """The canonical combination must be among the best matches."""
    rows = audit_candidates()
    canonical = next(r for r in rows if r.label == "n_A / K_rank")
    best_err = min(min(r.rel_err_eta_fit_pct, r.rel_err_primat_pct) for r in rows)
    canon_err = min(canonical.rel_err_eta_fit_pct, canonical.rel_err_primat_pct)
    assert math.isclose(canon_err, best_err, rel_tol=1e-9)


def test_audit_three_combinations_all_evaluate_to_3() -> None:
    """n_A/K_rank, n_R/N_BAM, F·R/K_pair all give exactly 3."""
    rows = audit_candidates()
    threes = [r for r in rows if math.isclose(r.value, 3.0, rel_tol=1e-12)]
    labels = {r.label for r in threes}
    assert "n_A / K_rank" in labels
    assert "n_R / N_BAM" in labels
    assert "F · R / K_pair" in labels


def test_audit_canonical_matches_target_within_tolerance() -> None:
    rows = audit_candidates(tolerance_pct=10.0)
    canon = next(r for r in rows if r.label == "n_A / K_rank")
    assert canon.matches_target is True


# ---------------------------------------------------------------------------
# substrate_li7_h helper
# ---------------------------------------------------------------------------


def test_substrate_li7_h_default_uses_primat() -> None:
    val = substrate_li7_h()
    expected = SBBN_LI7_PRIMAT_CENTRAL / SUBSTRATE_LI7_SUPPRESSION_FACTOR
    assert math.isclose(val, expected, rel_tol=1e-12)


def test_substrate_li7_h_eta_fit_normalisation() -> None:
    val = substrate_li7_h(SBBN_LI7_ETA_FIT)
    expected = SBBN_LI7_ETA_FIT / SUBSTRATE_LI7_SUPPRESSION_FACTOR
    assert math.isclose(val, expected, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def test_report_renders() -> None:
    text = report()
    assert len(text) > 200
    for kw in ("Substrate", "n_A", "K_rank", "Suppression", "Spite", "Mechanism", "audit"):
        # The report uses different cases — be lenient.
        assert kw.lower() in text.lower()


# ---------------------------------------------------------------------------
# Cross-check: integers tied to b3_constants
# ---------------------------------------------------------------------------


def test_n_A_canonical_is_15() -> None:
    assert n_A == 15


def test_K_rank_canonical_is_5() -> None:
    assert K_rank == 5


def test_R_canonical_is_3() -> None:
    """R = 3 alone gives 3 trivially — this is why we need a non-trivial derivation."""
    assert R == 3
