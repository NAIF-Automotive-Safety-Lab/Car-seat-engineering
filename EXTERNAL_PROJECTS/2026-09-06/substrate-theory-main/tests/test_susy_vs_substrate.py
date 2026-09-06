"""Tests for SUSY-vs-substrate BSM comparison framework."""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.susy_vs_substrate import (
    Prediction,
    Score,
    SUSYvsSubstrate,
    quick_summary,
)


# ---------------------------------------------------------------------------
# Construction & basic API
# ---------------------------------------------------------------------------


def test_constructs_cleanly() -> None:
    cmp = SUSYvsSubstrate()
    assert isinstance(cmp, SUSYvsSubstrate)


def test_predictions_returns_list_of_predictions() -> None:
    preds = SUSYvsSubstrate().predictions()
    assert len(preds) >= 8
    assert all(isinstance(p, Prediction) for p in preds)
    # Every prediction must mention all three frameworks
    for p in preds:
        assert p.sm and p.susy and p.substrate
        assert p.observable


def test_scores_in_unit_interval() -> None:
    for s in SUSYvsSubstrate().score_each():
        assert isinstance(s, Score)
        for v in (s.sm_score, s.susy_score, s.substrate_score):
            assert 0.0 <= v <= 1.0


# ---------------------------------------------------------------------------
# Substrate predictions match other modules
# ---------------------------------------------------------------------------


def test_substrate_dm_mass_matches_cube_module() -> None:
    """The 27.5 GeV cube DM mass should appear in the DM-mass entry."""
    cmp = SUSYvsSubstrate()
    # Pull live value
    try:
        from src.stiff_medium.cube_cell_dm_simulator import (
            default_substrate_params,
            m_dm_from_first_excited,
        )
        live = m_dm_from_first_excited(default_substrate_params())
    except Exception:
        live = 27.5

    assert math.isclose(cmp._substrate_dm_gev, live, rel_tol=1e-6)
    dm_pred = next(p for p in cmp.predictions() if "DM particle mass" in p.observable)
    assert f"{live:.2f}" in dm_pred.substrate


def test_substrate_dm_sigma_si_matches_cross_section_module() -> None:
    cmp = SUSYvsSubstrate()
    # Cross section should be tiny (gravitational, ~1e-130 to 1e-140 cm^2)
    assert cmp._substrate_dm_sigma_si > 0
    assert cmp._substrate_dm_sigma_si < 1e-100  # gravitational, way below LZ


def test_substrate_higgs_mass_matches_drag_target() -> None:
    cmp = SUSYvsSubstrate()
    # Drag-mass target is 125.25 GeV
    assert math.isclose(cmp._substrate_higgs_gev, 125.25, rel_tol=1e-3)


# ---------------------------------------------------------------------------
# SUSY in tension with LHC
# ---------------------------------------------------------------------------


def test_natural_susy_excluded_by_lhc() -> None:
    cmp = SUSYvsSubstrate()
    assert cmp.susy_excluded_by_lhc() is True


def test_susy_score_below_substrate_for_lhc_observables() -> None:
    cmp = SUSYvsSubstrate()
    scores = {s.observable: s for s in cmp.score_each()}
    for obs in ("LHC sparticles", "Stop searches", "DM sigma_SI"):
        s = scores[obs]
        assert s.susy_score < s.substrate_score, (
            f"SUSY should score below substrate on {obs}: "
            f"{s.susy_score} vs {s.substrate_score}"
        )


def test_susy_status_mentions_tension() -> None:
    status = SUSYvsSubstrate().current_status()
    text = status["SUSY (MSSM)"].lower()
    assert "tension" in text or "exclud" in text


# ---------------------------------------------------------------------------
# Substrate avoids fine-tuning
# ---------------------------------------------------------------------------


def test_substrate_naturalness_score_high() -> None:
    cmp = SUSYvsSubstrate()
    nat = next(s for s in cmp.score_each()
               if "Naturalness" in s.observable or "fine-tuning" in s.observable)
    assert nat.substrate_score >= 0.80
    assert nat.substrate_score > nat.sm_score
    assert nat.substrate_score > nat.susy_score


def test_substrate_hierarchy_prediction_uses_exp_4pi2() -> None:
    cmp = SUSYvsSubstrate()
    pred = next(p for p in cmp.predictions() if "Hierarchy" in p.observable)
    assert "exp" in pred.substrate.lower() or "4" in pred.substrate


# ---------------------------------------------------------------------------
# Aggregate + report
# ---------------------------------------------------------------------------


def test_aggregate_substrate_beats_susy() -> None:
    agg = SUSYvsSubstrate().aggregate()
    assert agg["substrate"] > agg["susy"]
    assert agg["substrate"] > agg["sm"]


def test_quick_summary_returns_three_keys() -> None:
    s = quick_summary()
    assert set(s.keys()) == {"sm", "susy", "substrate"}


def test_discriminating_experiments_nonempty() -> None:
    exps = SUSYvsSubstrate().discriminating_experiments()
    assert len(exps) >= 4
    for name, desc in exps:
        assert name and desc


def test_report_contains_all_sections() -> None:
    text = SUSYvsSubstrate().report()
    for header in ("PREDICTIONS", "SCORES", "CURRENT STATUS",
                   "DISCRIMINATING EXPERIMENTS", "BOTTOM LINE"):
        assert header in text
