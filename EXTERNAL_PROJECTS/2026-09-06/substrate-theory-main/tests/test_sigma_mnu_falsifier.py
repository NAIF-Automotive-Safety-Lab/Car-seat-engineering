"""Tests for the B3 substrate Sigma m_nu falsifier."""

from __future__ import annotations

import pytest

from src.stiff_medium.sigma_mnu_falsifier import (
    BOUNDS_CATALOG,
    MIN_SUM_NORMAL_HIERARCHY_MEV,
    SIGMA_MNU_PREDICTION_MEV,
    SigmaMNuFalsifier,
)


# ---------------------------------------------------------------------------
# (a) Substrate prediction matches the B3-anchored 60.5 meV value
# ---------------------------------------------------------------------------

def test_substrate_prediction_anchored_value():
    f = SigmaMNuFalsifier()
    assert f.prediction_meV == pytest.approx(60.5, abs=0.01)


def test_anchor_formula_reproduces_prediction_within_5pct():
    """The 15/16 -> oscillation chain should land near 60.5 meV."""
    val = SigmaMNuFalsifier.substrate_prediction_meV()
    assert val == pytest.approx(SIGMA_MNU_PREDICTION_MEV, rel=0.05)


# ---------------------------------------------------------------------------
# (b) Currently passes the LCDM cap (DESI DR2 + ACT + Planck, <= 64.2 meV)
# ---------------------------------------------------------------------------

def test_passes_lcdm_cap():
    f = SigmaMNuFalsifier()
    bound = next(b for b in BOUNDS_CATALOG if b.upper_meV == 64.2)
    report = f.evaluate_bound(bound)
    assert report.status in ("passes", "BRINK")
    assert report.margin_meV > 0
    # Within ~10% of bound -> should be flagged BRINK
    assert report.margin_fraction < 0.10


def test_passes_planck_loose_bound():
    f = SigmaMNuFalsifier()
    bound = next(b for b in BOUNDS_CATALOG if "Planck 2018" in b.name)
    report = f.evaluate_bound(bound)
    assert report.status == "passes"
    assert report.margin_meV > 50  # well clear


# ---------------------------------------------------------------------------
# (c) Fails strict free-streaming bound by ~14%
# ---------------------------------------------------------------------------

def test_fails_strict_free_streaming_by_about_14pct():
    f = SigmaMNuFalsifier()
    fs_bound = next(b for b in BOUNDS_CATALOG if "free-streaming" in b.name)
    report = f.evaluate_bound(fs_bound)
    assert report.status == "FAILS"
    # Prediction (60.5) exceeds bound (53) by (60.5-53)/53 = 14.15%
    overshoot = (f.prediction_meV - fs_bound.upper_meV) / fs_bound.upper_meV
    assert overshoot == pytest.approx(0.1415, abs=0.005)


# ---------------------------------------------------------------------------
# (d) DESI DR3 will be decisive
# ---------------------------------------------------------------------------

def test_desi_dr3_is_decisive():
    f = SigmaMNuFalsifier()
    decisive, rationale = f.predict_decisive_test()
    assert "DESI DR3" in decisive.name
    assert decisive.year <= 2027
    # Decision threshold must clear both the prediction and NH minimum
    threshold = max(SIGMA_MNU_PREDICTION_MEV, MIN_SUM_NORMAL_HIERARCHY_MEV)
    assert decisive.upper_meV < threshold
    assert "decision threshold" in rationale


# ---------------------------------------------------------------------------
# Auxiliary API surface
# ---------------------------------------------------------------------------

def test_falsification_threshold_equals_prediction():
    f = SigmaMNuFalsifier()
    assert f.falsification_threshold() == f.prediction_meV


def test_forecast_for_year_monotonic_improvement():
    f = SigmaMNuFalsifier()
    # 2024: best available is DESI DR1 (72 meV)
    b2024 = f.forecast_for_year(2024)
    # 2026: DR3 forecast (30 meV) becomes available
    b2026 = f.forecast_for_year(2026)
    # 2032: CMB-S4 (20 meV)
    b2032 = f.forecast_for_year(2032)
    assert b2024 is not None and b2026 is not None and b2032 is not None
    assert b2024.upper_meV >= b2026.upper_meV >= b2032.upper_meV


def test_forecast_for_year_returns_none_for_prehistory():
    f = SigmaMNuFalsifier()
    assert f.forecast_for_year(1900) is None


def test_combined_with_majorana_handles_missing_module():
    f = SigmaMNuFalsifier()
    out = f.combined_with_majorana()
    assert out["sigma_mnu_meV"] == pytest.approx(60.5, abs=0.01)
    # Module not present in repo -> consistency reports absence
    assert out["consistency"] in (
        "module absent",
        "ok",
        "module present, no m_bb function found",
    )


def test_status_report_runs_and_contains_key_sections():
    f = SigmaMNuFalsifier()
    text = f.status_report()
    assert "Substrate prediction" in text
    assert "Decisive test" in text
    assert "DESI DR2 strict free-streaming" in text
    # Print so the verifier sees current status when run with -s
    print("\n" + text)


def test_full_evaluation_reports_have_consistent_signs():
    f = SigmaMNuFalsifier()
    for r in f.evaluate_all():
        if r.margin_meV >= 0:
            assert r.status in ("passes", "BRINK")
        else:
            assert r.status == "FAILS"
