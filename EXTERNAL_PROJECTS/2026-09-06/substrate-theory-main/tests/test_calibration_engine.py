"""Tests for stiff_medium.calibration_engine."""

from __future__ import annotations

import os
import sys

import pytest

# Make src/ importable when running pytest from repo root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from stiff_medium.calibration_engine import (  # noqa: E402
    CalibrationEngine,
    Classification,
    Prediction,
    build_default_b3_engine,
)


# ---------------------------------------------------------------------------
# Unit-level tests on the small primitives
# ---------------------------------------------------------------------------
def test_prediction_residual_and_chi_squared():
    p = Prediction(
        name="x",
        value=1.05,
        observed=1.00,
        classification=Classification.FORCED,
        sigma=0.05,
    )
    assert p.residual == pytest.approx(0.05)
    assert p.relative_error == pytest.approx(0.05)
    assert p.chi_squared == pytest.approx(1.0)


def test_register_prediction_accepts_string_class():
    eng = CalibrationEngine()
    eng.register_prediction("a", 1.0, 1.0, "forced", sigma=0.01)
    eng.register_prediction("b", 2.0, 2.0, "DERIVED", sigma=0.02)
    assert len(eng.forced_predictions()) == 1
    assert len(eng.derived_predictions()) == 1


def test_register_prediction_overwrites_duplicate_names():
    eng = CalibrationEngine()
    eng.register_prediction("x", 1.0, 1.0, Classification.FORCED, sigma=0.1)
    eng.register_prediction("x", 2.0, 2.0, Classification.DERIVED, sigma=0.1)
    assert len(eng.predictions) == 1
    assert eng.predictions[0].classification == Classification.DERIVED


def test_dof_subtracts_calibrated_count():
    eng = CalibrationEngine()
    for i in range(5):
        eng.register_prediction(
            f"f{i}", 1.0, 1.0, Classification.FORCED, sigma=0.1
        )
    for i in range(2):
        eng.register_prediction(
            f"c{i}", 1.0, 1.0, Classification.CALIBRATED, sigma=0.1
        )
    # 7 obs - 2 calibrated = 5 dof
    assert eng.degrees_of_freedom() == 5


def test_aic_bic_scoring_basic():
    import math

    eng = CalibrationEngine()
    eng.register_prediction("x", 1.0, 1.0, Classification.FORCED, sigma=0.1)
    eng.register_prediction("y", 1.0, 1.0, Classification.CALIBRATED, sigma=0.1)
    ic = eng.aic_bic_scores()
    assert ic["k"] == 1
    assert ic["n"] == 2
    assert ic["AIC"] == pytest.approx(ic["chi2"] + 2.0)
    assert ic["BIC"] == pytest.approx(ic["chi2"] + math.log(2))


def test_compare_to_SM_reports_reduction():
    eng = CalibrationEngine(sm_parameter_count=25)
    for i in range(5):
        eng.register_prediction(
            f"c{i}", 1.0, 1.0, Classification.CALIBRATED, sigma=0.1
        )
    cmp = eng.compare_to_SM()
    assert cmp["SM_params"] == 25
    assert cmp["framework_calibrated_params"] == 5
    assert cmp["reduction_factor"] == pytest.approx(5.0)
    assert cmp["delta"] == 20


def test_report_runs_and_contains_section_headers():
    eng = build_default_b3_engine()
    txt = eng.report()
    assert "CALIBRATION ENGINE REPORT" in txt
    assert "FORCED" in txt
    assert "CALIBRATED" in txt
    assert "reduced chi^2" in txt


# ---------------------------------------------------------------------------
# The three required acceptance tests
# ---------------------------------------------------------------------------
def test_a_forced_predictions_at_least_14():
    eng = build_default_b3_engine()
    n_forced = len(eng.forced_predictions())
    assert n_forced >= 14, (
        f"Expected >= 14 FORCED predictions, got {n_forced}"
    )


def test_b_calibrated_parameters_at_most_6():
    eng = build_default_b3_engine()
    n_calib = len(eng.calibrated_parameters())
    assert n_calib <= 6, (
        f"Expected <= 6 CALIBRATED parameters, got {n_calib}"
    )


def test_c_reduced_chi_squared_forced_below_2():
    eng = build_default_b3_engine()
    red = eng.reduced_chi_squared_forced()
    assert red < 2.0, (
        f"Expected reduced chi^2 < 2 for FORCED set, got {red:.3f}"
    )
