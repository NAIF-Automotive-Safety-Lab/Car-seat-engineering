"""Tests for the topological-defect Poisson perturbation model."""

import math

import pytest

from src.stiff_medium.topological_defect_perturbations import (
    LAMBDA_QCD_GEV,
    N_M_DEFAULT,
    OBS_DELTA_RHO_OVER_RHO,
    TopologicalDefectPerturbations,
    summary,
)


@pytest.fixture
def tdp() -> TopologicalDefectPerturbations:
    return TopologicalDefectPerturbations()


def test_defaults(tdp):
    assert tdp.n_M == 268
    assert tdp.Lambda_QCD_GeV == pytest.approx(LAMBDA_QCD_GEV)
    assert tdp.observed_amplitude == OBS_DELTA_RHO_OVER_RHO


def test_freeze_out_amplitude_about_six_percent(tdp):
    A = tdp.freeze_out_amplitude()
    # 1/sqrt(268) = 0.06108...
    assert A == pytest.approx(1.0 / math.sqrt(268), rel=1e-12)
    assert 0.05 < A < 0.07


def test_freeze_out_independent_of_lambda_qcd():
    a = TopologicalDefectPerturbations(Lambda_QCD_GeV=0.2).freeze_out_amplitude()
    b = TopologicalDefectPerturbations(Lambda_QCD_GeV=0.5).freeze_out_amplitude()
    assert a == pytest.approx(b, rel=1e-12)


def test_suppression_factor_basic(tdp):
    assert tdp.suppression_factor(0.0) == pytest.approx(1.0)
    assert tdp.suppression_factor(1.0) == pytest.approx(math.exp(-1.0))
    assert tdp.suppression_factor(10.0) < tdp.suppression_factor(5.0)


def test_observable_amplitude_consistency(tdp):
    N = 5.0
    expected = tdp.freeze_out_amplitude() * math.exp(-N)
    assert tdp.observable_amplitude(N) == pytest.approx(expected, rel=1e-12)


def test_predicted_n_efolds_for_observed_around_8_to_9(tdp):
    """Required N* = ln(0.0611 / 1e-5) ~ 8.72."""
    N = tdp.predicted_n_efolds_for_observed()
    assert 8.0 < N < 9.0
    assert N == pytest.approx(math.log(tdp.freeze_out_amplitude() / 1e-5), rel=1e-12)


def test_eight_to_nine_efolds_reproduces_observed_amplitude(tdp):
    """Within an order of magnitude of 1e-5 for N in [8, 9]."""
    A_low = tdp.observable_amplitude(8.0)
    A_high = tdp.observable_amplitude(9.0)
    assert A_high < OBS_DELTA_RHO_OVER_RHO < A_low
    # both should be within an order of magnitude of observed
    assert 0.1 * OBS_DELTA_RHO_OVER_RHO < A_low < 10 * OBS_DELTA_RHO_OVER_RHO
    assert 0.1 * OBS_DELTA_RHO_OVER_RHO < A_high < 10 * OBS_DELTA_RHO_OVER_RHO


def test_round_trip_predicted_efolds(tdp):
    N_star = tdp.predicted_n_efolds_for_observed()
    assert tdp.observable_amplitude(N_star) == pytest.approx(
        OBS_DELTA_RHO_OVER_RHO, rel=1e-10
    )


def test_n_s_near_unity(tdp):
    assert tdp.predicted_n_s() == pytest.approx(1.0)
    # Planck observed
    assert abs(tdp.predicted_n_s() - tdp.n_s_observed()) < 0.05


def test_freeze_out_hubble_radius_femtometer_scale(tdp):
    R = tdp.freeze_out_hubble_radius_m()
    # hbar c / 0.218 GeV ~ 9e-16 m ~ 0.9 fm
    assert 1e-16 < R < 1e-14


def test_lambda_qcd_feature_scale_physics(tdp):
    """R_freeze ~ 0.9 fm; cosmic redshift z~9e11 -> sub-mm scale today.

    This is a HONEST result: the naive feature scale is far below galaxy-survey
    resolution, so the simple chain does NOT predict an Mpc bump. A
    Mpc-scale prediction would require additional defect coalescence /
    growth, not modeled here.
    """
    feat = tdp.Lambda_QCD_feature_scale()
    assert feat["R_today_m"] > feat["R_freeze_m"]
    # freeze-out radius is sub-fm to fm
    assert 1e-16 < feat["R_freeze_m"] < 1e-14
    # today: stretched by ~9e11 -> sub-mm
    assert 1e-5 < feat["R_today_m"] < 1e-1, f"R_today_m={feat['R_today_m']}"
    # In Mpc, this is many orders of magnitude below 1 Mpc
    assert feat["R_today_Mpc"] < 1e-20


def test_feature_scale_dominated_by_cosmic_expansion(tdp):
    """Cosmic redshift z_QCD ~ 9e11 dominates; e-folds are an amplitude knob."""
    a = tdp.Lambda_QCD_feature_scale(N_efolds=5.0)
    b = tdp.Lambda_QCD_feature_scale(N_efolds=10.0)
    # Scale today is set by cosmic redshift, not by amplitude e-folds
    assert a["R_today_Mpc"] == pytest.approx(b["R_today_Mpc"], rel=1e-12)
    # But amplitude stretch differs
    assert b["amplitude_stretch_factor"] > a["amplitude_stretch_factor"]
    # Cosmic redshift in expected ballpark
    assert 1e11 < a["z_QCD_to_today"] < 1e13


def test_compare_to_inflation_keys(tdp):
    cmp = tdp.compare_to_inflation()
    for key in [
        "amplitude_origin",
        "amplitude_value",
        "n_s_value",
        "tensor_to_scalar",
        "feature_at_Lambda_QCD",
        "free_parameters",
    ]:
        assert key in cmp
        assert "inflation" in cmp[key]
        assert "substrate" in cmp[key]


def test_gap_status_open_partial(tdp):
    g = tdp.gap_status()
    assert g["status"] == "open-partial"
    assert g["N_efolds_derived_from_substrate"] is None
    assert g["freeze_out_amplitude"] == pytest.approx(1.0 / math.sqrt(268), rel=1e-12)
    assert 8.0 < g["N_efolds_required"] < 9.0
    assert isinstance(g["what_works"], list) and len(g["what_works"]) >= 3
    assert isinstance(g["what_is_missing"], list) and len(g["what_is_missing"]) >= 2
    assert "falsifier" in g


def test_summary_smoke():
    s = summary()
    assert s["freeze_out_amplitude"] > 0.05
    assert "feature_scale" in s
    assert "gap_status" in s


def test_custom_n_M_changes_amplitude():
    """Sanity: doubling n_M reduces freeze amplitude by sqrt(2)."""
    a = TopologicalDefectPerturbations(n_M=268).freeze_out_amplitude()
    b = TopologicalDefectPerturbations(n_M=536).freeze_out_amplitude()
    assert a / b == pytest.approx(math.sqrt(2.0), rel=1e-12)
