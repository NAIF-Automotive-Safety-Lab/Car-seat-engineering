"""Tests for CMB B-mode polarization predictor."""

import math

import pytest

from src.stiff_medium.cmb_bmode_polarization import (
    BICEP_KECK_2021_LIMIT,
    CMB_S4_TARGET_SIGMA_R,
    LITEBIRD_TARGET_SIGMA_R,
    CMBBModePolarization,
)


@pytest.fixture
def model():
    return CMBBModePolarization()


# (a) substrate r = 0 exactly ------------------------------------------------
def test_substrate_r_is_zero(model):
    assert model.tensor_to_scalar_ratio_substrate() == 0.0


# (b) consistent with current BICEP/Keck upper limit -------------------------
def test_consistent_with_bicep_keck(model):
    comp = model.compare_to_bicep_keck()
    assert comp.experiment == "BICEP/Keck 2021"
    assert comp.consistent is True
    assert comp.substrate_r < BICEP_KECK_2021_LIMIT
    # Not yet decisive: BICEP/Keck cannot rule out small-r inflation
    assert comp.decisive is False


# (c) LiteBIRD and CMB-S4 will be decisive by ~2030 --------------------------
def test_litebird_is_decisive(model):
    comp = model.compare_to_litebird()
    assert comp.decisive is True
    assert comp.sensitivity_r == pytest.approx(LITEBIRD_TARGET_SIGMA_R)
    assert comp.sensitivity_r <= 1e-3
    assert comp.consistent is True


def test_cmb_s4_is_decisive(model):
    comp = model.compare_to_cmb_s4()
    assert comp.decisive is True
    assert comp.sensitivity_r == pytest.approx(CMB_S4_TARGET_SIGMA_R)
    assert comp.sensitivity_r <= 1e-3
    assert comp.consistent is True


def test_decisive_tests_below_current_limit(model):
    # Both experiments push at least an order of magnitude below BICEP/Keck
    for comp in (model.compare_to_litebird(), model.compare_to_cmb_s4()):
        assert comp.sensitivity_r * 10 <= BICEP_KECK_2021_LIMIT


# (d) lensing B-mode amplitude positive and peaks near l ~ 1000 --------------
def test_lensing_b_mode_positive_at_peak(model):
    val = model.lensing_b_mode(1000)
    assert val > 0
    assert math.isclose(val, 1.0e-4, rel_tol=1e-9)


def test_lensing_b_mode_peaks_near_1000(model):
    # value at l=1000 should exceed values far from the peak
    peak = model.lensing_b_mode(1000)
    low = model.lensing_b_mode(50)
    high = model.lensing_b_mode(3000)
    assert peak > low
    assert peak > high


def test_lensing_b_mode_zero_at_invalid_l(model):
    assert model.lensing_b_mode(0) == 0.0
    assert model.lensing_b_mode(-100) == 0.0


# Inflation B-mode behavior (sanity) -----------------------------------------
def test_inflation_b_mode_zero_for_r_zero(model):
    for l in (5, 80, 500, 1000):
        assert model.inflation_b_mode(l, 0.0) == 0.0


def test_inflation_b_mode_scales_with_r(model):
    a = model.inflation_b_mode(80, 0.01)
    b = model.inflation_b_mode(80, 0.001)
    c = model.inflation_b_mode(80, 0.0001)
    assert a > b > c > 0
    # linear in r
    assert math.isclose(a / b, 10.0, rel_tol=1e-6)
    assert math.isclose(b / c, 10.0, rel_tol=1e-6)


def test_inflation_recombination_bump_near_l_80(model):
    r = 0.01
    peak = model.inflation_b_mode(80, r)
    far = model.inflation_b_mode(400, r)
    assert peak > far


# Falsifier behavior ---------------------------------------------------------
def test_falsifier_consistent_with_no_detection(model):
    out = model.substrate_falsifier(r_observed=0.0005, r_uncertainty=0.001)
    assert out["substrate_falsified"] is False
    assert "CONSISTENT" in out["verdict"]


def test_falsifier_kills_substrate_on_clean_detection(model):
    # 5-sigma detection of r = 0.005 with sigma=0.001
    out = model.substrate_falsifier(r_observed=0.005, r_uncertainty=0.001)
    assert out["substrate_falsified"] is True
    assert out["significance_sigma"] == pytest.approx(5.0)
    assert "FALSIFIED" in out["verdict"]


# Report runs cleanly --------------------------------------------------------
def test_report_runs_and_mentions_key_facts(model):
    text = model.report()
    assert "r = 0.0000" in text
    assert "BICEP/Keck 2021" in text
    assert "LiteBIRD" in text
    assert "CMB-S4" in text
    assert "Falsifier" in text
