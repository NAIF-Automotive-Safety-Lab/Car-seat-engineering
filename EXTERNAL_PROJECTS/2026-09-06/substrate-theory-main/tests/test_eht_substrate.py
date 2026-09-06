"""Tests for src/stiff_medium/eht_substrate.py."""

from math import sqrt

import pytest

from src.stiff_medium.eht_substrate import EHTSubstrate


@pytest.fixture
def eht():
    return EHTSubstrate()


def test_M87_ring_within_5pct(eht):
    src = eht.M87_parameters()
    pred = eht.predicted_ring_diameter(src.mass_msun, src.distance_m, spin=0.0)
    err = abs(pred - 42.0) / 42.0
    assert err < 0.05, f"M87 prediction {pred:.2f} uas vs 42 uas (err={err:.2%})"


def test_SgrA_ring_within_5pct(eht):
    src = eht.SgrA_parameters()
    pred = eht.predicted_ring_diameter(src.mass_msun, src.distance_m, spin=0.0)
    err = abs(pred - 51.8) / 51.8
    assert err < 0.05, f"Sgr A* prediction {pred:.2f} uas vs 51.8 uas (err={err:.2%})"


def test_shadow_radius_is_critical_impact_parameter(eht):
    b_c = eht.shadow_radius(spin=0.0)
    assert b_c == pytest.approx(3.0 * sqrt(3.0), rel=1e-9)
    assert b_c == pytest.approx(5.196, abs=1e-3)


def test_shadow_diameter_units(eht):
    src = eht.M87_parameters()
    diam = eht.predicted_shadow_diameter(src.mass_msun, src.distance_m, spin=0.0)
    # Should be smaller than the observed ring (10.39 r_g vs 11 r_g).
    ring = eht.predicted_ring_diameter(src.mass_msun, src.distance_m, spin=0.0)
    assert diam < ring
    assert diam > 0.9 * ring  # within ~10%


def test_kerr_shadow_axial_ratio_high_spin(eht):
    res = eht.kerr_shadow(spin=0.99, inclination_deg=90.0)
    # Maximum distortion ~ 25% squash perpendicular to spin axis.
    assert 0.7 < res["axial_ratio"] < 1.0
    res_low = eht.kerr_shadow(spin=0.0, inclination_deg=17.0)
    assert res_low["axial_ratio"] == pytest.approx(1.0, abs=1e-6)


def test_brightness_asymmetry_positive(eht):
    res = eht.ring_brightness_asymmetry(spin=0.9, inclination_deg=17.0)
    assert res["brightness_ratio_approaching_to_receding"] > 1.0
    assert 0.0 < res["v_keplerian_over_c"] < 1.0


def test_substrate_horizon_mark_undetectable(eht):
    src = eht.M87_parameters()
    res = eht.substrate_horizon_mark(src)
    assert res["detectable"] is False
    assert res["l_Planck_over_r_g"] < 1e-30
    assert res["expected_deviation_uas"] < 1e-20


def test_compare_to_observations_structure(eht):
    comp = eht.compare_to_observations()
    assert set(comp.keys()) == {"M87*", "Sgr A*"}
    for row in comp.values():
        assert row["abs_pct_error"] < 5.0
        for key in (
            "mass_Msun",
            "distance_m",
            "r_g_uas",
            "predicted_ring_uas",
            "observed_ring_uas",
            "fractional_error",
        ):
            assert key in row


def test_report_runs_and_mentions_sources(eht):
    text = eht.report()
    assert "M87*" in text
    assert "Sgr A*" in text
    assert "r_g" in text


def test_M87_r_g_about_4_uas(eht):
    src = eht.M87_parameters()
    # Well-known: M87* one r_g ~ 3.8 uas.
    assert 3.5 < src.r_g_uas < 4.2


def test_SgrA_r_g_about_5_uas(eht):
    src = eht.SgrA_parameters()
    # Well-known: Sgr A* one r_g ~ 5.0 uas.
    assert 4.6 < src.r_g_uas < 5.4
