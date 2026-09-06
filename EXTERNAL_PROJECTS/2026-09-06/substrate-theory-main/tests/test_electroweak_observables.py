"""Tests for src/stiff_medium/electroweak_observables.py."""

from __future__ import annotations

import math

import pytest

from stiff_medium.electroweak_observables import (
    ElectroweakObservables,
    M_W_PDG,
    M_Z_PDG,
    SIN2_THETA_W_PDG,
    M_W_CDF2022,
    GAMMA_Z_INV_PDG,
    GAMMA_Z_NU_SM,
)


@pytest.fixture
def ew():
    return ElectroweakObservables()


def test_M_W_within_0p1_percent(ew):
    mw = ew.M_W()
    rel = abs(mw - M_W_PDG) / M_W_PDG
    assert rel < 1e-3, f"M_W = {mw}, rel error {rel:.3e}"


def test_M_Z_exact(ew):
    mz = ew.M_Z()
    # M_Z derived from M_W / cos θ_W; should agree with PDG to 0.05%.
    rel = abs(mz - M_Z_PDG) / M_Z_PDG
    assert rel < 5e-4, f"M_Z = {mz}, rel error {rel:.3e}"


def test_weinberg_angle_within_1_percent(ew):
    s2 = ew.weinberg_angle()
    rel = abs(s2 - SIN2_THETA_W_PDG) / SIN2_THETA_W_PDG
    assert rel < 1e-2, f"sin²θ_W = {s2}, rel error {rel:.3e}"


def test_rho_close_to_one(ew):
    rho = ew.rho_parameter()
    assert 0.999 < rho < 1.005


def test_Z_invisible_width_three_neutrinos(ew):
    g_inv = ew.Z_invisible_width(3)
    assert math.isclose(g_inv, 3 * GAMMA_Z_NU_SM, rel_tol=1e-9)
    # Effective N_nu from data should be near 3.
    n_eff = ew.n_neutrino_from_invisible_width()
    assert 2.9 < n_eff < 3.1


def test_CDF_anomaly_substrate_prediction(ew):
    res = ew.M_W_CDF_2022_anomaly()
    # Substrate-with-δT should be much closer to CDF than baseline is.
    assert res["M_W_substrate_with_delta_T_GeV"] > res["M_W_substrate_baseline_GeV"]
    assert res["agreement_with_CDF_pct"] < 0.2     # within 0.2% of CDF
    # Sanity: published CDF value is preserved in the report.
    assert math.isclose(res["M_W_CDF2022_GeV"], M_W_CDF2022)


def test_cross_section_at_Z_peak_is_large(ew):
    sigma_peak = ew.cross_section_e_e_hadrons(ew.M_Z())
    sigma_off  = ew.cross_section_e_e_hadrons(60.0)
    # Z peak should dominate over off-resonance.
    assert sigma_peak > 100.0 * sigma_off


def test_cross_section_increases_above_thresholds(ew):
    # R-ratio jumps when crossing charm and bottom thresholds.
    s_below_c = ew.cross_section_e_e_hadrons(2.0)
    s_above_c = ew.cross_section_e_e_hadrons(4.0)
    s_above_b = ew.cross_section_e_e_hadrons(15.0)
    # Per-event rate goes as R/s^2, so R-ratio jumps but absolute σ may
    # still fall.  Verify R ratio via σ * s^2 ordering.
    Rb_proxy = s_above_b * 15.0 ** 4
    Rc_proxy = s_above_c * 4.0 ** 4
    Ru_proxy = s_below_c * 2.0 ** 4
    assert Rc_proxy > Ru_proxy
    assert Rb_proxy > Rc_proxy


def test_compare_to_LEP_LHC_structure(ew):
    cmp = ew.compare_to_LEP_LHC()
    for key in ["M_W", "M_Z", "sin2_theta_W", "rho", "Z_invisible_width"]:
        assert key in cmp
    assert cmp["M_W"]["rel_error_pct"] < 0.1
    assert cmp["M_Z"]["rel_error_pct"] < 0.05


def test_report_contains_all_sections(ew):
    text = ew.report()
    for token in ("M_W", "M_Z", "sin2_theta_W", "rho",
                  "Z_invisible_width", "CDF-II"):
        assert token in text


def test_default_singleton_importable():
    from stiff_medium.electroweak_observables import default_ew
    assert isinstance(default_ew, ElectroweakObservables)
