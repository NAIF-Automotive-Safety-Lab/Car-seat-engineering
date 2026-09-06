"""Tests for src/stiff_medium/particle_interactions.py."""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.particle_interactions import (
    GEV2_TO_NB,
    M_E_GEV,
    ParticleInteractions,
)


@pytest.fixture
def pi_obj():
    return ParticleInteractions()


# ---------------------------------------------------------------------------
# (a) e+ e- -> mu+ mu- at sqrt(s) = 10 GeV
# ---------------------------------------------------------------------------
def test_ee_to_mumu_at_10_GeV(pi_obj):
    """Pointlike QED: sigma = 4 pi alpha^2 / (3 s)."""
    sqrt_s = 10.0
    s = sqrt_s ** 2
    sigma_nb = pi_obj.electron_positron_to_muons(s)

    expected_nb = (4.0 * math.pi * pi_obj.alpha ** 2) / (3.0 * s) * GEV2_TO_NB
    assert sigma_nb == pytest.approx(expected_nb, rel=1e-3)
    # known value: ~0.868 nb at sqrt(s) = 10 GeV
    assert sigma_nb == pytest.approx(0.868, rel=2e-2)


def test_ee_to_mumu_below_threshold(pi_obj):
    assert pi_obj.electron_positron_to_muons(0.01) == 0.0


def test_ee_to_mumu_scales_as_one_over_s(pi_obj):
    s1, s2 = 100.0, 400.0
    r = pi_obj.electron_positron_to_muons(s1) / pi_obj.electron_positron_to_muons(s2)
    assert r == pytest.approx(s2 / s1, rel=1e-3)


# ---------------------------------------------------------------------------
# (b) R-ratio above bb threshold
# ---------------------------------------------------------------------------
def test_R_ratio_above_bb_threshold(pi_obj):
    res = pi_obj.electron_positron_to_hadrons(20.0 ** 2)
    assert res["R"] == pytest.approx(11.0 / 3.0, rel=1e-12)
    assert "b" in res["active_flavors"]
    assert "t" not in res["active_flavors"]


def test_R_ratio_below_charm_threshold(pi_obj):
    # sqrt(s) ~ 2 GeV: only u,d,s active -> R = 3*(4/9 + 1/9 + 1/9) = 2
    res = pi_obj.electron_positron_to_hadrons(2.0 ** 2)
    assert res["R"] == pytest.approx(2.0, rel=1e-12)


def test_R_ratio_between_charm_and_bottom(pi_obj):
    # sqrt(s) ~ 5 GeV: u,d,s,c -> R = 3*(4+1+1+4)/9 = 10/3
    res = pi_obj.electron_positron_to_hadrons(5.0 ** 2)
    assert res["R"] == pytest.approx(10.0 / 3.0, rel=1e-12)


# ---------------------------------------------------------------------------
# (c) Compton edge formula
# ---------------------------------------------------------------------------
def test_compton_edge_cs137(pi_obj):
    """Cs-137: E_gamma = 661.7 keV -> Compton edge ~ 477.4 keV."""
    res = pi_obj.compton_scattering(0.6617e-3)             # GeV
    edge_keV = res["compton_edge_T_max_GeV"] * 1e6
    assert edge_keV == pytest.approx(477.3, rel=5e-3)


def test_compton_edge_general_formula(pi_obj):
    """T_max = E * 2x / (1 + 2x), with x = E / m_e."""
    E = 0.001                    # 1 MeV photon
    res = pi_obj.compton_scattering(E)
    x = E / M_E_GEV
    expected_Tmax = E * 2 * x / (1 + 2 * x)
    assert res["compton_edge_T_max_GeV"] == pytest.approx(expected_Tmax, rel=1e-9)


def test_compton_low_energy_thomson_limit(pi_obj):
    """At E << m_e, sigma_KN approaches sigma_T."""
    res = pi_obj.compton_scattering(1e-9)
    ratio = res["sigma_klein_nishina_mb"] / res["sigma_thomson_mb"]
    assert ratio == pytest.approx(1.0, rel=1e-4)


# ---------------------------------------------------------------------------
# remaining method smoke tests
# ---------------------------------------------------------------------------
def test_electron_proton_scattering_keys(pi_obj):
    res = pi_obj.electron_proton_scattering(4.0, math.radians(30.0))
    for k in ("E_prime_GeV", "Q2_GeV2", "form_factor_GE",
              "dsigma_dOmega_nb_per_sr"):
        assert k in res
    assert 0.0 < res["form_factor_GE"] < 1.0
    assert res["E_prime_GeV"] < 4.0
    assert res["dsigma_dOmega_nb_per_sr"] > 0.0


def test_proton_proton_total_grows_with_s(pi_obj):
    a = pi_obj.proton_proton_total(100.0 ** 2)
    b = pi_obj.proton_proton_total(7000.0 ** 2)
    c = pi_obj.proton_proton_total(13000.0 ** 2)
    assert a < b < c
    assert 90 < b < 110           # ~98 mb at 7 TeV
    assert 100 < c < 120          # ~110 mb at 13 TeV


def test_weak_cc_linear_in_E(pi_obj):
    s1 = pi_obj.weak_charged_current(1.0)
    s2 = pi_obj.weak_charged_current(10.0)
    assert s2 / s1 == pytest.approx(10.0, rel=1e-9)
    # rough magnitude check: ~10^-41 cm^2 per GeV (CC nu_e e elastic)
    assert 1e-42 < s1 < 1e-40


def test_dis_callan_gross(pi_obj):
    res = pi_obj.dis_structure_functions(0.3, 10.0)
    assert res["F1"] == pytest.approx(res["F2"] / (2.0 * 0.3), rel=1e-9)
    assert res["F2"] > 0.0


def test_pdf_basic_shape(pi_obj):
    p_low  = pi_obj.pdf_at_x_Q2(0.01)
    p_high = pi_obj.pdf_at_x_Q2(0.7)
    # sea/gluon dominate at small x
    assert p_low["g"] > p_high["g"]
    # valence u larger than valence d
    assert p_high["u"] > p_high["d"]


def test_compare_to_pdg_and_report(pi_obj):
    cmp = pi_obj.compare_to_PDG()
    assert isinstance(cmp, list) and len(cmp) >= 4
    for entry in cmp:
        assert {"observable", "prediction", "measured", "percent_dev"} <= set(entry)
        assert math.isfinite(entry["percent_dev"])
        assert abs(entry["percent_dev"]) < 25.0   # all within ~25%
    rep = pi_obj.report()
    assert "ParticleInteractions vs PDG" in rep
    assert "R-ratio" in rep


# ---------------------------------------------------------------------------
# input validation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "method, args",
    [
        ("electron_proton_scattering", (-1.0, 0.5)),
        ("electron_proton_scattering", (1.0, 0.0)),
        ("compton_scattering",         (-1.0,)),
        ("proton_proton_total",        (-1.0,)),
        ("weak_charged_current",       (-1.0,)),
        ("dis_structure_functions",    (1.5, 10.0)),
        ("dis_structure_functions",    (0.3, -1.0)),
        ("pdf_at_x_Q2",                (1.5,)),
    ],
)
def test_input_validation(pi_obj, method, args):
    with pytest.raises(ValueError):
        getattr(pi_obj, method)(*args)
