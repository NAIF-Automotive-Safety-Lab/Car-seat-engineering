"""Tests for HiggsBoson substrate model."""

import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from stiff_medium.higgs_boson_substrate import (  # noqa: E402
    HiggsBoson,
    M_H_LHC_GEV,
    V_EW_GEV,
    GAMMA_H_SM_MEV,
)


@pytest.fixture
def H():
    return HiggsBoson()


# (a) m_H ~ 125 GeV ----------------------------------------------------
def test_mass_substrate_near_125_GeV(H):
    m = H.mass_substrate()
    assert abs(m - 125.0) < 5.0, f"m_H = {m:.3f} GeV, expected ~125"


def test_mass_within_few_percent_of_LHC(H):
    m = H.mass_substrate()
    rel = abs(m - M_H_LHC_GEV) / M_H_LHC_GEV
    assert rel < 0.05, f"substrate prediction off by {rel*100:.2f}%"


# (b) BR(bb) > BR(WW) --------------------------------------------------
def test_bb_dominates_WW(H):
    br = H.branching_ratios()
    assert br["bb"] > br["WW"], "BR(bb) must exceed BR(WW)"


def test_branching_ratios_ordered(H):
    br = H.branching_ratios()
    # bb > WW > gg > tau tau > cc > ZZ
    assert br["bb"] > br["WW"] > br["gg"] > br["tau_tau"] > br["cc"] > br["ZZ"]


def test_branching_ratios_sum_close_to_unity(H):
    br = H.branching_ratios()
    s = sum(br.values())
    assert 0.95 < s < 1.01, f"BR sum = {s:.4f} not near 1"


def test_branching_ratios_all_positive(H):
    for ch, val in H.branching_ratios().items():
        assert val > 0, f"BR({ch}) non-positive"


# (c) total width ~ 4 MeV ---------------------------------------------
def test_width_about_4_MeV(H):
    Gamma = H.decay_width()
    assert 3.0 < Gamma < 5.0, f"Gamma = {Gamma:.3f} MeV, expected ~4"


def test_width_matches_SM(H):
    assert math.isclose(H.decay_width(), GAMMA_H_SM_MEV, rel_tol=0.01)


# VEV ------------------------------------------------------------------
def test_vev(H):
    v = H.vacuum_expectation_value()
    assert math.isclose(v, V_EW_GEV, rel_tol=1e-4)


# Self-coupling --------------------------------------------------------
def test_self_coupling_lambda(H):
    lam = H.self_coupling_lambda()
    assert 0.10 < lam < 0.16, f"lambda = {lam:.4f} outside expected range"


def test_lambda_relation(H):
    """m_H^2 = 2 lambda v^2."""
    m = H.mass_substrate()
    v = H.vacuum_expectation_value()
    lam = H.self_coupling_lambda()
    lhs = m ** 2
    rhs = 2.0 * lam * v ** 2
    assert math.isclose(lhs, rhs, rel_tol=1e-9)


# Production cross sections -------------------------------------------
def test_production_cross_sections(H):
    xs = H.production_cross_sections(13.0)
    for ch in ("ggF", "VBF", "WH", "ZH", "ttH"):
        assert ch in xs and xs[ch] > 0
    # ggF dominates
    assert xs["ggF"] > 10 * xs["VBF"]


def test_production_total_reasonable(H):
    xs = H.production_cross_sections(13.0)
    total = sum(xs.values())
    assert 40.0 < total < 70.0, f"total xs = {total:.2f} pb"


# Naturalness ----------------------------------------------------------
def test_naturalness_exp(H):
    nat = H.naturalness_substrate()
    assert nat["exp_4pi2_minus_1"] > 1e15
    # Match Planck/v_EW within an order of magnitude (it's a few percent)
    assert abs(nat["delta_pct"]) < 50.0


# Comparison + report --------------------------------------------------
def test_compare_to_LHC_keys(H):
    cmp = H.compare_to_LHC()
    for key in ("mass_GeV", "width_MeV", "lambda", "vev_GeV"):
        assert key in cmp


def test_report_runs(H):
    txt = H.report()
    assert "Higgs Boson" in txt
    assert "MASS" in txt
    assert "NATURALNESS" in txt
