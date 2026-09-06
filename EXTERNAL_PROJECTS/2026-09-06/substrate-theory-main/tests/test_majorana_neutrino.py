"""Tests for the B3 Majorana-neutrino / 0nubb predictor."""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.majorana_neutrino import (
    DEFAULT_M1_EV,
    DEFAULT_SIGMA_MNU_EV,
    DM21_SQ_EV2,
    DM31_SQ_EV2,
    EXPERIMENTS,
    MajoranaNeutrino,
    default_substrate_majorana,
)


def test_default_construction_runs():
    nu = default_substrate_majorana()
    assert nu.hierarchy == "NH"
    assert nu.masses_eV.shape == (3,)
    assert nu.U_e.shape == (3,)


def test_unitarity_of_Ue_row():
    """Sum_i |U_ei|^2 must equal 1 (PMNS row unitarity)."""
    nu = default_substrate_majorana()
    assert math.isclose(float(np.sum(nu.U_e)), 1.0, rel_tol=1e-12, abs_tol=1e-12)


def test_mass_ordering_NH():
    nu = default_substrate_majorana()
    m1, m2, m3 = nu.masses_eV
    assert m1 < m2 < m3


def test_mass_splittings_recovered_NH():
    nu = default_substrate_majorana()
    m1, m2, m3 = nu.masses_eV
    assert math.isclose(m2 * m2 - m1 * m1, DM21_SQ_EV2, rel_tol=1e-9)
    assert math.isclose(m3 * m3 - m1 * m1, DM31_SQ_EV2, rel_tol=1e-9)


def test_sigma_mnu_in_substrate_ballpark():
    """B3 predicts Sigma m_nu ~ 60.5 meV.  m_1=2.26 meV + splittings should
    give an internally consistent total in the same ballpark (within a few meV).
    """
    nu = default_substrate_majorana()
    sigma_meV = nu.sum_masses_meV()
    # Internal consistency: m_1 + sqrt(m1^2+Dm21) + sqrt(m1^2+Dm31) ~ 60.5 meV.
    assert 55.0 <= sigma_meV <= 65.0, f"Sigma m_nu = {sigma_meV} meV out of band"
    assert abs(sigma_meV - DEFAULT_SIGMA_MNU_EV * 1e3) < 5.0


def test_m_bb_zero_phase_positive():
    nu = default_substrate_majorana()
    val = nu.m_bb((0.0, 0.0))
    assert val > 0.0
    # constructive (all phases zero) should give the SUM, which equals max
    rng = nu.m_bb_range()
    assert math.isclose(val * 1e3, rng["max_meV"], rel_tol=1e-3)


def test_m_bb_range_in_audit_band():
    """Audit predicts m_bb in [0.055, 5.31] meV from a phase scan."""
    nu = default_substrate_majorana()
    rng = nu.m_bb_range(n_grid=361)
    assert 0.0 <= rng["min_meV"] < 1.0, f"min m_bb = {rng['min_meV']} meV"
    assert 1.0 < rng["max_meV"] < 6.0, f"max m_bb = {rng['max_meV']} meV"
    # full audit window
    assert rng["min_meV"] >= 0.0
    assert rng["max_meV"] <= 6.0


def test_kamland_zen_passes():
    """Current best limit (KamLAND-Zen 2024 < 36 meV) must NOT exclude B3."""
    nu = default_substrate_majorana()
    cmp = nu.compare_to_experiments()
    kz = cmp["KamLAND-Zen 2024"]
    assert kz["excluded"] is False
    assert "PASSES" in kz["status"]


def test_legend1000_or_nexo_will_discriminate():
    """At least one of LEGEND-1000 / nEXO must be tagged as a decisive test
    (its sensitivity band reaches B3 m_bb_max ~ 5 meV)."""
    nu = default_substrate_majorana()
    cmp = nu.compare_to_experiments()
    decisive = [
        name for name, row in cmp.items()
        if not row["is_limit"] and row["detectable"]
    ]
    assert any(n in decisive for n in ("LEGEND-1000", "nEXO")), (
        f"Expected LEGEND-1000 or nEXO to be decisive; got decisive={decisive}"
    )


def test_falsification_year_is_future():
    nu = default_substrate_majorana()
    f = nu.predict_falsification_year()
    assert f["decisive_year"] is not None
    assert f["decisive_year"] >= 2026
    assert f["experiment"] in {e.name for e in EXPERIMENTS if not e.is_limit}


def test_inverted_hierarchy_branch():
    nu = default_substrate_majorana()
    ih = nu.inverted_hierarchy_check(m3_lightest_eV=1.0e-3)
    # IH famously has a m_bb floor around ~15-20 meV.
    assert ih["IH_m_bb_min_meV"] > 5.0
    # IH is heavier overall than NH for tiny lightest mass.
    assert ih["IH_sum_mnu_meV"] > nu.sum_masses_meV() * 0.7


def test_phase_scan_min_below_max():
    nu = default_substrate_majorana()
    rng = nu.m_bb_range(n_grid=91)
    assert rng["min_meV"] < rng["max_meV"]


def test_report_renders():
    nu = default_substrate_majorana()
    txt = nu.report()
    assert "MAJORANA" in txt
    assert "m_bb" in txt
    assert "KamLAND-Zen" in txt


def test_invalid_hierarchy_raises():
    with pytest.raises(ValueError):
        MajoranaNeutrino(hierarchy="degenerate")


def test_alternate_m1_changes_sum():
    base = default_substrate_majorana()
    heavier = MajoranaNeutrino(m1_eV=5.0e-3)
    assert heavier.sum_masses_meV() > base.sum_masses_meV()
