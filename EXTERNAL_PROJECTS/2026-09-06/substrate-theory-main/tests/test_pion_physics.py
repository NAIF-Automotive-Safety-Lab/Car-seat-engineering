"""Tests for src.stiff_medium.pion_physics."""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.pion_physics import (
    PionPhysics,
    PionPrediction,
    FormFactorPoint,
    PDG_M_PI_CHARGED_MEV,
    PDG_F_PI_MEV,
    PDG_TAU_PI_CHARGED_NS,
)


@pytest.fixture
def pion() -> PionPhysics:
    return PionPhysics()


# ---------------------------------------------------------------------------
# Required acceptance criteria from spec
# ---------------------------------------------------------------------------

def test_m_pi_charged_within_1_percent(pion: PionPhysics) -> None:
    """m_π± within 1 % of the empirical 139.57 MeV."""
    m = pion.m_pi_charged()
    rel_err = abs(m - PDG_M_PI_CHARGED_MEV) / PDG_M_PI_CHARGED_MEV
    assert rel_err < 0.01, f"m_π± = {m:.3f} MeV, err = {rel_err:.4%}"


def test_f_pi_within_5_percent(pion: PionPhysics) -> None:
    """f_π within 5 % of the empirical 92.2 MeV (ChiPT)."""
    f = pion.f_pi()
    rel_err = abs(f - PDG_F_PI_MEV) / PDG_F_PI_MEV
    assert rel_err < 0.05, f"f_π = {f:.3f} MeV, err = {rel_err:.4%}"


def test_pion_lifetime_within_10_percent(pion: PionPhysics) -> None:
    """τ_π± within 10 % of the empirical 26.0 ns."""
    tau = pion.pion_lifetime()
    rel_err = abs(tau - PDG_TAU_PI_CHARGED_NS) / PDG_TAU_PI_CHARGED_NS
    assert rel_err < 0.10, f"τ = {tau:.3f} ns, err = {rel_err:.4%}"


# ---------------------------------------------------------------------------
# Other observables
# ---------------------------------------------------------------------------

def test_neutral_pion_lighter_than_charged(pion: PionPhysics) -> None:
    assert pion.m_pi_neutral() < pion.m_pi_charged()


def test_neutral_pion_close_to_PDG(pion: PionPhysics) -> None:
    """π0 within 1 % of 134.98 MeV."""
    m0 = pion.m_pi_neutral()
    assert abs(m0 - 134.9768) / 134.9768 < 0.01


def test_g_piNN_reasonable(pion: PionPhysics) -> None:
    g = pion.g_piNN()
    assert 12.0 < g < 15.0, f"g_πNN = {g:.3f}"


# ---------------------------------------------------------------------------
# Form factor
# ---------------------------------------------------------------------------

def test_form_factor_normalisation(pion: PionPhysics) -> None:
    """F_π(0) = 1 by charge conservation."""
    assert pion.pion_form_factor(0.0) == pytest.approx(1.0)


def test_form_factor_monotone_decreasing(pion: PionPhysics) -> None:
    fs = [pion.pion_form_factor(q) for q in (0.0, 0.5, 1.0, 5.0, 20.0)]
    for a, b in zip(fs, fs[1:]):
        assert b < a


def test_form_factor_negative_Q2_raises(pion: PionPhysics) -> None:
    with pytest.raises(ValueError):
        pion.pion_form_factor(-1.0)


def test_form_factor_table_returns_points(pion: PionPhysics) -> None:
    table = pion.form_factor_table()
    assert all(isinstance(p, FormFactorPoint) for p in table)
    assert len(table) >= 3


# ---------------------------------------------------------------------------
# GMOR & χPT
# ---------------------------------------------------------------------------

def test_chiral_perturbation_theory_keys(pion: PionPhysics) -> None:
    d = pion.chiral_perturbation_theory()
    for k in ("m_pi_GeV", "f_pi_GeV", "qbar_q_GeV3", "GMOR_lhs_GeV4",
              "GMOR_rhs_GeV4", "GMOR_ratio"):
        assert k in d


def test_GMOR_self_consistent(pion: PionPhysics) -> None:
    """LHS / RHS of GMOR ≈ 1 by construction (since m_π was derived from GMOR)."""
    d = pion.chiral_perturbation_theory()
    assert abs(d["GMOR_ratio"] - 1.0) < 1.0e-6


# ---------------------------------------------------------------------------
# PDG comparison & report
# ---------------------------------------------------------------------------

def test_compare_to_PDG_returns_predictions(pion: PionPhysics) -> None:
    rows = pion.compare_to_PDG()
    assert all(isinstance(r, PionPrediction) for r in rows)
    # At least the four headline observables.
    names = {r.name for r in rows}
    assert any("m_π±" in n for n in names)
    assert any("f_π" in n for n in names)
    assert any("τ_π±" in n for n in names)


def test_compare_all_within_15_percent(pion: PionPhysics) -> None:
    """No headline prediction is more than 15 % away from PDG."""
    for r in pion.compare_to_PDG():
        assert abs(r.fractional_error) < 0.15, f"{r.name}: {r.percent_error:+.2f}%"


def test_report_contains_key_sections(pion: PionPhysics) -> None:
    txt = pion.report()
    assert "PionPhysics" in txt
    assert "GMOR" in txt
    assert "F_π" in txt
    assert "τ_π±" in txt


# ---------------------------------------------------------------------------
# Construction with custom inputs
# ---------------------------------------------------------------------------

def test_custom_sigma_changes_predictions() -> None:
    a = PionPhysics(sigma_GeV2=0.180).f_pi()
    b = PionPhysics(sigma_GeV2=0.200).f_pi()
    assert b > a
