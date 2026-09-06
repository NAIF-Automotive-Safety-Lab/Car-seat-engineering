"""Tests for stiff_medium.bbn_substrate."""

from __future__ import annotations

import math

import pytest

from stiff_medium.bbn_substrate import (
    BBNObservable,
    BBNSubstrate,
    OBS_D_H,
    OBS_LI7_H,
    OBS_Y_P,
    SBBN_LI7_H,
    default_bbn,
)


@pytest.fixture
def bbn() -> BBNSubstrate:
    return BBNSubstrate()


# ----- primary observables ---------------------------------------------------


def test_Y_p_matches_observed(bbn: BBNSubstrate) -> None:
    yp = bbn.Y_p_He4()
    assert 0.24 < yp < 0.25
    assert abs(yp - OBS_Y_P) / OBS_Y_P < 0.02  # within 2%


def test_D_H_matches_observed(bbn: BBNSubstrate) -> None:
    dh = bbn.D_H_ratio()
    assert 1e-5 < dh < 5e-5
    assert abs(dh - OBS_D_H) / OBS_D_H < 0.05  # within 5%


def test_He3_H_in_band(bbn: BBNSubstrate) -> None:
    val = bbn.He3_H_ratio()
    assert 1e-6 < val < 5e-5


def test_eta_is_planck_value(bbn: BBNSubstrate) -> None:
    assert math.isclose(bbn.eta_baryon_to_photon(), 6.10e-10, rel_tol=1e-3)


def test_N_eff_is_standard(bbn: BBNSubstrate) -> None:
    assert math.isclose(bbn.N_eff(), 3.046, rel_tol=1e-3)


# ----- lithium puzzle --------------------------------------------------------


def test_Li7_substrate_below_sbbn(bbn: BBNSubstrate) -> None:
    """Substrate prediction must be lower than the over-predicting SBBN value."""
    assert bbn.Li7_H_ratio() < SBBN_LI7_H


def test_Li7_substrate_consistent_with_observed(bbn: BBNSubstrate) -> None:
    """Substrate-corrected Li-7/H should land within ~50% of the Spite plateau."""
    val = bbn.Li7_H_ratio()
    assert abs(val - OBS_LI7_H) / OBS_LI7_H < 0.5


def test_lithium_puzzle_resolution_dict(bbn: BBNSubstrate) -> None:
    res = bbn.lithium_puzzle_resolution()
    for key in (
        "sbbn_prediction",
        "observed",
        "sbbn_overprediction_factor",
        "substrate_suppression_factor",
        "substrate_prediction",
        "substrate_residual_factor",
        "mechanism",
        "puzzle_addressed",
    ):
        assert key in res
    # SBBN over-predicts by ~3
    assert 2.0 < res["sbbn_overprediction_factor"] < 4.0
    # substrate must address the puzzle
    assert res["puzzle_addressed"] is True


def test_suppression_factor_is_tunable() -> None:
    bbn = BBNSubstrate(li7_substrate_suppression=2.0)
    expected = SBBN_LI7_H / 2.0
    assert math.isclose(bbn.Li7_H_ratio(), expected, rel_tol=1e-9)


# ----- comparison + report ---------------------------------------------------


def test_compare_to_observed_table(bbn: BBNSubstrate) -> None:
    comp = bbn.compare_to_observed()
    assert set(comp) == {"Y_p", "D/H", "He-3/H", "Li-7/H"}
    for obs in comp.values():
        assert isinstance(obs, BBNObservable)
        assert obs.relative_error >= 0.0
    # Y_p, D/H, Li-7/H all within ~50% (Li-7 is the loose one)
    assert comp["Y_p"].relative_error < 0.05
    assert comp["D/H"].relative_error < 0.10
    assert comp["Li-7/H"].relative_error < 0.5


def test_report_contains_keywords(bbn: BBNSubstrate) -> None:
    text = bbn.report()
    for kw in ("BBN", "eta", "N_eff", "Y_p", "D/H", "Li-7", "Lithium", "substrate"):
        assert kw in text


def test_default_bbn_helper() -> None:
    bbn = default_bbn()
    assert isinstance(bbn, BBNSubstrate)
    assert math.isclose(bbn.N_eff(), 3.046, rel_tol=1e-3)
