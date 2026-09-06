"""Tests for src/stiff_medium/hubble_tension.py."""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.hubble_tension import (
    HubbleTension,
    MEASUREMENTS,
    SUBSTRATE_H0,
)


@pytest.fixture
def ht() -> HubbleTension:
    return HubbleTension()


def test_substrate_value() -> None:
    assert math.isclose(SUBSTRATE_H0, 71.92, rel_tol=1e-6)


def test_evaluate_returns_all_probes(ht: HubbleTension) -> None:
    rows = ht.evaluate_against_measurements()
    assert len(rows) == len(MEASUREMENTS)
    names = {r.name for r in rows}
    assert {"SH0ES", "Planck", "TRGB", "DESI+CMB",
            "ACT", "H0LiCOW", "Megamaser"} <= names


def test_within_2sigma_of_shoes(ht: HubbleTension) -> None:
    rows = {r.name: r for r in ht.evaluate_against_measurements()}
    assert rows["SH0ES"].n_sigma < 2.0, (
        f"substrate is {rows['SH0ES'].n_sigma:.2f}-sigma from SH0ES, expected <2"
    )


def test_three_sigma_or_more_from_planck(ht: HubbleTension) -> None:
    rows = {r.name: r for r in ht.evaluate_against_measurements()}
    # substrate at 71.92 vs Planck 67.40 +/- 0.50 -> ~6.4 sigma
    assert rows["Planck"].n_sigma >= 3.0, (
        f"substrate only {rows['Planck'].n_sigma:.2f}-sigma from Planck"
    )


def test_late_side_preferred_by_likelihood(ht: HubbleTension) -> None:
    combo = ht.combined_likelihood()
    assert combo["preferred_side"] == "late"
    assert combo["log_ratio_late_over_early"] > 0
    # late side should be very strongly preferred (Planck is ~6 sigma off)
    assert combo["log_ratio_late_over_early"] > 5.0


def test_sparc_cross_check_consistent(ht: HubbleTension) -> None:
    sparc = ht.cross_check_with_sparc()
    assert sparc["consistent_with_substrate"] is True
    assert sparc["sparc_vs_shoes_offset"] < sparc["sparc_vs_planck_offset"]
    assert sparc["sparc_side_preference"].startswith("late")


def test_early_universe_fix_specifies_delta(ht: HubbleTension) -> None:
    fix = ht.early_universe_modifications()
    assert "delta_h0_required" in fix
    # Planck = 67.4, substrate = 71.92 -> delta ~ +4.52
    assert "+4.5" in fix["delta_h0_required"]
    assert "early_dark_energy" in fix
    assert "extra_relativistic_dof" in fix


def test_late_universe_fix_specifies_delta(ht: HubbleTension) -> None:
    fix = ht.late_universe_modifications()
    # SH0ES = 73.04, substrate = 71.92 -> delta ~ -1.12
    assert "-1.1" in fix["delta_h0_required"]
    assert "cepheid_metallicity" in fix


def test_resolution_outlook_lists_experiments(ht: HubbleTension) -> None:
    outlook = ht.tension_resolution_outlook()
    for key in ("JWST_Cepheid_recalibration", "Euclid_weak_lensing",
                "DESI_DR3_BAO", "LISA_standard_sirens"):
        assert key in outlook
    assert "expected_resolution" in outlook


def test_report_contains_substrate_value(ht: HubbleTension) -> None:
    text = ht.report()
    assert "71.92" in text
    assert "SH0ES" in text and "Planck" in text
    assert "preferred side" in text.lower()
    # contains the table separators
    assert "Delta" in text and "n_sig" in text


def test_log_likelihood_signs(ht: HubbleTension) -> None:
    rows = {r.name: r for r in ht.evaluate_against_measurements()}
    # SH0ES likelihood should be much greater than Planck likelihood
    assert rows["SH0ES"].log_likelihood > rows["Planck"].log_likelihood
