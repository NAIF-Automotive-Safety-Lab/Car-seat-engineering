"""Tests for B3 substrate strong-CP module."""

from __future__ import annotations

import pytest

from src.stiff_medium.strong_cp_substrate import (
    NEDM_BOUND_PSI_2020,
    NEDM_SM_CKM_FLOOR,
    SM_THETA_BOUND,
    StrongCP,
    TheoryComparison,
)


@pytest.fixture
def cp() -> StrongCP:
    return StrongCP()


# ---------------------------------------------------------------------------
# (a) θ_QCD is forced to exactly zero by Möbius ℤ/2
# ---------------------------------------------------------------------------

def test_theta_substrate_prediction_is_exactly_zero(cp: StrongCP) -> None:
    theta = cp.theta_substrate_prediction()
    assert theta == 0.0
    # Trivially well below the implied SM bound:
    assert abs(theta) < SM_THETA_BOUND


def test_theta_zero_requires_K_pair_equals_two(cp: StrongCP) -> None:
    """Möbius mechanism is structural: only K_pair = 2 yields the ℤ/2."""
    cp.sheet_count = 3  # type: ignore[misc]
    with pytest.raises(RuntimeError):
        cp.theta_substrate_prediction()
    cp.sheet_count = 2  # restore


# ---------------------------------------------------------------------------
# (b) nEDM prediction sits below current PSI 2020 bound
# ---------------------------------------------------------------------------

def test_neutron_edm_strong_part_vanishes(cp: StrongCP) -> None:
    pred = cp.neutron_edm_prediction()
    assert pred["theta_QCD"] == 0.0
    assert pred["d_n_strong_e_cm"] == 0.0
    # Total survives only via CKM floor:
    assert pred["d_n_total_e_cm"] == NEDM_SM_CKM_FLOOR


def test_substrate_prediction_below_psi_2020_bound(cp: StrongCP) -> None:
    cmp = cp.compare_to_experiments()
    assert cmp["passes_bound"] is True
    # Margin should be enormous (~6 orders of magnitude)
    assert cmp["margin_below_bound"] > 1e5
    assert cmp["substrate_prediction_e_cm"] < NEDM_BOUND_PSI_2020


def test_compare_reports_psi_experiment(cp: StrongCP) -> None:
    cmp = cp.compare_to_experiments()
    assert "PSI" in cmp["experiment"]
    assert cmp["bound_e_cm"] == NEDM_BOUND_PSI_2020
    assert cmp["substrate_theta"] == 0.0


# ---------------------------------------------------------------------------
# (c) No axion needed in substrate
# ---------------------------------------------------------------------------

def test_axion_alternative_table_structure(cp: StrongCP) -> None:
    rows = cp.axion_alternative()
    assert len(rows) >= 5
    for row in rows:
        assert isinstance(row, TheoryComparison)
        assert row.aspect and row.standard_model and row.b3_substrate


def test_substrate_introduces_no_new_particle(cp: StrongCP) -> None:
    rows = cp.axion_alternative()
    new_particle_row = next(r for r in rows if "particle" in r.aspect.lower())
    assert "no" in new_particle_row.b3_substrate.lower()
    assert "axion" in new_particle_row.standard_model.lower()


def test_substrate_has_zero_free_parameters(cp: StrongCP) -> None:
    rows = cp.axion_alternative()
    fp_row = next(r for r in rows if "free parameter" in r.aspect.lower())
    assert "zero" in fp_row.b3_substrate.lower()


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def test_report_contains_key_phrases(cp: StrongCP) -> None:
    text = cp.report()
    assert "Möbius" in text
    assert "θ_QCD" in text
    assert "PSI" in text
    assert "axion" in text.lower()
    assert "Passes bound" in text and "True" in text
