"""Tests for the substrate vacuum-energy / cosmological-constant calculator."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.vacuum_energy import (
    VacuumEnergy,
    AlternativeProposal,
    RHO_LAMBDA_OBS_KG_M3,
    default,
)


# ---------------------------------------------------------------------------
# Core numerical predictions
# ---------------------------------------------------------------------------
def test_substrate_density_within_1pct_of_observed():
    v = default()
    pred = v.rho_Lambda_substrate()
    obs = v.rho_Lambda_observed()
    rel = abs(pred - obs) / obs
    assert rel < 0.01, f"substrate prediction off by {rel*100:.4f}% (>1%)"


def test_substrate_density_value_in_expected_range():
    v = default()
    pred = v.rho_Lambda_substrate()
    # Should be of order 1e-27 kg/m^3
    assert 1e-27 < pred < 1e-26


def test_observed_density_constant():
    assert default().rho_Lambda_observed() == pytest.approx(RHO_LAMBDA_OBS_KG_M3)


# ---------------------------------------------------------------------------
# Vacuum catastrophe
# ---------------------------------------------------------------------------
def test_vacuum_catastrophe_orders_around_120():
    v = default()
    orders = v.vacuum_catastrophe_orders()
    assert 115 < orders < 130, f"catastrophe orders = {orders} not ~120"


def test_substrate_beats_qft_by_at_least_115_orders():
    v = default()
    res = v.substrate_resolution()
    assert res["improvement_over_qft_orders"] > 115


def test_qft_naive_density_huge():
    v = default()
    assert v.rho_Lambda_qft_naive() > 1e90


# ---------------------------------------------------------------------------
# Equation of state
# ---------------------------------------------------------------------------
def test_w0_is_minus_one_exactly():
    assert default().dark_energy_equation_of_state() == -1.0


def test_w_of_z_is_constant_minus_one():
    v = default()
    z = np.linspace(0.0, 5.0, 50)
    w = v.evolution_with_redshift(z)
    assert isinstance(w, np.ndarray)
    assert np.all(w == -1.0)


def test_w_of_z_scalar():
    v = default()
    w0 = v.evolution_with_redshift(0.5)
    assert w0 == -1.0


def test_no_quintessence():
    info = default().dynamical_dark_energy_check()
    assert info["substrate_w0"] == -1.0
    assert info["substrate_wa"] == 0.0
    assert info["permits_quintessence"] is False
    assert any("DESI" in s for s in info["future_tests"])


# ---------------------------------------------------------------------------
# Alternatives
# ---------------------------------------------------------------------------
def test_compare_to_alternatives_includes_substrate():
    alts = default().compare_to_alternatives()
    names = [a.name for a in alts]
    assert any("substrate" in n.lower() for n in names)
    assert any("anthropic" in n.lower() or "multiverse" in n.lower() for n in names)
    assert any("qft" in n.lower() for n in names)


def test_substrate_has_zero_new_parameters():
    alts = default().compare_to_alternatives()
    sub = next(a for a in alts if "substrate" in a.name.lower())
    assert sub.new_parameters == 0
    assert sub.falsifiable is True
    assert sub.fine_tuning_orders < 1.0


def test_alternative_dataclass_immutable():
    a = AlternativeProposal(
        name="x", fine_tuning_orders=1.0, new_parameters=0,
        falsifiable=True, note="",
    )
    with pytest.raises(Exception):
        a.name = "y"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def test_report_contains_key_sections():
    text = default().report()
    assert "Vacuum Energy" in text
    assert "rho_Lambda" in text
    assert "Catastrophe" in text
    assert "w_0" in text


def test_substrate_resolution_dict_keys():
    res = default().substrate_resolution()
    for k in (
        "rho_substrate_kg_m3",
        "rho_observed_kg_m3",
        "relative_error",
        "relative_error_pct",
        "qft_catastrophe_orders",
        "improvement_over_qft_orders",
    ):
        assert k in res


# ---------------------------------------------------------------------------
# Constructor flexibility
# ---------------------------------------------------------------------------
def test_custom_lambda_qcd_changes_prediction():
    v1 = VacuumEnergy(lambda_qcd_MeV=200.0)
    v2 = VacuumEnergy(lambda_qcd_MeV=210.0)
    assert v2.rho_Lambda_substrate() > v1.rho_Lambda_substrate()


def test_c4_factor_scales_linearly():
    v1 = VacuumEnergy(c4_factor=1.0)
    v2 = VacuumEnergy(c4_factor=2.0)
    assert v2.rho_Lambda_substrate() == pytest.approx(
        2.0 * v1.rho_Lambda_substrate(), rel=1e-12
    )
