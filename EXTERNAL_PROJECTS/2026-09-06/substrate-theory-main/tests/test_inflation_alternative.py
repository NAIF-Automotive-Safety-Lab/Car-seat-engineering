"""Tests for inflation_alternative.InflationAlternative."""

from __future__ import annotations

import math

import pytest

from stiff_medium.inflation_alternative import (
    InflationAlternative,
    PLANCK_A_S,
    PLANCK_N_S,
    PLANCK_R_UPPER,
    B3_EXTRA_EFOLDS_NEEDED,
)


@pytest.fixture
def model() -> InflationAlternative:
    return InflationAlternative()


# ---------------------------------------------------------------------------
# (a) substrate r = 0
# ---------------------------------------------------------------------------
def test_substrate_tensor_to_scalar_is_exactly_zero(model: InflationAlternative):
    r = model.tensor_to_scalar_r()
    assert r.substrate == 0.0
    # the falsification map flags this as a sharp prediction
    fmap = model.falsification_summary()
    assert "r > 0" in fmap["substrate"]["primary"]


# ---------------------------------------------------------------------------
# (b) inflation r > 0 and within Planck/BICEP bound
# ---------------------------------------------------------------------------
def test_inflation_tensor_to_scalar_positive_and_bounded(model):
    r = model.tensor_to_scalar_r()
    assert isinstance(r.inflation, float)
    assert r.inflation > 0.0
    assert r.inflation <= PLANCK_R_UPPER


# ---------------------------------------------------------------------------
# (c) the proposed 8.7 e-folds of suppression close the amplitude gap
# ---------------------------------------------------------------------------
def test_8p7_efolds_suppression_matches_observed_A_s(model):
    """
    The 8.7-e-fold figure is the *target* number of e-folds whose
    derivation would close the full amplitude gap between the raw
    Poisson-defect amplitude (1/n_M) and the observed Planck A_s.
    With that many e-folds applied, the residual gap is under one dex.
    """
    res = model.density_perturbation_amplitude()
    gap_dex = abs(math.log10(PLANCK_A_S) - math.log10(res.substrate))
    assert gap_dex < 1.0, f"residual gap {gap_dex:.2f} dex too large"

    # And the class is using exactly 8.7 e-folds.
    assert model.SUBSTRATE_EXTRA_EFOLDS == pytest.approx(
        B3_EXTRA_EFOLDS_NEEDED
    )


def test_raw_substrate_amplitude_fails_by_many_orders(model):
    raw = 1.0 / model.SUBSTRATE_N_M
    gap_dex = math.log10(PLANCK_A_S) - math.log10(raw)
    # Raw 1/n_M sits ~17 dex below A_s; 8.7 e-folds of zeta suppression
    # (so 17.4 in A_s) overshoots slightly -- that's exactly the open
    # problem. Just check that the *raw* gap is huge.
    assert gap_dex > 10.0


# ---------------------------------------------------------------------------
# (d) both spectral indices within 5% of Planck central
# ---------------------------------------------------------------------------
def test_spectral_index_consistency_with_planck(model):
    res = model.spectral_index_n_s()
    for label, value in (("inflation", res.inflation),
                         ("substrate", res.substrate)):
        rel = abs(value - PLANCK_N_S) / PLANCK_N_S
        assert rel < 0.05, f"{label} n_s = {value} > 5% off Planck"


# ---------------------------------------------------------------------------
# Sanity: report renders, status board lists the open problem
# ---------------------------------------------------------------------------
def test_report_runs_and_contains_both_columns(model):
    out = model.report()
    assert "inflation" in out
    assert "substrate" in out
    assert "open problem" in out.lower()


def test_current_status_flags_amplitude_as_open(model):
    status = model.current_status()
    assert "inflation" in status and "substrate" in status
    assert "FAILS" in status["substrate"]["A_s"]
    assert "solved" in status["inflation"]["horizon"]


def test_all_comparisons_returns_seven_observables(model):
    results = model.all_comparisons()
    assert len(results) == 7
    names = {r.name for r in results}
    assert "A_s (scalar amplitude)" in names
    assert "r (tensor/scalar)" in names
    assert "reheating" in names
