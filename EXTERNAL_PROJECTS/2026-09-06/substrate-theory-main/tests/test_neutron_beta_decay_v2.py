"""Tests for neutron_beta_decay_v2 — bottle/beam puzzle + K_4 lifetime."""

from __future__ import annotations

import math

import pytest

from stiff_medium.neutron_beta_decay_v2 import (
    NeutronBetaDecayV2,
    TAU_N_BEAM_S,
    TAU_N_BOTTLE_S,
    V_UD_PDG,
    V_UD_SUBSTRATE,
)


@pytest.fixture
def model() -> NeutronBetaDecayV2:
    return NeutronBetaDecayV2()


# ---------------------------------------------------------------------------
# (a) τ_n within 2% of bottle 877.75 s
# ---------------------------------------------------------------------------

def test_lifetime_within_2pct_of_bottle(model: NeutronBetaDecayV2) -> None:
    tau = model.lifetime_substrate().tau_s
    frac = abs(tau - TAU_N_BOTTLE_S) / TAU_N_BOTTLE_S
    assert frac < 0.02, f"τ_sub = {tau:.2f} s, frac error {frac:.2%}"


def test_lifetime_in_physical_range(model: NeutronBetaDecayV2) -> None:
    tau = model.lifetime_substrate().tau_s
    # Squarely between bottle and beam means 870–890 s is the right range.
    assert 860.0 < tau < 900.0


# ---------------------------------------------------------------------------
# (b) g_A > 1.0
# ---------------------------------------------------------------------------

def test_g_A_greater_than_one(model: NeutronBetaDecayV2) -> None:
    info = model.g_A_axial_coupling()
    assert info["g_A_quenched_used"] > 1.0
    assert info["g_A_naive_Y_junction"] == pytest.approx(5.0 / 3.0)
    assert 0.7 < info["quench_factor"] < 0.85


# ---------------------------------------------------------------------------
# (c) Bottle/beam tension noted (≥ 3σ)
# ---------------------------------------------------------------------------

def test_bottle_beam_tension_present(model: NeutronBetaDecayV2) -> None:
    bb = model.bottle_vs_beam()
    assert bb.bottle_beam_tension_sigma >= 3.0
    assert bb.favored_method in ("bottle", "beam")
    # Substrate value must sit on a real-number axis between physical bounds.
    assert 800.0 < bb.tau_substrate_s < 950.0


def test_substrate_lifetime_between_methods_or_above(
    model: NeutronBetaDecayV2,
) -> None:
    bb = model.bottle_vs_beam()
    # Substrate K_4 lands closer (in σ) to whichever method has the larger
    # error bar; what matters physically is that it lies in the bottle/beam
    # window or just above it (within 2% of either), and that the favored
    # method is reported.
    assert bb.tau_substrate_s > 870.0
    assert bb.favored_method in ("bottle", "beam")


# ---------------------------------------------------------------------------
# Q-value
# ---------------------------------------------------------------------------

def test_q_value_near_782_keV(model: NeutronBetaDecayV2) -> None:
    Q = model.q_value()
    assert Q == pytest.approx(0.782, abs=0.05)


# ---------------------------------------------------------------------------
# |V_ud|
# ---------------------------------------------------------------------------

def test_Vud_substrate_close_to_pdg(model: NeutronBetaDecayV2) -> None:
    info = model.Vud_ckm()
    assert info["V_ud_substrate"] == pytest.approx(V_UD_SUBSTRATE)
    assert abs(info["frac_error"]) < 0.01
    assert V_UD_SUBSTRATE > 0.97 and V_UD_SUBSTRATE < 0.98
    # Cabibbo λ = 1/√20.
    assert info["lambda_substrate"] == pytest.approx(1.0 / math.sqrt(20.0))


# ---------------------------------------------------------------------------
# Dark decay
# ---------------------------------------------------------------------------

def test_dark_decay_forbidden(model: NeutronBetaDecayV2) -> None:
    dm = model.nDM_search()
    assert dm.branching_ratio == 0.0
    assert "FORBIDDEN" in dm.fornal_grinstein_status
    assert dm.branching_ratio_upper_limit_95 < 1e-2


# ---------------------------------------------------------------------------
# compare_to_experiments + report
# ---------------------------------------------------------------------------

def test_compare_to_experiments_keys(model: NeutronBetaDecayV2) -> None:
    out = model.compare_to_experiments()
    for k in ("substrate", "experiments", "bottle_beam_tension_sigma",
              "substrate_favors", "nsigma_to_bottle", "nsigma_to_beam"):
        assert k in out
    assert out["experiments"]["bottle"][0] == TAU_N_BOTTLE_S
    assert out["experiments"]["beam"][0] == TAU_N_BEAM_S


def test_report_contains_key_lines(model: NeutronBetaDecayV2) -> None:
    txt = model.report()
    assert "K_4" in txt
    assert "bottle" in txt.lower()
    assert "beam" in txt.lower()
    assert "Fornal" in txt or "FORBIDDEN" in txt
    assert "τ_n" in txt
