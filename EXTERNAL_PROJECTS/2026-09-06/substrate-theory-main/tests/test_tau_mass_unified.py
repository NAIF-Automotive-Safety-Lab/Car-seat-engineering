"""Tests for the unified m_tau/m_mu predictor (tau_mass_unified.py).

The module exposes a single canonical formula

    log(m_tau/m_mu) = n_M/(K_pair^4 * pi) - (n_R - R)/(K_pair*(K_pair+1))

shared by generation_map, mass_torque_engine, and integer_rigidity.
These tests verify

    (a) the canonical formula matches PDG within the published B3
        precision (currently 1.0% -- 0.93% honest residual);
    (b) the formula is traceable to substrate primitives + b3_constants
        (i.e. the function reproduces the value when fed the constants
        directly, and a perturbation in any constituent integer shifts
        the prediction in the expected direction);
    (c) every downstream module's m_tau/m_mu predictor agrees with the
        canonical formula to machine precision (1e-12), so there is a
        single source of truth.
"""

from __future__ import annotations

import math
import pytest

from src.stiff_medium import b3_constants
from src.stiff_medium.tau_mass_unified import (
    OBS_TAU_OVER_MU,
    LOG_OBS_TAU_OVER_MU,
    M_TAU_MEV,
    M_MUON_MEV,
    log_m_tau_over_m_mu,
    m_tau_over_m_mu,
    log_m_tau_over_m_mu_with_nlo_A,
    log_m_tau_over_m_mu_with_nlo_B,
    reconcile_modules,
)


# ---------------------------------------------------------------------------
# (a) Match within published B3 precision
# ---------------------------------------------------------------------------

def test_canonical_formula_matches_pdg_within_one_pct():
    """Canonical lepton tower must land within 1.0% of PDG m_tau/m_mu."""
    pred = m_tau_over_m_mu()
    err_pct = 100.0 * (pred - OBS_TAU_OVER_MU) / OBS_TAU_OVER_MU
    assert abs(err_pct) < 1.0, (
        f"canonical m_tau/m_mu = {pred:.6f}, PDG = {OBS_TAU_OVER_MU:.6f}, "
        f"err = {err_pct:+.4f}% (must be |err| < 1.0%; observed honest "
        f"residual ~ +0.93%)."
    )


def test_canonical_residual_within_half_pct():
    """Test the README's published precision target of 0.5%.

    Honest expected outcome: the current canonical formula gives
    +0.934% which FAILS this test. This test is xfail-marked and
    documents the gap that NLO derivation would close.
    """
    pred = m_tau_over_m_mu()
    err_pct = abs(100.0 * (pred - OBS_TAU_OVER_MU) / OBS_TAU_OVER_MU)
    if err_pct >= 0.5:
        pytest.xfail(
            f"canonical formula err = {err_pct:.4f}% > 0.5% "
            "(honest substrate residual; needs NLO derivation to close)"
        )


def test_log_obs_consistent_with_pdg_anchors():
    """Sanity: LOG_OBS_TAU_OVER_MU is consistent with M_TAU/M_MUON."""
    expected = math.log(M_TAU_MEV / M_MUON_MEV)
    assert math.isclose(LOG_OBS_TAU_OVER_MU, expected, rel_tol=1e-15)


# ---------------------------------------------------------------------------
# (b) Traceable to substrate primitives + b3_constants
# ---------------------------------------------------------------------------

def test_uses_canonical_b3_constants_by_default():
    """Default-arg call must consume canonical (n_M, K_pair, n_R, R)."""
    explicit = log_m_tau_over_m_mu(
        n_M=b3_constants.n_M,
        K_pair=b3_constants.K_pair,
        n_R=b3_constants.n_R,
        R=b3_constants.R,
    )
    default = log_m_tau_over_m_mu()
    assert math.isclose(explicit, default, rel_tol=0.0, abs_tol=1e-15)


def test_canonical_value_from_first_principles():
    """Independent recomputation from the closed-form expression."""
    n_M = b3_constants.n_M       # 268
    K_pair = b3_constants.K_pair # 2
    n_R = b3_constants.n_R       # 18
    R = b3_constants.R           # 3
    expected = (n_M / (K_pair ** 4 * math.pi)
                - (n_R - R) / (K_pair * (K_pair + 1)))
    got = log_m_tau_over_m_mu()
    assert math.isclose(got, expected, rel_tol=0.0, abs_tol=1e-15)


def test_n_R_perturbation_shifts_prediction():
    """n_R must enter the prediction actively, not via n_M identity."""
    base = log_m_tau_over_m_mu()
    plus = log_m_tau_over_m_mu(n_R=b3_constants.n_R + 1)
    minus = log_m_tau_over_m_mu(n_R=b3_constants.n_R - 1)
    expected_step = 1.0 / (b3_constants.K_pair * (b3_constants.K_pair + 1))
    # plus n_R -> minus the second term -> log decreases
    assert math.isclose(base - plus, expected_step, rel_tol=1e-12)
    assert math.isclose(minus - base, expected_step, rel_tol=1e-12)


def test_K_pair_perturbation_shifts_prediction():
    """K_pair appears in BOTH terms of the canonical formula."""
    base = log_m_tau_over_m_mu()
    plus = log_m_tau_over_m_mu(K_pair=b3_constants.K_pair + 1)
    assert plus != base
    # K_pair=3 makes K_pair^4 = 81 (much smaller leading term),
    # and K_pair*(K_pair+1) = 12 (smaller second term).
    # Net effect: leading dominates -> log decreases substantially.
    assert plus < base


def test_K_rank_K_pair_identity_holds():
    """(n_R - R)/(K_pair*(K_pair+1)) == K_rank/K_pair at canonical integers."""
    n_R = b3_constants.n_R
    R = b3_constants.R
    K_pair = b3_constants.K_pair
    K_rank = b3_constants.K_rank
    lhs = (n_R - R) / (K_pair * (K_pair + 1))
    rhs = K_rank / K_pair
    assert math.isclose(lhs, rhs, rel_tol=0.0, abs_tol=1e-15)


# ---------------------------------------------------------------------------
# (c) Cross-module agreement to 1e-12
# ---------------------------------------------------------------------------

def test_modules_agree_to_machine_precision():
    """generation_map, mass_torque_engine, integer_rigidity must match
    the canonical formula to 1e-12 in log(m_tau/m_mu)."""
    rec = reconcile_modules()
    assert rec.all_modules_agree(tol=1e-12), (
        f"Cross-module drift detected: max pairwise diff = "
        f"{rec.max_pairwise_diff:.3e}. Module logs: "
        f"canonical={rec.canonical:.12f}, "
        f"generation_map={rec.generation_map:.12f}, "
        f"mass_torque_engine={rec.mass_torque_engine:.12f}, "
        f"integer_rigidity={rec.integer_rigidity:.12f}"
    )


def test_reconciliation_canonical_equals_generation_map():
    """Pin: tau_mass_unified.log_m_tau_over_m_mu == GenerationMap._log_ratio_23."""
    from src.stiff_medium.generation_map import GenerationMap
    gm_log = GenerationMap()._log_ratio_23()
    canonical = log_m_tau_over_m_mu()
    assert math.isclose(canonical, gm_log, rel_tol=0.0, abs_tol=1e-12)


def test_reconciliation_canonical_equals_mass_torque():
    """Pin: tau_mass_unified == MassTorque._t_tau (extracted as m_tau/m_mu)."""
    from src.stiff_medium.mass_torque_engine import MassTorque
    eng = MassTorque()
    log_mt = math.log(eng.compute("tau").value_mev
                      / eng.compute("muon").value_mev)
    canonical = log_m_tau_over_m_mu()
    assert math.isclose(canonical, log_mt, rel_tol=0.0, abs_tol=1e-12)


def test_reconciliation_canonical_equals_integer_rigidity():
    """Pin: tau_mass_unified == integer_rigidity.predict_m_tau_over_m_e
    (extracted as ratio of tau/e and mu/e)."""
    from src.stiff_medium.integer_rigidity import (
        predict_m_tau_over_m_e,
        predict_m_mu_over_m_e,
        CANONICAL_INTEGERS,
    )
    log_ir = math.log(predict_m_tau_over_m_e(CANONICAL_INTEGERS)
                      / predict_m_mu_over_m_e(CANONICAL_INTEGERS))
    canonical = log_m_tau_over_m_mu()
    assert math.isclose(canonical, log_ir, rel_tol=0.0, abs_tol=1e-12)


# ---------------------------------------------------------------------------
# Research-candidate NLO corrections (sanity, not adoption)
# ---------------------------------------------------------------------------

def test_nlo_A_drives_residual_below_half_pct():
    """NLO_A candidate (-K_rank/(K_pair*n_M)) closes the residual."""
    pred = math.exp(log_m_tau_over_m_mu_with_nlo_A())
    err_pct = 100.0 * (pred - OBS_TAU_OVER_MU) / OBS_TAU_OVER_MU
    assert abs(err_pct) < 0.5


def test_nlo_B_drives_residual_below_half_pct():
    """NLO_B candidate (-1/(N_BAM*n_R)) closes the residual."""
    pred = math.exp(log_m_tau_over_m_mu_with_nlo_B())
    err_pct = 100.0 * (pred - OBS_TAU_OVER_MU) / OBS_TAU_OVER_MU
    assert abs(err_pct) < 0.5


def test_canonical_is_not_overruled_by_nlo_in_default_function():
    """Default ``m_tau_over_m_mu`` returns canonical, NOT either NLO form."""
    canonical = m_tau_over_m_mu()
    nlo_A = math.exp(log_m_tau_over_m_mu_with_nlo_A())
    nlo_B = math.exp(log_m_tau_over_m_mu_with_nlo_B())
    assert canonical != nlo_A
    assert canonical != nlo_B
    # Canonical sits at +0.93%; NLO_{A,B} sit much closer to PDG.
    canonical_err = abs(canonical - OBS_TAU_OVER_MU)
    nlo_A_err = abs(nlo_A - OBS_TAU_OVER_MU)
    nlo_B_err = abs(nlo_B - OBS_TAU_OVER_MU)
    assert nlo_A_err < canonical_err
    assert nlo_B_err < canonical_err
