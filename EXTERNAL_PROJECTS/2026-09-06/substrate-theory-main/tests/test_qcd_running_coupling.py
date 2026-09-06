"""Tests for QCDRunningCoupling."""

import numpy as np
import pytest

from src.stiff_medium.qcd_running_coupling import (
    ALPHA_EM,
    ALPHA_S_MZ_PDG,
    M_Z,
    QCDRunningCoupling,
)


@pytest.fixture
def qcd():
    return QCDRunningCoupling()


def test_alpha_s_at_M_Z_matches_PDG(qcd):
    a = qcd.alpha_s_at_M_Z()
    assert abs(a - ALPHA_S_MZ_PDG) / ALPHA_S_MZ_PDG < 0.05


def test_alpha_s_returns_anchor_at_M_Z_scale(qcd):
    a = qcd.alpha_s_at_scale(M_Z, loop_order=2)
    assert abs(a - ALPHA_S_MZ_PDG) / ALPHA_S_MZ_PDG < 1e-6


def test_substrate_anchor_is_16_alpha(qcd):
    a = qcd.alpha_s_substrate_anchor()
    assert a == pytest.approx(16.0 * ALPHA_EM, rel=1e-12)
    # 16*alpha ~ 0.1168 -- within ~1% of measured alpha_s(M_Z)
    assert abs(a - ALPHA_S_MZ_PDG) / ALPHA_S_MZ_PDG < 0.05


def test_lambda_qcd_in_expected_range(qcd):
    L1 = qcd.lambda_qcd_from_running(nf=5, loop_order=1)
    L2 = qcd.lambda_qcd_from_running(nf=5, loop_order=2)
    # Expect Lambda_QCD(nf=5) ~ 90-260 MeV depending on loop order/scheme
    assert 0.05 < L1 < 0.30, f"L1={L1}"
    assert 0.05 < L2 < 0.30, f"L2={L2}"


def test_running_monotonically_decreasing(qcd):
    mus = np.logspace(0, 4, 20)  # 1 GeV ... 10 TeV
    vals = [qcd.alpha_s_at_scale(float(mu), loop_order=2) for mu in mus]
    diffs = np.diff(vals)
    assert np.all(diffs < 0), "alpha_s must monotonically decrease with mu"


def test_asymptotic_freedom(qcd):
    af = qcd.asymptotic_freedom(mu_max_GeV=1.0e6, n=8)
    assert af["monotonically_decreasing"]
    # alpha_s at 10^6 GeV should be small (< 0.1)
    assert af["alpha_s_at_max"] < 0.1
    # and well below alpha_s(M_Z)
    assert af["alpha_s_at_max"] < ALPHA_S_MZ_PDG


def test_confinement_scale_in_GeV_range(qcd):
    mu = qcd.confinement_scale(alpha_s_target=1.0)
    # Confinement (alpha_s ~ 1) sits within an order of magnitude of Lambda_QCD
    assert 0.1 < mu < 2.0, f"confinement scale = {mu} GeV"


def test_alpha_s_grows_at_low_energy(qcd):
    # 1-loop: alpha_s grows toward Lambda_QCD
    a_high = qcd.alpha_s_at_scale(50.0, loop_order=1)
    a_low = qcd.alpha_s_at_scale(2.0, loop_order=1)
    assert a_low > a_high


def test_one_loop_vs_two_loop_close_at_high_E(qcd):
    a1 = qcd.alpha_s_at_scale(M_Z, loop_order=1)
    a2 = qcd.alpha_s_at_scale(M_Z, loop_order=2)
    # Both anchored to MZ, so equal at MZ
    assert a1 == pytest.approx(a2, rel=1e-9)
    # At 200 GeV they should still be close (< 5%)
    a1h = qcd.alpha_s_at_scale(200.0, loop_order=1)
    a2h = qcd.alpha_s_at_scale(200.0, loop_order=2)
    assert abs(a1h - a2h) / a2h < 0.05


def test_compare_to_PDG_within_5pct(qcd):
    cmp = qcd.compare_to_PDG()
    assert cmp["within_5pct_at_MZ"]
    assert abs(cmp["model_pct_diff"]) < 5.0


def test_report_contains_key_fields(qcd):
    r = qcd.report()
    assert "alpha_s(M_Z)" in r
    assert "Lambda_QCD" in r
    assert "Asymptotic freedom" in r
    assert "Confinement" in r


def test_invalid_inputs(qcd):
    with pytest.raises(ValueError):
        qcd.alpha_s_at_scale(-1.0)
    with pytest.raises(ValueError):
        qcd.alpha_s_at_scale(10.0, loop_order=3)
    with pytest.raises(ValueError):
        qcd.lambda_qcd_from_running(loop_order=5)
