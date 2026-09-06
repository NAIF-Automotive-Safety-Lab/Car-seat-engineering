"""Tests for stiff_medium.bbn_test (PDG/Planck comparison)."""

from __future__ import annotations

import math

import pytest

from stiff_medium.bbn_test import (
    BBNObservableFit,
    BBNTestResult,
    EtaScanResult,
    OBS_D_H_X1E5_CENTRAL,
    OBS_D_H_X1E5_SIGMA,
    OBS_HE3_H_X1E5_CENTRAL,
    OBS_LI7_H_X1E10_CENTRAL,
    OBS_LI7_H_X1E10_SIGMA,
    OBS_Y_P_CENTRAL,
    OBS_Y_P_SIGMA,
    default_test,
    report,
    run_bbn_test,
    scan_eta,
)


# ---------------------------------------------------------------------------
# Top-level BBNTestResult
# ---------------------------------------------------------------------------


def test_default_test_returns_result() -> None:
    res = default_test()
    assert isinstance(res, BBNTestResult)
    assert math.isclose(res.eta, 6.10e-10, rel_tol=1e-3)
    assert set(res.rows) == {"Y_p", "D/H", "He-3/H", "Li-7/H"}


def test_run_with_custom_eta_changes_d_h() -> None:
    res_low = run_bbn_test(eta=4.0e-10)
    res_hi = run_bbn_test(eta=8.0e-10)
    # D/H decreases with increasing eta (eta^{-1.6} fit).
    assert res_low.rows["D/H"].predicted > res_hi.rows["D/H"].predicted


def test_run_with_custom_eta_changes_li7() -> None:
    res_low = run_bbn_test(eta=4.0e-10)
    res_hi = run_bbn_test(eta=8.0e-10)
    # Li-7 increases with eta (eta^{2.1} fit).
    assert res_low.rows["Li-7/H"].predicted < res_hi.rows["Li-7/H"].predicted


# ---------------------------------------------------------------------------
# Per-observable n-sigma agreement
# ---------------------------------------------------------------------------


def test_Y_p_within_2_sigma() -> None:
    res = default_test()
    yp = res.rows["Y_p"]
    assert isinstance(yp, BBNObservableFit)
    assert yp.n_sigma <= 2.0
    assert yp.within_2s is True


def test_He3_H_within_2_sigma() -> None:
    res = default_test()
    he3 = res.rows["He-3/H"]
    assert he3.within_2s is True


def test_Li7_H_substrate_within_2_sigma() -> None:
    """Substrate-suppressed Li-7 must land within 2σ of the Spite plateau."""
    res = default_test()
    li7 = res.rows["Li-7/H"]
    assert li7.within_2s is True


def test_D_H_at_or_near_2_sigma() -> None:
    """D/H is the tightest constraint; substrate inherits SBBN value at ~2.4σ."""
    res = default_test()
    dh = res.rows["D/H"]
    # Allow up to ~3σ — this is a known SBBN tension that substrate inherits.
    assert dh.n_sigma < 3.0


def test_central_values_used() -> None:
    res = default_test()
    assert math.isclose(res.rows["Y_p"].observed, OBS_Y_P_CENTRAL)
    assert math.isclose(res.rows["Y_p"].sigma, OBS_Y_P_SIGMA)
    assert math.isclose(res.rows["D/H"].observed, OBS_D_H_X1E5_CENTRAL)
    assert math.isclose(res.rows["D/H"].sigma, OBS_D_H_X1E5_SIGMA)
    assert math.isclose(res.rows["He-3/H"].observed, OBS_HE3_H_X1E5_CENTRAL)
    assert math.isclose(res.rows["Li-7/H"].observed, OBS_LI7_H_X1E10_CENTRAL)
    assert math.isclose(res.rows["Li-7/H"].sigma, OBS_LI7_H_X1E10_SIGMA)


# ---------------------------------------------------------------------------
# Li-7 puzzle: substrate must beat SBBN
# ---------------------------------------------------------------------------


def test_sbbn_li7_overpredicts_by_factor_of_3() -> None:
    res = default_test()
    ratio = res.li7_sbbn_prediction_x1e10 / OBS_LI7_H_X1E10_CENTRAL
    assert 2.0 < ratio < 4.0


def test_substrate_li7_factor_3_below_sbbn() -> None:
    """Default suppression = 3.0 puts substrate at 1/3 of SBBN."""
    res = default_test()
    expected = res.li7_sbbn_prediction_x1e10 / 3.0
    assert math.isclose(res.li7_substrate_prediction_x1e10, expected, rel_tol=1e-9)


def test_substrate_beats_sbbn_on_li7() -> None:
    res = default_test()
    assert res.li7_substrate_beats_sbbn is True
    assert res.li7_substrate_n_sigma < res.li7_sbbn_n_sigma


def test_sbbn_li7_far_outside_2_sigma() -> None:
    """SBBN's ⁷Li/H is many σ from observation — that's the puzzle."""
    res = default_test()
    assert res.li7_sbbn_n_sigma > 5.0


def test_custom_suppression_changes_li7() -> None:
    res2 = run_bbn_test(li7_substrate_suppression=2.0)
    res3 = run_bbn_test(li7_substrate_suppression=3.0)
    res4 = run_bbn_test(li7_substrate_suppression=4.0)
    li7s = [r.li7_substrate_prediction_x1e10 for r in (res2, res3, res4)]
    # Monotonic decreasing with suppression factor.
    assert li7s[0] > li7s[1] > li7s[2]


# ---------------------------------------------------------------------------
# eta scan
# ---------------------------------------------------------------------------


def test_scan_eta_returns_grids() -> None:
    scan = scan_eta(n_pts=20)
    assert isinstance(scan, EtaScanResult)
    assert len(scan.eta_grid) == 20
    assert len(scan.D_H_x1e5) == 20
    assert len(scan.Li7_H_x1e10_sbbn) == 20
    assert len(scan.Li7_H_x1e10_sub) == 20


def test_scan_eta_d_h_monotonic_decreasing() -> None:
    scan = scan_eta(eta_min=3e-10, eta_max=1e-9, n_pts=15)
    dh = scan.D_H_x1e5
    for i in range(len(dh) - 1):
        assert dh[i] > dh[i + 1]


def test_scan_eta_li7_monotonic_increasing() -> None:
    scan = scan_eta(eta_min=3e-10, eta_max=1e-9, n_pts=15)
    li7 = scan.Li7_H_x1e10_sbbn
    for i in range(len(li7) - 1):
        assert li7[i] < li7[i + 1]


def test_scan_eta_Y_p_constant() -> None:
    scan = scan_eta(eta_min=3e-10, eta_max=1e-9, n_pts=10)
    yp = scan.Y_p
    for v in yp[1:]:
        assert math.isclose(v, yp[0], rel_tol=1e-9)


def test_scan_eta_substrate_li7_factor_of_3_below() -> None:
    scan = scan_eta(li7_substrate_suppression=3.0, n_pts=10)
    for sbbn, sub in zip(scan.Li7_H_x1e10_sbbn, scan.Li7_H_x1e10_sub):
        assert math.isclose(sub, sbbn / 3.0, rel_tol=1e-9)


def test_scan_eta_validation() -> None:
    with pytest.raises(ValueError):
        scan_eta(eta_min=1e-9, eta_max=3e-10)
    with pytest.raises(ValueError):
        scan_eta(n_pts=1)


# ---------------------------------------------------------------------------
# Verdict + report formatting
# ---------------------------------------------------------------------------


def test_verdict_mentions_li7() -> None:
    res = default_test()
    assert "Li-7" in res.verdict


def test_report_contains_observables() -> None:
    text = report(default_test())
    for kw in ("Y_p", "D/H", "³He/H", "⁷Li/H", "SBBN", "Substrate", "Verdict"):
        assert kw in text


def test_report_runs_without_error() -> None:
    text = report(default_test())
    assert len(text) > 200
