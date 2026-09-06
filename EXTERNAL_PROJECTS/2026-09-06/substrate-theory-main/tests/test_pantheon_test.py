"""Tests for stiff_medium.pantheon_test (Pantheon+ H_0 falsification)."""

from __future__ import annotations

import math

import pytest

from stiff_medium.pantheon_test import (
    C_KM_S,
    HubbleDiagramResidual,
    OMEGA_M_PANTHEON,
    PANTHEON_ALONE_H0,
    PANTHEON_BINNED,
    PantheonTestResult,
    PLANCK_H0,
    SHOES_PANTHEON_H0,
    SUBSTRATE_H0,
    TensionAgainstAnchor,
    chi2_from_residuals,
    comoving_distance_mpc,
    distance_modulus,
    evaluate_hubble_diagram,
    evaluate_tensions,
    hubble_parameter,
    luminosity_distance_mpc,
    report,
    run_pantheon_test,
    sigma_distance,
)


# ---------------------------------------------------------------------------
# Anchor constants (sanity)
# ---------------------------------------------------------------------------

def test_anchor_h0_ordering() -> None:
    """Planck (early) < substrate < SH0ES (late) — substrate is in the band."""
    assert PLANCK_H0 < SUBSTRATE_H0 < SHOES_PANTHEON_H0
    assert PLANCK_H0 < SUBSTRATE_H0 < PANTHEON_ALONE_H0


def test_substrate_h0_value() -> None:
    """Substrate prediction is 71.92 km/s/Mpc (Sigma m_nu chain)."""
    assert math.isclose(SUBSTRATE_H0, 71.92, abs_tol=1e-9)


# ---------------------------------------------------------------------------
# Sigma-distance helper
# ---------------------------------------------------------------------------

def test_sigma_distance_basic() -> None:
    # |10 - 5| / sqrt(1+1) = 5 / sqrt(2) ~ 3.536
    assert math.isclose(sigma_distance(10.0, 5.0, 1.0, 1.0),
                        5.0 / math.sqrt(2.0), rel_tol=1e-9)


def test_sigma_distance_zero_total() -> None:
    assert sigma_distance(1.0, 0.0, 0.0, 0.0) == float("inf")


# ---------------------------------------------------------------------------
# Tension table vs anchor probes
# ---------------------------------------------------------------------------

def test_evaluate_tensions_returns_three_anchors() -> None:
    rows = evaluate_tensions()
    names = {r.name for r in rows}
    assert names == {"SH0ES+Pantheon+", "Pantheon+ alone", "Planck"}
    assert all(isinstance(r, TensionAgainstAnchor) for r in rows)


def test_substrate_within_2_sigma_of_shoes() -> None:
    """Substrate H_0 = 71.92 should sit within 2 sigma of SH0ES (1.04 + 0.50)."""
    rows = evaluate_tensions()
    shoes = next(r for r in rows if r.name == "SH0ES+Pantheon+")
    assert shoes.n_sigma < 2.0
    assert shoes.delta_h0 < 0.0  # substrate is below SH0ES


def test_substrate_far_from_planck() -> None:
    """Substrate H_0 = 71.92 should sit > 5 sigma from Planck (0.50 + 0.50)."""
    rows = evaluate_tensions()
    planck = next(r for r in rows if r.name == "Planck")
    assert planck.n_sigma > 5.0
    assert planck.delta_h0 > 0.0  # substrate is above Planck


def test_pantheon_alone_close_to_shoes() -> None:
    """Pantheon+ alone H_0 ~ 73.4, similar to SH0ES, so substrate still within
    ~1.5 sigma of it."""
    rows = evaluate_tensions()
    pa = next(r for r in rows if r.name == "Pantheon+ alone")
    assert pa.n_sigma < 2.0


# ---------------------------------------------------------------------------
# Cosmology helpers
# ---------------------------------------------------------------------------

def test_hubble_parameter_at_z_zero_equals_h0() -> None:
    h0 = 70.0
    assert math.isclose(hubble_parameter(0.0, h0, 0.3), h0, rel_tol=1e-9)


def test_hubble_parameter_increases_with_z() -> None:
    h0 = 70.0
    assert hubble_parameter(0.5, h0, 0.3) > h0
    assert hubble_parameter(1.0, h0, 0.3) > hubble_parameter(0.5, h0, 0.3)


def test_comoving_distance_zero_at_z_zero() -> None:
    assert comoving_distance_mpc(0.0, 70.0, 0.3) == 0.0


def test_comoving_distance_low_z_hubble_law() -> None:
    """At small z, D_C ~ c*z/H_0 (Hubble's law)."""
    z = 0.001
    h0 = 70.0
    expected = C_KM_S * z / h0
    actual = comoving_distance_mpc(z, h0, 0.3)
    assert math.isclose(actual, expected, rel_tol=2e-3)


def test_luminosity_distance_factor() -> None:
    """D_L = (1+z) * D_C."""
    z = 0.5
    h0, om = 70.0, 0.3
    d_c = comoving_distance_mpc(z, h0, om)
    d_l = luminosity_distance_mpc(z, h0, om)
    assert math.isclose(d_l, (1.0 + z) * d_c, rel_tol=1e-9)


def test_distance_modulus_monotone_in_z() -> None:
    h0, om = 71.92, OMEGA_M_PANTHEON
    mus = [distance_modulus(z, h0, om) for z in (0.05, 0.1, 0.3, 0.6, 1.0)]
    for a, b in zip(mus, mus[1:]):
        assert b > a


def test_distance_modulus_decreases_with_h0() -> None:
    """Higher H_0 -> nearer universe -> smaller mu at fixed z."""
    z = 0.3
    om = OMEGA_M_PANTHEON
    mu_low_h0 = distance_modulus(z, 67.4, om)
    mu_high_h0 = distance_modulus(z, 73.04, om)
    assert mu_high_h0 < mu_low_h0


# ---------------------------------------------------------------------------
# Hubble-diagram residual sweep
# ---------------------------------------------------------------------------

def test_pantheon_binned_count() -> None:
    assert len(PANTHEON_BINNED) == 12


def test_pantheon_binned_z_monotone() -> None:
    zs = [b[0] for b in PANTHEON_BINNED]
    for a, b in zip(zs, zs[1:]):
        assert b > a


def test_evaluate_hubble_diagram_returns_per_bin_rows() -> None:
    rows = evaluate_hubble_diagram()
    assert len(rows) == len(PANTHEON_BINNED)
    assert all(isinstance(r, HubbleDiagramResidual) for r in rows)


def test_substrate_residuals_small_at_low_z() -> None:
    """Substrate H_0 disagrees with Pantheon+ best fit by ~1.5 km/s/Mpc;
    distance-modulus residuals at low z should be < 0.1 mag."""
    rows = evaluate_hubble_diagram()
    for r in rows[:6]:
        assert abs(r.residual_substrate) < 0.1


def test_planck_residuals_larger_than_substrate() -> None:
    """Planck H_0 = 67.40 disagrees with Pantheon+ best fit more than substrate."""
    rows = evaluate_hubble_diagram()
    for r in rows:
        assert abs(r.residual_planck) > abs(r.residual_substrate)


def test_chi2_substrate_smaller_than_planck() -> None:
    rows = evaluate_hubble_diagram()
    chi2_sub = chi2_from_residuals(rows, "substrate")
    chi2_planck = chi2_from_residuals(rows, "planck")
    assert chi2_sub < chi2_planck


def test_chi2_substrate_passes_per_bin_threshold() -> None:
    """Substrate per-bin chi^2 / N should be < 1 (good fit)."""
    rows = evaluate_hubble_diagram()
    chi2_sub = chi2_from_residuals(rows, "substrate")
    assert chi2_sub / len(rows) < 1.0


def test_chi2_unknown_which_raises() -> None:
    rows = evaluate_hubble_diagram()
    with pytest.raises(ValueError):
        chi2_from_residuals(rows, "asdf")


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------

def test_run_pantheon_test_returns_result() -> None:
    res = run_pantheon_test()
    assert isinstance(res, PantheonTestResult)


def test_substrate_lands_in_tension_band() -> None:
    res = run_pantheon_test()
    assert res.in_tension_band is True


def test_substrate_closer_to_late_side() -> None:
    res = run_pantheon_test()
    assert "late" in res.closer_to_side
    assert "SH0ES" in res.closer_to_side


def test_n_sigma_accessors() -> None:
    res = run_pantheon_test()
    assert res.n_sigma_vs_shoes < 2.0
    assert res.n_sigma_vs_planck > 5.0
    assert res.n_sigma_vs_pantheon_alone < 2.0


def test_verdict_mentions_in_tension_band() -> None:
    res = run_pantheon_test()
    assert "INSIDE" in res.verdict


def test_chi2_ranking_substrate_between_shoes_and_planck() -> None:
    """Substrate chi^2 should be > SH0ES (since central is Pantheon+ best fit)
    but << Planck.  Confirms substrate is the second-best fit of the three."""
    res = run_pantheon_test()
    assert res.chi2_shoes <= res.chi2_substrate
    assert res.chi2_substrate < res.chi2_planck


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def test_report_includes_substrate_label() -> None:
    text = report()
    assert "substrate prediction" in text
    assert "71.92" in text


def test_report_includes_per_anchor_lines() -> None:
    text = report()
    assert "SH0ES+Pantheon+" in text
    assert "Pantheon+ alone" in text
    assert "Planck" in text


def test_report_emits_chi2_block() -> None:
    text = report()
    assert "chi^2" in text
    assert "reduced chi^2" in text


def test_report_emits_verdict_line() -> None:
    text = report()
    assert "Verdict" in text
