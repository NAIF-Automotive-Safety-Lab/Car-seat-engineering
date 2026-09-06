"""Tests for sparc_test: substrate g_† = c·H_0/(2π) vs McGaugh+ 2016 SPARC."""

from __future__ import annotations

import math
import os

import numpy as np
import pytest

from stiff_medium.sparc_test import (
    HUBBLE_CANDIDATES,
    MCGAUGH,
    g_dagger_substrate,
    hubble_tension_verdict,
    predict_g_dagger_all,
    predict_mdar,
    predict_tully_fisher_slope,
    render_sparc_test,
    report,
)


# ---------------------------------------------------------------------------
# Empirical reference sanity
# ---------------------------------------------------------------------------


def test_mcgaugh_reference_constants_are_sane():
    """McGaugh+2016 reference should match the published numbers."""
    assert MCGAUGH.n_galaxies == 175
    assert math.isclose(MCGAUGH.g_dagger_central, 1.20e-10, rel_tol=1e-6)
    assert math.isclose(MCGAUGH.g_dagger_sigma, 0.02e-10, rel_tol=1e-6)
    assert MCGAUGH.tully_fisher_slope == pytest.approx(4.0, abs=1e-9)


def test_three_hubble_candidates_present():
    labels = {c.label for c in HUBBLE_CANDIDATES}
    assert labels == {"substrate", "SHOES", "Planck"}


# ---------------------------------------------------------------------------
# Substrate g_† closed form
# ---------------------------------------------------------------------------


def test_g_dagger_closed_form_matches_sparc_dynamics_engine():
    """The standalone closed form must match SPARCDynamics.g_dagger()."""
    from stiff_medium.sparc_dynamics import SPARCDynamics

    for h0 in (50.0, 67.4, 71.92, 73.04, 100.0):
        g_closed = g_dagger_substrate(h0)
        g_engine = SPARCDynamics(h0).g_dagger_from_substrate()
        assert math.isclose(g_closed, g_engine, rel_tol=1e-12), (
            f"H_0={h0}: closed-form {g_closed} vs engine {g_engine}"
        )


def test_g_dagger_scales_linearly_with_H0():
    """g_† = c·H_0/(2π) → doubling H_0 doubles g_†."""
    g1 = g_dagger_substrate(50.0)
    g2 = g_dagger_substrate(100.0)
    assert math.isclose(g2 / g1, 2.0, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# All three Hubble candidate predictions
# ---------------------------------------------------------------------------


def test_predict_g_dagger_all_returns_3_with_signs():
    preds = predict_g_dagger_all()
    assert len(preds) == 3
    # All three should land BELOW the McGaugh empirical (negative residual)
    # because c·H_0/(2π) ≈ 1.04–1.13 < 1.20.
    for p in preds:
        assert p.fractional_error < 0.0, p


def test_all_candidates_within_mcgaugh_systematic_envelope():
    """All 3 H_0 candidates should land within McGaugh's ±20% systematic."""
    preds = predict_g_dagger_all()
    for p in preds:
        assert p.passes_systematic, (
            f"{p.candidate.label} (H_0={p.candidate.H0_kms_per_mpc}): "
            f"{p.percent_error:+.2f}% (sys envelope = ±20%)"
        )


def test_shoes_residual_smaller_than_planck_residual():
    """The SH0ES side must beat the Planck side on |residual|."""
    by_label = {p.candidate.label: p for p in predict_g_dagger_all()}
    assert abs(by_label["SHOES"].fractional_error) < abs(
        by_label["Planck"].fractional_error
    )


def test_substrate_h0_gives_intermediate_residual():
    """B3 substrate H_0=71.92 lies between SH0ES and Planck residuals."""
    by_label = {p.candidate.label: p for p in predict_g_dagger_all()}
    sh = abs(by_label["SHOES"].fractional_error)
    pl = abs(by_label["Planck"].fractional_error)
    sub = abs(by_label["substrate"].fractional_error)
    assert sh <= sub <= pl, (
        f"SHOES={sh*100:.2f}%, substrate={sub*100:.2f}%, Planck={pl*100:.2f}%"
    )


def test_shoes_residual_within_8_percent():
    """SH0ES H_0 should give substrate g_† within 8% of McGaugh empirical."""
    by_label = {p.candidate.label: p for p in predict_g_dagger_all()}
    assert abs(by_label["SHOES"].fractional_error) < 0.08, by_label["SHOES"]


def test_planck_residual_at_least_10_percent():
    """Planck H_0 should be off by more than 10% — it loses to SH0ES."""
    by_label = {p.candidate.label: p for p in predict_g_dagger_all()}
    assert abs(by_label["Planck"].fractional_error) > 0.10, by_label["Planck"]


# ---------------------------------------------------------------------------
# Tully-Fisher slope = 1/4 (analytically forced)
# ---------------------------------------------------------------------------


def test_tully_fisher_slope_is_one_quarter_to_4_decimals():
    """Substrate forces V_flat ∝ M_b^(1/4) → slope must be 0.25 exactly."""
    tf = predict_tully_fisher_slope()
    assert tf.slope_substrate == pytest.approx(0.25, abs=1e-4)


def test_tully_fisher_slope_independent_of_H0():
    """The 1/4 slope is geometric — must hold for any H_0."""
    for h0 in (50.0, 67.4, 71.92, 73.04, 100.0):
        tf = predict_tully_fisher_slope(H0_kms_per_mpc=h0)
        assert tf.slope_substrate == pytest.approx(0.25, abs=1e-4)


def test_tully_fisher_slope_within_1_sigma_of_observed():
    """Substrate slope must agree with observed 4.0±0.1 (i.e. 0.25±0.025)."""
    tf = predict_tully_fisher_slope()
    assert tf.n_sigma < 1.0, tf


# ---------------------------------------------------------------------------
# MDAR / RAR limits
# ---------------------------------------------------------------------------


def test_mdar_returns_9_points_default():
    pts = predict_mdar()
    assert len(pts) == 9


def test_mdar_high_acceleration_limit_is_newtonian():
    """At g_bar >> g_†, g_obs should be ≈ g_bar."""
    pts = predict_mdar()
    # Highest-g_bar point (g_bar = 1e-8) should be in Newtonian regime
    p_high = pts[-1]
    ratio = p_high.g_obs_substrate / p_high.g_bar
    assert ratio == pytest.approx(1.0, abs=0.05), (
        f"g_obs/g_bar = {ratio:.4f} at log10(g_bar)={p_high.log10_g_bar}"
    )


def test_mdar_low_acceleration_limit_is_deep_mond():
    """At g_bar << g_†, g_obs ≈ √(g_bar · g_†) (factor √(g_†/g_bar) above)."""
    pts = predict_mdar()
    p_low = pts[0]            # log10(g_bar) = -13
    g_dag = g_dagger_substrate(73.04)
    g_mond = math.sqrt(p_low.g_bar * g_dag)
    # Should be within 20% of the deep-MOND asymptote
    rel = abs(p_low.g_obs_substrate - g_mond) / g_mond
    assert rel < 0.20, f"deep MOND check: rel = {rel:.3f}"


def test_mdar_monotone_increasing_in_g_bar():
    """g_obs(g_bar) must increase monotonically."""
    pts = predict_mdar(n_pts=50)
    g_obs = np.array([p.g_obs_substrate for p in pts])
    diffs = np.diff(g_obs)
    assert np.all(diffs > 0), "MDAR must be monotone increasing"


# ---------------------------------------------------------------------------
# Hubble tension cross-check
# ---------------------------------------------------------------------------


def test_hubble_verdict_prefers_shoes_side():
    """SPARC g_† constraint should prefer SH0ES over Planck."""
    v = hubble_tension_verdict()
    assert v["sparc_prefers_shoes_over_planck"] is True
    assert v["best_candidate"] == "SHOES"
    assert v["shoes_minus_planck_gap_pct"] > 0.0


def test_hubble_verdict_all_within_systematic():
    v = hubble_tension_verdict()
    assert v["all_within_systematic"] is True


def test_hubble_gap_at_least_5_percentage_points():
    """SH0ES must outperform Planck by at least 5 pct points on |residual|."""
    v = hubble_tension_verdict()
    assert v["shoes_minus_planck_gap_pct"] >= 5.0, v


# ---------------------------------------------------------------------------
# Report and visual artefacts
# ---------------------------------------------------------------------------


def test_report_runs_and_mentions_key_quantities():
    txt = report()
    assert "g_†" in txt
    assert "Tully-Fisher" in txt
    assert "MDAR" in txt or "Mass-discrepancy" in txt or "RAR" in txt or \
        "mass-discrepancy" in txt
    assert "SH0ES" in txt or "SHOES" in txt
    assert "Planck" in txt
    assert "substrate" in txt.lower()


def test_render_visual_creates_png(tmp_path):
    out = tmp_path / "129_sparc_test.png"
    path = render_sparc_test(out_path=str(out))
    assert os.path.isfile(path)
    assert os.path.getsize(path) > 5_000      # >5 KB sanity
