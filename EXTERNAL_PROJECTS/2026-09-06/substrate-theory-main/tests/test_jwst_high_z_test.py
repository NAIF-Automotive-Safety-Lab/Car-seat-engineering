"""Tests for src.stiff_medium.jwst_high_z_test (JWST z>10 comparison)."""

import math

import pytest

from src.stiff_medium.jwst_high_z_test import (
    SourceComparison,
    UVLFComparison,
    DensityComparison,
    JWSTHighZTestResult,
    JWST_OBSERVED_OVERABUNDANCE_DEX,
    M_UV_to_log_M_star,
    compare_per_source,
    compare_UV_LF,
    compare_stellar_mass_density,
    harikane_donnan_z10_LF,
    integrated_stellar_density,
    jwst_z_gt_10_sample,
    labbe_2023_density_points,
    make_lcdm_track,
    make_substrate_track,
    report,
    run_jwst_high_z_test,
)


# ---------------------------------------------------------------------------
# Reference data shape
# ---------------------------------------------------------------------------

def test_jwst_sample_includes_named_sources():
    names = {s.name for s in jwst_z_gt_10_sample()}
    # Must include the canonical confirmed z>10 detections from the spec
    for required in ("JADES-GS-z14-0", "GN-z11", "Maisie's Galaxy",
                     "CEERS-93316"):
        assert required in names


def test_jwst_sample_z_values_sane():
    sample = jwst_z_gt_10_sample()
    # Must include >= 4 z>10 sources for the test to be meaningful
    n_high_z = sum(1 for s in sample if s.z > 10.0)
    assert n_high_z >= 4
    # Highest-z must be z14 source (Carniani 2024)
    z_max_source = max(sample, key=lambda s: s.z)
    assert z_max_source.z >= 14.0


def test_uv_LF_points_in_z_range():
    pts = harikane_donnan_z10_LF()
    assert len(pts) >= 4
    for pt in pts:
        assert 9.0 <= pt.z <= 16.0
        assert -23.0 < pt.M_UV < -19.0
        assert -8.0 < pt.log_phi < -3.0


def test_labbe_density_points_high_mass():
    pts = labbe_2023_density_points()
    assert len(pts) >= 2
    for pt in pts:
        assert pt.z >= 7.0
        assert pt.log_M_star_min >= 9.0
        assert pt.log_rho_star > 0


# ---------------------------------------------------------------------------
# UV magnitude -> stellar mass
# ---------------------------------------------------------------------------

def test_M_UV_to_log_M_star_brighter_more_massive():
    """A brighter (more negative M_UV) galaxy must have larger M_*."""
    bright = M_UV_to_log_M_star(-22.0)
    faint = M_UV_to_log_M_star(-19.0)
    assert bright > faint
    # Anchor: M_UV = -21 -> log M ~ 9.0
    assert 8.5 < M_UV_to_log_M_star(-21.0) < 9.5


# ---------------------------------------------------------------------------
# Per-source comparison
# ---------------------------------------------------------------------------

def test_per_source_substrate_does_not_underpredict_at_fixed_Mhalo():
    """At fixed halo mass, substrate must predict >= LCDM (boost is positive)."""
    substrate = make_substrate_track()
    lcdm = make_lcdm_track()
    rows = compare_per_source(substrate, lcdm)
    assert len(rows) >= 4
    # Mean log excess should be non-negative (substrate boost adds halos)
    mean = sum(r.log_excess_substrate for r in rows) / len(rows)
    assert mean > -0.05, f"Substrate predicts FEWER than LCDM, mean = {mean}"


def test_per_source_finite_predictions():
    substrate = make_substrate_track()
    lcdm = make_lcdm_track()
    rows = compare_per_source(substrate, lcdm)
    for r in rows:
        assert math.isfinite(r.log_phi_substrate)
        assert math.isfinite(r.log_phi_lcdm)
        assert math.isfinite(r.log_excess_substrate)


# ---------------------------------------------------------------------------
# UV LF comparison
# ---------------------------------------------------------------------------

def test_compare_UV_LF_substrate_closer_than_LCDM():
    """Substrate should have lower |n-sigma| to obs than LCDM at most bins."""
    substrate = make_substrate_track()
    lcdm = make_lcdm_track()
    rows = compare_UV_LF(substrate, lcdm)
    assert len(rows) >= 4
    n_better = sum(1 for r in rows if abs(r.n_sigma_substrate) < abs(r.n_sigma_lcdm))
    # Substrate must be closer than LCDM in the majority of LF bins
    assert n_better >= len(rows) // 2 + 1, (
        f"Substrate not closer than LCDM in majority of LF bins "
        f"({n_better}/{len(rows)})"
    )


# ---------------------------------------------------------------------------
# Stellar mass density (Labbé+2023)
# ---------------------------------------------------------------------------

def test_integrated_density_decreases_with_M_min():
    """Density above 10^11 M_sun must be below density above 10^10 M_sun."""
    substrate = make_substrate_track()
    rho_low = integrated_stellar_density(substrate, z=10.0, log_M_star_min=10.0)
    rho_high = integrated_stellar_density(substrate, z=10.0, log_M_star_min=11.0)
    assert rho_low > rho_high > 0


def test_density_substrate_strictly_above_LCDM():
    """At every Labbé z, substrate density must exceed LCDM density."""
    substrate = make_substrate_track()
    lcdm = make_lcdm_track()
    rows = compare_stellar_mass_density(substrate, lcdm)
    assert len(rows) >= 2
    for r in rows:
        assert r.log_rho_substrate > r.log_rho_lcdm, (
            f"Substrate density NOT greater than LCDM at z={r.z}: "
            f"sub={r.log_rho_substrate}, lcdm={r.log_rho_lcdm}"
        )


# ---------------------------------------------------------------------------
# Top-level result + verdict
# ---------------------------------------------------------------------------

def test_run_returns_complete_result():
    r = run_jwst_high_z_test()
    assert isinstance(r, JWSTHighZTestResult)
    assert r.n_sources >= 4
    assert len(r.sources) == r.n_sources
    assert len(r.UV_LF) >= 4
    assert len(r.density) >= 2
    assert math.isfinite(r.mean_log_excess)
    assert math.isfinite(r.median_log_excess)


def test_total_log_boost_dex_in_observed_jwst_range():
    """Total substrate boost vs LCDM at fixed M_* must be in 0.3 - 1.5 dex.

    JWST observed over-abundance is ~ 0.5 - 1.0 dex (Boylan-Kolchin 2023).
    Substrate prediction must land in this same order-of-magnitude window
    or the framework fails the JWST high-z falsifier.
    """
    r = run_jwst_high_z_test()
    lo, hi = JWST_OBSERVED_OVERABUNDANCE_DEX
    # Allow slack on the high side (some tension is acceptable);
    # require boost > 0.3 dex (factor 2x) on the low side.
    assert 0.30 <= r.total_log_boost_dex <= hi + 0.5, (
        f"Substrate total boost {r.total_log_boost_dex:+.3f} dex "
        f"outside 0.30 - {hi+0.5:+.2f} dex window. "
        f"JWST observed overabundance is {lo:.2f} - {hi:.2f} dex."
    )


def test_eps_star_log_boost_positive():
    r = run_jwst_high_z_test()
    # Substrate eps_star (0.20) must exceed LCDM eps_star (0.06)
    assert r.eps_star_log_boost_dex > 0.0
    assert r.eps_star_substrate > r.eps_star_lcdm


def test_verdict_contains_match_keyword():
    """The verdict must classify the substrate prediction qualitatively."""
    r = run_jwst_high_z_test()
    assert isinstance(r.verdict, str)
    assert any(
        kw in r.verdict for kw in
        ("MATCH", "PARTIAL", "OVER-PREDICTED", "FAIL", "MARGINAL",
         "PATHOLOGICAL")
    )


def test_summary_is_serialisable():
    r = run_jwst_high_z_test()
    s = r.summary()
    for k, v in s.items():
        assert isinstance(v, float)
        assert math.isfinite(v)


def test_report_contains_key_sections():
    rep = report()
    for section in (
        "JWST high-z galaxies",
        "Per-source",
        "UV LF",
        "Cosmic stellar mass density",
        "Verdict",
        "eps_star",
    ):
        assert section in rep, f"Missing section in report: {section!r}"


# ---------------------------------------------------------------------------
# Cross-check against canonical Boylan-Kolchin tension expectation
# ---------------------------------------------------------------------------

def test_substrate_resolves_LCDM_tension_for_majority_LF_bins():
    """The substrate prediction should land closer to obs than LCDM at z>=10."""
    substrate = make_substrate_track()
    lcdm = make_lcdm_track()
    bins = compare_UV_LF(substrate, lcdm)
    sub_score = sum(abs(b.n_sigma_substrate) for b in bins)
    lcdm_score = sum(abs(b.n_sigma_lcdm) for b in bins)
    assert sub_score < lcdm_score, (
        f"Substrate sum |nσ| = {sub_score:.2f} not better than LCDM "
        f"{lcdm_score:.2f}"
    )
