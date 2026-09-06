"""Tests for ``src.stiff_medium.bcs_gap_ratio_test``.

Covers:
    * the substrate / BCS prediction constant 2π/e^γ
    * the dimensional ``measured_ratio`` arithmetic
    * deviation_percent sign and magnitude
    * the literature data table is well-formed and immutable
    * cross-check with the existing ``superconductivity_substrate``
      module so the two stay in sync (same constant)
    * ``run_test`` returns one row per material and respects the
      threshold flag
    * ``summary_stats`` aggregates correctly, including the elemental-
      only restriction
    * known per-material outcomes (BARE BCS):
        - Sn / Ta / In / Tl elemental weak-coupling materials match
          within 5%
        - Pb (strong-coupling Eliashberg) deviates above 20%
        - MgB_2 (multiband sigma+pi) deviates above 15%
    * Allen-Dynes strong-coupling correction:
        - rescues Pb to within 5%
        - improves Hg and Nb dramatically
        - MATERIAL_PHONON_PARAMS has the listed elements
        - predict_with_allen_dynes input validation
        - aggregate AD agreement: at least 7/10 within 5%, mean |dev| < 5%
    * MgB_2 multiband per-band breakdown matches Choi 2002 within 5%
    * the renderer writes a non-empty PNG to disk
"""

from __future__ import annotations

import math
import os
import tempfile

import pytest

from src.stiff_medium.bcs_gap_ratio_test import (
    ALLEN_DYNES_PREFACTOR,
    BCS_RATIO_PRED,
    EULER_GAMMA,
    K_B_MEV_K,
    MATERIAL_PHONON_PARAMS,
    MATERIALS,
    MaterialDatum,
    MgB2Bands,
    ResultRow,
    deviation_percent,
    measured_ratio,
    mgb2_multiband_prediction,
    predict_with_allen_dynes,
    render_bcs_gap_ratio_test,
    run_test,
    summary_stats,
)


# ---------------------------------------------------------------------------
# Constant: substrate / BCS prediction
# ---------------------------------------------------------------------------

def test_prediction_constant_value():
    """2 π / e^γ to 12 digits."""
    assert math.isclose(BCS_RATIO_PRED, 3.527753977724091, rel_tol=1e-12)


def test_prediction_definition_matches_2pi_over_eEgamma():
    """The constant equals 2π / exp(γ) by construction."""
    assert math.isclose(
        BCS_RATIO_PRED, 2.0 * math.pi / math.exp(EULER_GAMMA), rel_tol=1e-15
    )


def test_euler_gamma_value():
    """Euler-Mascheroni γ to 16 digits."""
    assert math.isclose(EULER_GAMMA, 0.5772156649015329, rel_tol=1e-15)


def test_kB_meV_per_K_value():
    """k_B in meV/K = 0.0861733... (textbook)."""
    assert math.isclose(K_B_MEV_K, 0.08617333262145178, rel_tol=1e-9)


def test_prediction_matches_superconductivity_substrate_module():
    """Cross-check: ``superconductivity_substrate`` exposes the same
    universal ratio.  Drift between the two would be a code-quality bug.
    """
    from src.stiff_medium.superconductivity_substrate import (
        BCS_UNIVERSAL_GAP_RATIO,
    )
    assert math.isclose(BCS_RATIO_PRED, BCS_UNIVERSAL_GAP_RATIO,
                        rel_tol=1e-15)


# ---------------------------------------------------------------------------
# Arithmetic: measured_ratio + deviation_percent
# ---------------------------------------------------------------------------

def test_measured_ratio_pb_value():
    """Pb: 2Δ=2.74 meV, T_c=7.20 K -> R = 2.74 / (k_B*7.20) ~ 4.416."""
    R = measured_ratio(7.20, 2.74)
    expected = 2.74 / (K_B_MEV_K * 7.20)
    assert math.isclose(R, expected, rel_tol=1e-12)
    assert math.isclose(R, 4.4156, rel_tol=1e-3)


def test_measured_ratio_dimensional_consistency():
    """If we double the gap, the ratio doubles; if we double T_c,
    the ratio halves.
    """
    R0 = measured_ratio(5.0, 1.0)
    assert math.isclose(measured_ratio(5.0, 2.0), 2.0 * R0, rel_tol=1e-12)
    assert math.isclose(measured_ratio(10.0, 1.0), 0.5 * R0, rel_tol=1e-12)


@pytest.mark.parametrize("Tc, gap", [
    (0.0, 1.0), (-1.0, 1.0), (1.0, 0.0), (1.0, -1.0),
])
def test_measured_ratio_rejects_nonpositive(Tc, gap):
    with pytest.raises(ValueError):
        measured_ratio(Tc, gap)


def test_deviation_percent_zero_when_match():
    assert deviation_percent(BCS_RATIO_PRED) == 0.0


def test_deviation_percent_sign():
    """R > pred -> positive deviation; R < pred -> negative."""
    assert deviation_percent(BCS_RATIO_PRED * 1.10) > 0
    assert deviation_percent(BCS_RATIO_PRED * 0.90) < 0


def test_deviation_percent_magnitude():
    """Exact 10% deviation reported as +10.0."""
    assert math.isclose(
        deviation_percent(BCS_RATIO_PRED * 1.10), 10.0, rel_tol=1e-12
    )


# ---------------------------------------------------------------------------
# MATERIALS table
# ---------------------------------------------------------------------------

def test_materials_table_has_ten_entries():
    assert len(MATERIALS) == 10


def test_materials_includes_expected_names():
    expected = {"Hg", "Pb", "Sn", "Al", "Nb", "V", "Ta", "In", "Tl", "MgB2"}
    actual = {m.name for m in MATERIALS}
    assert actual == expected


def test_materials_are_immutable():
    """MaterialDatum is a frozen dataclass."""
    m = MATERIALS[0]
    with pytest.raises((AttributeError, Exception)):  # FrozenInstanceError
        m.T_c_K = 999.9  # type: ignore[misc]


def test_materials_have_positive_values():
    for m in MATERIALS:
        assert m.T_c_K > 0.0
        assert m.gap_meV > 0.0
        assert m.klass in {"elemental", "multiband"}


# ---------------------------------------------------------------------------
# run_test
# ---------------------------------------------------------------------------

def test_run_test_returns_one_row_per_material():
    rows = run_test()
    assert len(rows) == len(MATERIALS)
    assert [r.name for r in rows] == [m.name for m in MATERIALS]


def test_run_test_within_5pct_flag_consistent():
    rows = run_test(threshold_pct=5.0)
    for r in rows:
        assert r.within_5pct == (abs(r.dev_pct) <= 5.0)


def test_run_test_threshold_strict_vs_loose():
    """Tightening threshold to 1% strictly reduces (or holds) the number
    of materials that pass; loosening to 30% accepts everyone here.
    """
    n_strict = sum(r.within_5pct for r in run_test(threshold_pct=1.0))
    n_default = sum(r.within_5pct for r in run_test(threshold_pct=5.0))
    n_loose = sum(r.within_5pct for r in run_test(threshold_pct=30.0))
    assert n_strict <= n_default <= n_loose
    assert n_loose == len(MATERIALS)


def test_run_test_R_pred_carries_through():
    """The R_pred field on each row equals the prediction passed in."""
    rows = run_test(R_pred=4.0)
    for r in rows:
        assert r.R_pred == 4.0


# ---------------------------------------------------------------------------
# Per-material expected verdicts
# ---------------------------------------------------------------------------

def _row(name: str) -> ResultRow:
    rows = run_test()
    for r in rows:
        if r.name == name:
            return r
    raise KeyError(name)


@pytest.mark.parametrize("name", ["Sn", "Ta", "In", "Tl"])
def test_weak_coupling_elementals_match_within_5pct(name):
    """Sn, Ta, In, Tl are textbook weak-coupling and should hit 3.528
    within 5%.
    """
    r = _row(name)
    assert r.within_5pct, (
        f"{name}: |dev|={abs(r.dev_pct):.2f}% > 5% threshold"
    )


def test_pb_strong_coupling_deviates_above_20pct():
    """Pb is the canonical strong-coupling material (Carbotte review,
    ratio ~4.4); should deviate by > 20%.
    """
    r = _row("Pb")
    assert r.dev_pct > 20.0, (
        f"Pb deviation {r.dev_pct:+.2f}% expected > +20% (strong coupling)"
    )


def test_mgb2_multiband_deviates_above_15pct():
    """MgB_2 is two-band sigma+pi; using the dominant sigma gap
    (~14 meV) gives > 15% above the single-band BCS line.
    """
    r = _row("MgB2")
    assert r.dev_pct > 15.0, (
        f"MgB2 deviation {r.dev_pct:+.2f}% expected > +15% (multiband)"
    )


def test_hg_mildly_strong_coupling():
    """Hg deviates above 5% but well below Pb (mildly strong-coupling)."""
    r = _row("Hg")
    assert r.dev_pct > 5.0
    assert r.dev_pct < _row("Pb").dev_pct


# ---------------------------------------------------------------------------
# summary_stats
# ---------------------------------------------------------------------------

def test_summary_stats_keys():
    stats = summary_stats(run_test())
    for key in ("n", "n_within_5pct", "mean_abs_dev",
                "max_abs_dev", "min_abs_dev", "elemental_only"):
        assert key in stats
    for key in ("n", "n_within_5pct", "mean_abs_dev",
                "max_abs_dev", "min_abs_dev"):
        assert key in stats["elemental_only"]


def test_summary_stats_consistency():
    rows = run_test()
    stats = summary_stats(rows)
    assert stats["n"] == len(rows)
    assert stats["n_within_5pct"] == sum(r.within_5pct for r in rows)
    assert stats["min_abs_dev"] <= stats["mean_abs_dev"] <= stats["max_abs_dev"]


def test_summary_stats_elemental_subset():
    rows = run_test()
    stats = summary_stats(rows)
    elem_count = sum(1 for r in rows if r.klass == "elemental")
    assert stats["elemental_only"]["n"] == elem_count
    assert stats["elemental_only"]["n"] == 9
    assert stats["n"] - stats["elemental_only"]["n"] == 1  # MgB2


def test_summary_stats_empty_raises():
    with pytest.raises(ValueError):
        summary_stats([])


# ---------------------------------------------------------------------------
# Allen-Dynes strong-coupling correction
# ---------------------------------------------------------------------------

def test_allen_dynes_prefactor_value():
    """The Allen-Dynes leading coefficient is 12.5 (Allen-Dynes 1975)."""
    assert ALLEN_DYNES_PREFACTOR == 12.5


def test_predict_with_allen_dynes_reduces_to_bcs_in_weak_coupling():
    """In the weak-coupling limit T_c << omega_log the correction
    vanishes (12.5 * (T_c/wlog)^2 * ln(wlog/(2 T_c)) -> 0).
    """
    R = predict_with_allen_dynes(lam=0.4, omega_log_K=10000.0, T_c_K=1.0)
    assert math.isclose(R, BCS_RATIO_PRED, rel_tol=1e-3)
    assert R > BCS_RATIO_PRED  # any T_c < wlog/2 gives positive correction


def test_predict_with_allen_dynes_pb():
    """Pb: lambda=1.55, omega_log=56 K, T_c=7.20 K
    -> AD ratio about 4.5 (matches Carbotte 1990 review for Pb).
    """
    lam, wlog = MATERIAL_PHONON_PARAMS["Pb"]
    R = predict_with_allen_dynes(lam, wlog, T_c_K=7.20)
    # Carbotte 1990 strong-coupling Pb ratio: ~ 4.3-4.5
    assert 4.3 < R < 4.6, f"Pb AD prediction {R:.3f} outside Carbotte band"


def test_predict_with_allen_dynes_returns_weak_when_out_of_range():
    """When omega_log <= 2 T_c, the leading-order AD form breaks
    down; the function falls back to the bare weak-coupling line.
    """
    R = predict_with_allen_dynes(lam=1.5, omega_log_K=10.0, T_c_K=10.0)
    assert R == BCS_RATIO_PRED
    R = predict_with_allen_dynes(lam=1.5, omega_log_K=15.0, T_c_K=10.0)
    assert R == BCS_RATIO_PRED


def test_predict_with_allen_dynes_monotone_in_Tc_over_wlog():
    """Holding lambda and omega_log fixed and increasing T_c (still <
    omega_log/2) should monotonically increase the AD-corrected ratio,
    because both (T_c/wlog)^2 grows and ln(wlog/(2 T_c)) is positive.
    """
    lam, wlog = 1.0, 200.0
    R_lo = predict_with_allen_dynes(lam, wlog, T_c_K=5.0)
    R_md = predict_with_allen_dynes(lam, wlog, T_c_K=20.0)
    R_hi = predict_with_allen_dynes(lam, wlog, T_c_K=50.0)
    assert R_lo < R_md < R_hi


@pytest.mark.parametrize("Tc, wlog, lam", [
    (-1.0, 100.0, 1.0),
    (0.0,  100.0, 1.0),
    (10.0,  -1.0, 1.0),
    (10.0,   0.0, 1.0),
    (10.0, 100.0, -0.1),
])
def test_predict_with_allen_dynes_input_validation(Tc, wlog, lam):
    with pytest.raises(ValueError):
        predict_with_allen_dynes(lam, wlog, Tc)


def test_phonon_params_includes_all_elementals():
    """Every elemental in the literature table has a (lambda, omega_log)
    entry in ``MATERIAL_PHONON_PARAMS``.  MgB_2 is also listed
    (single-band sigma fallback) but is handled multiband-correctly by
    ``mgb2_multiband_prediction``.
    """
    expected = {"Hg", "Pb", "Sn", "Al", "Nb", "V", "Ta", "In", "Tl", "MgB2"}
    assert set(MATERIAL_PHONON_PARAMS.keys()) == expected


def test_phonon_params_positive():
    for name, (lam, wlog) in MATERIAL_PHONON_PARAMS.items():
        assert lam > 0.0,  f"{name}: lambda must be positive"
        assert wlog > 0.0, f"{name}: omega_log must be positive"


def test_phonon_params_pb_strongest_coupling():
    """Among elementals Pb and Hg are known to be the most strongly
    coupled (lambda ~ 1.5+); Al and V are the weakest (lambda < 0.7).
    """
    pb_lam = MATERIAL_PHONON_PARAMS["Pb"][0]
    hg_lam = MATERIAL_PHONON_PARAMS["Hg"][0]
    al_lam = MATERIAL_PHONON_PARAMS["Al"][0]
    v_lam  = MATERIAL_PHONON_PARAMS["V"][0]
    assert pb_lam > 1.0
    assert hg_lam > 1.0
    assert al_lam < 0.7
    assert v_lam  < 0.7


# ---------------------------------------------------------------------------
# Per-material AD verdicts
# ---------------------------------------------------------------------------

def test_ad_correction_rescues_pb_within_5pct():
    """Pb deviates +25% from bare BCS but Allen-Dynes-corrected
    prediction agrees with measurement within 5%.
    """
    r = _row("Pb")
    assert abs(r.dev_pct_ad) < 5.0, (
        f"Pb AD-corrected deviation {r.dev_pct_ad:+.2f}% exceeds 5%"
    )
    # And the measured value is within 5% of the AD prediction (4.41 vs ~4.5)
    assert abs(r.R_meas - 4.416) < 0.01     # sanity on R_meas
    assert abs(r.R_pred_ad - r.R_meas) / r.R_meas * 100 < 5.0


def test_ad_correction_improves_hg_dramatically():
    """Hg bare deviation +12% drops to within 5% under AD."""
    r = _row("Hg")
    assert abs(r.dev_pct_ad) < 5.0
    assert abs(r.dev_pct_ad) < abs(r.dev_pct)


def test_ad_correction_improves_nb():
    """Nb bare +8.5% drops below 5% with AD (lambda=0.82)."""
    r = _row("Nb")
    assert abs(r.dev_pct_ad) < 5.0
    assert abs(r.dev_pct_ad) < abs(r.dev_pct)


def test_ad_does_not_help_weak_coupling_already_passing():
    """Sn / Ta / In / Tl already match within 5% bare; AD changes
    them only marginally (still within 5%).
    """
    for name in ("Sn", "Ta", "In", "Tl"):
        r = _row(name)
        assert abs(r.dev_pct_ad) < 5.0, (
            f"{name} AD-corrected exits 5% band: {r.dev_pct_ad:+.2f}%"
        )


def test_aggregate_ad_majority_within_5pct():
    """With Allen-Dynes corrections at least 7/10 materials sit
    within the 5 % band (vs 4/10 bare).  This is the headline
    'corrections improve agreement' claim of the test, but it is not
    9/10: Al and V are below the BCS line by ~6-8% which is NOT
    explained by strong coupling (signs go the wrong way), and the
    single-band AD formula does not fully account for the multiband
    sigma-channel rise in MgB_2 (it leaves +10% residual).
    """
    rows = run_test()
    n_bare = sum(1 for r in rows if r.within_5pct)
    n_ad   = sum(1 for r in rows if r.within_5pct_ad)
    assert n_bare == 4
    assert n_ad >= 7
    assert n_ad > n_bare


def test_aggregate_ad_mean_dev_below_5pct():
    """The Allen-Dynes-corrected mean |deviation| across all 10
    materials drops below 5% (from 8.4% bare).
    """
    stats = summary_stats(run_test())
    assert stats["mean_abs_dev"]    > 8.0     # bare ~ 8.4%
    assert stats["mean_abs_dev_ad"] < 5.0     # AD   ~ 3.8%
    assert stats["mean_abs_dev_ad"] < stats["mean_abs_dev"]


def test_aggregate_ad_max_dev_dropped_dramatically():
    """The Allen-Dynes correction drops the worst-case deviation
    from ~25% (Pb bare) to ~10% (MgB_2 single-band-AD residual).
    """
    stats = summary_stats(run_test())
    assert stats["max_abs_dev"]    > 20.0    # bare ~ 25.2% (Pb)
    assert stats["max_abs_dev_ad"] < 11.0    # AD   ~ 10.1% (MgB_2)


# ---------------------------------------------------------------------------
# MgB_2 multiband module
# ---------------------------------------------------------------------------

def test_mgb2_multiband_keys_present():
    out = mgb2_multiband_prediction()
    for key in (
        "R_sigma", "R_pi", "R_composite", "R_weighted_substrate",
        "R_pred_weak", "R_pred_dominant_AD",
        "dev_sigma_pct", "dev_composite_pct",
        "Choi_R_sigma_ref", "Choi_R_pi_ref",
    ):
        assert key in out, f"missing key {key!r}"


def test_mgb2_sigma_band_matches_choi_within_5pct():
    """Choi et al. 2002 Eliashberg gives R_sigma ~ 4.0; the substrate
    table uses 2 Delta_sigma = 14 meV, T_c = 39 K -> R_sigma = 4.17.
    Agreement with Choi within 5%.
    """
    out = mgb2_multiband_prediction()
    rel = abs(out["R_sigma"] - out["Choi_R_sigma_ref"]) / out["Choi_R_sigma_ref"]
    assert rel < 0.05, (
        f"R_sigma={out['R_sigma']:.3f} vs Choi {out['Choi_R_sigma_ref']:.2f} "
        f"({rel*100:.2f}% off)"
    )


def test_mgb2_pi_band_matches_choi_within_5pct():
    """Choi et al. 2002 Eliashberg gives R_pi ~ 1.2; substrate table
    uses 2 Delta_pi = 4 meV, T_c = 39 K -> R_pi = 1.19. Agreement
    with Choi within 5%.
    """
    out = mgb2_multiband_prediction()
    rel = abs(out["R_pi"] - out["Choi_R_pi_ref"]) / out["Choi_R_pi_ref"]
    assert rel < 0.05, (
        f"R_pi={out['R_pi']:.3f} vs Choi {out['Choi_R_pi_ref']:.2f} "
        f"({rel*100:.2f}% off)"
    )


def test_mgb2_dominant_ad_within_10pct_of_measured_sigma():
    """The single-band Allen-Dynes correction applied to the dominant
    sigma channel (lambda_sigma=0.96, omega_log_sigma~778 K) should
    bring the measured sigma-band ratio (4.17) within 10% --- still
    not within 5% because the multiband two-gap structure raises
    R_sigma above any single-band line.  This is the known-residual
    multiband signature.
    """
    out = mgb2_multiband_prediction()
    assert abs(out["dev_sigma_pct"]) < 11.0


def test_mgb2_bands_dataclass_defaults_match_literature():
    bands = MgB2Bands()
    assert bands.Tc_K == 39.0
    assert bands.Delta_sigma_meV == 7.0
    assert bands.Delta_pi_meV    == 2.0
    # full gaps
    assert 2.0 * bands.Delta_sigma_meV == 14.0
    assert 2.0 * bands.Delta_pi_meV    == 4.0
    # weights sum to 1
    assert 0.0 < bands.w_sigma < 1.0


def test_mgb2_weighted_substrate_equals_R_composite():
    """Algebraic identity: w_sigma * R_sigma + (1 - w_sigma) * R_pi
    equals R_composite (gap weighted average over same T_c).
    """
    out = mgb2_multiband_prediction()
    assert math.isclose(
        out["R_weighted_substrate"], out["R_composite"], rel_tol=1e-12
    )


# ---------------------------------------------------------------------------
# ResultRow / run_test wiring (new fields)
# ---------------------------------------------------------------------------

def test_resultrow_has_ad_fields():
    """ResultRow exposes R_pred_ad, dev_pct_ad, within_5pct_ad."""
    rows = run_test()
    r = rows[0]
    assert hasattr(r, "R_pred_ad")
    assert hasattr(r, "dev_pct_ad")
    assert hasattr(r, "within_5pct_ad")


def test_resultrow_within_5pct_ad_consistent_with_dev():
    rows = run_test(threshold_pct=5.0)
    for r in rows:
        assert r.within_5pct_ad == (abs(r.dev_pct_ad) <= 5.0)


def test_resultrow_R_pred_ad_equals_predict_function():
    """Each row's R_pred_ad should equal predict_with_allen_dynes
    called with the literature (lambda, omega_log) for that material.
    """
    rows = run_test()
    for r in rows:
        if r.name not in MATERIAL_PHONON_PARAMS:
            assert r.R_pred_ad == r.R_pred
            continue
        lam, wlog = MATERIAL_PHONON_PARAMS[r.name]
        expected = predict_with_allen_dynes(lam, wlog, r.T_c_K)
        assert math.isclose(r.R_pred_ad, expected, rel_tol=1e-12)


def test_summary_stats_has_ad_fields():
    stats = summary_stats(run_test())
    for key in ("n_within_5pct_ad", "mean_abs_dev_ad",
                "max_abs_dev_ad", "min_abs_dev_ad"):
        assert key in stats
    elem = stats["elemental_only"]
    for key in ("n_within_5pct_ad", "mean_abs_dev_ad",
                "max_abs_dev_ad", "min_abs_dev_ad"):
        assert key in elem


# ---------------------------------------------------------------------------
# Renderer
# ---------------------------------------------------------------------------

def test_render_bcs_gap_ratio_test_writes_png():
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "120_bcs_ratio_test.png")
        path = render_bcs_gap_ratio_test(out)
        assert path == out
        assert os.path.isfile(path)
        assert os.path.getsize(path) > 5_000  # non-trivial PNG
        # PNG signature: 0x89 'P' 'N' 'G'
        with open(path, "rb") as fh:
            head = fh.read(8)
        assert head[:4] == b"\x89PNG"
