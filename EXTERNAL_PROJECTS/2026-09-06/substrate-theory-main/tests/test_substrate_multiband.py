"""Tests for ``src.stiff_medium.substrate_multiband``.

Covers:
    * eps_unit = Lambda_QCD / n_R = 11.111 meV (the substrate
      band-gap unit, identical to the high-T_c bound generator)
    * gap_sigma_band derives Delta_sigma = 6.667 meV from K_4 face
      sp^2 loading 3 / K_rank
    * gap_pi_band   derives Delta_pi    = 1.852 meV from Q_3 axis
      Koide loading 1 / (F * R)
    * ratio_sigma_pi closed-form 3 * F * R / K_rank = 18 / 5 = 3.6
    * Both gaps within 10 % of measured values (7.0 and 2.0 meV)
    * Per-band BCS ratios within 10 % of Choi 2002 (4.0 and 1.2)
    * Algebraic identities (substrate prediction independent of which
      Tc anchor we feed for the BCS ratio computation; ratios scale
      inversely with Tc)
    * Bundle prediction has all expected keys, scalar types
    * MgB2SubstrateBands is frozen and reads canonical b3_constants
    * predict_with custom integers: shifting K_rank or n_R perturbs
      the predictions (rigidity check)
    * Cross-module consistency with bcs_gap_ratio_test
      mgb2_multiband_prediction (substrate keys appear and match
      direct calls)
"""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.b3_constants import (
    F,
    K_rank,
    LAMBDA_QCD_MEV,
    R,
    n_R,
)
from src.stiff_medium.substrate_multiband import (
    K_B_MEV_K,
    MGB2_DELTA_PI_MEV_EXP,
    MGB2_DELTA_SIGMA_MEV_EXP,
    MGB2_R_PI_REF,
    MGB2_R_SIGMA_REF,
    MGB2_TC_K_REF,
    MgB2SubstrateBands,
    bcs_band_ratio,
    eps_unit_meV,
    gap_pi_band,
    gap_sigma_band,
    predict_mgb2_two_gap,
    ratio_sigma_pi,
    substrate_multiband_summary,
)


# ---------------------------------------------------------------------------
# eps_unit = Lambda_QCD / n_R
# ---------------------------------------------------------------------------

def test_eps_unit_value():
    """eps_unit = 200 MeV / 18 = 11.111 meV."""
    assert math.isclose(eps_unit_meV(), 200.0 / 18.0, rel_tol=1e-12)
    assert math.isclose(eps_unit_meV(), 11.11111111111111, rel_tol=1e-9)


def test_eps_unit_matches_high_tc_ceiling():
    """eps_unit / k_B = T_c,max ~ 128.9 K (same ceiling as the
    high-T_c bound module).  This anchors the substrate-saturation
    cross-disciplinary tie-in.
    """
    Tc_max_K = eps_unit_meV() / K_B_MEV_K
    assert math.isclose(Tc_max_K, 128.94, rel_tol=1e-3)


# ---------------------------------------------------------------------------
# Sigma band: K_4 face-pair sp^2 loading
# ---------------------------------------------------------------------------

def test_gap_sigma_value():
    """Delta_sigma = (200/18) * (3/5) = 6.6667 meV."""
    assert math.isclose(gap_sigma_band(), (200.0 / 18.0) * (3.0 / 5.0),
                        rel_tol=1e-12)
    assert math.isclose(gap_sigma_band(), 6.666666666666666, rel_tol=1e-9)


def test_gap_sigma_within_10pct_of_measured():
    """Substrate prediction 6.667 meV vs measured 7.0 meV (4.8 % under)."""
    pred = gap_sigma_band()
    rel = abs(pred - MGB2_DELTA_SIGMA_MEV_EXP) / MGB2_DELTA_SIGMA_MEV_EXP
    assert rel < 0.10, f"Delta_sigma off by {rel*100:.2f}% (>10%)"


def test_gap_sigma_within_5pct_of_measured():
    """Substrate prediction within 5 % of measured 7 meV (4.76 %)."""
    pred = gap_sigma_band()
    rel = abs(pred - MGB2_DELTA_SIGMA_MEV_EXP) / MGB2_DELTA_SIGMA_MEV_EXP
    assert rel < 0.05, f"Delta_sigma off by {rel*100:.2f}% (>5%)"


def test_gap_sigma_uses_K_rank():
    """Shifting K_rank from 5 to 6 must shift the substrate sigma gap."""
    base = MgB2SubstrateBands()
    perturbed = MgB2SubstrateBands(
        K_rank=6,
        # keep all other integers at canonical values
    )
    assert gap_sigma_band(perturbed) != gap_sigma_band(base)
    # Specifically, K_rank goes up -> denominator up -> gap drops
    assert gap_sigma_band(perturbed) < gap_sigma_band(base)


def test_gap_sigma_uses_n_R():
    """Shifting n_R from 18 to 19 must shift the substrate sigma gap
    (eps_unit denominator).
    """
    base = MgB2SubstrateBands()
    perturbed = MgB2SubstrateBands(n_R=19)
    assert gap_sigma_band(perturbed) < gap_sigma_band(base)


# ---------------------------------------------------------------------------
# Pi band: Q_3 axis Koide loading
# ---------------------------------------------------------------------------

def test_gap_pi_value():
    """Delta_pi = (200/18) * (1 / (2*3)) = 1.8519 meV."""
    assert math.isclose(gap_pi_band(), (200.0 / 18.0) / 6.0, rel_tol=1e-12)
    assert math.isclose(gap_pi_band(), 1.8518518518518519, rel_tol=1e-9)


def test_gap_pi_within_10pct_of_measured():
    """Substrate prediction 1.852 meV vs measured 2.0 meV (7.4 % under)."""
    pred = gap_pi_band()
    rel = abs(pred - MGB2_DELTA_PI_MEV_EXP) / MGB2_DELTA_PI_MEV_EXP
    assert rel < 0.10, f"Delta_pi off by {rel*100:.2f}% (>10%)"


def test_gap_pi_uses_F_and_R():
    """Shifting F or R perturbs the substrate pi gap (Koide F*R)."""
    base = MgB2SubstrateBands()
    pf = MgB2SubstrateBands(F=3)        # F up -> pi gap shrinks
    pr = MgB2SubstrateBands(R=4)        # R up -> pi gap shrinks
    assert gap_pi_band(pf) < gap_pi_band(base)
    assert gap_pi_band(pr) < gap_pi_band(base)


# ---------------------------------------------------------------------------
# Sigma / pi ratio
# ---------------------------------------------------------------------------

def test_ratio_sigma_pi_closed_form():
    """ratio = 3 * F * R / K_rank = 18 / 5 = 3.6."""
    assert math.isclose(ratio_sigma_pi(), 18.0 / 5.0, rel_tol=1e-12)
    assert math.isclose(ratio_sigma_pi(), 3.6, rel_tol=1e-12)


def test_ratio_sigma_pi_matches_explicit_division():
    """Algebraic identity check: ratio = gap_sigma / gap_pi."""
    assert math.isclose(
        ratio_sigma_pi(),
        gap_sigma_band() / gap_pi_band(),
        rel_tol=1e-12,
    )


def test_ratio_sigma_pi_within_5pct_of_measured():
    """Measured 7/2 = 3.5; substrate 3.6; 2.86 % over."""
    measured = MGB2_DELTA_SIGMA_MEV_EXP / MGB2_DELTA_PI_MEV_EXP
    rel = abs(ratio_sigma_pi() - measured) / measured
    assert rel < 0.05, f"ratio off by {rel*100:.2f}% (>5%)"


# ---------------------------------------------------------------------------
# BCS ratios at T_c = 39 K
# ---------------------------------------------------------------------------

def test_bcs_ratio_sigma_within_10pct_of_choi():
    """R_sigma_substrate ~ 3.97 vs Choi 4.0 (0.8 % under)."""
    R_sigma = bcs_band_ratio(gap_sigma_band(), MGB2_TC_K_REF)
    rel = abs(R_sigma - MGB2_R_SIGMA_REF) / MGB2_R_SIGMA_REF
    assert rel < 0.10, f"R_sigma off by {rel*100:.2f}% (>10%)"


def test_bcs_ratio_sigma_within_2pct_of_choi():
    """Closer than 2 % to Choi (0.82 % actually)."""
    R_sigma = bcs_band_ratio(gap_sigma_band(), MGB2_TC_K_REF)
    rel = abs(R_sigma - MGB2_R_SIGMA_REF) / MGB2_R_SIGMA_REF
    assert rel < 0.02


def test_bcs_ratio_pi_within_10pct_of_choi():
    """R_pi_substrate ~ 1.10 vs Choi 1.2 (8.16 % under)."""
    R_pi = bcs_band_ratio(gap_pi_band(), MGB2_TC_K_REF)
    rel = abs(R_pi - MGB2_R_PI_REF) / MGB2_R_PI_REF
    assert rel < 0.10, f"R_pi off by {rel*100:.2f}% (>10%)"


def test_bcs_band_ratio_input_validation():
    with pytest.raises(ValueError):
        bcs_band_ratio(0.0, 39.0)
    with pytest.raises(ValueError):
        bcs_band_ratio(7.0, 0.0)
    with pytest.raises(ValueError):
        bcs_band_ratio(-1.0, 39.0)
    with pytest.raises(ValueError):
        bcs_band_ratio(7.0, -1.0)


def test_bcs_band_ratio_dimensional_consistency():
    """Doubling delta doubles R; doubling Tc halves R."""
    R0 = bcs_band_ratio(7.0, 39.0)
    assert math.isclose(bcs_band_ratio(14.0, 39.0), 2.0 * R0, rel_tol=1e-12)
    assert math.isclose(bcs_band_ratio(7.0, 78.0), 0.5 * R0, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Bundle prediction
# ---------------------------------------------------------------------------

def test_predict_mgb2_two_gap_keys():
    pred = predict_mgb2_two_gap()
    for key in (
        "Delta_sigma_MeV", "Delta_pi_MeV", "eps_unit_MeV",
        "R_sigma_pred", "R_pi_pred",
        "R_sigma_dev_pct", "R_pi_dev_pct",
        "Delta_sigma_dev_pct", "Delta_pi_dev_pct",
        "ratio_sigma_pi", "ratio_sigma_pi_exp",
        "ratio_sigma_pi_dev_pct",
    ):
        assert key in pred


def test_predict_mgb2_two_gap_values_match_individual_calls():
    pred = predict_mgb2_two_gap()
    bands = MgB2SubstrateBands()
    assert math.isclose(pred["Delta_sigma_MeV"], gap_sigma_band(bands),
                        rel_tol=1e-12)
    assert math.isclose(pred["Delta_pi_MeV"], gap_pi_band(bands),
                        rel_tol=1e-12)
    assert math.isclose(pred["eps_unit_MeV"], eps_unit_meV(bands),
                        rel_tol=1e-12)
    assert math.isclose(pred["ratio_sigma_pi"], ratio_sigma_pi(bands),
                        rel_tol=1e-12)


def test_predict_mgb2_two_gap_deviations_in_band():
    pred = predict_mgb2_two_gap()
    # All four headline deviations must sit within 10%
    assert abs(pred["Delta_sigma_dev_pct"]) < 10.0
    assert abs(pred["Delta_pi_dev_pct"])    < 10.0
    assert abs(pred["R_sigma_dev_pct"])     < 10.0
    assert abs(pred["R_pi_dev_pct"])        < 10.0
    assert abs(pred["ratio_sigma_pi_dev_pct"]) < 5.0


# ---------------------------------------------------------------------------
# Frozen dataclass + canonical defaults
# ---------------------------------------------------------------------------

def test_mgb2_substrate_bands_defaults():
    bands = MgB2SubstrateBands()
    assert bands.Lambda_QCD_MeV == LAMBDA_QCD_MEV == 200.0
    assert bands.n_R    == n_R    == 18
    assert bands.K_rank == K_rank == 5
    assert bands.F      == F      == 2
    assert bands.R      == R      == 3
    assert bands.sp2_bond_count == 3
    assert bands.pz_axis_count  == 1
    assert bands.Tc_K   == 39.0


def test_mgb2_substrate_bands_frozen():
    bands = MgB2SubstrateBands()
    with pytest.raises((AttributeError, Exception)):  # FrozenInstanceError
        bands.Lambda_QCD_MeV = 999.0  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Summary + integer rigidity
# ---------------------------------------------------------------------------

def test_substrate_multiband_summary_structure():
    s = substrate_multiband_summary()
    for key in ("anchors_and_integers", "experimental_anchors",
                "prediction", "headline_close_to_measurement_pct"):
        assert key in s
    head = s["headline_close_to_measurement_pct"]
    for key in ("Delta_sigma", "Delta_pi", "R_sigma", "R_pi", "ratio"):
        assert key in head


def test_summary_headlines_stay_under_10pct():
    s = substrate_multiband_summary()
    for key, val in s["headline_close_to_measurement_pct"].items():
        assert abs(val) < 10.0, (
            f"summary headline {key!r} = {val:+.2f}% breaches 10 %"
        )


def test_perturbing_n_R_changes_eps_unit_and_both_gaps():
    """Rigidity: a unit shift in n_R moves both gaps coherently."""
    base = MgB2SubstrateBands()
    bumped = MgB2SubstrateBands(n_R=20)
    assert eps_unit_meV(bumped) < eps_unit_meV(base)
    assert gap_sigma_band(bumped) < gap_sigma_band(base)
    assert gap_pi_band(bumped) < gap_pi_band(base)
    # Ratio is independent of n_R
    assert math.isclose(ratio_sigma_pi(bumped),
                        ratio_sigma_pi(base), rel_tol=1e-12)


def test_perturbing_lambda_qcd_changes_both_gaps_proportionally():
    base = MgB2SubstrateBands()
    bumped = MgB2SubstrateBands(Lambda_QCD_MeV=400.0)
    # Both gaps should double
    assert math.isclose(gap_sigma_band(bumped),
                        2.0 * gap_sigma_band(base), rel_tol=1e-12)
    assert math.isclose(gap_pi_band(bumped),
                        2.0 * gap_pi_band(base), rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Cross-module integration with bcs_gap_ratio_test
# ---------------------------------------------------------------------------

def test_bcs_gap_ratio_test_module_exposes_substrate_predictions():
    """``bcs_gap_ratio_test.mgb2_multiband_prediction`` should export
    the substrate-derived predictions (K_4 / Q_3 channel) so the
    substrate prediction is queryable from the existing API.
    """
    from src.stiff_medium.bcs_gap_ratio_test import mgb2_multiband_prediction
    out = mgb2_multiband_prediction()
    for key in (
        "Delta_sigma_substrate_pred_meV",
        "Delta_pi_substrate_pred_meV",
        "R_sigma_substrate_pred",
        "R_pi_substrate_pred",
        "dev_sigma_substrate_pct",
        "dev_pi_substrate_pct",
        "ratio_sigma_pi_substrate",
    ):
        assert key in out, f"missing substrate key {key!r}"


def test_bcs_gap_ratio_test_substrate_values_match_substrate_module():
    """Cross-check: bcs_gap_ratio_test substrate predictions are
    numerically identical to direct calls into substrate_multiband.
    """
    from src.stiff_medium.bcs_gap_ratio_test import mgb2_multiband_prediction
    out = mgb2_multiband_prediction()
    assert math.isclose(out["Delta_sigma_substrate_pred_meV"],
                        gap_sigma_band(), rel_tol=1e-12)
    assert math.isclose(out["Delta_pi_substrate_pred_meV"],
                        gap_pi_band(), rel_tol=1e-12)
    assert math.isclose(out["ratio_sigma_pi_substrate"],
                        ratio_sigma_pi(), rel_tol=1e-12)
