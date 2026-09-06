"""Tests for src/stiff_medium/debye_test.py.

Required claims:
    (a) Schema integrity:    every row in run_test() carries the expected
        keys with positive numeric values.
    (b) Sound-speed identities: Debye-averaged c_s satisfies the textbook
        identity 3/c_D³ = 1/c_L³ + 2/c_T³ to machine precision.
    (c) Atomic density formula: n_atoms = N_A · ρ / M to 1e-12 relative.
    (d) Anchor matches: diamond, copper, aluminum predictions reproduce the
        baseline values quoted in user_context (2228, 343, ~430 K) at the
        few-percent level.
    (e) Scaling: predicted Θ_D scales as expected with c_s and n.
    (f) Aggregate fit quality: at least 10/14 materials within 10%, and
        log-log Pearson r > 0.99 over the full set.
    (g) render_debye_test produces a non-empty PNG at the canonical path.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pytest

from src.stiff_medium.debye_test import (
    HBAR,
    KB,
    MATERIALS,
    Material,
    N_A,
    SIGMA_MAX_POISSON,
    debye_average_speed,
    debye_temperature,
    gamma_G_substrate,
    longitudinal_speed,
    material_sound_speeds,
    predict_theta_D,
    predict_with_quasiharmonic,
    quasi_harmonic_shift,
    render_debye_test,
    run_test,
    run_test_quasiharmonic,
    transverse_speed,
)


# --------------------------------------------------------------------------- #
# Claim (a): schema integrity                                                  #
# --------------------------------------------------------------------------- #


def test_run_test_schema():
    res = run_test()
    assert "rows" in res and "summary" in res
    rows = res["rows"]
    assert set(rows.keys()) == set(MATERIALS.keys())
    expected_keys = {
        "name", "lattice", "rho_mass", "B_GPa", "G_GPa", "n_atoms",
        "c_L", "c_T", "c_Debye",
        "theta_pred", "theta_meas", "rel_err",
    }
    for name, row in rows.items():
        assert set(row.keys()) == expected_keys, f"{name}: schema mismatch"
        # All numeric scalars (except 'name' and 'lattice') must be positive
        for k, v in row.items():
            if k in ("name", "lattice"):
                continue
            if k == "rel_err":
                continue   # may be slightly negative
            assert v > 0.0, f"{name}: {k} not positive ({v})"


def test_summary_contains_expected_keys():
    res = run_test()
    s = res["summary"]
    expected = {
        "n_materials", "mean_rel_err", "mean_abs_rel_err",
        "median_abs_rel_err", "max_abs_rel_err",
        "max_rel_err_material", "rms_rel_err",
        "pearson_r", "pearson_log_r",
        "loglog_slope", "loglog_intercept",
        "within_5pct", "within_10pct", "within_20pct",
    }
    assert expected.issubset(set(s.keys()))


# --------------------------------------------------------------------------- #
# Claim (b): Debye-averaging identity                                          #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", list(MATERIALS.keys()))
def test_debye_average_identity(name):
    """3/c_D³ must equal 1/c_L³ + 2/c_T³ to machine precision."""
    mat = MATERIALS[name]
    c_L, c_T, c_D = material_sound_speeds(mat)
    lhs = 3.0 / c_D ** 3
    rhs = 1.0 / c_L ** 3 + 2.0 / c_T ** 3
    assert math.isclose(lhs, rhs, rel_tol=1.0e-12), (
        f"{name}: Debye averaging broken (lhs={lhs}, rhs={rhs})"
    )


def test_longitudinal_faster_than_transverse():
    """Physically, c_L > c_T always (since 4G/3 > 0)."""
    for name, mat in MATERIALS.items():
        c_L, c_T, _ = material_sound_speeds(mat)
        assert c_L > c_T, f"{name}: c_L ({c_L}) <= c_T ({c_T})"


def test_debye_average_between_T_and_L():
    """c_T < c_D < c_L for every material."""
    for name, mat in MATERIALS.items():
        c_L, c_T, c_D = material_sound_speeds(mat)
        assert c_T < c_D < c_L, (
            f"{name}: c_D not between c_T and c_L "
            f"(c_T={c_T}, c_D={c_D}, c_L={c_L})"
        )


# --------------------------------------------------------------------------- #
# Claim (c): atomic-density identity                                           #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", list(MATERIALS.keys()))
def test_atomic_density_formula(name):
    mat = MATERIALS[name]
    expected = N_A * mat.rho_mass / mat.M_molar
    assert math.isclose(mat.n_atoms_per_m3, expected, rel_tol=1.0e-12), (
        f"{name}: n = N_A·ρ/M broken"
    )


# --------------------------------------------------------------------------- #
# Claim (d): canonical anchor materials match user_context targets             #
# --------------------------------------------------------------------------- #


def test_diamond_within_5pct():
    """Diamond Θ_D should match measured 2230 K within 5%."""
    row = predict_theta_D(MATERIALS["Diamond"])
    assert abs(row["rel_err"]) < 0.05, (
        f"Diamond Θ_D = {row['theta_pred']:.1f} K, measured 2230, "
        f"rel err {row['rel_err']:+.2%}"
    )


def test_copper_within_2pct():
    """Copper Θ_D should match measured 343 K within 2% (anchor calibration)."""
    row = predict_theta_D(MATERIALS["Copper"])
    assert abs(row["rel_err"]) < 0.02, (
        f"Copper Θ_D = {row['theta_pred']:.1f} K, measured 343, "
        f"rel err {row['rel_err']:+.2%}"
    )


def test_aluminum_within_5pct():
    """Aluminum Θ_D should match measured 428 K within 5%."""
    row = predict_theta_D(MATERIALS["Aluminum"])
    assert abs(row["rel_err"]) < 0.05, (
        f"Aluminum Θ_D = {row['theta_pred']:.1f} K, measured 428, "
        f"rel err {row['rel_err']:+.2%}"
    )


def test_silicon_within_2pct():
    """Silicon Θ_D should match measured 645 K within 2%."""
    row = predict_theta_D(MATERIALS["Silicon"])
    assert abs(row["rel_err"]) < 0.02, (
        f"Silicon Θ_D = {row['theta_pred']:.1f} K, measured 645, "
        f"rel err {row['rel_err']:+.2%}"
    )


# --------------------------------------------------------------------------- #
# Claim (e): scaling behaviour                                                 #
# --------------------------------------------------------------------------- #


def test_theta_D_proportional_to_c_s():
    """Doubling c_s must double Θ_D (with n fixed)."""
    n = 8.49e28        # Cu number density
    base = debye_temperature(2620.0, n)
    doubled = debye_temperature(5240.0, n)
    assert math.isclose(doubled / base, 2.0, rel_tol=1.0e-12)


def test_theta_D_scales_with_n_cube_root():
    """Θ_D ∝ n^{1/3}: 8x density => 2x Θ_D."""
    cs = 3000.0
    base = debye_temperature(cs, 1.0e28)
    octupled = debye_temperature(cs, 8.0e28)
    assert math.isclose(octupled / base, 2.0, rel_tol=1.0e-12)


# --------------------------------------------------------------------------- #
# Claim (f): aggregate fit quality                                             #
# --------------------------------------------------------------------------- #


def test_at_least_10_materials_within_10pct():
    """Substrate prediction must hit ≥10/14 materials within 10%."""
    s = run_test()["summary"]
    assert s["within_10pct"] >= 10, (
        f"Only {s['within_10pct']}/14 within 10%"
    )


def test_at_least_12_materials_within_20pct():
    s = run_test()["summary"]
    assert s["within_20pct"] >= 12


def test_log_log_pearson_above_0p99():
    """Geometric correlation must be ≥0.99 across 14× Θ_D dynamic range."""
    s = run_test()["summary"]
    assert s["pearson_log_r"] > 0.99, (
        f"log-log pearson r = {s['pearson_log_r']:.4f}"
    )


def test_loglog_slope_near_unity():
    """log(pred) = m·log(meas) + b: slope m must be within 0.10 of 1.0."""
    s = run_test()["summary"]
    assert abs(s["loglog_slope"] - 1.0) < 0.10, (
        f"log-log slope = {s['loglog_slope']:.3f}"
    )


def test_mean_abs_residual_under_10pct():
    s = run_test()["summary"]
    assert s["mean_abs_rel_err"] < 0.10, (
        f"mean |rel err| = {s['mean_abs_rel_err']:.2%}"
    )


def test_mean_rel_err_unbiased_within_8pct():
    """Mean (signed) residual reflects systematic bias; should be small."""
    s = run_test()["summary"]
    assert abs(s["mean_rel_err"]) < 0.08, (
        f"mean signed err = {s['mean_rel_err']:+.2%}"
    )


# --------------------------------------------------------------------------- #
# Claim (g): visual renders                                                    #
# --------------------------------------------------------------------------- #


def test_render_debye_test_writes_png(tmp_path):
    out = tmp_path / "124_debye_test.png"
    path = render_debye_test(out_path=str(out))
    assert os.path.exists(path)
    assert os.path.getsize(path) > 5_000, "PNG suspiciously small"


# --------------------------------------------------------------------------- #
# Sanity: physical constants                                                   #
# --------------------------------------------------------------------------- #


def test_physical_constants_correct():
    """Spot-check: HBAR, KB, N_A are CODATA 2019 SI values."""
    assert math.isclose(HBAR, 1.054571817e-34, rel_tol=1.0e-12)
    assert math.isclose(KB,   1.380649e-23,    rel_tol=1.0e-12)
    assert math.isclose(N_A,  6.02214076e23,   rel_tol=1.0e-12)


# --------------------------------------------------------------------------- #
# Quasi-harmonic / anharmonic correction tests                                #
# --------------------------------------------------------------------------- #
#
# These tests verify the quasi-harmonic substrate correction:
#
#   γ_G_substrate = (3/2) · (1 + ν) / (2 - 3ν),  ν ≤ σ_max = 1/2
#   Θ_D^qh        = Θ_D^harm · (1 + γ_G · α_V · T_eff),  T_eff = T_melt / 2
#
# The substrate γ_G is bounded above by 9/2 = 4.5 (when ν saturates at the
# substrate elastic limit σ_max = 1/2). This is THE substrate-derivable
# upper cap on Grüneisen anharmonicity.


def test_sigma_max_poisson_is_half():
    """Substrate elastic-limit cap on Poisson's ratio matches saturation.SIGMA_MAX."""
    assert math.isclose(SIGMA_MAX_POISSON, 0.5, rel_tol=1.0e-12)


def test_material_has_alpha_V_and_T_melt():
    """Every material must carry α_V and T_melt for the QH correction."""
    for name, mat in MATERIALS.items():
        assert mat.alpha_V > 0.0, f"{name}: alpha_V missing"
        assert mat.T_melt > 100.0, f"{name}: T_melt missing or unphysical"
        assert mat.gamma_G_lit > 0.0, f"{name}: gamma_G_lit missing"


def test_poisson_ratio_below_substrate_cap():
    """Every material's Poisson ratio must satisfy ν ≤ σ_max = 1/2 (substrate cap)."""
    for name, mat in MATERIALS.items():
        v = mat.poisson_ratio
        assert -1.0 < v < SIGMA_MAX_POISSON, (
            f"{name}: Poisson ratio v={v:.3f} violates substrate cap"
        )


def test_gamma_G_substrate_bounded_by_substrate_cap():
    """γ_G_substrate ≤ 9/2 ≈ 4.5 (from σ_max=1/2 cap on ν)."""
    for name, mat in MATERIALS.items():
        gG = gamma_G_substrate(mat)
        assert 0.5 < gG < 4.6, f"{name}: γ_G_sub={gG:.3f} outside substrate cap"


def test_gamma_G_substrate_close_to_literature():
    """γ_G_substrate(ν) should match literature γ_G to within ~50% for most metals.

    Belomestnykh-Tesleva is an isotropic-medium estimate; departures up to 50%
    are expected for Au (which has anomalously low G/B due to relativistic
    band effects) and Si/Ge (where covalent character softens γ_G below the
    elastic estimate). We require ≥10/14 within 50% as a sanity floor.
    """
    n_close = 0
    for name, mat in MATERIALS.items():
        gG_sub = gamma_G_substrate(mat)
        gG_lit = mat.gamma_G_lit
        rel = abs(gG_sub - gG_lit) / gG_lit
        if rel < 0.5:
            n_close += 1
    assert n_close >= 10, (
        f"Only {n_close}/14 materials have γ_G_sub within 50% of literature"
    )


def test_quasi_harmonic_shift_always_increases_theta():
    """The QH multiplier (1 + γ_G·α_V·T_melt/2) must be ≥ 1 for every material."""
    for name, mat in MATERIALS.items():
        s = quasi_harmonic_shift(mat)
        assert s >= 1.0, f"{name}: QH shift {s} < 1"
        # Realistic upper bound: γ_G ≤ 4.5, α_V ≤ 100ppm, T_melt ≤ 5000
        # ⇒ shift ≤ 1 + 4.5·1e-4·2500 = 2.13. Cap at 1.5 for sanity.
        assert s < 1.5, f"{name}: QH shift {s} unphysically large"


def test_quasi_harmonic_improves_lead_and_tin():
    """QH correction must reduce |residual| for Pb and Sn (the soft-mode failures)."""
    pb = predict_with_quasiharmonic(MATERIALS["Lead"])
    sn = predict_with_quasiharmonic(MATERIALS["Tin"])
    assert abs(pb["rel_err_qh"]) < abs(pb["rel_err_harm"]), (
        f"Pb QH correction did not improve: harm {pb['rel_err_harm']:+.2%}, "
        f"qh {pb['rel_err_qh']:+.2%}"
    )
    assert abs(sn["rel_err_qh"]) < abs(sn["rel_err_harm"]), (
        f"Sn QH correction did not improve: harm {sn['rel_err_harm']:+.2%}, "
        f"qh {sn['rel_err_qh']:+.2%}"
    )


def test_tin_qh_within_10pct():
    """Sn quasi-harmonic prediction must reach within 10% of measured 200 K.

    Soft β-Sn responds well to QH correction (Δ from -12% → ~-9%). This is
    the headline 'fix' of the QH machinery on the original failure set.
    """
    row = predict_with_quasiharmonic(MATERIALS["Tin"])
    assert abs(row["rel_err_qh"]) < 0.10, (
        f"Sn QH residual {row['rel_err_qh']:+.2%} not within 10% of 200 K"
    )


def test_lead_qh_substantially_improves_but_does_not_reach_10pct():
    """Pb cannot reach 10% with substrate QH alone — Kohn anomaly is required.

    HONEST SCOPING: substrate elastic continuum + Belomestnykh-Tesleva γ_G
    + Lindemann thermal scale T_melt/2 reduces Pb residual from -27% to
    roughly -21% (~6pp improvement). The remaining gap is NOT addressable
    by any continuum quasi-harmonic correction; it requires explicit
    treatment of the X-point Kohn anomaly (electron-phonon coupling at the
    Fermi surface). We therefore assert improvement, not within-10%.
    """
    row = predict_with_quasiharmonic(MATERIALS["Lead"])
    # Improves by at least 4 percentage points
    improvement_pp = abs(row["rel_err_harm"]) - abs(row["rel_err_qh"])
    assert improvement_pp > 0.04, (
        f"Pb improvement {improvement_pp:+.2%} smaller than expected 4pp"
    )
    # Final QH residual lies in (-25%, -15%) — Kohn-anomaly regime
    assert -0.25 < row["rel_err_qh"] < -0.15, (
        f"Pb QH residual {row['rel_err_qh']:+.2%} outside expected Kohn-anomaly band"
    )


def test_within_10pct_count_stable_or_improves():
    """Total #materials within 10% of measured must NOT decrease under QH."""
    s = run_test_quasiharmonic()["summary"]
    assert s["within_10pct_qh"] >= s["within_10pct_harm"], (
        f"QH degraded count: harm={s['within_10pct_harm']}, "
        f"qh={s['within_10pct_qh']}"
    )


def test_qh_does_not_destroy_well_fit_materials():
    """Materials at <2% under harmonic must stay <10% under QH (no over-correction)."""
    for name, mat in MATERIALS.items():
        row = predict_with_quasiharmonic(mat)
        if abs(row["rel_err_harm"]) < 0.02:
            assert abs(row["rel_err_qh"]) < 0.10, (
                f"{name}: harm OK ({row['rel_err_harm']:+.2%}) but QH "
                f"over-corrected to {row['rel_err_qh']:+.2%}"
            )


def test_qh_summary_schema():
    """run_test_quasiharmonic must expose harm/qh comparison and Pb/Sn focus."""
    s = run_test_quasiharmonic()["summary"]
    needed = {
        "n_materials", "T_eff_default",
        "mean_abs_rel_err_harm", "mean_abs_rel_err_qh",
        "max_abs_rel_err_harm", "max_abs_rel_err_qh",
        "within_10pct_harm", "within_10pct_qh", "improvement_pp",
        "Pb_harm", "Pb_qh", "Pb_meas",
        "Sn_harm", "Sn_qh", "Sn_meas",
    }
    assert needed.issubset(set(s.keys())), (
        f"missing summary keys: {needed - set(s.keys())}"
    )


def test_qh_row_schema():
    """Every quasi-harmonic row must carry the expected fields."""
    rows = run_test_quasiharmonic()["rows"]
    expected = {
        "name", "lattice", "v_poisson", "gamma_G_sub", "gamma_G_lit",
        "alpha_V", "T_melt", "T_eff", "shift",
        "theta_harm", "theta_qh", "theta_meas",
        "rel_err_harm", "rel_err_qh",
    }
    for name, row in rows.items():
        assert set(row.keys()) == expected, f"{name}: schema mismatch"
        assert row["theta_qh"] >= row["theta_harm"], (
            f"{name}: QH lowered Θ_D, expected stiffening"
        )
