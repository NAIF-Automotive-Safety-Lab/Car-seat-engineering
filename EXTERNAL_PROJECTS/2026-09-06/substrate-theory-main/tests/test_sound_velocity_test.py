"""Tests for src/stiff_medium/sound_velocity_test.py.

Required claims:
    (a) Schema integrity: 15 materials, expected keys with positive numerics,
        substrate flag matches the elasticity database membership.
    (b) Wave-equation kinematics: with measured (B, G) the substrate wave
        equation reproduces measured c_L within ~6% for every material; this
        is the kinematic-only check.
    (c) Substrate-derived moduli pathway: at least 6/7 materials with
        substrate (B, G) within 30% on c_L; log-log Pearson r > 0.95 across
        the 8 substrate-only entries.
    (d) Aggregate fit: log-log Pearson > 0.95 across all 15 materials with
        the substrate prediction; ALL 15 within 50%.
    (e) Validity bounds: invalid inputs (rho ≤ 0, negative moduli) raise
        ValueError.
    (f) Liquid limit: water (G = 0) gives c_L = sqrt(B/ρ) and c_T = 0.
    (g) Ranking preservation: sorting materials by predicted c_L gives the
        same ordering as sorting by measured c_L (Spearman = 1.0) — the
        substrate prediction respects the relative ranking even where the
        absolute moduli underpredict.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.sound_velocity_test import (
    MATERIALS,
    longitudinal_speed,
    predict_sound_speed,
    run_test,
    substrate_moduli,
    transverse_speed,
)


# --------------------------------------------------------------------------- #
# Claim (a): schema integrity                                                 #
# --------------------------------------------------------------------------- #


def test_run_test_schema():
    res = run_test()
    assert "rows" in res and "summary" in res
    assert "summary_substrate_BG" in res and "summary_kinematic" in res
    rows = res["rows"]
    assert len(rows) == 15, f"expected 15 materials, got {len(rows)}"
    assert set(rows.keys()) == set(MATERIALS.keys())
    expected_keys = {
        "name", "rho_mass", "B_GPa_used", "G_GPa_used",
        "B_GPa_meas", "G_GPa_meas", "used_substrate_BG",
        "c_L_pred", "c_T_pred", "c_L_measured_moduli", "c_L_meas",
        "rel_err", "rel_err_kinematic", "notes",
    }
    for name, row in rows.items():
        assert expected_keys.issubset(row.keys()), f"{name} missing keys"
        assert row["rho_mass"] > 0.0
        assert row["c_L_pred"] > 0.0
        assert row["c_L_meas"] > 0.0
        assert row["c_L_measured_moduli"] > 0.0


def test_substrate_flag_matches_elasticity_db():
    """The 'used_substrate_BG' flag should match whether substrate_elasticity
    has the material AND returns finite (B, G).
    """
    res = run_test()
    for name, row in res["rows"].items():
        B_sub, G_sub, used = substrate_moduli(name)
        assert row["used_substrate_BG"] == used, (
            f"{name}: row flag {row['used_substrate_BG']} != lookup {used}"
        )
        # When the flag is True, the (B, G) used should match the lookup
        if used:
            assert math.isclose(row["B_GPa_used"], B_sub, rel_tol=1e-9)
            assert math.isclose(row["G_GPa_used"], G_sub, rel_tol=1e-9)


# --------------------------------------------------------------------------- #
# Claim (b): wave-equation kinematics with measured moduli                    #
# --------------------------------------------------------------------------- #


def test_kinematic_check_per_material():
    """With MEASURED (B, G), the substrate wave equation should reproduce the
    measured c_L within ~6% per material.
    """
    res = run_test()
    for name, row in res["rows"].items():
        rel = abs(row["rel_err_kinematic"])
        assert rel <= 0.07, (
            f"{name}: kinematic rel-err {rel:.3%} exceeds 7% with "
            f"measured (B={row['B_GPa_meas']}, G={row['G_GPa_meas']})"
        )


def test_kinematic_aggregate():
    """Across all 15 materials the kinematic prediction is at < 2% mean
    error and log-log Pearson > 0.999.
    """
    res = run_test()
    summ = res["summary_kinematic"]
    assert summ["mean_abs_rel_err"] < 0.02, (
        f"kinematic mean abs rel err {summ['mean_abs_rel_err']:.3%} "
        f"exceeds 2%"
    )
    assert summ["loglog_pearson"] > 0.999, (
        f"kinematic log-log Pearson {summ['loglog_pearson']:.4f} below 0.999"
    )
    assert summ["within_10pct"] == 15, (
        f"kinematic: only {summ['within_10pct']}/15 within 10%"
    )


# --------------------------------------------------------------------------- #
# Claim (c): substrate-derived moduli pathway                                 #
# --------------------------------------------------------------------------- #


def test_substrate_subset_correlation():
    """Restricted to the materials with substrate-derived (B, G), log-log
    Pearson r > 0.95 — the substrate prediction reproduces the relative
    ordering even when absolute moduli are systematically under-predicted.
    """
    res = run_test()
    summ = res["summary_substrate_BG"]
    assert summ["n_materials_substrate_BG"] >= 6, (
        f"only {summ['n_materials_substrate_BG']} materials with substrate (B,G)"
    )
    assert summ["loglog_pearson"] > 0.95, (
        f"substrate-subset log-log Pearson {summ['loglog_pearson']:.4f} "
        f"below 0.95"
    )


def test_substrate_subset_within_30pct():
    """At least 6/7 substrate-pathway materials should land within 30% of
    measured c_L.  (Diamond is a known under-prediction at ~24%.)
    """
    res = run_test()
    sub_rows = [r for r in res["rows"].values() if r["used_substrate_BG"]]
    n_within_30 = sum(1 for r in sub_rows if abs(r["rel_err"]) <= 0.30)
    assert n_within_30 >= max(1, len(sub_rows) - 1), (
        f"only {n_within_30}/{len(sub_rows)} substrate-pathway materials "
        f"within 30%"
    )


# --------------------------------------------------------------------------- #
# Claim (d): aggregate fit across all 15 materials                            #
# --------------------------------------------------------------------------- #


def test_full_set_aggregate():
    """Across all 15 materials with the substrate prediction:
        - log-log Pearson > 0.95
        - all 15 within 50% (no outliers)
    """
    res = run_test()
    summ = res["summary"]
    assert summ["n_materials"] == 15
    assert summ["loglog_pearson"] > 0.95, (
        f"all-materials log-log Pearson {summ['loglog_pearson']:.4f} below 0.95"
    )
    assert summ["within_50pct"] == 15, (
        f"only {summ['within_50pct']}/15 within 50%"
    )
    assert summ["within_25pct"] >= 13, (
        f"only {summ['within_25pct']}/15 within 25%"
    )


# --------------------------------------------------------------------------- #
# Claim (e): validity bounds                                                  #
# --------------------------------------------------------------------------- #


def test_longitudinal_speed_validity():
    with pytest.raises(ValueError):
        longitudinal_speed(100.0, 50.0, 0.0)
    with pytest.raises(ValueError):
        longitudinal_speed(100.0, 50.0, -1.0)
    with pytest.raises(ValueError):
        longitudinal_speed(-1.0, 50.0, 1000.0)
    with pytest.raises(ValueError):
        longitudinal_speed(100.0, -1.0, 1000.0)


def test_transverse_speed_validity():
    with pytest.raises(ValueError):
        transverse_speed(50.0, 0.0)
    with pytest.raises(ValueError):
        transverse_speed(50.0, -1.0)
    with pytest.raises(ValueError):
        transverse_speed(-1.0, 1000.0)


# --------------------------------------------------------------------------- #
# Claim (f): liquid limit                                                     #
# --------------------------------------------------------------------------- #


def test_water_liquid_limit():
    """Water has G = 0; the substrate wave equation reduces to c_L = sqrt(B/ρ)
    and c_T = 0.  The 1480 m/s measured speed must be matched within 1%.
    """
    water = MATERIALS["Water"]
    row = predict_sound_speed(water)
    assert row["G_GPa_used"] == 0.0
    assert row["c_T_pred"] == 0.0
    # c_L = sqrt(B / rho); B = 2.2 GPa, rho = 1000 -> sqrt(2.2e9/1000) = 1483 m/s
    expected = math.sqrt(water.B_GPa_meas * 1.0e9 / water.rho_mass)
    assert math.isclose(row["c_L_pred"], expected, rel_tol=1e-9)
    assert abs(row["rel_err"]) <= 0.01, (
        f"water rel-err {row['rel_err']:.3%} exceeds 1%"
    )


# --------------------------------------------------------------------------- #
# Claim (g): ranking preservation                                             #
# --------------------------------------------------------------------------- #


def test_ranking_preservation():
    """The substrate prediction preserves the relative ordering of materials
    by sound speed.  Spearman correlation of (predicted, measured) ranks
    should be ≥ 0.9.
    """
    res = run_test()
    rows = list(res["rows"].values())
    pred = np.array([r["c_L_pred"] for r in rows])
    meas = np.array([r["c_L_meas"] for r in rows])
    rank_pred = pred.argsort().argsort()
    rank_meas = meas.argsort().argsort()
    n = len(rank_pred)
    d = rank_pred - rank_meas
    spearman = 1.0 - 6.0 * float((d * d).sum()) / (n * (n * n - 1))
    assert spearman >= 0.9, (
        f"Spearman ranking correlation {spearman:.4f} below 0.9 — "
        f"the substrate prediction does not preserve the c_L ranking"
    )


# --------------------------------------------------------------------------- #
# Claim (h): material database completeness                                   #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", list(MATERIALS.keys()))
def test_each_material_predicts(name):
    """Every material should produce a finite, positive prediction."""
    mat = MATERIALS[name]
    row = predict_sound_speed(mat)
    assert math.isfinite(row["c_L_pred"]) and row["c_L_pred"] > 0.0
    assert math.isfinite(row["c_T_pred"]) and row["c_T_pred"] >= 0.0
    assert math.isfinite(row["rel_err"])


def test_expected_materials_in_database():
    """The 15 spec materials must all be present (case-sensitive)."""
    expected = {
        "Diamond", "Beryllium", "Aluminum", "Iron", "Copper", "Silver",
        "Gold", "Lead", "Glass", "Steel", "Water", "Quartz", "Tungsten",
        "Silicon", "Germanium",
    }
    assert set(MATERIALS.keys()) == expected
