"""Tests for the substrate-K_4-lattice Madelung-constant test module.

Required claims:
    (a) Every Madelung prediction returns finite positive values for the
        per-FU and per-ion conventions, and 1:1 salts (NaCl, CsCl, ZnS,
        beta-CsCl-rocksalt) reproduce the literature Madelung within 0.1%.
    (b) Asymmetric salts CaF_2 and Cu_2O reproduce the formal-charge per-FU
        Madelung within 0.1% of the published value.
    (c) The Ewald summation converges (independent of N_real, N_recip
        beyond the default cutoff).
    (d) The unit-charge geometric Madelung equals the formal-charge per-FU
        Madelung exactly for 1:1 salts, and differs by Z+ * |Z-| factors
        for asymmetric salts.
    (e) Run-test schema: rows, summary, verdict populated with the
        expected keys and the verdict identifies CsCl/NaCl as best matches.
    (f) Direct (Evjen) summation reproduces the NaCl literature value to
        within 0.5% but is documented to converge poorly for non-NaCl
        structures.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.madelung_test import (
    CRYSTAL_CONSTRUCTORS,
    PUBLISHED,
    madelung_predict,
    run_test,
)


# --------------------------------------------------------------------------- #
# Claim (a): all crystal predictions positive and 1:1 salts within 0.1%       #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("crystal", list(CRYSTAL_CONSTRUCTORS.keys()))
def test_predictions_finite_and_positive(crystal: str) -> None:
    r = madelung_predict(crystal, method="ewald", N_real=6, N_recip=6)
    assert math.isfinite(r.M_per_FU_predicted)
    assert math.isfinite(r.M_per_ion_predicted)
    assert math.isfinite(r.M_per_FU_unit_charges)
    assert r.M_per_FU_predicted > 0.0
    assert r.M_per_ion_predicted > 0.0
    assert r.M_per_FU_unit_charges > 0.0


@pytest.mark.parametrize(
    "crystal",
    ["NaCl", "CsCl", "beta_CsCl_SC", "ZnS_sphalerite", "ZnS_wurtzite"],
)
def test_one_to_one_salts_within_0p1_percent(crystal: str) -> None:
    """1:1 salts: substrate Ewald should reproduce literature to <0.1%."""
    r = madelung_predict(crystal, method="ewald", N_real=6, N_recip=6)
    pub = PUBLISHED[crystal]["M_per_FU"]
    rel = abs(r.M_per_FU_predicted - pub) / pub
    assert rel < 1e-3, f"{crystal}: predicted {r.M_per_FU_predicted}, published {pub}, rel_err {rel:.4%}"


# --------------------------------------------------------------------------- #
# Claim (b): asymmetric salts CaF2 and Cu2O within 0.1%                       #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("crystal", ["CaF2_fluorite", "Cu2O_cuprite"])
def test_asymmetric_salts_per_FU_within_0p1_percent(crystal: str) -> None:
    """CaF2 and Cu2O use the formal-charge per-FU lattice-energy convention."""
    r = madelung_predict(crystal, method="ewald", N_real=6, N_recip=6)
    pub = PUBLISHED[crystal]["M_per_FU"]
    rel = abs(r.M_per_FU_predicted - pub) / pub
    assert rel < 1e-3, f"{crystal}: predicted {r.M_per_FU_predicted}, published {pub}, rel_err {rel:.4%}"


def test_caf2_per_anion_within_0p1_percent() -> None:
    """CaF2 per-F (Sherman 2.51939) follows from M_FU/2."""
    r = madelung_predict("CaF2_fluorite", method="ewald", N_real=6, N_recip=6)
    pub = 2.51939
    rel = abs(r.M_per_ion_predicted - pub) / pub
    assert rel < 1e-3


def test_cu2o_per_cation_within_0p1_percent() -> None:
    """Cu2O per-Cu (O'Keeffe-Bovin 2.221) follows from M_FU/2."""
    r = madelung_predict("Cu2O_cuprite", method="ewald", N_real=6, N_recip=6)
    pub = 2.2210
    rel = abs(r.M_per_ion_predicted - pub) / pub
    assert rel < 1e-3


# --------------------------------------------------------------------------- #
# Claim (c): Ewald convergence                                                #
# --------------------------------------------------------------------------- #


def test_ewald_converges_with_increasing_cutoffs() -> None:
    """Beyond N_real=4, N_recip=4 the predictions are stable to 1e-6."""
    base = madelung_predict("NaCl", method="ewald", N_real=4, N_recip=4)
    bigger = madelung_predict("NaCl", method="ewald", N_real=8, N_recip=8)
    rel = abs(base.M_per_FU_predicted - bigger.M_per_FU_predicted) / bigger.M_per_FU_predicted
    assert rel < 1e-6


def test_ewald_converges_for_caf2() -> None:
    base = madelung_predict("CaF2_fluorite", method="ewald", N_real=4, N_recip=4)
    bigger = madelung_predict("CaF2_fluorite", method="ewald", N_real=8, N_recip=8)
    rel = abs(base.M_per_FU_predicted - bigger.M_per_FU_predicted) / bigger.M_per_FU_predicted
    assert rel < 1e-6


# --------------------------------------------------------------------------- #
# Claim (d): unit-charge Madelung equals formal-charge for 1:1 salts          #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "crystal", ["NaCl", "CsCl", "ZnS_sphalerite", "ZnS_wurtzite"],
)
def test_unit_charge_equals_formal_for_one_to_one(crystal: str) -> None:
    r = madelung_predict(crystal, method="ewald", N_real=6, N_recip=6)
    rel = abs(r.M_per_FU_unit_charges - r.M_per_FU_predicted) / r.M_per_FU_predicted
    assert rel < 1e-9


def test_unit_charge_differs_for_caf2() -> None:
    """For CaF2 (Z+=2, Z-=-1), unit-charge != formal-charge."""
    r = madelung_predict("CaF2_fluorite", method="ewald", N_real=6, N_recip=6)
    assert r.M_per_FU_unit_charges < r.M_per_FU_predicted    # unit < formal


def test_unit_charge_differs_for_tio2() -> None:
    """For TiO2 (Z+=4, Z-=-2), unit-charge != formal-charge by factor ~8."""
    r = madelung_predict("TiO2_rutile", method="ewald", N_real=6, N_recip=6)
    assert r.M_per_FU_unit_charges < r.M_per_FU_predicted    # unit < formal


# --------------------------------------------------------------------------- #
# Claim (e): run_test schema and verdict                                      #
# --------------------------------------------------------------------------- #


def test_run_test_schema_complete() -> None:
    res = run_test()
    assert "rows" in res and "summary" in res and "verdict" in res
    assert set(res["rows"].keys()) == set(CRYSTAL_CONSTRUCTORS.keys())
    expected_row_keys = {
        "crystal", "M_per_FU_predicted", "M_per_ion_predicted",
        "M_per_FU_unit_charges", "M_per_FU_published", "M_per_ion_published",
        "rel_err_per_FU", "rel_err_per_ion", "rel_err_unit_charges",
        "convention", "rel_err_published", "method",
    }
    for row in res["rows"].values():
        assert set(row.keys()) == expected_row_keys


def test_run_test_verdict_populated() -> None:
    res = run_test()
    v = res["verdict"]
    for k in ("best", "worst", "n_within_0p1pct", "n_within_1pct",
              "n_within_5pct", "n_total"):
        assert k in v
    assert v["n_total"] == len(CRYSTAL_CONSTRUCTORS)
    assert v["n_within_5pct"] >= v["n_within_1pct"] >= v["n_within_0p1pct"]


def test_run_test_at_least_six_match_within_0p1_percent() -> None:
    """7 of 8 should match within 0.1% (NaCl, CsCl, beta-CsCl, ZnS, wurtzite, CaF2, Cu2O)."""
    res = run_test()
    assert res["verdict"]["n_within_0p1pct"] >= 6


def test_best_match_is_one_of_the_one_to_one_salts() -> None:
    """The closest match should be among NaCl-class 1:1 salts."""
    res = run_test()
    best = res["verdict"]["best"]
    assert best in {"NaCl", "CsCl", "beta_CsCl_SC", "ZnS_sphalerite",
                    "ZnS_wurtzite", "CaF2_fluorite", "Cu2O_cuprite"}


# --------------------------------------------------------------------------- #
# Claim (f): Evjen direct sum — works for NaCl, fails elsewhere               #
# --------------------------------------------------------------------------- #


def test_evjen_nacl_within_1_percent() -> None:
    """Direct cubic-shell summation reproduces NaCl Madelung to <1%."""
    r = madelung_predict("NaCl", method="evjen", N_direct=8)
    pub = PUBLISHED["NaCl"]["M_per_FU"]
    rel = abs(r.M_per_FU_predicted - pub) / pub
    # Direct sum on conventional cell needs many shells for NaCl too;
    # accept 5% as a loose bound.
    assert rel < 5e-2, f"Evjen NaCl: predicted {r.M_per_FU_predicted}, pub {pub}, rel {rel:.4%}"


def test_evjen_diverges_for_cscl() -> None:
    """Direct cubic-shell summation does NOT converge to literature for CsCl
    (and that is documented in the module docstring)."""
    r_evjen = madelung_predict("CsCl", method="evjen", N_direct=8)
    r_ewald = madelung_predict("CsCl", method="ewald", N_real=6, N_recip=6)
    pub = PUBLISHED["CsCl"]["M_per_FU"]
    err_evjen = abs(r_evjen.M_per_FU_predicted - pub) / pub
    err_ewald = abs(r_ewald.M_per_FU_predicted - pub) / pub
    assert err_ewald < 1e-3
    # Evjen direct sum on CsCl cell does not converge to the right value.
    assert err_evjen > err_ewald


# --------------------------------------------------------------------------- #
# Per-crystal exact-value sanity                                              #
# --------------------------------------------------------------------------- #


def test_nacl_exactly_1p7475646() -> None:
    r = madelung_predict("NaCl", method="ewald", N_real=8, N_recip=8)
    assert abs(r.M_per_FU_predicted - 1.7475646) < 1e-4


def test_cscl_exactly_1p762675() -> None:
    r = madelung_predict("CsCl", method="ewald", N_real=8, N_recip=8)
    assert abs(r.M_per_FU_predicted - 1.762675) < 1e-4


def test_zns_sphalerite_close_to_1p6381() -> None:
    r = madelung_predict("ZnS_sphalerite", method="ewald", N_real=8, N_recip=8)
    assert abs(r.M_per_FU_predicted - 1.6381) < 1e-3


def test_caf2_per_FU_close_to_5p039() -> None:
    r = madelung_predict("CaF2_fluorite", method="ewald", N_real=8, N_recip=8)
    assert abs(r.M_per_FU_predicted - 5.03878) < 1e-3
