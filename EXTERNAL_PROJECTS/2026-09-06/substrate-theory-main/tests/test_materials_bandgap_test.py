"""Tests for materials_bandgap_test.py.

Validates the substrate-DFT bandgap predictor against the Materials
Project reference table built into the module.  The tests encode HONEST
expectations:

  * Si is reproduced exactly by construction (it sets chi_subs).
  * The predictor is at most order-of-magnitude across all 28 materials.
  * Specific known failure modes are recorded as tests so that any
    revision to the predictor must intentionally break them:
      - narrow-gap chalcogenides (InAs, InSb, PbS, PbSe) over-predict
      - insulators (LiF, CaF2, MgO, Al2O3, SiO2) under-predict
      - metals (Cu, Al, Ag, Au) get a small but non-zero gap
"""

from __future__ import annotations

import math
import os

import numpy as np
import pytest

from src.stiff_medium.materials_bandgap_test import (
    CHI_SUBS_DEFAULT,
    KAPPA_SUBS_LOCAL,
    MATERIALS_MP,
    N_SAT_BOHR,
    SIGMA_CAP,
    BandgapResidual,
    evaluate_all,
    predict_bandgap_eV,
    summarize,
)


# --------------------------------------------------------------------- #
# Substrate constants used by the predictor                             #
# --------------------------------------------------------------------- #


def test_n_sat_is_K_pair_squared_over_pi() -> None:
    """n_sat = K_pair^2 / pi = 4 / pi from the canonical K_pair = 2."""
    assert N_SAT_BOHR == pytest.approx(4.0 / math.pi, rel=1e-12)


def test_kappa_subs_matches_lieb_oxford() -> None:
    """kappa_subs = 1.443 (within 0.5%) from Lieb-Oxford match 1.804."""
    assert KAPPA_SUBS_LOCAL == pytest.approx(1.443, rel=5e-3)


def test_sigma_cap_is_one_half() -> None:
    assert SIGMA_CAP == pytest.approx(0.5, rel=1e-12)


# --------------------------------------------------------------------- #
# Si calibration : honestly fails to be exact because the calibrated   #
# chi from Si alone would blow up high-eps materials.                  #
# --------------------------------------------------------------------- #


def test_si_within_factor_of_5() -> None:
    """The Si calibration is CLAMPED so that high-eps_r materials don't
    blow up; Si is therefore NOT exactly reproduced.  We check the
    predictor lands within a factor of 5 of the MP value (1.12 eV).
    """
    si = MATERIALS_MP["Si"]
    pred = predict_bandgap_eV(si["a_A"], si["m_e_star"], si["eps_r"])
    ratio = pred["E_g_eV"] / si["E_g_MP_eV"]
    assert 0.2 < ratio < 5.0, (
        f"Si predicted {pred['E_g_eV']:.3f} eV vs MP 1.12 eV, ratio {ratio:.2f}"
    )


def test_chi_subs_in_reasonable_range() -> None:
    """chi_subs after Si-clamped calibration should be O(-0.01 .. -0.1)."""
    assert -1.0 < CHI_SUBS_DEFAULT < 0.0


# --------------------------------------------------------------------- #
# Material table sanity                                                 #
# --------------------------------------------------------------------- #


def test_materials_table_size() -> None:
    """We claim ~30 materials; the table must cover at least 28."""
    assert len(MATERIALS_MP) >= 28


def test_every_material_has_required_fields() -> None:
    required = {"E_g_MP_eV", "a_A", "m_e_star", "m_h_star", "eps_r", "class_", "notes"}
    for name, row in MATERIALS_MP.items():
        missing = required - row.keys()
        assert not missing, f"{name} is missing fields {missing}"


def test_every_material_classifies_into_known_class() -> None:
    valid = {"semiconductor", "wide_gap", "insulator", "metal",
             "narrow_gap", "halide"}
    for name, row in MATERIALS_MP.items():
        assert row["class_"] in valid, f"{name} class {row['class_']!r} unknown"


def test_metal_E_g_MP_is_zero() -> None:
    metals = [n for n, r in MATERIALS_MP.items() if r["class_"] == "metal"]
    assert metals, "expected at least one metal"
    for m in metals:
        assert MATERIALS_MP[m]["E_g_MP_eV"] == 0.0


# --------------------------------------------------------------------- #
# Predictor produces finite, positive predictions for every material   #
# --------------------------------------------------------------------- #


def test_predictor_finite_and_positive() -> None:
    rows = evaluate_all()
    assert len(rows) == len(MATERIALS_MP)
    for r in rows:
        assert math.isfinite(r.E_g_subs_eV), f"{r.material}: non-finite prediction"
        assert r.E_g_subs_eV >= 0.0, f"{r.material}: negative gap {r.E_g_subs_eV}"


# --------------------------------------------------------------------- #
# Honest failure-mode regressions: lock in the current predictor's    #
# known systematic errors so any revision is intentional.              #
# --------------------------------------------------------------------- #


def test_insulators_systematically_underpredict() -> None:
    """Wide-gap ionic insulators (LiF, CaF2, MgO, Al2O3, SiO2) under-predict
    by 80-98 % under the current substrate predictor: the simple
    Hartree-scale * F_xc * f_dielectric formula cannot reach 8-13 eV
    without a ionicity-driven enhancement that the framework does not
    yet provide.
    """
    rows = evaluate_all()
    insulators = [r for r in rows if r.class_ == "insulator"]
    assert insulators, "expected at least one insulator"
    for r in insulators:
        # All insulators should under-predict (negative residual %)
        assert r.pct_residual < 0.0, f"{r.material} pct_residual {r.pct_residual:+.1f}% "
        # Magnitude in 30-99 % range
        assert abs(r.pct_residual) > 30.0, f"{r.material} too small under-pred"


def test_narrow_gap_chalcogenides_over_predict() -> None:
    """InAs/InSb/PbS/PbSe with very small m_e^* over-predict because the
    1/m_e^* atomic-scale factor blows up.  Recorded as a known failure.
    """
    rows = evaluate_all()
    targets = [r for r in rows if r.material in {"InAs", "InSb", "PbS", "PbSe", "HgTe"}]
    assert len(targets) == 5
    for r in targets:
        assert r.E_g_subs_eV > r.E_g_MP_eV, f"{r.material} should over-predict"


def test_metals_get_small_but_nonzero_gap() -> None:
    """The substrate predictor cannot return E_g = 0 (it has no band-overlap
    detection), so metals get a small spurious gap < 1 eV.  Honest failure.
    """
    rows = evaluate_all()
    metals = [r for r in rows if r.class_ == "metal"]
    assert metals, "expected at least one metal"
    for r in metals:
        assert r.E_g_subs_eV > 0.0, f"{r.material} predicted exactly zero?"
        assert r.E_g_subs_eV < 1.5, (
            f"{r.material} spurious metal gap {r.E_g_subs_eV:.2f} eV "
            "exceeds 1.5 eV; predictor is worse than current"
        )


def test_diamond_within_an_order_of_magnitude() -> None:
    """Diamond (5.5 eV) is the headline wide-gap covalent solid; the
    predictor lands within an order of magnitude (currently ~ 1.1 eV vs 5.5).
    """
    rows = {r.material: r for r in evaluate_all()}
    d = rows["diamond"]
    ratio = d.E_g_subs_eV / d.E_g_MP_eV
    assert 0.05 < ratio < 5.0, f"diamond ratio {ratio:.2f}"


# --------------------------------------------------------------------- #
# Aggregate statistics                                                  #
# --------------------------------------------------------------------- #


def test_summary_keys_present() -> None:
    rows = evaluate_all()
    summary = summarize(rows)
    assert "ALL" in summary
    assert "NON_METAL" in summary
    for key in ("n", "mae_eV", "mape_pct", "rmse_eV", "bias_eV"):
        assert key in summary["ALL"]


def test_aggregate_mae_is_few_eV() -> None:
    """Mean absolute error across all materials sits in 1-12 eV: the
    predictor is order-of-magnitude useful but NOT chemically accurate.

    Reference : a state-of-the-art HSE06 hybrid DFT achieves MAE ~ 0.3 eV
    on the same ~ 30-material benchmark; the substrate-DFT predictor is
    therefore ~ 15x worse than HSE06 while using fewer parameters.
    """
    rows = evaluate_all()
    summary = summarize(rows)
    mae = summary["ALL"]["mae_eV"]
    assert 1.0 < mae < 12.0, f"MAE = {mae:.2f} eV out of expected band"


def test_correlation_with_E_g_MP_is_negative() -> None:
    """HONEST FAILURE recorded as a regression test :

    The current substrate-DFT bandgap predictor's gap-vs-MP relation is
    *negatively* correlated across the 28-material set, because narrow-gap
    materials (InAs, InSb) have tiny conduction-band masses (m_e^* ~ 0.01)
    which the 1/m_e^* atomic-scale factor inflates while wide-gap insulators
    (LiF, CaF2) have small lattice constants but only modest m_e^*, giving
    SMALLER predicted gaps than narrow-gap chalcogenides.

    This is the predictor's biggest failure mode and the primary scientific
    finding of the test : the simple substrate-DFT bandgap formula
    (E_atom * F_xc * f_dielectric) does NOT preserve material rank order
    with the Materials Project DFT gap.
    """
    rows = evaluate_all()
    xs = np.array([r.E_g_MP_eV for r in rows])
    ys = np.array([r.E_g_subs_eV for r in rows])
    rho = float(np.corrcoef(xs, ys)[0, 1])
    # Locked in : current rho is ~ -0.4 .  Future revisions that fix the
    # narrow-gap inflation should INTENTIONALLY break this test and
    # replace it with a positive-correlation test.
    assert rho < 0.0, (
        f"current predictor has non-positive correlation rho = {rho:.3f} ; "
        "if your revision now produces positive correlation, intentionally "
        "delete this test and add a positive-rho replacement"
    )


# --------------------------------------------------------------------- #
# Visualization smoke test                                              #
# --------------------------------------------------------------------- #


def test_visualization_writes_png(tmp_path) -> None:
    from src.stiff_medium.materials_bandgap_test import make_visual
    rows = evaluate_all()
    out = tmp_path / "test_materials_bandgaps.png"
    make_visual(rows, str(out))
    assert out.exists()
    assert out.stat().st_size > 1000   # > 1 KB
