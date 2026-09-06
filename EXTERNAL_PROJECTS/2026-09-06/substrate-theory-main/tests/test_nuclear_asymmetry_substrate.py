"""Tests for nuclear_asymmetry_substrate -- substrate-derived a_sym."""

from __future__ import annotations

import pytest

from src.stiff_medium.b3_constants import (
    K_pair,
    K_rank,
    LAMBDA_QCD_MEV,
    N_BAM,
    R,
    n_A,
)
from src.stiff_medium.k4_face_pair_geometry import EPS_FACE_MEV
from src.stiff_medium.nuclear_asymmetry_substrate import (
    A_SYM_EMPIRICAL_MEV,
    AsymmetryDerivation,
    a_sym_koide_form,
    a_sym_lambda_form,
    a_sym_relative_error,
    a_sym_substrate,
    report,
)


# ---------------------------------------------------------------------------
# Numerical value of substrate a_sym
# ---------------------------------------------------------------------------

def test_a_sym_substrate_equals_eps_times_kpair_times_krank():
    """a_sym = ε_face × K_pair × K_rank from primitives."""
    expected = EPS_FACE_MEV * K_pair * K_rank
    assert abs(a_sym_substrate() - expected) < 1e-9


def test_a_sym_substrate_value_near_22p2_MeV():
    """Substrate a_sym ≈ 22.22 MeV (= 200/9)."""
    assert abs(a_sym_substrate() - 22.222) < 0.01


def test_a_sym_within_30_pct_of_empirical():
    """a_sym from substrate within 30 % of empirical 23 MeV (brief target)."""
    rel = abs(a_sym_relative_error())
    assert rel < 0.30, (
        f"substrate a_sym = {a_sym_substrate():.3f} MeV deviates from "
        f"empirical {A_SYM_EMPIRICAL_MEV} MeV by {100*rel:.2f}% > 30%"
    )


def test_a_sym_within_5_pct_of_empirical():
    """Stricter: substrate a_sym is within 5 % of empirical (it lands at 3.4 %)."""
    rel = abs(a_sym_relative_error())
    assert rel < 0.05


# ---------------------------------------------------------------------------
# Equivalent forms (numerically identical)
# ---------------------------------------------------------------------------

def test_lambda_form_equals_primitive_form():
    """Λ · K_pair·K_rank / (n_A·N_BAM) == ε_face · K_pair · K_rank.

    The two forms are mathematically identical (both = Λ_QCD/9), but
    EPS_FACE_MEV is rounded to 2.222 in ``k4_face_pair_geometry`` while
    the Λ-form keeps full float precision -- hence the 0.01 MeV tolerance.
    """
    assert abs(a_sym_lambda_form() - a_sym_substrate()) < 0.01


def test_koide_form_equals_primitive_form():
    """Λ / R² == 200 / 9 == ε_face · K_pair · K_rank (within ε_face roundoff)."""
    assert abs(a_sym_koide_form() - a_sym_substrate()) < 0.01


def test_lambda_form_equals_lambda_over_9():
    """The Λ-anchored form is exactly Λ_QCD / 9."""
    assert abs(a_sym_lambda_form() - LAMBDA_QCD_MEV / 9.0) < 1e-9


def test_koide_denominator_identity_holds():
    """K_pair · K_rank − 1 = R² = 9 holds exactly."""
    assert K_pair * K_rank - 1 == R * R == 9


# ---------------------------------------------------------------------------
# Report record
# ---------------------------------------------------------------------------

def test_report_returns_asymmetry_derivation_record():
    """report() returns an AsymmetryDerivation with substrate / empirical fields."""
    rep = report()
    assert isinstance(rep, AsymmetryDerivation)
    assert rep.a_sym_substrate_MeV == a_sym_substrate()
    assert rep.a_sym_empirical_MeV == A_SYM_EMPIRICAL_MEV
    assert "ε_face" in rep.derivation


def test_report_relative_error_signed():
    """err_pct sign reflects substrate < empirical -> negative."""
    rep = report()
    assert rep.err_pct < 0           # 22.22 < 23 -> negative
    assert abs(rep.err_pct) < 5.0    # 3.4 % match


def test_report_abs_error_about_0p8_MeV():
    """|22.22 - 23| ≈ 0.78 MeV."""
    rep = report()
    assert 0.5 < rep.err_abs_MeV < 1.0


# ---------------------------------------------------------------------------
# Substrate primitives consistency
# ---------------------------------------------------------------------------

def test_substrate_inputs_are_canonical():
    """Inputs (ε_face, K_pair, K_rank) match the canonical primitives."""
    assert abs(EPS_FACE_MEV - 2.222) < 0.01
    assert K_pair == 2
    assert K_rank == 5
    assert N_BAM == 6
    assert n_A == 15


def test_a_sym_is_positive():
    """Asymmetry coefficient is a positive penalty scale."""
    assert a_sym_substrate() > 0
