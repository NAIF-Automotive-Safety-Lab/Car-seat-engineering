"""Tests for atomic_spectroscopy_substrate."""

import math

import pytest

from stiff_medium.atomic_spectroscopy_substrate import (
    AtomicSpectroscopy,
    NIST,
)


@pytest.fixture(scope="module")
def spec():
    s = AtomicSpectroscopy()
    s.compare_to_NIST()
    return s


def _close(a, b, pct):
    return abs(a - b) / abs(b) * 100.0 < pct


# --- hydrogen series -------------------------------------------------------

def test_hydrogen_balmer_alpha(spec):
    val = spec.results["H_Balmer_alpha_nm"]
    assert _close(val, NIST["H_Balmer_alpha_nm"], 0.1), val


def test_hydrogen_lyman_alpha(spec):
    val = spec.results["H_Lyman_alpha_nm"]
    assert _close(val, NIST["H_Lyman_alpha_nm"], 0.1), val


def test_hydrogen_paschen_alpha(spec):
    val = spec.results["H_Paschen_alpha_nm"]
    assert _close(val, NIST["H_Paschen_alpha_nm"], 0.1), val


def test_hydrogen_all_within_1pct(spec):
    for k in (
        "H_Lyman_alpha_nm", "H_Lyman_beta_nm",
        "H_Balmer_alpha_nm", "H_Balmer_beta_nm",
        "H_Paschen_alpha_nm",
    ):
        ref_key = k  # NIST keys match
        assert _close(spec.results[k], NIST[ref_key], 1.0), k


# --- QED splittings --------------------------------------------------------

def test_lamb_shift(spec):
    val = spec.results["Lamb_2S2P_MHz"]
    assert _close(val, NIST["Lamb_shift_2S_2P_MHz"], 1.0), val


def test_hyperfine_21cm(spec):
    val = spec.results["Hyperfine_21cm_MHz"]
    # Fermi-contact leading form lands within ~0.5% of 1420.4 MHz
    assert _close(val, NIST["Hyperfine_21cm_MHz"], 1.0), val


# --- fine structure --------------------------------------------------------

def test_fine_structure_sign_and_scale():
    s = AtomicSpectroscopy()
    dE_2P32 = s.fine_structure(2, 1, 1.5)
    dE_2P12 = s.fine_structure(2, 1, 0.5)
    # 2P_{3/2} lies above 2P_{1/2}
    assert dE_2P32 > dE_2P12
    # alpha^2 suppression ~ 5e-5 of Rydberg
    from scipy import constants as C
    Ry = C.Rydberg * C.h * C.c
    assert abs(dE_2P32) < 1e-3 * Ry


# --- helium ---------------------------------------------------------------

def test_helium_singlet_triplet(spec):
    val = spec.results["He_1s2s_split_eV"]
    assert _close(val, NIST["He_1s2s_S_T_split_eV"], 1.0)


# --- alkalis --------------------------------------------------------------

@pytest.mark.parametrize("el,key", [
    ("Na", "Na_D2_nm"),
    ("K",  "K_D2_nm"),
    ("Rb", "Rb_D2_nm"),
    ("Cs", "Cs_D2_nm"),
])
def test_alkali_d_lines(spec, el, key):
    val = spec.results[key]
    assert _close(val, NIST[key], 1.0), (el, val, NIST[key])


# --- compare/report -------------------------------------------------------

def test_compare_returns_full_table(spec):
    tbl = spec.compare_to_NIST()
    assert "Hyperfine_21cm_MHz" in tbl
    for v in tbl.values():
        assert {"calc", "nist", "abs_err", "rel_err_pct"} <= set(v)


def test_report_string(spec):
    s = spec.report()
    assert "Substrate Atomic Spectroscopy" in s
    assert "worst" in s


def test_unknown_alkali_raises():
    with pytest.raises(ValueError):
        AtomicSpectroscopy().alkali_d_line("Fr")
