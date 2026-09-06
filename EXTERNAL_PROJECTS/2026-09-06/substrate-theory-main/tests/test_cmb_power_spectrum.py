"""Tests for CMBPowerSpectrum (substrate framework predictor)."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.cmb_power_spectrum import (
    CMBPowerSpectrum,
    PLANCK_TT_PEAKS,
    PLANCK_EE_FIRST_PEAK,
    T_CMB,
    planck_blackbody,
)


@pytest.fixture
def cmb():
    return CMBPowerSpectrum()


# --------------------------------------------------------------------- #
# (a) Temperature                                                       #
# --------------------------------------------------------------------- #

def test_temperature_matches_observed(cmb):
    assert cmb.temperature() == pytest.approx(2.7255, abs=1e-3)


def test_blackbody_peak_frequency():
    """Wien displacement gives ν_peak ≈ 160 GHz at T = 2.725 K."""
    nu = np.linspace(1e9, 1e12, 5000)
    B = planck_blackbody(nu, T=T_CMB)
    nu_peak_ghz = nu[np.argmax(B)] / 1e9
    assert 140.0 < nu_peak_ghz < 180.0


# --------------------------------------------------------------------- #
# (b) Acoustic peak positions                                           #
# --------------------------------------------------------------------- #

def test_first_four_peak_positions_within_5pct(cmb):
    pred = cmb.peak_positions()
    obs  = PLANCK_TT_PEAKS[:4]
    assert len(pred) == 4
    for p, o in zip(pred, obs):
        assert abs(p - o) / o < 0.05


# --------------------------------------------------------------------- #
# (c) Peak height ratio                                                 #
# --------------------------------------------------------------------- #

def test_second_over_first_ratio_within_30pct(cmb):
    h = cmb.peak_heights()
    observed_ratio = 2500.0 / 5750.0  # ≈ 0.435
    pred_ratio = h["second_over_first"]
    rel_err = abs(pred_ratio - observed_ratio) / observed_ratio
    assert rel_err < 0.30


def test_peak_heights_keys(cmb):
    h = cmb.peak_heights()
    for k in ("first", "second_over_first",
              "third_over_first", "fourth_over_first"):
        assert k in h
    assert h["first"] == 1.0


# --------------------------------------------------------------------- #
# (d) E-mode peak                                                       #
# --------------------------------------------------------------------- #

def test_e_mode_first_peak_location(cmb):
    """First EE peak should be near ℓ ≈ 140 (Planck observation)."""
    l_first = cmb.polarization_E_mode()
    assert abs(l_first - PLANCK_EE_FIRST_PEAK) / PLANCK_EE_FIRST_PEAK < 0.10


def test_e_mode_spectrum_has_peak_near_140(cmb):
    ells = np.arange(2, 1500)
    ee = cmb.polarization_E_mode(ells)
    l_peak = ells[np.argmax(ee)]
    assert 100 < l_peak < 200


# --------------------------------------------------------------------- #
# TT spectrum sanity                                                    #
# --------------------------------------------------------------------- #

def test_tt_spectrum_shape_and_first_peak(cmb):
    ells = np.arange(2, 2000)
    d_ell = cmb.tt_spectrum(ells)
    assert d_ell.shape == ells.shape
    assert np.all(d_ell >= 0)
    # First peak in 100 < ℓ < 350
    region = (ells > 100) & (ells < 350)
    l_peak = ells[region][np.argmax(d_ell[region])]
    assert 180 < l_peak < 260


def test_tt_spectrum_silk_damping(cmb):
    """High-ℓ tail must be suppressed below first peak amplitude."""
    ells = np.arange(2, 3000)
    d_ell = cmb.tt_spectrum(ells)
    first = d_ell[(ells > 180) & (ells < 260)].max()
    tail  = d_ell[ells > 2500].max()
    assert tail < 0.2 * first


# --------------------------------------------------------------------- #
# B-mode                                                                #
# --------------------------------------------------------------------- #

def test_b_mode_substrate_predicts_no_primordial(cmb):
    summary = cmb.polarization_B_mode()
    assert summary["r_tensor_substrate_prediction"] == 0.0
    assert summary["compatible_with_observations"] is True
    assert summary["r_tensor_observed_upper_bound"] > 0


def test_b_mode_spectrum_positive(cmb):
    ells = np.arange(2, 2000)
    bb = cmb.polarization_B_mode(ells)
    assert np.all(bb >= 0)


# --------------------------------------------------------------------- #
# Substrate signature                                                   #
# --------------------------------------------------------------------- #

def test_substrate_signature_keys(cmb):
    sig = cmb.substrate_horizon_flux_signature()
    for k in ("interpretation",
              "high_l_predictions",
              "low_l_feature_multipole",
              "low_l_suppression_factor",
              "lambda_qcd_GeV",
              "primordial_B_mode_prediction"):
        assert k in sig
    assert sig["primordial_B_mode_prediction"] == 0.0
    assert 0.0 < sig["low_l_suppression_factor"] <= 1.0


# --------------------------------------------------------------------- #
# Planck comparison                                                     #
# --------------------------------------------------------------------- #

def test_compare_to_planck_passes(cmb):
    cmp = cmb.compare_to_planck_2018()
    assert cmp["passes_within_5pct"] is True
    assert abs(cmp["T_fractional_residual"]) < 1e-3
    assert cmp["max_abs_peak_residual"] < 0.05
