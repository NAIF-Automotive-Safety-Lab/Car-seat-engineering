"""Tests for the paired CMB GEOMETRY + SIMULATOR module.

Covers the four required claims:
    (a) first acoustic peak l ~ 220
    (b) blackbody temperature T_CMB = 2.7255 K
    (c) E-mode amplitude (~5% of TT, first EE peak l ~ 140)
    (d) B-mode tensor ratio r ~ 0  (substrate prediction)
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.cmb_paired import (
    CMBGeometry,
    CMBSimulator,
    LENSING_B_PEAK_L,
    PLANCK_EE_FIRST_PEAK,
    PLANCK_TT_PEAKS,
    T_CMB,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom():
    return CMBGeometry(n_theta=60, n_phi=120, l_max=15)


@pytest.fixture
def sim(geom):
    return CMBSimulator(geometry=geom)


# --------------------------------------------------------------------- #
# (b) Blackbody temperature                                             #
# --------------------------------------------------------------------- #

def test_temperature_matches_2p7255(sim):
    assert sim.temperature() == pytest.approx(2.7255, abs=1e-4)


def test_geometry_horizon_temperature(geom):
    assert geom.horizon_temperature() == pytest.approx(T_CMB, abs=1e-4)


# --------------------------------------------------------------------- #
# (a) First acoustic peak                                               #
# --------------------------------------------------------------------- #

def test_first_acoustic_peak_near_220(sim):
    l_peak = sim.first_acoustic_peak()
    assert 195.0 < l_peak < 245.0   # within ~10% of 220


def test_first_five_peak_positions_match_planck(sim):
    peaks = sim.peak_positions()
    assert len(peaks) == 5
    for p, o in zip(peaks, PLANCK_TT_PEAKS[:5]):
        assert abs(p - o) / o < 0.05


def test_detect_peaks_returns_at_least_first_peak(sim):
    detected = sim.detect_peaks(5)
    assert len(detected) >= 1
    assert 195.0 < detected[0] < 245.0


def test_tt_spectrum_shape_positive(sim):
    ells = np.arange(2, 2000)
    d = sim.tt_spectrum(ells)
    assert d.shape == ells.shape
    assert np.all(d >= 0.0)


def test_tt_spectrum_silk_damped(sim):
    ells = np.arange(2, 3000)
    d = sim.tt_spectrum(ells)
    first = d[(ells > 180) & (ells < 260)].max()
    tail = d[ells > 2500].max()
    assert tail < 0.2 * first


# --------------------------------------------------------------------- #
# (c) E-mode polarization amplitude                                     #
# --------------------------------------------------------------------- #

def test_ee_first_peak_near_140(sim):
    l_peak = sim.first_ee_peak()
    assert abs(l_peak - PLANCK_EE_FIRST_PEAK) / PLANCK_EE_FIRST_PEAK < 0.20


def test_ee_amplitude_is_a_few_percent_of_tt(sim):
    r = sim.ee_amplitude_ratio()
    # Canonically EE peak D_l ~ 45 uK^2 vs TT ~ 5750 -> ~ 0.008
    # Allow a generous range to keep the test robust.
    assert 0.001 < r < 0.10


def test_ee_spectrum_positive(sim):
    ells = np.arange(2, 2000)
    ee = sim.ee_spectrum(ells)
    assert ee.shape == ells.shape
    assert np.all(ee >= 0.0)


# --------------------------------------------------------------------- #
# (d) B-mode: substrate predicts r = 0                                  #
# --------------------------------------------------------------------- #

def test_b_mode_r_tensor_is_zero(sim):
    assert sim.tensor_to_scalar_ratio() == pytest.approx(0.0, abs=1e-12)


def test_b_mode_summary_substrate_zero(sim):
    s = sim.b_mode_summary()
    assert s["r_tensor_substrate_prediction"] == 0.0
    assert s["compatible_with_observations"] is True
    assert s["r_tensor_observed_upper_bound"] > 0.0
    assert s["lensing_B_peak_l"] == LENSING_B_PEAK_L


def test_b_mode_spectrum_only_lensing_when_r_zero(sim):
    ells = np.arange(2, 2000)
    bb = sim.bb_spectrum(ells)
    # Lensing bump dominates around l ~ 1000
    l_peak = ells[np.argmax(bb)]
    assert 700 < l_peak < 1300
    # Low-l (l < 200, where primordial would dominate) is well below the
    # lensing peak when r = 0
    low = bb[ells < 200].max()
    peak = bb.max()
    assert low < 0.5 * peak


def test_b_mode_spectrum_with_nonzero_r_has_low_l_bump():
    sim_r = CMBSimulator(r_tensor=0.1)
    ells = np.arange(2, 2000)
    bb = sim_r.bb_spectrum(ells)
    # Now there is a low-l contribution at l ~ a few-tens
    assert bb[ells == 5][0] > 0.0


# --------------------------------------------------------------------- #
# Geometry: anisotropy on the 2-sphere                                  #
# --------------------------------------------------------------------- #

def test_geometry_grid_shape(geom):
    delta = geom.temperature_field_uK()
    assert delta.shape == (geom.n_theta, geom.n_phi)


def test_geometry_anisotropy_rms_close_to_target(geom):
    delta = geom.temperature_field_uK()
    rms = float(np.std(delta))
    # The toy generator targets ~ 110 uK; allow generous tolerance.
    assert 50.0 < rms < 250.0


def test_geometry_polarization_QU(geom):
    Q, U = geom.polarization_QU()
    assert Q.shape == (geom.n_theta, geom.n_phi)
    assert U.shape == (geom.n_theta, geom.n_phi)
    # Non-trivial polarization amplitude exists
    assert np.std(np.hypot(Q, U)) > 0.0


def test_geometry_substrate_signature(geom):
    sig = geom.substrate_signature()
    for key in ("interpretation", "no_inflation",
                "primordial_B_mode", "low_l_suppression_factor",
                "low_l_feature_multipole", "lambda_qcd_GeV"):
        assert key in sig
    assert sig["primordial_B_mode"] == 0.0
    assert sig["no_inflation"] is True


# --------------------------------------------------------------------- #
# Substrate low-l anomaly                                               #
# --------------------------------------------------------------------- #

def test_low_l_anomaly_keys(sim):
    a = sim.low_l_anomaly()
    for key in ("low_l_suppression_factor", "feature_multipole",
                "lambda_qcd_GeV", "matches_planck_anomaly"):
        assert key in a
    assert 0.0 < a["low_l_suppression_factor"] <= 1.0


def test_low_l_suppression_visible_in_tt(sim):
    """At l=2 the substrate suppression must reduce D_l below the
    un-suppressed plateau value (sw_plateau_uK2)."""
    d_at_2 = sim.tt_spectrum(np.array([2]))[0]
    # Without suppression the plateau at l=2 would be sw_plateau_uK2.
    # The substrate suppression (Lambda_QCD imprint) lowers this.
    assert d_at_2 < sim.sw_plateau_uK2
