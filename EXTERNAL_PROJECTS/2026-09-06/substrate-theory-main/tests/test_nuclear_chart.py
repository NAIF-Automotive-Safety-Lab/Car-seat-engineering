"""Tests for the substrate-flavoured nuclear chart calculator."""

from __future__ import annotations

import numpy as np
import pytest

from stiff_medium.nuclear_chart import (
    A_V_BARE_MEV,
    EPS_FACE_MEV,
    LAMBDA_QCD_MEV,
    MAGIC_NUMBERS,
    REFERENCE_BE,
    NuclearChart,
    calibrated_chart,
)


# ---------------------------------------------------------------------------
# Substrate constant sanity checks
# ---------------------------------------------------------------------------

def test_eps_face_value():
    assert EPS_FACE_MEV == pytest.approx(LAMBDA_QCD_MEV / 90.0, rel=1e-12)
    assert EPS_FACE_MEV == pytest.approx(2.2222222, rel=1e-5)


def test_a_v_bare_six_faces():
    assert A_V_BARE_MEV == pytest.approx(6.0 * EPS_FACE_MEV, rel=1e-12)


# ---------------------------------------------------------------------------
# Calibration must succeed and produce sane coefficients
# ---------------------------------------------------------------------------

def test_calibration_runs_and_yields_sane_coefficients():
    chart = calibrated_chart()
    c = chart.coeffs
    # Empirical SEMF ranges; allow generous bounds
    assert 14.0 < c.a_v < 17.0
    assert 12.0 < c.a_s < 22.0
    assert 0.5 < c.a_C < 1.0
    assert 18.0 < c.a_a < 28.0
    assert 5.0 < c.a_p < 18.0


# ---------------------------------------------------------------------------
# Required test cases from the spec
# ---------------------------------------------------------------------------

def test_deuteron_binding_energy_in_range():
    chart = calibrated_chart()
    be = chart.binding_energy(1, 1)
    # Spec: 1-3 MeV (observed 2.22 MeV).
    assert 1.0 <= be <= 3.0


def test_fe56_binding_per_nucleon_in_range():
    chart = calibrated_chart()
    bea = chart.binding_per_nucleon(26, 30)
    # Spec: 8.6-9.0 MeV/A (observed 8.79 MeV).
    assert 8.6 <= bea <= 9.0


def test_pb208_binding_total():
    chart = calibrated_chart()
    be = chart.binding_energy(82, 126)
    # Spec: 1620 +/- 50 MeV (observed 1636.4 MeV).
    assert 1570.0 <= be <= 1670.0


def test_magic_numbers_from_k4_shells_match_expected():
    chart = calibrated_chart()
    got = chart.magic_numbers_from_K4_shells()
    assert tuple(got) == MAGIC_NUMBERS


def test_valley_of_stability_matches_observation():
    chart = calibrated_chart()
    valley = dict(chart.valley_of_stability(A_max=210))
    # Spot-check well-known stable nuclei
    # A=4 -> Z=2 (alpha); A=12 -> Z=6 (C-12); A=56 -> Z=26 (Fe-56);
    # A=120 -> Z~50 (Sn-120); A=208 -> Z~82 (Pb-208)
    assert valley[4] == 2
    assert valley[12] == 6
    # A=40: Ar-40 (Z=18) and Ca-40 (Z=20) are both stable isobars
    assert abs(valley[40] - 20) <= 2
    assert abs(valley[56] - 26) <= 1
    assert abs(valley[120] - 50) <= 2
    assert abs(valley[208] - 82) <= 3


# ---------------------------------------------------------------------------
# Additional sanity tests
# ---------------------------------------------------------------------------

def test_alpha_binding_close_to_observed():
    chart = calibrated_chart()
    be = chart.binding_energy(2, 2)
    # Observed 28.30 MeV; SEMF struggles at light A, allow wider window.
    assert 20.0 <= be <= 35.0


def test_compare_to_ame_returns_residuals():
    chart = calibrated_chart()
    out = chart.compare_to_ame()
    # All reference nuclei present
    assert set(out.keys()) == set(REFERENCE_BE.keys())
    # Heavy nuclei should be predicted within ~5% relative
    rel_heavy = abs(out[(82, 126)]["rel_error"])
    assert rel_heavy < 0.05


def test_binding_energy_chart_shape_and_peak():
    chart = calibrated_chart()
    arr = chart.binding_energy_chart(Z_max=60, N_max=90)
    assert arr.shape == (61, 91)
    # Peak BE/A should be near Fe/Ni region (Z~26-28)
    peak_idx = np.unravel_index(np.argmax(arr), arr.shape)
    Z_peak, N_peak = peak_idx
    assert 22 <= Z_peak <= 32
    # And the peak value should be within sensible bounds
    assert 8.0 <= arr[peak_idx] <= 9.2


def test_mass_excess_monotone_for_isotopes():
    chart = calibrated_chart()
    # Mass excess should rise as you move far from the valley (Sn isotopes)
    me_stable = chart.mass_excess(50, 70)        # Sn-120 (stable)
    me_neutron_rich = chart.mass_excess(50, 90)  # Sn-140 (very n-rich)
    assert me_neutron_rich > me_stable


def test_pairing_signs():
    chart = calibrated_chart()
    # even-even should be > odd-A neighbours of same A-1 region (qualitative)
    be_ee = chart.binding_energy(20, 20)   # Ca-40, even-even
    be_oo = chart.binding_energy(21, 21)   # Sc-42, odd-odd  (smaller A=42 vs 40, careful)
    # Just check pairing term applies with correct sign for same nucleus
    assert chart._pairing(20, 20) > 0
    assert chart._pairing(21, 21) < 0
    assert chart._pairing(20, 21) == 0
