"""Tests for EM radiation simulator (substrate transverse modes)."""

import numpy as np
import pytest

from stiff_medium.em_radiation import (
    EMRadiation,
    default_electron_oscillator,
    EPS0,
    C,
    E_CHARGE,
    M_E,
    PI,
)


# ---------------------------------------------------------------------------
# (a) Larmor power for the test electron
# ---------------------------------------------------------------------------
def test_larmor_power_electron_10nm_1e15():
    em = default_electron_oscillator()
    a_peak = em.peak_acceleration()
    assert a_peak == pytest.approx(1.0e-8 * (1.0e15) ** 2)

    P_peak = EMRadiation.larmor_power(E_CHARGE, a_peak)
    # sanity: P > 0 and order of magnitude consistent with formula
    expected = (2.0 / 3.0) * (E_CHARGE ** 2 * a_peak ** 2) / (
        4.0 * PI * EPS0 * C ** 3
    )
    assert P_peak == pytest.approx(expected, rel=1e-12)
    assert P_peak > 0.0


def test_time_averaged_power_is_half_peak():
    em = default_electron_oscillator()
    a_peak = em.peak_acceleration()
    P_peak = EMRadiation.larmor_power(E_CHARGE, a_peak)
    assert em.time_averaged_power() == pytest.approx(0.5 * P_peak, rel=1e-12)


# ---------------------------------------------------------------------------
# (b) Dipole pattern
# ---------------------------------------------------------------------------
def test_dipole_pattern_max_at_pi_over_2_zero_on_axis():
    p_perp = EMRadiation.dipole_radiation_pattern(np.pi / 2)
    p_axis = EMRadiation.dipole_radiation_pattern(0.0)
    assert p_perp == pytest.approx(1.0, rel=1e-12)
    assert p_axis == pytest.approx(0.0, abs=1e-15)
    # ratio is the canonical sin^2 pattern signature
    assert p_perp > p_axis


def test_dipole_pattern_array():
    th = np.linspace(0, np.pi, 11)
    p = EMRadiation.dipole_radiation_pattern(th)
    assert p.shape == th.shape
    assert np.all(p >= 0)
    assert np.all(p <= 1.0 + 1e-12)


# ---------------------------------------------------------------------------
# (c) Substrate-field calculation reproduces Larmor
# ---------------------------------------------------------------------------
def test_substrate_calculation_agrees_with_larmor_to_5pct():
    em = default_electron_oscillator()
    out = em.verify_against_larmor(R=1.0, tol=0.05)
    assert out["passes"], (
        f"substrate {out['P_substrate']:.3e} vs "
        f"larmor {out['P_larmor']:.3e} (rel_err {out['rel_err']:.2%})"
    )
    # In fact they should agree to machine precision since the
    # angular integral gives 8 pi / 3 analytically.
    assert out["rel_err"] < 1e-6


def test_substrate_independent_of_observation_radius():
    em = default_electron_oscillator()
    P1 = em.substrate_radiation_power(R=1.0)
    P2 = em.substrate_radiation_power(R=137.0)
    assert P1 == pytest.approx(P2, rel=1e-10)


# ---------------------------------------------------------------------------
# (d) Cyclotron formula matches the textbook result
# ---------------------------------------------------------------------------
def test_cyclotron_matches_textbook():
    B = 1.0  # tesla
    v = 1.0e6  # m/s, non-relativistic
    P = EMRadiation.cyclotron_radiation(B, v)

    # Direct textbook form: P = q^4 B^2 v^2 / (6 pi eps0 m^2 c^3)
    P_book = (E_CHARGE ** 4 * B ** 2 * v ** 2) / (
        6.0 * PI * EPS0 * M_E ** 2 * C ** 3
    )
    assert P == pytest.approx(P_book, rel=1e-12)


# ---------------------------------------------------------------------------
# Bonus: quadrupole, synchrotron, retarded field, bremsstrahlung
# ---------------------------------------------------------------------------
def test_quadrupole_pattern_zero_on_axis_and_equator():
    Q = np.diag([1.0, 1.0, -2.0]) * 1e-30  # traceless
    P, pattern = EMRadiation.quadrupole_radiation(
        Q, omega=1e15, theta=np.array([0.0, np.pi / 4, np.pi / 2])
    )
    assert P > 0
    # sin^2(2 theta) zeros at 0 and pi/2, max at pi/4
    assert pattern[0] == pytest.approx(0.0, abs=1e-15)
    assert pattern[2] == pytest.approx(0.0, abs=1e-15)
    assert pattern[1] == pytest.approx(1.0, rel=1e-12)


def test_synchrotron_critical_freq_scales_as_gamma_cubed():
    B = 1.0
    w1 = EMRadiation.synchrotron_critical_freq(1.0, B)
    w10 = EMRadiation.synchrotron_critical_freq(10.0, B)
    assert w10 / w1 == pytest.approx(1000.0, rel=1e-12)
    # absolute value
    assert w1 == pytest.approx(1.5 * E_CHARGE * B / M_E, rel=1e-12)


def test_bremsstrahlung_scales_as_Z_squared_and_dv_squared():
    base = EMRadiation.bremsstrahlung(0.0, 1e6, Z=1.0)
    z2 = EMRadiation.bremsstrahlung(0.0, 1e6, Z=2.0)
    big_dv = EMRadiation.bremsstrahlung(0.0, 2e6, Z=1.0)
    assert z2 / base == pytest.approx(4.0, rel=1e-12)
    assert big_dv / base == pytest.approx(4.0, rel=1e-12)


def test_retarded_E_field_transverse_to_n():
    em = default_electron_oscillator()
    r_obs = np.array([1.0, 0.0, 0.5])
    E = em.retarded_E_field(r_obs, t=1e-16)
    n = r_obs / np.linalg.norm(r_obs)
    # Far-field radiation E must be perpendicular to n
    assert abs(np.dot(E, n)) < 1e-10 * np.linalg.norm(E)
