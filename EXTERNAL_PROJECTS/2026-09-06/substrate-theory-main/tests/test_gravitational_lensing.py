"""Tests for the gravitational lensing simulator."""

import numpy as np
import pytest
from scipy import constants as const

from src.stiff_medium.gravitational_lensing import (
    ARCSEC,
    KPC,
    MAS,
    MPC,
    M_SUN,
    GravitationalLensing,
    lmc_microlensing_einstein_radius,
    sun_grazing_deflection,
    to_arcsec,
    to_mas,
)


lens = GravitationalLensing()


# ---------- (a) Sun grazing deflection 1.75" ----------
def test_sun_grazing_deflection_arcsec():
    alpha = sun_grazing_deflection()
    assert 1.70 < alpha < 1.80, f"Sun deflection {alpha:.3f}\" off 1.75\""


def test_deflection_formula_matches_einstein():
    M = M_SUN
    b = 6.957e8  # solar radius
    alpha = lens.deflection_angle(M, b)
    expected = 4 * const.G * M / (b * const.c ** 2)
    assert np.isclose(alpha, expected, rtol=1e-12)


# ---------- (b) LMC microlensing Einstein radius ~ 1 mas ----------
def test_lmc_einstein_radius_mas():
    theta_E_mas = lmc_microlensing_einstein_radius(
        M_lens_solar=0.5, D_lens_kpc=10.0, D_source_kpc=50.0
    )
    # Standard textbook value is ~0.5-1 mas for 0.5 Msun midway to LMC
    assert 0.3 < theta_E_mas < 2.0, f"LMC theta_E={theta_E_mas:.3f} mas"


def test_einstein_radius_scales_as_sqrt_mass():
    D_L, D_S = 10 * KPC, 50 * KPC
    D_LS = D_S - D_L
    th1 = lens.einstein_radius(M_SUN, D_L, D_S, D_LS)
    th2 = lens.einstein_radius(4 * M_SUN, D_L, D_S, D_LS)
    assert np.isclose(th2 / th1, 2.0, rtol=1e-10)


# ---------- (c) Microlensing magnification at u=0.1 ~ 10 ----------
def test_microlensing_peak_magnification():
    A = lens.microlensing_magnification(0.1)
    # Exact: (0.01+2)/(0.1*sqrt(4.01)) = 2.01/0.2002... ~ 10.04
    assert 9.9 < A < 10.2, f"A(0.1)={A:.3f}"


def test_microlensing_large_u_unity():
    A = lens.microlensing_magnification(10.0)
    assert np.isclose(A, 1.0, atol=1e-2)


def test_microlensing_diverges_at_zero():
    A_small = lens.microlensing_magnification(1e-6)
    assert A_small > 1e5


# ---------- strong lensing two images ----------
def test_strong_lens_two_images_sum_and_product():
    theta_E = 1.0 * ARCSEC
    beta = 0.3 * ARCSEC
    tp, tm = lens.strong_lensing_images(beta, theta_E)
    # Vieta: tp+tm = beta, tp*tm = -theta_E^2
    assert np.isclose(tp + tm, beta, rtol=1e-12)
    assert np.isclose(tp * tm, -theta_E ** 2, rtol=1e-12)
    # Two images on opposite sides
    assert tp > 0 and tm < 0


def test_strong_lens_einstein_ring_when_aligned():
    theta_E = 1.0 * ARCSEC
    tp, tm = lens.strong_lensing_images(0.0, theta_E)
    assert np.isclose(tp, theta_E)
    assert np.isclose(tm, -theta_E)


# ---------- (d) Weak lensing shear ----------
def test_weak_lensing_shear_cluster():
    M = 1e14 * M_SUN
    R = 100 * KPC
    # use distance-based critical surface density: lens at z~0.3, source z~1
    D_L = 900 * MPC
    D_S = 1650 * MPC
    D_LS = 1100 * MPC
    gamma = lens.weak_lensing_shear(M, R, D_lens=D_L, D_source=D_S, D_LS=D_LS)
    # Point-mass shear at 100 kpc from 1e14 Msun is order unity — extended
    # NFW would be smaller, but the formula is mathematically correct.
    assert 1e-3 < gamma < 5.0, f"gamma_t = {gamma:.4f}"


def test_weak_lensing_shear_falls_as_inv_r2():
    M = 1e14 * M_SUN
    D_L, D_S, D_LS = 900 * MPC, 1650 * MPC, 1100 * MPC
    g1 = lens.weak_lensing_shear(M, 100 * KPC, D_lens=D_L, D_source=D_S, D_LS=D_LS)
    g2 = lens.weak_lensing_shear(M, 200 * KPC, D_lens=D_L, D_source=D_S, D_LS=D_LS)
    assert np.isclose(g1 / g2, 4.0, rtol=1e-10)


# ---------- (e) Time delay positive for trailing image ----------
def test_time_delay_positive_for_trailing_image():
    theta_E = 1.0 * ARCSEC
    beta = 0.3 * ARCSEC
    tp, tm = lens.strong_lensing_images(beta, theta_E)
    H_0 = 70e3 / MPC  # 70 km/s/Mpc in 1/s
    # tp is the brighter, leading image (arrives first); tm trails.
    # Convention: dt(theta1, theta2) returned with theta1=leading => positive
    dt = lens.time_delay(tp, tm, beta, theta_E, z_lens=0.5, H_0=H_0)
    assert dt > 0, f"dt = {dt} s should be positive for trailing image"


# ---------- cluster arc geometry sanity ----------
def test_cluster_arc_returns_dict():
    out = lens.cluster_lensing_arc(1e14 * M_SUN, 50 * KPC, 0.0)
    assert "theta_E_arcsec" in out
    assert out["theta_E_arcsec"] > 0


def test_schwarzschild_radius_solar():
    rs = lens.schwarzschild_radius(M_SUN)
    assert np.isclose(rs, 2953.0, atol=10.0)
