"""Tests for src.stiff_medium.saturation_horizon_geometry.

Verifies the four claims required by the saturation-horizon construction:

  (a) Cone is vertical (tilt = 0°) at σ = 0
  (b) Cone tilt is exactly 90° at σ = σ_max = 1/2
  (c) Schwarzschild radius r_s = 2 G M / c² is recovered as σ(r_s) = 1/2
  (d) Crack-tip stress is capped at σ_max (finite at every r > 0)

Plus a few sanity checks on monotonicity, potential growth, and the
HorizonSimulator photon trajectory.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.saturation_horizon_geometry import (
    SIGMA_MAX,
    G_NEWTON,
    C_LIGHT,
    SaturationHorizonGeometry,
    CrackTipGeometry,
    HorizonSimulator,
    cone_family,
    crack_stress_curve,
    potential_curve,
)


EPS_ANGLE = 1e-9
EPS_NUM = 1e-9


# ---------------------------------------------------------------------------
# (a) Cone vertical at σ = 0
# ---------------------------------------------------------------------------

def test_cone_vertical_at_zero_sigma():
    g = SaturationHorizonGeometry()
    tilt_rad = float(g.cone_tilt_angle(0.0))
    tilt_deg = float(g.cone_tilt_degrees(0.0))
    assert tilt_rad == pytest.approx(0.0, abs=EPS_ANGLE)
    assert tilt_deg == pytest.approx(0.0, abs=EPS_ANGLE)

    # Walls at σ=0 are degenerate (both lie along the t-axis)
    walls = g.cone_walls(0.0, r0=0.0, t0=0.0, half_height=1.0)
    # Right and left wall tips should both have r=0 (no spatial spread)
    assert walls["right"][1, 0] == pytest.approx(0.0, abs=EPS_NUM)
    assert walls["left"][1, 0] == pytest.approx(0.0, abs=EPS_NUM)
    # And ct = half_height (1.0)
    assert walls["right"][1, 1] == pytest.approx(1.0, abs=EPS_NUM)
    assert walls["left"][1, 1] == pytest.approx(1.0, abs=EPS_NUM)


# ---------------------------------------------------------------------------
# (b) Cone tilt is EXACTLY 90° at σ = 1/2
# ---------------------------------------------------------------------------

def test_cone_tilt_exactly_90_at_sigma_max():
    g = SaturationHorizonGeometry()
    tilt_rad = float(g.cone_tilt_angle(SIGMA_MAX))
    tilt_deg = float(g.cone_tilt_degrees(SIGMA_MAX))
    assert tilt_rad == pytest.approx(np.pi / 2, abs=EPS_ANGLE)
    assert tilt_deg == pytest.approx(90.0, abs=1e-6)

    # Walls at σ = σ_max should lie along the spatial (r) axis:
    # half-tilt = 45°, so each wall has equal r and ct components
    walls = g.cone_walls(SIGMA_MAX, r0=0.0, t0=0.0, half_height=1.0)
    # Right wall: dr = sin(45°), dct = cos(45°), each = √2 / 2
    assert walls["right"][1, 0] == pytest.approx(
        np.sin(np.pi / 4), abs=1e-6)
    assert walls["right"][1, 1] == pytest.approx(
        np.cos(np.pi / 4), abs=1e-6)
    # Left wall mirror
    assert walls["left"][1, 0] == pytest.approx(
        -np.sin(np.pi / 4), abs=1e-6)


def test_cone_tilt_45_degrees_at_sigma_eighth():
    """σ = σ_max / 4 = 1/8 should give 45° tilt by the 2·arctan(√(σ/σ_max))
    construction."""
    g = SaturationHorizonGeometry()
    sigma_eighth = SIGMA_MAX / 4.0   # = 1/8
    tilt = float(g.cone_tilt_degrees(sigma_eighth))
    # 2 · arctan(sqrt(1/4)) = 2 · arctan(1/2)
    expected = 2.0 * np.degrees(np.arctan(0.5))
    assert tilt == pytest.approx(expected, abs=1e-9)


def test_cone_tilt_monotone_increasing():
    g = SaturationHorizonGeometry()
    sigmas = np.linspace(0.0, SIGMA_MAX, 50)
    tilts = g.cone_tilt_degrees(sigmas)
    diffs = np.diff(tilts)
    assert np.all(diffs >= -1e-12)   # non-decreasing
    assert tilts[0] < tilts[-1]      # actually grows


# ---------------------------------------------------------------------------
# (c) Schwarzschild radius is recovered: σ(r_s) = 1/2
# ---------------------------------------------------------------------------

def test_schwarzschild_radius_recovered():
    g = SaturationHorizonGeometry()
    M = 1.989e30                      # solar mass (kg)
    r_s = g.schwarzschild_radius(M)
    expected = 2.0 * G_NEWTON * M / C_LIGHT ** 2
    assert r_s == pytest.approx(expected, rel=1e-12)
    # Solar Schwarzschild radius ≈ 2.95 km
    assert 2_900.0 < r_s < 3_000.0

    # σ(r_s) must equal σ_max = 1/2
    sigma_at_rs = float(g.sigma_at_radius(r_s, M))
    assert sigma_at_rs == pytest.approx(0.5, rel=1e-12)


def test_horizon_check_dictionary():
    g = SaturationHorizonGeometry()
    M = 1.0e30
    r_s = g.schwarzschild_radius(M)
    chk = g.horizon_check(M, r=r_s)
    assert chk["horizon_reproduced"] is True
    assert chk["sigma_at_rs"] == pytest.approx(0.5, abs=1e-12)
    assert chk["tilt_at_rs_deg"] == pytest.approx(90.0, abs=1e-6)


def test_clock_rate_vanishes_at_horizon():
    g = SaturationHorizonGeometry()
    f0 = float(g.clock_rate(0.0))
    f_h = float(g.clock_rate(SIGMA_MAX))
    assert f0 == pytest.approx(1.0, abs=1e-12)
    assert f_h == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# (d) Crack-tip stress capped at finite value
# ---------------------------------------------------------------------------

def test_crack_tip_stress_finite_everywhere():
    crack = CrackTipGeometry(a=1.0, sigma_inf=0.1)
    r = np.linspace(1e-6, 2.0, 500)
    stress = crack.stress_capped(r)
    # Every value finite
    assert np.all(np.isfinite(stress))
    # Bounded above by σ_max
    assert np.max(stress) <= SIGMA_MAX + 1e-12
    # Approaches σ_max near the tip
    assert stress[0] == pytest.approx(SIGMA_MAX, abs=1e-12)
    # And LEFM diverges to compare:
    classical = crack.stress_LEFM(np.array([1e-12]))
    assert classical[0] > SIGMA_MAX * 100   # massively over the cap


def test_crack_process_zone_radius_formula():
    crack = CrackTipGeometry(a=2.0, sigma_inf=0.1)
    r_p = crack.process_zone_radius()
    # r_p = (a/2) (σ_∞/σ_max)² = 1.0 · (0.1/0.5)² = 0.04
    assert r_p == pytest.approx(0.04, abs=1e-12)
    # And at r = r_p the LEFM stress equals σ_max exactly
    s_at_rp = float(crack.stress_LEFM(np.array([r_p]))[0])
    assert s_at_rp == pytest.approx(SIGMA_MAX, rel=1e-12)


def test_crack_peak_stress_is_sigma_max():
    crack = CrackTipGeometry(a=1.0, sigma_inf=0.05)
    assert crack.peak_stress() == pytest.approx(SIGMA_MAX)


def test_crack_curve_helper_outputs():
    data = crack_stress_curve(a=1.0, sigma_inf=0.1,
                              r_min=1e-4, r_max=2.0, n=200)
    for k in ("r", "stress_LEFM", "stress_capped",
              "process_zone_radius", "peak_stress"):
        assert k in data
    assert np.max(data["stress_capped"]) <= SIGMA_MAX + 1e-12


# ---------------------------------------------------------------------------
# Substrate potential V(σ) finite below cap, monotone, divergent at cap
# ---------------------------------------------------------------------------

def test_potential_zero_at_zero():
    g = SaturationHorizonGeometry()
    assert float(g.substrate_potential(0.0)) == pytest.approx(0.0, abs=1e-12)


def test_potential_monotone_increasing_in_sigma():
    sigma, V = potential_curve(n=200, K=1.0, sigma_clip=0.499)
    diffs = np.diff(V)
    assert np.all(diffs >= -1e-12)
    assert V[0] < V[-1]


def test_potential_capped_value_finite():
    g = SaturationHorizonGeometry()
    V_near_cap = g.potential_capped_value(sigma_clip=0.499, K=1.0)
    # Finite (not inf)
    assert np.isfinite(V_near_cap)
    # And large compared to V at small σ
    V_small = float(g.substrate_potential(0.05))
    assert V_near_cap > V_small * 10


# ---------------------------------------------------------------------------
# Cone-family helper produces n cones with consistent tilts
# ---------------------------------------------------------------------------

def test_cone_family_produces_correct_tilts():
    sigmas = np.array([0.0, 0.125, 0.25, 0.375, 0.5])
    cones = cone_family(sigmas, half_height=1.0)
    assert len(cones) == 5
    g = SaturationHorizonGeometry()
    for c, s in zip(cones, sigmas):
        expected_deg = float(g.cone_tilt_degrees(float(s)))
        assert c["tilt_degrees"] == pytest.approx(expected_deg, abs=1e-9)


# ---------------------------------------------------------------------------
# HorizonSimulator end-to-end: photon slows as it approaches horizon
# ---------------------------------------------------------------------------

def test_horizon_simulator_runs():
    # r_min = r_s exactly so the inner-boundary cone is at the horizon
    M = 1.0e30
    r_s = 2.0 * G_NEWTON * M / C_LIGHT ** 2
    sim = HorizonSimulator(M=M, r_min=r_s, r_max=50.0 * r_s,
                           n_r=100, n_t=200)
    res = sim.run()
    assert res.photon_path_r.size > 1
    # Photon must move outward (monotone non-decreasing radius)
    diffs = np.diff(res.photon_path_r)
    assert np.all(diffs >= -1e-3)
    # Cone tilt at r_s should be exactly 90° (horizon)
    assert res.cone_tilt_deg[0] == pytest.approx(90.0, abs=1e-6)
    # And tilt should monotonically decrease as r grows (outside horizon)
    assert res.cone_tilt_deg[-1] < res.cone_tilt_deg[0]


def test_horizon_simulator_metadata():
    sim = HorizonSimulator(M=1.0e30, r_min=5e3, r_max=5e4,
                           n_r=100, n_t=200)
    res = sim.run()
    md = res.metadata
    assert md["M_kg"] == 1.0e30
    assert md["r_schwarzschild_m"] == pytest.approx(
        2.0 * G_NEWTON * 1.0e30 / C_LIGHT ** 2, rel=1e-12)
    assert md["n_steps"] >= 2
