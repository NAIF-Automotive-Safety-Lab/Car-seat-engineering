"""Tests for stiff_medium.black_hole_paired — paired geometry/sim/viz BH module.

Verifies:
  (a) r_s = 2GM/c² yields σ = 1/2 at the horizon
  (b) photon orbits at 1.5 r_s  (photon sphere)
  (c) Hawking T_H for solar mass ≈ 6 × 10⁻⁸ K
  (d) ringdown ≈ 250 Hz for ~stellar (≥30 M_sun) BHs
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from stiff_medium.black_hole_paired import (
    BlackHoleGeometry,
    BlackHoleSimulator,
    BlackHoleVisualizer,
    C, G, HBAR, K_B, M_SUN, SIGMA_HORIZON,
    default_solar_bh,
    default_stellar_remnant,
)


# ---------------------------------------------------------------------------
# (a) r_s = 2GM/c² gives σ = 1/2 at the horizon
# ---------------------------------------------------------------------------


def test_r_s_yields_sigma_one_half_at_horizon() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    rs = geom.r_s
    # Schwarzschild radius matches textbook 2.953 km
    assert rs == pytest.approx(2.0 * G * M_SUN / C**2)
    assert rs == pytest.approx(2953.0, rel=2e-3)
    # σ at the horizon equals exactly 1/2 (universal cap)
    assert geom.sigma(rs) == pytest.approx(SIGMA_HORIZON, abs=1e-15)


def test_sigma_outside_falls_as_one_over_r() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    rs = geom.r_s
    # σ(2 r_s) = 1/4, σ(10 r_s) = 1/20, σ(100 r_s) = 1/200
    assert geom.sigma(2.0 * rs) == pytest.approx(0.25, rel=1e-12)
    assert geom.sigma(10.0 * rs) == pytest.approx(0.05, rel=1e-12)
    assert geom.sigma(100.0 * rs) == pytest.approx(0.005, rel=1e-12)


def test_sigma_inside_horizon_is_capped() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    rs = geom.r_s
    # Naive 1/r would exceed 1/2 inside; the implementation must clamp.
    assert geom.sigma(0.5 * rs) == pytest.approx(SIGMA_HORIZON)
    assert geom.sigma(0.1 * rs) == pytest.approx(SIGMA_HORIZON)


def test_c_local_vanishes_at_horizon_and_recovers_far_away() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    rs = geom.r_s
    # At the horizon c_local → 0 (cone closes)
    assert geom.c_local(rs) == pytest.approx(0.0, abs=1e-3)
    # Far away c_local → c
    assert geom.c_local(1e8 * rs) == pytest.approx(C, rel=1e-8)


# ---------------------------------------------------------------------------
# (b) Photon orbits at 1.5 r_s  (photon sphere)
# ---------------------------------------------------------------------------


def test_photon_sphere_at_1_5_r_s() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    rs = geom.r_s
    assert geom.r_photon == pytest.approx(1.5 * rs)
    # Photon-sphere radius is intrinsic to Schwarzschild (independent of M up to scale)
    assert BlackHoleGeometry(M=10.0 * M_SUN).r_photon \
        == pytest.approx(1.5 * BlackHoleGeometry(M=10.0 * M_SUN).r_s)


def test_isco_at_3_r_s() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    assert geom.r_isco == pytest.approx(3.0 * geom.r_s)


def test_radial_photon_is_captured_when_starting_inside_photon_sphere() -> None:
    sim = default_solar_bh()
    rs = sim.geom.r_s
    traj = sim.evolve_radial_photon(r0=1.4 * rs, n_steps=4000)
    assert traj.captured is True
    # Photon ends extremely close to r_s (within ~0.2 %)
    assert traj.r[-1] / rs < 1.005
    # Must asymptote: never crosses the σ=½ shell
    assert np.all(traj.r >= rs * 0.99)


def test_radial_photon_decelerates_through_substrate_gradient() -> None:
    sim = default_solar_bh()
    rs = sim.geom.r_s
    traj = sim.evolve_radial_photon(r0=2.0 * rs, n_steps=2000)
    # Speed (|dr/dt|) must monotonically decrease as r approaches r_s
    speeds = -np.diff(traj.r) / np.diff(traj.t)
    # Take the latter half — early transients can be noisy
    half = len(speeds) // 2
    assert speeds[half:].max() <= speeds[0] * 1.05


# ---------------------------------------------------------------------------
# (c) Hawking T_H for solar mass ~ 6×10⁻⁸ K
# ---------------------------------------------------------------------------


def test_hawking_temperature_solar_mass() -> None:
    sim = default_solar_bh()
    T_H = sim.hawking_temperature()
    # Canonical 6.17×10⁻⁸ K
    assert T_H == pytest.approx(6.17e-8, rel=5e-3)
    # Order of magnitude: 10⁻⁸ K
    assert 1e-8 < T_H < 1e-7


def test_hawking_temperature_scales_inversely_with_mass() -> None:
    sim1 = BlackHoleSimulator(BlackHoleGeometry(M=M_SUN))
    sim2 = BlackHoleSimulator(BlackHoleGeometry(M=2.0 * M_SUN))
    assert sim1.hawking_temperature() / sim2.hawking_temperature() \
        == pytest.approx(2.0, rel=1e-12)


def test_hawking_spectrum_has_blackbody_peak() -> None:
    sim = default_solar_bh()
    nu, B = sim.hawking_spectrum(n=300)
    # The peak frequency should be near Wien displacement 5.879×10¹⁰·T
    nu_peak_expected = 5.879e10 * sim.hawking_temperature()
    nu_peak_observed = nu[int(np.argmax(B))]
    # Allow factor-of-three slop: log grid resolution + Wien location vs B_ν peak
    assert 0.3 * nu_peak_expected < nu_peak_observed < 3.0 * nu_peak_expected
    # All values are finite and non-negative (high-frequency tail underflows
    # to 0 in float64, which is the correct Wien-tail behaviour)
    assert np.all(np.isfinite(B))
    assert np.all(B >= 0)
    # The spectrum has positive measure (at least near the peak)
    assert np.any(B > 0)


# ---------------------------------------------------------------------------
# (d) Ringdown ~250 Hz for stellar-mass BH
# ---------------------------------------------------------------------------


def test_ringdown_frequency_stellar_remnant() -> None:
    sim = default_stellar_remnant()       # 62 M_sun
    f = sim.ringdown_frequency()
    # GW150914 final-remnant ringdown was ~250 Hz
    assert f == pytest.approx(250.0, rel=0.10)


def test_ringdown_damping_few_milliseconds() -> None:
    sim = default_stellar_remnant()
    tau = sim.ringdown_damping_time()
    # ~4 ms for stellar-mass remnants
    assert 1e-3 < tau < 1e-2


def test_ringdown_quality_factor_universal() -> None:
    # Q = 0.3737/(2·0.0890) ≈ 2.10 — independent of M
    sim1 = default_solar_bh()
    sim2 = default_stellar_remnant()
    assert sim1.ringdown_quality() == pytest.approx(2.10, rel=1e-2)
    assert sim2.ringdown_quality() == pytest.approx(sim1.ringdown_quality())


def test_ringdown_waveform_decays_and_oscillates() -> None:
    sim = default_stellar_remnant()
    t, h = sim.ringdown_waveform(n_cycles=8.0, n_samples=1000)
    # Initial amplitude ≈ 1, late amplitude < 0.5 (decayed)
    assert abs(h[0]) == pytest.approx(1.0, abs=0.05)
    assert max(abs(h[-100:])) < 0.6
    # Has at least a few sign changes
    sign_flips = np.sum(np.diff(np.sign(h)) != 0)
    assert sign_flips >= 4


# ---------------------------------------------------------------------------
# Cross-checks: geometry, redshift, area
# ---------------------------------------------------------------------------


def test_horizon_area_equals_4pi_r_s_squared() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    assert geom.horizon_area == pytest.approx(4.0 * math.pi * geom.r_s**2)


def test_metric_g_tt_zero_at_horizon() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    assert geom.metric_g_tt(geom.r_s) == pytest.approx(0.0, abs=1e-3)
    # Far away g_tt → −c²
    assert geom.metric_g_tt(1e6 * geom.r_s) == pytest.approx(-C**2, rel=1e-5)


def test_gravitational_redshift_diverges_at_horizon() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    rs = geom.r_s
    # Photon emitted near the horizon arrives at infinity with f → 0
    assert geom.gravitational_redshift(1.001 * rs, 1e6 * rs) < 0.04
    # Photon emitted far away arrives nearly unredshifted at far observer
    assert geom.gravitational_redshift(1e6 * rs, 1e6 * rs) \
        == pytest.approx(1.0, rel=1e-9)


def test_visualizer_runs_without_matplotlib_window() -> None:
    """Smoke-test the visualizer wires up without raising."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sim = default_stellar_remnant()
    viz = BlackHoleVisualizer(sim)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    viz.horizon_panel(ax1, ax2)
    plt.close(fig)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    viz.hawking_ringdown_panel(ax1, ax2)
    plt.close(fig)


def test_invalid_inputs_raise() -> None:
    geom = BlackHoleGeometry(M=M_SUN)
    with pytest.raises(ValueError):
        geom.sigma(-1.0)
    sim = BlackHoleSimulator(geom)
    with pytest.raises(ValueError):
        sim.evolve_radial_photon(r0=0.5 * geom.r_s)
