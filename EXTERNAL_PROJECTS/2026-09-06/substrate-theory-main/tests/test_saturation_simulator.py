"""Tests for src.stiff_medium.saturation_simulator.

Verifies:
  - σ never exceeds 0.5 + ε in any of the four scenarios
  - c_local(σ) is strictly monotone-decreasing in σ
  - shock front relaxes to a finite thickness on the order of ξ
  - black-hole-like core has c_local → 0 at the centre
  - energy is conserved (or strictly decreasing under drag)
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.saturation_simulator import (
    SIGMA_MAX,
    apply_saturation_cap,
    c_local,
    cap_violation,
    energy_drift,
    measure_front_thickness,
    run_all_scenarios,
    scenario_black_hole,
    scenario_inverse_radius,
    scenario_linear_pulse,
    scenario_shock_front,
    speed_vs_sigma_curve,
    strain,
)


EPS = 1e-9


# ---------------------------------------------------------------------------
# Cap enforcement
# ---------------------------------------------------------------------------

def test_apply_saturation_cap_clips_to_half():
    u = np.array([-2.0, -0.6, -0.4, 0.0, 0.4, 0.6, 2.0])
    capped = apply_saturation_cap(u)
    assert np.all(np.abs(capped) <= SIGMA_MAX + EPS)
    # Below-cap values are untouched
    assert capped[2] == pytest.approx(-0.4)
    assert capped[4] == pytest.approx(0.4)


def test_strain_definition():
    u = np.array([-0.5, 0.0, 0.25, 0.5])
    assert np.allclose(strain(u), [0.5, 0.0, 0.25, 0.5])


# ---------------------------------------------------------------------------
# Wave-speed monotonicity
# ---------------------------------------------------------------------------

def test_c_local_monotone_decreasing():
    sigma, c = speed_vs_sigma_curve()
    diffs = np.diff(c)
    # Strictly non-increasing
    assert np.all(diffs <= 1e-12)
    # And actually decreases (not flat)
    assert c[0] > c[-1]


def test_c_local_vanishes_at_cap():
    assert c_local(np.array([SIGMA_MAX]))[0] == pytest.approx(0.0, abs=1e-12)
    assert c_local(np.array([0.0]))[0] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Scenario (a): linear pulse — cap respected
# ---------------------------------------------------------------------------

def test_linear_pulse_respects_cap():
    res = scenario_linear_pulse(t_max=4.0)
    assert cap_violation(res) <= SIGMA_MAX + 1e-6


# ---------------------------------------------------------------------------
# Scenario (b): inverse-radius singularity is regularised
# ---------------------------------------------------------------------------

def test_inverse_radius_capped():
    res = scenario_inverse_radius(t_max=1.0)
    sigma_max = cap_violation(res)
    assert sigma_max <= SIGMA_MAX + 1e-6
    # And the central peak is at the cap (the singular core saturates)
    assert sigma_max >= SIGMA_MAX - 0.05


# ---------------------------------------------------------------------------
# Scenario (c): shock front has finite thickness ~ ξ
# ---------------------------------------------------------------------------

def test_shock_front_finite_thickness():
    res = scenario_shock_front(t_max=3.0)
    assert cap_violation(res) <= SIGMA_MAX + 1e-6
    final_u = res.u_history[-1]
    width = measure_front_thickness(final_u, res.x)
    # Thickness should be O(ξ)=O(1).  Allow 0.05–10 ξ.
    assert 0.05 < width < 10.0


# ---------------------------------------------------------------------------
# Scenario (d): black-hole-like core has c_local → 0 at centre
# ---------------------------------------------------------------------------

def test_black_hole_horizon_like_core():
    res = scenario_black_hole(peak_sigma=0.499, t_max=0.05)
    assert cap_violation(res) <= SIGMA_MAX + 1e-6
    # At final time, find centre
    u_final = res.u_history[-1]
    sigma_final = np.abs(u_final) / 1.0
    c_field = c_local(sigma_final)
    # Edge has c ≈ 1; centre has c ≪ 1
    centre_idx = len(u_final) // 2
    assert c_field[centre_idx] < 0.4
    assert c_field[0] > 0.95
    assert c_field[-1] > 0.95


# ---------------------------------------------------------------------------
# Energy conservation / dissipation
# ---------------------------------------------------------------------------

def test_energy_conservation_no_drag():
    # The linear-pulse scenario uses γ = 0
    res = scenario_linear_pulse(t_max=4.0)
    de, rel = energy_drift(res)
    # Allow up to ~5% drift on this resolution / step
    assert abs(rel) < 0.05


def test_energy_decreases_with_drag():
    res = scenario_inverse_radius(t_max=1.0)
    e = res.energy
    # Energy should be (weakly) non-increasing once drag dominates
    assert e[-1] <= e[0] + 1e-6


# ---------------------------------------------------------------------------
# Bundle: full suite runs
# ---------------------------------------------------------------------------

def test_run_all_scenarios_completes():
    rep = run_all_scenarios()
    for name, v in rep["cap_violations"].items():
        assert v <= SIGMA_MAX + 1e-6, f"{name} violated cap with σ_max={v}"
    assert rep["shock_thickness"] > 0.0
