"""Tests for the 1D propagating EM field."""

import numpy as np

from stiff_medium.em_field import EMField1D


def test_field_starts_at_zero():
    field = EMField1D(x_min=-1.0, x_max=1.0, n_points=21, c=1.0)
    assert np.allclose(field.phi_curr, 0.0)
    assert np.allclose(field.phi_prev, 0.0)


def test_no_source_no_propagation():
    """With zero initial conditions and no source, field stays zero."""
    field = EMField1D(x_min=-1.0, x_max=1.0, n_points=21, c=1.0)
    for _ in range(50):
        field.step(dt=0.05)
    assert np.allclose(field.phi_curr, 0.0)


def test_source_creates_disturbance():
    """A point source at the origin should disturb the field."""
    field = EMField1D(x_min=-2.0, x_max=2.0, n_points=41, c=1.0)
    for _ in range(20):
        field.step(dt=0.05, sources={0.0: 1.0})
    # Field at origin should be non-zero
    assert abs(field.value_at(0.0)) > 1e-3


def test_disturbance_propagates_at_c():
    """A pulse at origin should reach a distant point at time = distance/c."""
    field = EMField1D(x_min=-1.0, x_max=10.0, n_points=221, c=1.0)
    # Single brief pulse at t=0
    field.step(dt=0.05, sources={0.0: 100.0})
    for _ in range(30):
        field.step(dt=0.05)
    # At t = 30*0.05 = 1.5, wave should have traveled ~1.5 units
    # Check field at x=1.5 vs x=5.0
    near = abs(field.value_at(1.5))
    far = abs(field.value_at(5.0))
    assert near > far, f"Wave should be near 1.5 (got {near}) more than near 5.0 (got {far})"


def test_cfl_violation_raises():
    """If c·dt/dx > 1, simulation should refuse (Courant-Friedrichs-Lewy)."""
    import pytest
    field = EMField1D(x_min=0.0, x_max=1.0, n_points=11, c=1.0)
    # dx = 0.1; c·dt/dx > 1 means dt > 0.1
    with pytest.raises(ValueError, match="CFL"):
        field.step(dt=0.5)


def test_resonant_absorber_gains_more_energy_than_non_resonant():
    """Source emitting at ω=2, two absorbers (one at ω=2, one at ω=5).
    The matched absorber should gain more energy."""
    DT = 0.025
    N_STEPS = 800
    omega_source = 2.0

    # Resonant case
    field = EMField1D(x_min=-5.0, x_max=15.0, n_points=401, c=1.0)
    abs_x = 8.0
    pos_r, vel_r, max_r = 0.0, 0.0, 0.0
    omega_r = 2.0
    for k in range(N_STEPS):
        t = k * DT
        src = 5.0 * omega_source**2 * np.sin(omega_source * t)
        field.step(DT, {0.0: src})
        force = field.gradient_at(abs_x) - omega_r**2 * pos_r - 0.05 * vel_r
        vel_r += force * DT
        pos_r += vel_r * DT
        max_r = max(max_r, abs(pos_r))

    # Non-resonant case
    field2 = EMField1D(x_min=-5.0, x_max=15.0, n_points=401, c=1.0)
    pos_n, vel_n, max_n = 0.0, 0.0, 0.0
    omega_n = 5.0
    for k in range(N_STEPS):
        t = k * DT
        src = 5.0 * omega_source**2 * np.sin(omega_source * t)
        field2.step(DT, {0.0: src})
        force = field2.gradient_at(abs_x) - omega_n**2 * pos_n - 0.05 * vel_n
        vel_n += force * DT
        pos_n += vel_n * DT
        max_n = max(max_n, abs(pos_n))

    assert max_r > 3 * max_n, (
        f"Resonant absorber should gain more energy ({max_r}) "
        f"than non-resonant ({max_n})"
    )
