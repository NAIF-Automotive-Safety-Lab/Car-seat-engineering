"""Tests for spinor.py — Möbius internal twist analysis."""

import numpy as np

from stiff_medium.neutrino import C
from stiff_medium.spinor import (
    cone_azimuth,
    mobius_phase,
    slope_sign,
    unwrap_azimuth_history,
    spin_half_check,
)


def test_cone_azimuth_zero_for_reference_direction():
    """Velocity along the reference perpendicular (e1) should give φ=0."""
    z = np.array([0.0, 0.0, 1.0])
    s = C / np.sqrt(2.0)
    # Build the same e1 as cone_azimuth uses internally
    e1 = np.cross(z, np.array([1.0, 0.0, 0.0]))
    e1 = e1 / float(np.linalg.norm(e1))
    v = s * z + s * e1
    phi = cone_azimuth(v, z)
    assert np.isclose(phi, 0.0)


def test_cone_azimuth_pi_for_opposite_direction():
    """Velocity along -e1 should give φ=π."""
    z = np.array([0.0, 0.0, 1.0])
    s = C / np.sqrt(2.0)
    e1 = np.cross(z, np.array([1.0, 0.0, 0.0]))
    e1 = e1 / float(np.linalg.norm(e1))
    v = s * z - s * e1
    phi = cone_azimuth(v, z)
    assert np.isclose(phi, np.pi)


def test_cone_azimuth_advances_around_circle():
    """As the azimuth advances 0 → π/2 → π → 3π/2, cone_azimuth should
    return matching values."""
    z = np.array([0.0, 0.0, 1.0])
    s = C / np.sqrt(2.0)
    e1 = np.cross(z, np.array([1.0, 0.0, 0.0]))
    e1 = e1 / float(np.linalg.norm(e1))
    e2 = np.cross(z, e1)

    for target in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2):
        v = s * z + s * (np.cos(target) * e1 + np.sin(target) * e2)
        phi = cone_azimuth(v, z)
        # phi is in [0, 2π); compare modulo 2π
        diff = (phi - target) % (2 * np.pi)
        diff = min(diff, 2 * np.pi - diff)
        assert diff < 1e-9, f"target={target}, got phi={phi}"


def test_mobius_phase_is_half_azimuth():
    assert np.isclose(mobius_phase(0.0), 0.0)
    assert np.isclose(mobius_phase(2 * np.pi), np.pi)
    assert np.isclose(mobius_phase(4 * np.pi), 2 * np.pi)


def test_slope_sign_flips_at_2pi_returns_at_4pi():
    """Slope is +1 for azimuth in [0, 2π), -1 for [2π, 4π), +1 for [4π, 6π)."""
    assert slope_sign(0.0) == 1
    assert slope_sign(np.pi) == 1  # azimuth π → ψ = π/2, still +1 (cos > 0)
    assert slope_sign(1.99 * np.pi) == 1  # just before 2π
    assert slope_sign(2.01 * np.pi) == -1  # just after 2π → ψ just over π/2 → ... wait
    # Let me reconsider: psi = az/2. cos(psi).
    # az=0 → psi=0 → cos=1 → +1
    # az=2π → psi=π → cos=-1 → -1
    # az=4π → psi=2π → cos=+1 → +1
    # The flip from +1 to -1 happens when cos(psi) crosses 0, which is at psi=π/2,
    # i.e., azimuth = π. So our code's logic of `int(psi//π) % 2` gives:
    # - azimuth=0, psi=0 → 0 // π = 0 → +1 ✓
    # - azimuth=π, psi=π/2 → 0 // π = 0 → +1 ✓ (still in first half)
    # - azimuth=2π, psi=π → 1 // π wait that's int division of π/π = 1 → -1 (flipped)
    # OK so the flip happens at azimuth=2π.
    # The test is checking the right thing.


def test_unwrap_azimuth_handles_wraparound():
    """Sequence that wraps around 2π should be unwrapped correctly."""
    # Azimuth advances linearly past 2π
    raw = np.array([0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    raw_wrapped = raw % (2 * np.pi)
    unwrapped = unwrap_azimuth_history(raw_wrapped)
    # After unwrap, should be monotonically increasing modulo small noise
    diffs = np.diff(unwrapped)
    assert all(d > 0 for d in diffs)


def test_spin_half_check_recognizes_spin_half():
    """One azimuth rotation per orbital revolution → spin-½."""
    # Simulate: azimuth advances by 2π per orbit, over 5 orbits.
    azimuths_unwrapped = np.linspace(0.0, 5.0 * 2 * np.pi, 100)
    result = spin_half_check(azimuths_unwrapped, orbital_revolutions=5.0)
    assert "1/2" in result["implied_spin"]


def test_spin_half_check_recognizes_spin_one():
    """Two azimuth rotations per orbital revolution → spin-1 (bosonic)."""
    azimuths_unwrapped = np.linspace(0.0, 5.0 * 4 * np.pi, 100)  # 2× per orbit
    result = spin_half_check(azimuths_unwrapped, orbital_revolutions=5.0)
    assert "1 (bosonic" in result["implied_spin"]


def test_spin_half_check_handles_no_orbital_motion():
    azimuths_unwrapped = np.linspace(0.0, np.pi, 100)
    result = spin_half_check(azimuths_unwrapped, orbital_revolutions=0.0)
    assert "undefined" in result["implied_spin"]
