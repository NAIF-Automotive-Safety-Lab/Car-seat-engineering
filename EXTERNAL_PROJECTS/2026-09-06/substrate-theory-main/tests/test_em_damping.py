"""Tests for EM radiation reaction."""

import numpy as np

from stiff_medium.atomic import (
    em_radiation_reaction,
    n_body_step_with_em_damping,
)


def test_em_force_zero_on_nucleus():
    positions = [
        np.array([0.0, 0.0, 0.0]),  # nucleus
        np.array([1.0, 0.0, 0.0]),  # electron
    ]
    velocities = [np.zeros(3), np.array([0.0, 1.0, 0.0])]
    charges = [+1.0, -1.0]
    bohr_radii = [0.0, 1.0]
    f = em_radiation_reaction(positions, velocities, 0, charges, bohr_radii, radiation_strength=1.0)
    assert np.allclose(f[0], 0.0), "Nucleus should feel no EM reaction force"


def test_em_force_zero_at_bohr_radius_with_no_radial_velocity():
    """Electron exactly at Bohr radius with purely tangential velocity:
    no radial drift → no damping force."""
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
    ]
    velocities = [np.zeros(3), np.array([0.0, 1.0, 0.0])]  # tangential only
    f = em_radiation_reaction(positions, velocities, 0, [+1.0, -1.0], [0.0, 1.0], 1.0)
    assert np.allclose(f[1], 0.0), "Electron at Bohr radius with tangential velocity should feel no EM"


def test_em_force_damps_outward_drift():
    """Electron drifting outward (radial velocity > 0) past Bohr radius:
    damping force pulls it back inward."""
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.5, 0.0, 0.0]),  # past r_bohr = 1.0
    ]
    velocities = [np.zeros(3), np.array([0.5, 0.0, 0.0])]  # outward radial velocity
    f = em_radiation_reaction(positions, velocities, 0, [+1.0, -1.0], [0.0, 1.0], 1.0)
    # Force on electron should be in -x direction (back toward nucleus)
    assert f[1][0] < 0, f"EM force should pull outward-drifting electron back inward; got {f[1]}"


def test_em_force_damps_inward_drift():
    """Electron drifting inward past Bohr radius: pushed back out."""
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.0, 0.0]),  # inside r_bohr = 1.0
    ]
    velocities = [np.zeros(3), np.array([-0.5, 0.0, 0.0])]  # inward radial velocity
    f = em_radiation_reaction(positions, velocities, 0, [+1.0, -1.0], [0.0, 1.0], 1.0)
    # Force should be in +x (push back outward)
    assert f[1][0] > 0, f"EM force should push inward-drifting electron outward; got {f[1]}"


def test_em_step_preserves_basic_structure():
    """Run with EM damping for 1000 steps; verify electrons stay bound."""
    a_n1 = 0.5
    pos = [
        np.array([0.0, 0.0, 0.0]),
        np.array([+a_n1, 0.0, 0.0]),
        np.array([-a_n1, 0.0, 0.0]),
    ]
    vel = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, +2.0, 0.0]),
        np.array([0.0, -2.0, 0.0]),
    ]
    masses = [1836.0, 1.0, 1.0]
    charges = [+2.0, -1.0, -1.0]
    spins = [0, 0, 1]
    bohr_radii = [0.0, a_n1, a_n1]

    for _ in range(1000):
        pos, vel = n_body_step_with_em_damping(
            pos, vel, masses, charges, spins, bohr_radii, nucleus_idx=0,
            dt=0.001,
            pauli_strength=0.5, pauli_radius=0.05,
            radiation_strength=2.0,
        )

    # Both electrons should still be bound
    d_a = float(np.linalg.norm(pos[1] - pos[0]))
    d_b = float(np.linalg.norm(pos[2] - pos[0]))
    assert d_a < 5 * a_n1, f"Electron A escaped: distance {d_a}"
    assert d_b < 5 * a_n1, f"Electron B escaped: distance {d_b}"
