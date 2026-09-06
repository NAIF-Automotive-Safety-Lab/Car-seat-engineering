"""Tests for Pauli-aware N-body force."""

import numpy as np

from stiff_medium.atomic import (
    n_body_force_with_pauli,
    n_body_step_with_pauli,
)


def test_pauli_force_zero_for_different_charges():
    """Pauli should only fire for SAME charge AND same spin. Different
    charges → no Pauli term, only Coulomb."""
    positions = [np.array([0.0, 0.0, 0.0]), np.array([0.05, 0.0, 0.0])]
    charges = [1.0, -1.0]  # different
    spins = [0, 0]
    f = n_body_force_with_pauli(positions, charges, spins, pauli_radius=0.1)
    f_no_pauli = n_body_force_with_pauli(positions, charges, [0, 1], pauli_radius=0.1)
    # Same magnitude either way (Pauli didn't fire because charges differ)
    assert np.isclose(np.linalg.norm(f[0]), np.linalg.norm(f_no_pauli[0]))


def test_pauli_force_zero_for_different_spins():
    """Same charge but different spins → no Pauli term."""
    positions = [np.array([0.0, 0.0, 0.0]), np.array([0.05, 0.0, 0.0])]
    charges = [-1.0, -1.0]  # same
    f_same_spin = n_body_force_with_pauli(positions, charges, [0, 0], pauli_radius=0.1)
    f_diff_spin = n_body_force_with_pauli(positions, charges, [0, 1], pauli_radius=0.1)
    # Same-spin should have stronger repulsion (Coulomb + Pauli)
    assert np.linalg.norm(f_same_spin[0]) > np.linalg.norm(f_diff_spin[0])


def test_pauli_force_repulsive():
    """Same-charge same-spin pair: Pauli adds to Coulomb repulsion (both push apart)."""
    positions = [np.array([0.0, 0.0, 0.0]), np.array([0.05, 0.0, 0.0])]
    charges = [-1.0, -1.0]
    spins = [0, 0]
    f = n_body_force_with_pauli(positions, charges, spins, pauli_radius=0.1, pauli_strength=1.0)
    # Force on particle 0 should point in −x (away from particle 1)
    assert f[0][0] < 0
    assert f[1][0] > 0


def test_pauli_force_falls_off_quickly():
    """At d ≫ pauli_radius, the Pauli term should be small compared to Coulomb."""
    # Close range
    positions_close = [np.array([0.0, 0.0, 0.0]), np.array([0.1, 0.0, 0.0])]
    f_close = n_body_force_with_pauli(positions_close, [-1.0, -1.0], [0, 0],
                                       pauli_strength=1.0, pauli_radius=0.1)
    f_close_no_pauli = n_body_force_with_pauli(positions_close, [-1.0, -1.0], [0, 1],
                                                pauli_strength=1.0, pauli_radius=0.1)
    pauli_contrib_close = np.linalg.norm(f_close[0]) - np.linalg.norm(f_close_no_pauli[0])

    # Far range
    positions_far = [np.array([0.0, 0.0, 0.0]), np.array([10.0, 0.0, 0.0])]
    f_far = n_body_force_with_pauli(positions_far, [-1.0, -1.0], [0, 0],
                                     pauli_strength=1.0, pauli_radius=0.1)
    f_far_no_pauli = n_body_force_with_pauli(positions_far, [-1.0, -1.0], [0, 1],
                                              pauli_strength=1.0, pauli_radius=0.1)
    pauli_contrib_far = np.linalg.norm(f_far[0]) - np.linalg.norm(f_far_no_pauli[0])

    # Pauli at far range should be much smaller than at close range
    assert pauli_contrib_far < pauli_contrib_close * 1e-6


def test_pauli_step_preserves_momentum_with_three_particles():
    """Total momentum still conserved with Pauli."""
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
    ]
    velocities = [np.zeros(3), np.zeros(3), np.zeros(3)]
    masses = [1.0, 1.0, 1.0]
    charges = [-1.0, -1.0, -1.0]
    spins = [0, 0, 1]

    initial_p = sum(m * v for m, v in zip(masses, velocities))
    for _ in range(100):
        positions, velocities = n_body_step_with_pauli(
            positions, velocities, masses, charges, spins, dt=0.001
        )
    final_p = sum(m * v for m, v in zip(masses, velocities))
    assert np.allclose(initial_p, final_p, atol=1e-9)
