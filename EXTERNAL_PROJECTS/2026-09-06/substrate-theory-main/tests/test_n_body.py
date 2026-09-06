"""Tests for N-body atomic dynamics."""

import numpy as np

from stiff_medium.atomic import n_body_force, n_body_newton_step


def test_n_body_force_zero_at_infinite_separation():
    """Two particles infinitely far apart feel no Coulomb force."""
    positions = [np.array([0.0, 0.0, 0.0]), np.array([1e10, 0.0, 0.0])]
    charges = [1.0, -1.0]
    forces = n_body_force(positions, charges)
    assert np.allclose(forces[0], 0.0, atol=1e-15)
    assert np.allclose(forces[1], 0.0, atol=1e-15)


def test_n_body_force_attractive_between_opposite_charges():
    """Particle with charge +1 should be pulled toward particle with charge −1."""
    positions = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
    charges = [1.0, -1.0]
    forces = n_body_force(positions, charges)
    # Force on particle 0 should point in +x (toward particle 1).
    assert forces[0][0] > 0
    # Force on particle 1 should point in −x (toward particle 0).
    assert forces[1][0] < 0
    # Magnitudes should be equal (Newton's third law).
    assert np.isclose(np.linalg.norm(forces[0]), np.linalg.norm(forces[1]))


def test_n_body_force_repulsive_between_same_charges():
    """Two same-charge particles should repel."""
    positions = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
    charges = [-1.0, -1.0]
    forces = n_body_force(positions, charges)
    # Force on particle 0 should point in −x (away from particle 1).
    assert forces[0][0] < 0
    # Force on particle 1 should point in +x (away from particle 0).
    assert forces[1][0] > 0


def test_n_body_force_three_particles():
    """Test that 3-body forces sum correctly."""
    # Symmetric: particle 0 at origin, particles 1 and 2 at +x and −x.
    # If charges are +1, −1, −1: forces on particle 0 should cancel by symmetry.
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([-1.0, 0.0, 0.0]),
    ]
    charges = [1.0, -1.0, -1.0]
    forces = n_body_force(positions, charges)
    # Force on particle 0 should be near zero (symmetric attractions cancel).
    assert np.allclose(forces[0], 0.0, atol=1e-12)


def test_n_body_step_preserves_total_momentum():
    """In a closed N-body system, total linear momentum is conserved."""
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.5, 0.0]),
        np.array([-0.5, -0.5, 0.0]),
    ]
    velocities = [
        np.array([0.1, 0.0, 0.0]),
        np.array([0.0, 0.2, 0.0]),
        np.array([-0.05, -0.1, 0.0]),
    ]
    masses = [1.0, 2.0, 0.5]
    charges = [1.0, -1.0, -1.0]

    initial_p = sum(m * v for m, v in zip(masses, velocities))

    for _ in range(100):
        positions, velocities = n_body_newton_step(
            positions, velocities, masses, charges, dt=0.001
        )

    final_p = sum(m * v for m, v in zip(masses, velocities))
    assert np.allclose(initial_p, final_p, atol=1e-9)


def test_helium_binding_2_electrons_stay_bound():
    """Two electrons binding to a Z=2 nucleus should stay bound for many steps."""
    M_E = 1.0
    M_NUCLEUS = 7294.3
    Z = 2

    # Set up at Bohr radius for Z=2: a_bohr = 1/Z = 0.5
    r_bohr = 1.0 / Z
    v_e = float(np.sqrt(Z / r_bohr))  # circular orbit velocity

    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([-r_bohr, 0.0, 0.0]),
        np.array([r_bohr, 0.0, 0.0]),
    ]
    velocities = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, v_e, 0.0]),
        np.array([0.0, -v_e, 0.0]),
    ]
    masses = [M_NUCLEUS, M_E, M_E]
    charges = [+2.0, -1.0, -1.0]

    for _ in range(10000):
        positions, velocities = n_body_newton_step(
            positions, velocities, masses, charges, dt=0.0001
        )

    d_a = float(np.linalg.norm(positions[1] - positions[0]))
    d_b = float(np.linalg.norm(positions[2] - positions[0]))
    # Bound = within 10 × Bohr radius
    assert d_a < 10 * r_bohr, f"Electron A escaped: distance {d_a} > {10 * r_bohr}"
    assert d_b < 10 * r_bohr, f"Electron B escaped: distance {d_b} > {10 * r_bohr}"
