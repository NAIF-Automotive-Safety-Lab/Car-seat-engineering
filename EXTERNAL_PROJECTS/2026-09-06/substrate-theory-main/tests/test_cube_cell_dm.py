"""Tests for cube_cell_dm_simulator."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.cube_cell_dm_simulator import (
    CubeCell,
    SubstrateParams,
    adjacency_matrix,
    default_substrate_params,
    em_cross_section,
    is_parity_conserving,
    leading_multipole_order,
    m_dm_from_first_excited,
    mode_breathing,
    mode_energy,
    mode_first_excited,
    mode_ground,
    octupole,
    parity_operator,
    parity_preserved_under_evolution,
    two_cell_gravity,
)


def test_cube_geometry_8_vertices_12_edges():
    cell = CubeCell()
    assert cell.vertices().shape == (8, 3)
    assert cell.edges().shape == (12, 2)


def test_no_triangular_faces():
    cell = CubeCell()
    assert cell.n_triangular_faces() == 0
    # all 6 faces are quadrilaterals
    assert all(len(f) == 4 for f in cell.faces())


def test_bipartite_parity_classes_4_and_4():
    parity = CubeCell.parity()
    assert int(np.sum(parity == 0)) == 4
    assert int(np.sum(parity == 1)) == 4


def test_ground_mode_zero_energy():
    cell = CubeCell()
    params = default_substrate_params()
    disp = mode_ground(cell, amp=0.1)
    E = mode_energy(cell, disp, params)
    assert E == pytest.approx(0.0, abs=1e-10)


def test_first_excited_positive_energy():
    cell = CubeCell()
    params = default_substrate_params()
    E = mode_energy(cell, mode_first_excited(cell, amp=1.0), params)
    assert E > 0


def test_breathing_positive_energy():
    cell = CubeCell()
    params = default_substrate_params()
    E = mode_energy(cell, mode_breathing(cell, amp=0.1), params)
    assert E > 0


def test_adjacency_commutes_with_parity():
    """Cube adjacency is bipartite — anti-commutes with P, but P^2=I, so
    H = A^2 commutes with P (both are parity-conserving observables)."""
    A = adjacency_matrix()
    H = A @ A
    assert is_parity_conserving(H)


def test_parity_conserved_under_unitary_evolution():
    """Even-power Hamiltonians preserve parity expectation."""
    A = adjacency_matrix()
    H = A @ A  # bipartite-block-diagonal
    P = parity_operator()
    # build a parity-eigenstate (eigenvector of P with eigenvalue +1)
    w, V = np.linalg.eigh(P)
    psi0 = V[:, np.argmax(w)].astype(np.complex128)
    t_values = np.linspace(0.0, 5.0, 11)
    assert parity_preserved_under_evolution(H, t_values, psi0)


def test_octupole_is_leading_multipole():
    cell = CubeCell(scale=1.0)
    assert leading_multipole_order(cell) == 3
    assert abs(octupole(cell)) > 1e-6


def test_m_dm_in_25_to_30_gev():
    params = default_substrate_params()
    m = m_dm_from_first_excited(params)
    assert 25.0 <= m <= 30.0


def test_em_cross_section_strongly_suppressed():
    cell = CubeCell(scale=1.0)
    sigma = em_cross_section(cell)
    # Long-wavelength octupole suppression at galactic v:
    # (m v / ℏ · a)^6 with a ~ 1 m is astronomically tiny — much smaller
    # than 10^-11 m^2 (the loose upper bound demanded by the spec).
    assert sigma < 1e-11


def test_two_cell_gravity_attractive_inverse_square():
    a = CubeCell(center=np.array([0.0, 0.0, 0.0]))
    b1 = CubeCell(center=np.array([1.0, 0.0, 0.0]))
    b2 = CubeCell(center=np.array([2.0, 0.0, 0.0]))
    V1 = two_cell_gravity(a, b1, m_dm_gev=27.5)
    V2 = two_cell_gravity(a, b2, m_dm_gev=27.5)
    assert V1 < 0 and V2 < 0
    # V ~ 1/r → V1 / V2 should equal 2
    assert V1 / V2 == pytest.approx(2.0, rel=1e-12)
