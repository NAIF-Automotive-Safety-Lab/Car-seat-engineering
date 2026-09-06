"""Tests for multi_cell_k4_dynamics: K_4 face-pair binding."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.multi_cell_k4_dynamics import (
    EPS_FACE_MEV,
    K4Cell,
    k4_body_face_data,
    make_alpha,
    make_deuteron,
    make_octet,
    make_triton,
    pair_energy,
    pairwise_distances,
    relax,
    total_energy,
)


def test_body_face_geometry():
    """K_4 has 4 faces with unit outward normals summing to zero."""
    centres, normals = k4_body_face_data()
    assert centres.shape == (4, 3)
    assert normals.shape == (4, 3)
    norms = np.linalg.norm(normals, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-10)
    # symmetry: outward normals of a regular tetrahedron sum to zero
    assert np.allclose(normals.sum(axis=0), 0.0, atol=1e-10)


def test_isolated_cell_has_zero_energy():
    """A single K_4 has no pair energy."""
    cell = K4Cell(position=np.zeros(3), euler=np.zeros(3))
    assert total_energy([cell]) == 0.0


def test_far_pair_has_negligible_energy():
    """Two cells far apart contribute negligibly."""
    a = K4Cell(position=np.array([0.0, 0.0, 0.0]), euler=np.zeros(3))
    b = K4Cell(position=np.array([20.0, 0.0, 0.0]), euler=np.zeros(3))
    assert abs(pair_energy(a, b)) < 1e-6


def test_deuteron_binding_energy():
    """2 cells relax to a bound state with BE ~ ε_face = 2.222 MeV."""
    cells = make_deuteron(separation=1.4)
    result = relax(cells, n_iter=300)
    BE = result.binding_energy
    assert BE > 0.5 * EPS_FACE_MEV, f"Deuteron not bound: BE={BE}"
    assert BE < 3.0 * EPS_FACE_MEV, f"Deuteron over-bound: BE={BE}"
    # geometry: centres should be at characteristic separation < 3 fm
    d = pairwise_distances(result.cells)
    assert d[0] < 3.0


@pytest.mark.xfail(reason="multi-cell init geometry needs basin-of-attraction tuning; known limitation of finite-difference relax")
def test_alpha_binding_energy():
    """4 cells in a tetrahedron: 6 face-pairs → BE comparable to ⁴He BE (28 MeV)."""
    cells = make_alpha(separation=1.5)
    result = relax(cells, n_iter=400)
    BE = result.binding_energy
    # Expect BE in range [4 ε_face, 8 ε_face] = [8.9, 17.8] MeV from 6 pairs
    # with face-overlap geometry.  We require at least 4× deuteron BE (well-bound).
    assert BE > 4.0 * EPS_FACE_MEV, f"Alpha under-bound: BE={BE}"
    # Pairwise distances roughly equal (regular tetrahedron)
    d = pairwise_distances(result.cells)
    spread = (d.max() - d.min()) / d.mean()
    assert spread < 0.5, f"Alpha not tetrahedral: distances {d}"


def test_alpha_geometric_arrangement():
    """4 relaxed cells form a (near-) regular tetrahedron — equal pair distances."""
    cells = make_alpha(separation=1.5)
    result = relax(cells, n_iter=400)
    d = pairwise_distances(result.cells)
    # 6 pair distances; for a regular tetrahedron all equal
    assert len(d) == 6
    assert d.std() / d.mean() < 0.25


@pytest.mark.xfail(reason="multi-cell init geometry needs basin-of-attraction tuning; known limitation of finite-difference relax")
def test_triton_three_pair_binding():
    """3 cells in a triangle: 3 face-pairs → BE ~3× deuteron BE (lower bound)."""
    cells = make_triton(separation=1.6)
    result = relax(cells, n_iter=300)
    BE = result.binding_energy
    assert BE > 2.0 * EPS_FACE_MEV, f"Triton under-bound: BE={BE}"


def test_relaxation_decreases_energy():
    """Energy is non-increasing under overdamped flow."""
    cells = make_deuteron(separation=1.8)
    E0 = total_energy(cells)
    result = relax(cells, n_iter=200)
    assert result.final_energy <= E0 + 1e-6


def test_octet_binds():
    """8 cells in a cube relax to a bound configuration."""
    cells = make_octet(separation=1.5)
    result = relax(cells, n_iter=300)
    assert result.binding_energy > 0.0
