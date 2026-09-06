"""Tests for k4_face_pair_geometry — explicit K_4 geometry + face-pair coupling."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.k4_face_pair_geometry import (
    DIHEDRAL_ANGLE_DEG,
    EPS_FACE_MEV,
    TETRAHEDRAL_ANGLE_DEG,
    XI_FM,
    K4FacePairCoupling,
    K4Geometry,
    rotation_aligning,
    rotation_matrix,
)


# ---------------------------------------------------------------------------
# K_4 geometry: combinatorics + invariants
# ---------------------------------------------------------------------------

def test_k4_combinatorics():
    """K_4 has 4 vertices, 6 edges, 4 faces (the complete graph on 4 nodes)."""
    g = K4Geometry()
    assert g.n_vertices == 4
    assert g.n_edges == 6
    assert g.n_faces == 4
    assert g.vertices.shape == (4, 3)
    assert g.edges.shape == (6, 2)
    assert g.faces.shape == (4, 3)


def test_k4_normals_unit():
    """All four face-normals have unit length."""
    g = K4Geometry()
    norms = np.linalg.norm(g.face_normals, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-12)


def test_k4_normals_sum_to_zero():
    """Outward face-normals of a regular tetrahedron sum to zero (4-fold symmetry)."""
    g = K4Geometry()
    s = g.normals_sum()
    assert np.allclose(s, 0.0, atol=1e-12), f"normals sum = {s}"


def test_tetrahedral_angle():
    """Vertex-centre-vertex angle is arccos(-1/3) ≈ 109.47°."""
    g = K4Geometry()
    assert abs(g.vertex_angle() - TETRAHEDRAL_ANGLE_DEG) < 1e-6


def test_normal_pair_angle():
    """Angle between two outward face-normals is arccos(-1/3) ≈ 109.47°.

    Outward face-normals point through opposite-vertex centroids, so
    they are antiparallel to scaled vertex vectors and their angle
    matches the vertex-vertex angle by symmetry.
    """
    g = K4Geometry()
    assert abs(g.normal_pair_angle() - TETRAHEDRAL_ANGLE_DEG) < 1e-6


def test_dihedral_angle():
    """Angle between two face PLANES is arccos(+1/3) ≈ 70.53°."""
    g = K4Geometry()
    assert abs(g.dihedral_angle() - DIHEDRAL_ANGLE_DEG) < 1e-6


def test_face_centroid_distance_matches_d_face():
    """Face centroid sits at distance d_face from the cell centre."""
    g = K4Geometry()
    distances = np.linalg.norm(g.face_centroids - g.centre, axis=1)
    assert np.allclose(distances, g.d_face, atol=1e-12)


def test_edges_connect_distinct_vertices():
    """All 6 edges are distinct unordered pairs of distinct vertices."""
    g = K4Geometry()
    pairs = {tuple(sorted(e)) for e in g.edges.tolist()}
    assert len(pairs) == 6
    for e in g.edges:
        assert e[0] != e[1]


def test_faces_are_triangles_opposite_vertices():
    """Each face is the set of 3 vertices opposite one excluded vertex."""
    g = K4Geometry()
    all_idx = set(range(4))
    for k, face in enumerate(g.faces):
        # face k should be opposite vertex k by our construction
        face_set = set(face.tolist())
        opposite = all_idx - face_set
        assert opposite == {k}


# ---------------------------------------------------------------------------
# Rotation helpers
# ---------------------------------------------------------------------------

def test_rotation_aligning_simple():
    """rotation_aligning(x̂, ŷ) takes (1,0,0) to (0,1,0)."""
    R = rotation_aligning(np.array([1.0, 0.0, 0.0]),
                          np.array([0.0, 1.0, 0.0]))
    out = R @ np.array([1.0, 0.0, 0.0])
    assert np.allclose(out, np.array([0.0, 1.0, 0.0]), atol=1e-12)


def test_rotation_aligning_identity():
    """Aligning a vector to itself is the identity."""
    v = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    R = rotation_aligning(v, v)
    assert np.allclose(R, np.eye(3), atol=1e-12)


def test_rotation_aligning_antiparallel():
    """Aligning v to -v is a 180° rotation."""
    v = np.array([0.0, 0.0, 1.0])
    R = rotation_aligning(v, -v)
    out = R @ v
    assert np.allclose(out, -v, atol=1e-12)


def test_rotation_matrix_orthogonal():
    """Rodrigues rotation is orthogonal with det = +1."""
    R = rotation_matrix(np.array([1.0, 2.0, 3.0]), 0.7)
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert abs(np.linalg.det(R) - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# Face-pair coupling: placement + matching
# ---------------------------------------------------------------------------

def test_place_two_cells_separation():
    """place_two_cells positions cells at ±d/2 along x."""
    coupling = K4FacePairCoupling()
    coupling.place_two_cells(d_centers=2.0)
    assert np.allclose(coupling.pos_a, np.array([-1.0, 0.0, 0.0]))
    assert np.allclose(coupling.pos_b, np.array([+1.0, 0.0, 0.0]))


def test_face_match_symmetric_placement():
    """For symmetric +x/-x placement the matched pair has anti-aligned normals.

    At d_centres = 1.4 fm with d_face = 0.8 fm the matched centroids are
    |1.4 - 2 × 0.8| = 0.2 fm apart; at d_centres = 1.6 fm they coincide.
    """
    coupling = K4FacePairCoupling()
    coupling.place_two_cells(d_centers=1.4)
    match = coupling.find_face_pair_match()
    # Normals should be very close to anti-aligned (dot ≈ -1)
    assert match.normal_dot < -0.99, f"matched dot = {match.normal_dot}"
    # Geometric prediction: |d_centres - 2 d_face| = |1.4 - 1.6| = 0.2 fm
    expected_d = abs(1.4 - 2 * coupling.cell_a.d_face)
    assert abs(match.centroid_distance - expected_d) < 0.02, (
        f"centroid dist = {match.centroid_distance}, expected {expected_d}")
    # And at d_centres = 1.6 fm the centroids coincide (within FP error)
    coupling.place_two_cells(d_centers=1.6)
    match2 = coupling.find_face_pair_match()
    assert match2.centroid_distance < 1e-6


def test_face_match_index_consistency():
    """By construction, matched face indices are (0, 0)."""
    coupling = K4FacePairCoupling()
    coupling.place_two_cells(d_centers=1.4)
    match = coupling.find_face_pair_match()
    assert match.face_a == 0
    assert match.face_b == 0


# ---------------------------------------------------------------------------
# Coupling energy: scan + minimum
# ---------------------------------------------------------------------------

def test_far_apart_zero_energy():
    """Cells at large separation have negligible coupling."""
    coupling = K4FacePairCoupling()
    E = coupling.coupling_energy(d=20.0, theta=0.0)
    assert abs(E) < 1e-6


def test_binding_energy_curve_has_minimum_near_xi():
    """U(d) curve has its minimum near d = ξ ≈ 1.4 fm."""
    coupling = K4FacePairCoupling()
    d_arr, E_arr = coupling.binding_energy_curve(
        d_range=np.linspace(0.6, 4.0, 80))
    d_min = d_arr[np.argmin(E_arr)]
    # d_min is the cell-centre separation; at the geometric optimum the
    # face centroids coincide → d_min = 2 * d_face = 1.6 fm.  We allow
    # some slack for the smooth gate.
    assert 1.0 < d_min < 2.0, f"d_min = {d_min} not in (1.0, 2.0) fm"


def test_minimum_energy_matches_eps_face():
    """At the relaxed minimum, BE = ε_face × 1 face = 2.222 MeV."""
    coupling = K4FacePairCoupling()
    d_min, theta_min, E_min = coupling.relax_to_minimum()
    BE = -E_min
    # We expect exactly one matched-face contribution ≈ ε_face;
    # the other 15 face-pair combinations are gated to zero by distance.
    assert abs(BE - EPS_FACE_MEV) < 0.05, f"BE = {BE} != {EPS_FACE_MEV}"


def test_relax_drives_to_anti_alignment():
    """After relaxation, the matched face-pair has dot ≈ -1."""
    coupling = K4FacePairCoupling()
    coupling.relax_to_minimum()
    match = coupling.find_face_pair_match()
    assert match.normal_dot < -0.99


def test_verify_deuteron_geometry():
    """verify_deuteron_geometry returns BE ≈ 2.222 MeV at d ≈ 1.4–1.6 fm."""
    coupling = K4FacePairCoupling()
    report = coupling.verify_deuteron_geometry(tol_be_mev=0.1, tol_d_fm=0.4)
    assert report['be_within_tol'], (
        f"BE = {report['binding_energy_MeV']} not near {EPS_FACE_MEV}")
    assert 1.0 < report['d_min_fm'] < 2.0
    assert report['matched_face_a'] == 0
    assert report['matched_face_b'] == 0


# ---------------------------------------------------------------------------
# Visualization smoke test (no display, tiny figure saved)
# ---------------------------------------------------------------------------

def test_make_visualization_smoke(tmp_path):
    """make_visualization runs without error and writes a PNG."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use('Agg')          # headless backend

    coupling = K4FacePairCoupling()
    coupling.place_two_cells(d_centers=1.4)
    out = tmp_path / "k4_pair.png"
    fig, ax = coupling.make_visualization(save_path=str(out))
    assert out.exists()
    assert out.stat().st_size > 0
    import matplotlib.pyplot as plt
    plt.close(fig)
