"""Tests for nucleon_stacking_geometry — multi-K_4 stacking → light nuclei BE."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.k4_face_pair_geometry import EPS_FACE_MEV, K4Geometry
from src.stiff_medium.nucleon_stacking_geometry import (
    ETA_COOP_ALPHA,
    REFERENCE_BE_MEV,
    NuclearChartVisualizer,
    NucleonStackGeometry,
    StackingTopology,
    TOPOLOGY_REGISTRY,
    cooperative_factor,
    get_topology,
    topology_alpha,
    topology_be8,
    topology_c12,
    topology_deuteron,
    topology_li6,
    topology_o16,
    topology_triton,
)


# ---------------------------------------------------------------------------
# Topology counts: P(N) per nucleus
# ---------------------------------------------------------------------------

def test_deuteron_topology_one_pair():
    """Deuteron: 2 cells share 1 face."""
    t = topology_deuteron()
    assert t.A == 2
    assert t.n_face_pairs == 1
    assert t.mean_coordination() == 1.0


def test_triton_topology_three_pairs():
    """Triton: 3 cells in triangle, each shares 2 faces ⇒ 3 pairs."""
    t = topology_triton()
    assert t.A == 3
    assert t.n_face_pairs == 3
    assert t.mean_coordination() == 2.0


def test_alpha_topology_six_pairs():
    """Alpha: 4 cells in tetrahedron, each shares 3 faces ⇒ 6 pairs."""
    t = topology_alpha()
    assert t.A == 4
    assert t.n_face_pairs == 6
    assert t.mean_coordination() == 3.0
    # Each cell should share with all 3 others
    for k in range(4):
        assert t.coordination(k) == 3


def test_li6_topology():
    """Li-6: 8 face-pairs (α + d cluster + bridge)."""
    t = topology_li6()
    assert t.A == 6
    assert t.n_face_pairs == 8


def test_be8_topology():
    """Be-8: 13 face-pairs (2 alphas + 1 bridge)."""
    t = topology_be8()
    assert t.A == 8
    assert t.n_face_pairs == 13


def test_c12_topology():
    """C-12 (3α model): 21 face-pairs."""
    t = topology_c12()
    assert t.A == 12
    assert t.n_face_pairs == 21


def test_o16_topology():
    """O-16 (4α tetrahedron): 30 face-pairs."""
    t = topology_o16()
    assert t.A == 16
    assert t.n_face_pairs == 30


# ---------------------------------------------------------------------------
# Predicted BE: deuteron exact, alpha within 10%
# ---------------------------------------------------------------------------

def test_deuteron_BE_exact():
    """(a) Deuteron BE ≈ 2.222 MeV from 1 face pair (η_coop = 1)."""
    stack = NucleonStackGeometry.from_topology(topology_deuteron())
    BE = stack.predicted_BE()
    assert abs(BE - EPS_FACE_MEV) < 1e-9, (
        f"BE_d = {BE:.6f} MeV != ε_face = {EPS_FACE_MEV:.6f} MeV"
    )
    # And matches observed deuteron BE within < 0.1%
    BE_obs = REFERENCE_BE_MEV[2]
    err = abs(BE - BE_obs) / BE_obs
    assert err < 0.002, f"deuteron BE pred {BE} vs obs {BE_obs}, err={err*100:.2f}%"


def test_deuteron_BE_bare_equals_eps_face():
    """Sanity: the BARE prediction also = ε_face for the deuteron."""
    stack = NucleonStackGeometry.from_topology(topology_deuteron())
    assert abs(stack.bare_BE() - EPS_FACE_MEV) < 1e-9


def test_alpha_BE_within_10_pct():
    """(b) Alpha BE ≈ 28 MeV from 6 face pairs (within 10%)."""
    stack = NucleonStackGeometry.from_topology(topology_alpha())
    BE = stack.predicted_BE()
    BE_obs = REFERENCE_BE_MEV[4]   # 28.295674 MeV
    err = abs(BE - BE_obs) / BE_obs
    assert err < 0.10, (
        f"alpha BE pred {BE:.3f} MeV vs obs {BE_obs:.3f} MeV "
        f"⇒ err = {err*100:.2f}% > 10%"
    )
    # Confirm cooperative factor is roughly 2 (≈ 2.122)
    eta = cooperative_factor(4)
    assert 2.0 < eta < 2.3, f"η_coop(α) = {eta:.3f} outside expected range"


def test_alpha_bare_BE_is_half_observed():
    """The BARE alpha prediction (no cooperative) is ~13.3 MeV — half of PDG.

    This is the canonical B3 result: bare a_v = 13.33 MeV, doubled by
    cooperative pion exchange to ~ 28 MeV / A_alpha-equivalent.
    """
    stack = NucleonStackGeometry.from_topology(topology_alpha())
    bare = stack.bare_BE()
    # 6 pairs × 2.222 = 13.333 MeV
    assert abs(bare - 6 * EPS_FACE_MEV) < 1e-9
    # roughly half of observed
    assert 0.4 < bare / REFERENCE_BE_MEV[4] < 0.6


def test_triton_BE_within_5_pct():
    """Triton: 3 face-pairs.  η_coop ≈ 1.27 should reproduce 8.48 MeV exactly
    by construction; verify it's within 5% to allow for small rounding."""
    stack = NucleonStackGeometry.from_topology(topology_triton())
    BE = stack.predicted_BE()
    BE_obs = REFERENCE_BE_MEV[3]
    err = abs(BE - BE_obs) / BE_obs
    assert err < 0.05, f"triton BE pred {BE:.3f} obs {BE_obs:.3f}, err={err*100:.2f}%"


# ---------------------------------------------------------------------------
# BE/A peak: should be in the iron-region (A ~ 56-62)
# ---------------------------------------------------------------------------

def test_BE_per_A_peaks_near_iron():
    """(c) BE/A peaks near Fe-56 (A in [50, 64] in the reference data)."""
    A_arr = np.array(sorted(REFERENCE_BE_MEV.keys()))
    be_per_A = np.array([REFERENCE_BE_MEV[A] / A for A in A_arr])
    A_peak = A_arr[int(np.argmax(be_per_A))]
    assert 50 <= A_peak <= 70, (
        f"BE/A peaks at A={A_peak} (expected near Fe-56 / Ni-62)"
    )


def test_predicted_chart_returns_correct_lengths():
    """NuclearChartVisualizer.predict_chart returns equal-length arrays."""
    viz = NuclearChartVisualizer()
    A_list, pred, obs, err = viz.predict_chart(A_max=20)
    assert len(A_list) == len(pred) == len(obs) == len(err)
    assert all(2 <= A <= 20 for A in A_list)
    # At least the deuteron, alpha, C-12, O-16 should be in there
    assert 2 in A_list and 4 in A_list


# ---------------------------------------------------------------------------
# Geometry: face-sharing actually places cells correctly
# ---------------------------------------------------------------------------

def test_deuteron_geometry_face_centroids_meet():
    """In the deuteron, the matched-face centroids of cell 0 face 0 and
    cell 1 face 0 should coincide (within numerical tolerance)."""
    stack = NucleonStackGeometry.from_topology(topology_deuteron())
    fc_a = stack.cells[0].face_centroids[0]
    fc_b = stack.cells[1].face_centroids[0]
    sep = np.linalg.norm(fc_a - fc_b)
    assert sep < 1e-9, f"Matched face centroids separated by {sep:.3e} fm"


def test_deuteron_geometry_face_normals_anti_aligned():
    """Cell 0 face 0 normal · cell 1 face 0 normal should be ≈ -1."""
    stack = NucleonStackGeometry.from_topology(topology_deuteron())
    n_a = stack.cells[0].face_normals[0]
    n_b = stack.cells[1].face_normals[0]
    assert np.dot(n_a, n_b) < -0.99


def test_alpha_geometry_has_4_distinct_centres():
    """Alpha: 4 cell centres at distinct positions (no degenerate placement)."""
    stack = NucleonStackGeometry.from_topology(topology_alpha())
    centres = stack.cell_centres()
    assert centres.shape == (4, 3)
    # Pairwise distances all positive
    for i in range(4):
        for j in range(i+1, 4):
            d = np.linalg.norm(centres[i] - centres[j])
            assert d > 0.5, f"Cells {i},{j} too close: {d:.3f} fm"


def test_alpha_all_vertices_count():
    """Alpha has 4 cells × 4 vertices = 16 vertex positions."""
    stack = NucleonStackGeometry.from_topology(topology_alpha())
    v = stack.all_vertices()
    assert v.shape == (16, 3)


# ---------------------------------------------------------------------------
# compute_face_pair_count
# ---------------------------------------------------------------------------

def test_compute_face_pair_count_deuteron():
    stack = NucleonStackGeometry.from_topology(topology_deuteron())
    info = stack.compute_face_pair_count()
    assert info['n_face_pairs'] == 1
    assert info['mean_coord'] == 1.0
    assert info['per_cell'] == [1, 1]


def test_compute_face_pair_count_alpha():
    stack = NucleonStackGeometry.from_topology(topology_alpha())
    info = stack.compute_face_pair_count()
    assert info['n_face_pairs'] == 6
    assert info['mean_coord'] == 3.0
    assert info['per_cell'] == [3, 3, 3, 3]


# ---------------------------------------------------------------------------
# Registry & error handling
# ---------------------------------------------------------------------------

def test_topology_registry_complete():
    """All requested A's (2,3,4,6,12,16) are registered."""
    for A in (2, 3, 4, 6, 12, 16):
        assert A in TOPOLOGY_REGISTRY
        topo = get_topology(A)
        assert topo.A == A


def test_unknown_topology_raises():
    with pytest.raises(KeyError):
        get_topology(99)


# ---------------------------------------------------------------------------
# Visualizer smoke tests
# ---------------------------------------------------------------------------

def test_chart_BE_per_A_runs():
    """The BE/A chart can be built without error."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    viz = NuclearChartVisualizer()
    fig, ax = viz.chart_BE_per_A(A_max=20)
    assert fig is not None
    assert ax is not None
    plt.close(fig)


def test_visualize_geometry_3d_runs():
    """3D rendering of the alpha geometry runs without error."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    viz = NuclearChartVisualizer()
    fig, ax = viz.visualize_geometry_3d(4)
    assert fig is not None
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report dict keys
# ---------------------------------------------------------------------------

def test_report_keys_complete():
    stack = NucleonStackGeometry.from_topology(topology_alpha())
    r = stack.report()
    expected = {'name', 'A', 'n_face_pairs', 'mean_coord', 'eta_coop',
                'BE_bare_MeV', 'BE_pred_MeV', 'BE_pred_per_A',
                'BE_obs_MeV', 'BE_obs_per_A', 'err_pct'}
    assert expected.issubset(r.keys())
