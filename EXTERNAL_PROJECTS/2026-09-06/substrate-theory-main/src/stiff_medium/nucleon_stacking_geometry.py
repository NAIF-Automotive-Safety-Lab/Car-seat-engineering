"""Multi-nucleon K_4 cell stacking — explicit geometry of light nuclei.

This module makes the *multi-nucleon* assembly of K_4 cells explicit, runnable
and visualisable.  Where ``k4_face_pair_geometry.py`` models a single matched
face pair (the deuteron), this module models the full face-sharing graph of
N ≥ 2 cells: deuteron (2), triton (3, triangle), alpha (4, tetrahedron),
Li-6 (6, octahedron), C-12 (12, icosahedron), O-16 (16, FCC slab).

THE PHYSICS
-----------
A nucleon ↔ one K_4 cell.  Each cell has 4 face-binding handles (the 4
triangular faces).  When two cells share a face, they bind by ε_face = 2.222
MeV — the deuteron BE, fixed by Λ_QCD/(n_A · N_BAM) = 200/90.

For N face-shared cells we count the face-pair graph:

    P(N) = number of (i,j) cell-pairs sharing a face (i < j)

The *bare* binding energy is

    BE_bare(N) = P(N) · ε_face                                    (1)

This already reproduces the deuteron exactly: P(2)=1 ⇒ 2.222 MeV.

In larger stacks each shared face sits inside a *cooperative* environment
(neighbouring shared faces enhance the local stress-tensor coupling — the
same physics behind the 15.5 MeV/A volume term in the SEMF being roughly
2× the bare a_v = 13.33 MeV the geometry produces).  We absorb this into
a single dimensionless cooperative factor

    BE(N) = η_coop(N) · P(N) · ε_face                             (2)

where the leading-order η_coop(N) for tightly-packed stacks is fixed by
the alpha binding energy: η_coop(α) ≈ 28.30/13.33 ≈ 2.12.  We use
η_coop = 1 for the deuteron (no neighbours to cooperate with) and
η_coop ≈ 2.12 for N ≥ 4 (saturated coordination).  The triton sits in
between (1.5, intermediate coordination).

SHAPES (face-sharing graphs)
----------------------------
- N=2 deuteron: 2 cells share 1 face — line.
- N=3 triton:   3 cells in a triangle (each cell shares 2 faces with the
                other two) — 3 face-pairs.
- N=4 alpha:    4 cells in a tetrahedron (each shares 3 faces) — 6 face-pairs.
- N=6 Li-6:     two alpha-like clusters joined; ~9 face-pairs (octahedral).
- N=12 C-12:    3 triangle stacks (3α model) ≈ 18 face-pairs.
- N=16 O-16:    4 alpha stacks ≈ 24 face-pairs.

All cells are placed by *propagating* face-matching: starting from one
seed cell, each subsequent cell is positioned so that one of its faces
is centroid-aligned and normal-anti-aligned with a free face of an
already-placed cell.  This gives an explicit 3D point cloud of all 4N
vertices, suitable for plotting.

UNITS
-----
Length: fm.  Energy: MeV.

REFERENCES
----------
- k4_face_pair_geometry.py     — single-pair geometry (deuteron).
- multi_cell_k4_dynamics.py    — full Newton dynamics for N cells.
- nuclear_chart.py             — SEMF-style BE/A predictions.
- nbam_from_k4_edges.py        — derivation ε_face = Λ_QCD/(n_A · N_BAM).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from src.stiff_medium.k4_face_pair_geometry import (
    EPS_FACE_MEV,
    K4Geometry,
    rotation_aligning,
    rotation_matrix,
)


# ---------------------------------------------------------------------------
# Cooperative-binding factor η_coop(N)
# ---------------------------------------------------------------------------
# Calibrated so the alpha (N=4, P=6) reproduces its experimental BE.
# η_coop * P * ε_face = 28.295 MeV  ⇒  η_coop = 28.295 / (6 * 2.222) ≈ 2.122
ETA_COOP_ALPHA: float = 28.295674 / (6.0 * EPS_FACE_MEV)  # ≈ 2.122
ETA_COOP_TRITON: float = 8.482 / (3.0 * EPS_FACE_MEV)     # ≈ 1.272
ETA_COOP_DEUTERON: float = 1.0                              # by construction


# Reference (AME2020 / PDG) total binding energies for the chart.
REFERENCE_BE_MEV: Dict[int, float] = {
    2:    2.224573,    # deuteron (¹H + n)
    3:    8.481798,    # triton (³H)
    4:   28.295674,    # alpha (⁴He)
    6:   31.994561,    # Li-6
    8:   56.498,       # Be-8 (≈ 2 × alpha, near-degenerate)
   12:   92.161728,    # C-12
   16:  127.619337,    # O-16
   20:  160.644859,    # Ne-20
   24:  198.257,       # Mg-24
   28:  236.536,       # Si-28
   32:  271.781,       # S-32
   40:  342.052,       # Ca-40
   56:  492.253892,    # Fe-56
   62:  545.259,       # Ni-62 (peak BE/A ≈ 8.7945 MeV)
}


# ---------------------------------------------------------------------------
# Face-sharing graph topology per nucleus
# ---------------------------------------------------------------------------
# (N_cells, list of (cell_a, face_a, cell_b, face_b)) pairs that share a face.
# These shapes are the "platonic" face-sharing graphs used by the geometry
# class — they fix P(N) (the face-pair count) without ambiguity.

@dataclass
class StackingTopology:
    """Combinatorial face-sharing topology for an N-cell stack.

    Attributes
    ----------
    name : str
        Human-readable label ("deuteron", "alpha", ...).
    A : int
        Mass number = number of K_4 cells.
    pairs : list of (i, fi, j, fj)
        Face-sharing pairs.  Cell i face fi shares with cell j face fj.
        Stored with i < j.  This list defines P(N) = len(pairs).
    """

    name: str
    A: int
    pairs: List[Tuple[int, int, int, int]]

    @property
    def n_face_pairs(self) -> int:
        return len(self.pairs)

    def coordination(self, cell: int) -> int:
        """Number of faces of `cell` that are shared with another cell."""
        return sum(1 for (i, _fi, j, _fj) in self.pairs if i == cell or j == cell)

    def mean_coordination(self) -> float:
        """Average shared-face count per cell.  Equals 2·P/N."""
        return 2.0 * self.n_face_pairs / self.A


def topology_deuteron() -> StackingTopology:
    """2 cells share 1 face."""
    return StackingTopology("deuteron", 2, [(0, 0, 1, 0)])


def topology_triton() -> StackingTopology:
    """3 cells in a triangle: each cell shares 2 faces ⇒ 3 face-pairs.

    Each cell has 4 faces; the triangle uses 2 shared faces per cell,
    leaving 2 free faces per cell.
    """
    pairs = [
        (0, 0, 1, 0),
        (1, 1, 2, 0),
        (0, 1, 2, 1),
    ]
    return StackingTopology("triton", 3, pairs)


def topology_alpha() -> StackingTopology:
    """4 cells in a tetrahedron: each cell shares 3 faces ⇒ 6 face-pairs.

    All 4 cells are mutually face-sharing.  Each cell saturates 3 of its 4
    faces; the 4th face is free (outward-pointing).
    """
    pairs = [
        (0, 0, 1, 0),
        (0, 1, 2, 0),
        (0, 2, 3, 0),
        (1, 1, 2, 1),
        (1, 2, 3, 1),
        (2, 2, 3, 2),
    ]
    return StackingTopology("alpha", 4, pairs)


def topology_li6() -> StackingTopology:
    """6 cells = 1 alpha cluster (cells 0..3, 6 pairs) + 1 deuteron
    (cells 4,5, 1 pair) + 1 bridge (alpha-tip face 3 ↔ deuteron face 1)
    = 8 face-pairs total.

    Reflects the empirical α + d cluster structure of Li-6.
    """
    pairs = [
        (0, 0, 1, 0), (0, 1, 2, 0), (0, 2, 3, 0),
        (1, 1, 2, 1), (1, 2, 3, 1), (2, 2, 3, 2),     # alpha cluster
        (4, 0, 5, 0),                                   # deuteron
        (3, 3, 4, 1),                                   # α–d bridge
    ]
    return StackingTopology("Li-6", 6, pairs)


def topology_be8() -> StackingTopology:
    """8 cells = 2 alpha clusters touching on 1 face ⇒ 13 face-pairs."""
    pairs = []
    # Cluster A: 0..3
    for (i, fi, j, fj) in [
        (0, 0, 1, 0), (0, 1, 2, 0), (0, 2, 3, 0),
        (1, 1, 2, 1), (1, 2, 3, 1), (2, 2, 3, 2),
    ]:
        pairs.append((i, fi, j, fj))
    # Cluster B: 4..7
    for (i, fi, j, fj) in [
        (4, 0, 5, 0), (4, 1, 6, 0), (4, 2, 7, 0),
        (5, 1, 6, 1), (5, 2, 7, 1), (6, 2, 7, 2),
    ]:
        pairs.append((i, fi, j, fj))
    # Inter-cluster bridge (cell 3 face 3 shares with cell 4 face 3)
    pairs.append((3, 3, 4, 3))
    return StackingTopology("Be-8", 8, pairs)


def topology_c12() -> StackingTopology:
    """12 cells = 3 alpha clusters in a triangle ⇒ 21 face-pairs.

    Hoyle-state-like 3α geometry.  3 alpha clusters (6 pairs each = 18)
    + 3 inter-cluster bridges = 21.  Each cell on the bridge has its
    shared face count bumped from 3 to 4 (saturated).
    """
    pairs = []
    for k in range(3):
        base = 4 * k
        for (i, fi, j, fj) in [
            (0, 0, 1, 0), (0, 1, 2, 0), (0, 2, 3, 0),
            (1, 1, 2, 1), (1, 2, 3, 1), (2, 2, 3, 2),
        ]:
            pairs.append((base + i, fi, base + j, fj))
    # Inter-cluster bridges: cluster k face 3 of one tip cell to next cluster
    pairs.append((3, 3, 4, 3))     # cluster 0 → 1
    pairs.append((7, 3, 8, 3))     # cluster 1 → 2
    pairs.append((11, 3, 0, 3))    # cluster 2 → 0 (closes the triangle)
    return StackingTopology("C-12", 12, pairs)


def topology_o16() -> StackingTopology:
    """16 cells = 4 alpha clusters in a tetrahedron ⇒ 30 face-pairs.

    Each alpha cluster has 6 internal pairs (24 total).  6 inter-cluster
    bridges (one per pair of clusters in the tetrahedron) ⇒ 30.
    """
    pairs = []
    for k in range(4):
        base = 4 * k
        for (i, fi, j, fj) in [
            (0, 0, 1, 0), (0, 1, 2, 0), (0, 2, 3, 0),
            (1, 1, 2, 1), (1, 2, 3, 1), (2, 2, 3, 2),
        ]:
            pairs.append((base + i, fi, base + j, fj))
    # Six inter-cluster bridges (clusters in tetrahedral arrangement)
    bridges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for (ka, kb) in bridges:
        # Use the 'tip' cell (index 3) of each cluster, face 3
        pairs.append((4 * ka + 3, 3, 4 * kb + 3, 3))
    return StackingTopology("O-16", 16, pairs)


# Registry of named topologies.
TOPOLOGY_REGISTRY = {
    2:  topology_deuteron,
    3:  topology_triton,
    4:  topology_alpha,
    6:  topology_li6,
    8:  topology_be8,
    12: topology_c12,
    16: topology_o16,
}


def get_topology(A: int) -> StackingTopology:
    """Return the canonical face-sharing topology for mass number A."""
    if A in TOPOLOGY_REGISTRY:
        return TOPOLOGY_REGISTRY[A]()
    raise KeyError(f"No registered topology for A={A}.  "
                   f"Known: {sorted(TOPOLOGY_REGISTRY.keys())}")


# ---------------------------------------------------------------------------
# Alpha-cluster decomposition for light nuclei
# ---------------------------------------------------------------------------
#
# Several light nuclei are EMPIRICALLY known to consist of pre-formed
# alpha-cluster substructures plus loose (deuteron, triton, neutron)
# valence pieces.  The bare ``η_coop * P * ε_face`` close-packed model is
# a poor fit for these cases because it assumes saturated coordination
# throughout, while the cluster geometry has BOTH tight intra-cluster
# bonds (alpha core, η ≈ 2.12) AND weak inter-cluster bridges (deuteron
# regime, η ≈ 1.0).
#
# For each cluster nucleus we record:
#   * subclusters    : list of dict with name, A, count
#   * n_aa_bridges_saturated : int — number of α-α bridges that sit in a
#                                     saturated environment (3-α triangle
#                                     or higher), worth η_alpha · ε_face each
#   * n_aa_bridges_unsaturated : int — number of α-α bridges that are loose
#                                       (Be-8 like), worth ε_face each
#   * n_other_bridges : int — number of bridges between alpha and any
#                              non-alpha cluster, worth ε_face each
#
# The cluster-aware BE is
#
#     BE = sum_subcluster(BE_subcluster)
#        + n_aa_bridges_saturated   · η_alpha · ε_face
#        + (n_aa_bridges_unsaturated + n_other_bridges) · ε_face
#
# where BE_alpha = 28.295674 MeV, BE_triton = 8.481798 MeV, BE_d = 2.222 MeV,
# and BE_n = 0 (free neutron has no binding partner).
#
# The bridge counts are minimal spanning-tree style: enough to connect
# all subclusters, plus α-α-α triangular closures when 3+ alphas allow
# saturation.  This reproduces the empirical α-cluster picture (Ikeda
# diagram) and matches AME2020 binding energies within ~5%.

@dataclass(frozen=True)
class SubCluster:
    """One pre-formed nuclear sub-cluster.

    Attributes
    ----------
    name : str
        Label, e.g. "alpha", "triton", "deuteron", "neutron".
    A : int
        Mass number of this subcluster.
    BE_MeV : float
        Binding energy of the isolated subcluster (MeV, positive).
    """
    name: str
    A: int
    BE_MeV: float


# Canonical subclusters (BE values from REFERENCE_BE_MEV).
SUB_ALPHA   = SubCluster(name="alpha",    A=4, BE_MeV=28.295674)
SUB_TRITON  = SubCluster(name="triton",   A=3, BE_MeV=8.481798)
SUB_DEUTERON = SubCluster(name="deuteron", A=2, BE_MeV=2.224573)
SUB_NEUTRON = SubCluster(name="neutron",  A=1, BE_MeV=0.0)


@dataclass(frozen=True)
class ClusterTopology:
    """Alpha-cluster decomposition for a light nucleus.

    Bridge bookkeeping:
      * n_aa_bridges_saturated  : α-α bonds in a saturated triangle/tetra
        configuration (3+ alphas mutually bonded), worth η_alpha · ε_face.
      * n_aa_bridges_unsaturated : isolated α-α bonds (Be-8 like),
        worth ε_face.
      * n_other_bridges          : α to triton/deuteron/neutron bonds,
        worth ε_face.
    """
    name: str
    Z: int
    A: int
    subclusters: Tuple[SubCluster, ...]
    n_aa_bridges_saturated: int = 0
    n_aa_bridges_unsaturated: int = 0
    n_other_bridges: int = 0

    @property
    def n_bridges_total(self) -> int:
        return (self.n_aa_bridges_saturated
                + self.n_aa_bridges_unsaturated
                + self.n_other_bridges)

    @property
    def effective_face_pair_count(self) -> int:
        """Equivalent face-pair count summing all internal + bridge pairs.

        For each subcluster of mass A_s we count the canonical face-pair
        count from TOPOLOGY_REGISTRY (alpha=6, triton=3, deuteron=1,
        neutron=0), then add bridges as 1 face-pair each.

        This is the "P" used in cluster-mode reporting; the actual BE is
        computed via :func:`cluster_aware_BE_MeV` below.
        """
        face_pairs_per_sub = {"alpha": 6, "triton": 3, "deuteron": 1, "neutron": 0}
        sub_count = sum(face_pairs_per_sub[s.name] for s in self.subclusters)
        return sub_count + self.n_bridges_total


# ---- explicit cluster topologies ---------------------------------------

def topology_7li_cluster() -> ClusterTopology:
    """⁷Li = α + ³H, joined by 1 α–³H bridge (η = 1, deuteron-like).

    BE = BE(α) + BE(³H) + 1 · ε_face = 28.30 + 8.48 + 2.22 = 39.00 MeV
    AME2020 = 39.245 MeV  (residual −0.6 %)
    """
    return ClusterTopology(
        name="7Li", Z=3, A=7,
        subclusters=(SUB_ALPHA, SUB_TRITON),
        n_other_bridges=1,
    )


def topology_9be_cluster() -> ClusterTopology:
    """⁹Be = 2α + n, joined by 1 α–n bridge (η = 1).

    The two alphas form a Be-8-like core (Be-8 itself is unbound by
    92 keV, so we treat the α-α direct bond as "absent" — it adds no
    extra binding) and the loose neutron sits between them with one
    face-pair shared.

    BE = 2·BE(α) + 1 · ε_face = 56.59 + 2.22 = 58.81 MeV
    AME2020 = 58.165 MeV  (residual +1.1 %)
    """
    return ClusterTopology(
        name="9Be", Z=4, A=9,
        subclusters=(SUB_ALPHA, SUB_ALPHA, SUB_NEUTRON),
        n_other_bridges=1,
    )


def topology_10b_cluster() -> ClusterTopology:
    """¹⁰B = 2α + d, joined by 1 α–α (unsaturated, η = 1) and 2 α–d bridges.

    BE = 2·BE(α) + BE(d) + 1 · ε_face (α–α) + 2 · ε_face (α–d)
       = 56.59 + 2.22 + 2.22 + 4.44 = 65.48 MeV
    AME2020 = 64.751 MeV  (residual +1.1 %)
    """
    return ClusterTopology(
        name="10B", Z=5, A=10,
        subclusters=(SUB_ALPHA, SUB_ALPHA, SUB_DEUTERON),
        n_aa_bridges_unsaturated=1,
        n_other_bridges=2,
    )


def topology_14n_cluster() -> ClusterTopology:
    """¹⁴N = 3α + d, with the 3 alphas forming a saturated triangle
    (η_alpha applies to the 3 α-α bonds) plus 1 α–d bridge (η = 1).

    BE = 3·BE(α) + BE(d) + 3 · η_alpha · ε_face + 1 · ε_face
       = 84.89 + 2.22 + 14.14 + 2.22 = 103.48 MeV
    AME2020 = 104.659 MeV  (residual −1.1 %)
    """
    return ClusterTopology(
        name="14N", Z=7, A=14,
        subclusters=(SUB_ALPHA, SUB_ALPHA, SUB_ALPHA, SUB_DEUTERON),
        n_aa_bridges_saturated=3,
        n_other_bridges=1,
    )


CLUSTER_TOPOLOGIES: Dict[Tuple[int, int], ClusterTopology] = {
    (3, 7):  topology_7li_cluster(),
    (4, 9):  topology_9be_cluster(),
    (5, 10): topology_10b_cluster(),
    (7, 14): topology_14n_cluster(),
}


def get_cluster_topology(Z: int, A: int) -> Optional[ClusterTopology]:
    """Return the registered cluster decomposition for (Z, A), or None."""
    return CLUSTER_TOPOLOGIES.get((Z, A))


def cluster_face_pair_count(Z: int, A: int) -> Optional[int]:
    """Equivalent face-pair count P for a cluster nucleus, or None.

    Returns the sum of (sub-cluster internal face-pairs) + (inter-cluster
    bridge face-pairs).  Used for cluster-aware reporting; the actual BE
    is computed via :func:`cluster_aware_BE_MeV`.
    """
    topo = get_cluster_topology(Z, A)
    if topo is None:
        return None
    return topo.effective_face_pair_count


def cluster_aware_BE_MeV(Z: int, A: int) -> Optional[float]:
    """Cluster-aware binding energy for (Z, A), or None if no topology.

        BE = sum_subcluster(BE_subcluster)
           + n_aa_saturated · η_alpha · ε_face
           + (n_aa_unsaturated + n_other) · ε_face

    The α-α saturated bridges (η_alpha) reflect 3+ mutually-bonded alphas
    (triangle / tetrahedron); unsaturated α-α bridges (η = 1) reflect
    isolated Be-8-like pairs.  The "other" bridges (α–triton, α–deuteron,
    α–neutron) all use η = 1 because the non-alpha partner has no
    cooperative environment of its own.
    """
    topo = get_cluster_topology(Z, A)
    if topo is None:
        return None
    be_sub = sum(s.BE_MeV for s in topo.subclusters)
    be_aa_sat = topo.n_aa_bridges_saturated * ETA_COOP_ALPHA * EPS_FACE_MEV
    be_other = ((topo.n_aa_bridges_unsaturated + topo.n_other_bridges)
                * 1.0 * EPS_FACE_MEV)
    return be_sub + be_aa_sat + be_other


def cooperative_factor(A: int) -> float:
    """Empirical cooperative-binding factor η_coop(A).

    Calibrated so:
      A=2  (deuteron) ⇒ η = 1   (no cooperative neighbours)
      A=3  (triton)   ⇒ η ≈ 1.27 (intermediate)
      A=4  (alpha)    ⇒ η ≈ 2.12 (saturated coordination)
      A≥6 (cluster)  ⇒ scales with mean coordination, capped near alpha value
    """
    if A == 2:
        return ETA_COOP_DEUTERON
    if A == 3:
        return ETA_COOP_TRITON
    if A == 4:
        return ETA_COOP_ALPHA
    # For A>=6, η_coop SATURATES at the alpha value.  Additional binding
    # in cluster nuclei (C-12, O-16, ...) comes from inter-cluster face
    # pairs (already counted in the topology), not from enhanced per-face
    # coupling.  The empirical BE/A then falls naturally at A>=6 because
    # the inter-cluster "bridge" pairs have lower effective η than the
    # tight intra-cluster pairs — captured here by the simple cap.
    return ETA_COOP_ALPHA


# ---------------------------------------------------------------------------
# Geometry: place N K_4 cells in face-sharing world-frame configuration
# ---------------------------------------------------------------------------

@dataclass
class NucleonStackGeometry:
    """Place N K_4 cells in a face-sharing 3D configuration.

    Workflow:
        1. ``stack = NucleonStackGeometry.from_topology(topology_alpha())``
        2. ``stack.cells`` — list of N K4Geometry instances in world frame
        3. ``stack.predicted_BE()`` — N · c / 2 · ε_face  (with η_coop)
        4. ``stack.compute_face_pair_count()`` — coordination per cell
    """

    topology: StackingTopology
    cells: List[K4Geometry] = field(default_factory=list)
    eps_face: float = EPS_FACE_MEV
    use_cooperative: bool = True

    # ----- construction -----
    @classmethod
    def from_topology(cls, topology: StackingTopology,
                      eps_face: float = EPS_FACE_MEV,
                      use_cooperative: bool = True
                      ) -> "NucleonStackGeometry":
        """Place N K_4 cells in 3D so the topology's face-sharing pairs are
        actually centroid-aligned in world coordinates."""
        stack = cls(topology=topology, eps_face=eps_face,
                    use_cooperative=use_cooperative)
        stack._build_geometry()
        return stack

    def _build_geometry(self) -> None:
        """Iteratively place cells.  Cell 0 sits at the origin in body
        frame; each subsequent cell is placed face-matched against an
        already-placed neighbour.

        We use a *first-pair-from-list* strategy: for cell j, find the
        first (i, fi, j, fj) pair in the topology where i < j, and place
        cell j so that its face fj is centroid-matched and
        normal-anti-aligned with cell i's face fi.

        This gives a unique world-frame placement (modulo torsion which
        is fixed to 0).  Geometric residuals (other shared faces in the
        topology that DON'T quite line up after this seeded placement)
        are fine for visualisation — what matters for predicted_BE is
        the topological count, not the world-frame fit.
        """
        N = self.topology.A
        d_centres = 2.0 * K4Geometry().d_face   # face-to-face touching
        self.cells = []

        # Cell 0: identity placement
        seed = K4Geometry()
        # Wrap into world-frame transformed copy (identity rot, zero shift)
        self.cells.append(seed.transformed(np.zeros(3), np.eye(3)))

        for j in range(1, N):
            placed = False
            for (i, fi, jj, fj) in self.topology.pairs:
                if jj == j and i < j:
                    self._place_cell_against(j, i, fi, fj, d_centres)
                    placed = True
                    break
                if i == j and jj < j:
                    # Symmetric: this pair has cell j first; place against jj
                    self._place_cell_against(j, jj, fj, fi, d_centres)
                    placed = True
                    break
            if not placed:
                # Cell unconstrained by listed pairs — place along +x as fallback
                offset = np.array([d_centres * (j + 1), 0.0, 0.0])
                fresh = K4Geometry()
                self.cells.append(fresh.transformed(offset, np.eye(3)))

    def _place_cell_against(self, j: int, i: int,
                            face_i: int, face_j: int,
                            d_centres: float) -> None:
        """Place cell j so that its face_j is normal-anti-aligned with
        cell i's face_i and centroid-touching."""
        # World-frame normal of cell i's face_i (already placed)
        n_i_world = self.cells[i].face_normals[face_i]
        c_i_world = self.cells[i].face_centroids[face_i]

        # Body-frame normal of fresh cell, face face_j: must align to -n_i_world
        body = K4Geometry()
        body_n_j = body.face_normals[face_j]
        target_normal = -n_i_world
        rot_j = rotation_aligning(body_n_j, target_normal)

        # After rotation, body face centroid moves to (rot_j @ body_c_j).
        # We want this to coincide with c_i_world; therefore translate by
        # (c_i_world - rot_j @ body_c_j).
        body_c_j = body.face_centroids[face_j]
        translation = c_i_world - rot_j @ body_c_j

        cell_j = body.transformed(translation, rot_j)
        # Append (or replace if we've already placed)
        if j < len(self.cells):
            self.cells[j] = cell_j
        else:
            # Pad with placeholders if topology pairs out of order
            while len(self.cells) < j:
                self.cells.append(K4Geometry().transformed(
                    np.array([d_centres * (len(self.cells) + 1), 0, 0]),
                    np.eye(3)))
            self.cells.append(cell_j)

    # ----- counts -----
    def compute_face_pair_count(self) -> Dict[str, float]:
        """Coordination data for the stack.

        Returns
        -------
        dict with:
            n_face_pairs : int   — total P(N) = len(pairs)
            mean_coord   : float — 2·P/N (avg shared faces per cell)
            per_cell     : list of int  — shared-face count per cell index
        """
        topo = self.topology
        per_cell = [topo.coordination(k) for k in range(topo.A)]
        return {
            'n_face_pairs': topo.n_face_pairs,
            'mean_coord': topo.mean_coordination(),
            'per_cell': per_cell,
        }

    # ----- binding energy -----
    def predicted_BE(self) -> float:
        """Predicted total binding energy (MeV).

            BE = η_coop(A) · P(N) · ε_face

        where P(N) = ∑ (shared faces per cell) / 2 = n_face_pairs.
        """
        P = self.topology.n_face_pairs
        eta = cooperative_factor(self.topology.A) if self.use_cooperative else 1.0
        return eta * P * self.eps_face

    def predicted_BE_per_A(self) -> float:
        """BE / A (MeV per nucleon)."""
        return self.predicted_BE() / max(self.topology.A, 1)

    def bare_BE(self) -> float:
        """Bare prediction (no cooperative factor): P · ε_face."""
        return self.topology.n_face_pairs * self.eps_face

    def report(self) -> Dict[str, float]:
        """Full diagnostic dict for printing or testing."""
        topo = self.topology
        ref = REFERENCE_BE_MEV.get(topo.A, float('nan'))
        be = self.predicted_BE()
        return {
            'name': topo.name,
            'A': topo.A,
            'n_face_pairs': topo.n_face_pairs,
            'mean_coord': topo.mean_coordination(),
            'eta_coop': cooperative_factor(topo.A),
            'BE_bare_MeV': self.bare_BE(),
            'BE_pred_MeV': be,
            'BE_pred_per_A': be / topo.A,
            'BE_obs_MeV': ref,
            'BE_obs_per_A': ref / topo.A if not np.isnan(ref) else float('nan'),
            'err_pct': 100.0 * (be - ref) / ref if not np.isnan(ref) else float('nan'),
        }

    # ----- world-frame point cloud -----
    def all_vertices(self) -> NDArray[np.float64]:
        """All 4N vertex positions stacked into a (4N, 3) array."""
        return np.vstack([c.vertices for c in self.cells])

    def all_face_centroids(self) -> NDArray[np.float64]:
        """All 4N face centroids stacked into a (4N, 3) array."""
        return np.vstack([c.face_centroids for c in self.cells])

    def cell_centres(self) -> NDArray[np.float64]:
        """Cell centroids: (N, 3) array."""
        return np.array([c.centre for c in self.cells])


# ---------------------------------------------------------------------------
# Nuclear chart visualizer
# ---------------------------------------------------------------------------

@dataclass
class NuclearChartVisualizer:
    """Plots BE/A vs A and 3D stacking geometry for the predicted nuclei."""

    topologies: Dict[int, StackingTopology] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.topologies:
            self.topologies = {A: get_topology(A)
                               for A in TOPOLOGY_REGISTRY.keys()}

    def predict_chart(self, A_max: int = 20
                      ) -> Tuple[List[int], List[float], List[float], List[float]]:
        """Return (A_list, BE_pred_per_A, BE_obs_per_A, err_pct).

        Includes only A's that have a registered topology (and ≤ A_max).
        Plus extrapolated entries for the BE/A peak comparison.
        """
        A_list, pred, obs, err = [], [], [], []
        for A in sorted(self.topologies.keys()):
            if A > A_max:
                continue
            stack = NucleonStackGeometry.from_topology(self.topologies[A])
            r = stack.report()
            A_list.append(A)
            pred.append(r['BE_pred_per_A'])
            obs.append(r['BE_obs_per_A'])
            err.append(r['err_pct'])
        return A_list, pred, obs, err

    def chart_BE_per_A(self, ax=None, A_max: int = 20):
        """Plot BE/A vs A — predicted (cooperative model) vs observed PDG.

        Also overlays the empirical SEMF/BE-curve (extrapolated from the
        Ni-62 peak ≈ 8.79 MeV/A).
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = ax.figure

        # Predicted (this module)
        A_list, pred, obs, err = self.predict_chart(A_max=A_max)
        A_arr = np.array(A_list)
        pred_arr = np.array(pred)
        obs_arr = np.array(obs)

        # Reference observed BE/A across the full PDG range up to A_max
        # (extends past A_max to show the Fe-56 peak)
        obs_full_A = sorted([A for A in REFERENCE_BE_MEV if A <= max(A_max, 62)])
        obs_full_BE_per_A = [REFERENCE_BE_MEV[A] / A for A in obs_full_A]

        ax.plot(obs_full_A, obs_full_BE_per_A, 'ko-', label='PDG / AME2020',
                markersize=7, linewidth=1, alpha=0.85)
        ax.plot(A_arr, pred_arr, 'rs--', label='K_4 stacking (cooperative)',
                markersize=10, linewidth=1.5, alpha=0.9)

        # Mark Fe-56 peak
        if 56 in REFERENCE_BE_MEV:
            be56 = REFERENCE_BE_MEV[56] / 56
            ax.axvline(56, color='gray', linestyle=':', alpha=0.5)
            ax.text(56, be56 + 0.2, 'Fe-56 peak', ha='center', fontsize=9)

        ax.set_xlabel('Mass number A')
        ax.set_ylabel('BE / A  [MeV / nucleon]')
        ax.set_title('Nuclear binding energy per nucleon: K_4 stacking vs PDG\n'
                     f"(ε_face = {EPS_FACE_MEV:.3f} MeV from Λ_QCD/(n_A·N_BAM))")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 10)
        return fig, ax

    def visualize_geometry_3d(self, A: int, ax=None,
                              show_face_normals: bool = False,
                              highlight_shared: bool = True):
        """3D plot of N face-shared K_4 cells for nucleus A."""
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        topo = get_topology(A)
        stack = NucleonStackGeometry.from_topology(topo)

        if ax is None:
            fig = plt.figure(figsize=(9, 8))
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        # Per-cell colors
        cmap = plt.get_cmap('tab10')
        colors = [cmap(i % 10) for i in range(topo.A)]

        for k, cell in enumerate(stack.cells):
            v = cell.vertices
            color = colors[k]
            ax.scatter(v[:, 0], v[:, 1], v[:, 2], c=[color], s=30,
                       edgecolors='black', linewidths=0.5,
                       label=f'cell {k}')
            for (a, b) in cell.edges:
                ax.plot(*zip(v[a], v[b]), color=color, alpha=0.55, lw=1)
            polys = [v[f] for f in cell.faces]
            ax.add_collection3d(Poly3DCollection(
                polys, facecolor=color, edgecolor=color, alpha=0.18))

            if show_face_normals:
                fc, fn = cell.face_centroids, cell.face_normals
                for q in range(4):
                    ax.quiver(fc[q, 0], fc[q, 1], fc[q, 2],
                              fn[q, 0]*0.3, fn[q, 1]*0.3, fn[q, 2]*0.3,
                              color=color, alpha=0.5,
                              arrow_length_ratio=0.3)

        # Highlight shared faces (red)
        if highlight_shared:
            for (i, fi, j, fj) in topo.pairs:
                cell_i = stack.cells[i]
                vi = cell_i.vertices[cell_i.faces[fi]]
                ax.add_collection3d(Poly3DCollection(
                    [vi], facecolor='red', edgecolor='darkred', alpha=0.4))

        ax.set_xlabel('x [fm]'); ax.set_ylabel('y [fm]'); ax.set_zlabel('z [fm]')
        report = stack.report()
        ax.set_title(
            f"{topo.name}: A={topo.A}, P={topo.n_face_pairs} face-pairs\n"
            f"BE_pred = {report['BE_pred_MeV']:.2f} MeV "
            f"(obs = {report['BE_obs_MeV']:.2f}, err = {report['err_pct']:+.1f}%)"
        )
        if topo.A <= 6:
            ax.legend(fontsize=8, loc='upper left')
        return fig, ax


# ---------------------------------------------------------------------------
# Convenience: end-to-end demo
# ---------------------------------------------------------------------------

def demo() -> None:
    """Print a small chart of all registered nuclei."""
    print(f"{'name':10} {'A':>3} {'P':>3} {'coord':>6} "
          f"{'eta':>6} {'BE_pred':>10} {'BE_obs':>10} {'err%':>7}")
    print("-" * 72)
    for A in sorted(TOPOLOGY_REGISTRY.keys()):
        stack = NucleonStackGeometry.from_topology(get_topology(A))
        r = stack.report()
        print(f"{r['name']:10} {r['A']:>3d} {r['n_face_pairs']:>3d} "
              f"{r['mean_coord']:>6.2f} {r['eta_coop']:>6.2f} "
              f"{r['BE_pred_MeV']:>10.3f} {r['BE_obs_MeV']:>10.3f} "
              f"{r['err_pct']:>+7.2f}")


if __name__ == "__main__":
    demo()
