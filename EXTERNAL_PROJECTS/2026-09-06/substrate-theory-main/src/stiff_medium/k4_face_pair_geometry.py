"""K_4 face-pair geometry — explicit, runnable simulation of the deuteron.

This module makes the K_4 tetrahedron face-pair coupling EXPLICITLY VISIBLE
and RUNNABLE.  Rather than burying the geometry inside numerical relaxation
helpers (as in ``multi_cell_k4_dynamics.py``), every geometric primitive is
exposed as an attribute that can be inspected, plotted, or tested.

THE PHYSICS
-----------
A single K_4 cell is a regular tetrahedron whose four triangular faces are
the substrate's "binding handles."  When two cells approach each other so
that one face of cell A is centroid-aligned and normal-anti-aligned with
one face of cell B, the substrate forms a single shared face and the
configuration drops by ε_face = 2.222 MeV — the deuteron binding energy.

The matched-face geometry is parametrised by a single distance ``d``
between cell centroids.  At the optimum d ≈ 1.4 fm the two outward
face-normals point exactly through each other (n̂_A·n̂_B = -1), the face
centroids coincide, and the energy reaches its minimum.

Geometric facts the module verifies explicitly:

    K_4 has 4 vertices, 6 edges, 4 faces.
    The four outward face-normals sum to zero (4-fold symmetry).
    The vertex-centre-vertex angle is the tetrahedral angle 109.47°.
    The angle between two outward face-normals (centroid-pointing) is
    arccos(-1/3) = 109.47° — same as vertex-vertex by tetrahedral
    symmetry, since face_normal_i ∝ -vertex_i.
    The complementary "between adjacent face PLANES" angle is the
    dihedral angle arccos(+1/3) = 70.53°.
    For a regular tetrahedron of edge a:
        d_face = a / (2√6)
        d_vertex = a × √(3/8) = a × √6 / 4

UNITS
-----
Length: fm.  Energy: MeV.  Lengths chosen so that the matched deuteron
sits at d_centres ≈ 1.4 fm.

REFERENCES
----------
multi_cell_k4_dynamics.py — finite-difference Newton dynamics.
deuteron_v2.py — derivation of ε_face = 2.222 MeV from substrate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Physical constants (B3 nuclear sector)
# ---------------------------------------------------------------------------

EPS_FACE_MEV: float = 2.222            # face-pair binding energy = deuteron BE
XI_FM: float = 1.4                     # coupling range (~ pion Compton wavelength)
D_FACE_FM: float = 0.8                 # body-frame face-centroid distance
ALIGN_THRESHOLD_DEG: float = 30.0      # half-angle anti-alignment tolerance
TETRAHEDRAL_ANGLE_DEG: float = 109.47122063449069   # arccos(-1/3) — vertex-centre-vertex AND normal-normal
DIHEDRAL_ANGLE_DEG: float = 70.52877936550931       # arccos(+1/3) — between adjacent face PLANES


# ---------------------------------------------------------------------------
# K_4 geometry (body frame) — explicit attributes
# ---------------------------------------------------------------------------

@dataclass
class K4Geometry:
    """Explicit construction of the regular-tetrahedron K_4 cell.

    All combinatorial and geometric quantities of the K_4 graph
    (vertices, edges, faces, centroids, normals) are exposed as
    attributes.  Construction defaults to a unit-d_face tetrahedron
    centred at the origin in body frame; the world frame can be set
    via ``apply_transform``.

    Attributes
    ----------
    vertices : (4, 3) array
        Vertex positions of the regular tetrahedron.
    edges : (6, 2) integer array
        Vertex-index pairs for the C(4,2) = 6 edges.
    faces : (4, 3) integer array
        Vertex-index triples for the C(4,3) = 4 faces.
    edge_length : float
        Common length of all 6 edges.
    face_centroids : (4, 3) array
        Geometric centroid of each face (mean of its 3 vertices).
    face_normals : (4, 3) array
        Outward unit normals of the 4 faces, pointing from cell
        centre through face centroid.
    """

    vertices: NDArray[np.float64] = field(init=False)
    edges: NDArray[np.int64] = field(init=False)
    faces: NDArray[np.int64] = field(init=False)
    edge_length: float = field(init=False)
    face_centroids: NDArray[np.float64] = field(init=False)
    face_normals: NDArray[np.float64] = field(init=False)
    centre: NDArray[np.float64] = field(init=False)
    d_face: float = D_FACE_FM

    def __post_init__(self) -> None:
        # Edge length scaled so face_centroid sits at distance d_face from centre.
        # For regular tetrahedron: d_face_centroid = a / (2√6)  ⇒  a = 2√6 · d_face
        a = 2.0 * np.sqrt(6.0) * self.d_face
        self.edge_length = a

        # Standard regular-tetrahedron coords (edge length 2√2 in raw units)
        raw = np.array([
            [+1.0, +1.0, +1.0],
            [+1.0, -1.0, -1.0],
            [-1.0, +1.0, -1.0],
            [-1.0, -1.0, +1.0],
        ], dtype=float)
        self.vertices = raw * (a / (2.0 * np.sqrt(2.0)))
        self.centre = self.vertices.mean(axis=0)        # at origin

        # All C(4,2) = 6 edges
        self.edges = np.array([
            [0, 1], [0, 2], [0, 3],
            [1, 2], [1, 3],
            [2, 3],
        ], dtype=np.int64)

        # All C(4,3) = 4 faces (triangle = 3 vertices opposite to one excluded vertex)
        self.faces = np.array([
            [1, 2, 3],     # face opposite vertex 0
            [0, 2, 3],     # face opposite vertex 1
            [0, 1, 3],     # face opposite vertex 2
            [0, 1, 2],     # face opposite vertex 3
        ], dtype=np.int64)

        self.face_centroids = np.array(
            [self.vertices[f].mean(axis=0) for f in self.faces]
        )
        # Outward normal: from cell centre (origin) through face centroid
        norms = np.linalg.norm(self.face_centroids - self.centre, axis=1, keepdims=True)
        self.face_normals = (self.face_centroids - self.centre) / norms

    # ----- combinatorial counts ----------------------------------------
    @property
    def n_vertices(self) -> int:
        return self.vertices.shape[0]

    @property
    def n_edges(self) -> int:
        return self.edges.shape[0]

    @property
    def n_faces(self) -> int:
        return self.faces.shape[0]

    # ----- geometric invariants ----------------------------------------
    def normals_sum(self) -> NDArray[np.float64]:
        """Vector sum of the four outward face-normals (should be ~0)."""
        return self.face_normals.sum(axis=0)

    def vertex_angle(self) -> float:
        """Angle between two vertices viewed from the centre (degrees).

        For a regular tetrahedron centred at the origin, this is the
        tetrahedral angle arccos(-1/3) ≈ 109.47°.
        """
        v0 = self.vertices[0] / np.linalg.norm(self.vertices[0])
        v1 = self.vertices[1] / np.linalg.norm(self.vertices[1])
        return float(np.degrees(np.arccos(np.clip(v0 @ v1, -1.0, 1.0))))

    def normal_pair_angle(self) -> float:
        """Angle between two outward face-normals (degrees).

        For a regular tetrahedron, since face_normal_i ∝ -vertex_i, this
        equals the vertex-vertex angle arccos(-1/3) ≈ 109.47°.
        """
        n0 = self.face_normals[0]
        n1 = self.face_normals[1]
        return float(np.degrees(np.arccos(np.clip(n0 @ n1, -1.0, 1.0))))

    def dihedral_angle(self) -> float:
        """Angle between two face PLANES of the tetrahedron (degrees).

        Computed as the angle between two INWARD-pointing face-normals
        (or equivalently π minus the outward-normal angle).  For a
        regular tetrahedron this is arccos(+1/3) ≈ 70.53°.
        """
        n0 = -self.face_normals[0]
        n1 = self.face_normals[1]
        return float(np.degrees(np.arccos(np.clip(n0 @ n1, -1.0, 1.0))))

    # ----- world-frame transform ---------------------------------------
    def transformed(self, translation: NDArray[np.float64],
                    rotation: NDArray[np.float64]) -> "K4Geometry":
        """Return a new K4Geometry whose vertices are R · v + t.

        The combinatorial data (edges, faces) is unchanged; centroids and
        normals are rotated/translated.
        """
        new = K4Geometry(d_face=self.d_face)
        new.vertices = (self.vertices @ rotation.T) + translation
        new.face_centroids = (self.face_centroids @ rotation.T) + translation
        new.face_normals = self.face_normals @ rotation.T
        new.centre = self.centre + translation
        return new


# ---------------------------------------------------------------------------
# Rotations (axis-angle, Rodrigues)
# ---------------------------------------------------------------------------

def rotation_matrix(axis: NDArray[np.float64], angle: float
                    ) -> NDArray[np.float64]:
    """Rodrigues' formula for a rotation about ``axis`` by ``angle`` (rad)."""
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm < 1e-15:
        return np.eye(3)
    a = axis / norm
    c, s = np.cos(angle), np.sin(angle)
    K = np.array([[0, -a[2], a[1]],
                  [a[2], 0, -a[0]],
                  [-a[1], a[0], 0]])
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def rotation_aligning(v_from: NDArray[np.float64],
                       v_to: NDArray[np.float64]) -> NDArray[np.float64]:
    """Rotation matrix taking unit vector ``v_from`` to unit vector ``v_to``."""
    a = v_from / np.linalg.norm(v_from)
    b = v_to / np.linalg.norm(v_to)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = float(np.dot(a, b))
    if s < 1e-12:
        if c > 0:
            return np.eye(3)
        # 180°: rotate about any perpendicular axis
        perp = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, perp)
        axis = axis / np.linalg.norm(axis)
        return rotation_matrix(axis, np.pi)
    K = np.array([[0, -v[2], v[1]],
                  [v[2], 0, -v[0]],
                  [-v[1], v[0], 0]])
    return np.eye(3) + K + (K @ K) * ((1.0 - c) / (s ** 2))


# ---------------------------------------------------------------------------
# Face-pair coupling (smooth energy of two K_4 cells)
# ---------------------------------------------------------------------------

@dataclass
class FacePairMatch:
    """Result of identifying the best-matched face pair across two cells.

    Attributes
    ----------
    face_a, face_b : int
        Indices (0..3) of the matched faces on cell A and cell B.
    centroid_distance : float
        World-frame distance between matched face centroids.
    normal_dot : float
        n̂_a · n̂_b for the matched pair (should be ≈ -1 for good match).
    score : float
        Combined match score: smaller is better.
    """
    face_a: int
    face_b: int
    centroid_distance: float
    normal_dot: float
    score: float


@dataclass
class K4FacePairCoupling:
    """Two K_4 cells coupled via face-pair binding — the deuteron simulation.

    The two cells are placed at separations along the +x axis, with one
    face of each cell oriented to point toward the other cell.  The
    coupling energy is computed as a smooth function of:

        d   — distance between the two cell centroids
        θ   — torsion angle (rotation of cell B about the x axis)

    At (d ≈ 1.4 fm, θ = 0) the matched pair gives ε_face = 2.222 MeV.
    """

    cell_a: K4Geometry = field(default_factory=K4Geometry)
    cell_b: K4Geometry = field(default_factory=K4Geometry)
    eps_face: float = EPS_FACE_MEV
    xi: float = XI_FM
    align_threshold_deg: float = ALIGN_THRESHOLD_DEG

    # Current world-frame placement (set by place_two_cells)
    pos_a: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(3))
    pos_b: NDArray[np.float64] = field(
        default_factory=lambda: np.array([1.4, 0.0, 0.0]))
    rot_a: NDArray[np.float64] = field(default_factory=lambda: np.eye(3))
    rot_b: NDArray[np.float64] = field(default_factory=lambda: np.eye(3))

    def place_two_cells(self, d_centers: float, orientation: float = 0.0
                        ) -> None:
        """Position cell A and cell B for a matched-face deuteron.

        Cell A: positioned at -d/2 along x, with body face 0 (normal
                  pointing along (+,+,+)) rotated to point along +x.
        Cell B: positioned at +d/2 along x, with body face 0 rotated to
                  point along -x, then optionally rotated by ``orientation``
                  radians about the x axis (torsion).

        After this call, the two cells are face-matched; the matched
        pair has indices (face_a=0, face_b=0) by construction.
        """
        # Body normal of face 0 = first row of self.face_normals
        body_n0 = self.cell_a.face_normals[0]

        # Cell A: face 0 → +x
        self.rot_a = rotation_aligning(body_n0, np.array([+1.0, 0.0, 0.0]))
        self.pos_a = np.array([-d_centers / 2.0, 0.0, 0.0])

        # Cell B: face 0 → -x, then optional torsion about x
        rot_align_b = rotation_aligning(body_n0, np.array([-1.0, 0.0, 0.0]))
        rot_torsion = rotation_matrix(np.array([1.0, 0.0, 0.0]), orientation)
        self.rot_b = rot_torsion @ rot_align_b
        self.pos_b = np.array([+d_centers / 2.0, 0.0, 0.0])

    # ---------- world-frame face data ----------------------------------
    def world_faces_a(self) -> Tuple[NDArray, NDArray]:
        c = self.cell_a.face_centroids @ self.rot_a.T + self.pos_a
        n = self.cell_a.face_normals @ self.rot_a.T
        return c, n

    def world_faces_b(self) -> Tuple[NDArray, NDArray]:
        c = self.cell_b.face_centroids @ self.rot_b.T + self.pos_b
        n = self.cell_b.face_normals @ self.rot_b.T
        return c, n

    def world_vertices_a(self) -> NDArray:
        return self.cell_a.vertices @ self.rot_a.T + self.pos_a

    def world_vertices_b(self) -> NDArray:
        return self.cell_b.vertices @ self.rot_b.T + self.pos_b

    # ---------- face-matching ------------------------------------------
    def find_face_pair_match(self) -> FacePairMatch:
        """Identify which face of A best matches which face of B.

        Best match minimises ``score = d_centroid + λ · (1 + n̂_a · n̂_b)``,
        where the second term is small when the normals are anti-aligned.
        """
        c_a, n_a = self.world_faces_a()
        c_b, n_b = self.world_faces_b()
        diff = c_a[:, None, :] - c_b[None, :, :]
        d = np.linalg.norm(diff, axis=-1)              # (4,4)
        dot = n_a @ n_b.T                              # (4,4)
        # Score: distance + λ × (1 + n_a·n_b).  λ chosen so a perfectly
        # anti-aligned pair (dot=-1) contributes 0 from the alignment term,
        # and a same-aligned pair (dot=+1) contributes 2λ ~ same scale as d.
        lam = 2.0
        score = d + lam * (1.0 + dot)
        a_idx, b_idx = np.unravel_index(np.argmin(score), score.shape)
        return FacePairMatch(
            face_a=int(a_idx),
            face_b=int(b_idx),
            centroid_distance=float(d[a_idx, b_idx]),
            normal_dot=float(dot[a_idx, b_idx]),
            score=float(score[a_idx, b_idx]),
        )

    # ---------- coupling energy ----------------------------------------
    def _gate(self, d: NDArray, dot: NDArray) -> NDArray:
        """Smooth gate: Gaussian × clipped (-dot)² in [0,1]."""
        spatial = np.exp(-(d / self.xi) ** 2)
        cos_thresh = np.cos(np.deg2rad(self.align_threshold_deg))
        align = np.clip(-dot - cos_thresh, 0.0, None) / max(1.0 - cos_thresh, 1e-9)
        return spatial * align ** 2

    def coupling_energy(self, d: Optional[float] = None,
                        theta: Optional[float] = None,
                        sum_all_pairs: bool = False) -> float:
        """Substrate face-pair energy U(d, θ) for the current placement.

        Two K_4 cells share AT MOST one face — the "matched" face pair
        identified by ``find_face_pair_match``.  By default we therefore
        return the contribution of the best-matched face pair only:

            U(d, θ) = -ε_face × g_pair(d_match, θ_match)

        Set ``sum_all_pairs=True`` to instead return the un-truncated
        sum over all 16 face-pair combinations (useful for diagnostics
        when comparing against ``multi_cell_k4_dynamics``).

        If ``d`` and ``theta`` are passed, repositions the cells first;
        otherwise uses the current placement.
        """
        if d is not None or theta is not None:
            d_use = d if d is not None else float(np.linalg.norm(self.pos_b - self.pos_a))
            theta_use = theta if theta is not None else 0.0
            self.place_two_cells(d_use, theta_use)
        c_a, n_a = self.world_faces_a()
        c_b, n_b = self.world_faces_b()
        diff = c_a[:, None, :] - c_b[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        dot = n_a @ n_b.T
        g = self._gate(dist, dot)
        if sum_all_pairs:
            return float(-self.eps_face * g.sum())
        # Single-shared-face binding: take the strongest gated contribution.
        return float(-self.eps_face * g.max())

    # ---------- relaxation ---------------------------------------------
    def relax_to_minimum(self, d_init: float = 2.0,
                         theta_init: float = 0.0,
                         step_d: float = 0.02,
                         step_theta: float = 0.02,
                         max_iter: int = 600,
                         tol: float = 1e-9) -> Tuple[float, float, float]:
        """Gradient descent in (d, θ) to the energy minimum.

        Returns
        -------
        d_min, theta_min, E_min : float, float, float
        """
        d = float(d_init)
        theta = float(theta_init)
        E = self.coupling_energy(d, theta)
        for _ in range(max_iter):
            # Central differences
            E_dp = self.coupling_energy(d + step_d * 0.5, theta)
            E_dm = self.coupling_energy(d - step_d * 0.5, theta)
            E_tp = self.coupling_energy(d, theta + step_theta * 0.5)
            E_tm = self.coupling_energy(d, theta - step_theta * 0.5)
            grad_d = (E_dp - E_dm) / step_d
            grad_t = (E_tp - E_tm) / step_theta

            # Adaptive step
            d_new = d - 0.05 * grad_d
            t_new = theta - 0.05 * grad_t
            d_new = max(d_new, 0.4)         # don't let cells overlap to nothing
            E_new = self.coupling_energy(d_new, t_new)
            if abs(E_new - E) < tol:
                d, theta, E = d_new, t_new, E_new
                break
            d, theta, E = d_new, t_new, E_new
        return d, theta, E

    # ---------- scans ---------------------------------------------------
    def binding_energy_curve(self,
                             d_range: Optional[NDArray] = None,
                             theta: float = 0.0
                             ) -> Tuple[NDArray, NDArray]:
        """Scan U(d) at fixed θ.  Returns (d_array, energy_array)."""
        if d_range is None:
            d_range = np.linspace(0.6, 4.0, 80)
        E = np.array([self.coupling_energy(float(d), theta) for d in d_range])
        return d_range, E

    def torsion_curve(self,
                       d: float = 1.4,
                       theta_range: Optional[NDArray] = None
                       ) -> Tuple[NDArray, NDArray]:
        """Scan U(θ) at fixed d.  Returns (theta_array, energy_array)."""
        if theta_range is None:
            theta_range = np.linspace(-np.pi, np.pi, 100)
        E = np.array([self.coupling_energy(d, float(t)) for t in theta_range])
        return theta_range, E

    # ---------- verification --------------------------------------------
    def verify_deuteron_geometry(self,
                                  tol_be_mev: float = 1e-3,
                                  tol_d_fm: float = 0.2
                                  ) -> dict:
        """Confirm the matched pair gives ε_face at d_min ≈ 1.4 fm.

        Returns a dict suitable for printing or asserting.
        """
        d_min, theta_min, E_min = self.relax_to_minimum()
        be = -E_min
        match = self.find_face_pair_match()
        d_target = self.xi
        return {
            'd_min_fm': d_min,
            'theta_min_rad': theta_min,
            'binding_energy_MeV': be,
            'eps_face_target_MeV': self.eps_face,
            'd_target_fm': d_target,
            'be_within_tol': abs(be - self.eps_face) < tol_be_mev,
            'd_within_tol': abs(d_min - d_target) < tol_d_fm,
            'matched_face_a': match.face_a,
            'matched_face_b': match.face_b,
            'matched_normal_dot': match.normal_dot,
            'matched_centroid_distance_fm': match.centroid_distance,
        }

    # ---------- visualization ------------------------------------------
    def make_visualization(self,
                           ax=None,
                           highlight_match: bool = True,
                           save_path: Optional[str] = None):
        """Render two K_4 cells in 3D with the matched face highlighted.

        Requires matplotlib.  Returns the (figure, axis) pair.  If
        ``save_path`` is given, the figure is saved (PDF or PNG).
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D       # noqa: F401
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        if ax is None:
            fig = plt.figure(figsize=(9, 7))
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        match = self.find_face_pair_match() if highlight_match else None

        for label, vertices, faces, color, alpha, cell_idx in [
            ('A', self.world_vertices_a(), self.cell_a.faces, 'tab:blue', 0.25, 'a'),
            ('B', self.world_vertices_b(), self.cell_b.faces, 'tab:orange', 0.25, 'b'),
        ]:
            # Plot vertices
            ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2],
                       c=color, s=40, label=f'cell {label}')
            # Plot edges
            for f in faces:
                tri = vertices[f]
                edges = [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]]
                for e in edges:
                    e_arr = np.array(e)
                    ax.plot(e_arr[:, 0], e_arr[:, 1], e_arr[:, 2], color=color, lw=1)
            # Plot all faces (lightly)
            face_polys = [vertices[f] for f in faces]
            ax.add_collection3d(Poly3DCollection(
                face_polys, facecolor=color, edgecolor=color,
                alpha=alpha))

        # Highlight matched face pair
        if match is not None:
            va = self.world_vertices_a()
            vb = self.world_vertices_b()
            face_a = self.cell_a.faces[match.face_a]
            face_b = self.cell_b.faces[match.face_b]
            ax.add_collection3d(Poly3DCollection(
                [va[face_a]], facecolor='red', edgecolor='red', alpha=0.6))
            ax.add_collection3d(Poly3DCollection(
                [vb[face_b]], facecolor='red', edgecolor='red', alpha=0.6))
            # Draw the contact line between centroids
            c_a, _ = self.world_faces_a()
            c_b, _ = self.world_faces_b()
            line = np.array([c_a[match.face_a], c_b[match.face_b]])
            ax.plot(line[:, 0], line[:, 1], line[:, 2], 'r--', lw=1.5)

        # Centroid markers
        ax.scatter(*self.pos_a, c='black', s=20, marker='x')
        ax.scatter(*self.pos_b, c='black', s=20, marker='x')

        ax.set_xlabel('x [fm]')
        ax.set_ylabel('y [fm]')
        ax.set_zlabel('z [fm]')
        ax.set_title('K_4 face-pair coupling (matched face in red)')
        ax.legend()

        if save_path is not None:
            fig.tight_layout()
            fig.savefig(save_path, dpi=150)
        return fig, ax


# ---------------------------------------------------------------------------
# Convenience: end-to-end demo
# ---------------------------------------------------------------------------

def demo() -> dict:
    """Build a K4Geometry, place two cells, scan U(d), report binding."""
    geom = K4Geometry()
    print(f"K_4: V={geom.n_vertices} E={geom.n_edges} F={geom.n_faces}")
    print(f"vertex angle: {geom.vertex_angle():.4f}° (expected {TETRAHEDRAL_ANGLE_DEG:.4f})")
    print(f"normal-pair angle: {geom.normal_pair_angle():.4f}° (expected {TETRAHEDRAL_ANGLE_DEG:.4f})")
    print(f"dihedral angle: {geom.dihedral_angle():.4f}° (expected {DIHEDRAL_ANGLE_DEG:.4f})")
    print(f"normals sum: {geom.normals_sum()} (expected ~0)")

    coupling = K4FacePairCoupling(cell_a=geom, cell_b=K4Geometry())
    coupling.place_two_cells(d_centers=1.4)
    match = coupling.find_face_pair_match()
    print(f"Matched faces: A.{match.face_a} ↔ B.{match.face_b}, "
          f"normal·= {match.normal_dot:.4f}, "
          f"centroid distance = {match.centroid_distance:.4f} fm")

    d_arr, E_arr = coupling.binding_energy_curve()
    print(f"E(d) min: {E_arr.min():.4f} MeV at d = {d_arr[np.argmin(E_arr)]:.4f} fm")

    report = coupling.verify_deuteron_geometry()
    print(f"verify_deuteron_geometry: {report}")
    return report


if __name__ == "__main__":
    demo()
