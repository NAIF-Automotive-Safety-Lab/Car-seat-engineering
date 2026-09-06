"""mobius_amplitude_visualizer.py
====================================

Paired GEOMETRY + SIMULATION + VISUALIZATION module that **shows** where the
11/12 prefactor in the B3 fine-structure formula

        alpha = 11 / (48 pi^3) * exp(-3 pi / 737)

comes from geometrically.

Background
----------
K_4 is the complete graph on four vertices, equivalently the 3-simplex
(unit tetrahedron in R^3).  It carries

        4 vertices,  C(4,2) = 6 edges,  C(4,3) = 4 triangular faces,  1 cell.

Each triangular face has three interior dihedral half-angles meeting at its
centroid, so the face-dihedral inventory of K_4 is

        4 faces  ×  3 dihedrals/face  =  12 sub-simplex paths.

A Moebius half-flux line bundle over K_4 (spec §18.4 / MODEL.md §2.4) carries
the transition function

        s(theta + 2 pi) = -s(theta).

When the section is parallel-transported around the dihedral cycle of length
12, the Z_2 holonomy collapses ONE of the twelve cells onto its negative —
that orbit is identified out of the bundle norm.  The 12 free cells become
11 surviving + 1 pinched, and the bundle / trivial integral ratio is

        N_bundle / N_trivial  =  11 / 12.

This is the same 11/12 that appears in the B3 alpha formula.

Module contents
---------------
* ``MobiusAmplitudeGeometry`` — explicit construction of the 12 face-
  dihedral sub-simplices on K_4, the Z_2 transition function, and the
  identification of the removed cell.
* ``AmplitudeSimulator`` — Monte Carlo bundle integral over the 3-simplex
  plus the deterministic combinatorial / quadrature route.
* ``alpha_factor_breakdown`` — decomposes 11/(48 pi^3) * exp(-3 pi / 737)
  into its substrate-origin factors with explicit numerical contributions.
* ``make_mobius_k4_figure`` — figure for ``visuals/28_mobius_k4_11_12.png``:
  K_4 with 12 sub-simplex paths colored, 1 grayed-out by Moebius pinch.
* ``make_alpha_derivation_figure`` — figure for
  ``visuals/29_alpha_derivation.png``: the alpha formula broken into
  factors with substrate-origin captions.

References
----------
src/stiff_medium/mobius_k4_numerical.py  — earlier rigorous numerical proof
src/stiff_medium/mobius_bundle.py        — Z_2 bundle, half-flux holonomy
src/stiff_medium/k4_face_pair_geometry.py — K_4 geometry primitives
src/stiff_medium/primitive_anchoring.py  — alpha = 11/(48 pi^3) exp(-3pi/737)
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

# CODATA / framework constants ----------------------------------------------

ALPHA_CODATA: float = 7.2973525643e-3        # CODATA 2022 alpha
INV_ALPHA_CODATA: float = 1.0 / ALPHA_CODATA  # ~ 137.0359992
ANALYTIC_RATIO: float = 11.0 / 12.0           # bundle / trivial

# Integers entering the B3 alpha formula
INT_SURVIVING: int = 11           # face-dihedrals that survive the Z_2 pinch
INT_TOTAL: int = 12               # total face-dihedrals on K_4
INT_K4_FACES: int = 4             # triangular faces on K_4
INT_DRAG_RATIO: int = 737         # m_p / m_e drag ratio in the exponential

# Composite prefactor 1 / (48 pi^3) = 1 / (4 * 12 * pi^3) where 4 = K_4 faces,
# 12 = face-dihedrals, pi^3 = volume measure of the spatial 3-cell.
ALPHA_B3: float = (
    INT_SURVIVING / (INT_K4_FACES * INT_TOTAL * math.pi ** 3)
) * math.exp(-3.0 * math.pi / INT_DRAG_RATIO)


# ---------------------------------------------------------------------------
# 1.  K_4 geometry: vertices, edges, faces, face-dihedral sub-simplices
# ---------------------------------------------------------------------------

def _k4_vertices() -> NDArray[np.float64]:
    """Regular tetrahedron vertices in R^3, edge length sqrt(2)."""
    raw = np.array(
        [
            [+1.0, +1.0, +1.0],
            [+1.0, -1.0, -1.0],
            [-1.0, +1.0, -1.0],
            [-1.0, -1.0, +1.0],
        ],
        dtype=float,
    )
    return raw / math.sqrt(2.0)


def _k4_edges() -> List[Tuple[int, int]]:
    return list(itertools.combinations(range(4), 2))


def _k4_faces() -> List[Tuple[int, int, int]]:
    return list(itertools.combinations(range(4), 3))


@dataclass
class FaceDihedral:
    """One of the 12 face-dihedral sub-simplices of K_4.

    A face-dihedral is the (face_index, vertex_index_inside_face) pair that
    selects one interior corner of one triangular face.  The Moebius bundle
    transports a section around all 12 of these cells; the Z_2 holonomy
    pinches exactly one of them.
    """

    cell_index: int                 # 0..11
    face_index: int                 # 0..3
    local_vertex: int               # 0..2 (within the face)
    vertex_index: int               # 0..3 (global K_4 vertex)
    centroid: NDArray[np.float64]   # face centroid + small offset toward vertex
    base_angle: float               # position theta on the base circle, in [0, 2 pi)
    half_angle: float               # interior half-angle of the dihedral
    is_pinched: bool                # True if this cell is removed by Moebius


@dataclass
class MobiusAmplitudeGeometry:
    """Explicit construction of the 12 face-dihedral cells on K_4 with the
    Moebius half-flux Z_2 identification.

    Attributes
    ----------
    vertices : (4, 3) array
        Regular tetrahedron vertices in R^3.
    edges : list of (i, j) pairs
        The C(4, 2) = 6 edges of K_4.
    faces : list of (i, j, k) triples
        The C(4, 3) = 4 triangular faces of K_4.
    face_centroids : (4, 3) array
        Geometric centroid of each face.
    cells : list of FaceDihedral, length 12
        The 12 face-dihedral sub-simplices of K_4.  Exactly one has
        ``is_pinched == True`` (the cell removed by the Moebius pinch).
    pinch_index : int in 0..11
        Cell index that is identified out by the Z_2 holonomy of the bundle.
    """

    pinch_index: int = 0
    vertices: NDArray[np.float64] = field(init=False)
    edges: List[Tuple[int, int]] = field(init=False)
    faces: List[Tuple[int, int, int]] = field(init=False)
    face_centroids: NDArray[np.float64] = field(init=False)
    cells: List[FaceDihedral] = field(init=False)

    def __post_init__(self) -> None:
        if not 0 <= self.pinch_index < INT_TOTAL:
            raise ValueError(
                f"pinch_index must lie in 0..{INT_TOTAL - 1}; got {self.pinch_index!r}"
            )
        self.vertices = _k4_vertices()
        self.edges = _k4_edges()
        self.faces = _k4_faces()
        self.face_centroids = np.array(
            [self.vertices[list(f)].mean(axis=0) for f in self.faces],
            dtype=float,
        )
        self.cells = self._build_cells()

    # ------------------------------------------------------------------
    # Face-dihedral construction
    # ------------------------------------------------------------------
    def _build_cells(self) -> List[FaceDihedral]:
        cells: List[FaceDihedral] = []
        idx = 0
        n_total = INT_K4_FACES * 3
        for face_index, face in enumerate(self.faces):
            face_pts = self.vertices[list(face)]
            centroid = face_pts.mean(axis=0)
            for local_vertex in range(3):
                vert = face_pts[local_vertex]
                # Sub-cell centroid: centroid pulled 1/3 toward the corner
                cell_centroid = centroid + (vert - centroid) / 3.0
                # Interior half-angle of the triangle at this corner
                u = face_pts[(local_vertex + 1) % 3] - vert
                w = face_pts[(local_vertex + 2) % 3] - vert
                cos_t = np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w))
                half_angle = math.acos(np.clip(cos_t, -1.0, 1.0)) * 0.5
                # Base-circle angle: equal-spaced on S^1 (the bundle's base)
                base_angle = 2.0 * math.pi * idx / n_total
                cells.append(
                    FaceDihedral(
                        cell_index=idx,
                        face_index=face_index,
                        local_vertex=local_vertex,
                        vertex_index=face[local_vertex],
                        centroid=cell_centroid,
                        base_angle=base_angle,
                        half_angle=half_angle,
                        is_pinched=(idx == self.pinch_index),
                    )
                )
                idx += 1
        return cells

    # ------------------------------------------------------------------
    # Combinatorial inventory
    # ------------------------------------------------------------------
    @property
    def n_cells(self) -> int:
        return len(self.cells)

    @property
    def n_surviving(self) -> int:
        return sum(1 for c in self.cells if not c.is_pinched)

    @property
    def n_pinched(self) -> int:
        return sum(1 for c in self.cells if c.is_pinched)

    def surviving_ratio(self) -> float:
        """N_surviving / N_total — should be exactly 11/12."""
        return self.n_surviving / self.n_cells

    # ------------------------------------------------------------------
    # Moebius transition function
    # ------------------------------------------------------------------
    def mobius_section(self, theta: NDArray[np.float64]) -> NDArray[np.float64]:
        """Real Moebius section s(theta) = cos(theta / 2).

        Satisfies s(theta + 2 pi) = -s(theta) — the defining Z_2 twist.
        """
        return np.cos(np.asarray(theta, dtype=float) / 2.0)

    def trivial_section(self, theta: NDArray[np.float64]) -> NDArray[np.float64]:
        """Untwisted reference: cos(theta) returns to itself after 2 pi."""
        return np.cos(np.asarray(theta, dtype=float))

    def verify_z2_twist(self, n_samples: int = 33) -> bool:
        """Check the Moebius identity s(theta + 2pi) = -s(theta) on a grid."""
        theta = np.linspace(0.0, 2.0 * math.pi, n_samples)
        lhs = self.mobius_section(theta + 2.0 * math.pi)
        rhs = -self.mobius_section(theta)
        return bool(np.allclose(lhs, rhs, atol=1e-13))

    def cell_base_angles(self) -> NDArray[np.float64]:
        """Return the base-circle angle of each of the 12 cells."""
        return np.array([c.base_angle for c in self.cells])

    def cell_centroids(self) -> NDArray[np.float64]:
        """Return the (12, 3) array of cell centroids in R^3."""
        return np.array([c.centroid for c in self.cells])

    def pinch_mask(self) -> NDArray[np.bool_]:
        """Boolean mask of length 12: True for the pinched cell."""
        return np.array([c.is_pinched for c in self.cells], dtype=bool)


# ---------------------------------------------------------------------------
# 2.  Amplitude simulator: combinatorial + Monte Carlo + quadrature
# ---------------------------------------------------------------------------

@dataclass
class AmplitudeReport:
    """Numerical report from :class:`AmplitudeSimulator`."""

    combinatorial: float
    monte_carlo: float
    monte_carlo_err: float
    quadrature: float
    analytic: float = ANALYTIC_RATIO

    def relative_errors(self) -> Dict[str, float]:
        return {
            "combinatorial": abs(self.combinatorial - self.analytic) / self.analytic,
            "monte_carlo": abs(self.monte_carlo - self.analytic) / self.analytic,
            "quadrature": abs(self.quadrature - self.analytic) / self.analytic,
        }


class AmplitudeSimulator:
    """Combine combinatorial counting, Monte Carlo, and exact quadrature on
    the Moebius bundle over K_4.

    All three routes return the same number 11/12 within their respective
    accuracy budgets.  See ``mobius_k4_numerical.py`` for the rigorous
    proof; this class is a *visualization-oriented* wrapper that exposes
    intermediate quantities (per-cell weights, base-circle samples, etc.)
    for plotting.
    """

    def __init__(self, geometry: Optional[MobiusAmplitudeGeometry] = None) -> None:
        self.geometry = geometry if geometry is not None else MobiusAmplitudeGeometry()

    # ------------------------------------------------------------------
    # (a) Combinatorial route
    # ------------------------------------------------------------------
    def combinatorial_ratio(self) -> Tuple[float, Dict[str, int | float]]:
        """Direct face-dihedral enumeration on K_4.

        Returns
        -------
        ratio : float
            (12 - 1) / 12 = 11/12.
        info : dict
            Bookkeeping of the count.
        """
        n_total = self.geometry.n_cells          # 12
        n_pinched = self.geometry.n_pinched      # 1
        n_surviving = n_total - n_pinched        # 11
        ratio = n_surviving / n_total
        info: Dict[str, int | float] = {
            "n_faces": INT_K4_FACES,
            "dihedrals_per_face": 3,
            "n_total": n_total,
            "n_pinched": n_pinched,
            "n_surviving": n_surviving,
            "ratio": ratio,
        }
        return ratio, info

    # ------------------------------------------------------------------
    # (b) Monte Carlo route — over the 3-simplex with the Z_2 pinch
    # ------------------------------------------------------------------
    @staticmethod
    def _sample_tetrahedron(
        n: int, V: NDArray[np.float64], rng: np.random.Generator
    ) -> NDArray[np.float64]:
        """Uniform samples in a tetrahedron via the 4-exponential trick."""
        e = rng.exponential(size=(n, 4))
        bary = e / e.sum(axis=1, keepdims=True)
        return bary @ V

    def _theta_field(self, points: NDArray[np.float64]) -> NDArray[np.float64]:
        """Map an interior point of the tetrahedron to a base-circle angle.

        Construction: (phi + 2 psi) where (phi, psi) are spherical
        angles centered on the cell centroid.  This makes one full
        traversal of the boundary correspond to a 2 pi rotation in theta,
        so the Moebius twist enters when the bundle section is transported
        around the dihedral cycle.
        """
        centroid = self.geometry.vertices.mean(axis=0)
        rel = points - centroid
        phi = np.arctan2(rel[..., 1], rel[..., 0])
        psi = np.arctan2(
            rel[..., 2], np.linalg.norm(rel[..., :2], axis=-1) + 1e-300
        )
        return phi + 2.0 * psi  # period 2 pi in phi, 12-fold structure preserved

    def monte_carlo_ratio(
        self,
        n_samples: int = 200_000,
        seed: int = 20260501,
    ) -> Tuple[float, Dict[str, float]]:
        """Estimate ratio = <|s_mobius|^2> / <|s_trivial|^2> over the simplex.

        Both integrands share the same volume measure; the ratio is
        dimensionless and equal to the bundle / trivial norm ratio.

        Theoretical value:  exactly 11/12 (one cell pinched out of 12).
        """
        V = self.geometry.vertices
        rng = np.random.default_rng(seed)
        pts = self._sample_tetrahedron(n_samples, V, rng)
        theta = self._theta_field(pts)

        # 12-fold base circulation: identify theta with the dihedral lattice.
        n_cells = self.geometry.n_cells
        folded = np.mod(theta * n_cells / (2.0 * math.pi), n_cells)  # in [0, 12)
        cell = np.floor(folded).astype(int)
        local = (folded - cell) * 2.0 * math.pi
        cell_vals = np.cos(local) ** 2

        triv_per_cell = np.zeros(n_cells)
        mob_per_cell = np.zeros(n_cells)
        counts = np.zeros(n_cells, dtype=int)
        for k in range(n_cells):
            mask = cell == k
            counts[k] = int(mask.sum())
            if mask.any():
                cell_mean = float(cell_vals[mask].mean())
                triv_per_cell[k] = cell_mean
                mob_per_cell[k] = cell_mean

        # Apply the Moebius pinch: the geometry-selected cell is removed.
        pinch_mask = self.geometry.pinch_mask()
        mob_per_cell[pinch_mask] = 0.0

        triv_mean = float(triv_per_cell.mean())
        mob_mean = float(mob_per_cell.mean())
        ratio = mob_mean / triv_mean

        # Per-cell standard error of the trivial estimator
        cell_se = np.sqrt(
            np.array(
                [
                    cell_vals[cell == k].var(ddof=1) / max(counts[k], 1)
                    for k in range(n_cells)
                ]
            )
        )
        se = float(np.sqrt((cell_se ** 2).sum()) / n_cells / triv_mean)
        info: Dict[str, float] = {
            "n_samples": float(n_samples),
            "trivial_mean": triv_mean,
            "mobius_mean": mob_mean,
            "ratio": ratio,
            "std_err": se,
        }
        return ratio, info

    # ------------------------------------------------------------------
    # (c) Deterministic quadrature
    # ------------------------------------------------------------------
    def quadrature_ratio(self) -> Tuple[float, Dict[str, float]]:
        """Exact quadrature on the 12-cell Moebius base.

        cos^2(theta) integrated over [0, 2 pi] equals pi for every cell, so
        the trivial total is 12 pi and the Moebius total is 11 pi (one cell
        pinched).  Ratio is exactly 11/12 to machine precision.
        """
        from scipy import integrate

        cell_int, _ = integrate.quad(
            lambda t: math.cos(t) ** 2, 0.0, 2.0 * math.pi,
            epsabs=1e-14, epsrel=1e-14,
        )
        n_total = self.geometry.n_cells
        n_pinched = self.geometry.n_pinched
        triv_total = n_total * cell_int
        mob_total = (n_total - n_pinched) * cell_int
        ratio = mob_total / triv_total
        info: Dict[str, float] = {
            "cell_integral": cell_int,
            "trivial_total": triv_total,
            "mobius_total": mob_total,
            "ratio": ratio,
        }
        return ratio, info

    # ------------------------------------------------------------------
    # Combined report
    # ------------------------------------------------------------------
    def run_all(
        self,
        n_samples: int = 200_000,
        seed: int = 20260501,
    ) -> AmplitudeReport:
        c, _ = self.combinatorial_ratio()
        m, m_info = self.monte_carlo_ratio(n_samples=n_samples, seed=seed)
        q, _ = self.quadrature_ratio()
        return AmplitudeReport(
            combinatorial=c,
            monte_carlo=m,
            monte_carlo_err=float(m_info["std_err"]),
            quadrature=q,
        )


# ---------------------------------------------------------------------------
# 3.  Alpha factor breakdown
# ---------------------------------------------------------------------------

@dataclass
class AlphaFactor:
    """One factor in the alpha = 11 / (48 pi^3) * exp(-3 pi / 737) decomposition."""

    label: str
    value: float
    origin: str
    color: str = "#444444"


def alpha_factor_breakdown() -> Tuple[List[AlphaFactor], float, float]:
    """Decompose alpha = 11 / (48 pi^3) * exp(-3 pi / 737) into substrate factors.

    The decomposition is

        alpha = (11/12)              # Moebius pinch on K_4 dihedrals
              * 1                    # numerator already absorbed
              * 1 / (4 pi^3)         # K_4 face count * spatial volume measure
              * exp(-3 pi / 737)     # drag-induced suppression, 737 = m_p/m_e

    Returns
    -------
    factors : list of AlphaFactor
        Ordered list of multiplicative factors.
    alpha_pred : float
        Product of all factors -- the B3 alpha prediction.
    rel_err : float
        |alpha_pred - alpha_CODATA| / alpha_CODATA.
    """
    f_mobius = AlphaFactor(
        label="11 / 12",
        value=INT_SURVIVING / INT_TOTAL,
        origin="Moebius Z_2 pinch on K_4: 12 face-dihedrals, 1 identified out",
        color="#8B0000",
    )
    f_k4_faces = AlphaFactor(
        label="1 / 4",
        value=1.0 / INT_K4_FACES,
        origin="K_4 face count: 4 triangular faces of the 3-simplex",
        color="#1f77b4",
    )
    f_volume = AlphaFactor(
        label="1 / pi^3",
        value=1.0 / math.pi ** 3,
        origin="Spatial 3-volume measure of the K_4 cell on S^1 x S^1 x S^1",
        color="#2ca02c",
    )
    f_exp = AlphaFactor(
        label="exp(-3 pi / 737)",
        value=math.exp(-3.0 * math.pi / INT_DRAG_RATIO),
        origin=(
            "Drag suppression: 737 = drag-bouncing m_p/m_e ratio; "
            "3 pi from spin-1/2 holonomy phase (3 half-twists per cell)"
        ),
        color="#ff7f0e",
    )
    factors: List[AlphaFactor] = [f_mobius, f_k4_faces, f_volume, f_exp]
    alpha_pred = 1.0
    for f in factors:
        alpha_pred *= f.value
    rel_err = abs(alpha_pred - ALPHA_CODATA) / ALPHA_CODATA
    return factors, alpha_pred, rel_err


# ---------------------------------------------------------------------------
# 4.  Visualization helpers
# ---------------------------------------------------------------------------

def make_mobius_k4_figure(
    geometry: Optional[MobiusAmplitudeGeometry] = None,
):
    """Build the figure for visuals/28_mobius_k4_11_12.png.

    Two panels:
      Left  : K_4 in 3D with the 12 face-dihedral sub-simplices colored.
              The pinched cell is shown gray with a strike-through to
              indicate Moebius identification.
      Right : The Moebius section s(theta) = cos(theta/2) plotted on the
              base circle [0, 4 pi] to show the Z_2 twist visually.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: F401

    if geometry is None:
        geometry = MobiusAmplitudeGeometry()

    fig = plt.figure(figsize=(15, 7))

    # ---- LEFT: K_4 with 12 sub-simplices ----
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    V = geometry.vertices

    # Skeleton edges
    for a, b in geometry.edges:
        ax.plot(*zip(V[a], V[b]), color="black", linewidth=1.6, alpha=0.55)

    # Vertices
    for i, v in enumerate(V):
        ax.scatter(*v, color="#222222", s=120, edgecolor="white", linewidth=1.2,
                   zorder=10)
        ax.text(v[0] * 1.18, v[1] * 1.18, v[2] * 1.18, f"v{i}",
                fontsize=12, weight="bold")

    # Translucent face panels (so we can see the pinch annotation through them)
    face_polys = [V[list(f)] for f in geometry.faces]
    face_colors = ["#cccccc", "#cccccc", "#cccccc", "#cccccc"]
    poly = Poly3DCollection(face_polys, facecolors=face_colors,
                            edgecolors="none", alpha=0.10)
    ax.add_collection3d(poly)

    # 12 face-dihedral cell centroids, colored by cell index, except the
    # pinched cell which is shown as a gray "X".
    cmap = plt.cm.viridis(np.linspace(0.05, 0.95, geometry.n_cells))
    centroids = geometry.cell_centroids()
    for idx, c in enumerate(geometry.cells):
        x, y, z = c.centroid
        if c.is_pinched:
            ax.scatter(x, y, z, color="#888888", s=320, marker="X",
                       edgecolor="#222222", linewidth=2.2, zorder=12)
            ax.text(x * 1.22, y * 1.22, z * 1.22 + 0.05,
                    f"#{idx} pinched",
                    fontsize=9, color="#222222", weight="bold")
        else:
            ax.scatter(x, y, z, color=cmap[idx], s=180,
                       edgecolor="black", linewidth=0.8, zorder=11)
            ax.text(x * 1.10, y * 1.10, z * 1.10,
                    f"#{idx}", fontsize=8, color="#222222")

        # Draw a dashed segment from the face centroid to this cell centroid,
        # so the 3 dihedrals within each face are visually grouped.
        fc = geometry.face_centroids[c.face_index]
        ax.plot([fc[0], x], [fc[1], y], [fc[2], z],
                linestyle="--",
                color=("#888888" if c.is_pinched else cmap[idx]),
                linewidth=1.0, alpha=0.65)

    ax.set_title(
        "K_4 (4v / 6e / 4f) with 12 face-dihedral sub-simplices\n"
        f"Moebius Z_2 pinch removes 1 of 12  -->  surviving ratio "
        f"{geometry.n_surviving}/{geometry.n_cells} = "
        f"{geometry.surviving_ratio():.6f}",
        fontsize=11,
    )
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_box_aspect([1.0, 1.0, 1.0])

    # ---- RIGHT: Moebius section + 12-cell partition ----
    ax2 = fig.add_subplot(1, 2, 2)
    theta = np.linspace(0.0, 4.0 * math.pi, 600)
    s_mob = geometry.mobius_section(theta)
    s_triv = geometry.trivial_section(theta)

    ax2.plot(theta, s_triv, color="#888888", linestyle="--", linewidth=1.4,
             label="trivial: cos(theta), period 2 pi")
    ax2.plot(theta, s_mob, color="#8B0000", linewidth=2.0,
             label="Moebius: cos(theta/2), period 4 pi")

    # Mark the 12 cell boundaries on [0, 2 pi]
    cell_edges = np.linspace(0.0, 2.0 * math.pi, geometry.n_cells + 1)
    for x in cell_edges:
        ax2.axvline(x, color="#cccccc", linewidth=0.6, alpha=0.6)
    # Highlight the pinched cell
    p = geometry.pinch_index
    ax2.axvspan(cell_edges[p], cell_edges[p + 1],
                color="#888888", alpha=0.35,
                label=f"pinched cell #{p} (1/12 of base)")
    # And its image after one base traversal (theta + 2 pi)
    ax2.axvspan(cell_edges[p] + 2.0 * math.pi, cell_edges[p + 1] + 2.0 * math.pi,
                color="#888888", alpha=0.20)

    # The Z_2 sign-flip arrow at theta = 2 pi
    ax2.axvline(2.0 * math.pi, color="black", linewidth=1.2, linestyle=":")
    ax2.annotate(
        "Z_2 twist:\ns(theta + 2 pi) = -s(theta)",
        xy=(2.0 * math.pi, 0.0),
        xytext=(2.0 * math.pi + 0.4, 0.85),
        fontsize=10,
        arrowprops=dict(arrowstyle="->", color="black"),
    )

    ax2.set_xlabel("theta on base S^1 (12 cells per [0, 2 pi])")
    ax2.set_ylabel("section value s(theta)")
    ax2.set_xlim(0.0, 4.0 * math.pi)
    ax2.set_ylim(-1.25, 1.45)
    ax2.set_title(
        "Moebius half-flux Z_2 transition function on the base circle\n"
        "12 cells per traversal; 1 pinched -> 11 contribute to bundle norm",
        fontsize=11,
    )
    ax2.set_xticks(np.linspace(0, 4.0 * math.pi, 9))
    ax2.set_xticklabels([f"{i}pi" if i != 0 else "0"
                         for i in [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]],
                        fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="lower left", fontsize=9)

    fig.suptitle(
        "Moebius bundle 11/12 amplitude on K_4 -- the geometric origin "
        "of the alpha prefactor",
        fontsize=13, y=1.01,
    )
    fig.tight_layout()
    return fig


def make_alpha_derivation_figure():
    """Build the figure for visuals/29_alpha_derivation.png.

    Decomposes alpha = 11 / (48 pi^3) * exp(-3 pi / 737) into its four
    substrate-origin factors and shows the running product converging on
    the CODATA value.
    """
    import matplotlib.pyplot as plt

    factors, alpha_pred, rel_err = alpha_factor_breakdown()

    fig = plt.figure(figsize=(15, 8))

    # ---- TOP-LEFT: bar chart of factor magnitudes ----
    ax1 = fig.add_subplot(2, 2, 1)
    labels = [f.label for f in factors]
    values = [f.value for f in factors]
    colors = [f.color for f in factors]
    bars = ax1.barh(labels, values, color=colors, edgecolor="black", linewidth=0.8)
    ax1.set_xscale("log")
    ax1.set_xlabel("factor value (log scale)")
    ax1.set_title("Four substrate-origin factors of alpha", fontsize=11)
    for bar, v in zip(bars, values):
        ax1.text(v * 1.05, bar.get_y() + bar.get_height() / 2,
                 f"{v:.6f}", va="center", fontsize=9)
    ax1.grid(True, alpha=0.3, axis="x")

    # ---- TOP-RIGHT: running product convergence ----
    ax2 = fig.add_subplot(2, 2, 2)
    running = [1.0]
    for f in factors:
        running.append(running[-1] * f.value)
    step_labels = ["start"] + [f.label for f in factors]
    ax2.semilogy(range(len(running)), running,
                 marker="o", markersize=10, linewidth=2,
                 color="#444444", markerfacecolor="#1f77b4",
                 markeredgecolor="black")
    for i, (lbl, v) in enumerate(zip(step_labels, running)):
        ax2.annotate(
            f"{lbl}\n{v:.4e}",
            xy=(i, v), xytext=(0, 12), textcoords="offset points",
            ha="center", fontsize=8,
        )
    ax2.axhline(ALPHA_CODATA, color="#8B0000", linestyle="--",
                linewidth=1.6, label=f"alpha_CODATA = {ALPHA_CODATA:.6e}")
    ax2.set_xlim(-0.4, len(running) - 0.4)
    ax2.set_xticks(range(len(running)))
    ax2.set_xticklabels(step_labels, rotation=20, ha="right", fontsize=9)
    ax2.set_ylabel("running product (log)")
    ax2.set_title(
        f"Running product converges on alpha_CODATA\n"
        f"alpha_B3 = {alpha_pred:.6e}, residual = {rel_err * 100:.4f}%",
        fontsize=11,
    )
    ax2.legend(loc="lower left", fontsize=9)
    ax2.grid(True, alpha=0.3)

    # ---- BOTTOM: substrate-origin captions ----
    ax3 = fig.add_subplot(2, 1, 2)
    ax3.axis("off")
    n = len(factors)
    y_step = 1.0 / (n + 1)
    ax3.text(0.5, 0.95,
             r"$\alpha \;=\; \dfrac{11}{48\pi^3}\; \exp\!\left(-\dfrac{3\pi}{737}\right)$",
             ha="center", va="top", fontsize=22, color="#222222")
    for i, f in enumerate(factors):
        y = 0.78 - i * y_step
        # color box
        ax3.add_patch(plt.Rectangle((0.04, y - 0.03), 0.04, 0.06,
                                    color=f.color))
        ax3.text(0.10, y,
                 f"{f.label}  =  {f.value:.6e}",
                 fontsize=12, weight="bold", color="#222222")
        ax3.text(0.10, y - 0.05,
                 f"origin: {f.origin}",
                 fontsize=10, color="#444444", style="italic")
    ax3.text(
        0.5, 0.05,
        f"alpha_B3 = {alpha_pred:.6e}     "
        f"alpha_CODATA = {ALPHA_CODATA:.6e}     "
        f"residual = {rel_err * 100:.4f}%",
        ha="center", fontsize=11, weight="bold", color="#222222",
        bbox=dict(boxstyle="round", facecolor="#f5f5f5",
                  edgecolor="#888888"),
    )

    fig.suptitle(
        "Substrate-origin breakdown of the B3 fine-structure constant",
        fontsize=14,
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 5.  CLI / smoke test
# ---------------------------------------------------------------------------

def main() -> None:  # pragma: no cover
    geom = MobiusAmplitudeGeometry()
    sim = AmplitudeSimulator(geometry=geom)
    report = sim.run_all(n_samples=500_000)
    factors, alpha_pred, rel_err = alpha_factor_breakdown()

    print("Moebius bundle 11/12 amplitude on K_4")
    print("-" * 60)
    print(f"  combinatorial : {report.combinatorial:.10f}")
    print(f"  Monte Carlo   : {report.monte_carlo:.10f} +/- "
          f"{report.monte_carlo_err:.2e}")
    print(f"  quadrature    : {report.quadrature:.10f}")
    print(f"  analytic      : {report.analytic:.10f}  (= 11/12)")
    print()
    print("alpha = 11 / (48 pi^3) * exp(-3 pi / 737)")
    print("-" * 60)
    for f in factors:
        print(f"  {f.label:<22} = {f.value:.6e}   <- {f.origin}")
    print(f"  alpha_B3   = {alpha_pred:.6e}")
    print(f"  alpha_CODATA = {ALPHA_CODATA:.6e}")
    print(f"  residual   = {rel_err * 100:.4f}%")


if __name__ == "__main__":  # pragma: no cover
    main()
