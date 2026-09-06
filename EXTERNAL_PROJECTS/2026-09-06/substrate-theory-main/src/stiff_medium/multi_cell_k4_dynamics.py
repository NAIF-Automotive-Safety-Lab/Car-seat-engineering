"""Multi-cell K_4 tetrahedron dynamics — face-pair coupling for nuclear binding.

This module simulates N tetrahedral (K_4) substrate cells in a 3D box,
interacting via face-pair coupling.  In the B3 framework a single K_4
cell is a nucleon; two K_4 sharing a face form a deuteron (BE = ε_face ≈
2.222 MeV); four K_4 in a tetrahedral arrangement form an alpha particle
(BE ≈ 28 MeV ≈ 6 ε_face × ~2 face-overlap correction → see results).

THE PHYSICS
-----------
Each K_4 cell is parametrised by:

    r_i ∈ ℝ³                position of the cell centre
    R_i ∈ SO(3)             orientation, given by Euler angles (α, β, γ)

A K_4 has four triangular faces with outward normals n̂_a^(i) (a = 1..4),
each face centred on the body-frame point f_a (centroid of the
opposite-vertex triangle, at distance d_face from the centre).

Face-pair coupling.  Two cells i, j contribute -ε_face per pair of
faces (a from i, b from j) when:

    1. The face centres are close:        ||(r_i + R_i f_a) - (r_j + R_j f_b)|| < ξ
    2. The face normals are anti-aligned: n̂_a^(i) · n̂_b^(j) < -cos(30°)

Smoothing.  We use a smooth Gaussian × cos² gating so the energy is
differentiable and Newtonian dynamics is well-defined:

    g_pair(d, θ) = exp(-(d/ξ)²) × max(0, -n̂_a · n̂_b)²

    E_pair = -ε_face × Σ_{a,b} g_pair(d_ab, θ_ab)

The total potential energy is U = Σ_{i<j} E_pair(i, j).

Time evolution.  Newtonian translation only (orientations are minimised
quasi-statically each step via gradient descent on U; rotational inertia
of the K_4 is small compared to substrate drag in the regime of
interest):

    m_cell × r̈_i = -∇_{r_i} U  -  γ ṙ_i

with m_cell = m_nucleon = 938.272 MeV/c² and substrate drag γ chosen
to overdamp the motion to equilibrium.  The Euler-angle update is a
simple gradient step on U with the same γ.

The binding energy of an N-cell configuration is

    BE(N) = -U_min(N)

evaluated at the relaxed equilibrium.

UNITS
-----
Energy: MeV.  Length: fm.  Time: fm/c (≈ 3.3 × 10⁻²⁴ s).
m_nucleon = 938.272 MeV/c²; in (MeV)(fm/c)²/fm² this is 938.272.

REFERENCES
----------
Spec §face-pair coupling, B3 deuteron derivation (deuteron_v2.py),
alpha tetrahedron (alpha_tetrahedron.py).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple, List, Optional

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Physical constants (B3 nuclear sector)
# ---------------------------------------------------------------------------

EPS_FACE_MEV: float = 2.222            # face-pair binding (deuteron BE)
M_NUCLEON_MEV: float = 938.272         # nucleon mass in MeV/c²
XI_FM: float = 1.4                     # face-coupling range (~ pion Compton λ̄)
D_FACE_FM: float = 0.8                 # body-frame face-centre distance from K_4 centre
ALIGN_THRESHOLD_DEG: float = 30.0      # anti-alignment tolerance


# ---------------------------------------------------------------------------
# K_4 geometry (body frame)
# ---------------------------------------------------------------------------

def k4_body_vertices() -> NDArray[np.float64]:
    """Four vertices of a regular tetrahedron centred at the origin.

    Edge length is normalised so that the face-centre distance to the
    centre is D_FACE_FM.  For a regular tetrahedron of edge length a:

        d_face (centre → face-centroid) = a / (2√6)

    so a = 2√6 × D_FACE_FM.
    """
    a = 2.0 * np.sqrt(6.0) * D_FACE_FM
    # standard regular-tetrahedron coords, scaled
    raw = np.array([
        [+1.0, +1.0, +1.0],
        [+1.0, -1.0, -1.0],
        [-1.0, +1.0, -1.0],
        [-1.0, -1.0, +1.0],
    ], dtype=float)
    # raw edge length is 2√2; scale so edge = a
    return raw * (a / (2.0 * np.sqrt(2.0)))


def k4_body_face_data() -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Body-frame face centroids and outward-normal unit vectors.

    Each face is the triangle opposite one vertex; its centroid is
    -v_a / 3 (mean of the other three vertices for a centroid-zero
    tetrahedron), and its outward normal points from the centre through
    the face centroid.

    Returns:
        (centres, normals): each shape (4, 3).
    """
    V = k4_body_vertices()
    centres = np.zeros((4, 3))
    normals = np.zeros((4, 3))
    for a in range(4):
        # face = three vertices excluding vertex a
        face_verts = np.delete(V, a, axis=0)
        c = face_verts.mean(axis=0)
        centres[a] = c
        n = c / np.linalg.norm(c)        # outward (away from origin)
        normals[a] = n
    return centres, normals


# ---------------------------------------------------------------------------
# Rotations (Euler angles, ZYZ convention)
# ---------------------------------------------------------------------------

def euler_to_matrix(alpha: float, beta: float, gamma: float
                    ) -> NDArray[np.float64]:
    """ZYZ Euler-angle rotation matrix R = Rz(α) Ry(β) Rz(γ)."""
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    Rz1 = np.array([[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]])
    Ry  = np.array([[cb, 0.0, sb], [0.0, 1.0, 0.0], [-sb, 0.0, cb]])
    Rz2 = np.array([[cg, -sg, 0.0], [sg, cg, 0.0], [0.0, 0.0, 1.0]])
    return Rz1 @ Ry @ Rz2


# ---------------------------------------------------------------------------
# Cell state
# ---------------------------------------------------------------------------

@dataclass
class K4Cell:
    """A single K_4 tetrahedron cell."""
    position: NDArray[np.float64]            # shape (3,), in fm
    euler: NDArray[np.float64]               # shape (3,) (α, β, γ) in radians
    velocity: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(3))

    def world_face_data(self) -> Tuple[NDArray[np.float64],
                                        NDArray[np.float64]]:
        """Return (centres, normals) for this cell's 4 faces in world coords."""
        body_c, body_n = k4_body_face_data()
        R = euler_to_matrix(*self.euler)
        centres = self.position + body_c @ R.T
        normals = body_n @ R.T
        return centres, normals


# ---------------------------------------------------------------------------
# Face-pair coupling energy
# ---------------------------------------------------------------------------

def _gate(d: NDArray[np.float64], dot: NDArray[np.float64],
          xi: float = XI_FM) -> NDArray[np.float64]:
    """Smooth gate: Gaussian in distance × clipped (-dot)² in alignment.

    d:   pairwise face-centre distances  (M,)
    dot: pairwise face-normal dot products (M,)

    Returns g ∈ [0, 1].
    """
    spatial = np.exp(-(d / xi) ** 2)
    # anti-alignment: want n_a · n_b ≈ -1.  Use clipped (-dot) ∈ [0,1] then square.
    cos_thresh = np.cos(np.deg2rad(ALIGN_THRESHOLD_DEG))
    align = np.clip(-dot - cos_thresh, 0.0, None) / max(1.0 - cos_thresh, 1e-9)
    return spatial * align ** 2


def pair_energy(cell_i: K4Cell, cell_j: K4Cell,
                eps: float = EPS_FACE_MEV,
                xi: float = XI_FM) -> float:
    """Face-pair coupling energy between two cells.

    Sums over all 16 face-pair combinations the smoothed contribution
    -ε × g_pair(d, θ).  In an idealized "shared face" geometry, exactly
    one pair contributes ≈ -ε; geometry constrains the rest.
    """
    c_i, n_i = cell_i.world_face_data()
    c_j, n_j = cell_j.world_face_data()
    # pairwise (4, 4)
    diff = c_i[:, None, :] - c_j[None, :, :]
    d = np.linalg.norm(diff, axis=-1)
    dot = n_i @ n_j.T
    g = _gate(d, dot, xi=xi)
    return float(-eps * g.sum())


def total_energy(cells: List[K4Cell], eps: float = EPS_FACE_MEV,
                 xi: float = XI_FM) -> float:
    """Total potential energy U = Σ_{i<j} E_pair."""
    U = 0.0
    n = len(cells)
    for i in range(n):
        for j in range(i + 1, n):
            U += pair_energy(cells[i], cells[j], eps=eps, xi=xi)
    return U


# ---------------------------------------------------------------------------
# Numerical gradients (finite-difference)
# ---------------------------------------------------------------------------

def _grad_position(cells: List[K4Cell], idx: int,
                   eps_fd: float = 1e-4) -> NDArray[np.float64]:
    """∂U/∂r_idx via central differences."""
    g = np.zeros(3)
    for k in range(3):
        cells[idx].position[k] += eps_fd
        Up = total_energy(cells)
        cells[idx].position[k] -= 2.0 * eps_fd
        Um = total_energy(cells)
        cells[idx].position[k] += eps_fd
        g[k] = (Up - Um) / (2.0 * eps_fd)
    return g


def _grad_euler(cells: List[K4Cell], idx: int,
                eps_fd: float = 1e-4) -> NDArray[np.float64]:
    """∂U/∂(α,β,γ)_idx via central differences."""
    g = np.zeros(3)
    for k in range(3):
        cells[idx].euler[k] += eps_fd
        Up = total_energy(cells)
        cells[idx].euler[k] -= 2.0 * eps_fd
        Um = total_energy(cells)
        cells[idx].euler[k] += eps_fd
        g[k] = (Up - Um) / (2.0 * eps_fd)
    return g


# ---------------------------------------------------------------------------
# Time evolution
# ---------------------------------------------------------------------------

@dataclass
class RelaxResult:
    """Output of a relaxation run."""
    cells: List[K4Cell]
    energy_history: NDArray[np.float64]
    final_energy: float
    binding_energy: float           # = -final_energy (ε_face units MeV)
    n_iter: int


def relax(cells: List[K4Cell],
          dt: float = 0.01,
          gamma: float = 200.0,
          n_iter: int = 4000,
          tol: float = 1e-7,
          orient_step: float = 0.01,
          m_cell: float = M_NUCLEON_MEV,
          eps: float = EPS_FACE_MEV,
          xi: float = XI_FM,
          ) -> RelaxResult:
    """Overdamped Newtonian relaxation to the energy minimum.

    Translational EOM:  m r̈ = -∇U − γ ṙ  (semi-implicit Euler)
    Orientational EOM:  α̇ = -orient_step × ∂U/∂α   (gradient flow)

    Args:
        cells:       Initial configuration (modified in place is OK).
        dt:          Time step in fm/c.
        gamma:       Drag coefficient (large → overdamped).
        n_iter:      Maximum iterations.
        tol:         Stop when |ΔU| < tol over 20 iterations.
        orient_step: Gradient-descent step for Euler angles.
        m_cell:      Effective mass per cell.
        eps:         Face-pair coupling strength.
        xi:          Coupling range.

    Returns:
        RelaxResult with the relaxed cells and energy history.
    """
    n = len(cells)
    history = []
    prev = total_energy(cells, eps=eps, xi=xi)
    history.append(prev)

    for it in range(1, n_iter + 1):
        # Compute all gradients first (so we update cells synchronously)
        grad_r = [np.zeros(3) for _ in range(n)]
        grad_e = [np.zeros(3) for _ in range(n)]
        for i in range(n):
            grad_r[i] = _grad_position(cells, i)
            grad_e[i] = _grad_euler(cells, i)

        for i in range(n):
            # semi-implicit damped Newton
            v = cells[i].velocity
            a = (-grad_r[i] - gamma * v) / m_cell
            cells[i].velocity = v + dt * a
            cells[i].position = cells[i].position + dt * cells[i].velocity
            # Euler-angle gradient flow
            cells[i].euler = cells[i].euler - orient_step * grad_e[i]

        E = total_energy(cells, eps=eps, xi=xi)
        history.append(E)

        # convergence: smoothed change
        if it > 30:
            window = np.array(history[-20:])
            if np.max(window) - np.min(window) < tol:
                break
        prev = E

    final_E = history[-1]
    return RelaxResult(
        cells=cells,
        energy_history=np.array(history),
        final_energy=final_E,
        binding_energy=-final_E,
        n_iter=len(history) - 1,
    )


# ---------------------------------------------------------------------------
# Initial-condition factories
# ---------------------------------------------------------------------------

def _orient_facing(face_idx: int, target_dir: NDArray[np.float64]
                   ) -> NDArray[np.float64]:
    """Choose Euler angles so that body-face `face_idx` points along target_dir.

    For simplicity, we align via a lookup: the four body normals are
    fixed; we find Euler angles that rotate normal[face_idx] onto
    target_dir.  Uses Rodrigues' formula then converts to ZYZ Euler.
    """
    _, body_n = k4_body_face_data()
    n0 = body_n[face_idx]
    t = target_dir / np.linalg.norm(target_dir)
    # rotation that takes n0 → t
    v = np.cross(n0, t)
    s = np.linalg.norm(v)
    c = np.dot(n0, t)
    if s < 1e-9:
        if c > 0:
            R = np.eye(3)
        else:
            # 180° around any perpendicular axis
            perp = np.array([1.0, 0.0, 0.0]) if abs(n0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            axis = np.cross(n0, perp)
            axis /= np.linalg.norm(axis)
            K = np.array([[0, -axis[2], axis[1]],
                          [axis[2], 0, -axis[0]],
                          [-axis[1], axis[0], 0]])
            R = np.eye(3) + 2.0 * K @ K
    else:
        K = np.array([[0, -v[2], v[1]],
                      [v[2], 0, -v[0]],
                      [-v[1], v[0], 0]])
        R = np.eye(3) + K + K @ K * ((1.0 - c) / (s ** 2))
    # Convert R back to ZYZ Euler
    beta = np.arccos(np.clip(R[2, 2], -1.0, 1.0))
    if np.sin(beta) > 1e-6:
        alpha = np.arctan2(R[1, 2], R[0, 2])
        gamma = np.arctan2(R[2, 1], -R[2, 0])
    else:
        alpha = np.arctan2(R[1, 0], R[0, 0])
        gamma = 0.0
    return np.array([alpha, beta, gamma])


def make_deuteron(separation: float = 1.4) -> List[K4Cell]:
    """Two K_4 cells facing each other along x — initial deuteron config."""
    # cell 0 at -separation/2, face 0 (its body normal[0]) pointing +x
    # cell 1 at +separation/2, face 0 pointing -x
    r0 = np.array([-separation / 2.0, 0.0, 0.0])
    r1 = np.array([+separation / 2.0, 0.0, 0.0])
    e0 = _orient_facing(0, np.array([+1.0, 0.0, 0.0]))
    e1 = _orient_facing(0, np.array([-1.0, 0.0, 0.0]))
    return [K4Cell(position=r0.copy(), euler=e0),
            K4Cell(position=r1.copy(), euler=e1)]


def make_triton(separation: float = 1.6) -> List[K4Cell]:
    """Three K_4 cells in an equilateral triangle (3 face-pairs)."""
    cells = []
    for k in range(3):
        ang = 2.0 * np.pi * k / 3.0
        r = separation * np.array([np.cos(ang), np.sin(ang), 0.0])
        # face roughly pointing toward centre
        toward = -r / np.linalg.norm(r)
        e = _orient_facing(0, toward)
        cells.append(K4Cell(position=r, euler=e))
    return cells


def make_alpha(separation: float = 1.5) -> List[K4Cell]:
    """Four K_4 cells at the vertices of a tetrahedron (6 face-pairs)."""
    # outer-tetrahedron vertices
    outer = np.array([
        [+1.0, +1.0, +1.0],
        [+1.0, -1.0, -1.0],
        [-1.0, +1.0, -1.0],
        [-1.0, -1.0, +1.0],
    ], dtype=float)
    outer *= separation / np.linalg.norm(outer[0])
    cells = []
    for r in outer:
        toward = -r / np.linalg.norm(r)
        e = _orient_facing(0, toward)
        cells.append(K4Cell(position=r.copy(), euler=e))
    return cells


def make_octet(separation: float = 1.5) -> List[K4Cell]:
    """Eight K_4 cells at the vertices of a cube — toy 8-cell nucleus."""
    cells = []
    for sx in (-1, +1):
        for sy in (-1, +1):
            for sz in (-1, +1):
                r = separation * np.array([sx, sy, sz], dtype=float) / np.sqrt(3.0)
                toward = -r / np.linalg.norm(r)
                e = _orient_facing(0, toward)
                cells.append(K4Cell(position=r, euler=e))
    return cells


def make_hexadecet(separation: float = 1.5) -> List[K4Cell]:
    """Sixteen K_4 cells: two interlocked octets.  Crude initial guess."""
    cells = make_octet(separation=separation)
    # Add a second shell at radius 2.4
    for sx in (-1, +1):
        for sy in (-1, +1):
            for sz in (-1, +1):
                r = 2.4 * separation * np.array([sx, sy, sz], dtype=float) / np.sqrt(3.0)
                toward = -r / np.linalg.norm(r)
                e = _orient_facing(0, toward)
                cells.append(K4Cell(position=r, euler=e))
    return cells


# ---------------------------------------------------------------------------
# Convenience: report
# ---------------------------------------------------------------------------

def binding_energy_report(N_list: Tuple[int, ...] = (2, 3, 4, 8, 16),
                           experimental_be: dict = None
                           ) -> dict:
    """Relax each configuration and report binding energies.

    Args:
        N_list:           Cluster sizes to evaluate.
        experimental_be:  Optional dict {N: BE_exp_MeV} for comparison.

    Returns:
        dict mapping N → {'BE_calc', 'BE_exp', 'rel_error'}.
    """
    if experimental_be is None:
        experimental_be = {2: 2.224, 3: 8.482, 4: 28.296,
                           8: 56.5, 16: 127.6}      # ²H ³H ⁴He ⁸Be(approx) ¹⁶O

    factories = {
        2:  make_deuteron,
        3:  make_triton,
        4:  make_alpha,
        8:  make_octet,
        16: make_hexadecet,
    }
    out = {}
    for N in N_list:
        cells = factories[N]()
        result = relax(cells, n_iter=2000)
        be_calc = result.binding_energy
        be_exp = experimental_be.get(N, None)
        rel = abs(be_calc - be_exp) / be_exp if be_exp else None
        out[N] = {
            'BE_calc': be_calc,
            'BE_exp':  be_exp,
            'rel_error': rel,
            'n_iter':  result.n_iter,
        }
    return out


def configuration_radius(cells: List[K4Cell]) -> float:
    """RMS radius of cell centres about their centroid."""
    P = np.array([c.position for c in cells])
    com = P.mean(axis=0)
    return float(np.sqrt(((P - com) ** 2).sum(axis=1).mean()))


def pairwise_distances(cells: List[K4Cell]) -> NDArray[np.float64]:
    """All pairwise centre-to-centre distances (length n*(n-1)/2)."""
    P = np.array([c.position for c in cells])
    n = len(cells)
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            out.append(np.linalg.norm(P[i] - P[j]))
    return np.array(out)
