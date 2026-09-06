"""Phonon dispersion calculator for substrate field linear excitations.

Linearizes substrate Lagrangian L = (ρ/2)(∂_t φ)² - (K/2)(∇φ)² - (ξ/2) φ²
around vacuum on a discrete lattice. Diagonalizes the dynamical matrix
D(k) at each k in the first Brillouin zone to extract acoustic and
optical phonon branches.

Physical correspondence (substrate -> condensed-matter):
    K  = elastic stiffness  -> spring constant scaling
    ρ  = inertial density   -> on-site mass
    ξ  = optical mass term  -> on-site restoring spring (gives optical gap)

Long-wavelength limit recovers c_s = sqrt(K/ρ).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
from scipy.linalg import eigh

# ---------------------------------------------------------------------------
# Lattice geometry definitions
# ---------------------------------------------------------------------------

# Each lattice gives:
#   a:        primitive lattice vectors (rows)
#   basis:    fractional coords of basis atoms (n_atoms × dim)
#   neighbors: list of (delta_R, weight) — Cartesian connection vectors
#              between basis atom pairs (i, j, dR_cartesian, weight)
#   high_sym: dict of high-symmetry k-point names -> Cartesian k-vectors

LatticeData = Dict[str, object]


def _chain_1d(a: float = 1.0, n_basis: int = 1) -> LatticeData:
    a_vec = np.array([[a]])
    if n_basis == 1:
        basis = np.array([[0.0]])
        bonds = [(0, 0, np.array([+a])), (0, 0, np.array([-a]))]
    else:
        basis = np.array([[0.0], [0.5 * a]])
        bonds = [
            (0, 1, np.array([+0.5 * a])),
            (0, 1, np.array([-0.5 * a])),
            (1, 0, np.array([+0.5 * a])),
            (1, 0, np.array([-0.5 * a])),
        ]
    high_sym = {"G": np.array([0.0]), "X": np.array([np.pi / a])}
    return dict(a=a_vec, basis=basis, bonds=bonds, high_sym=high_sym, dim=1)


def _square_2d(a: float = 1.0) -> LatticeData:
    a_vec = a * np.eye(2)
    basis = np.array([[0.0, 0.0]])
    bonds = [
        (0, 0, np.array([+a, 0.0])),
        (0, 0, np.array([-a, 0.0])),
        (0, 0, np.array([0.0, +a])),
        (0, 0, np.array([0.0, -a])),
    ]
    high_sym = {
        "G": np.array([0.0, 0.0]),
        "X": np.array([np.pi / a, 0.0]),
        "M": np.array([np.pi / a, np.pi / a]),
    }
    return dict(a=a_vec, basis=basis, bonds=bonds, high_sym=high_sym, dim=2)


def _hex_2d(a: float = 1.0) -> LatticeData:
    a1 = a * np.array([1.0, 0.0])
    a2 = a * np.array([0.5, np.sqrt(3) / 2])
    a_vec = np.vstack([a1, a2])
    basis = np.array([[0.0, 0.0]])
    angs = np.deg2rad([0, 60, 120, 180, 240, 300])
    bonds = [(0, 0, a * np.array([np.cos(t), np.sin(t)])) for t in angs]
    high_sym = {
        "G": np.array([0.0, 0.0]),
        "M": np.array([2 * np.pi / a, 0.0]) / np.sqrt(3),
        "K": np.array([4 * np.pi / (3 * a), 0.0]),
    }
    return dict(a=a_vec, basis=basis, bonds=bonds, high_sym=high_sym, dim=2)


def _fcc_3d(a: float = 1.0) -> LatticeData:
    a_vec = (a / 2) * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    basis = np.array([[0.0, 0.0, 0.0]])
    nn = (a / 2) * np.array(
        [
            [1, 1, 0], [-1, -1, 0], [1, -1, 0], [-1, 1, 0],
            [1, 0, 1], [-1, 0, -1], [1, 0, -1], [-1, 0, 1],
            [0, 1, 1], [0, -1, -1], [0, 1, -1], [0, -1, 1],
        ],
        dtype=float,
    )
    bonds = [(0, 0, v) for v in nn]
    high_sym = {
        "G": np.zeros(3),
        "X": (2 * np.pi / a) * np.array([1.0, 0.0, 0.0]),
        "L": (np.pi / a) * np.array([1.0, 1.0, 1.0]),
        "K": (2 * np.pi / a) * np.array([0.75, 0.75, 0.0]),
    }
    return dict(a=a_vec, basis=basis, bonds=bonds, high_sym=high_sym, dim=3)


def _bcc_3d(a: float = 1.0) -> LatticeData:
    a_vec = (a / 2) * np.array(
        [[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=float
    )
    basis = np.array([[0.0, 0.0, 0.0]])
    nn = (a / 2) * np.array(
        [
            [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1],
            [-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1],
        ],
        dtype=float,
    )
    bonds = [(0, 0, v) for v in nn]
    high_sym = {
        "G": np.zeros(3),
        "H": (2 * np.pi / a) * np.array([1.0, 0.0, 0.0]),
        "N": (np.pi / a) * np.array([1.0, 1.0, 0.0]),
        "P": (np.pi / a) * np.array([1.0, 1.0, 1.0]),
    }
    return dict(a=a_vec, basis=basis, bonds=bonds, high_sym=high_sym, dim=3)


def _diamond_3d(a: float = 1.0) -> LatticeData:
    """Diamond: FCC with two-atom basis at (0,0,0) and (a/4)(1,1,1)."""
    a_vec = (a / 2) * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    d = a / 4
    basis = np.array([[0.0, 0.0, 0.0], [d, d, d]])
    # Each atom 0 connects to 4 atom-1 neighbors at tetrahedral angles
    tet = (a / 4) * np.array(
        [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
    )
    bonds = []
    for v in tet:
        bonds.append((0, 1, v))
        bonds.append((1, 0, -v))
    high_sym = {
        "G": np.zeros(3),
        "X": (2 * np.pi / a) * np.array([1.0, 0.0, 0.0]),
        "L": (np.pi / a) * np.array([1.0, 1.0, 1.0]),
        "K": (2 * np.pi / a) * np.array([0.75, 0.75, 0.0]),
    }
    return dict(a=a_vec, basis=basis, bonds=bonds, high_sym=high_sym, dim=3)


_LATTICE_FACTORIES: Dict[str, Callable[..., LatticeData]] = {
    "chain": _chain_1d,
    "chain_diatomic": lambda a=1.0: _chain_1d(a=a, n_basis=2),
    "square": _square_2d,
    "hex": _hex_2d,
    "fcc": _fcc_3d,
    "bcc": _bcc_3d,
    "diamond": _diamond_3d,
}


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------


@dataclass
class PhononDispersion:
    """Phonon dispersion on a substrate-field lattice.

    Parameters
    ----------
    K : float
        Substrate elastic stiffness (units of energy/length²).
    rho : float
        Substrate inertial density (units of mass/length^d).
    xi : float
        Optical (on-site) mass term (gives optical-mode gap).
    lattice : str
        One of: 'chain', 'chain_diatomic', 'square', 'hex',
        'fcc', 'bcc', 'diamond'.
    a : float
        Lattice constant.
    masses : sequence of float, optional
        Per-basis-atom masses (default: all = rho).
    """

    K: float
    rho: float
    xi: float = 0.0
    lattice: str = "fcc"
    a: float = 1.0
    masses: Sequence[float] = field(default=None)

    def __post_init__(self):
        if self.lattice not in _LATTICE_FACTORIES:
            raise ValueError(f"Unknown lattice '{self.lattice}'")
        self._geo: LatticeData = _LATTICE_FACTORIES[self.lattice](a=self.a)
        self.dim: int = int(self._geo["dim"])
        self.basis: np.ndarray = np.asarray(self._geo["basis"])
        self.n_atoms: int = self.basis.shape[0]
        self.bonds: List[Tuple[int, int, np.ndarray]] = list(self._geo["bonds"])
        self.high_sym: Dict[str, np.ndarray] = dict(self._geo["high_sym"])
        if self.masses is None:
            self.masses = np.full(self.n_atoms, self.rho)
        else:
            self.masses = np.asarray(self.masses, dtype=float)
            if self.masses.shape[0] != self.n_atoms:
                raise ValueError("masses length must equal n_atoms")
        # n DOF per atom = dim (longitudinal + transverse displacements).
        self.dof_per_atom: int = self.dim
        self.n_modes: int = self.n_atoms * self.dof_per_atom
        # Degeneracies on bonds for proper k=0 acoustic-mode test
        # (sum of all bond force constants on each atom -> "self" term).

    # -- core dynamical matrix ------------------------------------------------

    def dynamical_matrix(self, k: np.ndarray) -> np.ndarray:
        """Build D(k): (n_modes × n_modes) Hermitian matrix.

        Uses harmonic-spring force-constant model:
            Φ_{iα,jβ}(R) = -K * (e_R)_α (e_R)_β     for bond (i,j) at R
        Self-term enforces acoustic sum rule:
            Φ_{iα,iα}(0) = -Σ_{j,R≠} Φ_{ij,αβ}(R)
        Plus optical on-site term ξ δ_{ij} δ_{αβ}.
        Returns D = M^{-1/2} Φ(k) M^{-1/2}.
        """
        k = np.atleast_1d(np.asarray(k, dtype=float))
        n = self.n_modes
        d = self.dof_per_atom
        Phi = np.zeros((n, n), dtype=complex)

        # Off-diagonal hops: Σ_R Φ_{ij}(R) e^{i k·R}
        for (i, j, R) in self.bonds:
            r = np.linalg.norm(R)
            if r < 1e-14:
                continue
            ehat = R / r
            phase = np.exp(1j * np.dot(k, R))
            block = -self.K * np.outer(ehat, ehat) * phase
            Phi[i * d:(i + 1) * d, j * d:(j + 1) * d] += block

        # Self / on-site term: enforce acoustic sum rule
        # Self_{i} = -Σ_{j,R} Φ_{ij}(R)|_{phase=1}  + ξ I
        for i in range(self.n_atoms):
            self_block = np.zeros((d, d), dtype=complex)
            for (ii, jj, R) in self.bonds:
                if ii != i:
                    continue
                r = np.linalg.norm(R)
                if r < 1e-14:
                    continue
                ehat = R / r
                self_block += self.K * np.outer(ehat, ehat)
            self_block += self.xi * np.eye(d)
            Phi[i * d:(i + 1) * d, i * d:(i + 1) * d] += self_block

        # Mass-weight: D = M^{-1/2} Φ M^{-1/2}
        m_vec = np.repeat(self.masses, d)
        m_inv_sqrt = 1.0 / np.sqrt(m_vec)
        D = (m_inv_sqrt[:, None] * Phi) * m_inv_sqrt[None, :]
        # Symmetrize numerically
        D = 0.5 * (D + D.conj().T)
        return D

    def frequencies(self, k: np.ndarray) -> np.ndarray:
        """Eigenfrequencies ω(k) (sorted ascending). ω = sqrt(max(λ, 0))."""
        D = self.dynamical_matrix(k)
        evals = eigh(D, eigvals_only=True)
        evals = np.where(evals < 0, 0.0, evals)
        return np.sqrt(evals)

    def modes(self, k: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Return (frequencies, eigenvectors). Columns are mode vectors."""
        D = self.dynamical_matrix(k)
        evals, evecs = eigh(D)
        evals = np.where(evals < 0, 0.0, evals)
        return np.sqrt(evals), evecs

    # -- API methods ----------------------------------------------------------

    def dispersion_curve(
        self, path: Sequence[str], n_points: int = 50
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ω(k) along a k-path through high-symmetry points.

        Parameters
        ----------
        path : sequence of point names (e.g. ['G', 'X', 'M', 'G'])
        n_points : k-points per segment.

        Returns
        -------
        x      : (N,) cumulative k-distance
        omega  : (N, n_modes) frequencies
        ticks  : (n_segments+1,) x-positions of high-symmetry points
        """
        pts = [np.atleast_1d(self.high_sym[p]).astype(float) for p in path]
        ticks = [0.0]
        xs: List[float] = []
        ks: List[np.ndarray] = []
        cum = 0.0
        for i in range(len(pts) - 1):
            k0, k1 = pts[i], pts[i + 1]
            seg_len = float(np.linalg.norm(k1 - k0))
            for s in range(n_points):
                t = s / n_points
                ks.append(k0 + t * (k1 - k0))
                xs.append(cum + t * seg_len)
            cum += seg_len
            ticks.append(cum)
        ks.append(pts[-1])
        xs.append(cum)
        omega = np.array([self.frequencies(k) for k in ks])
        return np.asarray(xs), omega, np.asarray(ticks)

    def density_of_states(
        self,
        omega_range: Tuple[float, float] | None = None,
        n_bins: int = 100,
        n_kgrid: int = 12,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Phonon DOS via Monkhorst-Pack-like uniform grid + histogram.

        Returns (omega_centers, dos). DOS is normalized so that
        sum(dos) * dω == n_modes (total number of branches).
        """
        # Uniform grid in fractional reciprocal-lattice coordinates [-0.5, 0.5)
        axes = [np.linspace(-0.5, 0.5, n_kgrid, endpoint=False)] * self.dim
        mesh = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(
            -1, self.dim
        )
        # Reciprocal lattice
        a_mat = np.asarray(self._geo["a"])
        b_mat = 2 * np.pi * np.linalg.inv(a_mat).T  # rows: b_i
        all_freqs: List[float] = []
        for frac in mesh:
            k = frac @ b_mat
            all_freqs.extend(self.frequencies(k).tolist())
        all_freqs = np.asarray(all_freqs)
        if omega_range is None:
            omega_range = (0.0, all_freqs.max() * 1.01)
        hist, edges = np.histogram(all_freqs, bins=n_bins, range=omega_range)
        centers = 0.5 * (edges[:-1] + edges[1:])
        n_k = mesh.shape[0]
        domega = edges[1] - edges[0]
        # Normalize so ∫ dos dω = n_modes (total branches per unit cell)
        dos = hist.astype(float) / (n_k * domega)
        return centers, dos

    def sound_speed(self, direction: np.ndarray | None = None,
                    k_small: float = 1e-3) -> float:
        """Recover acoustic sound speed c_s = ω/|k| at small k.

        Returns the average over the lowest dim acoustic branches.
        """
        if direction is None:
            direction = np.zeros(self.dim)
            direction[0] = 1.0
        direction = np.asarray(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)
        k = k_small * direction
        omega = self.frequencies(k)
        # Acoustic = lowest dim branches
        acoustic = omega[: self.dim]
        return float(np.mean(acoustic) / k_small)

    def debye_temperature(self, hbar: float = 1.0,
                          kB: float = 1.0,
                          n_kgrid: int = 12) -> float:
        """Debye temperature Θ_D = (ħ/k_B) ω_D where ω_D from total mode count.

        Uses average sound speed via Debye formula:
            ω_D = c_s * (6 π² N/V)^{1/3}   (3D)
        For 2D / 1D, uses dimension-appropriate formula. Returns Θ_D.
        """
        # Average sound speed over several Cartesian directions
        if self.dim == 3:
            dirs = np.array(
                [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                 [1, 1, 0], [1, 0, 1], [0, 1, 1],
                 [1, 1, 1]], dtype=float
            )
        elif self.dim == 2:
            dirs = np.array([[1, 0], [0, 1], [1, 1]], dtype=float)
        else:
            dirs = np.array([[1.0]])
        cs_vals = [self.sound_speed(d) for d in dirs]
        c_s = float(np.mean(cs_vals))

        # Number density n = n_atoms / V_cell
        a_mat = np.asarray(self._geo["a"])
        if self.dim == 3:
            V = abs(np.linalg.det(a_mat))
            n_density = self.n_atoms / V
            omega_D = c_s * (6 * np.pi**2 * n_density) ** (1.0 / 3.0)
        elif self.dim == 2:
            V = abs(np.linalg.det(a_mat))
            n_density = self.n_atoms / V
            omega_D = c_s * (4 * np.pi * n_density) ** 0.5
        else:
            V = abs(a_mat[0, 0])
            n_density = self.n_atoms / V
            omega_D = c_s * np.pi * n_density
        return float(hbar * omega_D / kB)

    def acoustic_branches(self, k: np.ndarray) -> np.ndarray:
        """Return the dim lowest (acoustic) branches at given k."""
        omega = self.frequencies(k)
        return omega[: self.dim]

    def optical_branches(self, k: np.ndarray) -> np.ndarray:
        """Return the optical branches: 3(n_atoms-1) for 3D, etc."""
        omega = self.frequencies(k)
        if self.n_atoms == 1:
            return np.array([])
        return omega[self.dim:]


# ---------------------------------------------------------------------------
# Convenience: physical-units check for real materials
# ---------------------------------------------------------------------------

def debye_T_from_real_material(
    c_s_avg_m_per_s: float,
    n_atoms_per_m3: float,
    hbar: float = 1.054571817e-34,
    kB: float = 1.380649e-23,
) -> float:
    """Standard Debye-temperature formula in SI for cross-checking."""
    omega_D = c_s_avg_m_per_s * (6 * np.pi**2 * n_atoms_per_m3) ** (1 / 3)
    return hbar * omega_D / kB


def estimate_debye_diamond() -> float:
    """Diamond cross-check: c_s ≈ 12000 m/s, n ≈ 1.76e29 atoms/m³.

    Returns Θ_D in kelvin. Expected ≈ 2230 K.
    """
    return debye_T_from_real_material(12000.0, 1.76e29)


def estimate_debye_aluminum() -> float:
    """Aluminum: c_s ≈ 5100 m/s, n ≈ 6.02e28 atoms/m³. Expected ≈ 428 K."""
    return debye_T_from_real_material(5100.0, 6.02e28)


def estimate_debye_copper() -> float:
    """Copper: Debye-averaged c_s ≈ 2620 m/s, n ≈ 8.49e28 atoms/m³.
    Expected ≈ 343 K. (Longitudinal c_L ≈ 4760, transverse c_T ≈ 2325;
    Debye average is dominated by transverse.)"""
    return debye_T_from_real_material(2620.0, 8.49e28)


def estimate_debye_aluminum_v2() -> float:
    """Aluminum Debye-averaged c_s ≈ 3450 m/s, n ≈ 6.02e28 atoms/m³."""
    return debye_T_from_real_material(3450.0, 6.02e28)
