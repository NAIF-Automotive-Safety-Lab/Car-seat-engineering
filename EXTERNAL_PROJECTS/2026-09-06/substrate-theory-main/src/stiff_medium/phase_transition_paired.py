"""Paired GEOMETRY + SIM for the cosmological substrate de-saturation epoch.

This module models the (B3) cosmological phase transition in which the
saturated substrate (σ ≈ 1/2 in the pre-Big-Bang era) coherently
de-saturates through the nucleation of "low-σ" bubbles.  The bubbles
expand, collide, and percolate; once a percolating cluster spans the
horizon, the saturation cap is broken globally and ordinary expansion
begins.

PAIRED PATTERN
--------------
``PhaseTransitionGeometry``
    The static lattice picture.  A 2-D Ising-like substrate with sites
    σ_i ∈ {+1 (saturated), -1 (de-saturated)}.  Provides nearest-neighbour
    coupling, a critical temperature derived from the Onsager solution
    on the square lattice,

        T_c (J=1) = 2 / log(1 + √2) ≈ 2.269185

    and clustering / percolation utilities (Hoshen-Kopelman flood-fill).
    The geometry also exposes ``critical_sigma`` ≡ the order-parameter
    threshold |⟨σ⟩| < σ_c at which the percolating cluster appears
    (here we use σ_c = 0 since the Ising transition is symmetric, but
    we expose it as a parameter to allow biased / cosmological variants).

``PhaseTransitionSimulator``
    The dynamics.  Metropolis-Hastings (default) or Glauber dynamics on
    the geometry.  Records:

      * the spin grid history (snapshots),
      * the order parameter ⟨σ⟩(t),
      * the largest cluster size L(t),
      * the bubble (cluster of -1 sites) size distribution,
      * the nucleation rate Γ(T): bubbles per unit area per sweep.

Both objects share the lattice; the geometry is the *kinematic* layer
(neighbours, clusters, T_c) and the simulator is the *dynamic* layer
(MC moves, history).

PHYSICAL CORRESPONDENCE (B3 framework)
--------------------------------------
* Hot phase  T ≫ T_c, ⟨σ⟩ ≈ 0:  random-walk substrate, no coherence.
* Symmetric saturated phase  T = 0, ⟨σ⟩ = +1:  the cap is locked,
  no expansion.  This is the pre-Big-Bang state of the substrate
  cosmology (cosmology_simulator.py).
* The cosmological de-saturation epoch sits *just below* T_c on the
  cooling side: bubbles of -1 (low-σ) nucleate inside the +1 sea,
  grow, collide, and percolate.  Once percolation occurs, the
  saturation cap is broken on the horizon scale and a(t) starts
  growing — this is the matching condition with the de-saturation
  rate in cosmology_simulator.SubstrateParams.

The output is intentionally pedagogical / dimensionless: J = 1 sets the
temperature scale and the lattice unit is the substrate coherence length
ξ. No new free parameter is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Constants — Onsager critical temperature for the 2-D square Ising model
# ---------------------------------------------------------------------------

J_COUPLING_DEFAULT: float = 1.0
T_CRITICAL_ONSAGER: float = 2.0 / np.log(1.0 + np.sqrt(2.0))  # ≈ 2.269185


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

@dataclass
class PhaseTransitionGeometry:
    """2-D Ising-like substrate with σ ∈ {+1, -1} as order parameter.

    Parameters
    ----------
    L : int
        Linear lattice size (the lattice has L * L sites).
    J : float
        Nearest-neighbour exchange coupling. ``J = 1`` is the canonical
        choice; the critical temperature scales as T_c ∝ J.
    h : float
        External bias (analogous to a chemical potential favouring the
        de-saturated phase). ``h = 0`` keeps the symmetric Ising point.
    sigma_critical : float
        Order-parameter threshold (in the ⟨σ⟩ sense) below which the
        ``-1`` cluster is expected to percolate.  Defaults to ``0`` —
        the symmetric Ising transition.
    seed : int
        RNG seed used to initialise the lattice.

    Notes
    -----
    The lattice is square with periodic boundary conditions.  Each cell
    has four neighbours (von Neumann).  This is the geometric layer:
    no dynamics are stored here, only the static structure plus
    cluster / percolation queries.
    """

    L: int = 48
    J: float = J_COUPLING_DEFAULT
    h: float = 0.0
    sigma_critical: float = 0.0
    seed: int = 0

    grid: NDArray[np.int8] = field(init=False)

    def __post_init__(self) -> None:
        rng = np.random.default_rng(self.seed)
        # Start saturated (+1) by default; simulator will optionally seed
        # the lattice differently (e.g. random or biased).
        self.grid = np.ones((self.L, self.L), dtype=np.int8)
        # Tiny initial noise to break perfect symmetry without changing ⟨σ⟩.
        flips = rng.random((self.L, self.L)) < 0.001
        self.grid[flips] = -1

    # ------------------------------------------------------------------
    # Static / kinematic queries
    # ------------------------------------------------------------------
    @property
    def critical_temperature(self) -> float:
        """Onsager T_c for the 2-D square Ising model in units of J."""
        return self.J * T_CRITICAL_ONSAGER

    @property
    def critical_sigma(self) -> float:
        """Order-parameter threshold for the percolating bubble cluster."""
        return self.sigma_critical

    def neighbours(self, i: int, j: int) -> List[Tuple[int, int]]:
        """Return the 4 von-Neumann neighbours with periodic BCs."""
        L = self.L
        return [
            ((i + 1) % L, j),
            ((i - 1) % L, j),
            (i, (j + 1) % L),
            (i, (j - 1) % L),
        ]

    def order_parameter(self, grid: Optional[NDArray] = None) -> float:
        """Mean spin ⟨σ⟩ on the supplied grid (defaults to self.grid)."""
        g = self.grid if grid is None else grid
        return float(np.mean(g))

    def energy(self, grid: Optional[NDArray] = None) -> float:
        """Ising Hamiltonian H = -J Σ σ_i σ_j - h Σ σ_i (per spin)."""
        g = self.grid if grid is None else grid
        right = np.roll(g, -1, axis=1)
        down = np.roll(g, -1, axis=0)
        bond = -self.J * (g * right + g * down).sum()
        field_term = -self.h * g.sum()
        return float((bond + field_term) / g.size)

    # ------------------------------------------------------------------
    # Cluster / percolation utilities (Hoshen-Kopelman style flood-fill)
    # ------------------------------------------------------------------
    def label_clusters(
        self,
        grid: Optional[NDArray] = None,
        target: int = -1,
    ) -> Tuple[NDArray[np.int32], Dict[int, int]]:
        """Label connected ``target``-sites and return (labels, sizes).

        Parameters
        ----------
        grid : ndarray, optional
            Spin grid; defaults to ``self.grid``.
        target : int
            Spin value identifying the bubble interior (``-1`` for
            de-saturated phase by default).

        Returns
        -------
        labels : ndarray of int
            Same shape as the grid; each cell carries its cluster id
            (0 means "not in target phase").
        sizes : dict
            Mapping ``cluster_id -> cluster_size``.
        """
        g = self.grid if grid is None else grid
        L = g.shape[0]
        labels = np.zeros((L, L), dtype=np.int32)
        sizes: Dict[int, int] = {}
        next_id = 1
        for i in range(L):
            for j in range(L):
                if g[i, j] != target or labels[i, j] != 0:
                    continue
                # iterative flood-fill (BFS)
                stack = [(i, j)]
                labels[i, j] = next_id
                sizes[next_id] = 0
                while stack:
                    ci, cj = stack.pop()
                    sizes[next_id] += 1
                    for ni, nj in self.neighbours(ci, cj):
                        if g[ni, nj] == target and labels[ni, nj] == 0:
                            labels[ni, nj] = next_id
                            stack.append((ni, nj))
                next_id += 1
        return labels, sizes

    def largest_cluster_size(
        self, grid: Optional[NDArray] = None, target: int = -1
    ) -> int:
        """Largest contiguous bubble (target-spin) cluster."""
        _, sizes = self.label_clusters(grid=grid, target=target)
        return max(sizes.values()) if sizes else 0

    def has_percolating_cluster(
        self, grid: Optional[NDArray] = None, target: int = -1
    ) -> bool:
        """True iff any ``target`` cluster spans both lattice axes.

        Spans = the cluster touches all rows OR all columns (we check
        both directions, periodic BCs make the distinction subtle but
        for the snapshot we use the simpler "touches both extremes").
        """
        labels, sizes = self.label_clusters(grid=grid, target=target)
        if not sizes:
            return False
        for cid in sizes:
            mask = labels == cid
            spans_rows = mask.any(axis=1).all()
            spans_cols = mask.any(axis=0).all()
            if spans_rows or spans_cols:
                return True
        return False


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

@dataclass
class PhaseTransitionSimulator:
    """Metropolis or Glauber dynamics on the phase-transition lattice.

    Parameters
    ----------
    geometry : PhaseTransitionGeometry
        The lattice substrate (provides spin field and neighbour rules).
    T : float
        Bath temperature in units of J (so T = T_critical at the
        Onsager point ≈ 2.269 J).
    n_sweeps : int
        Number of full lattice sweeps (one sweep = L² site updates).
    dynamics : str
        ``"metropolis"`` or ``"glauber"``.
    snapshot_every : int
        Record full grid snapshots every N sweeps.
    seed_grid : str
        ``"saturated"`` starts σ = +1 everywhere (B3 pre-BB state),
        ``"random"`` initialises uniformly in {±1},
        ``"hot"`` is identical to ``"random"``.
    seed : int
        RNG seed for the dynamics.
    bubble_threshold : int
        Minimum cluster size considered a "bubble" (for nucleation
        statistics we ignore single-site flips).
    """

    geometry: PhaseTransitionGeometry
    T: float = 2.10
    n_sweeps: int = 200
    dynamics: str = "metropolis"
    snapshot_every: int = 10
    seed_grid: str = "saturated"
    seed: int = 0
    bubble_threshold: int = 2

    history: Dict[str, List] = field(init=False, default_factory=dict)

    def __post_init__(self) -> None:
        if self.dynamics not in ("metropolis", "glauber"):
            raise ValueError(f"unknown dynamics {self.dynamics!r}")
        self.history = {
            "sweep": [],
            "order_parameter": [],
            "energy": [],
            "largest_cluster": [],
            "n_bubbles": [],
            "snapshots": [],
            "snapshot_sweeps": [],
        }
        self._rng = np.random.default_rng(self.seed)
        self._init_lattice()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------
    def _init_lattice(self) -> None:
        L = self.geometry.L
        if self.seed_grid in ("saturated", "+1"):
            self.geometry.grid = np.ones((L, L), dtype=np.int8)
        elif self.seed_grid in ("random", "hot"):
            self.geometry.grid = np.where(
                self._rng.random((L, L)) < 0.5, 1, -1
            ).astype(np.int8)
        elif self.seed_grid in ("desaturated", "-1"):
            self.geometry.grid = -np.ones((L, L), dtype=np.int8)
        else:
            raise ValueError(f"unknown seed_grid {self.seed_grid!r}")

    # ------------------------------------------------------------------
    # Local energy change for a flip
    # ------------------------------------------------------------------
    def _local_dE(self, i: int, j: int) -> float:
        g = self.geometry.grid
        s = g[i, j]
        L = self.geometry.L
        nbsum = (
            g[(i + 1) % L, j]
            + g[(i - 1) % L, j]
            + g[i, (j + 1) % L]
            + g[i, (j - 1) % L]
        )
        # ΔE for flipping σ_i: H = -J σ_i Σ σ_j - h σ_i  →  ΔE = 2 σ_i (J Σ σ_j + h)
        return 2.0 * s * (self.geometry.J * nbsum + self.geometry.h)

    # ------------------------------------------------------------------
    # MC step (one sweep = L² random-site updates)
    # ------------------------------------------------------------------
    def _sweep(self) -> None:
        L = self.geometry.L
        T = max(self.T, 1e-9)
        beta = 1.0 / T
        # Random site order
        idx = self._rng.integers(0, L, size=(L * L, 2))
        rand = self._rng.random(L * L)
        for k in range(L * L):
            i, j = int(idx[k, 0]), int(idx[k, 1])
            dE = self._local_dE(i, j)
            if self.dynamics == "metropolis":
                if dE <= 0.0 or rand[k] < np.exp(-beta * dE):
                    self.geometry.grid[i, j] = -self.geometry.grid[i, j]
            else:  # Glauber
                p_flip = 1.0 / (1.0 + np.exp(beta * dE))
                if rand[k] < p_flip:
                    self.geometry.grid[i, j] = -self.geometry.grid[i, j]

    # ------------------------------------------------------------------
    # Public driver
    # ------------------------------------------------------------------
    def run(self) -> Dict[str, List]:
        for s in range(self.n_sweeps):
            self._sweep()
            sigma_mean = self.geometry.order_parameter()
            energy = self.geometry.energy()
            _, sizes = self.geometry.label_clusters(target=-1)
            largest = max(sizes.values()) if sizes else 0
            n_bubbles = sum(
                1 for v in sizes.values() if v >= self.bubble_threshold
            )
            self.history["sweep"].append(s)
            self.history["order_parameter"].append(sigma_mean)
            self.history["energy"].append(energy)
            self.history["largest_cluster"].append(largest)
            self.history["n_bubbles"].append(n_bubbles)
            if (s % self.snapshot_every) == 0 or s == self.n_sweeps - 1:
                self.history["snapshots"].append(self.geometry.grid.copy())
                self.history["snapshot_sweeps"].append(s)
        return self.history

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    def bubble_size_distribution(self) -> NDArray[np.int64]:
        """Cluster sizes of the de-saturated phase in the *current* grid."""
        _, sizes = self.geometry.label_clusters(target=-1)
        return np.array(sorted(sizes.values()), dtype=np.int64)

    def percolated(self) -> bool:
        return self.geometry.has_percolating_cluster(target=-1)

    def nucleation_rate(self) -> float:
        """Bubbles per unit area per sweep, averaged over the run."""
        if not self.history["n_bubbles"]:
            return 0.0
        area = float(self.geometry.L ** 2)
        return float(np.mean(self.history["n_bubbles"])) / area


# ---------------------------------------------------------------------------
# Convenience: temperature-scan utility used by both tests and the visualiser
# ---------------------------------------------------------------------------

def nucleation_rate_vs_T(
    temperatures: NDArray[np.float64],
    L: int = 32,
    n_sweeps: int = 80,
    bubble_threshold: int = 2,
    seed: int = 1,
) -> Dict[str, NDArray[np.float64]]:
    """Run a temperature scan starting from the saturated state.

    Returns a dict with ``T``, ``rate`` (mean bubbles / area / sweep),
    ``order_parameter`` (final ⟨σ⟩), and ``largest_cluster_frac``.
    """
    Ts = np.asarray(temperatures, dtype=float)
    rates = np.zeros_like(Ts)
    sigmas = np.zeros_like(Ts)
    largest = np.zeros_like(Ts)
    for k, T in enumerate(Ts):
        geom = PhaseTransitionGeometry(L=L, seed=seed + k)
        sim = PhaseTransitionSimulator(
            geometry=geom,
            T=float(T),
            n_sweeps=n_sweeps,
            snapshot_every=max(n_sweeps, 1),
            seed_grid="saturated",
            seed=seed + k,
            bubble_threshold=bubble_threshold,
        )
        sim.run()
        rates[k] = sim.nucleation_rate()
        sigmas[k] = geom.order_parameter()
        largest[k] = sim.history["largest_cluster"][-1] / (L * L)
    return {
        "T": Ts,
        "rate": rates,
        "order_parameter": sigmas,
        "largest_cluster_frac": largest,
    }


__all__ = [
    "T_CRITICAL_ONSAGER",
    "PhaseTransitionGeometry",
    "PhaseTransitionSimulator",
    "nucleation_rate_vs_T",
]
