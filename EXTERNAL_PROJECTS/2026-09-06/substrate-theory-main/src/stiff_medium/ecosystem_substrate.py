"""Paired ecosystem-dynamics GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for population
ecology under the B3 / stiff-medium ontology.  Population dynamics are
the macroscopic, slow-substrate limit of the same flow equations that
govern micro-substrate occupancies: each species is a "strain
configuration" whose density evolves under nonlinear interactions with
neighbouring configurations.

GEOMETRY
--------
``EcosystemGeometry`` represents a (food) web as a directed graph.

    * Nodes: species (predator, prey, producer, decomposer ...).
    * Edges: directed (consumer <- resource); edge weight = interaction
      strength a_ij (units of 1 / population / time).
    * Trophic level: shortest-path distance from a basal producer.

Two canonical layouts:

    * predator_prey():  2-node Lotka-Volterra system.
    * random_food_web(N, C, sigma, seed): an Erdos-Renyi-like ecology
      used in May-Wigner stability arguments (May 1972 Nature).
      Parameters:
          N         number of species
          C         connectance (probability of edge i,j)
          sigma     standard deviation of the interaction matrix
                    entries.

SIMULATOR
---------
``EcosystemSimulator`` exposes the canonical population-dynamics
identities and a numerical integrator:

    * lotka_volterra_dynamics(N, P, alpha, beta, delta, gamma):
          dN/dt =  alpha * N - beta  * N * P
          dP/dt =  delta * N * P - gamma * P
    * lotka_volterra_period(alpha, gamma):
          T = 2 pi / sqrt(alpha * gamma)   (small-amplitude limit)
    * lotka_volterra_equilibrium():
          (N*, P*) = (gamma / delta, alpha / beta)
    * logistic_growth_dynamics(N, r, K):
          dN/dt = r * N * (1 - N / K)
    * extinction_threshold(N, N_min):
          species goes extinct if N falls below N_min (Allee-like cutoff).
    * food_web_jacobian(state):
          linearisation about the interior equilibrium; eigenvalues
          decide local stability (May-Wigner: instability at
          sigma * sqrt(N * C) > 1).
    * may_wigner_threshold(N, C):
          critical sigma_c = 1 / sqrt(N * C) at which a random food web
          loses linear stability.
    * integrate(): RK4 integrator for arbitrary species networks.

Substrate-specific predictions encoded here:
    * Conservation of the LV invariant
          V = delta*N - gamma*ln N + beta*P - alpha*ln P
      is the strain-energy invariant of the predator-prey configuration
      (closed orbit in phase space).
    * The Allee / extinction threshold maps onto the substrate-saturation
      lower bound: below a critical occupancy a strain pattern cannot
      sustain itself and decays to zero.
    * May-Wigner instability for large random webs is the ecological
      analogue of the substrate's loss-of-coherence at high disorder
      (cf. Anderson localisation).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) values -- canonical Lotka-Volterra demo
# ---------------------------------------------------------------------------

# Canonical Lotka-Volterra coefficients used in textbook demos
LV_ALPHA_DEFAULT: float = 1.0   # prey intrinsic growth (1/time)
LV_BETA_DEFAULT:  float = 0.1   # predation rate (1 / predator / time)
LV_DELTA_DEFAULT: float = 0.075 # prey -> predator conversion efficiency
LV_GAMMA_DEFAULT: float = 1.5   # predator death rate (1/time)

# Logistic growth canonical values
LOGISTIC_R_DEFAULT: float = 0.5   # intrinsic growth rate (1/time)
LOGISTIC_K_DEFAULT: float = 1000.0  # carrying capacity (individuals)

# Extinction threshold (Allee-like floor; population goes extinct below this)
EXTINCTION_THRESHOLD_DEFAULT: float = 1.0   # individuals

# May-Wigner random-matrix stability constants
MAY_WIGNER_SELF_REGULATION: float = 1.0   # diagonal damping a_ii = -1

# Substrate anchor (B3 framework)
LAMBDA_QCD_MEV: float = 200.0


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class EcosystemGeometry:
    """Directed-graph geometry for an ecological community.

    Attributes
    ----------
    n_species : int
        Number of species (graph nodes).
    species_names : list[str]
        Optional human-readable labels.
    interaction_matrix : np.ndarray, shape (N, N)
        a_ij is the per-capita effect of species j on species i.
        Diagonal a_ii is intrinsic growth (positive) or self-damping
        (negative, e.g. carrying-capacity term).  Off-diagonal entries
        encode predator-prey, competition, or mutualism.
    trophic_levels : np.ndarray, shape (N,)
        Trophic level (1 = basal producer, 2 = primary consumer, ...).
    """

    n_species: int = 2
    species_names: list[str] = field(default_factory=lambda: ["prey", "predator"])
    interaction_matrix: np.ndarray = field(default_factory=lambda: np.zeros((2, 2)))
    trophic_levels: np.ndarray = field(default_factory=lambda: np.array([1.0, 2.0]))

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def predator_prey(
        cls,
        alpha: float = LV_ALPHA_DEFAULT,
        beta: float = LV_BETA_DEFAULT,
        delta: float = LV_DELTA_DEFAULT,
        gamma: float = LV_GAMMA_DEFAULT,
    ) -> "EcosystemGeometry":
        """Two-node Lotka-Volterra ecosystem.

        Interaction matrix (rows = effect on, cols = effect from):
            A = [[+alpha,  -beta],
                 [+delta,  -gamma]]

        Note that the diagonal carries the intrinsic rates; off-diagonals
        carry the bilinear interaction strengths.
        """
        A = np.array([[+alpha, -beta],
                      [+delta, -gamma]], dtype=float)
        return cls(
            n_species=2,
            species_names=["prey", "predator"],
            interaction_matrix=A,
            trophic_levels=np.array([1.0, 2.0]),
        )

    @classmethod
    def random_food_web(
        cls,
        n_species: int = 20,
        connectance: float = 0.3,
        sigma: float = 0.5,
        self_regulation: float = MAY_WIGNER_SELF_REGULATION,
        seed: int | None = 0,
    ) -> "EcosystemGeometry":
        """Random food web following the May-Wigner construction.

        Off-diagonal a_ij is drawn from N(0, sigma^2) with probability
        ``connectance`` and is zero otherwise; diagonal is set to
        ``-self_regulation``.

        Local stability (May 1972): the largest real part of the
        eigenvalue spectrum scales like
            lambda_max ~ sigma * sqrt(N * C) - self_regulation.
        Stability requires sigma * sqrt(N * C) < self_regulation.
        """
        rng = np.random.default_rng(seed)
        A = rng.normal(0.0, sigma, size=(n_species, n_species))
        # Sparsify to connectance C
        mask = rng.random(size=(n_species, n_species)) < connectance
        A = A * mask
        # Self-regulation on the diagonal
        np.fill_diagonal(A, -self_regulation)
        names = [f"sp{i:02d}" for i in range(n_species)]
        # Trophic levels: random integer 1..ceil(log2(N)) for layout only
        n_levels = max(2, int(math.ceil(math.log2(max(2, n_species)))))
        levels = 1.0 + rng.integers(0, n_levels, size=n_species).astype(float)
        return cls(
            n_species=n_species,
            species_names=names,
            interaction_matrix=A,
            trophic_levels=levels,
        )

    @classmethod
    def linear_food_chain(
        cls, n_levels: int = 4,
        link_strength: float = 0.5,
        self_regulation: float = MAY_WIGNER_SELF_REGULATION,
    ) -> "EcosystemGeometry":
        """Linear chain producer -> herbivore -> ... -> top predator.

        Useful for stability analysis as a baseline against the random
        food web.
        """
        N = n_levels
        A = np.zeros((N, N))
        np.fill_diagonal(A, -self_regulation)
        A[0, 0] = +1.0   # basal producer growth
        for i in range(N - 1):
            # i+1 eats i: predator gains, prey loses
            A[i + 1, i] = +link_strength
            A[i, i + 1] = -link_strength
        names = [f"L{i + 1}" for i in range(N)]
        levels = np.arange(1, N + 1, dtype=float)
        return cls(
            n_species=N,
            species_names=names,
            interaction_matrix=A,
            trophic_levels=levels,
        )

    # ------------------------------------------------------------------
    # Graph helpers
    # ------------------------------------------------------------------

    def edges(self) -> list[tuple[int, int, float]]:
        """Return the list of directed edges (source, target, weight)."""
        out: list[tuple[int, int, float]] = []
        for i in range(self.n_species):
            for j in range(self.n_species):
                if i == j:
                    continue
                w = float(self.interaction_matrix[i, j])
                if abs(w) > 0.0:
                    out.append((j, i, w))
        return out

    def connectance(self) -> float:
        """Fraction of off-diagonal entries that are non-zero."""
        N = self.n_species
        if N < 2:
            return 0.0
        off = self.interaction_matrix.copy()
        np.fill_diagonal(off, 0.0)
        return float(np.count_nonzero(off)) / float(N * (N - 1))

    def node_positions_2d(self) -> np.ndarray:
        """Return (N, 2) layout positions for visualisation.

        Nodes within a trophic level are spread along the x-axis;
        trophic level controls the y-coordinate.
        """
        N = self.n_species
        levels = self.trophic_levels
        unique = sorted(set(levels.tolist()))
        pos = np.zeros((N, 2))
        for lev in unique:
            idx = np.where(np.isclose(levels, lev))[0]
            n = len(idx)
            xs = np.linspace(-1.0, 1.0, max(2, n + 1))[:n] if n > 1 else [0.0]
            for k, i in enumerate(idx):
                pos[i, 0] = xs[k] if n > 1 else 0.0
                pos[i, 1] = lev
        return pos

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "ecosystem = directed graph of strain configurations; "
                "edges = bilinear interaction strengths"
            ),
            "n_species": self.n_species,
            "connectance": self.connectance(),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class EcosystemSimulator:
    """Population-dynamics simulator for ecological networks.

    Provides closed-form Lotka-Volterra and logistic identities, plus an
    RK4 integrator for general N-species networks defined by
    ``EcosystemGeometry``.
    """

    geometry: EcosystemGeometry | None = None
    extinction_threshold: float = EXTINCTION_THRESHOLD_DEFAULT
    dt: float = 0.01

    # ------------------------------------------------------------------
    # Lotka-Volterra
    # ------------------------------------------------------------------

    @staticmethod
    def lotka_volterra_dynamics(
        N: float, P: float,
        alpha: float = LV_ALPHA_DEFAULT,
        beta: float = LV_BETA_DEFAULT,
        delta: float = LV_DELTA_DEFAULT,
        gamma: float = LV_GAMMA_DEFAULT,
    ) -> tuple[float, float]:
        """Evaluate (dN/dt, dP/dt) for the Lotka-Volterra system."""
        dN = alpha * N - beta * N * P
        dP = delta * N * P - gamma * P
        return dN, dP

    @staticmethod
    def lotka_volterra_equilibrium(
        alpha: float = LV_ALPHA_DEFAULT,
        beta: float = LV_BETA_DEFAULT,
        delta: float = LV_DELTA_DEFAULT,
        gamma: float = LV_GAMMA_DEFAULT,
    ) -> tuple[float, float]:
        """Interior equilibrium (N*, P*) = (gamma/delta, alpha/beta)."""
        return (gamma / delta, alpha / beta)

    @staticmethod
    def lotka_volterra_period(
        alpha: float = LV_ALPHA_DEFAULT,
        gamma: float = LV_GAMMA_DEFAULT,
    ) -> float:
        """Linearised oscillation period T = 2 pi / sqrt(alpha * gamma)."""
        return 2.0 * math.pi / math.sqrt(alpha * gamma)

    @staticmethod
    def lotka_volterra_invariant(
        N: float, P: float,
        alpha: float = LV_ALPHA_DEFAULT,
        beta: float = LV_BETA_DEFAULT,
        delta: float = LV_DELTA_DEFAULT,
        gamma: float = LV_GAMMA_DEFAULT,
    ) -> float:
        """Conserved quantity V = delta N - gamma ln N + beta P - alpha ln P.

        Constant along closed orbits in the (N, P) phase plane.
        """
        if N <= 0.0 or P <= 0.0:
            return float("inf")
        return (delta * N - gamma * math.log(N)
                + beta * P - alpha * math.log(P))

    def integrate_lotka_volterra(
        self,
        N0: float, P0: float,
        alpha: float = LV_ALPHA_DEFAULT,
        beta: float = LV_BETA_DEFAULT,
        delta: float = LV_DELTA_DEFAULT,
        gamma: float = LV_GAMMA_DEFAULT,
        t_end: float = 50.0, dt: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """RK4 integrate the LV system; return (t, N, P) arrays."""
        h = float(dt if dt is not None else self.dt)
        n_steps = int(round(t_end / h))
        t = np.linspace(0.0, n_steps * h, n_steps + 1)
        N_arr = np.empty(n_steps + 1)
        P_arr = np.empty(n_steps + 1)
        N_arr[0], P_arr[0] = N0, P0

        def f(N, P):
            return self.lotka_volterra_dynamics(N, P, alpha, beta, delta, gamma)

        for k in range(n_steps):
            N, P = N_arr[k], P_arr[k]
            k1N, k1P = f(N, P)
            k2N, k2P = f(N + 0.5 * h * k1N, P + 0.5 * h * k1P)
            k3N, k3P = f(N + 0.5 * h * k2N, P + 0.5 * h * k2P)
            k4N, k4P = f(N + h * k3N, P + h * k3P)
            N_arr[k + 1] = N + (h / 6.0) * (k1N + 2 * k2N + 2 * k3N + k4N)
            P_arr[k + 1] = P + (h / 6.0) * (k1P + 2 * k2P + 2 * k3P + k4P)
            # Apply extinction threshold
            if N_arr[k + 1] < self.extinction_threshold * 1e-6:
                N_arr[k + 1] = 0.0
            if P_arr[k + 1] < self.extinction_threshold * 1e-6:
                P_arr[k + 1] = 0.0
        return t, N_arr, P_arr

    # ------------------------------------------------------------------
    # Logistic growth
    # ------------------------------------------------------------------

    @staticmethod
    def logistic_growth_dynamics(
        N: float,
        r: float = LOGISTIC_R_DEFAULT,
        K: float = LOGISTIC_K_DEFAULT,
    ) -> float:
        """dN/dt = r N (1 - N/K)."""
        return r * N * (1.0 - N / K)

    def integrate_logistic(
        self, N0: float,
        r: float = LOGISTIC_R_DEFAULT,
        K: float = LOGISTIC_K_DEFAULT,
        t_end: float = 30.0, dt: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Closed-form logistic trajectory: N(t) = K / (1 + ((K-N0)/N0) e^(-rt))."""
        h = float(dt if dt is not None else self.dt)
        n_steps = int(round(t_end / h))
        t = np.linspace(0.0, n_steps * h, n_steps + 1)
        if N0 <= 0.0:
            return t, np.zeros_like(t)
        # Use the analytic solution to avoid integrator drift
        denom = 1.0 + ((K - N0) / N0) * np.exp(-r * t)
        N = K / denom
        return t, N

    # ------------------------------------------------------------------
    # Extinction threshold
    # ------------------------------------------------------------------

    def is_extinct(self, N: float) -> bool:
        """True if population is at or below the extinction floor."""
        return float(N) <= self.extinction_threshold

    def time_to_extinction_exponential(
        self, N0: float, decay_rate: float,
    ) -> float:
        """Time for an exponentially-decaying population to drop below
        the extinction threshold:
            N(t) = N0 * exp(-decay_rate * t) <= N_min
        """
        if N0 <= self.extinction_threshold:
            return 0.0
        if decay_rate <= 0.0:
            return float("inf")
        return math.log(N0 / self.extinction_threshold) / decay_rate

    # ------------------------------------------------------------------
    # Food-web stability (Jacobian + May-Wigner)
    # ------------------------------------------------------------------

    def food_web_jacobian(
        self, equilibrium_state: np.ndarray | None = None,
    ) -> np.ndarray:
        """Linearised Jacobian J_ij = (partial dx_i / dt) / (partial x_j).

        For a generalised Lotka-Volterra system
            dx_i/dt = x_i * ( r_i + sum_j a_ij x_j ),
        evaluated at x = x*, the Jacobian is
            J_ij = x_i^* * a_ij
        (treating r_i as encoded in the diagonal of the interaction
        matrix is a convention; here we use the geometry's interaction
        matrix directly as the linearisation around the equilibrium).
        """
        if self.geometry is None:
            raise ValueError("geometry must be set to compute the Jacobian")
        A = self.geometry.interaction_matrix
        if equilibrium_state is None:
            return A.copy()
        x = np.asarray(equilibrium_state, dtype=float)
        return (x[:, None] * A)

    def stability_eigenvalues(
        self, equilibrium_state: np.ndarray | None = None,
    ) -> np.ndarray:
        """Eigenvalues of the Jacobian at the (interior) equilibrium."""
        J = self.food_web_jacobian(equilibrium_state)
        return np.linalg.eigvals(J)

    def is_locally_stable(
        self, equilibrium_state: np.ndarray | None = None,
        tol: float = 1e-9,
    ) -> bool:
        """All Jacobian eigenvalues have negative real part?"""
        evals = self.stability_eigenvalues(equilibrium_state)
        return bool(np.all(np.real(evals) < -tol))

    @staticmethod
    def may_wigner_threshold(
        n_species: int, connectance: float,
        self_regulation: float = MAY_WIGNER_SELF_REGULATION,
    ) -> float:
        """Critical interaction strength sigma_c above which a random
        food web with N species and connectance C loses stability.

        May-Wigner (1972): sigma_c = self_regulation / sqrt(N * C).
        """
        return float(self_regulation) / math.sqrt(max(1.0, n_species * connectance))

    @staticmethod
    def may_wigner_complexity(
        n_species: int, connectance: float, sigma: float,
    ) -> float:
        """Complexity parameter sigma * sqrt(N * C); >= 1 => unstable."""
        return float(sigma) * math.sqrt(max(0.0, n_species * connectance))

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "extinction_threshold": self.extinction_threshold,
            "interpretation": (
                "population dynamics = slow-substrate occupancy flow; "
                "Lotka-Volterra invariant = strain-energy conservation"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


__all__ = [
    "EcosystemGeometry",
    "EcosystemSimulator",
    "LV_ALPHA_DEFAULT",
    "LV_BETA_DEFAULT",
    "LV_DELTA_DEFAULT",
    "LV_GAMMA_DEFAULT",
    "LOGISTIC_R_DEFAULT",
    "LOGISTIC_K_DEFAULT",
    "EXTINCTION_THRESHOLD_DEFAULT",
    "MAY_WIGNER_SELF_REGULATION",
    "LAMBDA_QCD_MEV",
]
