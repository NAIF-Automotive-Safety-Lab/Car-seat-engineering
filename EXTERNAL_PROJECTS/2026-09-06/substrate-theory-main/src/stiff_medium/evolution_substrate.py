"""Paired evolutionary-dynamics GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for population genetics
and adaptive evolution under the B3 / stiff-medium ontology.

GEOMETRY
--------
``EvolutionGeometry`` represents a population in allele-frequency space
together with a fitness landscape over phenotype/genotype coordinates.

The natural geometric objects are

    * the simplex S^{n-1} of allele frequencies p = (p_1, ..., p_n)
      with sum_i p_i = 1, p_i >= 0
    * a smooth fitness map W : R^d -> R giving the Wrightian fitness
      landscape (here a sum of Gaussian peaks plus an optional ridge)
    * a per-locus mutation rate mu (per-generation)

Under the substrate ontology the fitness landscape is the local strain
relaxation of the substrate around a phenotypic configuration: peaks are
metastable strain minima of the developmental field, valleys are
energetically forbidden corridors.

SIMULATOR
---------
``EvolutionSimulator`` implements the canonical population-genetic models:

(1) Hardy-Weinberg equilibrium (large N, no selection, no mutation):
        p^2 : 2pq : q^2  for genotype frequencies AA : Aa : aa

(2) Wright-Fisher drift in finite populations of size N (haploid, or
    diploid 2N gene copies):  next-generation count of allele A is
        n' ~ Binomial(2N, p)
    Drift causes |Delta p| ~ sqrt(p(1-p)/(2N)) per generation.

(3) Selection coefficient s, fixation probability:
    For a beneficial mutation with selection coefficient s arising
    once in a population of size N, Kimura (1962) gave

        P_fix(s, N) = (1 - exp(-2s))      / (1 - exp(-4Ns))         (diploid)
                   ~ 2s                    for small s, large N
                   ~ 1/(2N)                for s -> 0  (neutral)

    For s > 0 and 4Ns >> 1 the haploid limit P_fix -> 2s recovers
    Haldane (1927).

(4) Speciation via reproductive isolation:
    Two daughter populations diverge in allele frequencies due to drift
    and divergent selection on different fitness peaks.  Reproductive
    isolation grows monotonically with genetic distance D = sum_i |p_i^A - p_i^B|.
    A simple Bateson-Dobzhansky-Muller model gives

        I(t) = 1 - exp(-k * D(t))

    where k > 0 controls the rate at which incompatibilities accumulate.
    When I crosses 0.5 the populations are conventionally called species.

Substrate-specific predictions encoded here:
    * Drift dominates selection when 4Ns < 1  (Ohta 1973 nearly-neutral
      threshold) -- this is a topological condition on the substrate's
      local strain noise floor.
    * The number of accessible adaptive walks scales with the topology
      of the fitness landscape: a single Gaussian peak admits a unique
      gradient ascent path from any basin point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) values
# ---------------------------------------------------------------------------

# Mutation-rate scale per nucleotide per replication (eukaryotic genome)
DEFAULT_MUTATION_RATE: float = 1.0e-8

# Default population size (illustrative drift regime)
DEFAULT_POPULATION_SIZE: int = 1000

# Haldane (1927) approximation: P_fix ~ 2s  for s > 0 and 4Ns >> 1
HALDANE_FIXATION_FACTOR: float = 2.0

# Ohta (1973) nearly-neutral threshold: drift dominates when 4Ns < 1
NEARLY_NEUTRAL_THRESHOLD: float = 1.0   # in units of 4Ns

# Reproductive isolation crossover (Bateson-Dobzhansky-Muller)
SPECIATION_THRESHOLD: float = 0.5       # I = 0.5 marks species boundary


# ---------------------------------------------------------------------------
# Hardy-Weinberg helpers
# ---------------------------------------------------------------------------

def hardy_weinberg_genotype_frequencies(p: float) -> Tuple[float, float, float]:
    """Return (AA, Aa, aa) genotype frequencies under Hardy-Weinberg.

    For a biallelic locus with allele-A frequency p and allele-a frequency
    q = 1-p, the equilibrium genotype frequencies are p^2 : 2pq : q^2.
    """
    if not (0.0 <= p <= 1.0):
        raise ValueError(f"allele frequency must be in [0,1], got {p}")
    q = 1.0 - p
    return (p * p, 2.0 * p * q, q * q)


def kimura_fixation_probability(s: float, N: int) -> float:
    """Kimura (1962) diploid fixation probability for a single new mutant.

    P_fix(s, N) = (1 - exp(-2s)) / (1 - exp(-4Ns))

    For s -> 0 this approaches 1/(2N) (neutral fixation by drift).
    For 4Ns >> 1 it approaches 2s (Haldane 1927).
    """
    if N <= 0:
        raise ValueError(f"N must be positive, got {N}")
    if abs(s) < 1e-12:
        return 1.0 / (2.0 * N)
    num = 1.0 - math.exp(-2.0 * s)
    den = 1.0 - math.exp(-4.0 * N * s)
    if abs(den) < 1e-18:
        # 4Ns ~ 0 limit -> neutral
        return 1.0 / (2.0 * N)
    return num / den


def haldane_fixation_probability(s: float) -> float:
    """Haldane (1927) limiting form: P_fix ~ 2s for small s > 0."""
    if s <= 0.0:
        return 0.0
    return HALDANE_FIXATION_FACTOR * s


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class EvolutionGeometry:
    """Population on a fitness landscape in allele-frequency space.

    Parameters
    ----------
    n_alleles : int
        Number of allelic types (>=2).  The simplex of allele frequencies
        is S^{n-1}.
    allele_frequencies : np.ndarray, shape (n_alleles,)
        Initial allele frequencies, sum to 1.  Defaults to uniform.
    mutation_rate : float
        Per-locus per-generation mutation rate.
    landscape_peaks : sequence of (centre, height, width) tuples
        Each peak is a 2D Gaussian on the (x,y) phenotype plane:
            W(x,y) += height * exp(- ||(x,y) - centre||^2 / (2 width^2))
    """

    n_alleles: int = 2
    allele_frequencies: Optional[np.ndarray] = None
    mutation_rate: float = DEFAULT_MUTATION_RATE
    landscape_peaks: Sequence[Tuple[Tuple[float, float], float, float]] = field(
        default_factory=lambda: [
            ((0.0, 0.0), 1.0, 1.0),         # central modest peak
            ((3.0, 2.5), 1.5, 0.8),         # high adaptive peak
            ((-2.5, 2.0), 0.9, 0.9),        # secondary peak
            ((2.0, -2.5), 0.7, 1.0),        # low side peak
        ]
    )

    def __post_init__(self) -> None:
        if self.n_alleles < 2:
            raise ValueError(f"n_alleles must be >= 2, got {self.n_alleles}")
        if self.mutation_rate < 0.0:
            raise ValueError(f"mutation_rate must be >= 0, got {self.mutation_rate}")
        if self.allele_frequencies is None:
            self.allele_frequencies = np.full(
                self.n_alleles, 1.0 / self.n_alleles, dtype=float,
            )
        else:
            self.allele_frequencies = np.asarray(
                self.allele_frequencies, dtype=float,
            )
            if self.allele_frequencies.shape != (self.n_alleles,):
                raise ValueError(
                    f"allele_frequencies must have shape ({self.n_alleles},)"
                )
            if not np.isclose(self.allele_frequencies.sum(), 1.0, atol=1e-9):
                raise ValueError("allele_frequencies must sum to 1")
            if (self.allele_frequencies < -1e-12).any():
                raise ValueError("allele_frequencies must be non-negative")

    # --------------------------- landscape API -------------------------- #

    def fitness(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Evaluate the Wrightian fitness landscape W(x,y).

        Parameters
        ----------
        x, y : np.ndarray
            Coordinates (broadcastable).

        Returns
        -------
        np.ndarray
            Fitness values, same shape as broadcast(x,y).
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        W = np.zeros(np.broadcast_shapes(x.shape, y.shape), dtype=float)
        for (cx, cy), h, w in self.landscape_peaks:
            r2 = (x - cx) ** 2 + (y - cy) ** 2
            W += h * np.exp(-0.5 * r2 / (w * w))
        return W

    def fitness_gradient(self, x: float, y: float) -> Tuple[float, float]:
        """Return (dW/dx, dW/dy) at a single point (x,y)."""
        gx = 0.0
        gy = 0.0
        for (cx, cy), h, w in self.landscape_peaks:
            dx = x - cx
            dy = y - cy
            r2 = dx * dx + dy * dy
            g = h * math.exp(-0.5 * r2 / (w * w)) / (w * w)
            gx += -dx * g
            gy += -dy * g
        return (gx, gy)

    def landscape_grid(
        self,
        x_range: Tuple[float, float] = (-5.0, 5.0),
        y_range: Tuple[float, float] = (-5.0, 5.0),
        n: int = 80,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (X, Y, W) meshgrid of the fitness landscape."""
        xs = np.linspace(x_range[0], x_range[1], n)
        ys = np.linspace(y_range[0], y_range[1], n)
        X, Y = np.meshgrid(xs, ys)
        W = self.fitness(X, Y)
        return X, Y, W

    def n_landscape_peaks(self) -> int:
        """Return the number of Gaussian peaks defining the landscape."""
        return len(self.landscape_peaks)

    def simplex_dimension(self) -> int:
        """Dimension of the allele-frequency simplex S^{n-1}."""
        return self.n_alleles - 1


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class EvolutionSimulator:
    """Canonical population-genetic simulator on an EvolutionGeometry.

    Parameters
    ----------
    geometry : EvolutionGeometry
        The fitness landscape and initial allele frequencies.
    population_size : int
        Diploid population size N (so 2N gene copies per generation).
    selection_coefficient : float
        Selection coefficient s on allele A (relative fitness 1+s vs 1).
    """

    geometry: EvolutionGeometry = field(default_factory=EvolutionGeometry)
    population_size: int = DEFAULT_POPULATION_SIZE
    selection_coefficient: float = 0.0

    def __post_init__(self) -> None:
        if self.population_size <= 0:
            raise ValueError(
                f"population_size must be > 0, got {self.population_size}"
            )

    # ------------------------ Hardy-Weinberg ---------------------------- #

    def hardy_weinberg(self, p: float) -> Tuple[float, float, float]:
        """Return (AA, Aa, aa) genotype frequencies at HW equilibrium."""
        return hardy_weinberg_genotype_frequencies(p)

    def is_in_hw_equilibrium(
        self,
        genotype_freqs: Tuple[float, float, float],
        atol: float = 1e-3,
    ) -> bool:
        """Test whether observed (AA, Aa, aa) frequencies sit at HW.

        Recovers p = AA + Aa/2 from the genotype counts and compares the
        expected (p^2, 2pq, q^2) profile within tolerance `atol`.
        """
        f_aa, f_ab, f_bb = genotype_freqs
        total = f_aa + f_ab + f_bb
        if not math.isclose(total, 1.0, abs_tol=1e-6):
            raise ValueError(
                f"genotype frequencies must sum to 1, got {total}"
            )
        p = f_aa + 0.5 * f_ab
        expected = hardy_weinberg_genotype_frequencies(p)
        return all(
            math.isclose(o, e, abs_tol=atol)
            for o, e in zip(genotype_freqs, expected)
        )

    # ------------------------ Wright-Fisher drift ----------------------- #

    def wright_fisher_step(
        self,
        p: float,
        rng: np.random.Generator,
    ) -> float:
        """Single Wright-Fisher generation: n' ~ Binomial(2N, p)/(2N).

        Optionally biased by selection_coefficient s on allele A:
        sampling probability becomes p(1+s)/((1+s)p + (1-p)).
        """
        if self.selection_coefficient != 0.0:
            s = self.selection_coefficient
            wbar = (1.0 + s) * p + (1.0 - p)
            p_eff = (1.0 + s) * p / wbar
        else:
            p_eff = p
        # mutation: A <-> a at rate mu, symmetric per-locus
        mu = self.geometry.mutation_rate
        if mu > 0.0:
            p_eff = (1.0 - mu) * p_eff + mu * (1.0 - p_eff)
        n_copies = 2 * self.population_size
        n_A = rng.binomial(n_copies, p_eff)
        return n_A / n_copies

    def wright_fisher_trajectory(
        self,
        p0: float,
        n_generations: int,
        seed: Optional[int] = None,
    ) -> np.ndarray:
        """Return allele-frequency trajectory over `n_generations` generations."""
        rng = np.random.default_rng(seed)
        traj = np.empty(n_generations + 1, dtype=float)
        traj[0] = p0
        p = p0
        for t in range(1, n_generations + 1):
            p = self.wright_fisher_step(p, rng)
            traj[t] = p
        return traj

    def drift_variance_per_generation(self, p: float) -> float:
        """Wright-Fisher per-generation variance of p:  p(1-p)/(2N)."""
        return p * (1.0 - p) / (2.0 * self.population_size)

    # ------------------------ Selection / fixation ---------------------- #

    def fixation_probability(
        self,
        s: Optional[float] = None,
        N: Optional[int] = None,
    ) -> float:
        """Kimura (1962) fixation probability of a single new mutant."""
        s_use = self.selection_coefficient if s is None else s
        N_use = self.population_size if N is None else N
        return kimura_fixation_probability(s_use, N_use)

    def haldane_approximation(self, s: Optional[float] = None) -> float:
        """Haldane limit P_fix ~ 2s for small positive s."""
        s_use = self.selection_coefficient if s is None else s
        return haldane_fixation_probability(s_use)

    def empirical_fixation_probability(
        self,
        s: float,
        N: int,
        n_trials: int = 2000,
        max_generations: int = 20000,
        seed: Optional[int] = None,
    ) -> float:
        """Estimate fixation probability of a single mutant by simulation.

        Each trial seeds one mutant copy in 2N gene copies (p0 = 1/(2N))
        and Wright-Fisher-evolves until p hits 0 or 1.
        """
        rng = np.random.default_rng(seed)
        original_s = self.selection_coefficient
        original_N = self.population_size
        self.selection_coefficient = s
        self.population_size = N
        try:
            n_fixed = 0
            p0 = 1.0 / (2.0 * N)
            for _ in range(n_trials):
                p = p0
                for _ in range(max_generations):
                    p = self.wright_fisher_step(p, rng)
                    if p <= 0.0 or p >= 1.0:
                        break
                if p >= 1.0:
                    n_fixed += 1
            return n_fixed / n_trials
        finally:
            self.selection_coefficient = original_s
            self.population_size = original_N

    def drift_dominates(
        self, s: Optional[float] = None, N: Optional[int] = None,
    ) -> bool:
        """Ohta (1973) nearly-neutral test: drift dominates when 4Ns < 1."""
        s_use = self.selection_coefficient if s is None else s
        N_use = self.population_size if N is None else N
        return 4.0 * N_use * abs(s_use) < NEARLY_NEUTRAL_THRESHOLD

    # ------------------------ Adaptive walks ---------------------------- #

    def adaptive_walk(
        self,
        start: Tuple[float, float],
        step_size: float = 0.05,
        n_steps: int = 400,
        noise: float = 0.02,
        seed: Optional[int] = None,
    ) -> np.ndarray:
        """Hill-climbing trajectory with stochastic noise on the landscape.

        Returns
        -------
        np.ndarray, shape (n_steps+1, 2)
            Sequence of (x, y) positions.
        """
        rng = np.random.default_rng(seed)
        path = np.empty((n_steps + 1, 2), dtype=float)
        path[0] = start
        x, y = start
        for t in range(1, n_steps + 1):
            gx, gy = self.geometry.fitness_gradient(x, y)
            gnorm = math.hypot(gx, gy)
            if gnorm > 0.0:
                gx /= gnorm
                gy /= gnorm
            x += step_size * gx + noise * rng.standard_normal()
            y += step_size * gy + noise * rng.standard_normal()
            path[t] = (x, y)
        return path

    # ------------------------ Speciation -------------------------------- #

    def genetic_distance(self, p_A: np.ndarray, p_B: np.ndarray) -> float:
        """L1 distance between two allele-frequency vectors."""
        return float(np.sum(np.abs(np.asarray(p_A) - np.asarray(p_B))))

    def reproductive_isolation(
        self, distance: float, k: float = 4.0,
    ) -> float:
        """Bateson-Dobzhansky-Muller isolation: I = 1 - exp(-k * D).

        I = 0   -> fully interfertile
        I = 1   -> fully isolated (distinct species)
        Crosses 0.5 when D = ln(2) / k.
        """
        if distance < 0.0:
            raise ValueError(f"distance must be >= 0, got {distance}")
        return 1.0 - math.exp(-k * distance)

    def is_speciated(self, isolation: float) -> bool:
        """Return True if reproductive isolation has crossed the species boundary."""
        return isolation >= SPECIATION_THRESHOLD

    def divergence_simulation(
        self,
        n_generations: int = 1500,
        n_loci: int = 50,
        s_pop_A: float = +0.005,
        s_pop_B: float = -0.005,
        k: float = 4.0,
        seed: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Simulate two daughter populations diverging under drift + divergent selection.

        Both populations start from p = 0.5 at every locus.  Population A
        feels selection coefficient ``s_pop_A`` on allele A, population B
        feels ``s_pop_B``.  Mutation and Wright-Fisher drift act on top.

        Returns
        -------
        gens : np.ndarray, shape (n_generations+1,)
        p_A  : np.ndarray, shape (n_generations+1, n_loci)
        p_B  : np.ndarray, shape (n_generations+1, n_loci)
        I    : np.ndarray, shape (n_generations+1,)   reproductive isolation
        """
        rng = np.random.default_rng(seed)
        N = self.population_size
        mu = self.geometry.mutation_rate
        n_copies = 2 * N
        gens = np.arange(n_generations + 1)
        p_A = np.empty((n_generations + 1, n_loci), dtype=float)
        p_B = np.empty((n_generations + 1, n_loci), dtype=float)
        p_A[0] = 0.5
        p_B[0] = 0.5
        for t in range(1, n_generations + 1):
            for which, p_array, s in (
                (p_A, p_A[t - 1], s_pop_A),
                (p_B, p_B[t - 1], s_pop_B),
            ):
                p_eff = p_array.copy()
                if s != 0.0:
                    wbar = (1.0 + s) * p_eff + (1.0 - p_eff)
                    p_eff = (1.0 + s) * p_eff / wbar
                if mu > 0.0:
                    p_eff = (1.0 - mu) * p_eff + mu * (1.0 - p_eff)
                p_eff = np.clip(p_eff, 0.0, 1.0)
                n_succ = rng.binomial(n_copies, p_eff)
                which[t] = n_succ / n_copies
        # reproductive isolation as a function of mean L1 distance
        D = np.mean(np.abs(p_A - p_B), axis=1) * n_loci  # total L1 over loci
        I = 1.0 - np.exp(-k * D / n_loci)                # normalise per locus
        return gens, p_A, p_B, I
