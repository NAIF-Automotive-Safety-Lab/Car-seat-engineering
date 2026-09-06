"""Tests for the paired evolutionary-dynamics GEOMETRY + SIMULATOR module.

Covers the three required claims:
    (a) Hardy-Weinberg equilibrium for unselected loci
    (b) Fixation probability ~ 2s for small s, large N (Haldane / Kimura)
    (c) Drift dominates selection when 4Ns < 1 (Ohta nearly-neutral regime)

Plus geometric and population-genetic sanity:
    * allele frequencies live on the simplex (sum = 1, non-negative)
    * fitness landscape is a sum of Gaussian peaks (multiple maxima)
    * Wright-Fisher per-generation variance = p(1-p)/(2N)
    * adaptive walks ascend the gradient
    * reproductive isolation grows monotonically with genetic distance
    * speciation crossing at I = 0.5
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.evolution_substrate import (
    DEFAULT_MUTATION_RATE,
    DEFAULT_POPULATION_SIZE,
    HALDANE_FIXATION_FACTOR,
    NEARLY_NEUTRAL_THRESHOLD,
    SPECIATION_THRESHOLD,
    EvolutionGeometry,
    EvolutionSimulator,
    haldane_fixation_probability,
    hardy_weinberg_genotype_frequencies,
    kimura_fixation_probability,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom_default() -> EvolutionGeometry:
    return EvolutionGeometry()


@pytest.fixture
def sim_neutral() -> EvolutionSimulator:
    """No selection, finite population for drift studies."""
    return EvolutionSimulator(
        geometry=EvolutionGeometry(mutation_rate=0.0),
        population_size=200,
        selection_coefficient=0.0,
    )


@pytest.fixture
def sim_selected() -> EvolutionSimulator:
    """Beneficial mutation, large population (Haldane regime)."""
    return EvolutionSimulator(
        geometry=EvolutionGeometry(mutation_rate=0.0),
        population_size=10000,
        selection_coefficient=0.01,
    )


# --------------------------------------------------------------------- #
# Geometry sanity                                                       #
# --------------------------------------------------------------------- #

def test_geometry_default_allele_frequencies_uniform(geom_default):
    """Default allele frequencies form the uniform distribution on S^{n-1}."""
    p = geom_default.allele_frequencies
    assert p.shape == (2,)
    assert np.allclose(p, 0.5)
    assert math.isclose(float(p.sum()), 1.0)


def test_geometry_simplex_dimension():
    """The allele-frequency simplex S^{n-1} has dimension n-1."""
    g3 = EvolutionGeometry(n_alleles=3)
    g4 = EvolutionGeometry(n_alleles=4)
    assert g3.simplex_dimension() == 2
    assert g4.simplex_dimension() == 3


def test_geometry_rejects_invalid_frequencies():
    """allele_frequencies must sum to 1 and be non-negative."""
    with pytest.raises(ValueError):
        EvolutionGeometry(n_alleles=2, allele_frequencies=np.array([0.6, 0.6]))
    with pytest.raises(ValueError):
        EvolutionGeometry(n_alleles=2, allele_frequencies=np.array([1.5, -0.5]))


def test_geometry_rejects_negative_mutation_rate():
    with pytest.raises(ValueError):
        EvolutionGeometry(mutation_rate=-1e-3)


def test_geometry_landscape_has_multiple_peaks(geom_default):
    """The default landscape carries >= 2 Gaussian peaks (rugged)."""
    assert geom_default.n_landscape_peaks() >= 2


def test_geometry_fitness_is_positive_at_peak(geom_default):
    """Fitness at the highest peak's centre exceeds fitness far away."""
    peaks = geom_default.landscape_peaks
    centre = max(peaks, key=lambda p: p[1])[0]
    W_centre = float(geom_default.fitness(np.array(centre[0]), np.array(centre[1])))
    W_far = float(geom_default.fitness(np.array(20.0), np.array(20.0)))
    assert W_centre > W_far


def test_geometry_landscape_grid_shape(geom_default):
    X, Y, W = geom_default.landscape_grid(n=32)
    assert X.shape == (32, 32)
    assert Y.shape == (32, 32)
    assert W.shape == (32, 32)
    assert np.all(W >= 0.0)


# --------------------------------------------------------------------- #
# (a) Hardy-Weinberg equilibrium                                        #
# --------------------------------------------------------------------- #

def test_hw_genotype_frequencies_p_half():
    """p = 0.5 -> AA : Aa : aa = 0.25 : 0.50 : 0.25."""
    f = hardy_weinberg_genotype_frequencies(0.5)
    assert f == pytest.approx((0.25, 0.50, 0.25))
    assert sum(f) == pytest.approx(1.0)


def test_hw_genotype_frequencies_extremes():
    """p = 0 -> all aa;  p = 1 -> all AA."""
    assert hardy_weinberg_genotype_frequencies(0.0) == pytest.approx((0.0, 0.0, 1.0))
    assert hardy_weinberg_genotype_frequencies(1.0) == pytest.approx((1.0, 0.0, 0.0))


def test_hw_genotype_frequencies_sum_to_one():
    """For arbitrary p in (0,1), p^2 + 2pq + q^2 = (p+q)^2 = 1."""
    for p in [0.1, 0.27, 0.42, 0.63, 0.85]:
        f = hardy_weinberg_genotype_frequencies(p)
        assert sum(f) == pytest.approx(1.0, abs=1e-12)


def test_unselected_locus_relaxes_to_hw(sim_neutral):
    """A randomly mating, unselected, non-mutating locus sits at HW.

    We check the equilibrium predicate directly: the textbook HW
    profile (p^2, 2pq, q^2) satisfies the simulator's HW test.
    """
    for p in [0.2, 0.5, 0.7]:
        gen_freqs = hardy_weinberg_genotype_frequencies(p)
        assert sim_neutral.is_in_hw_equilibrium(gen_freqs, atol=1e-9)


def test_non_hw_genotype_frequencies_detected(sim_neutral):
    """Out-of-equilibrium genotype frequencies are flagged as non-HW.

    p = 0.5 -> HW expects (0.25, 0.50, 0.25).  An excess homozygote
    profile (0.40, 0.20, 0.40) is far from HW and must fail the test.
    """
    bad = (0.40, 0.20, 0.40)
    assert sim_neutral.is_in_hw_equilibrium(bad, atol=1e-3) is False


def test_hw_population_under_pure_drift_stays_near_starting_p():
    """Under pure drift (no selection, no mutation), the average allele
    frequency over many independent populations is conserved E[p_t] = p_0."""
    sim = EvolutionSimulator(
        geometry=EvolutionGeometry(mutation_rate=0.0),
        population_size=500,
        selection_coefficient=0.0,
    )
    rng = np.random.default_rng(0)
    p0 = 0.5
    n_pops = 200
    p_final = []
    for k in range(n_pops):
        p = p0
        for _ in range(20):
            p = sim.wright_fisher_step(p, rng)
        p_final.append(p)
    mean_p = float(np.mean(p_final))
    # E[p_t] = p_0 to within sampling error ~ sigma / sqrt(n_pops)
    sigma = math.sqrt(p0 * (1.0 - p0) / (2.0 * sim.population_size))
    se = sigma * math.sqrt(20) / math.sqrt(n_pops)
    assert abs(mean_p - p0) < 5.0 * se + 0.05


# --------------------------------------------------------------------- #
# (b) Fixation probability ~ 2s for small s, large N                    #
# --------------------------------------------------------------------- #

def test_kimura_fixation_neutral_limit():
    """At s = 0, Kimura's formula collapses to 1/(2N)."""
    N = 1000
    p_fix = kimura_fixation_probability(0.0, N)
    assert p_fix == pytest.approx(1.0 / (2.0 * N), rel=1e-9)


def test_kimura_fixation_large_4Ns_recovers_haldane():
    """For small s and 4Ns >> 1, P_fix -> 2s (Haldane 1927).

    The exact Kimura formula gives (1 - exp(-2s)) / (1 - exp(-4Ns)) and
    expands as 2s * (1 - s + ...) for small s.  We require both that the
    Haldane shortcut returns exactly 2s (it does, by definition) and that
    Kimura agrees with 2s at the leading-order O(s) expected accuracy.
    """
    s = 1.0e-4   # small enough that the O(s) correction is negligible
    N = 1_000_000
    assert 4.0 * N * s > 100.0
    p_fix_kim = kimura_fixation_probability(s, N)
    p_fix_hal = haldane_fixation_probability(s)
    # leading-order agreement: 1 - exp(-2s) = 2s - 2 s^2 + ... ; relative
    # error is ~ s, i.e. 1e-4 here, so a 1e-3 tolerance is comfortable
    assert p_fix_kim == pytest.approx(p_fix_hal, rel=1.0e-3)
    assert p_fix_hal == pytest.approx(HALDANE_FIXATION_FACTOR * s, rel=1e-12)


def test_kimura_fixation_intermediate_4Ns():
    """For 4Ns ~ a few, Kimura sits between neutral 1/(2N) and Haldane 2s."""
    s = 0.001
    N = 1000
    p_fix = kimura_fixation_probability(s, N)
    p_neutral = 1.0 / (2.0 * N)
    p_haldane = HALDANE_FIXATION_FACTOR * s
    assert p_neutral < p_fix < p_haldane * 1.05


def test_simulator_fixation_probability_consistent(sim_selected):
    """The simulator's analytic call agrees with the standalone Kimura formula."""
    p_sim = sim_selected.fixation_probability()
    p_kim = kimura_fixation_probability(
        sim_selected.selection_coefficient,
        sim_selected.population_size,
    )
    assert p_sim == pytest.approx(p_kim, rel=1e-12)


def test_simulator_haldane_approximation_two_s(sim_selected):
    """The simulator's Haldane shortcut returns exactly 2s."""
    p = sim_selected.haldane_approximation()
    assert p == pytest.approx(2.0 * sim_selected.selection_coefficient, rel=1e-12)


def test_empirical_fixation_probability_matches_2s_for_strong_selection():
    """Stochastic Wright-Fisher trials reproduce the Haldane 2s prediction
    for a strongly-selected mutant in a large population.

    We use s = 0.05 (4Ns = 200 >> 1) and average over 4000 trials, so the
    expected count of fixations is ~ 4000 * 2 * 0.05 = 400, with binomial
    standard deviation ~ sqrt(400 * 0.9) ~ 19.  Our 50% relative tolerance
    gives a comfortable >= 5-sigma envelope.
    """
    sim = EvolutionSimulator(
        geometry=EvolutionGeometry(mutation_rate=0.0),
        population_size=1000,
        selection_coefficient=0.05,
    )
    p_est = sim.empirical_fixation_probability(
        s=0.05, N=1000, n_trials=4000, max_generations=10000, seed=7,
    )
    expected = 2.0 * 0.05
    assert p_est == pytest.approx(expected, rel=0.50)


# --------------------------------------------------------------------- #
# (c) Drift dominates at small N (Ohta 1973: 4Ns < 1)                   #
# --------------------------------------------------------------------- #

def test_drift_dominates_in_small_population():
    """Tiny N + tiny s -> 4Ns < 1 -> drift wins."""
    sim = EvolutionSimulator(
        geometry=EvolutionGeometry(),
        population_size=20,
        selection_coefficient=0.001,
    )
    assert 4.0 * sim.population_size * abs(sim.selection_coefficient) < 1.0
    assert sim.drift_dominates() is True


def test_selection_dominates_in_large_population():
    """Large N + non-trivial s -> 4Ns >> 1 -> selection wins."""
    sim = EvolutionSimulator(
        geometry=EvolutionGeometry(),
        population_size=100000,
        selection_coefficient=0.01,
    )
    assert 4.0 * sim.population_size * abs(sim.selection_coefficient) > 1.0
    assert sim.drift_dominates() is False


def test_drift_variance_per_generation_formula(sim_neutral):
    """Wright-Fisher per-generation variance is p(1-p)/(2N)."""
    p = 0.4
    expected = p * (1.0 - p) / (2.0 * sim_neutral.population_size)
    assert sim_neutral.drift_variance_per_generation(p) == pytest.approx(expected)


def test_drift_variance_empirical_matches_theory():
    """Sampled variance of (Delta p)^2 over many one-step trajectories
    matches p(1-p)/(2N).  Uses N=200, p0=0.5."""
    sim = EvolutionSimulator(
        geometry=EvolutionGeometry(mutation_rate=0.0),
        population_size=200,
        selection_coefficient=0.0,
    )
    rng = np.random.default_rng(123)
    p0 = 0.5
    n_trials = 5000
    deltas = np.empty(n_trials)
    for k in range(n_trials):
        deltas[k] = sim.wright_fisher_step(p0, rng) - p0
    emp_var = float(np.var(deltas))
    th_var = sim.drift_variance_per_generation(p0)
    assert emp_var == pytest.approx(th_var, rel=0.10)


def test_drift_dominates_at_small_N_traj():
    """At small N (= 30), an unselected population's allele frequency
    wanders by O(1) within ~ 4N generations.  We require the trajectory
    to leave the 0.4-0.6 starting band by then."""
    sim = EvolutionSimulator(
        geometry=EvolutionGeometry(mutation_rate=0.0),
        population_size=30,
        selection_coefficient=0.0,
    )
    traj = sim.wright_fisher_trajectory(p0=0.5, n_generations=120, seed=11)
    # std of accumulated drift over 4N generations is ~ sqrt(4N * p(1-p)/(2N))
    # = sqrt(2 * 0.25) ~ 0.71  -> the trajectory must wander widely
    assert (traj.max() - traj.min()) > 0.2


# --------------------------------------------------------------------- #
# Adaptive walks                                                        #
# --------------------------------------------------------------------- #

def test_adaptive_walk_ascends_landscape(sim_neutral):
    """Hill-climbing (with weak noise) ends at higher fitness than it started."""
    start = (-3.0, -3.0)
    path = sim_neutral.adaptive_walk(
        start=start, step_size=0.05, n_steps=300, noise=0.0, seed=42,
    )
    W_start = float(sim_neutral.geometry.fitness(np.array(start[0]), np.array(start[1])))
    W_end = float(sim_neutral.geometry.fitness(
        np.array(path[-1, 0]), np.array(path[-1, 1])
    ))
    assert W_end > W_start


def test_adaptive_walk_returns_correct_shape(sim_neutral):
    path = sim_neutral.adaptive_walk(start=(0.0, 0.0), n_steps=50, seed=1)
    assert path.shape == (51, 2)


# --------------------------------------------------------------------- #
# Speciation                                                            #
# --------------------------------------------------------------------- #

def test_genetic_distance_zero_for_identical_populations(sim_neutral):
    p = np.array([0.4, 0.3, 0.3])
    assert sim_neutral.genetic_distance(p, p) == pytest.approx(0.0)


def test_genetic_distance_l1(sim_neutral):
    p_A = np.array([1.0, 0.0])
    p_B = np.array([0.0, 1.0])
    assert sim_neutral.genetic_distance(p_A, p_B) == pytest.approx(2.0)


def test_reproductive_isolation_monotonic(sim_neutral):
    """I(D) = 1 - exp(-k D) is monotonically increasing in D."""
    Ds = np.linspace(0.0, 5.0, 50)
    Is = np.array([sim_neutral.reproductive_isolation(d) for d in Ds])
    assert np.all(np.diff(Is) >= -1e-12)
    assert Is[0] == pytest.approx(0.0)
    assert Is[-1] > 0.99


def test_speciation_threshold_crossing(sim_neutral):
    """I crosses 0.5 at D = ln(2)/k (k = 4 default)."""
    k = 4.0
    D_cross = math.log(2.0) / k
    I_at = sim_neutral.reproductive_isolation(D_cross, k=k)
    assert I_at == pytest.approx(SPECIATION_THRESHOLD, abs=1e-9)
    assert sim_neutral.is_speciated(I_at + 1e-6)
    assert not sim_neutral.is_speciated(I_at - 1e-6)


def test_divergence_simulation_runs_and_increases_isolation():
    """Two diverging populations under opposite selection accumulate
    reproductive isolation over time."""
    sim = EvolutionSimulator(
        geometry=EvolutionGeometry(mutation_rate=1e-4),
        population_size=200,
    )
    gens, p_A, p_B, I = sim.divergence_simulation(
        n_generations=600, n_loci=20,
        s_pop_A=+0.02, s_pop_B=-0.02,
        seed=99,
    )
    assert gens.shape == (601,)
    assert p_A.shape == (601, 20)
    assert p_B.shape == (601, 20)
    assert I.shape == (601,)
    assert I[0] == pytest.approx(0.0, abs=0.05)
    # divergent selection -> isolation grows
    assert I[-1] > I[0] + 0.1
