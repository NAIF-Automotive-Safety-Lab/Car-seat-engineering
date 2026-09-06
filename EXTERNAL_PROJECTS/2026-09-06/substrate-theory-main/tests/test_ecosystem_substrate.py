"""Tests for the paired ecosystem-dynamics GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) Lotka-Volterra small-amplitude period T = 2 pi / sqrt(alpha gamma)
    (b) Lotka-Volterra interior equilibrium (gamma/delta, alpha/beta)
    (c) Logistic growth saturates at carrying capacity K
    (d) Extinction threshold (Allee floor) drives sub-floor populations to 0
    (e) May-Wigner instability for large random food webs
        (sigma * sqrt(N C) > 1 implies a positive-real-part eigenvalue)
    (f) Linear food chain has all-negative eigenvalues at the chosen
        link strength (locally stable equilibrium baseline).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.ecosystem_substrate import (
    EXTINCTION_THRESHOLD_DEFAULT,
    LAMBDA_QCD_MEV,
    LOGISTIC_K_DEFAULT,
    LOGISTIC_R_DEFAULT,
    LV_ALPHA_DEFAULT,
    LV_BETA_DEFAULT,
    LV_DELTA_DEFAULT,
    LV_GAMMA_DEFAULT,
    MAY_WIGNER_SELF_REGULATION,
    EcosystemGeometry,
    EcosystemSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def predator_prey_geom() -> EcosystemGeometry:
    return EcosystemGeometry.predator_prey()


@pytest.fixture
def sim() -> EcosystemSimulator:
    return EcosystemSimulator()


# --------------------------------------------------------------------- #
# (a) Lotka-Volterra oscillation period                                 #
# --------------------------------------------------------------------- #

def test_lv_period_formula(sim):
    """Linearised period must equal 2 pi / sqrt(alpha gamma)."""
    alpha = 1.0
    gamma = 1.5
    T = sim.lotka_volterra_period(alpha=alpha, gamma=gamma)
    expected = 2.0 * math.pi / math.sqrt(alpha * gamma)
    assert abs(T - expected) < 1e-9


def test_lv_period_scales_with_root_alpha_gamma(sim):
    """Halving alpha*gamma should multiply the period by sqrt(2)."""
    T1 = sim.lotka_volterra_period(alpha=1.0, gamma=1.0)
    T2 = sim.lotka_volterra_period(alpha=0.5, gamma=1.0)
    assert abs(T2 / T1 - math.sqrt(2.0)) < 1e-9


def test_lv_period_matches_simulation(sim):
    """RK4-integrated trajectory should oscillate with the predicted period."""
    alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5
    N_eq = gamma / delta
    P_eq = alpha / beta
    # Start CLOSE to equilibrium so the small-amplitude period is accurate
    N0 = N_eq * 1.02
    P0 = P_eq * 1.00
    t, N, P = sim.integrate_lotka_volterra(
        N0=N0, P0=P0,
        alpha=alpha, beta=beta, delta=delta, gamma=gamma,
        t_end=80.0, dt=0.001,
    )
    # Find peaks in N to estimate the period
    dN = np.diff(N)
    sign_changes = np.where((dN[:-1] > 0) & (dN[1:] <= 0))[0] + 1
    assert len(sign_changes) >= 3
    peak_times = t[sign_changes]
    measured_T = float(np.mean(np.diff(peak_times)))
    expected_T = 2.0 * math.pi / math.sqrt(alpha * gamma)
    rel_err = abs(measured_T - expected_T) / expected_T
    assert rel_err < 0.05, (
        f"measured T = {measured_T:.3f} vs expected {expected_T:.3f}"
    )


# --------------------------------------------------------------------- #
# (b) Lotka-Volterra interior equilibrium                               #
# --------------------------------------------------------------------- #

def test_lv_equilibrium_formula(sim):
    """Equilibrium = (gamma/delta, alpha/beta)."""
    alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5
    N_eq, P_eq = sim.lotka_volterra_equilibrium(alpha, beta, delta, gamma)
    assert abs(N_eq - gamma / delta) < 1e-12
    assert abs(P_eq - alpha / beta) < 1e-12


def test_lv_equilibrium_is_fixed_point(sim):
    """Derivatives vanish at the equilibrium."""
    alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5
    N_eq, P_eq = sim.lotka_volterra_equilibrium(alpha, beta, delta, gamma)
    dN, dP = sim.lotka_volterra_dynamics(N_eq, P_eq,
                                         alpha, beta, delta, gamma)
    assert abs(dN) < 1e-12
    assert abs(dP) < 1e-12


def test_lv_invariant_conserved_by_rk4(sim):
    """V(N, P) should be conserved along an LV orbit (within RK4 drift)."""
    alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5
    N0, P0 = 30.0, 12.0
    t, N, P = sim.integrate_lotka_volterra(
        N0=N0, P0=P0,
        alpha=alpha, beta=beta, delta=delta, gamma=gamma,
        t_end=20.0, dt=0.001,
    )
    V0 = sim.lotka_volterra_invariant(N[0], P[0],
                                       alpha, beta, delta, gamma)
    V_traj = np.array([
        sim.lotka_volterra_invariant(N[k], P[k],
                                      alpha, beta, delta, gamma)
        for k in range(0, len(t), 200)
    ])
    drift = float(np.max(np.abs(V_traj - V0)) / abs(V0))
    assert drift < 0.01, f"LV invariant drift = {drift:.4f}"


# --------------------------------------------------------------------- #
# (c) Logistic growth saturates at K                                    #
# --------------------------------------------------------------------- #

def test_logistic_saturates_at_K(sim):
    """N(t) -> K as t -> infty for r > 0, N0 > 0."""
    K = 1000.0
    r = 0.5
    t, N = sim.integrate_logistic(N0=10.0, r=r, K=K, t_end=60.0, dt=0.05)
    assert abs(N[-1] - K) / K < 1e-3


def test_logistic_growth_rate_at_zero(sim):
    """dN/dt = r N at small N (linearisation)."""
    r = 0.5
    K = 1000.0
    Nsmall = 1.0
    rate = sim.logistic_growth_dynamics(Nsmall, r=r, K=K)
    expected = r * Nsmall * (1.0 - Nsmall / K)
    assert abs(rate - expected) < 1e-12


def test_logistic_dN_dt_zero_at_K(sim):
    """No growth at carrying capacity."""
    rate = sim.logistic_growth_dynamics(LOGISTIC_K_DEFAULT,
                                        r=LOGISTIC_R_DEFAULT,
                                        K=LOGISTIC_K_DEFAULT)
    assert abs(rate) < 1e-9


# --------------------------------------------------------------------- #
# (d) Extinction threshold                                              #
# --------------------------------------------------------------------- #

def test_population_below_threshold_is_extinct(sim):
    """Below the floor the species is declared extinct."""
    assert sim.is_extinct(0.5) is True


def test_population_above_threshold_survives(sim):
    """Above the floor the species survives."""
    assert sim.is_extinct(10.0) is False


def test_time_to_extinction_exponential(sim):
    """An exponentially-decaying population reaches the floor in
    t = ln(N0 / N_min) / decay_rate."""
    decay = 0.1
    N0 = 100.0
    t_ext = sim.time_to_extinction_exponential(N0=N0, decay_rate=decay)
    expected = math.log(N0 / sim.extinction_threshold) / decay
    assert abs(t_ext - expected) < 1e-9
    assert t_ext > 0.0


def test_lv_decay_drives_extinction_with_no_prey(sim):
    """With no prey, predator decays exponentially and goes extinct."""
    # Prey absent -> dP/dt = -gamma P
    t, N, P = sim.integrate_lotka_volterra(
        N0=0.0, P0=50.0,
        t_end=20.0, dt=0.01,
    )
    assert P[-1] < P[0]
    # Eventually drops to the extinction floor
    assert P[-1] < 5.0


# --------------------------------------------------------------------- #
# (e) May-Wigner instability for large random food webs                 #
# --------------------------------------------------------------------- #

def test_may_wigner_threshold_formula(sim):
    """sigma_c = 1 / sqrt(N C) for self_regulation = 1."""
    N, C = 50, 0.4
    sigma_c = sim.may_wigner_threshold(N, C, self_regulation=1.0)
    expected = 1.0 / math.sqrt(N * C)
    assert abs(sigma_c - expected) < 1e-12


def test_small_random_web_is_stable(sim):
    """A random web well below the May-Wigner threshold is stable."""
    N, C = 20, 0.2
    sigma_c = sim.may_wigner_threshold(N, C)
    sigma = 0.3 * sigma_c   # comfortably below threshold
    geom = EcosystemGeometry.random_food_web(
        n_species=N, connectance=C, sigma=sigma, seed=1,
    )
    sim_local = EcosystemSimulator(geometry=geom)
    evals = sim_local.stability_eigenvalues()
    # Largest real part should be negative (stable)
    assert float(np.max(np.real(evals))) < 0.0


def test_large_random_web_is_unstable(sim):
    """At sigma well ABOVE the May-Wigner threshold the largest real
    eigenvalue is positive => unstable equilibrium."""
    N, C = 100, 0.5
    sigma_c = sim.may_wigner_threshold(N, C)
    sigma = 3.0 * sigma_c   # well above threshold
    geom = EcosystemGeometry.random_food_web(
        n_species=N, connectance=C, sigma=sigma, seed=2,
    )
    sim_local = EcosystemSimulator(geometry=geom)
    evals = sim_local.stability_eigenvalues()
    assert float(np.max(np.real(evals))) > 0.0


def test_may_wigner_ensemble_transition(sim):
    """Across a small ensemble, the FRACTION of unstable webs increases
    from ~0 (sigma << sigma_c) to ~1 (sigma >> sigma_c)."""
    N, C = 30, 0.3
    sigma_c = sim.may_wigner_threshold(N, C)
    n_trials = 12

    def fraction_unstable(sigma):
        unstable = 0
        for s in range(n_trials):
            geom = EcosystemGeometry.random_food_web(
                n_species=N, connectance=C, sigma=sigma, seed=100 + s,
            )
            sim_local = EcosystemSimulator(geometry=geom)
            evals = sim_local.stability_eigenvalues()
            if float(np.max(np.real(evals))) > 0.0:
                unstable += 1
        return unstable / n_trials

    f_low = fraction_unstable(0.3 * sigma_c)
    f_high = fraction_unstable(3.0 * sigma_c)
    assert f_low < 0.2
    assert f_high > 0.8


# --------------------------------------------------------------------- #
# (f) Linear food chain baseline stability                              #
# --------------------------------------------------------------------- #

def test_linear_food_chain_geometry_shape():
    """Linear chain has trophic levels 1..N and N x N matrix."""
    geom = EcosystemGeometry.linear_food_chain(n_levels=4)
    assert geom.n_species == 4
    assert geom.interaction_matrix.shape == (4, 4)
    assert np.allclose(geom.trophic_levels, np.arange(1, 5))


def test_linear_food_chain_has_predator_prey_edges():
    """Off-diagonal links exist for adjacent trophic levels."""
    geom = EcosystemGeometry.linear_food_chain(n_levels=4, link_strength=0.5)
    edges = geom.edges()
    # Expect an upward (predator gains) and downward (prey loses) link
    pairs = {(s, t) for (s, t, _) in edges}
    for i in range(3):
        assert (i, i + 1) in pairs
        assert (i + 1, i) in pairs


# --------------------------------------------------------------------- #
# Geometry helpers                                                      #
# --------------------------------------------------------------------- #

def test_predator_prey_geometry_matrix_signs(predator_prey_geom):
    """A[0,0] = +alpha (prey growth); A[1,1] = -gamma (predator death)."""
    A = predator_prey_geom.interaction_matrix
    assert A[0, 0] > 0.0
    assert A[1, 1] < 0.0
    assert A[0, 1] < 0.0  # predation hurts prey
    assert A[1, 0] > 0.0  # prey feeds predator


def test_random_food_web_connectance_close_to_target():
    """Empirical connectance of a random web matches the requested C."""
    rng_seed = 123
    N = 100
    C_target = 0.25
    geom = EcosystemGeometry.random_food_web(
        n_species=N, connectance=C_target, sigma=0.1, seed=rng_seed,
    )
    measured = geom.connectance()
    # Allow some sampling error (1-sigma ~ sqrt(C(1-C)/N(N-1)))
    assert abs(measured - C_target) < 0.05


def test_node_positions_have_distinct_y_per_level():
    """Trophic-level y-coordinates are distinct across levels."""
    geom = EcosystemGeometry.linear_food_chain(n_levels=4)
    pos = geom.node_positions_2d()
    ys = pos[:, 1]
    assert len(set(ys.tolist())) == 4


# --------------------------------------------------------------------- #
# Substrate metadata                                                    #
# --------------------------------------------------------------------- #

def test_geometry_signature_has_lambda_qcd(predator_prey_geom):
    sig = predator_prey_geom.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
    assert sig["n_species"] == 2


def test_simulator_signature_records_threshold(sim):
    sig = sim.substrate_signature()
    assert sig["extinction_threshold"] == EXTINCTION_THRESHOLD_DEFAULT
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
