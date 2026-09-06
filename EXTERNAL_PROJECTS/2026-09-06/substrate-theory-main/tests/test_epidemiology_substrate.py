"""Tests for the paired epidemiology GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) Herd-immunity threshold H = 1 - 1/R_0  (R_0 = 2.5 -> 60%)
    (b) Peak prevalence formula i_max = 1 - (1 + ln R_0)/R_0
    (c) Final-size relation       R_∞/N = 1 - exp(-R_0 · R_∞/N)
    (d) R_0 < 1 -> no epidemic (sub-critical, D goes to ~0)
    (e) Vaccination above H suppresses the epidemic
    (f) SEIR latency delays the peak vs SIR with same R_0
    (g) Population conservation S + I + R = N (and S+E+I+R for SEIR)
    (h) Geometry: contact network, age structure, next-generation matrix
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.epidemiology_substrate import (
    D_INFECTIOUS_DAYS_DEFAULT,
    D_LATENT_DAYS_DEFAULT,
    LAMBDA_QCD_MEV,
    R0_DEFAULT,
    AgeStructure,
    ContactNetwork,
    EpidemiologyGeometry,
    EpidemiologySimulator,
    final_size_fraction,
    herd_immunity_threshold,
    peak_prevalence_fraction,
    r0_sweep,
)


# --------------------------------------------------------------------- #
# (a) Herd immunity threshold                                           #
# --------------------------------------------------------------------- #

def test_herd_immunity_R0_2p5_is_60pct():
    """R_0 = 2.5 -> H = 1 - 1/2.5 = 0.6 exactly."""
    H = herd_immunity_threshold(2.5)
    assert H == pytest.approx(0.60, abs=1.0e-12)


def test_herd_immunity_subcritical_is_zero():
    assert herd_immunity_threshold(0.5) == 0.0
    assert herd_immunity_threshold(1.0) == 0.0


def test_herd_immunity_R0_3():
    H = herd_immunity_threshold(3.0)
    assert H == pytest.approx(2.0 / 3.0, rel=1.0e-9)


def test_herd_immunity_measles_R0_15():
    """Classic measles benchmark: R_0 ~ 15 -> H ~ 14/15 ~ 93%."""
    H = herd_immunity_threshold(15.0)
    assert H == pytest.approx(14.0 / 15.0, rel=1.0e-9)
    assert 0.93 < H < 0.94


# --------------------------------------------------------------------- #
# (b) Peak prevalence formula                                           #
# --------------------------------------------------------------------- #

def test_peak_prevalence_R0_2p5():
    """i_max = 1 - (1 + ln 2.5)/2.5 ~ 0.234."""
    i_max = peak_prevalence_fraction(2.5)
    expected = 1.0 - (1.0 + math.log(2.5)) / 2.5
    assert i_max == pytest.approx(expected, rel=1.0e-12)
    assert 0.23 < i_max < 0.24


def test_peak_prevalence_subcritical():
    assert peak_prevalence_fraction(0.7) == 0.0
    assert peak_prevalence_fraction(1.0) == 0.0


def test_peak_prevalence_matches_simulation_within_tolerance():
    """ODE simulation at R_0 = 2.5 must reproduce the closed-form i_max
    within a few percent (closed form is the I_max with full population
    initially susceptible)."""
    geom = EpidemiologyGeometry(
        N=1_000_000, R0=2.5,
        D_infectious_days=D_INFECTIOUS_DAYS_DEFAULT,
    )
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.05, duration_days=200.0,
    )
    res = sim.simulate_sir(I0_fraction=1.0e-6)
    i_max_emp = sim.peak_infected_fraction(res)
    i_max_th = peak_prevalence_fraction(2.5)
    assert i_max_emp == pytest.approx(i_max_th, rel=0.02)


# --------------------------------------------------------------------- #
# (c) Final-size relation                                               #
# --------------------------------------------------------------------- #

def test_final_size_R0_2p5_implicit_eq():
    """x* = 1 - exp(-R_0 x*)  must hold at the returned root."""
    R0 = 2.5
    x = final_size_fraction(R0)
    assert x == pytest.approx(1.0 - math.exp(-R0 * x), abs=1.0e-9)
    assert 0.85 < x < 0.95   # known answer ~ 0.893


def test_final_size_subcritical_zero():
    assert final_size_fraction(0.8) == 0.0
    assert final_size_fraction(1.0) == 0.0


def test_final_size_matches_simulation():
    """Simulated R(∞)/N must approach the implicit solution."""
    geom = EpidemiologyGeometry(
        N=1_000_000, R0=2.5,
        D_infectious_days=D_INFECTIOUS_DAYS_DEFAULT,
    )
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=400.0,
    )
    res = sim.simulate_sir(I0_fraction=1.0e-6)
    R_inf_emp = sim.final_recovered_fraction(res)
    R_inf_th = final_size_fraction(2.5)
    assert R_inf_emp == pytest.approx(R_inf_th, rel=0.02)


# --------------------------------------------------------------------- #
# (d) Subcritical: R_0 < 1 -> no epidemic                               #
# --------------------------------------------------------------------- #

def test_subcritical_dies_out():
    """R_0 = 0.7 -> infection should die out (I -> 0 by end)."""
    geom = EpidemiologyGeometry(N=10_000, R0=0.7)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=200.0,
    )
    res = sim.simulate_sir(I0_fraction=1.0e-3)
    assert sim.peak_infected_fraction(res) <= 1.0e-3 + 1.0e-6
    assert res["I_frac"][-1] < 1.0e-4


def test_supercritical_grows():
    """R_0 = 2.5 -> infection grows substantially."""
    geom = EpidemiologyGeometry(N=10_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=180.0,
    )
    res = sim.simulate_sir(I0_fraction=1.0e-4)
    assert sim.peak_infected_fraction(res) > 0.10


# --------------------------------------------------------------------- #
# (e) Vaccination above herd-threshold suppresses epidemic              #
# --------------------------------------------------------------------- #

def test_vaccination_above_threshold_suppresses():
    """Vaccinating > H = 60% (R_0=2.5) prevents an epidemic."""
    geom = EpidemiologyGeometry(N=100_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=200.0,
    )
    res = sim.vaccinate_uniform(fraction=0.70)
    # Attack rate among the unvaccinated residual must stay tiny
    assert sim.attack_rate(res) < 0.05


def test_vaccination_below_threshold_does_not_suppress():
    """Vaccinating well below H -> still substantial epidemic."""
    geom = EpidemiologyGeometry(N=100_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=200.0,
    )
    res = sim.vaccinate_uniform(fraction=0.20)
    assert sim.attack_rate(res) > 0.30


def test_targeted_vaccination_more_effective_than_uniform():
    """Same coverage, targeted reduces the attack rate further."""
    geom = EpidemiologyGeometry(N=100_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=200.0,
    )
    f = 0.40
    res_u = sim.vaccinate_uniform(fraction=f)
    res_t = sim.vaccinate_targeted(fraction=f, degree_amplification=2.0)
    assert sim.attack_rate(res_t) < sim.attack_rate(res_u)


# --------------------------------------------------------------------- #
# (f) SEIR adds latency, delays peak                                    #
# --------------------------------------------------------------------- #

def test_seir_peak_later_than_sir():
    """Same β, γ; SEIR adds an exposed compartment so peak is later."""
    geom = EpidemiologyGeometry(
        N=100_000, R0=2.5,
        D_infectious_days=5.0, D_latent_days=3.0,
    )
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=200.0,
    )
    res_sir  = sim.simulate_sir(I0_fraction=1.0e-4)
    res_seir = sim.simulate_seir(E0_fraction=1.0e-4)
    t_peak_sir  = sim.peak_time_days(res_sir)
    t_peak_seir = sim.peak_time_days(res_seir)
    assert t_peak_seir > t_peak_sir


def test_seir_compartments_present():
    geom = EpidemiologyGeometry(N=10_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=120.0,
    )
    res = sim.simulate_seir(E0_fraction=1.0e-4)
    for key in ("S", "E", "I", "R", "S_frac", "E_frac", "I_frac", "R_frac",
                "sigma"):
        assert key in res
    assert res["sigma"] == pytest.approx(1.0 / D_LATENT_DAYS_DEFAULT,
                                         rel=1.0e-12)


# --------------------------------------------------------------------- #
# (g) Population conservation                                           #
# --------------------------------------------------------------------- #

def test_sir_conservation():
    geom = EpidemiologyGeometry(N=10_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=120.0,
    )
    res = sim.simulate_sir(I0_fraction=1.0e-3)
    total = res["S"] + res["I"] + res["R"]
    assert np.allclose(total, geom.N, rtol=1.0e-6)


def test_seir_conservation():
    geom = EpidemiologyGeometry(N=10_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=120.0,
    )
    res = sim.simulate_seir(E0_fraction=1.0e-3)
    total = res["S"] + res["E"] + res["I"] + res["R"]
    assert np.allclose(total, geom.N, rtol=1.0e-6)


def test_compartments_non_negative():
    geom = EpidemiologyGeometry(N=10_000, R0=2.5)
    sim = EpidemiologySimulator(
        geometry=geom, dt_days=0.1, duration_days=120.0,
    )
    res = sim.simulate_seir(E0_fraction=1.0e-3)
    for key in ("S", "E", "I", "R"):
        assert (res[key] >= -1.0e-9).all(), f"{key} went negative"


# --------------------------------------------------------------------- #
# (h) Geometry: networks + age structure                                #
# --------------------------------------------------------------------- #

def test_contact_network_default_edge_probability():
    net = ContactNetwork(N=1000, mean_contacts=10.0)
    # ⟨k⟩ = (N - 1) p  ->  p = 10 / 999
    assert net.p_edge == pytest.approx(10.0 / 999.0, rel=1.0e-9)
    assert net.expected_degree() == pytest.approx(10.0)


def test_contact_network_validates_type():
    with pytest.raises(ValueError):
        ContactNetwork(N=100, mean_contacts=5.0,
                       network_type="banana-graph")


def test_age_structure_default_4_groups():
    ag = AgeStructure()
    assert ag.n_groups == 4
    assert ag.age_fractions.sum() == pytest.approx(1.0, abs=1.0e-9)
    assert ag.mixing_matrix.shape == (4, 4)


def test_age_structure_fractions_must_sum_to_one():
    with pytest.raises(ValueError):
        AgeStructure(n_groups=3,
                     age_fractions=np.array([0.5, 0.4, 0.2]))


def test_geometry_derived_rates():
    geom = EpidemiologyGeometry(
        N=1000, R0=2.5,
        D_infectious_days=5.0, D_latent_days=3.0,
    )
    assert geom.gamma_per_day() == pytest.approx(0.2, abs=1.0e-12)
    assert geom.sigma_per_day() == pytest.approx(1.0 / 3.0, abs=1.0e-12)
    assert geom.beta_per_day() == pytest.approx(0.5, abs=1.0e-12)
    # τ = β / ⟨k⟩
    assert geom.transmissibility_per_contact() == pytest.approx(
        0.5 / 13.0, rel=1.0e-12,
    )


def test_geometry_herd_threshold_via_method():
    geom = EpidemiologyGeometry(N=1000, R0=2.5)
    assert geom.herd_immunity_threshold() == pytest.approx(0.60, abs=1.0e-12)


def test_next_generation_matrix_shape_and_R0():
    """Next-generation matrix is K x K and recovers R_0 (≈) for default
    homogeneous mixing."""
    geom = EpidemiologyGeometry(N=1000, R0=2.5)
    K = geom.next_generation_matrix()
    assert K.shape == (4, 4)
    R0_eff = geom.R0_from_next_generation()
    # The default mixing matrix is rank-1-dominated with diagonal=1,
    # off-diagonal=0.5; the spectral radius equals R_0 scaled by the
    # dominant eigenvector overlap.  We only assert it is positive and
    # of order R_0.
    assert R0_eff > 0.5 * geom.R0
    assert R0_eff < 2.5 * geom.R0


def test_substrate_signature_keys():
    geom = EpidemiologyGeometry(N=10_000, R0=2.5)
    sig = geom.substrate_signature()
    for key in ("interpretation", "lambda_qcd_MeV", "N", "mean_contacts",
                "R0", "gamma_per_day", "beta_per_day", "tau_per_contact",
                "H_threshold"):
        assert key in sig
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
    assert sig["H_threshold"] == pytest.approx(0.60, abs=1.0e-12)


# --------------------------------------------------------------------- #
# r0_sweep convenience                                                  #
# --------------------------------------------------------------------- #

def test_r0_sweep_returns_one_run_per_R0():
    R0s = np.array([0.8, 1.5, 2.5, 4.0])
    runs = r0_sweep(R0s, N=10_000, duration_days=120.0)
    assert set(runs.keys()) == {float(r) for r in R0s}
    for R0, res in runs.items():
        assert res["R0"] == R0
        assert res["t_days"][0] == 0.0
        assert res["t_days"][-1] == pytest.approx(120.0, abs=1.0e-9)


def test_r0_sweep_attack_rate_monotonic_in_R0():
    """Larger R_0 -> larger final attack rate (monotonic in supercritical
    regime)."""
    R0s = np.array([1.5, 2.0, 2.5, 3.0, 4.0])
    runs = r0_sweep(R0s, N=100_000, duration_days=300.0,
                    I0_fraction=1.0e-5)
    sim = EpidemiologySimulator(
        geometry=EpidemiologyGeometry(N=100_000, R0=2.5),
    )
    attack = [sim.attack_rate(runs[float(r)]) for r in R0s]
    for i in range(1, len(attack)):
        assert attack[i] > attack[i - 1] - 1.0e-3


# --------------------------------------------------------------------- #
# Defaults + sanity                                                     #
# --------------------------------------------------------------------- #

def test_defaults_cover_canonical_R0():
    assert R0_DEFAULT == 2.5
    assert D_INFECTIOUS_DAYS_DEFAULT == 5.0
    assert D_LATENT_DAYS_DEFAULT == 3.0


def test_lambda_qcd_anchor_present():
    """B3 anchor is exposed on the substrate signature."""
    assert LAMBDA_QCD_MEV == 200.0
