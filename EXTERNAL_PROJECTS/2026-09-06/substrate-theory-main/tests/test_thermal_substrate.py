"""Tests for thermal substrate equilibrium simulator."""
import numpy as np
import pytest

from src.stiff_medium.thermal_substrate import (
    ThermalSubstrate,
    SIGMA_SB_TRUE,
)


@pytest.fixture(scope="module")
def equilibrated_low_T():
    sub = ThermalSubstrate(N=12, dx=1.0, rho=1.0, K=1.0, gamma=1.0,
                           dt=0.05, seed=1)
    sub.equilibrate(T=1.0, n_steps=2500)
    return sub


def test_equipartition_kinetic(equilibrated_low_T):
    """<1/2 rho v^2> per dof = 1/2 k_B T (k_B=1, so KE per dof = T/2)."""
    sub = equilibrated_low_T
    KE = sub.kinetic_energy_density()  # 1/2 rho <v^2>
    expected = 0.5 * 1.0  # T = 1
    assert abs(KE - expected) / expected < 0.05, (
        f"equipartition KE off: got {KE}, expected {expected}")


def test_equipartition_kinetic_temp_match(equilibrated_low_T):
    """Kinetic temperature 2*KE/k_B should match set T within 5%."""
    sub = equilibrated_low_T
    info = sub.equipartition_check()
    Tk = info["kinetic_temperature"]
    assert abs(Tk - 1.0) < 0.05, f"kinetic T = {Tk}"


def test_energy_density_linear_in_T():
    """Classical lattice field -> u(T) linear in T (equipartition).

    At fixed N, dx the energy per unit volume should grow ~ k_B T.
    Slope check at a couple of temperatures.
    """
    sub = ThermalSubstrate(N=10, dx=1.0, rho=1.0, K=1.0, gamma=1.0,
                           dt=0.05, seed=2)
    Ts = np.array([0.5, 1.0, 2.0])
    us = []
    for T in Ts:
        us.append(sub.energy_density(float(T), n_equil=1200, n_avg=120))
    us = np.array(us)
    # ratios should track T-ratios within 10%
    r1 = us[1] / us[0]; r2 = us[2] / us[0]
    assert abs(r1 - 2.0) < 0.2, f"u(2T)/u(T) = {r1}"
    assert abs(r2 - 4.0) < 0.5, f"u(4T)/u(T) = {r2}"


def test_t4_planck_prediction():
    """Derived Planck u_Pl ~ T^4 (analytic, given c_s, hbar=1).

    Verifies the closed-form scaling exposed by the simulator's t4_scaling().
    """
    sub = ThermalSubstrate(N=8, dx=1.0, rho=1.0, K=1.0, gamma=1.0, dt=0.05, seed=3)
    Ts = np.array([0.5, 1.0, 2.0])
    info = sub.t4_scaling(Ts, n_equil=600, n_avg=60)
    u_pl = info["u_planck_sim"]
    # Planck T^4 ratios exact within rounding
    assert abs(u_pl[1] / u_pl[0] - 16.0) < 1e-6
    assert abs(u_pl[2] / u_pl[0] - 256.0) < 1e-6


def test_stefan_boltzmann_extraction():
    """sigma_SB extracted should match 5.67e-8 W/m^2/K^4 within 5%."""
    sub = ThermalSubstrate(N=10, dx=1.0, rho=1.0, K=1.0, gamma=1.0,
                           dt=0.05, seed=4)
    sigma = sub.stefan_boltzmann_constant(T_dimless=1.0,
                                          n_equil=1500, n_avg=200)
    rel = abs(sigma - SIGMA_SB_TRUE) / SIGMA_SB_TRUE
    assert rel < 0.05, f"sigma_SB = {sigma:.3e}, true = {SIGMA_SB_TRUE:.3e}, rel = {rel:.3%}"


def test_heat_capacity_classical_limit():
    """C_V per dof should approach k_B = 1 (sim units) in classical limit."""
    sub = ThermalSubstrate(N=8, dx=1.0, rho=1.0, K=1.0, gamma=1.0,
                           dt=0.05, seed=5)
    sub.equilibrate(T=1.0, n_steps=2000)
    C = sub.heat_capacity(samples=300, sample_every=4, T=1.0)
    # Each scalar dof contributes 1 k_B (1/2 from KE, 1/2 from PE).
    # Stochastic estimate is noisy at modest sample size; allow generous band.
    assert 0.4 < C < 1.6, f"C_V per dof = {C}"


def test_mode_occupation_classical_equipartition():
    """<E_k> for each mode ~ k_B T in classical limit."""
    sub = ThermalSubstrate(N=8, dx=1.0, rho=1.0, K=1.0, gamma=1.0,
                           dt=0.05, seed=6)
    sub.equilibrate(T=1.0, n_steps=1500)
    info = sub.mode_occupation(n_avg=80, sample_every=4, T=1.0)
    Ek = info["E_k"]
    # exclude k=0 (zero mode under-resolved)
    Ek_phys = Ek[1:]
    mean = float(np.mean(Ek_phys))
    assert 0.5 < mean < 1.7, (
        f"<E_k> mean = {mean}, expected ~ k_B T = 1")
