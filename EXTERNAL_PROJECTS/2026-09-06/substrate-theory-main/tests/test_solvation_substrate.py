"""Tests for solvation_substrate."""

import math

import numpy as np
import pytest

from stiff_medium.solvation_substrate import (
    DEBYE_TAU,
    DIELECTRIC_CONSTANTS,
    DIELECTRIC_INFINITY,
    EXPERIMENTAL_HYDRATION_KJ,
    ION_CHARGE,
    ION_RADII_A,
    ION_RADII_BORN_A,
    SolvationGeometry,
    SolvationSimulator,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _close_pct(a, b, pct):
    return abs(a - b) / abs(b) * 100.0 < pct


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def test_geometry_unknown_ion_raises():
    with pytest.raises(ValueError):
        SolvationGeometry(ion="Xx")


def test_geometry_shell_count_and_positions():
    g = SolvationGeometry(ion="Na+", n_shell=6, shell_count=1)
    pos = g.shell_positions(0)
    assert pos.shape == (6, 3)
    # all on one sphere
    radii = np.linalg.norm(pos - g.solute_position, axis=1)
    assert np.allclose(radii, radii[0])
    # radius matches r_ion + offset
    assert math.isclose(radii[0], g.shell_radius_A(0), rel_tol=1e-9)


def test_geometry_orientations_point_correctly():
    # Cation: dipole points away from ion (-r̂)
    g_cat = SolvationGeometry(ion="Na+", n_shell=6)
    pos = g_cat.shell_positions(0) - g_cat.solute_position
    rhat = pos / np.linalg.norm(pos, axis=1, keepdims=True)
    mu = g_cat.shell_orientations(0)
    # mu = -rhat for cation
    np.testing.assert_allclose(mu, -rhat, atol=1e-12)

    # Anion: dipole points toward the ion (+r̂)
    g_an = SolvationGeometry(ion="Cl-", n_shell=6)
    pos = g_an.shell_positions(0) - g_an.solute_position
    rhat = pos / np.linalg.norm(pos, axis=1, keepdims=True)
    mu = g_an.shell_orientations(0)
    np.testing.assert_allclose(mu, rhat, atol=1e-12)


def test_geometry_multi_shell_radii_increase():
    g = SolvationGeometry(ion="Mg2+", n_shell=6, shell_count=3)
    r0 = g.shell_radius_A(0)
    r1 = g.shell_radius_A(1)
    r2 = g.shell_radius_A(2)
    assert r2 > r1 > r0
    # All-shell stacking gives N = n_shell * shell_count points
    all_pos = g.all_shell_positions()
    assert all_pos.shape == (6 * 3, 3)


def test_geometry_ion_dipole_energy_negative():
    # Both cations and anions get favourable (negative) energy
    for ion in ("Na+", "Cl-", "Mg2+"):
        g = SolvationGeometry(ion=ion, n_shell=6, shell_count=1)
        E = g.ion_dipole_energy_kJmol()
        assert E < 0.0, ion


# ---------------------------------------------------------------------------
# Born solvation formula
# ---------------------------------------------------------------------------

def test_born_vacuum_to_vacuum_zero():
    # Born ΔG vanishes when initial == final permittivity
    e = SolvationSimulator.born_solvation_kJmol(1, 1.5, eps_r=1.0)
    assert math.isclose(e, 0.0, abs_tol=1e-12)


def test_born_scales_with_z_squared():
    e1 = SolvationSimulator.born_solvation_kJmol(1, 1.0, 78.5)
    e2 = SolvationSimulator.born_solvation_kJmol(2, 1.0, 78.5)
    # ΔG ∝ Z² so e2 / e1 = 4
    assert math.isclose(e2 / e1, 4.0, rel_tol=1e-6)


def test_born_scales_inversely_with_radius():
    e_small = SolvationSimulator.born_solvation_kJmol(1, 1.0, 78.5)
    e_big = SolvationSimulator.born_solvation_kJmol(1, 2.0, 78.5)
    # |ΔG| ∝ 1/r so smaller cavity is more favourable
    assert abs(e_small) > abs(e_big)
    assert math.isclose(e_small / e_big, 2.0, rel_tol=1e-6)


def test_born_negative_in_water():
    e = SolvationSimulator.born_solvation_kJmol(1, 1.5, 78.5)
    assert e < 0.0


# ---------------------------------------------------------------------------
# Hydration energies vs experiment (within 5 %)
# ---------------------------------------------------------------------------

def test_Na_hydration_within_5pct():
    sim = SolvationSimulator()
    e = sim.hydration_energy_kJmol("Na+", "water")
    assert _close_pct(e, EXPERIMENTAL_HYDRATION_KJ["Na+"], 5.0), e


def test_Cl_hydration_within_5pct():
    sim = SolvationSimulator()
    e = sim.hydration_energy_kJmol("Cl-", "water")
    assert _close_pct(e, EXPERIMENTAL_HYDRATION_KJ["Cl-"], 5.0), e


def test_Mg_hydration_within_5pct():
    sim = SolvationSimulator()
    e = sim.hydration_energy_kJmol("Mg2+", "water")
    assert _close_pct(e, EXPERIMENTAL_HYDRATION_KJ["Mg2+"], 5.0), e


@pytest.mark.parametrize(
    "ion", ["Na+", "K+", "Li+", "Cl-", "F-", "Mg2+", "Ca2+"]
)
def test_all_ions_within_5pct(ion):
    sim = SolvationSimulator()
    e = sim.hydration_energy_kJmol(ion, "water")
    assert _close_pct(e, EXPERIMENTAL_HYDRATION_KJ[ion], 5.0), (ion, e)


def test_hydration_unknown_ion_raises():
    sim = SolvationSimulator()
    with pytest.raises(ValueError):
        sim.hydration_energy_kJmol("Xx", "water")


def test_hydration_unknown_solvent_raises():
    sim = SolvationSimulator()
    with pytest.raises(ValueError):
        sim.hydration_energy_kJmol("Na+", "milk")


# ---------------------------------------------------------------------------
# Dielectric constants
# ---------------------------------------------------------------------------

def test_water_dielectric_constant():
    sim = SolvationSimulator()
    assert sim.dielectric_constant("water") == 78.5


def test_methanol_dielectric_constant():
    sim = SolvationSimulator()
    assert sim.dielectric_constant("methanol") == 32.7


def test_dmso_dielectric_constant():
    sim = SolvationSimulator()
    assert sim.dielectric_constant("DMSO") == 47.0


def test_hexane_dielectric_constant():
    sim = SolvationSimulator()
    assert sim.dielectric_constant("hexane") == 1.9


def test_dielectric_unknown_raises():
    sim = SolvationSimulator()
    with pytest.raises(ValueError):
        sim.dielectric_constant("benzaldehyde")


# ---------------------------------------------------------------------------
# Debye relaxation
# ---------------------------------------------------------------------------

def test_debye_low_freq_limit():
    # ω → 0  ⇒  ε(ω) → ε_s
    eps = SolvationSimulator.debye_eps(0.0, 78.5, 5.2, 8.27e-12)
    assert math.isclose(eps.real, 78.5, rel_tol=1e-9)
    assert math.isclose(eps.imag, 0.0, abs_tol=1e-9)


def test_debye_high_freq_limit():
    # ω → ∞  ⇒  ε(ω) → ε_∞
    eps = SolvationSimulator.debye_eps(1e20, 78.5, 5.2, 8.27e-12)
    assert math.isclose(eps.real, 5.2, rel_tol=1e-3)


def test_debye_loss_peak_at_omega_tau_one():
    # ε″ peaks at ω·τ = 1; check the peak location of a sweep.
    sim = SolvationSimulator()
    spec = sim.water_dielectric_spectrum()
    f = spec["f_Hz"]
    eps_imag = spec["eps_imag"]
    omega = 2.0 * math.pi * f
    tau = DEBYE_TAU["water"]
    omega_tau = omega * tau
    # peak should be at omega*tau ~ 1
    i_peak = int(np.argmax(eps_imag))
    # within a factor of 2 (log-spaced grid)
    assert 0.5 < omega_tau[i_peak] < 2.0


def test_debye_real_part_drops_through_relaxation():
    # ε'(low f) ≈ ε_s, ε'(high f) ≈ ε_∞
    sim = SolvationSimulator()
    spec = sim.water_dielectric_spectrum()
    eps_real = spec["eps_real"]
    assert math.isclose(eps_real[0], DIELECTRIC_CONSTANTS["water"], rel_tol=1e-3)
    assert eps_real[-1] < DIELECTRIC_CONSTANTS["water"] / 4.0
    assert math.isclose(eps_real[-1], DIELECTRIC_INFINITY["water"], rel_tol=0.1)


# ---------------------------------------------------------------------------
# Compare table & report
# ---------------------------------------------------------------------------

def test_compare_hydration_returns_full_table():
    sim = SolvationSimulator()
    table = sim.compare_hydration(("Na+", "Cl-", "Mg2+"))
    assert {"Na+", "Cl-", "Mg2+"} <= set(table)
    for row in table.values():
        assert {"calc_kJmol", "exp_kJmol", "abs_err", "rel_err_pct"} <= set(row)


def test_report_string():
    s = SolvationSimulator().report()
    assert "Substrate Solvation" in s
    assert "Na+" in s
    assert "ε_r" in s


def test_ion_charge_table_consistent():
    # Every ion radius entry must have a charge
    assert set(ION_RADII_A) == set(ION_CHARGE) == set(ION_RADII_BORN_A)
