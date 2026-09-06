"""Tests for BoundStateSpectrum solver."""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from stiff_medium.bound_state_spectrum import (  # noqa: E402
    BoundStateSpectrum,
    SubstrateParams,
    RYDBERG_eV,
)


@pytest.fixture
def solver():
    return BoundStateSpectrum(SubstrateParams())


def test_hydrogen_ground_state(solver):
    eigs, analytic = solver.solve_hydrogen(n_states=1, l=0,
                                           r_max_bohr=80.0, N=8000)
    # ground state at -13.6 eV +/- 1%
    assert abs(eigs[0] - (-RYDBERG_eV)) / RYDBERG_eV < 0.01, (
        f"H ground state {eigs[0]:.4f} eV not within 1% of -13.6 eV")


def test_hydrogen_excited_p_d_f(solver):
    # 2p (l=1, n=2): -13.6/4
    eigs_p, _ = solver.solve_hydrogen(n_states=1, l=1, r_max_bohr=120.0, N=8000)
    assert abs(eigs_p[0] - (-RYDBERG_eV / 4)) / (RYDBERG_eV / 4) < 0.02
    # 3d (l=2, n=3): -13.6/9
    eigs_d, _ = solver.solve_hydrogen(n_states=1, l=2, r_max_bohr=200.0, N=8000)
    assert abs(eigs_d[0] - (-RYDBERG_eV / 9)) / (RYDBERG_eV / 9) < 0.03
    # 4f (l=3, n=4): -13.6/16
    eigs_f, _ = solver.solve_hydrogen(n_states=1, l=3, r_max_bohr=300.0, N=8000)
    assert abs(eigs_f[0] - (-RYDBERG_eV / 16)) / (RYDBERG_eV / 16) < 0.05


def test_harmonic_oscillator(solver):
    n = 4
    vals, analytic = solver.solve_harmonic(n_states=n, omega=1.0, mass=1.0,
                                           hbar=1.0, r_max=15.0, N=8000)
    # 3D radial l=0 ladder: hbar omega (2 n_r + 3/2)
    for i in range(n):
        assert abs(vals[i] - analytic[i]) / analytic[i] < 0.02, (
            f"HO level {i}: {vals[i]} vs {analytic[i]}")


def test_finite_well_count(solver):
    # Choose well so theory predicts a known number of bound l=0 states.
    # z0 = a sqrt(2 m V0)/hbar; for a=1, m=1, hbar=1, V0=50 -> z0=10 -> 4 bound states
    bound, expected = solver.solve_finite_well(V0=50.0, a=1.0, mass=1.0,
                                               hbar=1.0, r_max=25.0, N=8000, k=10)
    # numerical Dirichlet box adds at most one boundary-shifted state
    assert abs(len(bound) - expected) <= 1, (
        f"Found {len(bound)} bound states, expected ~{expected}")


def test_deuteron_K4(solver):
    BE_pred, BE_emp, bound = solver.solve_substrate_K4(
        n_states=3, V0_MeV=35.0, a_fm=2.1, r_max_fm=40.0, N=8000)
    # Expect order-of-magnitude agreement with empirical 2.2245 MeV.
    # Standard square-well fit hits ~2.2 MeV; allow 50% tolerance to absorb
    # grid/discretisation drift and the coarse square-well approximation.
    assert BE_pred > 0.0, "no bound state found"
    assert abs(BE_pred - BE_emp) / BE_emp < 0.5, (
        f"deuteron BE_pred={BE_pred:.3f} MeV vs empirical {BE_emp} MeV")


def test_drag_derived_electron_mass():
    """m_e from drag-cone-bouncing matches 0.511 MeV within 5%."""
    sub = SubstrateParams.electron_anchor(K=1.0, rho=1.0, xi=1.0)
    m_drag = sub.drag_mass_SI()
    M_E_kg = 9.1093837015e-31
    rel_err = abs(m_drag - M_E_kg) / M_E_kg
    assert rel_err < 0.05, (
        f"drag-derived m_e = {m_drag:.4e} kg, expected {M_E_kg:.4e}, "
        f"rel_err={rel_err:.2%}")
    # Static method should agree with the bound method.
    m_static = BoundStateSpectrum.effective_mass_from_drag(
        gamma=sub.gamma, K=sub.K, rho=sub.rho, xi=sub.xi)
    assert abs(m_static - m_drag) / m_drag < 1e-12


def test_hydrogen_with_drag_mass(solver):
    """Hydrogen ground state = -13.6 eV using drag-derived m_e."""
    sub = SubstrateParams.electron_anchor()
    eigs, analytic, m_drag = solver.solve_hydrogen_substrate(
        sub=sub, n_states=1, l=0, r_max_bohr=80.0, N=8000)
    assert abs(eigs[0] - (-RYDBERG_eV)) / RYDBERG_eV < 0.01, (
        f"H ground state with drag-derived m_e: {eigs[0]:.4f} eV "
        f"vs -13.6 eV")
    info = solver.verify_drag_mass_consistency(sub=sub)
    assert info["passes"], info


def test_mass_independent_of_drag_anchor(solver):
    """Different self-consistent (K, rho, xi, gamma) anchors give same m_e."""
    anchors = [
        SubstrateParams.electron_anchor(K=1.0, rho=1.0, xi=1.0),
        SubstrateParams.electron_anchor(K=4.0, rho=1.0, xi=1.0),     # c_s x2
        SubstrateParams.electron_anchor(K=1.0, rho=2.5, xi=0.7),
        SubstrateParams.electron_anchor(K=9.0, rho=4.0, xi=2.3),
    ]
    masses = [s.drag_mass_SI() for s in anchors]
    m0 = masses[0]
    for m in masses[1:]:
        assert abs(m - m0) / m0 < 1e-10, (
            f"anchor-mass inconsistency: {m} vs {m0}")
    # Hydrogen ground state should also be invariant.
    energies = []
    for s in anchors:
        eigs, _, _ = solver.solve_hydrogen_substrate(
            sub=s, n_states=1, l=0, r_max_bohr=80.0, N=8000)
        energies.append(eigs[0])
    e0 = energies[0]
    for e in energies[1:]:
        assert abs(e - e0) / abs(e0) < 1e-6, (
            f"hydrogen E_gs differs across anchors: {e} vs {e0}")


def test_report_runs(solver):
    solver.solve_harmonic(n_states=3, omega=1.0, mass=1.0, hbar=1.0,
                          r_max=12.0, N=2000)
    text = solver.report()
    assert "BoundStateSpectrum" in text
    assert "Harmonic" in text
