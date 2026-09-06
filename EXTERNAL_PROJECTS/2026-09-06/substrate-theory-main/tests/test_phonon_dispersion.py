"""Tests for phonon_dispersion.py."""

import numpy as np
import pytest

from src.stiff_medium.phonon_dispersion import (
    PhononDispersion,
    estimate_debye_diamond,
    estimate_debye_aluminum,
    estimate_debye_copper,
)


def test_k0_acoustic_zero_modes_3d():
    """At k=0 with no optical mass, acoustic phonons must have ω=0."""
    pd = PhononDispersion(K=1.0, rho=1.0, xi=0.0, lattice="fcc", a=1.0)
    omega = pd.frequencies(np.zeros(3))
    # Three acoustic modes at k=0 must vanish
    assert np.all(omega[:3] < 1e-8), f"Expected ω[:3]=0, got {omega[:3]}"


def test_k0_diamond_acoustic_zero_optical_finite():
    """Diamond: 3 acoustic + 3 optical. Acoustic at k=0 must vanish;
    optical at k=0 must be finite."""
    pd = PhononDispersion(K=1.0, rho=1.0, xi=0.0, lattice="diamond", a=1.0)
    omega = pd.frequencies(np.zeros(3))
    assert pd.n_modes == 6
    assert np.all(omega[:3] < 1e-7), f"Acoustic non-zero: {omega[:3]}"
    # Optical branches must be finite even with ξ=0 (different masses or
    # diatomic basis with different bonding gives optical gap)
    assert np.all(omega[3:] > 1e-3), f"Optical zero: {omega[3:]}"


def test_sound_speed_recovery_chain():
    """1D chain: c_s = sqrt(K/ρ) at long wavelength."""
    K, rho = 4.0, 1.0
    pd = PhononDispersion(K=K, rho=rho, xi=0.0, lattice="chain", a=1.0)
    cs = pd.sound_speed(k_small=1e-4)
    expected = np.sqrt(K / rho)
    assert abs(cs - expected) / expected < 0.02, f"cs={cs}, expected≈{expected}"


def test_sound_speed_recovery_square():
    """2D square lattice: longitudinal c_s = sqrt(K/ρ)."""
    K, rho = 2.0, 0.5
    pd = PhononDispersion(K=K, rho=rho, xi=0.0, lattice="square", a=1.0)
    cs = pd.sound_speed(direction=np.array([1.0, 0.0]), k_small=1e-4)
    expected = np.sqrt(K / rho)
    # Longitudinal mode along x has c_L = sqrt(K/ρ); average over L+T includes
    # transverse zero — so check the longitudinal eigenvalue specifically.
    omega = pd.frequencies(1e-4 * np.array([1.0, 0.0]))
    cL = omega.max() / 1e-4
    assert abs(cL - expected) / expected < 0.02, f"cL={cL}, expected≈{expected}"


def test_optical_branches_diatomic_chain():
    """1D diatomic chain has 1 acoustic + 1 optical branch."""
    pd = PhononDispersion(
        K=1.0, rho=1.0, xi=0.0, lattice="chain_diatomic", a=1.0,
        masses=[1.0, 2.0],
    )
    # n_atoms=2, dim=1, so n_modes=2, optical = 1
    opt = pd.optical_branches(np.array([0.1]))
    ac = pd.acoustic_branches(np.array([0.1]))
    assert len(ac) == 1
    assert len(opt) == 1
    # Optical at k=0 should be > 0; acoustic should be small but nonzero
    omega0 = pd.frequencies(np.zeros(1))
    assert omega0[0] < 1e-8, f"acoustic@k=0 should vanish: {omega0[0]}"
    assert omega0[1] > 0.5, f"optical@k=0 should be finite: {omega0[1]}"


def test_dispersion_curve_path_shape():
    pd = PhononDispersion(K=1.0, rho=1.0, lattice="fcc", a=1.0)
    x, omega, ticks = pd.dispersion_curve(["G", "X", "L", "G"], n_points=20)
    assert omega.shape[1] == 3
    assert omega.shape[0] == x.shape[0]
    assert len(ticks) == 4


def test_dos_integrates_to_total_modes():
    """∫ dos dω should equal n_modes (per unit cell, per k-point)."""
    pd = PhononDispersion(K=1.0, rho=1.0, lattice="fcc", a=1.0)
    centers, dos = pd.density_of_states(n_bins=80, n_kgrid=8)
    domega = centers[1] - centers[0]
    total = dos.sum() * domega
    # Approximation: depends on bin width; should match n_modes within ~few%
    assert abs(total - pd.n_modes) / pd.n_modes < 0.1, (
        f"DOS total = {total}, expected ≈ {pd.n_modes}"
    )


def test_debye_temperature_diamond_cross_check():
    """Real-material Debye T for diamond should be ~2230 K within 30%."""
    theta = estimate_debye_diamond()
    assert 1500 < theta < 3000, f"diamond Θ_D = {theta} K, expected ≈ 2230"


def test_debye_temperature_aluminum_cross_check():
    theta = estimate_debye_aluminum()
    assert 300 < theta < 600, f"Al Θ_D = {theta} K, expected ≈ 428"


def test_debye_temperature_copper_cross_check():
    theta = estimate_debye_copper()
    assert 240 < theta < 480, f"Cu Θ_D = {theta} K, expected ≈ 343"


def test_substrate_debye_scales_with_sound_speed():
    """Doubling K (so c_s increases by sqrt(2)) should scale Θ_D by sqrt(2)."""
    pd1 = PhononDispersion(K=1.0, rho=1.0, lattice="fcc", a=1.0)
    pd2 = PhononDispersion(K=2.0, rho=1.0, lattice="fcc", a=1.0)
    t1 = pd1.debye_temperature()
    t2 = pd2.debye_temperature()
    ratio = t2 / t1
    assert abs(ratio - np.sqrt(2)) < 0.05, f"ratio={ratio}, expected √2"


def test_optical_gap_from_xi():
    """Adding ξ > 0 opens a gap in all branches at k=0."""
    pd = PhononDispersion(K=1.0, rho=1.0, xi=0.5, lattice="fcc", a=1.0)
    omega = pd.frequencies(np.zeros(3))
    # All modes should be lifted by sqrt(ξ/ρ)
    expected_gap = np.sqrt(0.5 / 1.0)
    assert np.all(np.abs(omega - expected_gap) < 1e-6), (
        f"With ξ=0.5, expected ω ≈ {expected_gap}, got {omega}"
    )


def test_hermitian_dynamical_matrix():
    pd = PhononDispersion(K=1.0, rho=1.0, lattice="diamond", a=1.0)
    k = np.array([0.3, 0.1, -0.2])
    D = pd.dynamical_matrix(k)
    assert np.allclose(D, D.conj().T, atol=1e-12)
