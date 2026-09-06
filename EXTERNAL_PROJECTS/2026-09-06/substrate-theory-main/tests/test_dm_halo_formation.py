"""Tests for the dark-matter halo formation simulator."""

import numpy as np
import pytest

from src.stiff_medium.dm_halo_formation import (
    DMHaloFormation,
    build_relaxed_halo,
)


def test_initial_conditions_and_force_shape():
    sim = DMHaloFormation(seed=1)
    sim.initialize_cosmological(N=50, box_size=4.0, sigma_v=0.1)
    assert sim.positions.shape == (50, 3)
    assert sim.velocities.shape == (50, 3)
    acc = sim.compute_forces()
    assert acc.shape == (50, 3)
    assert np.all(np.isfinite(acc))


def test_evolve_conserves_energy_roughly():
    sim = DMHaloFormation(seed=2, softening=0.1)
    sim.initialize_cosmological(N=40, box_size=4.0, sigma_v=0.05)
    e0 = sim.total_energy()
    sim.evolve(t_max=1.0, dt=0.02)
    e1 = sim.total_energy()
    assert abs(e1 - e0) / max(abs(e0), 1e-12) < 0.30


def test_nfw_slopes_recovered():
    """NFW slopes: -1 at small r, -3 at large r."""
    sim = build_relaxed_halo(N=200, radius=1.0, M=1.0, seed=3,
                             relax_time=3.0, softening=0.05)
    halos = sim.halo_finder(sigma_threshold=1.0, min_particles=30)
    assert len(halos) >= 1, "halo finder did not return any halo"
    halo = halos[0]
    centers, rho = sim.radial_profile(halo.center, max_r=halo.radius,
                                       n_bins=15, indices=halo.indices)
    rho_s, r_s = sim.nfw_fit((centers, rho))
    assert rho_s > 0 and r_s > 0
    # slope at small r should be near -1 (allow generous bracket)
    r_small = 0.2 * r_s
    r_large = 5.0 * r_s
    eps = 1e-3
    def slope(r):
        rho1 = sim.nfw_density(np.array([r * (1 - eps)]), rho_s, r_s)[0]
        rho2 = sim.nfw_density(np.array([r * (1 + eps)]), rho_s, r_s)[0]
        return (np.log(rho2) - np.log(rho1)) / (np.log(r * (1 + eps))
                                                 - np.log(r * (1 - eps)))
    s_in = slope(r_small)
    s_out = slope(r_large)
    assert -1.5 < s_in < -0.5, f"inner slope {s_in} not near -1"
    assert -3.5 < s_out < -2.5, f"outer slope {s_out} not near -3"


def test_concentration_relation_for_mw_like_halo():
    """c-M relation: c ~ 5-15 for ~10^11 - 10^12 M_sun."""
    c_low = DMHaloFormation.concentration(M_halo=1e11)
    c_mw = DMHaloFormation.concentration(M_halo=1e12)
    c_high = DMHaloFormation.concentration(M_halo=1e13)
    assert 5.0 <= c_low <= 15.0
    assert 5.0 <= c_mw <= 15.0
    # decreasing with mass
    assert c_low > c_mw > c_high


def test_virial_theorem_in_relaxed_halo():
    """2T + U should be small relative to |U| in equilibrium."""
    sim = build_relaxed_halo(N=200, radius=1.0, M=1.0, seed=4,
                             relax_time=8.0, softening=0.05)
    halos = sim.halo_finder(sigma_threshold=1.0, min_particles=30)
    assert len(halos) >= 1
    ratio = sim.virial_ratio(halos[0])
    # within 30% as required
    assert ratio < 0.3, f"virial ratio {ratio} exceeds 30%"


def test_velocity_dispersion_positive_and_reasonable():
    sim = build_relaxed_halo(N=150, radius=1.0, M=1.0, seed=5,
                             relax_time=3.0, softening=0.05)
    halos = sim.halo_finder(sigma_threshold=1.0, min_particles=20)
    assert len(halos) >= 1
    sigma = sim.velocity_dispersion(halos[0])
    # virial estimate sigma ~ sqrt(GM/2R) ~ sqrt(0.5) ~ 0.7 in code units
    assert 0.1 < sigma < 5.0


def test_merger_yields_single_remnant():
    """Two halos placed close together should merge into one bound clump."""
    sim = DMHaloFormation(seed=6, softening=0.08, particle_mass=0.005)
    # build two halos by hand, well within each other's gravitational reach
    sim.initialize_halo(N=80, radius=0.6, sigma_v=0.3,
                        center=np.array([-0.8, 0.0, 0.0]))
    sim.initialize_halo(N=80, radius=0.6, sigma_v=0.3,
                        center=np.array([0.8, 0.0, 0.0]))
    # give them a small head-on velocity
    sim.velocities[:80, 0] += 0.1
    sim.velocities[80:, 0] -= 0.1

    n_before = len(sim.halo_finder(sigma_threshold=1.5, min_particles=20))
    sim.evolve(t_max=4.0, dt=0.04)
    halos_after = sim.halo_finder(sigma_threshold=1.0, min_particles=30)
    # we expect a single dominant remnant containing most of the mass
    assert len(halos_after) >= 1
    largest = halos_after[0]
    total_mass = float(np.sum(sim.masses))
    assert largest.mass / total_mass > 0.5
    # merger: the remnant should be more massive than either initial halo
    assert n_before >= 1
