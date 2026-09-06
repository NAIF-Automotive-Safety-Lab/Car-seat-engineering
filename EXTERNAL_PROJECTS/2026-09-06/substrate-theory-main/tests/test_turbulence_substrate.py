"""Tests for the paired fluid-turbulence GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) Kolmogorov inertial-range slope is -5/3 (analytic and synthetic field)
    (b) universal Kolmogorov constant C_K = 1.5
    (c) ε = u_rms³ / L recovers the large-eddy dissipation estimate
    (d) Reynolds number Re = U L / ν, regime transitions at Re_crit ≈ 2300
    (e) Strouhal number St = f D / U with St_cylinder ≈ 0.21 ± 0.02
    (f) substrate viscosity bridge ν = γ ξ²/2 reproduces molecular viscosity
    (g) Kármán vortex street has the right wavelength λ = D / St
    (h) vorticity finite-difference operator: ω = ∂_x v - ∂_y u
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.turbulence_substrate import (
    KOLMOGOROV_CONSTANT,
    KOLMOGOROV_EXPONENT,
    NU_AIR_M2_S,
    NU_WATER_M2_S,
    RE_CRIT_PIPE,
    RE_CRIT_PIPE_LOWER,
    RE_CRIT_PIPE_UPPER,
    STROUHAL_CYLINDER,
    TurbulenceGeometry,
    TurbulenceSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> TurbulenceGeometry:
    return TurbulenceGeometry(Lx=1.0, Ly=1.0, Nx=128, Ny=128)


@pytest.fixture
def sim(geom: TurbulenceGeometry) -> TurbulenceSimulator:
    return TurbulenceSimulator(geometry=geom, nu_m2_s=NU_AIR_M2_S)


# --------------------------------------------------------------------- #
# (a) Kolmogorov inertial-range slope -5/3                              #
# --------------------------------------------------------------------- #

def test_kolmogorov_exponent_is_minus_five_thirds():
    """The canonical inertial-range exponent must be -5/3."""
    assert math.isclose(KOLMOGOROV_EXPONENT, -5.0 / 3.0, rel_tol=1e-15)


def test_kolmogorov_spectrum_analytic_slope(sim):
    """E(k) = C_K ε^(2/3) k^(-5/3) — fitted log-log slope must be -5/3."""
    k = np.logspace(0.0, 4.0, 50)
    E = sim.kolmogorov_spectrum(k, epsilon_w_kg=1.0)
    slope = sim.fit_inertial_slope(k, E, k_min=1.0e1, k_max=1.0e3)
    assert abs(slope - (-5.0 / 3.0)) < 1e-3, (
        f"analytic slope {slope:.4f} (expected -5/3)"
    )


def test_kolmogorov_constant_is_one_point_five():
    """Universal Kolmogorov constant C_K = 1.5 (Sreenivasan 1995)."""
    assert math.isclose(KOLMOGOROV_CONSTANT, 1.5)


def test_synthetic_field_shell_spectrum_slope(sim):
    """The synthetic 2-D field with k^(-5/3) shell amplitude must give a
    fitted spectral slope close to -5/3 over the inertial decade."""
    field = sim.synthetic_kolmogorov_field(
        epsilon_w_kg=1.0, L_m=0.5, seed=42,
    )
    k_bin, E_bin = sim.shell_average_spectrum(
        field["u"], field["v"], field["kx"], field["ky"], n_bins=48,
    )
    # Fit inside the inertial decade
    k_min = 4.0 * 2.0 * math.pi / sim.geometry.Lx
    k_max = 0.4 * math.pi * sim.geometry.Nx / sim.geometry.Lx
    slope = sim.fit_inertial_slope(k_bin, E_bin, k_min=k_min, k_max=k_max)
    assert abs(slope - (-5.0 / 3.0)) < 0.30, (
        f"synthetic slope {slope:.3f} (expected -5/3 ± 0.30)"
    )


# --------------------------------------------------------------------- #
# (b) Energy dissipation rate (large-eddy estimate)                     #
# --------------------------------------------------------------------- #

def test_dissipation_rate_units(sim):
    """ε has units m²/s³; for u_rms = 1 m/s, L = 1 m, ε = 1 m²/s³."""
    eps = sim.energy_dissipation_rate(u_rms=1.0, L_m=1.0)
    assert math.isclose(eps, 1.0)


def test_dissipation_rate_scales_as_u_cubed(sim):
    """ε ∝ u_rms³ for fixed L."""
    eps1 = sim.energy_dissipation_rate(u_rms=1.0, L_m=1.0)
    eps2 = sim.energy_dissipation_rate(u_rms=2.0, L_m=1.0)
    assert math.isclose(eps2 / eps1, 8.0, rel_tol=1e-12)


# --------------------------------------------------------------------- #
# (c) Reynolds number transition                                        #
# --------------------------------------------------------------------- #

def test_reynolds_number_definition(sim):
    """Re = U L / ν at standard air viscosity."""
    Re = sim.reynolds_number(U_m_s=1.0, L_m=NU_AIR_M2_S)
    assert math.isclose(Re, 1.0, rel_tol=1e-12)


def test_re_crit_pipe_value():
    """Pipe-flow critical Reynolds number is the textbook 2300."""
    assert math.isclose(RE_CRIT_PIPE, 2300.0)
    assert RE_CRIT_PIPE_LOWER <= RE_CRIT_PIPE <= RE_CRIT_PIPE_UPPER


def test_laminar_below_2300(sim):
    """Re = 1500 in the lower laminar branch returns 'laminar'."""
    assert sim.laminar_or_turbulent(1500.0) == "laminar"


def test_transitional_at_3000(sim):
    """Re = 3000 returns 'transitional'."""
    assert sim.laminar_or_turbulent(3000.0) == "transitional"


def test_turbulent_above_4000(sim):
    """Re = 1e5 returns 'turbulent'."""
    assert sim.laminar_or_turbulent(1.0e5) == "turbulent"


def test_water_pipe_transition_velocity(sim):
    """For water (ν = 1e-6 m²/s) in a 2-cm pipe, the laminar–turbulent
    transition occurs at U ≈ Re_crit · ν / D ≈ 0.115 m/s."""
    D = 0.02  # 2-cm pipe
    U_crit = RE_CRIT_PIPE * NU_WATER_M2_S / D
    assert abs(U_crit - 0.115) < 0.005


# --------------------------------------------------------------------- #
# (d) Strouhal number for vortex shedding                               #
# --------------------------------------------------------------------- #

def test_strouhal_cylinder_value():
    """Canonical St = 0.21 for a circular cylinder in Re ∼ 100 – 10⁵."""
    assert abs(STROUHAL_CYLINDER - 0.21) < 1e-12


def test_strouhal_definition(sim):
    """St = f D / U."""
    St = sim.strouhal_number(shedding_freq_hz=2.1, diameter_m=0.1, U_m_s=1.0)
    assert math.isclose(St, 0.21)


def test_shedding_frequency_inverse(sim):
    """Round-trip: f = St U / D."""
    f = sim.shedding_frequency(U_m_s=1.0, diameter_m=0.1)
    assert math.isclose(f, 2.1)


def test_vortex_field_wavelength(sim):
    """Generated Kármán street has wavelength λ = U / f = D / St."""
    field = sim.vortex_shedding_field(U_m_s=1.0, diameter_m=0.1)
    lam = field["wavelength_m"]
    assert abs(lam - 0.1 / STROUHAL_CYLINDER) < 1e-9


def test_vortex_field_alternating_sign(sim):
    """Successive vortices in the wake alternate sign — verified via the
    integrated vorticity in two adjacent half-wavelength windows."""
    field = sim.vortex_shedding_field(
        U_m_s=1.0, diameter_m=0.1, n_pairs=4,
    )
    omega = field["vorticity"]
    geom = sim.geometry
    X, Y = geom.meshgrid()
    cyl_x, cyl_y = field["cylinder_xy"]
    lam = field["wavelength_m"]
    # First window (top vortex)
    mask_top = (
        (X > cyl_x + 0.4 * 0.1)
        & (X < cyl_x + 0.4 * 0.1 + lam)
        & (Y > cyl_y)
    )
    mask_bot = (
        (X > cyl_x + 0.4 * 0.1)
        & (X < cyl_x + 0.4 * 0.1 + lam)
        & (Y < cyl_y)
    )
    sum_top = float(omega[mask_top].sum())
    sum_bot = float(omega[mask_bot].sum())
    assert sum_top > 0
    assert sum_bot < 0


# --------------------------------------------------------------------- #
# (e) Substrate-viscosity bridge ν = γ ξ²/2                              #
# --------------------------------------------------------------------- #

def test_substrate_viscosity_bridge(sim):
    """ν = (1/2) γ ξ² recovers the air kinematic viscosity for sensible
    molecular ξ (mean free path) and γ (collision rate).

    Air at STP: mean free path ξ ≈ 6.8e-8 m, collision rate γ ≈ 6.5e9 Hz
    -> ν_predicted = 1/2 · 6.5e9 · (6.8e-8)² ≈ 1.5e-5 m²/s.
    """
    nu = TurbulenceSimulator.viscosity_from_substrate(
        gamma_hz=6.5e9, xi_m=6.8e-8,
    )
    # Match air viscosity within a factor of 2 (collision-rate uncertainty)
    assert 0.5 * NU_AIR_M2_S < nu < 2.0 * NU_AIR_M2_S


def test_viscosity_scales_quadratic_in_xi():
    """ν ∝ ξ² for fixed γ."""
    nu1 = TurbulenceSimulator.viscosity_from_substrate(1.0e10, 1.0e-8)
    nu2 = TurbulenceSimulator.viscosity_from_substrate(1.0e10, 2.0e-8)
    assert math.isclose(nu2 / nu1, 4.0, rel_tol=1e-12)


# --------------------------------------------------------------------- #
# (f) Geometry — vorticity operator + grid                              #
# --------------------------------------------------------------------- #

def test_geometry_grid_size(geom):
    """Mesh has the requested resolution."""
    X, Y = geom.meshgrid()
    assert X.shape == (geom.Ny, geom.Nx)
    assert Y.shape == (geom.Ny, geom.Nx)


def test_vorticity_solid_body_rotation(geom):
    """Solid-body rotation u = -y, v = +x has uniform ω = 2."""
    X, Y = geom.meshgrid()
    cx = 0.5 * geom.Lx
    cy = 0.5 * geom.Ly
    u = -(Y - cy)
    v = +(X - cx)
    dx, dy = geom.cell_size()
    omega = TurbulenceGeometry.vorticity_2d(u, v, dx, dy)
    interior = omega[2:-2, 2:-2]
    assert abs(interior.mean() - 2.0) < 1.0e-3


def test_kolmogorov_eta_scale(sim):
    """η = (ν³/ε)^(1/4) — at ν=1e-5, ε=1: η ≈ 1.78e-4 m."""
    eta = TurbulenceGeometry.kolmogorov_scale_m(
        nu_m2_s=1.0e-5, epsilon_w_kg=1.0,
    )
    assert abs(eta - (1.0e-5) ** 0.75) < 1e-9


def test_cascade_band_edges_ordering(sim):
    """k_L < k_eta — large eddies have smaller wavenumber than dissipative
    eddies."""
    band = sim.geometry.cascade_band_edges(
        nu_m2_s=1.0e-5, epsilon_w_kg=1.0, u_rms=1.0,
    )
    assert band["k_L"] < band["k_eta"]
