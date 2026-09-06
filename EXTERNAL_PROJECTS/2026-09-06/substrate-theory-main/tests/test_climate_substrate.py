"""Tests for the paired climate-energy-balance GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) substrate-derived sigma_SB matches CODATA value
    (b) Earth no-atmosphere effective T = 255 K (radiative equilibrium)
    (c) Earth with atmosphere T_surface = 288 K (greenhouse warming)
    (d) Venus runaway greenhouse drives surface T to 737 K
    (e) Mars thin atmosphere gives surface T ~ 210 K
    (f) ECS in IPCC AR6 likely band 1.5 - 4.5 K
    (g) CO2 doubling forcing ~ 3.7 W/m^2
    (h) 3-cell circulation: Hadley/Ferrel/Polar surface fractions sum to 1
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.climate_substrate import (
    CH4_PREINDUSTRIAL_PPB,
    CH4_PRESENT_PPB,
    CO2_PREINDUSTRIAL_PPM,
    CO2_PRESENT_PPM,
    ECS_CENTRAL_K,
    ECS_LOWER_K,
    ECS_UPPER_K,
    N2O_PREINDUSTRIAL_PPB,
    N2O_PRESENT_PPB,
    PLANETS,
    SIGMA_SB,
    SOLAR_CONSTANT_1AU_W_M2,
    ClimateGeometry,
    ClimateSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> ClimateGeometry:
    return ClimateGeometry(n_atm_layers=10)


@pytest.fixture
def sim() -> ClimateSimulator:
    return ClimateSimulator()


# --------------------------------------------------------------------- #
# (a) Substrate-derived Stefan-Boltzmann constant                       #
# --------------------------------------------------------------------- #

def test_sigma_SB_matches_codata():
    """sigma_SB = pi^2 k_B^4 / (60 hbar^3 c^2) must match CODATA value."""
    codata_sigma = 5.670374419e-8
    rel = abs(SIGMA_SB - codata_sigma) / codata_sigma
    assert rel < 1.0e-6, f"substrate sigma_SB {SIGMA_SB} vs CODATA {codata_sigma}"


def test_solar_constant_at_1AU_close_to_1361():
    """Solar luminosity / 4 pi (1 AU)^2 should give ~1361 W/m^2."""
    assert abs(SOLAR_CONSTANT_1AU_W_M2 - 1361.0) / 1361.0 < 0.01


# --------------------------------------------------------------------- #
# (b) Earth no-atmosphere effective T = 255 K                            #
# --------------------------------------------------------------------- #

def test_earth_no_atm_eq_temp_255K(sim):
    """Earth with no atmosphere should have T_eff ~ 255 K."""
    T = sim.planet_equilibrium_K("earth", with_atmosphere=False)
    # Standard textbook value 254-255 K
    assert abs(T - 255.0) < 2.0, f"Earth no-atm T = {T:.2f} K (expect ~255)"


def test_earth_no_atm_uses_correct_formula(sim):
    """Bare T_eff = ((1-A) S / (4 sigma))^(1/4)."""
    T = sim.equilibrium_temperature_K(1361.0, 0.306, emissivity=1.0)
    expected = ((1.0 - 0.306) * 1361.0 / (4.0 * SIGMA_SB)) ** 0.25
    assert abs(T - expected) < 1.0e-9


# --------------------------------------------------------------------- #
# (c) Earth with atmosphere T_surface = 288 K (greenhouse)              #
# --------------------------------------------------------------------- #

def test_earth_with_atm_surface_288K(sim):
    """Earth with greenhouse atmosphere should give T_surface ~ 288 K."""
    T = sim.planet_equilibrium_K("earth", with_atmosphere=True)
    assert abs(T - 288.0) < 2.0, f"Earth surface T = {T:.2f} K (expect 288)"


def test_greenhouse_raises_temperature(sim):
    """Adding atmosphere must raise T (greenhouse trapping)."""
    T_no_atm = sim.planet_equilibrium_K("earth", with_atmosphere=False)
    T_with   = sim.planet_equilibrium_K("earth", with_atmosphere=True)
    assert T_with > T_no_atm, "greenhouse should warm the surface"
    delta = T_with - T_no_atm
    # Standard greenhouse warming of Earth is ~33 K
    assert 25.0 < delta < 40.0, f"Earth greenhouse delta T = {delta:.1f} K"


# --------------------------------------------------------------------- #
# (d) Venus runaway greenhouse                                           #
# --------------------------------------------------------------------- #

def test_venus_runaway_737K(sim):
    """Venus surface temperature must be ~737 K (Venera lander)."""
    T = sim.planet_equilibrium_K("venus", with_atmosphere=True)
    assert abs(T - 737.0) / 737.0 < 0.05, (
        f"Venus surface T = {T:.1f} K (expect 737)"
    )


def test_venus_no_atm_close_to_earth_no_atm(sim):
    """Without atmosphere, Venus T_eff is comparable to Earth (high albedo)."""
    T = sim.planet_equilibrium_K("venus", with_atmosphere=False)
    # Venus high albedo (0.77) compensates higher solar flux -> ~232 K
    assert 220.0 < T < 250.0, f"Venus no-atm T = {T:.2f} K"


def test_venus_atm_warming_500plus(sim):
    """Venus atmosphere adds at least 500 K of greenhouse warming."""
    T_no = sim.planet_equilibrium_K("venus", with_atmosphere=False)
    T_w = sim.planet_equilibrium_K("venus", with_atmosphere=True)
    assert (T_w - T_no) > 500.0, "Venus runaway should add >500 K"


# --------------------------------------------------------------------- #
# (e) Mars thin atmosphere gives ~210 K                                  #
# --------------------------------------------------------------------- #

def test_mars_surface_210K(sim):
    """Mars surface T ~ 210 K (Viking/Curiosity measurement)."""
    T = sim.planet_equilibrium_K("mars", with_atmosphere=True)
    assert abs(T - 210.0) < 5.0, f"Mars surface T = {T:.2f} K (expect ~210)"


def test_mars_thin_atm_minimal_greenhouse(sim):
    """Mars thin atm contributes < 5 K of greenhouse warming."""
    T_no = sim.planet_equilibrium_K("mars", with_atmosphere=False)
    T_w = sim.planet_equilibrium_K("mars", with_atmosphere=True)
    assert (T_w - T_no) < 5.0, "Mars thin atm should add < 5 K"


# --------------------------------------------------------------------- #
# (f) ECS in IPCC AR6 likely band                                        #
# --------------------------------------------------------------------- #

def test_ecs_central_value(sim):
    """Central ECS should be 3.0 K per CO2 doubling."""
    assert sim.climate_sensitivity_ECS_K() == ECS_CENTRAL_K


def test_ecs_band_in_ipcc_range(sim):
    """ECS likely range 1.5 - 4.5 K per IPCC AR6."""
    low, high = sim.climate_sensitivity_band_K()
    assert low == ECS_LOWER_K and high == ECS_UPPER_K
    assert 1.0 < low < 2.0, "AR6 lower bound should be 1.5 K"
    assert 4.0 < high < 5.0, "AR6 upper bound should be 4.5 K"


def test_temperature_response_doubling(sim):
    """Temperature response at 2x preindustrial CO2 = ECS."""
    dT = sim.temperature_response_K(2.0 * CO2_PREINDUSTRIAL_PPM)
    assert abs(dT - ECS_CENTRAL_K) < 1.0e-9


def test_temperature_response_quadrupling(sim):
    """Temperature response at 4x preindustrial CO2 = 2 * ECS."""
    dT = sim.temperature_response_K(4.0 * CO2_PREINDUSTRIAL_PPM)
    assert abs(dT - 2.0 * ECS_CENTRAL_K) < 1.0e-9


# --------------------------------------------------------------------- #
# (g) CO2 doubling forcing ~ 3.7 W/m^2                                  #
# --------------------------------------------------------------------- #

def test_co2_doubling_forcing(sim):
    """RF for doubled CO2 should be ~ 3.7 W/m^2 (Myhre 1998 / IPCC)."""
    rf = sim.co2_forcing_w_m2(2.0 * CO2_PREINDUSTRIAL_PPM)
    # 5.35 * ln(2) = 3.71 W/m^2
    assert abs(rf - 3.71) < 0.05, f"CO2 doubling RF = {rf:.2f} W/m^2"


def test_present_day_co2_forcing_positive(sim):
    """Present CO2 (420 ppm) vs pre-industrial (280 ppm) should give RF > 2 W/m^2."""
    rf = sim.co2_forcing_w_m2(CO2_PRESENT_PPM)
    assert rf > 2.0, f"Present CO2 RF should exceed 2 W/m^2, got {rf:.2f}"


def test_ch4_forcing_positive_present(sim):
    """Present CH4 (1900 ppb) vs pre-industrial (700 ppb) should give RF > 0.4 W/m^2."""
    rf = sim.ch4_forcing_w_m2(CH4_PRESENT_PPB)
    assert rf > 0.4, f"Present CH4 RF should exceed 0.4 W/m^2, got {rf:.2f}"


def test_n2o_forcing_positive_present(sim):
    """Present N2O (335 ppb) vs pre-industrial (270 ppb) should give RF > 0.1 W/m^2."""
    rf = sim.n2o_forcing_w_m2(N2O_PRESENT_PPB)
    assert rf > 0.1, f"Present N2O RF should exceed 0.1 W/m^2, got {rf:.2f}"


def test_aerosol_forcing_negative(sim):
    """Aerosol direct + indirect ERF must be net negative (cooling)."""
    rf = sim.aerosol_forcing_w_m2(0.12)
    assert rf < 0, "aerosols should be net cooling"


def test_total_present_day_forcing(sim):
    """Total present-day anthropogenic RF should be in the 2-3 W/m^2 range."""
    rf = sim.total_forcing_w_m2(
        co2_ppm=CO2_PRESENT_PPM,
        ch4_ppb=CH4_PRESENT_PPB,
        n2o_ppb=N2O_PRESENT_PPB,
        aerosol_optical_depth=0.12,
    )
    assert 1.0 < rf < 4.0, f"Total RF = {rf:.2f} W/m^2 (expect ~2-3)"


# --------------------------------------------------------------------- #
# (h) 3-cell circulation geometry                                       #
# --------------------------------------------------------------------- #

def test_three_cells_present(geom):
    """All three meridional cells (Hadley/Ferrel/Polar) must be defined."""
    cells = geom.cell_latitude_bounds()
    assert set(cells.keys()) == {"hadley", "ferrel", "polar"}


def test_cell_bounds_partition_hemisphere(geom):
    """Cells must tile 0 to 90 deg without gap or overlap."""
    cells = geom.cell_latitude_bounds()
    assert cells["hadley"] == (0.0, 30.0)
    assert cells["ferrel"] == (30.0, 60.0)
    assert cells["polar"]  == (60.0, 90.0)


def test_cell_area_fractions_sum_to_one(geom):
    """Sum of single-hemisphere cell area fractions = sin(90 deg) = 1."""
    total = sum(
        geom.cell_area_fraction(c) for c in ("hadley", "ferrel", "polar")
    )
    assert abs(total - 1.0) < 1.0e-9


def test_hadley_largest_fraction(geom):
    """Hadley (tropical) cell has the largest surface fraction."""
    h = geom.cell_area_fraction("hadley")
    f = geom.cell_area_fraction("ferrel")
    p = geom.cell_area_fraction("polar")
    assert h > f > p, f"expected Hadley > Ferrel > Polar; got {h}/{f}/{p}"


# --------------------------------------------------------------------- #
# Vertical atmosphere column                                            #
# --------------------------------------------------------------------- #

def test_temperature_decreases_in_troposphere(geom):
    """T(z) must decrease with z below the tropopause (lapse rate)."""
    T = geom.temperature_profile_K()
    # First few layers (in troposphere) should be monotonically decreasing
    z = geom.altitude_grid_km()
    trop_mask = z < 11.0
    T_trop = T[trop_mask]
    diffs = np.diff(T_trop)
    assert np.all(diffs < 0.0), "tropospheric T must decrease with altitude"


def test_pressure_decreases_with_altitude(geom):
    """p(z) must decrease monotonically with altitude (hydrostatic)."""
    p = geom.pressure_profile_Pa()
    diffs = np.diff(p)
    assert np.all(diffs < 0.0), "pressure must decrease with altitude"


def test_surface_pressure_1bar(geom):
    """Surface pressure of the column should be ~1 bar = 1.01325e5 Pa."""
    p = geom.pressure_profile_Pa()
    # First layer is centred at half-spacing above ground
    z = geom.altitude_grid_km()
    p0 = 1.01325e5 * math.exp(-z[0] / 8.5)
    assert abs(p[0] - p0) < 1.0e3


# --------------------------------------------------------------------- #
# Stefan-Boltzmann emission                                             #
# --------------------------------------------------------------------- #

def test_emission_at_288K_is_390(sim):
    """sigma * 288^4 ~ 390 W/m^2 (Earth surface upward flux)."""
    F = sim.stefan_boltzmann_emission_w_m2(288.0)
    assert abs(F - 390.0) < 5.0, f"flux = {F:.1f} W/m^2 (expect ~390)"


def test_emission_scales_as_T4(sim):
    """Doubling T should multiply emission by 16."""
    F1 = sim.stefan_boltzmann_emission_w_m2(150.0)
    F2 = sim.stefan_boltzmann_emission_w_m2(300.0)
    assert abs(F2 / F1 - 16.0) < 1.0e-9
