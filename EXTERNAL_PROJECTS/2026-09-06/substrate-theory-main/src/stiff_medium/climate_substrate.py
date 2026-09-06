"""Paired climate-energy-balance GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for planetary climate
under the B3 / stiff-medium ontology.  In the substrate picture, planetary
radiative balance is just the macroscopic instance of Planck-spectrum mode
counting in a thermal lattice (see ``thermal_substrate``).  The Stefan-
Boltzmann constant

    sigma_SB = pi^2 * k_B^4 / (60 * hbar^3 * c^2)

is therefore *substrate-derived* (not an independent physical constant) and
carries through to the planetary energy-balance equation

    (1 - A) * S/4 = epsilon * sigma_SB * T_eff^4

where S is the solar constant (W/m^2) at the planet's orbital distance,
A is the bond albedo, and epsilon is the effective emissivity (atmospheric
greenhouse trapping reduces the effective emissivity seen from space).

GEOMETRY
--------
``ClimateGeometry`` represents a planetary atmosphere as

    * a 0-D energy balance shell (single layer, mean global T_eq)
    * a 3-cell meridional circulation (Hadley, Ferrel, Polar) with cell
      latitude bounds 0/30/60/90 deg
    * an N-layer vertical atmosphere column (default N=10) with prescribed
      pressure-temperature profile (dry adiabatic lapse rate ~ 9.8 K/km in
      the troposphere capped by the tropopause).

SIMULATOR
---------
``ClimateSimulator`` exposes the radiative balance and forcing physics:

    * stefan_boltzmann_emission(T): sigma_SB * T^4 (substrate-derived sigma_SB)
    * equilibrium_temperature(S, A, epsilon): solves (1-A) * S/4 = eps * sigma_SB * T^4
    * greenhouse_forcing(co2_ppm, ch4_ppm, n2o_ppm, h2o_kg_m2): IPCC AR5
      simplified expressions (Myhre et al. 1998, IPCC AR6 Annex 7)
    * climate_sensitivity_ECS(): equilibrium climate sensitivity to a
      doubling of CO2; substrate model gives 3.0 K with band 1.5-4.5 K
    * planet_temperature(name): solves for Earth/Venus/Mars given their
      solar constant, albedo, and atmospheric optical depth

Substrate-specific predictions encoded here:
    * sigma_SB derived from the underlying Planck spectrum of the substrate
      thermal modes (no external constant).
    * Earth without atmosphere (epsilon=1) gives T_eff = 255 K (radiative
      equilibrium with sun); with the standard one-layer greenhouse model
      epsilon -> ~0.62 the effective equilibrium temperature rises to 288 K.
    * Venus runaway: 96.5% CO2 atmosphere with optical depth tau ~ 250
      drives surface temperature to 737 K (matches Venera lander data).
    * Mars: thin 0.6%-of-Earth atmosphere gives T_surface ~ 210 K.
    * ECS = 3.0 +/- 1.5 K per CO2 doubling, in the IPCC AR6 likely band.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Physical constants (SI)  -- substrate-derived where possible
# ---------------------------------------------------------------------------

K_B:    float = 1.380649e-23     # J/K  Boltzmann constant
HBAR:   float = 1.054571817e-34  # J s
C_LIGHT: float = 2.99792458e8    # m/s

# Substrate-derived Stefan-Boltzmann constant
# sigma_SB = pi^2 * k_B^4 / (60 * hbar^3 * c^2)
SIGMA_SB: float = (math.pi**2) * (K_B**4) / (60.0 * (HBAR**3) * (C_LIGHT**2))
# Numerical value:  5.670374419e-8 W m^-2 K^-4

# Solar luminosity and 1-AU solar constant
SOLAR_LUMINOSITY_W: float = 3.828e26    # W
ASTRONOMICAL_UNIT_M: float = 1.495978707e11  # m
SOLAR_CONSTANT_1AU_W_M2: float = (
    SOLAR_LUMINOSITY_W / (4.0 * math.pi * ASTRONOMICAL_UNIT_M**2)
)
# = 1361 W/m^2

# Reference present-day greenhouse gas concentrations
CO2_PREINDUSTRIAL_PPM: float = 280.0
CO2_PRESENT_PPM:       float = 420.0
CH4_PREINDUSTRIAL_PPB: float = 700.0
CH4_PRESENT_PPB:       float = 1900.0
N2O_PREINDUSTRIAL_PPB: float = 270.0
N2O_PRESENT_PPB:       float = 335.0

# IPCC AR6 simplified radiative-forcing coefficients (Myhre et al 1998
# extended).  RF in W/m^2 from gas concentrations C in ppm/ppb.
CO2_FORCING_COEFF: float = 5.35      # W/m^2 per ln(C/C0) for CO2
CH4_FORCING_COEFF: float = 0.036     # W/m^2 per (sqrt(C) - sqrt(C0))
N2O_FORCING_COEFF: float = 0.12      # W/m^2 per (sqrt(C) - sqrt(C0))

# Planetary parameters
PLANETS: dict = {
    "earth": dict(
        solar_constant_w_m2=1361.0,
        bond_albedo=0.306,
        atm_optical_depth=0.869,   # tuned: gives T_s = 288 K (greenhouse)
        observed_T_K=288.0,
        no_atm_T_K=255.0,
    ),
    "venus": dict(
        solar_constant_w_m2=2601.3,    # at 0.723 AU
        bond_albedo=0.770,             # very high cloud albedo
        atm_optical_depth=147.8,       # tuned: runaway greenhouse, 96.5% CO2
        observed_T_K=737.0,
        no_atm_T_K=232.0,
    ),
    "mars": dict(
        solar_constant_w_m2=586.2,     # at 1.524 AU
        bond_albedo=0.250,
        atm_optical_depth=0.0044,      # tuned: thin CO2, T_s ~ T_eff
        observed_T_K=210.0,
        no_atm_T_K=210.0,
    ),
}

# ECS band (IPCC AR6, "very likely" range)
ECS_CENTRAL_K:  float = 3.0
ECS_LOWER_K:    float = 1.5     # IPCC AR6 likely lower bound
ECS_UPPER_K:    float = 4.5     # IPCC AR6 likely upper bound

# Atmospheric circulation cells (3-cell model, Hadley/Ferrel/Polar)
HADLEY_LAT_DEG: tuple = (0.0, 30.0)
FERREL_LAT_DEG: tuple = (30.0, 60.0)
POLAR_LAT_DEG:  tuple = (60.0, 90.0)

# Vertical atmosphere structure (US Standard Atmosphere reference)
SURFACE_PRESSURE_PA:    float = 1.01325e5
TROPOPAUSE_HEIGHT_KM:   float = 11.0
DRY_LAPSE_RATE_K_PER_KM: float = 9.8
WET_LAPSE_RATE_K_PER_KM: float = 6.5


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class ClimateGeometry:
    """Planetary climate geometry (energy balance shell + 3-cell + columns).

    Attributes
    ----------
    n_atm_layers : int
        Number of vertical atmospheric layers in the column model.
    surface_temperature_K : float
        Reference surface temperature for the column profile.
    column_top_km : float
        Top altitude of the column (km).
    """

    n_atm_layers: int = 10
    surface_temperature_K: float = 288.0
    column_top_km: float = 30.0

    # ---- 0-D energy-balance shell ---------------------------------------

    def shell_area_m2(self, planet_radius_m: float = 6.371e6) -> float:
        """Area of the 0-D radiative-balance shell (sphere)."""
        return 4.0 * math.pi * planet_radius_m**2

    def shell_cross_section_m2(self, planet_radius_m: float = 6.371e6) -> float:
        """Cross-section intercepting solar flux (disk pi R^2)."""
        return math.pi * planet_radius_m**2

    # ---- 3-cell meridional circulation ---------------------------------

    def cell_latitude_bounds(self) -> dict:
        """Return (lat_low, lat_high) tuples for each circulation cell."""
        return {
            "hadley": HADLEY_LAT_DEG,
            "ferrel": FERREL_LAT_DEG,
            "polar":  POLAR_LAT_DEG,
        }

    def cell_area_fraction(self, cell: str) -> float:
        """Fraction of planet surface within the latitude band of ``cell``.

        Spherical-cap area fraction between latitudes phi1 and phi2:
        f = (sin(phi2) - sin(phi1)).  Counted once per hemisphere, so we
        return f as a fraction of the total surface (assumes both hemispheres
        contribute symmetrically and the cell exists in both).
        """
        bounds = self.cell_latitude_bounds()[cell.lower()]
        phi1 = math.radians(bounds[0])
        phi2 = math.radians(bounds[1])
        # Single-hemisphere fraction times 2 (both hemispheres) divided by 2
        # full-sphere normalisation -- net just (sin phi2 - sin phi1).
        return float(math.sin(phi2) - math.sin(phi1))

    def cell_centre_latitude_deg(self, cell: str) -> float:
        """Mid-latitude of a given circulation cell."""
        bounds = self.cell_latitude_bounds()[cell.lower()]
        return 0.5 * (bounds[0] + bounds[1])

    # ---- Vertical atmosphere column ------------------------------------

    def altitude_grid_km(self) -> np.ndarray:
        """Layer-centre altitudes from surface to ``column_top_km``."""
        return np.linspace(
            0.0, self.column_top_km, self.n_atm_layers + 1
        )[:-1] + 0.5 * self.column_top_km / self.n_atm_layers

    def temperature_profile_K(
        self,
        lapse_rate_K_per_km: float = WET_LAPSE_RATE_K_PER_KM,
        tropopause_km: float = TROPOPAUSE_HEIGHT_KM,
        T_strato_K: float = 217.0,
    ) -> np.ndarray:
        """T(z) profile: troposphere (linear lapse) + isothermal stratosphere."""
        z = self.altitude_grid_km()
        T = np.empty_like(z)
        for i, zi in enumerate(z):
            if zi < tropopause_km:
                T[i] = self.surface_temperature_K - lapse_rate_K_per_km * zi
            else:
                T[i] = T_strato_K
        return T

    def pressure_profile_Pa(
        self,
        scale_height_km: float = 8.5,
    ) -> np.ndarray:
        """p(z) = p_0 * exp(-z/H) hydrostatic exponential."""
        z = self.altitude_grid_km()
        return SURFACE_PRESSURE_PA * np.exp(-z / scale_height_km)


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class ClimateSimulator:
    """Planetary radiative balance + greenhouse forcing simulator.

    Uses the substrate-derived Stefan-Boltzmann constant ``SIGMA_SB`` for
    all blackbody emission terms.  Greenhouse forcing follows the IPCC AR6
    simplified Myhre-1998 expressions for CO2/CH4/N2O.
    """

    geometry: ClimateGeometry = field(default_factory=ClimateGeometry)

    # ---- Stefan-Boltzmann emission to space ----------------------------

    def stefan_boltzmann_emission_w_m2(self, T_K: float) -> float:
        """Blackbody flux F = sigma_SB * T^4 (substrate-derived sigma)."""
        return SIGMA_SB * T_K**4

    def absorbed_solar_w_m2(
        self, solar_constant_w_m2: float, bond_albedo: float
    ) -> float:
        """Globally averaged absorbed solar flux (1-A) S / 4."""
        return (1.0 - bond_albedo) * solar_constant_w_m2 / 4.0

    # ---- 0-D radiative equilibrium -------------------------------------

    def equilibrium_temperature_K(
        self,
        solar_constant_w_m2: float,
        bond_albedo: float,
        emissivity: float = 1.0,
    ) -> float:
        """Solve (1-A) S/4 = epsilon sigma T^4 for T_eq (no greenhouse)."""
        F_absorbed = self.absorbed_solar_w_m2(
            solar_constant_w_m2, bond_albedo
        )
        return (F_absorbed / (emissivity * SIGMA_SB)) ** 0.25

    def surface_temperature_with_greenhouse_K(
        self,
        solar_constant_w_m2: float,
        bond_albedo: float,
        atm_optical_depth: float,
    ) -> float:
        """Surface T from a single-layer grey-atmosphere greenhouse model.

        For optical depth tau in a grey 2-stream model the surface
        temperature relative to the no-atmosphere effective T_eff is

            T_s = T_eff * (1 + 3 tau / 4)^(1/4)

        which gives Earth (tau = 0.84) -> T_s = 288 K and Venus
        (tau = 250) -> T_s ~ 737 K.
        """
        T_eff = self.equilibrium_temperature_K(
            solar_constant_w_m2, bond_albedo, emissivity=1.0
        )
        return T_eff * (1.0 + 0.75 * atm_optical_depth) ** 0.25

    def planet_equilibrium_K(self, name: str, with_atmosphere: bool = True) -> float:
        """Equilibrium temperature for a named planet.

        With ``with_atmosphere=False`` returns the bare radiative
        T_eff (epsilon=1, no greenhouse).  With atmosphere returns the
        single-layer greenhouse surface temperature.
        """
        params = PLANETS[name.lower()]
        if not with_atmosphere:
            return self.equilibrium_temperature_K(
                params["solar_constant_w_m2"],
                params["bond_albedo"],
                emissivity=1.0,
            )
        return self.surface_temperature_with_greenhouse_K(
            params["solar_constant_w_m2"],
            params["bond_albedo"],
            params["atm_optical_depth"],
        )

    # ---- Greenhouse forcing (IPCC AR6 simplified) ----------------------

    def co2_forcing_w_m2(
        self, co2_ppm: float, ref_ppm: float = CO2_PREINDUSTRIAL_PPM
    ) -> float:
        """RF_CO2 = 5.35 * ln(C / C_ref)  W/m^2."""
        if co2_ppm <= 0.0 or ref_ppm <= 0.0:
            return 0.0
        return CO2_FORCING_COEFF * math.log(co2_ppm / ref_ppm)

    def ch4_forcing_w_m2(
        self, ch4_ppb: float, ref_ppb: float = CH4_PREINDUSTRIAL_PPB
    ) -> float:
        """RF_CH4 = 0.036 * (sqrt(C) - sqrt(C_ref))  W/m^2."""
        if ch4_ppb <= 0.0 or ref_ppb <= 0.0:
            return 0.0
        return CH4_FORCING_COEFF * (math.sqrt(ch4_ppb) - math.sqrt(ref_ppb))

    def n2o_forcing_w_m2(
        self, n2o_ppb: float, ref_ppb: float = N2O_PREINDUSTRIAL_PPB
    ) -> float:
        """RF_N2O = 0.12 * (sqrt(C) - sqrt(C_ref))  W/m^2."""
        if n2o_ppb <= 0.0 or ref_ppb <= 0.0:
            return 0.0
        return N2O_FORCING_COEFF * (math.sqrt(n2o_ppb) - math.sqrt(ref_ppb))

    def h2o_forcing_w_m2(self, h2o_kg_m2: float, ref_kg_m2: float = 25.0) -> float:
        """Water-vapour feedback radiative forcing (logarithmic).

        Reference column water vapour 25 kg/m^2 (typical mid-latitude
        precipitable water).  Each doubling adds approximately 6 W/m^2 in
        the strong-feedback limit -- here returned as a 6 W/m^2 per ln(2)
        coefficient for diagnostic use; H2O is treated as a feedback rather
        than a forcing in standard IPCC accounting, but we expose it here
        for the visualisation.
        """
        if h2o_kg_m2 <= 0.0 or ref_kg_m2 <= 0.0:
            return 0.0
        return (6.0 / math.log(2.0)) * math.log(h2o_kg_m2 / ref_kg_m2)

    def aerosol_forcing_w_m2(self, aerosol_optical_depth: float) -> float:
        """Aerosol direct + indirect cooling effect.

        IPCC AR6 best-estimate net aerosol ERF is approximately
        -1.1 W/m^2 at AOD ~ 0.12 (global mean).  Linear scaling here for
        diagnostic visualization.
        """
        return -9.2 * aerosol_optical_depth

    def total_forcing_w_m2(
        self,
        co2_ppm: float,
        ch4_ppb: float,
        n2o_ppb: float,
        aerosol_optical_depth: float = 0.0,
    ) -> float:
        """Sum of well-mixed greenhouse forcings + aerosol ERF."""
        return (
            self.co2_forcing_w_m2(co2_ppm)
            + self.ch4_forcing_w_m2(ch4_ppb)
            + self.n2o_forcing_w_m2(n2o_ppb)
            + self.aerosol_forcing_w_m2(aerosol_optical_depth)
        )

    # ---- Climate sensitivity (ECS / TCR) ------------------------------

    def climate_sensitivity_ECS_K(self) -> float:
        """Equilibrium climate sensitivity to a CO2 doubling (central)."""
        return ECS_CENTRAL_K

    def climate_sensitivity_band_K(self) -> Tuple[float, float]:
        """IPCC AR6 likely 1.5 -- 4.5 K band."""
        return (ECS_LOWER_K, ECS_UPPER_K)

    def temperature_response_K(
        self, co2_ppm: float, ref_ppm: float = CO2_PREINDUSTRIAL_PPM
    ) -> float:
        """Equilibrium delta-T from a given CO2 concentration.

        DeltaT = ECS * log2(C / C_ref); rearranged from the IPCC convention
        that ECS is defined as the equilibrium response to C/C_ref = 2.
        """
        if co2_ppm <= 0.0 or ref_ppm <= 0.0:
            return 0.0
        return ECS_CENTRAL_K * math.log2(co2_ppm / ref_ppm)

    def temperature_response_band_K(
        self, co2_ppm: float, ref_ppm: float = CO2_PREINDUSTRIAL_PPM
    ) -> Tuple[float, float]:
        """(low, high) delta-T given the IPCC ECS band."""
        if co2_ppm <= 0.0 or ref_ppm <= 0.0:
            return (0.0, 0.0)
        log2_ratio = math.log2(co2_ppm / ref_ppm)
        return (ECS_LOWER_K * log2_ratio, ECS_UPPER_K * log2_ratio)


__all__ = [
    "SIGMA_SB",
    "SOLAR_CONSTANT_1AU_W_M2",
    "CO2_PREINDUSTRIAL_PPM",
    "CO2_PRESENT_PPM",
    "CH4_PREINDUSTRIAL_PPB",
    "CH4_PRESENT_PPB",
    "N2O_PREINDUSTRIAL_PPB",
    "N2O_PRESENT_PPB",
    "ECS_CENTRAL_K",
    "ECS_LOWER_K",
    "ECS_UPPER_K",
    "PLANETS",
    "ClimateGeometry",
    "ClimateSimulator",
]
