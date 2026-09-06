"""
Cube DM Cross-Section
=====================

Closed-form spin-independent (SI) and spin-dependent (SD) cross-sections for
the B3-substrate cube cell dark matter candidate (m_DM = 27.5 GeV).

The cube cell is electromagnetically inert at leading dipole/quadrupole order;
its EM interaction is octupole-suppressed (~alpha^2 (k a)^6 a^2 ~ 10^-105 m^2,
handled by cube_cell_dm_simulator). The remaining channel for nucleon
scattering is *gravitational* exchange, which gives a direct-detection
cross-section many orders of magnitude below current and projected sensitivity.

Closed form
-----------
For non-relativistic 2-to-2 elastic gravitational scattering of a DM particle
(mass m_chi) on a nucleon (mass m_N), the t-channel single-graviton exchange
amplitude in the Born approximation gives a Newtonian-like differential cross
section.  Integrating over solid angle with an IR cutoff at the de Broglie
wavelength of the relative motion yields the order-of-magnitude estimate

    sigma_grav  ~  ( 4 pi G m_chi m_N / v^2 )^2

where v is the relative velocity (~ 230 km/s in the Milky Way halo).
This is the "gravitational Rutherford" cross section used as the canonical
floor for any DM candidate that couples only via gravity.

Author: B3 framework, audit_07 follow-up
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import math

from scipy import constants as sc


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
G_NEWTON = sc.G                       # m^3 kg^-1 s^-2
C_LIGHT = sc.c                        # m / s
HBAR = sc.hbar                        # J s
EV = sc.e                             # J / eV
GEV = 1.0e9 * EV                      # J / GeV
M_PROTON = sc.proton_mass             # kg
M_NEUTRON = sc.neutron_mass           # kg
M_NUCLEON = 0.5 * (M_PROTON + M_NEUTRON)
AMU = sc.atomic_mass                  # kg

# Astrophysics
V_HALO = 230.0e3                      # m/s, local DM relative velocity
RHO_DM_LOCAL = 0.4 * GEV / C_LIGHT**2 / 1.0e-6  # 0.4 GeV/cc -> kg / m^3
# 0.4 GeV/c^2 / cm^3 = 0.4 * GEV/c^2 * 1e6 m^-3
RHO_DM_LOCAL = (0.4 * GEV / C_LIGHT**2) * 1.0e6  # kg / m^3

# Direct-detection limits at ~30 GeV (m^2 = cm^2 * 1e-4)
LIMIT_LZ_2024_CM2 = 2.0e-48           # cm^2  (LZ 2024, ~30 GeV WIMP)
LIMIT_PANDAX_2024_CM2 = 3.0e-48       # cm^2  (PandaX-4T 2024, ~27 GeV)
LIMIT_XENONNT_PROJ_CM2 = 1.0e-48      # cm^2  (XENON-nT projected 2027)
NEUTRINO_FLOOR_CM2 = 1.0e-49          # cm^2  (coherent nu-N scattering floor)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _gev_to_kg(m_gev: float) -> float:
    return m_gev * GEV / C_LIGHT**2


def _reduced_mass(m1: float, m2: float) -> float:
    return m1 * m2 / (m1 + m2)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------
@dataclass
class CubeDMCrossSection:
    """
    Spin-independent / spin-dependent nucleon cross-section for the
    27.5 GeV cube DM candidate, gravity-only channel.
    """

    m_dm_gev: float = 27.5
    v_rel: float = V_HALO            # m/s
    rho_dm: float = RHO_DM_LOCAL     # kg/m^3 local halo density
    target_A: int = 131              # default xenon-131
    target_amu: float = 131.293      # xenon atomic mass

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------
    @property
    def m_dm_kg(self) -> float:
        return _gev_to_kg(self.m_dm_gev)

    @property
    def mu_nucleon(self) -> float:
        """Reduced mass DM <-> nucleon (kg)."""
        return _reduced_mass(self.m_dm_kg, M_NUCLEON)

    @property
    def mu_nucleus(self) -> float:
        """Reduced mass DM <-> target nucleus (kg)."""
        m_nuc = self.target_amu * AMU
        return _reduced_mass(self.m_dm_kg, m_nuc)

    # ------------------------------------------------------------------
    # Cross-sections
    # ------------------------------------------------------------------
    def sigma_si_nucleon_m2(self) -> float:
        """
        Spin-independent gravitational cross-section per nucleon, in m^2.

        sigma ~ ( 4 pi G m_DM m_N / v^2 )^2
        """
        prefactor = 4.0 * math.pi * G_NEWTON * self.m_dm_kg * M_NUCLEON
        return (prefactor / self.v_rel**2) ** 2

    def sigma_si_nucleon_cm2(self) -> float:
        return self.sigma_si_nucleon_m2() * 1.0e4

    def sigma_sd_nucleon_cm2(self) -> float:
        """
        Spin-dependent piece. Gravity is universal/spin-blind, so the
        spin-dependent component is suppressed by an additional factor
        ~ (v/c)^2 from spin-orbit-like couplings.
        """
        return self.sigma_si_nucleon_cm2() * (self.v_rel / C_LIGHT) ** 2

    def sigma_si_nucleus_cm2(self) -> float:
        """
        Coherent nuclear cross-section (~A^2 enhancement on the
        per-nucleon piece, with reduced-mass scaling).
        """
        per_nuc = self.sigma_si_nucleon_cm2()
        ratio = (self.mu_nucleus / self.mu_nucleon) ** 2
        return self.target_A**2 * ratio * per_nuc

    # ------------------------------------------------------------------
    # Comparisons / consequences
    # ------------------------------------------------------------------
    def compare_to_limits(self) -> Dict[str, float]:
        """
        Ratio sigma_predicted / sigma_limit (per nucleon, cm^2).

        Values << 1 mean the prediction is well below the limit (consistent
        with observed null results).
        """
        s = self.sigma_si_nucleon_cm2()
        return {
            "sigma_si_nucleon_cm2": s,
            "ratio_to_LZ_2024":      s / LIMIT_LZ_2024_CM2,
            "ratio_to_PandaX_2024":  s / LIMIT_PANDAX_2024_CM2,
            "ratio_to_XENONnT_2027": s / LIMIT_XENONNT_PROJ_CM2,
            "ratio_to_nu_floor":     s / NEUTRINO_FLOOR_CM2,
        }

    def event_rate(self,
                   detector_mass_kg: float,
                   exposure_yr: float,
                   target_amu: float | None = None,
                   A: int | None = None) -> float:
        """
        Expected events / (kg yr) * detector_mass * exposure.

        R = N_T * Phi_DM * sigma_nucleus
        with
            N_T   = detector_mass / (target_amu * AMU)   (number of nuclei)
            Phi   = (rho_DM / m_DM) * v_rel              (DM number flux)
        """
        if target_amu is None:
            target_amu = self.target_amu
        if A is None:
            A = self.target_A

        # Temporarily swap target if user passed override
        old_a, old_amu = self.target_A, self.target_amu
        self.target_A, self.target_amu = A, target_amu
        sigma_m2 = self.sigma_si_nucleus_cm2() * 1.0e-4
        self.target_A, self.target_amu = old_a, old_amu

        n_targets = detector_mass_kg / (target_amu * AMU)
        n_dm = self.rho_dm / self.m_dm_kg              # m^-3
        flux = n_dm * self.v_rel                       # m^-2 s^-1
        rate_per_s = n_targets * flux * sigma_m2
        seconds = exposure_yr * sc.year
        return rate_per_s * seconds

    def discovery_horizon(self,
                          target_events: float = 1.0,
                          exposure_yr: float = 1.0) -> float:
        """
        Detector mass (kg) needed to expect `target_events` over `exposure_yr`.
        """
        unit_rate = self.event_rate(detector_mass_kg=1.0,
                                    exposure_yr=exposure_yr)
        if unit_rate <= 0.0:
            return math.inf
        return target_events / unit_rate

    # ------------------------------------------------------------------
    # Sanity / framework consistency
    # ------------------------------------------------------------------
    def verify_null_consistency(self) -> Dict[str, object]:
        """
        Confirm the substrate prediction is well below current limits, so
        LZ/PandaX nulls are trivially consistent with the cube-cell DM model.
        """
        cmp = self.compare_to_limits()
        worst = max(cmp["ratio_to_LZ_2024"],
                    cmp["ratio_to_PandaX_2024"],
                    cmp["ratio_to_XENONnT_2027"])
        consistent = worst < 1.0
        return {
            "sigma_si_nucleon_cm2": cmp["sigma_si_nucleon_cm2"],
            "below_LZ_by_factor":   1.0 / cmp["ratio_to_LZ_2024"],
            "below_nu_floor_by":    1.0 / cmp["ratio_to_nu_floor"],
            "consistent_with_nulls": consistent,
        }

    def compare_to_wimp_scenarios(self) -> Dict[str, Dict[str, float]]:
        """
        Order-of-magnitude contrast between substrate cube DM (gravity-only)
        and standard WIMP benchmarks at the same mass (~27 GeV).

        Typical SUSY neutralino at this mass: sigma_SI ~ 1e-46 to 1e-45 cm^2
        (already excluded by current xenon experiments).
        """
        sub = self.sigma_si_nucleon_cm2()
        wimp_typical = 1.0e-45
        wimp_pessimistic = 1.0e-48
        return {
            "substrate_cube_DM": {
                "sigma_si_cm2": sub,
                "channel": "gravity-only",
                "detectable_now": sub > LIMIT_LZ_2024_CM2,
            },
            "SUSY_neutralino_typical": {
                "sigma_si_cm2": wimp_typical,
                "channel": "weak (Z/h exchange)",
                "ratio_to_substrate": wimp_typical / sub,
            },
            "SUSY_neutralino_just_below_limit": {
                "sigma_si_cm2": wimp_pessimistic,
                "channel": "weak (suppressed)",
                "ratio_to_substrate": wimp_pessimistic / sub,
            },
            "neutrino_floor": {
                "sigma_si_cm2": NEUTRINO_FLOOR_CM2,
                "ratio_to_substrate": NEUTRINO_FLOOR_CM2 / sub,
            },
        }

    # ------------------------------------------------------------------
    # Pretty summary
    # ------------------------------------------------------------------
    def summary(self) -> str:
        s = self.sigma_si_nucleon_cm2()
        sd = self.sigma_sd_nucleon_cm2()
        nuc = self.sigma_si_nucleus_cm2()
        rate_lz = self.event_rate(detector_mass_kg=7e3, exposure_yr=1.0)
        horizon = self.discovery_horizon()
        return (
            f"Cube DM gravitational cross-section\n"
            f"  m_DM            = {self.m_dm_gev:.2f} GeV\n"
            f"  sigma_SI/nuc    = {s:.3e} cm^2\n"
            f"  sigma_SD/nuc    = {sd:.3e} cm^2\n"
            f"  sigma_SI/Xe-131 = {nuc:.3e} cm^2 (coherent A^2)\n"
            f"  LZ ratio        = {s/LIMIT_LZ_2024_CM2:.3e}\n"
            f"  expected events at LZ (7 t * yr) = {rate_lz:.3e}\n"
            f"  discovery horizon (1 ev/yr)      = {horizon:.3e} kg\n"
        )


__all__ = ["CubeDMCrossSection"]
