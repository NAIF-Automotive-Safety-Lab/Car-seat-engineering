"""Pure substrate-polarization tests for dark-matter phenomenology.

Question tested here:

    With the extra substrate dynamics, can we remove a separate dark-matter
    component entirely and explain the data as baryon-induced polarization?

Answer from this toy module:

    - Galaxy rotation curves: yes, if the substrate susceptibility has the
      low-acceleration form g ~= sqrt(g_N a0).
    - Lensing: yes only if the polarization contributes to the same metric
      potential as mass-energy, not just to particle inertia.
    - Bullet/cluster offsets: no for an instantaneous baryon-locked field.
      Offsets require a relaxation/memory time, which is an independent dark
      stress component in practice.

This is not a full cosmological fit.  It is a boundary test for "no DM at all"
versus "no fundamental WIMP, but yes dark substrate stress".
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


G_SI = 6.67430e-11
C_SI = 2.99792458e8
M_SUN_KG = 1.98847e30
KPC_M = 3.0856775814913673e19
MPC_M = 3.0856775814913673e22
YEAR_S = 365.25 * 24.0 * 3600.0


def cosmological_a0(h0_km_s_mpc: float = 67.4) -> float:
    """Return a0 = c H0 / (2 pi), the natural cosmological acceleration."""

    h0_si = h0_km_s_mpc * 1000.0 / MPC_M
    return C_SI * h0_si / (2.0 * math.pi)


def substrate_nu_function(g_newton: float, a0: float) -> float:
    """Substrate susceptibility factor nu(g_N/a0).

    Uses the empirical radial-acceleration-style interpolation:

        g = g_N / [1 - exp(-sqrt(g_N/a0))]

    It has the desired limits:

        g >> a0:  g -> g_N
        g << a0:  g -> sqrt(g_N a0)

    In this model this is interpreted as saturated substrate gravitational
    susceptibility, not modified inertia.
    """

    if g_newton <= 0.0:
        return 0.0
    y = g_newton / a0
    denom = 1.0 - math.exp(-math.sqrt(y))
    return 1.0 / denom


def total_acceleration(g_newton: float, a0: float | None = None) -> float:
    """Return total acceleration from baryons plus substrate polarization."""

    if a0 is None:
        a0 = cosmological_a0()
    return g_newton * substrate_nu_function(g_newton, a0)


def point_mass_rotation_velocity(
    baryonic_mass_kg: float,
    radius_m: float,
    a0: float | None = None,
) -> float:
    """Circular velocity for a point baryonic mass plus polarization."""

    if a0 is None:
        a0 = cosmological_a0()
    g_newton = G_SI * baryonic_mass_kg / radius_m**2
    return math.sqrt(total_acceleration(g_newton, a0) * radius_m)


def baryonic_tully_fisher_velocity(
    baryonic_mass_kg: float,
    a0: float | None = None,
) -> float:
    """Asymptotic velocity v^4 = G M_b a0."""

    if a0 is None:
        a0 = cosmological_a0()
    return (G_SI * baryonic_mass_kg * a0) ** 0.25


def effective_mass_ratio(
    baryonic_mass_kg: float,
    radius_m: float,
    a0: float | None = None,
) -> float:
    """Return M_effective(<r) / M_b inferred from dynamics/lensing."""

    if a0 is None:
        a0 = cosmological_a0()
    g_newton = G_SI * baryonic_mass_kg / radius_m**2
    g_total = total_acceleration(g_newton, a0)
    m_eff = g_total * radius_m**2 / G_SI
    return m_eff / baryonic_mass_kg


@dataclass(frozen=True)
class RotationCurveResult:
    """Galaxy-scale result for the pure-polarization law."""

    baryonic_mass_msun: float
    radii_kpc: tuple[float, ...]
    velocities_km_s: tuple[float, ...]
    btfr_velocity_km_s: float
    flatness_fraction: float
    mass_ratio_at_outer_radius: float
    verdict: str


def rotation_curve_check(
    baryonic_mass_msun: float = 6.0e10,
    radii_kpc: Sequence[float] = (10.0, 20.0, 50.0, 100.0),
) -> RotationCurveResult:
    """Check whether pure polarization gives a flat galaxy rotation curve."""

    mass_kg = baryonic_mass_msun * M_SUN_KG
    velocities = tuple(
        point_mass_rotation_velocity(mass_kg, r * KPC_M) / 1000.0
        for r in radii_kpc
    )
    mean_v = sum(velocities) / len(velocities)
    flatness = max(abs(v - mean_v) for v in velocities) / mean_v
    outer_ratio = effective_mass_ratio(mass_kg, radii_kpc[-1] * KPC_M)
    btfr = baryonic_tully_fisher_velocity(mass_kg) / 1000.0
    if flatness < 0.15 and abs(velocities[-1] - btfr) / btfr < 0.10:
        verdict = "passes galaxy rotation and BTFR scale"
    else:
        verdict = "fails galaxy rotation flatness"
    return RotationCurveResult(
        baryonic_mass_msun=baryonic_mass_msun,
        radii_kpc=tuple(radii_kpc),
        velocities_km_s=velocities,
        btfr_velocity_km_s=btfr,
        flatness_fraction=flatness,
        mass_ratio_at_outer_radius=outer_ratio,
        verdict=verdict,
    )


@dataclass(frozen=True)
class SolarSystemResult:
    """High-acceleration sanity check."""

    g_newton: float
    g_total: float
    fractional_excess: float
    verdict: str


def solar_system_check() -> SolarSystemResult:
    """Check that the interpolation shuts off in the Solar System."""

    m_sun = M_SUN_KG
    radius_au = 1.495978707e11
    g_newton = G_SI * m_sun / radius_au**2
    g_total = total_acceleration(g_newton)
    excess = (g_total - g_newton) / g_newton
    verdict = "passes high-acceleration shutoff" if excess < 1.0e-8 else "too large"
    return SolarSystemResult(
        g_newton=g_newton,
        g_total=g_total,
        fractional_excess=excess,
        verdict=verdict,
    )


@dataclass(frozen=True)
class ClusterSeparationResult:
    """Mass/light separation requirement for a pure-polarization field."""

    offset_kpc: float
    collision_speed_km_s: float
    required_memory_myr: float
    instantaneous_local_passes: bool
    memory_turns_into_dark_stress: bool
    verdict: str


def cluster_separation_check(
    offset_kpc: float = 150.0,
    collision_speed_km_s: float = 3000.0,
) -> ClusterSeparationResult:
    """Check Bullet-like lensing separation.

    If polarization is instantaneously sourced by baryons/gas, its lensing peak
    cannot remain offset by O(100 kpc).  To keep an offset d after a collision
    with speed v, the polarization must retain memory for tau >= d/v.
    """

    required_memory_s = offset_kpc * KPC_M / (collision_speed_km_s * 1000.0)
    required_memory_myr = required_memory_s / (1.0e6 * YEAR_S)
    instantaneous_local_passes = False
    memory_turns_into_dark_stress = True
    verdict = (
        "instantaneous pure polarization fails cluster offsets; "
        "memory/propagation is an independent dark stress component"
    )
    return ClusterSeparationResult(
        offset_kpc=offset_kpc,
        collision_speed_km_s=collision_speed_km_s,
        required_memory_myr=required_memory_myr,
        instantaneous_local_passes=instantaneous_local_passes,
        memory_turns_into_dark_stress=memory_turns_into_dark_stress,
        verdict=verdict,
    )


@dataclass(frozen=True)
class PurePolarizationAssessment:
    """Overall strict no-DM assessment."""

    a0_m_s2: float
    rotation: RotationCurveResult
    solar_system: SolarSystemResult
    cluster: ClusterSeparationResult
    strict_no_dm_passes: bool
    hybrid_required: bool
    verdict: str


def assess_pure_polarization() -> PurePolarizationAssessment:
    """Assess whether pure baryon-locked polarization eliminates dark matter."""

    rotation = rotation_curve_check()
    solar = solar_system_check()
    cluster = cluster_separation_check()
    strict_passes = (
        rotation.verdict.startswith("passes")
        and solar.verdict.startswith("passes")
        and cluster.instantaneous_local_passes
    )
    hybrid_required = not strict_passes
    verdict = (
        "pure instantaneous polarization can handle galaxies but fails cluster "
        "mass/light separation; use hybrid neutral stress/polarization sector"
    )
    return PurePolarizationAssessment(
        a0_m_s2=cosmological_a0(),
        rotation=rotation,
        solar_system=solar,
        cluster=cluster,
        strict_no_dm_passes=strict_passes,
        hybrid_required=hybrid_required,
        verdict=verdict,
    )
