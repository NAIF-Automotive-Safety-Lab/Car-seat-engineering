"""Memory-clock trials for the hybrid dark-stress sector.

The dark-stress parameter closure found a sharp dimensionless relaxation count:

    N_relax = 4*pi^2 + 3*pi ~= 48.903

The open piece is the dimensional clock tau_clock.  This module tests whether
tau_clock can come from the hybrid dark-stress dynamics rather than being
inserted as "1 Myr" by hand.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .dark_stress_hybrid import hybrid_cluster_check, hybrid_self_interaction
from .substrate_polarization_dm import KPC_M, M_SUN_KG, YEAR_S


G_SI = 6.67430e-11
MYR_S = 1.0e6 * YEAR_S


def relaxation_count() -> float:
    """Return N_relax = 4*pi^2 + 3*pi."""

    return 4.0 * math.pi**2 + 3.0 * math.pi


@dataclass(frozen=True)
class CoherenceCrossingClock:
    """Clock from crossing one polarization coherence length."""

    ell_pol_kpc: float
    velocity_km_s: float
    tau_clock_myr: float
    tau_pol_myr: float
    required_tau_myr: float
    error_pct: float
    mechanism: str
    verdict: str


def coherence_crossing_clock(
    ell_pol_kpc: float = 1.0,
    velocity_km_s: float = 1000.0,
) -> CoherenceCrossingClock:
    """tau_clock = ell_pol / v_dark."""

    tau_clock = ell_pol_kpc * KPC_M / (velocity_km_s * 1000.0) / MYR_S
    required = hybrid_cluster_check().required_memory_myr
    tau_pol = relaxation_count() * tau_clock
    err = (tau_pol / required - 1.0) * 100.0
    verdict = (
        "viable if ell_pol and v_dark derive from cluster-scale substrate modes"
        if abs(err) < 5.0
        else "wrong memory scale"
    )
    return CoherenceCrossingClock(
        ell_pol_kpc=ell_pol_kpc,
        velocity_km_s=velocity_km_s,
        tau_clock_myr=tau_clock,
        tau_pol_myr=tau_pol,
        required_tau_myr=required,
        error_pct=err,
        mechanism="tau_clock = ell_pol / v_dark, with N_relax = 4*pi^2 + 3*pi",
        verdict=verdict,
    )


@dataclass(frozen=True)
class SelfInteractionClock:
    """Clock from neutral-stress self-interaction mean-free time."""

    density_kg_m3: float
    velocity_km_s: float
    sigma_over_m_m2_kg: float
    tau_clock_myr: float
    tau_pol_myr: float
    required_tau_myr: float
    density_required_kg_m3: float
    density_required_msun_kpc3: float
    verdict: str


def self_interaction_clock(
    density_kg_m3: float = 1.0e-22,
    velocity_km_s: float = 3000.0,
) -> SelfInteractionClock:
    """tau_clock = 1 / (rho * sigma/m * v).

    This tests whether the memory time could be a collisional relaxation time.
    """

    sidm = hybrid_self_interaction()
    sigma_over_m_m2_kg = sidm.sigma_over_m_cm2_g * 0.1
    v_m_s = velocity_km_s * 1000.0
    tau_clock_s = 1.0 / (density_kg_m3 * sigma_over_m_m2_kg * v_m_s)
    tau_clock_myr = tau_clock_s / MYR_S
    required = hybrid_cluster_check().required_memory_myr
    tau_pol = relaxation_count() * tau_clock_myr

    # If the collisional time itself is to supply the required polarization
    # memory, rho_req follows from tau_required = 1/(rho sigma/m v).
    density_required = 1.0 / (required * MYR_S * sigma_over_m_m2_kg * v_m_s)
    density_required_msun_kpc3 = density_required * KPC_M**3 / M_SUN_KG

    if tau_clock_myr > 100.0:
        verdict = "too slow at typical cluster densities; not the memory clock"
    elif density_required_msun_kpc3 > 1.0e8:
        verdict = "requires dense core; not generic"
    else:
        verdict = "possible collisional clock"
    return SelfInteractionClock(
        density_kg_m3=density_kg_m3,
        velocity_km_s=velocity_km_s,
        sigma_over_m_m2_kg=sigma_over_m_m2_kg,
        tau_clock_myr=tau_clock_myr,
        tau_pol_myr=tau_pol,
        required_tau_myr=required,
        density_required_kg_m3=density_required,
        density_required_msun_kpc3=density_required_msun_kpc3,
        verdict=verdict,
    )


@dataclass(frozen=True)
class FreeFallClock:
    """Clock from gravitational dynamical/free-fall time."""

    density_kg_m3: float
    tau_ff_myr: float
    required_density_kg_m3: float
    required_density_msun_kpc3: float
    verdict: str


def free_fall_clock(density_kg_m3: float = 1.0e-21) -> FreeFallClock:
    """Free-fall time t_ff = sqrt(3*pi/(32*G*rho))."""

    tau_ff_s = math.sqrt(3.0 * math.pi / (32.0 * G_SI * density_kg_m3))
    tau_ff_myr = tau_ff_s / MYR_S
    required = hybrid_cluster_check().required_memory_myr
    req_s = required * MYR_S
    required_density = 3.0 * math.pi / (32.0 * G_SI * req_s**2)
    required_density_msun_kpc3 = required_density * KPC_M**3 / M_SUN_KG
    verdict = (
        "possible if polarization coherence is set by cluster-core dynamical density"
        if 1.0e7 <= required_density_msun_kpc3 <= 1.0e8
        else "density scale not plausible"
    )
    return FreeFallClock(
        density_kg_m3=density_kg_m3,
        tau_ff_myr=tau_ff_myr,
        required_density_kg_m3=required_density,
        required_density_msun_kpc3=required_density_msun_kpc3,
        verdict=verdict,
    )


@dataclass(frozen=True)
class MemoryClockAssessment:
    """Bundle of memory-clock trials."""

    relaxation_count: float
    coherence_clock: CoherenceCrossingClock
    self_interaction_clock: SelfInteractionClock
    free_fall_clock: FreeFallClock
    verdict: str


def assess_memory_clock() -> MemoryClockAssessment:
    """Assess candidate dimensional clocks for tau_pol."""

    coherence = coherence_crossing_clock()
    collisional = self_interaction_clock()
    freefall = free_fall_clock()
    if abs(coherence.error_pct) < 5.0 and collisional.verdict.startswith("too slow"):
        verdict = (
            "best clock is a kpc-scale polarization coherence crossing; "
            "self-interaction mean-free time is too slow"
        )
    else:
        verdict = "memory-clock mechanism remains unresolved"
    return MemoryClockAssessment(
        relaxation_count=relaxation_count(),
        coherence_clock=coherence,
        self_interaction_clock=collisional,
        free_fall_clock=freefall,
        verdict=verdict,
    )
