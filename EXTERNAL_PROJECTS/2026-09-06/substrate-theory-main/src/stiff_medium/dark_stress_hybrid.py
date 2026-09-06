"""Hybrid neutral-kink / substrate-polarization dark stress.

This module continues the strict no-DM test.  Pure instantaneous polarization
can fit galaxy rotation curves but fails cluster mass/light separation.  The
minimal surviving model is therefore a hybrid:

    - a mobile neutral-kink stress component that can separate from gas;
    - a substrate polarization component that gives the galaxy acceleration law;
    - both source the same gravitational potential for lensing.

The aim here is to quantify how much of the dark stress must be mobile.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .substrate_polarization_dm import (
    KPC_M,
    YEAR_S,
    cluster_separation_check,
    rotation_curve_check,
)


ALPHA_EM = 1.0 / 137.035999084
XI_QCD_FM = 0.2
GEV_TO_G = 1.78266192e-24
FM_TO_CM = 1.0e-13


@dataclass(frozen=True)
class DarkStressEquations:
    """Symbolic form of the hybrid dark-stress equations."""

    poisson: str
    kink_continuity: str
    polarization_relaxation: str
    quasi_static_closure: str
    interpretation: str


def coupled_dark_stress_equations() -> DarkStressEquations:
    """Return the proposed coupled equations in compact form."""

    return DarkStressEquations(
        poisson="laplacian(Phi) = 4*pi*G*(rho_b + rho_kink + rho_pol)",
        kink_continuity="d_t rho_kink + div(rho_kink*v_kink) = collision/self-interaction terms",
        polarization_relaxation=(
            "tau_pol*d_t rho_pol + rho_pol - ell_pol^2*laplacian(rho_pol) "
            "= chi(g_N)*(rho_b + epsilon_k*rho_kink)"
        ),
        quasi_static_closure=(
            "M_dark_eff(<r)/M_b = nu(g_N/a0) - 1, split into mobile kink "
            "and locked polarization fractions"
        ),
        interpretation=(
            "No external WIMP sector is required, but rho_kink + rho_pol is "
            "a real dark substrate-stress source for gravity and lensing."
        ),
    )


@dataclass(frozen=True)
class DarkStressSplit:
    """Split an effective dark ratio into mobile and polarization parts."""

    effective_dark_to_baryon: float
    mobile_fraction_of_dark: float
    mobile_to_baryon: float
    polarization_to_baryon: float
    total_to_baryon: float


def split_dark_stress(
    effective_dark_to_baryon: float,
    mobile_fraction_of_dark: float,
) -> DarkStressSplit:
    """Split dark stress between mobile kink and locked polarization pieces."""

    if effective_dark_to_baryon < 0.0:
        raise ValueError("effective_dark_to_baryon must be non-negative")
    if not (0.0 <= mobile_fraction_of_dark <= 1.0):
        raise ValueError("mobile_fraction_of_dark must be in [0, 1]")
    mobile = mobile_fraction_of_dark * effective_dark_to_baryon
    pol = (1.0 - mobile_fraction_of_dark) * effective_dark_to_baryon
    return DarkStressSplit(
        effective_dark_to_baryon=effective_dark_to_baryon,
        mobile_fraction_of_dark=mobile_fraction_of_dark,
        mobile_to_baryon=mobile,
        polarization_to_baryon=pol,
        total_to_baryon=1.0 + effective_dark_to_baryon,
    )


@dataclass(frozen=True)
class HybridClusterResult:
    """Cluster mass/light separation test for the hybrid model."""

    cosmic_dark_to_baryon: float
    mobile_fraction_of_dark: float
    mobile_fraction_of_total_lensing: float
    minimum_mobile_fraction_of_dark: float
    required_memory_myr: float
    polarization_memory_myr: float
    memory_offset_kpc: float
    passes_mobile_lensing: bool
    passes_memory_offset: bool
    verdict: str


def minimum_mobile_fraction_for_lensing(
    dark_to_baryon: float = 5.36,
    required_mobile_total_fraction: float = 0.70,
) -> float:
    """Mobile fraction of dark stress needed to dominate total lensing mass."""

    if dark_to_baryon <= 0.0:
        return float("inf")
    return required_mobile_total_fraction * (1.0 + dark_to_baryon) / dark_to_baryon


def hybrid_cluster_check(
    cosmic_dark_to_baryon: float = 5.36,
    mobile_fraction_of_dark: float = 0.85,
    required_mobile_total_fraction: float = 0.70,
    polarization_memory_myr: float = 60.0,
    offset_kpc: float = 150.0,
    collision_speed_km_s: float = 3000.0,
) -> HybridClusterResult:
    """Check whether hybrid stress can support Bullet-like separation."""

    pure = cluster_separation_check(offset_kpc, collision_speed_km_s)
    mobile_total = (
        mobile_fraction_of_dark
        * cosmic_dark_to_baryon
        / (1.0 + cosmic_dark_to_baryon)
    )
    min_mobile = minimum_mobile_fraction_for_lensing(
        cosmic_dark_to_baryon,
        required_mobile_total_fraction,
    )
    memory_offset_kpc = (
        collision_speed_km_s
        * 1000.0
        * polarization_memory_myr
        * 1.0e6
        * YEAR_S
        / KPC_M
    )
    passes_mobile = mobile_total >= required_mobile_total_fraction
    passes_memory = polarization_memory_myr >= pure.required_memory_myr
    if passes_mobile and passes_memory:
        verdict = "passes cluster offset as hybrid dark substrate stress"
    elif passes_mobile:
        verdict = "mobile kink stress passes; polarization memory can be shorter"
    elif passes_memory:
        verdict = "memory passes but mobile fraction is too small"
    else:
        verdict = "fails cluster separation"
    return HybridClusterResult(
        cosmic_dark_to_baryon=cosmic_dark_to_baryon,
        mobile_fraction_of_dark=mobile_fraction_of_dark,
        mobile_fraction_of_total_lensing=mobile_total,
        minimum_mobile_fraction_of_dark=min_mobile,
        required_memory_myr=pure.required_memory_myr,
        polarization_memory_myr=polarization_memory_myr,
        memory_offset_kpc=memory_offset_kpc,
        passes_mobile_lensing=passes_mobile,
        passes_memory_offset=passes_memory,
        verdict=verdict,
    )


@dataclass(frozen=True)
class HybridGalaxyResult:
    """Galaxy rotation plus dark-stress decomposition."""

    rotation_verdict: str
    outer_effective_dark_to_baryon: float
    outer_mobile_to_baryon: float
    outer_polarization_to_baryon: float
    outer_mobile_fraction_of_total_lensing: float
    verdict: str


def hybrid_galaxy_check(mobile_fraction_of_dark: float = 0.85) -> HybridGalaxyResult:
    """Decompose the galaxy effective halo into mobile and polarization parts."""

    rotation = rotation_curve_check()
    dark_ratio = rotation.mass_ratio_at_outer_radius - 1.0
    split = split_dark_stress(dark_ratio, mobile_fraction_of_dark)
    mobile_total = split.mobile_to_baryon / split.total_to_baryon
    verdict = (
        "passes galaxy rotation with hybrid decomposition"
        if rotation.verdict.startswith("passes")
        else "fails galaxy rotation"
    )
    return HybridGalaxyResult(
        rotation_verdict=rotation.verdict,
        outer_effective_dark_to_baryon=dark_ratio,
        outer_mobile_to_baryon=split.mobile_to_baryon,
        outer_polarization_to_baryon=split.polarization_to_baryon,
        outer_mobile_fraction_of_total_lensing=mobile_total,
        verdict=verdict,
    )


@dataclass(frozen=True)
class HybridSelfInteractionResult:
    """Self-interaction estimate for the neutral polarization halo."""

    radius_fm: float
    mass_gev: float
    sigma_over_m_cm2_g: float
    verdict: str


def hybrid_self_interaction(
    mass_gev: float = 48.6,
    radius_fm: float | None = None,
) -> HybridSelfInteractionResult:
    """Use R = xi_QCD/alpha unless a radius is supplied."""

    if radius_fm is None:
        radius_fm = XI_QCD_FM / ALPHA_EM
    sigma_cm2 = math.pi * (radius_fm * FM_TO_CM) ** 2
    sigma_over_m = sigma_cm2 / (mass_gev * GEV_TO_G)
    if 0.1 <= sigma_over_m <= 1.0:
        verdict = "self-interacting but not over-collisional"
    elif sigma_over_m < 0.1:
        verdict = "mostly collisionless"
    else:
        verdict = "likely over-collisional"
    return HybridSelfInteractionResult(
        radius_fm=radius_fm,
        mass_gev=mass_gev,
        sigma_over_m_cm2_g=sigma_over_m,
        verdict=verdict,
    )


@dataclass(frozen=True)
class HybridDarkStressAssessment:
    """Overall assessment for the hybrid dark-stress route."""

    equations: DarkStressEquations
    galaxy: HybridGalaxyResult
    cluster: HybridClusterResult
    self_interaction: HybridSelfInteractionResult
    strict_no_dm: bool
    fundamental_wimp_needed: bool
    substrate_dark_stress_needed: bool
    verdict: str


def assess_hybrid_dark_stress(
    mobile_fraction_of_dark: float = 0.85,
    polarization_memory_myr: float = 60.0,
) -> HybridDarkStressAssessment:
    """Assess the neutral-kink / polarization hybrid route."""

    galaxy = hybrid_galaxy_check(mobile_fraction_of_dark)
    cluster = hybrid_cluster_check(
        mobile_fraction_of_dark=mobile_fraction_of_dark,
        polarization_memory_myr=polarization_memory_myr,
    )
    self_interaction = hybrid_self_interaction()
    passes = (
        galaxy.verdict.startswith("passes")
        and cluster.verdict.startswith("passes")
        and self_interaction.verdict == "self-interacting but not over-collisional"
    )
    verdict = (
        "hybrid route passes toy galaxy/cluster/self-interaction tests; derive dynamics next"
        if passes
        else "hybrid route still fails at least one toy test"
    )
    return HybridDarkStressAssessment(
        equations=coupled_dark_stress_equations(),
        galaxy=galaxy,
        cluster=cluster,
        self_interaction=self_interaction,
        strict_no_dm=False,
        fundamental_wimp_needed=False,
        substrate_dark_stress_needed=True,
        verdict=verdict,
    )
