"""Cluster-scale dynamics check for the hybrid dark-stress sector.

The parameter and scale closures make sharp claims:

    Omega_dark/Omega_b ~= 5.350
    f_mobile ~= 0.8408
    tau_pol ~= 48.77 Myr

This module asks whether those numbers form a coherent cluster picture.  The
key distinction is that the cluster-scale lensing separation must be carried by
the mobile neutral-kink component; the slower polarization response supplies
local memory and galaxy-scale susceptibility, not long-range transport across
the whole Bullet-like offset.
"""

from __future__ import annotations

from dataclasses import dataclass

from .dark_stress_hybrid import hybrid_cluster_check
from .dark_stress_parameter_closure import (
    dark_to_baryon_phase_space_closure,
    mobile_fraction_phase_closure,
)
from .dark_stress_scale_closure import MYR_S, assess_dark_stress_scale_closure
from .substrate_polarization_dm import KPC_M


@dataclass(frozen=True)
class ClusterDynamicsAssessment:
    """Cluster dynamics implied by the current dark-stress closures."""

    dark_to_baryon: float
    mobile_fraction_of_dark: float
    locked_fraction_of_dark: float
    mobile_to_baryon: float
    locked_to_baryon: float
    mobile_fraction_of_total_lensing: float
    mobile_peak_to_gas_locked_ratio: float
    target_offset_kpc: float
    predicted_memory_offset_kpc: float
    offset_error_pct: float
    tau_pol_myr: float
    ell_pol_kpc: float
    v_dark_km_s: float
    stress_horizon_kpc: float
    stress_horizon_fraction_of_offset: float
    coherence_steps_during_memory: float
    collision_coherence_steps: float
    centroid_offset_kpc: float
    passes_mobile_lensing: bool
    passes_offset_memory: bool
    polarization_alone_spans_offset: bool
    verdict: str


def assess_cluster_dynamics(
    target_offset_kpc: float = 150.0,
    collision_speed_km_s: float = 3000.0,
    required_mobile_total_fraction: float = 0.70,
) -> ClusterDynamicsAssessment:
    """Assess whether the dark closures survive a Bullet-like cluster audit."""

    abundance = dark_to_baryon_phase_space_closure()
    mobile = mobile_fraction_phase_closure(abundance.predicted_dark_to_baryon)
    scale = assess_dark_stress_scale_closure()

    dark_to_baryon = abundance.predicted_dark_to_baryon
    mobile_fraction = mobile.mobile_fraction
    locked_fraction = mobile.locked_fraction
    mobile_to_baryon = dark_to_baryon * mobile_fraction
    locked_to_baryon = dark_to_baryon * locked_fraction
    total_lensing_to_baryon = 1.0 + dark_to_baryon
    mobile_total = mobile_to_baryon / total_lensing_to_baryon
    peak_ratio = mobile_to_baryon / (1.0 + locked_to_baryon)

    tau_pol_myr = scale.tau_pol_myr
    predicted_offset = (
        collision_speed_km_s * 1000.0 * tau_pol_myr * MYR_S / KPC_M
    )
    offset_error = (predicted_offset / target_offset_kpc - 1.0) * 100.0

    stress_horizon = scale.speed.v_dark_km_s * 1000.0 * tau_pol_myr * MYR_S / KPC_M
    horizon_fraction = stress_horizon / predicted_offset
    coherence_steps = stress_horizon / scale.coherence.ell_pol_kpc
    collision_steps = predicted_offset / scale.coherence.ell_pol_kpc

    # If gas/locked polarization sit near the collision center while mobile
    # stress follows the galaxies, this is the total lensing centroid offset.
    centroid_offset = mobile_total * predicted_offset

    cluster = hybrid_cluster_check(
        cosmic_dark_to_baryon=dark_to_baryon,
        mobile_fraction_of_dark=mobile_fraction,
        required_mobile_total_fraction=required_mobile_total_fraction,
        polarization_memory_myr=tau_pol_myr,
        offset_kpc=target_offset_kpc,
        collision_speed_km_s=collision_speed_km_s,
    )

    passes_offset = abs(offset_error) < 1.0
    polarization_spans_offset = stress_horizon >= predicted_offset
    passes_mobile = (
        cluster.passes_mobile_lensing
        and mobile_total >= required_mobile_total_fraction
        and peak_ratio > 1.0
    )

    if passes_mobile and passes_offset and not polarization_spans_offset:
        verdict = (
            "hybrid cluster picture is coherent: mobile kink stress carries "
            "the large lensing separation, while polarization supplies local "
            "memory; derive the coupled transport equations next"
        )
    elif polarization_spans_offset:
        verdict = "polarization propagation is too large; model collapses toward pure modified gravity"
    elif not passes_mobile:
        verdict = "mobile neutral stress is too small to dominate cluster lensing"
    else:
        verdict = "memory clock misses the required cluster offset"

    return ClusterDynamicsAssessment(
        dark_to_baryon=dark_to_baryon,
        mobile_fraction_of_dark=mobile_fraction,
        locked_fraction_of_dark=locked_fraction,
        mobile_to_baryon=mobile_to_baryon,
        locked_to_baryon=locked_to_baryon,
        mobile_fraction_of_total_lensing=mobile_total,
        mobile_peak_to_gas_locked_ratio=peak_ratio,
        target_offset_kpc=target_offset_kpc,
        predicted_memory_offset_kpc=predicted_offset,
        offset_error_pct=offset_error,
        tau_pol_myr=tau_pol_myr,
        ell_pol_kpc=scale.coherence.ell_pol_kpc,
        v_dark_km_s=scale.speed.v_dark_km_s,
        stress_horizon_kpc=stress_horizon,
        stress_horizon_fraction_of_offset=horizon_fraction,
        coherence_steps_during_memory=coherence_steps,
        collision_coherence_steps=collision_steps,
        centroid_offset_kpc=centroid_offset,
        passes_mobile_lensing=passes_mobile,
        passes_offset_memory=passes_offset,
        polarization_alone_spans_offset=polarization_spans_offset,
        verdict=verdict,
    )
