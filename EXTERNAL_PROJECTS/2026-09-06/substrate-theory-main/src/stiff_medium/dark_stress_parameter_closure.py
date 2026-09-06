"""Parameter-closure candidates for the hybrid dark-stress sector.

The hybrid model in :mod:`stiff_medium.dark_stress_hybrid` needs four numbers:

    - Omega_dark / Omega_baryon
    - mobile fraction of dark stress
    - polarization memory time
    - neutral halo radius

This module tests compact substrate-geometric candidates for those numbers.
The first two are genuinely dimensionless and therefore meaningful closure
targets.  The memory time still needs a derived substrate clock, so it is
reported as a target-shape/mnemonic rather than a closed derivation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .dark_stress_hybrid import (
    ALPHA_EM,
    XI_QCD_FM,
    hybrid_cluster_check,
    hybrid_self_interaction,
)


OBS_DARK_TO_BARYON = 5.36


@dataclass(frozen=True)
class DarkAbundanceClosure:
    """Candidate for Omega_dark / Omega_baryon."""

    formula: str
    predicted_dark_to_baryon: float
    observed_dark_to_baryon: float
    error_pct: float
    mechanism: str
    verdict: str


def dark_to_baryon_phase_space_closure(
    observed_dark_to_baryon: float = OBS_DARK_TO_BARYON,
) -> DarkAbundanceClosure:
    """Omega_dark/Omega_b ~= (2*pi - 1) * (1 + 1/(8*pi^2)).

    Physical reading:
        - one bright/baryonic orientation is removed from a full U(1) angular
          phase-space measure, leaving 2*pi - 1 dark neutral orientations;
        - a one-loop substrate correction 1/(8*pi^2) appears repeatedly in
          the lepton and dark-sector trials.
    """

    predicted = (2.0 * math.pi - 1.0) * (1.0 + 1.0 / (8.0 * math.pi**2))
    err = (predicted / observed_dark_to_baryon - 1.0) * 100.0
    verdict = (
        "subpercent abundance candidate; derive phase-space measure"
        if abs(err) < 1.0
        else "miss"
    )
    return DarkAbundanceClosure(
        formula="(2*pi - 1) * (1 + 1/(8*pi^2))",
        predicted_dark_to_baryon=predicted,
        observed_dark_to_baryon=observed_dark_to_baryon,
        error_pct=err,
        mechanism=(
            "full U(1) neutral orientation phase space minus one bright "
            "orientation, with one-loop substrate correction"
        ),
        verdict=verdict,
    )


@dataclass(frozen=True)
class MobileFractionClosure:
    """Candidate for the mobile/locked split of dark substrate stress."""

    formula: str
    mobile_fraction: float
    locked_fraction: float
    minimum_required_mobile_fraction: float
    margin: float
    mechanism: str
    verdict: str


def mobile_fraction_phase_closure(dark_to_baryon: float | None = None) -> MobileFractionClosure:
    """f_mobile = 1 - 1/(2*pi), with one zero mode locked as polarization."""

    if dark_to_baryon is None:
        dark_to_baryon = dark_to_baryon_phase_space_closure().predicted_dark_to_baryon
    mobile = 1.0 - 1.0 / (2.0 * math.pi)
    locked = 1.0 - mobile
    cluster = hybrid_cluster_check(cosmic_dark_to_baryon=dark_to_baryon)
    min_req = cluster.minimum_mobile_fraction_of_dark
    margin = mobile - min_req
    verdict = (
        "passes cluster mobile-fraction threshold"
        if margin >= 0.0
        else "too locked; fails cluster threshold"
    )
    return MobileFractionClosure(
        formula="1 - 1/(2*pi)",
        mobile_fraction=mobile,
        locked_fraction=locked,
        minimum_required_mobile_fraction=min_req,
        margin=margin,
        mechanism=(
            "one U(1) angular zero mode remains locked to polarization; "
            "the remaining phase-space measure is mobile neutral stress"
        ),
        verdict=verdict,
    )


@dataclass(frozen=True)
class MemoryTimeClosure:
    """Candidate for tau_pol in Myr."""

    formula: str
    predicted_tau_myr: float
    required_tau_myr: float
    error_pct: float
    mechanism: str
    verdict: str


def polarization_memory_time_candidate(required_tau_myr: float | None = None) -> MemoryTimeClosure:
    """tau_pol ~= 4*pi^2 + 3*pi Myr.

    This is not closed because the Myr clock is not derived here.  It is a
    compact target showing that the required Bullet-like memory scale is an
    O(4*pi^2) relaxation count plus a 3-sector closure correction.
    """

    if required_tau_myr is None:
        required_tau_myr = hybrid_cluster_check().required_memory_myr
    predicted = 4.0 * math.pi**2 + 3.0 * math.pi
    err = (predicted / required_tau_myr - 1.0) * 100.0
    return MemoryTimeClosure(
        formula="(4*pi^2 + 3*pi) * tau_clock, with tau_clock provisionally 1 Myr",
        predicted_tau_myr=predicted,
        required_tau_myr=required_tau_myr,
        error_pct=err,
        mechanism=(
            "orientation-vortex relaxation count plus three-sector closure; "
            "the substrate clock tau_clock is still open"
        ),
        verdict="numerically sharp but dimensionful clock is not derived",
    )


@dataclass(frozen=True)
class HaloRadiusClosure:
    """Candidate for the neutral halo radius and self-interaction."""

    formula: str
    radius_fm: float
    sigma_over_m_cm2_g: float
    mechanism: str
    verdict: str


def halo_radius_alpha_closure() -> HaloRadiusClosure:
    """R_halo = xi_QCD / alpha."""

    radius = XI_QCD_FM / ALPHA_EM
    sidm = hybrid_self_interaction(radius_fm=radius)
    return HaloRadiusClosure(
        formula="xi_QCD / alpha",
        radius_fm=radius,
        sigma_over_m_cm2_g=sidm.sigma_over_m_cm2_g,
        mechanism="neutral core dressed by charge-symmetric polarization halo",
        verdict=sidm.verdict,
    )


@dataclass(frozen=True)
class DarkStressParameterClosure:
    """Bundle of the dark-sector parameter-closure candidates."""

    abundance: DarkAbundanceClosure
    mobile_fraction: MobileFractionClosure
    memory_time: MemoryTimeClosure
    halo_radius: HaloRadiusClosure
    verdict: str


def assess_dark_stress_parameter_closure() -> DarkStressParameterClosure:
    """Assess the compact parameter-closure candidates."""

    abundance = dark_to_baryon_phase_space_closure()
    mobile = mobile_fraction_phase_closure(abundance.predicted_dark_to_baryon)
    memory = polarization_memory_time_candidate()
    halo = halo_radius_alpha_closure()
    if (
        abs(abundance.error_pct) < 1.0
        and mobile.margin > 0.0
        and 0.1 <= halo.sigma_over_m_cm2_g <= 1.0
    ):
        verdict = (
            "dimensionless dark-sector closures are promising; memory clock "
            "and dynamical equations remain open"
        )
    else:
        verdict = "dark-sector closure fails at least one dimensionless gate"
    return DarkStressParameterClosure(
        abundance=abundance,
        mobile_fraction=mobile,
        memory_time=memory,
        halo_radius=halo,
        verdict=verdict,
    )
