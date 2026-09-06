"""Electromagnetic darkness gates for hybrid dark substrate stress.

The observational fact is not merely that dark matter is dim.  To be invisible
to EM instruments, a candidate must fail all ordinary EM visibility gates:

    - no photon emission from a charge-asymmetric channel;
    - no detector-resonant absorption band;
    - no appreciable reflection/scattering by EM fields.

In this substrate model, the mobile component is a heavy neutral-kink stress
with no charge-asymmetric channel, while the locked component is a coherent
long-wavelength polarization mode.  This module keeps that distinction explicit.
"""

from __future__ import annotations

from dataclasses import dataclass

from .dark_stress_hybrid import GEV_TO_G, hybrid_self_interaction
from .dark_stress_scale_closure import MYR_S, assess_dark_stress_scale_closure
from .substrate_polarization_dm import C_SI, KPC_M


H_PLANCK_J_S = 6.62607015e-34
EV_J = 1.602176634e-19
GEV_J = 1.0e9 * EV_J
FM_M = 1.0e-15
CM_M = 1.0e-2


@dataclass(frozen=True)
class EMVisibilityGate:
    """EM visibility channels for a dark-stress component."""

    emits_photons: bool
    absorbs_detector_photons: bool
    reflects_em_fields: bool
    gravitationally_visible: bool
    reason: str

    @property
    def em_dark(self) -> bool:
        return not (
            self.emits_photons
            or self.absorbs_detector_photons
            or self.reflects_em_fields
        )


@dataclass(frozen=True)
class HeavyNeutralKinkDarkness:
    """EM-darkness diagnostics for the mobile neutral-kink component."""

    mass_gev: float
    mass_kg: float
    halo_radius_fm: float
    compton_wavelength_fm: float
    halo_de_broglie_wavelength_fm: float
    compton_frequency_hz: float
    shear_mode_frequency_hz: float
    shear_mode_energy_ev: float
    gate: EMVisibilityGate
    verdict: str


@dataclass(frozen=True)
class CoherentPolarizationDarkness:
    """EM-darkness diagnostics for the locked polarization component."""

    coherence_length_kpc: float
    memory_time_myr: float
    memory_frequency_hz: float
    memory_quantum_energy_ev: float
    light_wavelength_for_memory_frequency_kpc: float
    gate: EMVisibilityGate
    verdict: str


@dataclass(frozen=True)
class EMDarknessAssessment:
    """Combined EM-darkness assessment for hybrid dark stress."""

    heavy_neutral_kink: HeavyNeutralKinkDarkness
    coherent_polarization: CoherentPolarizationDarkness
    ordinary_em_detection_expected: bool
    gravitational_detection_expected: bool
    verdict: str


def _mass_gev_to_kg(mass_gev: float) -> float:
    return mass_gev * GEV_TO_G / 1000.0


def assess_heavy_neutral_kink_darkness(
    halo_speed_km_s: float = 220.0,
) -> HeavyNeutralKinkDarkness:
    """Assess EM visibility of the mobile neutral-kink component."""

    sidm = hybrid_self_interaction()
    mass_kg = _mass_gev_to_kg(sidm.mass_gev)
    compton_m = H_PLANCK_J_S / (mass_kg * C_SI)
    de_broglie_m = H_PLANCK_J_S / (mass_kg * halo_speed_km_s * 1000.0)
    compton_frequency = sidm.mass_gev * GEV_J / H_PLANCK_J_S
    scale = assess_dark_stress_scale_closure()
    radius_m = sidm.radius_fm * FM_M
    shear_frequency = scale.speed.v_dark_km_s * 1000.0 / radius_m
    shear_energy_ev = H_PLANCK_J_S * shear_frequency / EV_J
    gate = EMVisibilityGate(
        emits_photons=False,
        absorbs_detector_photons=False,
        reflects_em_fields=False,
        gravitationally_visible=True,
        reason=(
            "neutral-kink stress carries no charge-asymmetric EM channel; "
            "its characteristic internal scales are gamma/nuclear, not broad "
            "radio/optical detector-resonant bands"
        ),
    )
    verdict = (
        "heavy neutral mobile stress: gravitationally visible and self-interacting, "
        "but EM-dark because it has no emission/absorption/reflection channel"
    )
    return HeavyNeutralKinkDarkness(
        mass_gev=sidm.mass_gev,
        mass_kg=mass_kg,
        halo_radius_fm=sidm.radius_fm,
        compton_wavelength_fm=compton_m / FM_M,
        halo_de_broglie_wavelength_fm=de_broglie_m / FM_M,
        compton_frequency_hz=compton_frequency,
        shear_mode_frequency_hz=shear_frequency,
        shear_mode_energy_ev=shear_energy_ev,
        gate=gate,
        verdict=verdict,
    )


def assess_coherent_polarization_darkness() -> CoherentPolarizationDarkness:
    """Assess EM visibility of the locked polarization component."""

    scale = assess_dark_stress_scale_closure()
    memory_s = scale.tau_pol_myr * MYR_S
    frequency = 1.0 / memory_s
    energy_ev = H_PLANCK_J_S * frequency / EV_J
    light_wavelength_kpc = C_SI / frequency / KPC_M
    gate = EMVisibilityGate(
        emits_photons=False,
        absorbs_detector_photons=False,
        reflects_em_fields=False,
        gravitationally_visible=True,
        reason=(
            "locked substrate polarization is a coherent stress field with "
            "cluster-scale wavelength and ultra-low frequency, not a photon source"
        ),
    )
    verdict = (
        "ultra-long-wavelength coherent polarization: visible through gravity, "
        "not through ordinary EM resonance"
    )
    return CoherentPolarizationDarkness(
        coherence_length_kpc=scale.coherence.ell_pol_kpc,
        memory_time_myr=scale.tau_pol_myr,
        memory_frequency_hz=frequency,
        memory_quantum_energy_ev=energy_ev,
        light_wavelength_for_memory_frequency_kpc=light_wavelength_kpc,
        gate=gate,
        verdict=verdict,
    )


def assess_em_darkness() -> EMDarknessAssessment:
    """Assess whether hybrid dark stress satisfies EM invisibility constraints."""

    heavy = assess_heavy_neutral_kink_darkness()
    coherent = assess_coherent_polarization_darkness()
    ordinary_em_detection = not (heavy.gate.em_dark and coherent.gate.em_dark)
    gravity = heavy.gate.gravitationally_visible and coherent.gate.gravitationally_visible
    verdict = (
        "hybrid dark stress satisfies the EM-darkness gate: mobile stress is "
        "heavy neutral matter-like stress, while locked polarization is an "
        "ultra-light coherent substrate mode"
        if not ordinary_em_detection and gravity
        else "hybrid dark stress leaks into an ordinary EM visibility channel"
    )
    return EMDarknessAssessment(
        heavy_neutral_kink=heavy,
        coherent_polarization=coherent,
        ordinary_em_detection_expected=ordinary_em_detection,
        gravitational_detection_expected=gravity,
        verdict=verdict,
    )
