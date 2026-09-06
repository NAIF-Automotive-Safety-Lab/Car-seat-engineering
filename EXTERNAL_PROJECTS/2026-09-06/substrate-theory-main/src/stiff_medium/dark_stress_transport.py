"""Finite-speed 1D transport toy model for hybrid dark stress.

This is not a full cluster simulation.  It is a consistency harness for the
coupled picture left open by the dark-sector audits:

    - baryonic gas is near the collision center after a Bullet-like event;
    - mobile neutral-kink stress follows the galaxies to the offset peak;
    - locked polarization can only spread within v_dark*tau_pol.

The goal is to prevent the model from hiding the cluster offset inside an
instantaneous polarization field.  If the total lensing peak is offset, it must
be because the mobile component is large enough.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .dark_stress_cluster_dynamics import assess_cluster_dynamics


@dataclass(frozen=True)
class TransportProfiles:
    """One-dimensional post-collision stress profiles."""

    x_kpc: np.ndarray
    baryon: np.ndarray
    mobile_kink: np.ndarray
    locked_polarization: np.ndarray
    total_lensing: np.ndarray
    offset_kpc: float
    polarization_horizon_kpc: float


@dataclass(frozen=True)
class TransportAssessment:
    """Summary of the 1D transport consistency check."""

    baryon_mass: float
    mobile_mass: float
    locked_mass: float
    total_lensing_mass: float
    mobile_peak_kpc: float
    baryon_peak_kpc: float
    polarization_peak_kpc: float
    total_peak_kpc: float
    total_peak_offset_error_kpc: float
    central_total_density: float
    mobile_peak_total_density: float
    mobile_to_central_peak_ratio: float
    polarization_density_at_mobile_peak: float
    polarization_leakage_fraction: float
    verdict: str


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    """Compatibility wrapper for NumPy trapezoid integration."""

    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def _normalize(profile: np.ndarray, x: np.ndarray, mass: float) -> np.ndarray:
    area = _trapz(profile, x)
    if area <= 0.0:
        raise ValueError("profile area must be positive")
    return profile * (mass / area)


def gaussian_profile(
    x_kpc: np.ndarray,
    center_kpc: float,
    sigma_kpc: float,
    mass: float,
) -> np.ndarray:
    """Area-normalized Gaussian profile."""

    if sigma_kpc <= 0.0:
        raise ValueError("sigma_kpc must be positive")
    profile = np.exp(-0.5 * ((x_kpc - center_kpc) / sigma_kpc) ** 2)
    return _normalize(profile, x_kpc, mass)


def compact_polarization_profile(
    x_kpc: np.ndarray,
    center_kpc: float,
    horizon_kpc: float,
    mass: float,
) -> np.ndarray:
    """Compact triangular finite-speed polarization profile."""

    if horizon_kpc <= 0.0:
        raise ValueError("horizon_kpc must be positive")
    distance = np.abs(x_kpc - center_kpc)
    profile = np.maximum(1.0 - distance / horizon_kpc, 0.0)
    return _normalize(profile, x_kpc, mass)


def post_collision_profiles(
    x_min_kpc: float = -250.0,
    x_max_kpc: float = 250.0,
    dx_kpc: float = 0.5,
    baryon_sigma_kpc: float = 25.0,
    mobile_sigma_kpc: float = 25.0,
) -> TransportProfiles:
    """Build post-collision baryon, mobile-kink, and polarization profiles."""

    if dx_kpc <= 0.0:
        raise ValueError("dx_kpc must be positive")
    dynamics = assess_cluster_dynamics()
    x = np.arange(x_min_kpc, x_max_kpc + 0.5 * dx_kpc, dx_kpc)
    baryon = gaussian_profile(x, 0.0, baryon_sigma_kpc, 1.0)
    mobile = gaussian_profile(
        x,
        dynamics.predicted_memory_offset_kpc,
        mobile_sigma_kpc,
        dynamics.mobile_to_baryon,
    )
    locked = compact_polarization_profile(
        x,
        0.0,
        dynamics.stress_horizon_kpc,
        dynamics.locked_to_baryon,
    )
    return TransportProfiles(
        x_kpc=x,
        baryon=baryon,
        mobile_kink=mobile,
        locked_polarization=locked,
        total_lensing=baryon + mobile + locked,
        offset_kpc=dynamics.predicted_memory_offset_kpc,
        polarization_horizon_kpc=dynamics.stress_horizon_kpc,
    )


def profile_mass(profile: np.ndarray, x_kpc: np.ndarray) -> float:
    """Return profile area."""

    return _trapz(profile, x_kpc)


def peak_location(profile: np.ndarray, x_kpc: np.ndarray) -> float:
    """Return x coordinate of profile maximum."""

    return float(x_kpc[int(np.argmax(profile))])


def assess_transport_profiles() -> TransportAssessment:
    """Assess whether finite-speed profiles preserve the cluster-lensing offset."""

    profiles = post_collision_profiles()
    x = profiles.x_kpc

    baryon_mass = profile_mass(profiles.baryon, x)
    mobile_mass = profile_mass(profiles.mobile_kink, x)
    locked_mass = profile_mass(profiles.locked_polarization, x)
    total_mass = profile_mass(profiles.total_lensing, x)

    mobile_peak = peak_location(profiles.mobile_kink, x)
    baryon_peak = peak_location(profiles.baryon, x)
    pol_peak = peak_location(profiles.locked_polarization, x)
    total_peak = peak_location(profiles.total_lensing, x)
    peak_error = total_peak - profiles.offset_kpc

    central_idx = int(np.argmin(np.abs(x)))
    mobile_idx = int(np.argmin(np.abs(x - mobile_peak)))
    central_density = float(profiles.total_lensing[central_idx])
    mobile_density = float(profiles.total_lensing[mobile_idx])
    peak_ratio = mobile_density / central_density
    pol_at_mobile = float(profiles.locked_polarization[mobile_idx])

    outside_horizon = np.abs(x) > profiles.polarization_horizon_kpc
    leakage = profile_mass(profiles.locked_polarization[outside_horizon], x[outside_horizon])
    leakage_fraction = leakage / locked_mass if locked_mass > 0.0 else 0.0

    if abs(peak_error) < 2.0 and peak_ratio > 2.0 and leakage_fraction < 1.0e-6:
        verdict = (
            "finite-speed hybrid transport keeps the lensing maximum on the "
            "mobile neutral-stress peak; polarization remains local"
        )
    elif leakage_fraction >= 1.0e-6:
        verdict = "polarization leaks outside its finite-speed horizon"
    elif peak_ratio <= 2.0:
        verdict = "mobile peak does not dominate central gas plus polarization"
    else:
        verdict = "total lensing peak misses the mobile offset"

    return TransportAssessment(
        baryon_mass=baryon_mass,
        mobile_mass=mobile_mass,
        locked_mass=locked_mass,
        total_lensing_mass=total_mass,
        mobile_peak_kpc=mobile_peak,
        baryon_peak_kpc=baryon_peak,
        polarization_peak_kpc=pol_peak,
        total_peak_kpc=total_peak,
        total_peak_offset_error_kpc=peak_error,
        central_total_density=central_density,
        mobile_peak_total_density=mobile_density,
        mobile_to_central_peak_ratio=peak_ratio,
        polarization_density_at_mobile_peak=pol_at_mobile,
        polarization_leakage_fraction=leakage_fraction,
        verdict=verdict,
    )
