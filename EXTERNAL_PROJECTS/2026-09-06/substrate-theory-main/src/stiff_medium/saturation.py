"""Saturation barrier — the §18.39 σ_max = ½ elastic limit.

The substrate has a maximum strain σ_max = ½. Beyond this, the medium
cannot deform elastically.

This single mechanism explains:
- Black hole horizons (where local σ → ½)
- Big Bang initial state (universe-wide σ = ½)
- No singularities (σ is capped)
- Cosmological constant bound (vacuum strain bounded)
- CMB de-saturation phase transition
"""

import numpy as np
from typing import Tuple


SIGMA_MAX = 0.5  # universal substrate elastic limit


def is_saturated(sigma: float) -> bool:
    """Check if local strain has reached saturation."""
    return sigma >= SIGMA_MAX


def saturation_barrier_potential(phi: float, phi_max: float, K: float = 1.0) -> float:
    """V(φ) = -K log(1 - (φ/φ_max)²) / 2

    Diverges as φ → φ_max, providing the saturation barrier.
    Zero at φ = 0, grows quadratically (K φ²/φ_max² approximately) for small φ.
    """
    if abs(phi) >= phi_max:
        return float('inf')
    y = phi / phi_max
    return -K * np.log(1 - y**2) / 2


def saturation_force(phi: float, phi_max: float, K: float = 1.0) -> float:
    """F = -dV/dφ for saturation-limit potential.

    F = -K φ / (φ_max² - φ²)
    """
    if abs(phi) >= phi_max:
        return -np.sign(phi) * float('inf')
    return -K * phi / (phi_max**2 - phi**2)


def schwarzschild_horizon_strain(M_kg: float, r_m: float,
                                  G: float = 6.674e-11,
                                  c: float = 2.998e8) -> float:
    """Compute σ at given r from a point mass M.

    σ(r) = GM / (r c²)

    For r = r_s = 2GM/c²: σ = ½ (universal saturation)
    """
    return G * M_kg / (r_m * c**2)


def is_inside_horizon(M_kg: float, r_m: float,
                       G: float = 6.674e-11,
                       c: float = 2.998e8) -> bool:
    """Check if (M, r) point is inside Schwarzschild horizon."""
    return schwarzschild_horizon_strain(M_kg, r_m, G, c) >= SIGMA_MAX


def schwarzschild_radius(M_kg: float, G: float = 6.674e-11,
                          c: float = 2.998e8) -> float:
    """r_s = 2GM/c² where σ = ½."""
    return 2 * G * M_kg / c**2


def time_dilation_factor(sigma: float) -> float:
    """Clock rate ratio dt/dt_∞ = √(1 - 2σ).

    σ = 0: no dilation (far from gravity)
    σ → ½: dilation factor → 0 (clocks freeze at horizon)
    """
    if sigma >= SIGMA_MAX:
        return 0.0
    return np.sqrt(1 - 2 * sigma)


def saturated_region_check_3d(rho_field: np.ndarray, dx: float,
                              K: float = 1.0,
                              G: float = 6.674e-11,
                              c: float = 2.998e8) -> np.ndarray:
    """Given a 3D mass density field, return 3D bool field of saturated regions.

    For each point, compute integrated mass within local radius and check
    if σ(local r) ≥ σ_max.
    """
    # Simplified: compute σ as ρ × dx² × G / c² × radial-factor
    # Actual computation would need spherical integral; this is order-of-magnitude
    sigma_field = G * rho_field * dx**3 / (dx * c**2)
    return sigma_field >= SIGMA_MAX


def saturation_bounded_strain(K: float, phi_max_squared: float) -> float:
    """Maximum vacuum strain energy density given saturation limit.

    From V(φ_max → 0.99 φ_max): ε ≈ K log(0.01) / 2 ≈ -K × 2.3
    But sigma_max = 0.5 is the actual cap on physical σ.

    Returns max strain energy density allowed: K × σ_max² / 2
    """
    return 0.5 * K * SIGMA_MAX**2
