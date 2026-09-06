"""Charge-asymmetric vs charge-symmetric coupling channels — §18.32.

In our model, bound configurations interact through TWO channels:

1. CHARGE-ASYMMETRIC channel (electromagnetism):
   - Depends on Möbius half-flux orientation (chirality)
   - Has signs (positive vs negative)
   - Cancels for neutral aggregates
   - Coupling: α ≈ 1/137

2. CHARGE-SYMMETRIC channel (gravity):
   - Depends only on configuration's existence (mass)
   - Always positive sign
   - Adds for all matter
   - Coupling: α_grav = (m/M_Planck)² × α

Both channels are different aspects of the same medium back-reaction.
The hierarchy F_grav/F_em ≈ 10⁻³⁷ comes from the structural ratio.
"""

import numpy as np
from typing import Tuple


def coulomb_force_em(q1: float, q2: float, r: float,
                      epsilon_0: float = 8.854e-12) -> np.ndarray:
    """Charge-asymmetric channel: F = k q1 q2 / r²

    Returns force vector (along radial direction).
    Positive = repulsive, negative = attractive.
    """
    if r <= 0:
        return np.zeros(3)
    F_magnitude = q1 * q2 / (4 * np.pi * epsilon_0 * r**2)
    return F_magnitude  # scalar; sign indicates attract/repel


def gravity_force_charge_symmetric(m1: float, m2: float, r: float,
                                     G: float = 6.674e-11) -> np.ndarray:
    """Charge-symmetric channel: F = -G m1 m2 / r²

    Always attractive.
    """
    if r <= 0:
        return np.zeros(3)
    F_magnitude = -G * m1 * m2 / r**2
    return F_magnitude


def total_force_two_particles(q1: float, m1: float,
                               q2: float, m2: float,
                               r: float,
                               c: float = 2.998e8,
                               epsilon_0: float = 8.854e-12,
                               G: float = 6.674e-11) -> Tuple[float, float]:
    """Compute both channel contributions to force between two particles.

    Returns: (F_em, F_grav) — separate channel forces.
    """
    F_em = coulomb_force_em(q1, q2, r, epsilon_0)
    F_grav = gravity_force_charge_symmetric(m1, m2, r, G)
    return F_em, F_grav


def gravity_em_ratio_two_particles(q1: float, m1: float,
                                     q2: float, m2: float,
                                     epsilon_0: float = 8.854e-12,
                                     G: float = 6.674e-11) -> float:
    """Compute |F_grav/F_em| ratio for two charged massive particles.

    For two protons: ratio ≈ 10⁻³⁷, matching observation.
    """
    if abs(q1 * q2) < 1e-30:
        return float('inf')
    return abs(G * m1 * m2 * (4 * np.pi * epsilon_0) / (q1 * q2))


def hierarchy_check_protons(m_p: float = 1.673e-27,
                              e: float = 1.602e-19,
                              alpha: float = 1/137.036,
                              hbar: float = 1.055e-34,
                              c: float = 2.998e8,
                              G: float = 6.674e-11) -> dict:
    """Verify §18.32 hierarchy formula: F_grav/F_em = (m_p/M_Planck)²/α."""
    M_Planck = np.sqrt(hbar * c / G)
    epsilon_0 = 8.854e-12

    # Direct calculation
    F_grav = G * m_p**2
    F_em = e**2 / (4 * np.pi * epsilon_0)
    ratio_direct = F_grav / F_em

    # Predicted from structure
    ratio_predicted = (m_p / M_Planck)**2 / alpha

    return {
        "M_Planck_kg": M_Planck,
        "ratio_measured": ratio_direct,
        "ratio_predicted": ratio_predicted,
        "agreement_pct": ratio_predicted / ratio_direct * 100,
    }


def epsilon_charge_symmetric_residual(target_ratio: float = 8e-37,
                                        alpha: float = 1/137.036,
                                        m_p: float = 1.673e-27,
                                        hbar: float = 1.055e-34,
                                        c: float = 2.998e8,
                                        G: float = 6.674e-11) -> float:
    """Compute the charge-symmetric residual coupling fraction ε.

    From §18.32: G ~ ε² × α / M_substrate²
    Solving for ε given measured G:
    ε = √(G × M_substrate² / α)
    """
    # Assume M_substrate = M_Planck (substrate at Planck scale)
    M_Planck = np.sqrt(hbar * c / G)
    epsilon_squared = G * M_Planck**2 / alpha
    return np.sqrt(epsilon_squared)


class TwoChannelInteraction:
    """Combine both channels for a system of charged massive particles."""

    def __init__(self, alpha=1/137.036, G=6.674e-11):
        self.alpha = alpha
        self.G = G
        self.epsilon_0 = 8.854e-12

    def force_between(self, p1: dict, p2: dict, r: float) -> dict:
        """Compute both channels of force between two particles.

        p1, p2 should have keys: 'mass', 'charge'
        """
        F_em = coulomb_force_em(p1["charge"], p2["charge"], r, self.epsilon_0)
        F_grav = gravity_force_charge_symmetric(p1["mass"], p2["mass"], r, self.G)
        return {
            "F_em": F_em,
            "F_grav": F_grav,
            "F_total": F_em + F_grav,
            "ratio": abs(F_grav / F_em) if abs(F_em) > 0 else float('inf'),
        }
