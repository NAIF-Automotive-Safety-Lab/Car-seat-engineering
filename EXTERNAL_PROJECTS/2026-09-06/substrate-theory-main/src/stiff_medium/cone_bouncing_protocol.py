"""cone_bouncing_protocol.py
================================

Canonical reference for the cone-bouncing mass mechanism (MODEL.md §2.5).

Background
----------
Three modules in the codebase currently use *different* expressions for the
cone-bouncing frequency ω_b that sets a particle's rest mass via

    m c² = ℏ ω_b .

Audit (audit_03_drag_mass.md) flagged the inconsistency:

    A.  drag_mass_generator.py
        ω_b = (c/ξ) · √(1 + α · γ̃),     γ̃ = γ ξ / √(K ρ)
        — treats drag as a sub-leading multiplicative dressing of the
          Compton frequency c/ξ.  Reduces to ω₀ = c/ξ at γ → 0.

    B.  mass_torque_engine.py
        ω_b = γ · ξ · ρ · √K · f_int
        — *linear* in γ.  Vanishes at γ = 0 (no Compton baseline).
          Dimensionally only consistent in a non-SI "natural-baseline"
          unit system where K = ρ = ξ = 1.

    D.  lattice_substrate_2d.py
        ω_b² = ω₀² + (γ / 2ρ)²
        — exact damped-oscillator frequency-shift identity from the
          substrate Lagrangian.  Reduces to ω₀ at γ → 0.

This module derives ω_b directly from the substrate Lagrangian

    L = ½ ρ (∂_t u)² − ½ K |∇u|² − V(u) − γ u · ∂_t u

and provides ``omega_b_canonical`` as the single source of truth.

Derivation
----------
For a localized bound mode of size ξ, the Euler–Lagrange equation in the
small-amplitude (V''(0) > 0) limit reduces to

    ρ ü + γ u̇ + K k² u = 0,        k ≡ 1/ξ.

This is a damped harmonic oscillator with bare frequency

    ω₀² = K k² / ρ = (K / ρ) · (1/ξ²) = (c/ξ)²,    c ≡ √(K/ρ),

and damping rate Γ = γ / (2ρ).  The *peak* (resonance) frequency of the
amplitude response — i.e. the observable cone-bouncing frequency — is

    ω_b² = ω₀² + Γ²
         = (c/ξ)² + (γ / 2ρ)² .                       (CANONICAL)

This is **Formula D**.  Notice:

    ω_b = (c/ξ) · √( 1 + ( γ ξ / (2 ρ c) )² )
        = (c/ξ) · √( 1 + (γ̃/2)² · (1/(Kρ)/(Kρ)) ... )

Re-expressing in the dimensionless γ̃ ≡ γξ/√(Kρ):

    γ̃ / 2 = γ ξ / (2 √(K ρ)) = (γ / 2ρ) · (ξ √ρ / √K) = (γ / 2ρ) / (c/ξ)

so

    ω_b = (c/ξ) · √( 1 + (γ̃/2)² ).                     (canonical, γ̃ form)

This shows Formula A is the *linearization* of canonical at small γ̃ if
α = (1/2)·(γ̃/2) — i.e. A's α drift is the missing quadratic.  Formula B
is recovered only as the dominant-drag (γ̃ ≫ 1) limit *and only* in
unit-baseline natural units; it cannot be the canonical form because it
has no Compton baseline.

Anchor
------
At the electron Compton anchor ξ = ξ_e ≡ ℏ/(m_e c) and γ = 0:

    ω_b = c / ξ_e = m_e c² / ℏ
    m   = ℏ ω_b / c² = m_e          ✓

Drag (γ > 0) only *adds* mass — never subtracts.  This matches the
"smaller Q ⇒ larger mass" statement in MODEL.md §2.5.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any

import numpy as np
from scipy import constants as sc


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

C       = sc.c                # m/s
HBAR    = sc.hbar             # J·s
M_E     = sc.m_e              # kg
EV      = sc.e                # J/eV
MEV     = 1.0e6 * EV          # J/MeV

# Electron Compton wavelength λ_C / (2π) = ℏ / (m_e c)
XI_E    = HBAR / (M_E * C)    # m  ≈ 3.8616e-13


# ---------------------------------------------------------------------------
# Canonical formula
# ---------------------------------------------------------------------------

def omega_b_canonical(
    K: float,
    rho: float,
    xi: float,
    gamma: float,
    config: Optional[Dict[str, Any]] = None,
) -> float:
    """Canonical cone-bouncing angular frequency.

    Derived from L = ½ρ u̇² − ½K|∇u|² − V(u) − γ u u̇ as the resonance
    frequency of the damped bound mode of size ξ:

        ω_b² = (c/ξ)² + (γ / 2ρ)²,     c ≡ √(K/ρ).

    Parameters
    ----------
    K, rho : float
        Substrate stiffness and density.  Must be > 0.  Together set the
        substrate wave speed c_eff = √(K/ρ).  At the canonical anchor we
        choose K, ρ so that c_eff = c (speed of light); this is the
        substrate-rest-frame convention.
    xi : float
        Bound-state size (Compton wavelength of the particle).
    gamma : float
        Drag coefficient (≥ 0).  Adds mass via Γ = γ/(2ρ).
    config : dict, optional
        Override any of the above with keys ``K``, ``rho``, ``xi``,
        ``gamma``.

    Returns
    -------
    float
        ω_b in rad/s (when inputs are in SI).
    """
    if config:
        K     = float(config.get("K",     K))
        rho   = float(config.get("rho",   rho))
        xi    = float(config.get("xi",    xi))
        gamma = float(config.get("gamma", gamma))

    if K <= 0 or rho <= 0 or xi <= 0:
        raise ValueError("K, rho, xi must be strictly positive.")
    if gamma < 0:
        raise ValueError("gamma must be non-negative (drag adds mass).")

    c_eff   = np.sqrt(K / rho)
    omega_0 = c_eff / xi
    Gamma   = gamma / (2.0 * rho)
    return float(np.sqrt(omega_0 * omega_0 + Gamma * Gamma))


def mass_from_omega(omega_b: float) -> float:
    """Rest mass m = ℏ ω_b / c²  (kg, SI)."""
    return float(HBAR * omega_b / (C * C))


def dimensionless_drag_to_mass(mass_kg: float, gamma: float) -> float:
    """The audit-required O(1) sanity ratio  m c² / (ℏ γ).

    For drag-dominated states this approaches 1/2 (Γ = γ/2ρ, ρ chosen so
    c_eff = c, ξ small enough that ω₀ ≪ Γ).  For Compton-dominated
    states (γ → 0) it diverges — by design — because the mass is set
    entirely by the geometric ξ, not by γ.  We return the value as-is;
    callers test it in the regime they're interested in.
    """
    if gamma <= 0:
        return float("inf")
    return float(mass_kg * C * C / (HBAR * gamma))


# ---------------------------------------------------------------------------
# Cross-validation against the three legacy formulas
# ---------------------------------------------------------------------------

def omega_A(K: float, rho: float, xi: float, gamma: float,
            alpha: float = 1.0) -> float:
    """Formula A (drag_mass_generator):  ω = (c/ξ)·√(1 + α γ̃)."""
    c_eff       = np.sqrt(K / rho)
    gamma_tilde = gamma * xi / np.sqrt(K * rho)
    return float((c_eff / xi) * np.sqrt(1.0 + alpha * gamma_tilde))


def omega_B(K: float, rho: float, xi: float, gamma: float,
            f_int: float = 1.0) -> float:
    """Formula B (mass_torque_engine):  ω = γ ξ ρ √K · f_int."""
    return float(gamma * xi * rho * np.sqrt(K) * f_int)


def omega_D(K: float, rho: float, xi: float, gamma: float) -> float:
    """Formula D (lattice_substrate_2d):  ω² = ω₀² + (γ/2ρ)²."""
    omega_0 = np.sqrt(K / rho) / xi
    shift   = gamma / (2.0 * rho)
    return float(np.sqrt(omega_0 * omega_0 + shift * shift))


@dataclass
class CrossCompareRow:
    label: str
    omega_b: float
    mass_kg: float
    rel_err_vs_canonical: float


def cross_compare(K: float, rho: float, xi: float, gamma: float,
                  alpha: float = 1.0, f_int: float = 1.0) -> Dict[str, CrossCompareRow]:
    """Tabulate ω_b from {canonical, A, B, D} for given primitives.

    Returns a dict label → CrossCompareRow.  Useful for the test report.
    """
    can = omega_b_canonical(K, rho, xi, gamma)

    def row(lbl: str, omega: float) -> CrossCompareRow:
        rel = abs(omega - can) / can if can else float("inf")
        return CrossCompareRow(lbl, omega, mass_from_omega(omega), rel)

    return {
        "canonical": row("canonical (ω₀² + Γ², MODEL §2.5)", can),
        "A":         row("A  drag_mass_generator   √(1+αγ̃)", omega_A(K, rho, xi, gamma, alpha)),
        "B":         row("B  mass_torque_engine    γξρ√K·f", omega_B(K, rho, xi, gamma, f_int)),
        "D":         row("D  lattice_substrate_2d  ω₀²+Γ²", omega_D(K, rho, xi, gamma)),
    }


# ---------------------------------------------------------------------------
# Anchor verification utility
# ---------------------------------------------------------------------------

def electron_anchor_primitives() -> Dict[str, float]:
    """Canonical primitives that put the electron exactly on anchor.

    Choose K, ρ so that c_eff = √(K/ρ) = c, and ξ = ξ_e = ℏ/(m_e c).
    With γ = 0 this gives ω_b = m_e c² / ℏ exactly.
    """
    rho = 1.0           # arbitrary units; only ratio K/ρ matters
    K   = rho * C * C   # ⇒ c_eff = c
    return {"K": K, "rho": rho, "xi": XI_E, "gamma": 0.0}


def predicted_electron_mass_mev(primitives: Optional[Dict[str, float]] = None) -> float:
    """Convenience: predicted electron mass in MeV/c² at given primitives."""
    p = primitives or electron_anchor_primitives()
    omega = omega_b_canonical(**p)
    m_kg  = mass_from_omega(omega)
    return float(m_kg * C * C / MEV)


# ---------------------------------------------------------------------------
# Documentation report
# ---------------------------------------------------------------------------

CANONICAL_VERDICT = """\
Canonical formula:  ω_b² = (c/ξ)² + (γ/2ρ)²

Most consistent legacy module:  D  (lattice_substrate_2d)
  — identical to canonical; derivable directly from the Euler–Lagrange
    equation of the substrate Lagrangian.

Formula A (drag_mass_generator):
  agrees with canonical at γ = 0 (both → c/ξ); deviates at finite γ
  because A is linear-in-γ̃ under the radical while canonical is
  quadratic-in-Γ.  A is a small-drag *linearization*, not exact.

Formula B (mass_torque_engine):
  disagrees with canonical at γ = 0 (B → 0, canonical → c/ξ).  Cannot
  reproduce the electron mass without the Compton baseline.  Only valid
  as a drag-dominated limit in unit-baseline natural units.
"""


if __name__ == "__main__":
    p = electron_anchor_primitives()
    omega = omega_b_canonical(**p)
    m_mev = predicted_electron_mass_mev(p)
    print(f"ω_b (anchor)  = {omega:.6e} rad/s")
    print(f"m_e predicted = {m_mev:.6f} MeV/c² (target 0.510999)")
    print()
    print(CANONICAL_VERDICT)
