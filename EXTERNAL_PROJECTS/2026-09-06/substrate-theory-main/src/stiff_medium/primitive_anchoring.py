"""Substrate primitive anchoring — fits (K, ρ, ξ, γ) to known constants.

The B3 / stiff-medium ansatz posits four substrate primitives:
    K   — bulk modulus (Pa = J/m³)
    ρ   — mass density (kg/m³)
    ξ   — coherence length (m)
    γ   — inverse coherence time (1/s)

These are anchored via two universal identities:
    c    = sqrt(K/ρ)              wave speed
    ℏ    = K ξ⁴ / c               action quantum (MODEL.md §1.4)

Two more relations are needed to close the system; they come from the
choice of length anchor (ξ ~ λ_C of some reference particle, or ℓ_Planck)
and a temporal anchor for γ (γ ~ c/ξ by default — single-cell light crossing).

This module compares three anchor choices for self-consistency:
    1. ELECTRON_COMPTON   ξ = ℏ/(m_e c)
    2. PROTON_COMPTON     ξ = ℏ/(m_p c)
    3. PLANCK_LENGTH      ξ = sqrt(ℏ G / c³)

For each anchor, derived quantities are computed and compared to:
    α               = 11/(48 π³) · exp(-3π/737)   [B3 fine-structure formula]
    hierarchy       = exp(4π² - 1)                 [B3 hierarchy magnitude]
    m_p / m_e       = 1836.15267343                [PDG]
    m_kink          = 8 ℏ / (c ξ)                  [sine-Gordon, MeV scale at λ_C]

The "best" anchor is the one whose derived constants reproduce observed
values within the tightest residuals.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np
import scipy.constants as sc

# -----------------------------------------------------------------------------
# Physical constants (CODATA via scipy.constants)
# -----------------------------------------------------------------------------
C_LIGHT: float = sc.c                        # m/s
HBAR: float = sc.hbar                        # J·s
G_NEWTON: float = sc.G                       # m³/(kg·s²)
M_E: float = sc.m_e                          # kg
M_P: float = sc.m_p                          # kg
ALPHA_OBS: float = sc.alpha                  # ≈ 7.2973525693e-3
E_CHARGE: float = sc.e                       # C
EPS0: float = sc.epsilon_0                   # F/m

# Reference-particle Compton wavelengths (REDUCED: λ̄_C = ℏ / (m c))
LAMBDA_C_E: float = HBAR / (M_E * C_LIGHT)   # ≈ 3.8616e-13 m
LAMBDA_C_P: float = HBAR / (M_P * C_LIGHT)   # ≈ 2.1031e-16 m
ELL_PLANCK: float = math.sqrt(HBAR * G_NEWTON / C_LIGHT ** 3)

# B3 derived targets
ALPHA_B3: float = 11.0 / (48.0 * math.pi ** 3) * math.exp(-3.0 * math.pi / 737.0)
HIERARCHY_B3: float = math.exp(4.0 * math.pi ** 2 - 1.0)
MP_OVER_ME: float = M_P / M_E

ANCHOR_CHOICES: Tuple[str, ...] = (
    "electron_compton",
    "proton_compton",
    "planck_length",
    "electron_drag",
)

# Map for drag-mass anchoring: name -> rest mass in kg
PARTICLE_MASS_KG: Dict[str, float] = {
    "electron": M_E,
    "proton": M_P,
}


# -----------------------------------------------------------------------------
# Result container
# -----------------------------------------------------------------------------
@dataclass
class AnchorResult:
    """Holds primitives + cross-check residuals for one anchor choice."""

    anchor: str
    xi: float
    K: float
    rho: float
    gamma: float
    derived: Dict[str, float] = field(default_factory=dict)
    observed: Dict[str, float] = field(default_factory=dict)
    residuals: Dict[str, float] = field(default_factory=dict)

    def total_residual(self) -> float:
        """Sum of |relative residual| across all cross-checks."""
        return float(sum(abs(r) for r in self.residuals.values()))


# -----------------------------------------------------------------------------
# Core class
# -----------------------------------------------------------------------------
class PrimitiveAnchoring:
    """Solve and cross-check substrate primitives (K, ρ, ξ, γ) per anchor."""

    def __init__(self, anchor: str = "electron_compton") -> None:
        if anchor not in ANCHOR_CHOICES:
            raise ValueError(
                f"unknown anchor {anchor!r}; choose one of {ANCHOR_CHOICES}"
            )
        self.anchor = anchor
        self.xi: float = 0.0
        self.K: float = 0.0
        self.rho: float = 0.0
        self.gamma: float = 0.0
        self._solve()

    # -- solving ---------------------------------------------------------------
    def _xi_from_anchor(self) -> float:
        if self.anchor == "electron_compton":
            return LAMBDA_C_E
        if self.anchor == "proton_compton":
            return LAMBDA_C_P
        if self.anchor == "planck_length":
            return ELL_PLANCK
        if self.anchor == "electron_drag":
            # Same length scale as electron_compton; γ is reset by drag-mass
            # solver below.
            return LAMBDA_C_E
        raise AssertionError(self.anchor)

    def _solve(self) -> None:
        """Solve for (K, ρ, ξ, γ) given the anchor + (c, ℏ).

        Identities used:
            ξ      from anchor
            K      = ℏ c / ξ⁴      from   ℏ = K ξ⁴ / c
            ρ      = K / c²         from   c = sqrt(K/ρ)
            γ      = c / ξ          single-cell light-crossing rate (default)

        For the ``electron_drag`` anchor γ is *not* the light-crossing rate
        but is solved from the drag-as-mass-generator identity
            m_e c² = ℏ γ
        which fixes γ to the electron's Compton angular frequency.
        """
        self.xi = self._xi_from_anchor()
        self.K = HBAR * C_LIGHT / self.xi ** 4
        self.rho = self.K / C_LIGHT ** 2
        self.gamma = C_LIGHT / self.xi
        if self.anchor == "electron_drag":
            self.gamma = self.solve_with_drag_mass_anchor("electron")

    # -- drag-mass machinery ---------------------------------------------------
    def solve_with_drag_mass_anchor(self, particle: str = "electron") -> float:
        """Solve γ from the drag-as-mass-generator identity.

        Given (K, ρ, ξ) already fixed by the length anchor + (c, ℏ), the
        drag rate γ is determined by requiring the cone-bouncing rest-mass
        relation
            m c² = ℏ γ
        to reproduce the observed mass of ``particle``.

        Args:
            particle: Name in ``PARTICLE_MASS_KG`` (e.g. ``"electron"``).

        Returns:
            γ in 1/s such that ℏ γ = m c².
        """
        if particle not in PARTICLE_MASS_KG:
            raise ValueError(
                f"unknown particle {particle!r}; choose one of "
                f"{tuple(PARTICLE_MASS_KG)}"
            )
        m = PARTICLE_MASS_KG[particle]
        return m * C_LIGHT ** 2 / HBAR

    def drag_to_mass_ratio(self, particle: str = "electron") -> float:
        """Dimensionless ratio m c² / (ℏ γ) for the current anchor.

        If drag-as-mass holds and γ is the correct mass-generation primitive
        for ``particle``, this ratio is O(1) — exactly 1 in the canonical
        identity m c² = ℏ γ.

        Args:
            particle: Name in ``PARTICLE_MASS_KG``.
        """
        if particle not in PARTICLE_MASS_KG:
            raise ValueError(
                f"unknown particle {particle!r}; choose one of "
                f"{tuple(PARTICLE_MASS_KG)}"
            )
        m = PARTICLE_MASS_KG[particle]
        return m * C_LIGHT ** 2 / (HBAR * self.gamma)

    def derived_drag_mass_kg(self) -> float:
        """Rest mass implied by the current γ via m = ℏ γ / c²."""
        return HBAR * self.gamma / C_LIGHT ** 2

    def derived_drag_mass_MeV(self) -> float:
        """Drag-derived rest mass in MeV/c² for the current γ."""
        m_kg = self.derived_drag_mass_kg()
        return m_kg * C_LIGHT ** 2 / (1.0e6 * E_CHARGE)

    def cross_check_drag_alpha(self) -> Dict[str, float]:
        """Cross-check: with γ from drag, recompute α = 11/(48 π³)·exp(-3π/737).

        The B3 α formula is a pure-number identity that doesn't depend on
        γ directly, but here we verify that fixing γ via drag-mass leaves
        the α reconstruction (via the ℏc = K ξ⁴ identity) unchanged and
        consistent with the B3 closed form to the same residual as before.

        Returns:
            Dict with ``alpha_drag`` (the B3 closed-form value),
            ``alpha_from_primitives`` (via Coulomb identity), and
            ``residual_vs_B3`` (relative).
        """
        alpha_b3 = ALPHA_B3
        alpha_prim = self.derived_alpha()
        residual = (alpha_prim - alpha_b3) / alpha_b3
        return {
            "alpha_drag": alpha_b3,
            "alpha_from_primitives": alpha_prim,
            "residual_vs_B3": residual,
        }

    # -- derived cross-checks --------------------------------------------------
    def derived_alpha(self) -> float:
        """α from primitives via Coulomb identity α = e² / (4π ε_0 ℏ c).

        Reconstructed using K ξ⁴ = ℏ c → α = e² μ_0 c / (2 h) (CODATA).
        Returns the value implied by the substrate consistency.
        """
        # α = e² / (4π ε_0 ℏ c) — independent of anchor by construction
        # but expressed via primitives: ℏ c = K ξ⁴, so
        return E_CHARGE ** 2 / (4.0 * math.pi * EPS0 * self.K * self.xi ** 4)

    def derived_kink_mass_MeV(self) -> float:
        """Sine-Gordon kink mass m_kink = 8 ℏ / (c ξ) in MeV/c²."""
        m_kink_kg = 8.0 * HBAR / (C_LIGHT * self.xi)
        return m_kink_kg * C_LIGHT ** 2 / (1.0e6 * E_CHARGE)

    def derived_kink_over_electron(self) -> float:
        """Kink mass / electron mass — should equal 8 if ξ = λ_C(e)."""
        m_kink_kg = 8.0 * HBAR / (C_LIGHT * self.xi)
        return m_kink_kg / M_E

    def derived_kink_over_proton(self) -> float:
        return (8.0 * HBAR / (C_LIGHT * self.xi)) / M_P

    def derived_hierarchy(self) -> float:
        """Length-scale hierarchy ξ / ℓ_Planck — compare to exp(4π²-1).

        At ξ = λ_C(e) this is ≈ 2.4e22, compared to exp(4π²-1) ≈ 6.7e16.
        Ratio ≈ 3.5e5 — captures gravitational hierarchy magnitude only
        for the Planck-anchored case (ratio = 1 exactly).
        """
        return self.xi / ELL_PLANCK

    def derived_mp_over_me_via_kink(self) -> float:
        """If ξ = λ_C(e), kink = 8 m_e; ratio of m_p to that times 8."""
        m_kink_kg = 8.0 * HBAR / (C_LIGHT * self.xi)
        # m_p / m_kink × 8 reproduces m_p/m_e if ξ = λ_C(e)
        return 8.0 * M_P / m_kink_kg

    # -- cross-check report ---------------------------------------------------
    def compute_cross_checks(self) -> AnchorResult:
        derived = {
            "alpha": self.derived_alpha(),
            "m_kink_MeV": self.derived_kink_mass_MeV(),
            "kink_over_electron": self.derived_kink_over_electron(),
            "mp_over_me_implied": self.derived_mp_over_me_via_kink(),
            "xi_over_lPlanck": self.derived_hierarchy(),
            "K_Pa": self.K,
            "rho_kg_m3": self.rho,
            "gamma_Hz": self.gamma,
            "drag_mass_electron_MeV": self.derived_drag_mass_MeV(),
            "drag_to_mass_ratio_e": self.drag_to_mass_ratio("electron"),
            "drag_to_mass_ratio_p": self.drag_to_mass_ratio("proton"),
        }
        observed = {
            "alpha": ALPHA_OBS,
            "alpha_B3_formula": ALPHA_B3,
            "mp_over_me": MP_OVER_ME,
            "hierarchy_B3": HIERARCHY_B3,
            "kink_over_electron_target": 8.0,
        }
        residuals = {
            "alpha_vs_CODATA": (derived["alpha"] - ALPHA_OBS) / ALPHA_OBS,
            "alpha_vs_B3_formula": (derived["alpha"] - ALPHA_B3) / ALPHA_B3,
            "kink_8me_check": (
                derived["kink_over_electron"] - 8.0
            ) / 8.0,
            "mp_me_consistency": (
                derived["mp_over_me_implied"] - MP_OVER_ME
            ) / MP_OVER_ME,
            "hierarchy_B3_match": (
                derived["xi_over_lPlanck"] - HIERARCHY_B3
            ) / HIERARCHY_B3,
        }
        return AnchorResult(
            anchor=self.anchor,
            xi=self.xi,
            K=self.K,
            rho=self.rho,
            gamma=self.gamma,
            derived=derived,
            observed=observed,
            residuals=residuals,
        )

    # -- pretty printing ------------------------------------------------------
    def report_cross_checks(self) -> str:
        """Return a formatted multi-line table; also prints to stdout."""
        result = self.compute_cross_checks()
        lines: List[str] = []
        lines.append(f"=== Anchor: {self.anchor} ===")
        lines.append(f"  ξ      = {self.xi:.4e} m")
        lines.append(f"  K      = {self.K:.4e} Pa")
        lines.append(f"  ρ      = {self.rho:.4e} kg/m^3")
        lines.append(f"  γ      = {self.gamma:.4e} 1/s")
        lines.append("")
        lines.append("  Derived quantities:")
        for k, v in result.derived.items():
            lines.append(f"    {k:24s} = {v:.6e}")
        lines.append("")
        drag = self.cross_check_drag_alpha()
        lines.append("")
        lines.append("  Drag-mass cross-check:")
        lines.append(
            f"    m_e via drag (MeV)       = {self.derived_drag_mass_MeV():.6e}"
        )
        lines.append(
            f"    m c²/(ℏγ)  electron      = {self.drag_to_mass_ratio('electron'):.6e}"
        )
        lines.append(
            f"    m c²/(ℏγ)  proton        = {self.drag_to_mass_ratio('proton'):.6e}"
        )
        lines.append(
            f"    α (drag-formula)         = {drag['alpha_drag']:.6e}"
        )
        lines.append(
            f"    α residual vs B3 formula = {drag['residual_vs_B3']:+.3e}"
        )
        lines.append("")
        lines.append("  Residuals (relative):")
        for k, v in result.residuals.items():
            lines.append(f"    {k:24s} = {v:+.3e}")
        lines.append(f"  total |Δ|              = {result.total_residual():.3e}")
        text = "\n".join(lines)
        print(text)
        return text


# -----------------------------------------------------------------------------
# Multi-anchor comparison
# -----------------------------------------------------------------------------
def compare_anchors() -> Dict[str, AnchorResult]:
    """Run all anchors and return their AnchorResults keyed by anchor name."""
    out: Dict[str, AnchorResult] = {}
    for anchor in ANCHOR_CHOICES:
        out[anchor] = PrimitiveAnchoring(anchor).compute_cross_checks()
    return out


def best_anchor(results: Dict[str, AnchorResult] | None = None) -> str:
    """Return the anchor name with smallest total residual on the
    cross-checks where the anchor is *expected* to be self-consistent.

    The "kink_8me_check" and "mp_me_consistency" residuals are most
    sensitive to the anchor choice; we score by their combined |residual|.
    """
    if results is None:
        results = compare_anchors()
    scores: Dict[str, float] = {}
    for name, res in results.items():
        # electron_drag is a γ-anchor, not a length-anchor; its kink/mp checks
        # are degenerate with electron_compton — exclude from length-anchor
        # competition.
        if name == "electron_drag":
            continue
        scores[name] = abs(res.residuals["kink_8me_check"]) + abs(
            res.residuals["mp_me_consistency"]
        )
    return min(scores, key=scores.get)


def report_all_anchors() -> str:
    """Full multi-anchor report; identifies the best self-consistent anchor."""
    lines: List[str] = []
    results = compare_anchors()
    for anchor in ANCHOR_CHOICES:
        lines.append(PrimitiveAnchoring(anchor).report_cross_checks())
        lines.append("")
    winner = best_anchor(results)
    lines.append(f">>> Best self-consistent anchor: {winner}")
    text = "\n".join(lines)
    return text


if __name__ == "__main__":  # pragma: no cover
    report_all_anchors()
