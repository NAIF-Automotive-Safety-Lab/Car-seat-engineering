"""Substrate primitives — §18.2, §18.21, §18.22, §18.31, §18.32, §18.46.

Derives specific numerical values for K, ρ, ξ, ε_0, φ_max from the
constraint system imposed by observed physical constants.

Constraint system (in SI units throughout):
    (1)  c   = sqrt(K/ρ)                         [wave speed]
    (2)  hbar = K ξ⁴ / c                          [action quantum, §18.46.1]
    (3)  m_kink c² = 8 hbar c / ξ                 [sine-Gordon kink, natural units]
    (4)  α   = e² / (4π K ξ⁴) = e² / (4π hbar c) [redundant once (2) is used]
    (5)  Λ   = 8π G ρ_Λ / c²                      [determines ε_0 separately]
    (6)  σ_max = 0.5  →  φ_max                    [saturation cutoff, §18.39]

Key derivation (combining (2) and (3)):
    K = hbar c / ξ⁴     from (2)
    m_kink = 8 hbar / (c ξ)  from (3) with (2) substituted
           = 8 × (hbar / (m_e c)) / ξ × m_e   … if ξ = λ_C(electron)
    → when ξ = λ_C = hbar/(m_e c):  m_kink = 8 m_e  ≈ 4.1 MeV

    This is the CORRECT sine-Gordon kink mass in 1+1D:
        m_kink = 8 m_0    where m_0 = hbar/(c ξ) is the meson mass.

Note on the spec §18.21 formula "m_ν = 8 ρ ξ":
    The spec writes m_kink = 8 ρ ξ with ρ = K/c², claiming this gives 27 GeV.
    Dimensionally: [kg/m³][m] = kg/m²  — this is a surface mass density, not a mass.
    Numerically: 8 × 1.58e7 kg/m³ × 3.86e-13 m = 4.88e-5 kg/m²,
    and 4.88e-5 kg × c² / (1.602e-10 J/GeV) ≈ 2.7e22 GeV, NOT 27 GeV.
    The correct SI formula is m_kink = 8 hbar / (c ξ) = 8 m_e ≈ 4.1 MeV when ξ = λ_C.
    See KINK_MASS_NOTE in solve_from_constraints() for the full analysis.

Multi-scale finding:
    A single ξ can satisfy ALL three constraints (1)–(3) simultaneously.
    With ξ = λ_C(electron) = hbar/(m_e c):
        K   = hbar c / ξ⁴  ≈ 1.42e24 J/m³
        ρ   = K / c²        ≈ 1.58e7 kg/m³
        m_kink = 8 m_e ≈ 4.1 MeV  (kink is NOT a heavy carrier at this ξ)

    The multi-scale issue arises when trying to use a SINGLE ξ for:
        • hbar at Planck scale:  ξ_P = l_Planck ≈ 1.6e-35 m
        • kink = electron:       ξ_A = λ_C(e)  ≈ 3.9e-13 m
        • SM neutrino mass:      ξ_ν ~ hbar c / (m_ν c²) ≈ 2e-6 m
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Observed physical constants (CODATA 2022 / PDG 2024)
# ---------------------------------------------------------------------------

C_SI: Final[float] = 2.99792458e8         # m/s, exact
HBAR_SI: Final[float] = 1.054571817e-34   # J·s
M_E_SI: Final[float] = 9.1093837015e-31   # kg
E_CHARGE: Final[float] = 1.602176634e-19  # C, exact
ALPHA: Final[float] = 7.2973525693e-3     # ≈ 1/137.036
G_SI: Final[float] = 6.67430e-11          # m³ kg⁻¹ s⁻²
LAMBDA_OBS: Final[float] = 1.089e-52      # m⁻², observed cosmological constant

# Coulomb coupling  e² / (4π ε₀) [J·m]
COULOMB_COUPLING: Final[float] = (
    E_CHARGE**2 / (4 * math.pi * 8.8541878128e-12)
)

# Cosmological dark-energy density: ρ_Λ = Λ c² / (8π G)
RHO_LAMBDA_OBS: Final[float] = (
    LAMBDA_OBS * C_SI**2 / (8 * math.pi * G_SI)
)

# Planck length l_P = sqrt(hbar G / c³)
L_PLANCK: Final[float] = math.sqrt(HBAR_SI * G_SI / C_SI**3)

# Planck mass M_P = sqrt(hbar c / G)
M_PLANCK: Final[float] = math.sqrt(HBAR_SI * C_SI / G_SI)


# ---------------------------------------------------------------------------
# Data class for the substrate primitives
# ---------------------------------------------------------------------------

@dataclass
class SubstratePrimitives:
    """Container for the five substrate primitives K, ρ, ξ, ε_0, φ_max.

    All values are in SI units unless stated otherwise.

    Attributes:
        K:         Stiffness modulus  [J/m³]  = [Pa].
        rho:       Effective density  [kg/m³].
        xi:        Natural length scale ξ  [m].
        epsilon_0: Vacuum energy offset ε₀ [J/m³].
        phi_max:   Saturation field amplitude φ_max [m, natural strain unit].
        xi_label:  Human-readable description of which ξ this is.
        notes:     Free-form derivation notes.
    """

    K: float
    rho: float
    xi: float
    epsilon_0: float
    phi_max: float
    xi_label: str = ""
    notes: list[str] = field(default_factory=list)

    # ------------------------------------------------------------------
    # Verification methods — all check SI-consistent formulas
    # ------------------------------------------------------------------

    def verify_c(self) -> dict[str, float | bool]:
        """Verify c = sqrt(K/ρ).

        Returns:
            Dict with 'c_derived', 'c_observed', 'relative_error', 'consistent'.
        """
        c_derived = math.sqrt(self.K / self.rho)
        rel_err = abs(c_derived - C_SI) / C_SI
        return {
            "c_derived": c_derived,
            "c_observed": C_SI,
            "relative_error": rel_err,
            "consistent": rel_err < 1e-6,
        }

    def verify_hbar(self) -> dict[str, float | bool]:
        """Verify hbar = K ξ⁴ / c.

        This relation is the foundation of §18.46.1: the substrate's natural
        action quantum at length scale ξ.

        Returns:
            Dict with 'hbar_derived', 'hbar_observed', 'relative_error', 'consistent'.
        """
        hbar_derived = self.K * self.xi**4 / C_SI
        rel_err = abs(hbar_derived - HBAR_SI) / HBAR_SI
        return {
            "hbar_derived": hbar_derived,
            "hbar_observed": HBAR_SI,
            "relative_error": rel_err,
            "consistent": rel_err < 1e-6,
        }

    def verify_alpha(self) -> dict[str, float | bool]:
        """Verify α = e² / (4π K ξ⁴).

        Note: substituting K ξ⁴ = hbar c (from verify_hbar) gives
        α = e² / (4π hbar c), the standard textbook identity independent of ξ.
        So this constraint is REDUNDANT once the hbar constraint holds.

        Returns:
            Dict with 'alpha_derived', 'alpha_observed', 'relative_error',
            'consistent', 'is_redundant'.
        """
        Kxi4 = self.K * self.xi**4
        alpha_derived = COULOMB_COUPLING / Kxi4
        rel_err = abs(alpha_derived - ALPHA) / ALPHA
        # Verify redundancy: alpha using hbar c directly
        alpha_from_hbar = COULOMB_COUPLING / (HBAR_SI * C_SI)
        return {
            "alpha_derived": alpha_derived,
            "alpha_from_hbar_c": alpha_from_hbar,
            "alpha_observed": ALPHA,
            "relative_error": rel_err,
            "consistent": rel_err < 1e-6,
            "is_redundant": abs(alpha_from_hbar / ALPHA - 1) < 1e-6,
            "note": (
                "alpha = e²/(4π hbar c) — always true regardless of ξ. "
                "The alpha constraint gives NO new information about ξ."
            ),
        }

    def verify_electron_mass(self) -> dict[str, float | bool]:
        """Verify m_kink = 8 hbar / (c ξ)  (sine-Gordon kink mass).

        The CORRECT SI formula for the 1+1D sine-Gordon kink mass is:
            m_kink c² = 8 hbar c / ξ   [in Joules]
            m_kink    = 8 hbar / (c ξ)  [in kg]

        This equals 8 m_0 where m_0 = hbar/(cξ) is the small-oscillation
        (meson) mass scale of the substrate at length scale ξ.

        When ξ = λ_C(electron) = hbar/(m_e c), this gives m_kink = 8 m_e.

        Note on spec §18.46.2 formula "m_e c² = 8K/ξ":
            Dimensionally [J/m³]/[m] = J/m⁴ — NOT energy.
            The formula is inconsistent in SI. The correct version is
            m_e c² = 8 hbar c / ξ (which equals 8K ξ³/ρ = 8Kξ/c² × c² ...
            reduced to 8 hbar c/ξ via Kξ⁴ = hbar c).

        Returns:
            Dict with mass values, ratio to observed m_e, 'consistent'.
        """
        # Correct formula: m_kink = 8 hbar / (c ξ)
        m_kink = 8.0 * HBAR_SI / (C_SI * self.xi)  # kg
        m_kink_c2_J = m_kink * C_SI**2              # J
        m_kink_MeV = m_kink_c2_J / 1.602176634e-13  # MeV

        rel_err_vs_me = abs(m_kink - M_E_SI) / M_E_SI
        ratio_to_me = m_kink / M_E_SI

        return {
            "m_kink_kg": m_kink,
            "m_kink_MeV": m_kink_MeV,
            "me_observed_kg": M_E_SI,
            "ratio_kink_to_me": ratio_to_me,
            "relative_error_vs_me": rel_err_vs_me,
            "consistent_if_kink_is_me": rel_err_vs_me < 1e-3,
            "note": (
                f"Kink mass = {ratio_to_me:.4f} × m_e. "
                "When ξ = λ_C(electron), kink = 8 m_e ≈ 4.1 MeV. "
                "The kink is not the electron but 8× heavier."
            ),
        }

    # ------------------------------------------------------------------
    # Multi-scale analysis
    # ------------------------------------------------------------------

    def multi_scale_analysis(self) -> dict[str, object]:
        """Assess whether a single ξ can satisfy all constraints simultaneously.

        Computes the ξ implied by each physical scale and reports agreement.

        Key constraint algebra:
            From (2): K = hbar c / ξ⁴
            From (3): m_kink = 8 hbar / (c ξ)
            For m_kink = m_e: ξ = 8 hbar / (m_e c) = 8 λ_C(electron)
            For kink to be the electron: ξ = 8 λ_C
            For ξ = λ_C: m_kink = 8 m_e (kink is 8× heavier than electron)

        The α constraint is always satisfied (redundant).

        Returns:
            Dict with per-constraint ξ values, ratios, and honest verdict.
        """
        # ξ from Compton wavelength identification (§18.12, §18.21)
        xi_compton = HBAR_SI / (M_E_SI * C_SI)

        # ξ such that kink mass = m_e exactly (from 8ℏ/(cξ) = m_e c²)
        # → ξ = 8ℏ/(m_e c) = 8 λ_C
        xi_kink_equals_me = 8.0 * HBAR_SI / (M_E_SI * C_SI)

        # ξ_Planck = l_Planck  (where ℏ = Kξ⁴/c with Planck-scale stiffness)
        xi_planck = L_PLANCK

        # ξ for SM neutrino (m_ν ~ 0.1 eV): ξ_ν = ℏc/(m_ν c²)
        m_nu_J = 0.1 * 1.602176634e-19  # 0.1 eV in J
        xi_nu = HBAR_SI * C_SI / m_nu_J

        # Kink mass at ξ = λ_C (this object's ξ)
        m_kink_at_this_xi = 8.0 * HBAR_SI / (C_SI * self.xi)
        m_kink_MeV = m_kink_at_this_xi * C_SI**2 / 1.602176634e-13
        kink_to_me = m_kink_at_this_xi / M_E_SI

        # Is this ξ consistent with kink = electron?
        kink_equals_me = abs(kink_to_me - 1.0) < 0.01

        # Can ξ = λ_C satisfy ℏ and m_kink = m_e simultaneously?
        # No: at ξ = λ_C, m_kink = 8 m_e (factor 8 off).
        # At ξ = 8λ_C, m_kink = m_e exactly (but ξ ≠ λ_C).
        single_xi_ok = kink_equals_me

        if kink_equals_me:
            verdict = (
                "CONSISTENT: this ξ gives m_kink = m_e exactly. "
                "(ξ = 8 λ_C(electron))"
            )
        else:
            verdict = (
                f"MULTI-SCALE: at ξ = λ_C(electron), kink mass = {kink_to_me:.3f} × m_e "
                f"= {m_kink_MeV:.3f} MeV. "
                "To get kink = m_e exactly, need ξ = 8 λ_C. "
                "Different physics lives at different ξ scales."
            )

        return {
            "xi_compton_m": xi_compton,
            "xi_for_kink_equals_me_m": xi_kink_equals_me,
            "xi_planck_m": xi_planck,
            "xi_neutrino_m": xi_nu,
            "m_kink_at_this_xi_kg": m_kink_at_this_xi,
            "m_kink_at_this_xi_MeV": m_kink_MeV,
            "kink_to_me_ratio": kink_to_me,
            "single_xi_satisfies_kink_equals_me": single_xi_ok,
            "verdict": verdict,
        }


# ---------------------------------------------------------------------------
# Factory: solve constraint system from observables
# ---------------------------------------------------------------------------

#: Note on the spec §18.21 "m_ν = 8ρξ ≈ 27 GeV" claim.
#:
#: The spec §18.21 writes: "m_ν = 8 × 1.58×10⁷ × 3.86×10⁻¹³ kg ≈ 4.88×10⁻⁵ kg ≈ 27 GeV/c²"
#: The dimensional check: 8 × [kg/m³] × [m] = kg/m² (surface mass density, NOT mass).
#: And numerically: 4.88e-5 kg × c² / (1.602e-10 J/GeV) ≈ 2.7×10²² GeV, NOT 27 GeV.
#:
#: The CORRECT 1+1D sine-Gordon kink mass in SI units is:
#:     m_kink = 8 ℏ / (c ξ)  [kg]
#: At ξ = λ_C(electron) = ℏ/(m_e c):
#:     m_kink = 8 ℏ / (c × ℏ/(m_e c)) = 8 m_e ≈ 4.1 MeV
#:
#: The spec's "27 GeV" figure does not follow from its stated formula.
#: The kink at ξ = λ_C is ~4.1 MeV, not 27 GeV.
#: The multi-scale issue remains: this kink (4.1 MeV) is heavier than m_e (511 keV)
#: by factor 8, not equal to it. The kink is NOT the electron at this ξ.
KINK_MASS_NOTE: str = (
    "Spec §18.21 formula m_kink = 8ρξ is dimensionally inconsistent (gives kg/m², not kg). "
    "Correct formula: m_kink = 8ℏ/(cξ). At ξ=λ_C(e): m_kink = 8m_e ≈ 4.1 MeV. "
    "The kink is ~8× heavier than the electron at ξ=λ_C, NOT 27 GeV."
)


def solve_from_constraints() -> tuple[SubstratePrimitives, SubstratePrimitives, dict]:
    """Determine K, ρ, ξ from observed physical constants.

    Two solutions are constructed based on two interpretations of ξ:

    Solution A — ξ = λ_C(electron) (§18.21 primary):
        ξ_A  = ℏ / (m_e c)        [electron Compton wavelength]
        K_A  = ℏ c / ξ_A⁴         [from ℏ = K ξ⁴/c]
        ρ_A  = K_A / c²            [from c = sqrt(K/ρ)]
        m_kink = 8ℏ/(cξ_A) = 8m_e ≈ 4.1 MeV  [kink is 8× electron mass]

    Solution B — ξ chosen so kink mass = m_e exactly:
        m_kink = 8ℏ/(cξ) = m_e  →  ξ_B = 8ℏ/(m_e c) = 8 λ_C
        K_B  = ℏ c / ξ_B⁴
        ρ_B  = K_B / c²
        m_kink = m_e by construction; ξ_B = 8 λ_C ≈ 3.09e-12 m

    ε_0 determination (for observed Λ, §18.46.2):
        ρ_Λ = Λ c² / (8π G)   [dark energy density]
        ε_0 = ρ_Λ c²           [vacuum energy offset in J/m³]

    φ_max determination (§18.39, σ_max = 0.5):
        φ_max = ξ  (natural choice: saturation at kink amplitude scale)

    Returns:
        Tuple of (solution_A, solution_B, analysis_dict).
    """
    # ------------------------------------------------------------------ #
    # SOLUTION A: ξ = λ_C(electron) — the §18.21 primary route          #
    # ------------------------------------------------------------------ #
    xi_A = HBAR_SI / (M_E_SI * C_SI)         # 3.862e-13 m
    K_A = HBAR_SI * C_SI / xi_A**4           # J/m³
    rho_A = K_A / C_SI**2                    # kg/m³
    m_kink_A_kg = 8.0 * HBAR_SI / (C_SI * xi_A)     # = 8 m_e
    m_kink_A_MeV = m_kink_A_kg * C_SI**2 / 1.602176634e-13

    # Vacuum energy offset from observed Λ
    epsilon_0 = RHO_LAMBDA_OBS * C_SI**2     # J/m³

    # φ_max = ξ (natural saturation scale: strain amplitude = length scale of kink)
    phi_max_A = xi_A

    sol_A = SubstratePrimitives(
        K=K_A,
        rho=rho_A,
        xi=xi_A,
        epsilon_0=epsilon_0,
        phi_max=phi_max_A,
        xi_label="xi = lambda_C(electron) = hbar/(m_e c) = 3.862e-13 m",
        notes=[
            "ξ = λ_C(electron): primary identification of §18.21.",
            "ℏ = Kξ⁴/c is exactly satisfied by construction.",
            "α = e²/(4πKξ⁴) = e²/(4πℏc) = α_fine — redundant, no ξ info.",
            f"Kink mass = 8ℏ/(cξ) = 8m_e ≈ {m_kink_A_MeV:.2f} MeV (NOT the electron).",
            "ρ_A ≈ 1.58e7 kg/m³  (white-dwarf density range).",
            "K_A ≈ 1.42e24 J/m³  (very stiff substrate).",
            "Kink is 8× heavier than electron — a distinct heavier excitation.",
        ],
    )

    # ------------------------------------------------------------------ #
    # SOLUTION B: ξ = 8λ_C so kink mass = m_e exactly                  #
    # ------------------------------------------------------------------ #
    # 8ℏ/(cξ) = m_e  →  ξ_B = 8ℏ/(m_e c) = 8 λ_C
    xi_B = 8.0 * HBAR_SI / (M_E_SI * C_SI)  # = 8 xi_A
    K_B = HBAR_SI * C_SI / xi_B**4
    rho_B = K_B / C_SI**2
    phi_max_B = xi_B

    sol_B = SubstratePrimitives(
        K=K_B,
        rho=rho_B,
        xi=xi_B,
        epsilon_0=epsilon_0,
        phi_max=phi_max_B,
        xi_label=f"xi = 8*lambda_C = 8*hbar/(m_e c) = {xi_B:.3e} m",
        notes=[
            "ξ = 8λ_C: chosen so kink mass = m_e exactly.",
            "Kink mass = 8ℏ/(cξ) = 8ℏ/(c×8λ_C) = m_e by construction.",
            "ℏ = Kξ⁴/c also exactly satisfied.",
            "α constraint still redundant (independent of ξ).",
            "But ξ ≠ λ_C: the kink and Compton scales are distinct.",
            "ρ_B ≈ 3.8e18 kg/m³  (neutron-star to black-hole density range).",
            "K_B much smaller; substrate softer at larger ξ.",
        ],
    )

    # ------------------------------------------------------------------ #
    # Constraint analysis                                                 #
    # ------------------------------------------------------------------ #
    # Verify α redundancy: α = e²/(4πℏc) exactly
    alpha_check = COULOMB_COUPLING / (HBAR_SI * C_SI)

    analysis: dict[str, object] = {
        "constraints": {
            "c = sqrt(K/rho)": "INDEPENDENT — fixes ρ once K is known",
            "hbar = K xi^4 / c": "INDEPENDENT — fixes K once ξ is known",
            "m_kink = 8 hbar/(c xi)": "INDEPENDENT — fixes ξ once m_kink target set",
            "alpha = e^2/(4pi K xi^4)": (
                "REDUNDANT: substituting K xi^4 = hbar c gives "
                f"alpha = e^2/(4pi hbar c) = {alpha_check:.6e} "
                f"= alpha_obs ({ALPHA:.6e}). No information about ξ."
            ),
            "Lambda = 8pi G rho_Lambda / c^2": "INDEPENDENT — determines ε_0 separately",
            "sigma_max = 0.5 -> phi_max": "SEMI-INDEPENDENT — phi_max = ξ (natural choice)",
        },
        "independent_constraints": 4,    # c, ℏ, m_kink target, Λ
        "redundant_constraints": 1,       # α
        "solution_A": {
            "xi_m": xi_A,
            "K_Jpm3": K_A,
            "rho_kgpm3": rho_A,
            "kink_mass_kg": m_kink_A_kg,
            "kink_mass_MeV": m_kink_A_MeV,
            "kink_to_me_ratio": m_kink_A_kg / M_E_SI,
            "interpretation": "kink = 8 m_e ≈ 4.1 MeV (distinct from electron)",
        },
        "solution_B": {
            "xi_m": xi_B,
            "K_Jpm3": K_B,
            "rho_kgpm3": rho_B,
            "kink_mass_kg": M_E_SI,
            "kink_mass_MeV": 0.51099895,
            "kink_to_me_ratio": 1.0,
            "interpretation": "kink = m_e exactly by construction (ξ = 8λ_C)",
        },
        "xi_ratio_B_over_A": xi_B / xi_A,  # = 8.0 exactly
        "epsilon_0_Jpm3": epsilon_0,
        "rho_Lambda_kgpm3": RHO_LAMBDA_OBS,
        "alpha_redundancy_check": alpha_check,
        "alpha_observed": ALPHA,
        "kink_mass_note": KINK_MASS_NOTE,
    }

    return sol_A, sol_B, analysis
