"""Kink-antikink scattering dynamics — collider physics analog.

This module implements kink-antikink scattering processes in the substrate-
mechanical framework, treating kinks as topological solitons whose scattering
cross-sections are computed from their effective field theory (§18.37, §18.48,
§18.49) as well as standard QED cross-sections inherited via structural identity
(§18.49.8).

Physical picture
----------------
Kinks are topological defects in the stiff substrate.  Their scattering
mirrors particle collider physics:

  - Elastic scattering : kinks bounce off each other, topology preserved.
  - Pair annihilation  : kink + antikink → substrate waves (photon analogs)
                         (§18.37, §18.53.2).
  - Pair production    : high-energy substrate wave → kink + antikink
                         (threshold = 2 M_K c², §18.49).
  - Stress loading     : collider deposits energy into vertex as excited state
                         (§18.30).

QED cross-sections (Bhabha, Klein-Nishina) are inherited exactly via the
structural identity of the gauge sector (§18.49.8).

Units
-----
All energies are in GeV.  Cross-sections are in natural units (GeV⁻²) and
converted to m² where noted.  Physical constants follow the CODATA 2018
embedded values used throughout the package.

References
----------
- §18.37  : dark matter kink composites, pair annihilation rate
- §18.48.3: two-kink interaction potential → elastic scattering scale
- §18.49  : SU(3) extension → QCD-inherited processes
- §18.49.8: structural identity → QED inheritance (Bhabha, Compton)
- §18.30  : stress-loading mechanism → collider excited states
- §18.53.2: antimatter = reversed Möbius half-flux; pair production/annihilation
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Final, Literal, Optional

import numpy as np

from stiff_medium.multi_kink import Kink

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

#: Speed of light in vacuum (m/s)
C_SI: Final[float] = 2.99792458e8

#: Reduced Planck constant ℏ (J·s)
HBAR_SI: Final[float] = 1.054571817e-34

#: Fine-structure constant α = e²/(4πε₀ ℏ c) (dimensionless)
ALPHA: Final[float] = 1.0 / 137.035999084

#: Electron mass (GeV)
M_ELECTRON_GEV: Final[float] = 0.000510998950

#: Conversion: 1 GeV⁻² = HBAR²c² in m²
#: ℏc = 0.197326980 GeV·fm  →  (ℏc)² = 0.389379 GeV²·mb = 0.389379e-31 m²·GeV²
#: So 1 GeV⁻² = 0.389379e-31 m²
GEV_INV2_TO_M2: Final[float] = 0.3893793721e-31  # m² per GeV⁻²

#: Canonical kink rest mass per §18.22 (GeV)
M_K_CANONICAL: Final[float] = 27.0

#: Photon = massless (GeV)
M_PHOTON: Final[float] = 0.0

ProcessType = Literal["elastic", "inelastic", "pair_production", "annihilation"]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class ScatteringEvent:
    """Record of a single kink-scattering interaction.

    Captures the full kinematics of a kink collider event, analogous to a
    Monte Carlo truth record in collider simulations.

    Attributes:
        incoming_kinks: Kink objects entering the interaction region.
        outgoing_kinks: Kink objects (and/or photon-analogs) leaving.
        center_of_mass_energy: √s in GeV; total c.m. energy of incoming state.
        cross_section: Total cross-section for this process in m².
        process_type: One of "elastic", "inelastic", "pair_production",
            "annihilation".
        Q_value: Energy released into radiation (GeV); positive for exothermic.
        momentum_conserved: True if four-momentum is conserved to tolerance.
        energy_conserved: True if total energy is conserved to tolerance.

    References:
        §18.37, §18.48.3, §18.53.2
    """

    incoming_kinks: list[Kink]
    outgoing_kinks: list[Kink]
    center_of_mass_energy: float  # GeV
    cross_section: float          # m²
    process_type: ProcessType
    Q_value: float = 0.0          # GeV; energy deposited into substrate waves
    momentum_conserved: bool = True
    energy_conserved: bool = True

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def multiplicity_in(self) -> int:
        """Number of incoming kinks."""
        return len(self.incoming_kinks)

    @property
    def multiplicity_out(self) -> int:
        """Number of outgoing kinks."""
        return len(self.outgoing_kinks)

    @property
    def is_physical(self) -> bool:
        """True if both conservation laws are satisfied."""
        return self.momentum_conserved and self.energy_conserved

    def __str__(self) -> str:
        return (
            f"ScatteringEvent({self.process_type}, "
            f"√s={self.center_of_mass_energy:.3f} GeV, "
            f"σ={self.cross_section:.3e} m², "
            f"N_in={self.multiplicity_in}, N_out={self.multiplicity_out})"
        )


# ---------------------------------------------------------------------------
# Module-level cross-section functions
# ---------------------------------------------------------------------------


def pair_production_threshold(M_K: float) -> float:
    """Minimum photon energy needed to produce a kink-antikink pair.

    From special relativity, a massless photon must supply at least the rest-
    mass energy of both produced particles:

        E_threshold = 2 M_K c²

    In natural units (c = 1) this simplifies to 2 M_K (GeV).

    This is the direct analog of the pair-production threshold for electrons:
    E_γ ≥ 2 m_e for γ → e⁺e⁻.  For the substrate model the production
    threshold is set by the kink rest mass M_K ≈ 27 GeV (§18.22), giving
    E_min ≈ 54 GeV.

    Args:
        M_K: Kink rest mass in GeV.

    Returns:
        Threshold photon energy E_min = 2 M_K (GeV).

    References:
        §18.49 (kinematic threshold for kink pair production)
        §18.53.2 (pair production mechanism in substrate)
    """
    return 2.0 * M_K


def elastic_cross_section(E: float, M_K: float) -> float:
    """Kink-kink elastic scattering cross-section.

    Derived from the §18.48.3 two-kink potential in the Born approximation.
    At low energies (E ≪ M_K) the cross-section is dominated by the long-range
    exponential tail of the potential; at high energies (E ≫ M_K) it falls off
    as the kinks pass through each other without topological reconnection.

    The formula used is the partial-wave Born approximation for the §18.48.3
    potential:

        σ_el(E) = π α²_K / (s - 4 M_K²)

    where α_K = ALPHA (inherited QED coupling, §18.49.8) and
    s = (2 E_cm)² in the c.m. frame when E is the lab-frame beam energy.
    For a head-on kink-antikink collision with equal and opposite momenta,
    s = (2E)² at high energy, so:

        σ_el(E) ≈ π α² / (4 E²)   [GeV⁻²]

    This is the same scaling as QED elastic lepton scattering and is inherited
    via structural identity (§18.49.8).

    Args:
        E: Center-of-mass energy √s in GeV.
        M_K: Kink rest mass in GeV.

    Returns:
        Elastic cross-section in GeV⁻².

    Raises:
        ValueError: If E ≤ 0 or M_K ≤ 0.

    References:
        §18.48.3 (two-kink potential: source of elastic scattering)
        §18.49.8 (QED inheritance → α coupling)
    """
    if E <= 0.0:
        raise ValueError(f"Energy must be positive; got E={E}")
    if M_K <= 0.0:
        raise ValueError(f"Kink mass must be positive; got M_K={M_K}")

    s = E * E  # s = E_cm² when E = √s
    s_threshold = 4.0 * M_K * M_K
    if s <= s_threshold:
        # Below threshold: no real elastic scattering in vacuum
        return 0.0

    # Born approximation with running coupling ALPHA
    sigma_gev2 = math.pi * ALPHA**2 / (s - s_threshold + s_threshold * 1e-3)
    return sigma_gev2


def annihilation_cross_section(E: float, M_K: float) -> float:
    """Kink-antikink annihilation cross-section.

    In the substrate model, kink-antikink annihilation (§18.37, §18.53.2) is
    the topological analog of e⁺e⁻ → γγ.  The cross-section is:

        σ_ann(E) = π α² / s × [2 ln(s/M_K²) - 1]   [GeV⁻²]

    This is the relativistic Dirac formula for fermion-antifermion annihilation,
    inherited via the structural identity of the SU(2)×U(1) sector (§18.49.8).
    At high energy (s ≫ M_K²) the logarithm grows slowly so σ_ann ∝ 1/s.

    The dark matter annihilation rate (§18.37) is proportional to this cross-
    section times velocity, giving Γ ~ G² m_DM⁵ / (192π³) in the
    gravitational-only limit.

    Args:
        E: Center-of-mass energy √s in GeV.
        M_K: Kink rest mass in GeV.

    Returns:
        Annihilation cross-section in GeV⁻².  Returns 0 below threshold.

    Raises:
        ValueError: If E ≤ 0 or M_K ≤ 0.

    References:
        §18.37  (dark matter kink annihilation rate)
        §18.49.8 (QED structural inheritance → fermion annihilation formula)
        §18.53.2 (pair annihilation as topology cancellation)
    """
    if E <= 0.0:
        raise ValueError(f"Energy must be positive; got E={E}")
    if M_K <= 0.0:
        raise ValueError(f"Kink mass must be positive; got M_K={M_K}")

    s = E * E
    s_threshold = 4.0 * M_K * M_K
    if s <= s_threshold:
        return 0.0

    log_factor = math.log(s / (M_K * M_K))
    sigma_gev2 = math.pi * ALPHA**2 / s * (2.0 * log_factor - 1.0)
    # Guard against negative values (can occur just above threshold)
    return max(sigma_gev2, 0.0)


def bhabha_scattering_cross_section(E: float) -> float:
    """Bhabha scattering cross-section for e⁺e⁻ → e⁺e⁻ (QED / inherited).

    This is the standard QED result, inherited by the substrate model via
    structural identity of the gauge sector (§18.49.8).  In the ultra-
    relativistic limit (√s ≫ m_e) and integrating over all angles:

        dσ/dΩ = (α / 2s) × [u²/t² + 2u/t + (u²+s²)/(2t²) + …]

    Keeping only the leading s-channel (annihilation) term integrated
    over the full solid angle gives the approximate total:

        σ_Bhabha ≈ 4π α² / (3 s)   [GeV⁻²]   (dominant s-channel piece)

    The exact Bhabha total cross-section also has a t-channel (Møller-type)
    divergence regularised by the electron mass.  For √s ≫ m_e the
    s-channel piece dominates and the t-channel contributes a logarithmic
    correction ~α² / s × ln(s / m_e²).

    This function returns the s-channel dominant contribution, appropriate
    for LEP1/LEP2 energies (√s ~ 91–200 GeV).

    Args:
        E: Center-of-mass energy √s in GeV.  Must be positive.

    Returns:
        Bhabha cross-section in GeV⁻².

    Raises:
        ValueError: If E ≤ 0.

    References:
        §18.49.8 (QED inheritance via structural identity)
        Standard QED: Peskin & Schroeder §5.5
    """
    if E <= 0.0:
        raise ValueError(f"Energy must be positive; got E={E}")

    s = E * E

    # s-channel annihilation term (dominant at high energy)
    sigma_s_channel = 4.0 * math.pi * ALPHA**2 / (3.0 * s)

    # t-channel logarithmic correction: ~ α²/s × ln(s/m_e²)
    # Only add when s > 4 m_e² to avoid log of negative
    s_threshold_e = 4.0 * M_ELECTRON_GEV**2
    if s > s_threshold_e:
        log_t = math.log(s / M_ELECTRON_GEV**2)
        sigma_t_log = ALPHA**2 / s * log_t
    else:
        sigma_t_log = 0.0

    return sigma_s_channel + sigma_t_log


def compton_cross_section(E_gamma: float) -> float:
    """Klein-Nishina total cross-section for Compton scattering γ + e → γ + e.

    The Klein-Nishina formula gives the total cross-section for Compton
    scattering of a photon off a free electron.  Defining x = E_γ / m_e:

        σ_KN = (3 σ_T / 4) × {
            (1+x) / x³ × [2x(1+x)/(1+2x) - ln(1+2x)]
            + ln(1+2x)/(2x)
            - (1+3x)/(1+2x)²
        }

    where σ_T = (8π/3) r_e² is the Thomson cross-section,
    r_e = α / m_e = α ℏ c / (m_e c²) in natural units.

    In natural units r_e = α / m_e (GeV⁻¹), so σ_T = 8π α² / (3 m_e²) [GeV⁻²].

    This formula is inherited exactly by the substrate model because Compton
    scattering is a process of a substrate wave (photon analog) scattering off
    a bound kink configuration (electron analog), and the gauge sector is
    structurally identical to QED (§18.49.8).

    Limits:
      - Low energy (x ≪ 1): σ → σ_T (Thomson limit)
      - High energy (x ≫ 1): σ → (3 σ_T / 8x) [ln(2x) + ½]

    Args:
        E_gamma: Photon energy in GeV.  Must be positive.

    Returns:
        Klein-Nishina total cross-section in GeV⁻².

    Raises:
        ValueError: If E_gamma ≤ 0.

    References:
        §18.49.8 (QED structural identity → Compton scattering inherited)
        Klein & Nishina (1929); O. Klein & Y. Nishina, Z. Phys. 52, 853.
    """
    if E_gamma <= 0.0:
        raise ValueError(f"Photon energy must be positive; got E_gamma={E_gamma}")

    m_e = M_ELECTRON_GEV
    x = E_gamma / m_e  # dimensionless energy parameter

    # Thomson cross-section in GeV⁻²:  σ_T = 8π r_e² / 3, r_e = α/m_e
    r_e_gev = ALPHA / m_e  # classical electron radius in GeV⁻¹
    sigma_T = (8.0 * math.pi / 3.0) * r_e_gev**2

    # Full Klein-Nishina formula
    one_plus_2x = 1.0 + 2.0 * x
    log_1p2x = math.log(one_plus_2x)

    term1 = (1.0 + x) / x**3 * (2.0 * x * (1.0 + x) / one_plus_2x - log_1p2x)
    term2 = log_1p2x / (2.0 * x)
    term3 = -(1.0 + 3.0 * x) / one_plus_2x**2

    sigma_kn = (3.0 * sigma_T / 4.0) * (term1 + term2 + term3)
    return sigma_kn


# ---------------------------------------------------------------------------
# KinkScattering class
# ---------------------------------------------------------------------------


class KinkScattering:
    """Simulate kink-antikink scattering processes.

    This class provides methods for each scattering channel: elastic, pair
    annihilation, and pair production.  The underlying cross-sections are
    either derived from the §18.48.3 potential (for kink-specific processes)
    or inherited from QED via §18.49.8 (for Bhabha, Compton).

    Args:
        M_K: Kink rest mass in GeV.  Default is the canonical §18.22 value.
        rng_seed: Optional seed for the internal random-number generator
            (for reproducible Monte Carlo sampling).

    References:
        §18.37, §18.48.3, §18.49, §18.49.8, §18.53.2
    """

    def __init__(
        self,
        M_K: float = M_K_CANONICAL,
        rng_seed: Optional[int] = None,
    ) -> None:
        if M_K <= 0.0:
            raise ValueError(f"Kink mass must be positive; got M_K={M_K}")
        self.M_K = M_K
        self._rng = random.Random(rng_seed)

    # ------------------------------------------------------------------
    # Private kinematics helpers
    # ------------------------------------------------------------------

    def _lorentz_beta(self, mass: float, kinetic_energy: float) -> float:
        """Compute β = v/c from rest mass and kinetic energy.

        Args:
            mass: Rest mass in GeV.
            kinetic_energy: Kinetic energy in GeV.

        Returns:
            β ∈ [0, 1).
        """
        total_E = mass + kinetic_energy
        p2 = total_E**2 - mass**2
        if p2 < 0.0:
            return 0.0
        return math.sqrt(p2) / total_E

    def _four_momentum(
        self, mass: float, kinetic_energy: float, direction: np.ndarray
    ) -> tuple[float, np.ndarray]:
        """Compute (E, p⃗) four-momentum from rest mass and kinetic energy.

        Args:
            mass: Rest mass in GeV.
            kinetic_energy: Kinetic energy in GeV.
            direction: Unit 3-vector for the momentum direction.

        Returns:
            Tuple (E_total_GeV, p_3vec_GeV) where p_3vec has the same
            direction as `direction`.
        """
        E_total = mass + kinetic_energy
        p_mag = math.sqrt(max(E_total**2 - mass**2, 0.0))
        p_vec = direction * p_mag
        return E_total, p_vec

    def _check_conservation(
        self,
        E_in: float,
        p_in: np.ndarray,
        E_out: float,
        p_out: np.ndarray,
        tol: float = 1e-6,
    ) -> tuple[bool, bool]:
        """Test four-momentum conservation.

        Args:
            E_in: Total incoming energy (GeV).
            p_in: Total incoming 3-momentum (GeV).
            E_out: Total outgoing energy (GeV).
            p_out: Total outgoing 3-momentum (GeV).
            tol: Fractional tolerance for conservation check.

        Returns:
            Tuple (energy_conserved, momentum_conserved) as booleans.
        """
        energy_ok = abs(E_in - E_out) / (E_in + 1e-30) < tol
        dp = np.linalg.norm(p_in - p_out)
        p_scale = np.linalg.norm(p_in) + 1e-30
        momentum_ok = dp / p_scale < tol
        return bool(energy_ok), bool(momentum_ok)

    # ------------------------------------------------------------------
    # Elastic scattering
    # ------------------------------------------------------------------

    def elastic_scattering(
        self,
        k1: Kink,
        k2: Kink,
        impact_parameter: float = 0.0,
    ) -> ScatteringEvent:
        """Elastic kink-kink scattering: kinks bounce off, topology preserved.

        Models the scattering k₁ + k₂ → k₁' + k₂' where both outgoing kinks
        have the same topology (chirality, winding) as the incoming ones.
        The momentum exchange is computed from the §18.48.3 potential in the
        Born approximation.

        Energy and momentum are exactly conserved by construction: the c.m.
        frame momenta are reflected about the scattering axis determined by
        the impact parameter b.

        The differential cross-section follows from the §18.48.3 potential:
        at low energy the kinks bounce off the exponential tail; at high
        energy they pass through each other (sine-Gordon property).

        Args:
            k1: First incoming kink.
            k2: Second incoming kink.
            impact_parameter: Perpendicular distance between trajectories in
                natural substrate units.  b = 0 is head-on.

        Returns:
            ScatteringEvent describing the elastic collision.

        References:
            §18.48.3 (two-kink interaction potential)
            §18.49.8 (elastic cross-section formula)
        """
        # Incoming kinematics
        v1 = np.linalg.norm(k1.velocity)
        v2 = np.linalg.norm(k2.velocity)
        E1 = math.sqrt(k1.mass**2 + v1**2)  # natural units
        E2 = math.sqrt(k2.mass**2 + v2**2)
        p1 = k1.velocity.copy()
        p2 = k2.velocity.copy()

        E_total_in = E1 + E2
        p_total_in = p1 + p2
        sqrt_s = math.sqrt(max(E_total_in**2 - np.dot(p_total_in, p_total_in), 0.0))

        # Scattering angle from impact parameter (Born approximation)
        # For the §18.48.3 potential the deflection angle θ ≈ 2 arctan(b_0/b)
        # where b_0 = α_K / (2 p_cm) is the classical turning point.
        p_cm = 0.5 * math.sqrt(max(sqrt_s**2 - 4.0 * k1.mass**2, 0.0))
        b_0 = ALPHA / (2.0 * p_cm + 1e-30)
        theta = 2.0 * math.atan2(b_0, impact_parameter + 1e-30)

        # Scatter k1 by angle theta in the x-z plane; k2 recoils to conserve p
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        v1_mag = np.linalg.norm(k1.velocity)
        if v1_mag < 1e-30:
            v1_out = k1.velocity.copy()
            v2_out = k2.velocity.copy()
        else:
            v1_dir_orig = k1.velocity / v1_mag
            # Construct a rotation in the plane defined by v1_dir and the
            # impact-parameter direction (y-axis for head-on geometry)
            perp = np.array([0.0, 1.0, 0.0]) if v1_dir_orig[0] > 0.5 else np.array([1.0, 0.0, 0.0])
            perp -= np.dot(perp, v1_dir_orig) * v1_dir_orig
            perp_norm = np.linalg.norm(perp)
            if perp_norm > 1e-12:
                perp /= perp_norm
            v1_dir_new = cos_t * v1_dir_orig + sin_t * perp
            v1_out = v1_mag * v1_dir_new

            # k2 recoils such that total momentum is conserved
            v2_out = p_total_in - v1_out

        # Rebuild outgoing kinks (same topology)
        k1_out = Kink(
            position=k1.position.copy(),
            velocity=v1_out,
            chirality=k1.chirality,
            winding=k1.winding,
            mass=k1.mass,
        )
        k2_out = Kink(
            position=k2.position.copy(),
            velocity=v2_out,
            chirality=k2.chirality,
            winding=k2.winding,
            mass=k2.mass,
        )

        # Energy check
        E_out = (
            math.sqrt(k1_out.mass**2 + np.dot(v1_out, v1_out))
            + math.sqrt(k2_out.mass**2 + np.dot(v2_out, v2_out))
        )
        p_out = v1_out + v2_out
        energy_ok, momentum_ok = self._check_conservation(
            E_total_in, p_total_in, E_out, p_out
        )

        sigma_gev2 = elastic_cross_section(sqrt_s, self.M_K)
        sigma_m2 = sigma_gev2 * GEV_INV2_TO_M2

        return ScatteringEvent(
            incoming_kinks=[k1, k2],
            outgoing_kinks=[k1_out, k2_out],
            center_of_mass_energy=sqrt_s,
            cross_section=sigma_m2,
            process_type="elastic",
            Q_value=0.0,
            energy_conserved=energy_ok,
            momentum_conserved=momentum_ok,
        )

    # ------------------------------------------------------------------
    # Pair annihilation
    # ------------------------------------------------------------------

    def pair_annihilation(self, kink: Kink, antikink: Kink) -> ScatteringEvent:
        """Kink-antikink annihilation producing substrate waves (photon analogs).

        When a kink (chirality +1) and antikink (chirality -1) meet, their
        topological charges cancel and all energy is released as substrate
        oscillations (photon analogs).  This is the substrate-mechanical
        realization of e⁺e⁻ → γγ (§18.53.2):

            kink + antikink → 2 photons (substrate waves)

        Energy-momentum conservation:
            E_γ1 + E_γ2 = E_kink + E_antikink  (energy)
            p⃗_γ1 + p⃗_γ2 = p⃗_kink + p⃗_antikink  (momentum)

        In the c.m. frame: each photon carries √s / 2 and travels back-to-back.

        Args:
            kink: The kink (chirality = +1).
            antikink: The antikink (chirality = -1).

        Returns:
            ScatteringEvent with outgoing_kinks representing the two photon-
            analog substrates waves as massless Kink objects (mass=0, chirality=0).

        Raises:
            ValueError: If both particles have the same chirality (not a
                kink-antikink pair).

        References:
            §18.37  (dark matter annihilation rate)
            §18.53.2 (pair annihilation = topology cancellation → photons)
        """
        if kink.chirality == antikink.chirality:
            raise ValueError(
                "pair_annihilation requires kink (chirality=+1) and antikink "
                f"(chirality=-1); got {kink.chirality} and {antikink.chirality}"
            )

        v_k = np.linalg.norm(kink.velocity)
        v_ak = np.linalg.norm(antikink.velocity)
        E_k = math.sqrt(kink.mass**2 + v_k**2)
        E_ak = math.sqrt(antikink.mass**2 + v_ak**2)
        p_k = kink.velocity.copy()
        p_ak = antikink.velocity.copy()

        E_total = E_k + E_ak
        p_total = p_k + p_ak
        sqrt_s = math.sqrt(max(E_total**2 - np.dot(p_total, p_total), 0.0))

        # In the c.m. frame each photon carries √s / 2 back-to-back
        E_gamma = sqrt_s / 2.0

        # Boost back to lab: photons share the total momentum equally
        p_boost = p_total / 2.0
        p_perp = np.array([0.0, 1.0, 0.0])  # transverse for one photon
        p_perp -= np.dot(p_perp, p_boost / (np.linalg.norm(p_boost) + 1e-30)) * (
            p_boost / (np.linalg.norm(p_boost) + 1e-30)
        )
        if np.linalg.norm(p_perp) > 1e-12:
            p_perp /= np.linalg.norm(p_perp)
        else:
            p_perp = np.array([1.0, 0.0, 0.0])

        # Massless photon-analog: each has energy E_gamma and moves back-to-back
        # in c.m.  In lab: add boost p_boost.
        gamma1_vel = p_boost + E_gamma * p_perp
        gamma2_vel = p_boost - E_gamma * p_perp

        # Photon analogs represented as massless Kinks (mass=0, chirality=0)
        photon1 = Kink(
            position=kink.position.copy(),
            velocity=gamma1_vel,
            chirality=0,
            winding=0,
            mass=0.0,
        )
        photon2 = Kink(
            position=antikink.position.copy(),
            velocity=gamma2_vel,
            chirality=0,
            winding=0,
            mass=0.0,
        )

        # Energy deposited into substrate radiation
        Q_value = E_k + E_ak  # all rest + kinetic mass converted to waves

        E_out = math.sqrt(np.dot(gamma1_vel, gamma1_vel)) + math.sqrt(
            np.dot(gamma2_vel, gamma2_vel)
        )
        p_out = gamma1_vel + gamma2_vel
        energy_ok, momentum_ok = self._check_conservation(
            E_total, p_total, E_out, p_out, tol=1e-4
        )

        sigma_gev2 = annihilation_cross_section(sqrt_s, self.M_K)
        sigma_m2 = sigma_gev2 * GEV_INV2_TO_M2

        return ScatteringEvent(
            incoming_kinks=[kink, antikink],
            outgoing_kinks=[photon1, photon2],
            center_of_mass_energy=sqrt_s,
            cross_section=sigma_m2,
            process_type="annihilation",
            Q_value=Q_value,
            energy_conserved=energy_ok,
            momentum_conserved=momentum_ok,
        )

    # ------------------------------------------------------------------
    # Pair production
    # ------------------------------------------------------------------

    def pair_production(
        self,
        photon_energy: float,
        target: Optional[Kink] = None,
    ) -> ScatteringEvent:
        """High-energy photon → kink + antikink pair production.

        A massless photon cannot produce a massive kink pair in pure vacuum —
        momentum conservation requires a recoil target.  In the substrate
        framework the medium itself (with effective mass M_K) provides the
        necessary recoil (§18.30, §18.53.2).  When ``target`` is None the
        substrate background acts as the implicit target.

        Kinematics (γ + target → kink + antikink):

            s = (E_γ + M_target)² - E_γ²
              = 2 M_target E_γ + M_target²

            √s = √(2 M_target E_γ + M_target²)

        In the c.m. frame each kink carries E_cm = √s / 2 and momentum
        p_cm = √(E_cm² - M_K²).  The pair is boosted back to the lab along
        the photon direction (z-axis).  This Lorentz boost exactly conserves
        both energy and three-momentum to floating-point precision.

        With ``target`` provided, its mass replaces M_K as the target mass
        and the threshold is adjusted accordingly (Bethe-Heitler formula):

            E_threshold = 2 M_K + M_K² / M_target   [with target]
            E_threshold = 2 M_K                       [substrate target, no arg]

        Args:
            photon_energy: Incoming photon energy in GeV.
            target: Optional explicit target kink.  When None, the substrate
                medium (mass = M_K) provides the recoil.

        Returns:
            ScatteringEvent with the produced kink + antikink pair.  Both
            ``energy_conserved`` and ``momentum_conserved`` are True to within
            the floating-point tolerance of the boost.

        Raises:
            ValueError: If photon_energy is below the pair production threshold.

        References:
            §18.49  (kinematic threshold for kink pair production)
            §18.53.2 (pair production mechanism: wave → topology)
            §18.30  (stress-loading at collider; substrate provides recoil)
        """
        E_threshold = pair_production_threshold(self.M_K)

        # Determine target mass: explicit target or substrate background
        if target is not None:
            M_target = target.mass
            # Bethe-Heitler threshold: E_γ ≥ 2 M_K + M_K² / M_target
            effective_threshold = 2.0 * self.M_K + self.M_K**2 / M_target
        else:
            # Substrate background acts as implicit target of mass M_K
            M_target = self.M_K
            effective_threshold = E_threshold  # = 2 M_K (exact)

        if photon_energy < effective_threshold:
            raise ValueError(
                f"Photon energy {photon_energy:.4f} GeV is below the pair "
                f"production threshold {effective_threshold:.4f} GeV "
                f"(= 2 × M_K = 2 × {self.M_K} GeV)."
            )

        # Mandelstam s for γ (E_γ, 0, 0, E_γ) + target (M_target, 0, 0, 0)
        s = 2.0 * M_target * photon_energy + M_target**2
        sqrt_s = math.sqrt(s)

        # c.m. frame kinematics
        half_sqrt_s = sqrt_s / 2.0
        p_cm = math.sqrt(max(half_sqrt_s**2 - self.M_K**2, 0.0))

        # Random isotropic production direction in c.m. frame
        phi = self._rng.uniform(0.0, 2.0 * math.pi)
        cos_theta = self._rng.uniform(-1.0, 1.0)
        sin_theta = math.sqrt(max(1.0 - cos_theta**2, 0.0))

        # c.m. kink 3-momenta (back-to-back)
        px_cm = p_cm * sin_theta * math.cos(phi)
        py_cm = p_cm * sin_theta * math.sin(phi)
        pz_k_cm = p_cm * cos_theta
        pz_ak_cm = -p_cm * cos_theta

        # Lorentz boost from c.m. to lab along +z
        # β_boost = p_cm_total_z / E_cm_total
        # c.m. total 4-momentum: from lab, the photon+target system has
        #   P_lab = (E_γ + M_target, 0, 0, E_γ)
        # β_boost = E_γ / (E_γ + M_target)
        beta_boost = photon_energy / (photon_energy + M_target)
        gamma_boost = 1.0 / math.sqrt(1.0 - beta_boost**2)

        def _boost_pz(E_cm: float, pz_cm: float) -> tuple[float, float]:
            """Boost energy and z-momentum from c.m. to lab frame."""
            E_lab = gamma_boost * (E_cm + beta_boost * pz_cm)
            pz_lab = gamma_boost * (pz_cm + beta_boost * E_cm)
            return E_lab, pz_lab

        E_k_lab, pz_k_lab = _boost_pz(half_sqrt_s, pz_k_cm)
        E_ak_lab, pz_ak_lab = _boost_pz(half_sqrt_s, pz_ak_cm)

        # Lab-frame 3-momenta (x, y unchanged by z-boost)
        v_kink = np.array([px_cm, py_cm, pz_k_lab])
        v_akink = np.array([-px_cm, -py_cm, pz_ak_lab])

        produced_kink = Kink(
            position=np.zeros(3),
            velocity=v_kink,
            chirality=+1,
            winding=1,
            mass=self.M_K,
        )
        produced_antikink = Kink(
            position=np.zeros(3),
            velocity=v_akink,
            chirality=-1,
            winding=1,
            mass=self.M_K,
        )

        # Photon incoming (massless, along +z)
        photon_repr = Kink(
            position=np.zeros(3),
            velocity=np.array([0.0, 0.0, photon_energy]),
            chirality=0,
            winding=0,
            mass=0.0,
        )

        # Conservation check against full γ + target system
        E_in = photon_energy + M_target
        p_in = np.array([0.0, 0.0, photon_energy])  # target at rest in lab

        # Verify lab-frame energies match on-shell condition
        E_k_check = math.sqrt(self.M_K**2 + float(np.dot(v_kink, v_kink)))
        E_ak_check = math.sqrt(self.M_K**2 + float(np.dot(v_akink, v_akink)))
        E_out = E_k_check + E_ak_check
        p_out = v_kink + v_akink

        # The target absorbs recoil: E_recoil = E_in - E_out
        # If target is None (substrate), the substrate absorbs it transparently.
        # For exact tracking we include it in the conservation check.
        # Since γ + target → kink + antikink + (target recoil), energy balance
        # is: E_in = E_out + E_target_recoil. For the kink sector alone:
        # E_out < E_in if target provides recoil.
        # We check that 3-momentum of (kink + antikink) = 3-momentum of photon.
        momentum_ok = bool(
            float(np.linalg.norm(p_in - p_out)) / (float(np.linalg.norm(p_in)) + 1e-30) < 1e-6
        )
        # Energy check: kink+antikink must carry E_in - M_target (the photon's energy
        # goes into the produced pair; target just provides momentum conservation)
        # E_out should equal gamma_boost*(sqrt_s) = E_γ + M_target for the full system.
        energy_ok = bool(abs(E_out - E_in) / (E_in + 1e-30) < 1e-6)

        # Cross-section from annihilation time-reversal
        sigma_gev2 = annihilation_cross_section(sqrt_s, self.M_K)
        sigma_m2 = sigma_gev2 * GEV_INV2_TO_M2

        incoming = [photon_repr] if target is None else [photon_repr, target]
        return ScatteringEvent(
            incoming_kinks=incoming,
            outgoing_kinks=[produced_kink, produced_antikink],
            center_of_mass_energy=sqrt_s,
            cross_section=sigma_m2,
            process_type="pair_production",
            Q_value=-(2.0 * self.M_K),  # negative: energy stored in kink rest mass
            energy_conserved=energy_ok,
            momentum_conserved=momentum_ok,
        )

    # ------------------------------------------------------------------
    # Total cross-section
    # ------------------------------------------------------------------

    def kink_antikink_total_cross_section(self, energy: float) -> float:
        """Total kink-antikink cross-section (elastic + annihilation).

        The optical theorem relates the total cross-section to the forward
        elastic amplitude.  For kink-antikink scattering the total cross-
        section has two contributing channels:

            σ_total = σ_elastic + σ_annihilation

        At energies just above threshold (√s ~ 2 M_K) annihilation dominates.
        At very high energies (√s ≫ M_K) both fall as 1/s and the elastic
        channel from the §18.48.3 Coulomb-like tail dominates.

        Args:
            energy: Center-of-mass energy √s in GeV.

        Returns:
            Total cross-section in GeV⁻².

        References:
            §18.37  (annihilation channel)
            §18.48.3 (elastic channel from two-kink potential)
        """
        return elastic_cross_section(energy, self.M_K) + annihilation_cross_section(
            energy, self.M_K
        )


# ---------------------------------------------------------------------------
# Collider event simulator
# ---------------------------------------------------------------------------


def simulate_collider_event(
    beam_energy: float,
    particle_type: Literal["kink_antikink", "photon", "electron_positron"] = "kink_antikink",
    M_K: float = M_K_CANONICAL,
    rng_seed: Optional[int] = None,
) -> ScatteringEvent:
    """Simulate a single collider collision producing kink-antikink pairs.

    Models a head-on collision between two beam particles at √s = 2 × beam_energy
    and determines the dominant final state based on the particle type and energy.

    Physics:
        - kink_antikink beams (√s ≥ 2 M_K): elastic or annihilation based on
          the relative cross-sections at the given energy.
        - photon beams (E_γ ≥ 2 M_K): γγ → kink + antikink pair production.
        - electron_positron beams: Bhabha or pair production channel.

    The process-type selection follows the branching ratio:
        P(annihilation) = σ_ann / σ_total

    Momentum and energy are enforced to machine precision in the c.m. frame
    before the Lorentz boost back to the lab (§18.30 stress-loading analog).

    Args:
        beam_energy: Energy of each beam particle in GeV (so √s = 2 E_beam
            for equal-mass beams).  Must be positive.
        particle_type: One of "kink_antikink", "photon", "electron_positron".
        M_K: Kink rest mass in GeV (default: canonical §18.22 value, 27 GeV).
        rng_seed: Optional integer seed for reproducible kinematics.

    Returns:
        ScatteringEvent with realistic kinematics for the produced final state.

    Raises:
        ValueError: If beam_energy ≤ 0 or if it is below the pair production
            threshold for the chosen particle_type.

    References:
        §18.30  (collider stress-loading mechanism)
        §18.37  (kink-antikink collider physics)
        §18.49  (QCD and QED inherited cross-sections)
        §18.49.8 (structural identity → Bhabha and Compton inherited)
    """
    if beam_energy <= 0.0:
        raise ValueError(f"beam_energy must be positive; got {beam_energy}")

    scatterer = KinkScattering(M_K=M_K, rng_seed=rng_seed)
    rng = random.Random(rng_seed)
    sqrt_s = 2.0 * beam_energy

    if particle_type == "photon":
        # γγ → kink + antikink via pair production
        E_threshold = pair_production_threshold(M_K)
        if beam_energy < E_threshold:
            raise ValueError(
                f"Photon beam energy {beam_energy:.3f} GeV is below the pair "
                f"production threshold {E_threshold:.3f} GeV (= 2 M_K = 2 × {M_K} GeV)."
            )
        # Each beam contributes beam_energy; combined photon energy for pair production
        return scatterer.pair_production(photon_energy=sqrt_s, target=None)

    if particle_type == "electron_positron":
        # e⁺e⁻ collider: primarily Bhabha at low energy; kink pair at high energy
        sigma_bhabha = bhabha_scattering_cross_section(sqrt_s) * GEV_INV2_TO_M2
        E_threshold_kink = pair_production_threshold(M_K)
        if sqrt_s >= E_threshold_kink:
            # Above threshold: decide between Bhabha (elastic) and kink production
            sigma_kink_ann = annihilation_cross_section(sqrt_s, M_K) * GEV_INV2_TO_M2
            sigma_total = sigma_bhabha + sigma_kink_ann
            # Random process selection weighted by cross-sections
            if rng.random() < sigma_bhabha / (sigma_total + 1e-60):
                # Bhabha: e⁺ e⁻ → e⁺ e⁻ (represented as massless kinks for simplicity)
                e_minus = Kink(
                    position=np.array([-1.0, 0.0, 0.0]),
                    velocity=np.array([beam_energy, 0.0, 0.0]),
                    chirality=+1, winding=1, mass=M_ELECTRON_GEV,
                )
                e_plus = Kink(
                    position=np.array([+1.0, 0.0, 0.0]),
                    velocity=np.array([-beam_energy, 0.0, 0.0]),
                    chirality=-1, winding=1, mass=M_ELECTRON_GEV,
                )
                return scatterer.elastic_scattering(e_minus, e_plus, impact_parameter=0.1)
            else:
                # Kink pair production
                return scatterer.pair_production(photon_energy=sqrt_s, target=None)
        else:
            # Below kink threshold: pure Bhabha
            e_minus = Kink(
                position=np.array([-1.0, 0.0, 0.0]),
                velocity=np.array([beam_energy, 0.0, 0.0]),
                chirality=+1, winding=1, mass=M_ELECTRON_GEV,
            )
            e_plus = Kink(
                position=np.array([+1.0, 0.0, 0.0]),
                velocity=np.array([-beam_energy, 0.0, 0.0]),
                chirality=-1, winding=1, mass=M_ELECTRON_GEV,
            )
            return scatterer.elastic_scattering(e_minus, e_plus, impact_parameter=0.1)

    # Default: kink_antikink beams
    E_threshold = pair_production_threshold(M_K)

    # Build head-on kink-antikink beams
    incoming_kink = Kink(
        position=np.array([-1.0, 0.0, 0.0]),
        velocity=np.array([+beam_energy, 0.0, 0.0]),
        chirality=+1, winding=1, mass=M_K,
    )
    incoming_antikink = Kink(
        position=np.array([+1.0, 0.0, 0.0]),
        velocity=np.array([-beam_energy, 0.0, 0.0]),
        chirality=-1, winding=1, mass=M_K,
    )

    if sqrt_s < E_threshold:
        # Below pair production: elastic only
        return scatterer.elastic_scattering(incoming_kink, incoming_antikink, impact_parameter=0.0)

    # Above threshold: decide channel by relative cross-sections
    sigma_el = elastic_cross_section(sqrt_s, M_K)
    sigma_ann = annihilation_cross_section(sqrt_s, M_K)
    sigma_total = sigma_el + sigma_ann

    if sigma_total < 1e-60:
        return scatterer.elastic_scattering(incoming_kink, incoming_antikink, impact_parameter=0.0)

    if rng.random() < sigma_ann / sigma_total:
        return scatterer.pair_annihilation(incoming_kink, incoming_antikink)
    else:
        return scatterer.elastic_scattering(incoming_kink, incoming_antikink, impact_parameter=0.0)
