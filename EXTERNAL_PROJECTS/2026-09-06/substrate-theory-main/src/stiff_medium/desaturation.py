"""Substrate de-saturation phase transition — §18.44.

The CMB is interpreted as the latent heat released when the universe-wide
substrate transitions from the fully saturated state (σ = ½, the Big Bang
initial condition) to the ordinary elastic regime (σ < ½, where structure
can form).

Per §18.44: this is a first-order phase transition in the medium.
    - Before: σ = ½ everywhere, plastic-saturated, no observers, no clocks.
    - During: phase boundary sweeps; latent energy releases as radiation.
    - After: σ < ½, kinks crystallise, atoms form, structure begins.

The energy density released is:
    ε_latent = ½ K (σ_initial² - σ_final²)

where K = c⁷ / (ℏ G²) is the Planck stress — the substrate stiffness at
the natural unit set by quantum gravity.

After the transition the radiation redshifts as a⁻⁴.  Matching today's
observed T_CMB = 2.7255 K constrains the expansion factor since the
transition, placing the event near the Planck epoch.

The bubble-nucleation rate governing how the boundary propagates is:
    Γ ∝ exp(-S / ℏ),    S = 27π² ε_latent ξ⁴ / (2 Δε²)
in the thin-wall approximation (Coleman-Callan bounce action), where Δε
is the energy difference across the phase boundary.

References:
    - Spec §18.44 (CMB as de-saturation phase transition)
    - Spec §18.39 (σ_max = ½ saturation barrier)
    - Spec §18.40 (cosmological initial condition σ = ½)
    - Coleman (1977), Phys. Rev. D 15, 2929 (bounce action)
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

#: Speed of light in vacuum (m s⁻¹).
c: Final[float] = 2.998e8

#: Reduced Planck constant (J s).
hbar: Final[float] = 1.055e-34

#: Newton's gravitational constant (m³ kg⁻¹ s⁻²).
G: Final[float] = 6.674e-11

#: Boltzmann constant (J K⁻¹).
k_B: Final[float] = 1.381e-23

#: Stefan-Boltzmann constant (W m⁻² K⁻⁴).
sigma_SB: Final[float] = 5.670e-8

#: Radiation constant a = 4σ_SB / c (J m⁻³ K⁻⁴).
a_rad: Final[float] = 4.0 * sigma_SB / c

#: Planck stress K_Planck = c⁷ / (ℏ G²)  (Pa = J m⁻³).
#: This is the substrate stiffness at the quantum-gravity scale (§18.21).
K_PLANCK: Final[float] = c**7 / (hbar * G**2)

#: Universal elastic saturation limit (§18.39).
SIGMA_MAX: Final[float] = 0.5


# ---------------------------------------------------------------------------
# Dataclass
# ---------------------------------------------------------------------------


@dataclass
class DesaturationTransition:
    """First-order phase transition of the substrate from σ=½ to σ<½ (§18.44).

    Models the de-saturation event that releases the CMB latent heat.  The
    substrate begins in the fully saturated plastic regime (σ_initial = ½)
    and evolves toward some final strain σ_final once the phase boundary has
    swept through.

    Attributes:
        sigma_initial: Strain before the transition.  Default ½ (saturated,
            the Big Bang initial condition per §18.40 and §18.42).
        sigma_final: Strain after the transition.  Default 0.0 (empty
            vacuum; the extreme case).  Realistic values are 0 < σ_final ≪ ½.
        K: Substrate stiffness modulus K (Pa).  Default is the Planck stress
            c⁷ / (ℏ G²), the natural scale at which the substrate's elastic
            limit is reached (§18.21).
        phi_max: Saturation field cutoff φ_max (dimensionless).  Per §18.39
            the saturation barrier potential V(φ) diverges at |φ| = φ_max.
            Default 1.0.
    """

    sigma_initial: float = SIGMA_MAX        # saturated: σ = ½
    sigma_final: float = 0.0               # empty vacuum: σ = 0
    K: float = field(default_factory=lambda: K_PLANCK)   # Planck stress
    phi_max: float = 1.0                   # saturation field cutoff

    # ------------------------------------------------------------------
    # Methods
    # ------------------------------------------------------------------

    def latent_heat_density(self) -> float:
        """Energy density released per unit volume during the phase transition.

        Per §18.44, the latent heat is the elastic strain energy stored in
        the saturated medium that is liberated when σ drops from σ_initial
        to σ_final:

            ε_latent = ½ K (σ_initial² - σ_final²)

        This is exactly analogous to the latent heat of a first-order
        thermodynamic phase transition (water → ice), but for the substrate
        transitioning from its plastic-saturated to ordinary-elastic regime.

        Returns:
            Latent heat energy density in J m⁻³.  Always non-negative for
            σ_initial ≥ σ_final ≥ 0.
        """
        return 0.5 * self.K * (self.sigma_initial**2 - self.sigma_final**2)

    def evolve(self, dt: float, current_sigma: float) -> float:
        """Advance the transition by one time step dt.

        The substrate strain relaxes toward σ_final at a rate set by the
        elastic relaxation timescale τ = 1 / (K / (ρ c)):

            dσ/dt = -(σ - σ_final) / τ

        For a Planck-scale K this relaxation is essentially instantaneous
        on cosmological timescales.  Here we use a phenomenological form:

            σ_next = σ_final + (σ - σ_final) exp(-dt / τ_relax)

        where τ_relax = ℏ / (K_PLANCK × ξ³) with ξ the Planck length.

        Args:
            dt: Time step in seconds.
            current_sigma: Current local strain σ.

        Returns:
            Updated strain σ_next after the time step.
        """
        xi_planck = np.sqrt(hbar * G / c**3)           # Planck length ~ 1.6e-35 m
        tau_relax = hbar / (self.K * xi_planck**3)     # relaxation timescale (s)
        decay = np.exp(-dt / tau_relax)
        return self.sigma_final + (current_sigma - self.sigma_final) * decay

    def is_complete(self) -> bool:
        """Return True if the phase transition has reached equilibrium.

        The transition is considered complete when σ_initial has already
        relaxed to within a small tolerance of σ_final, or equivalently
        when the remaining latent heat is negligible compared with the
        initial latent heat.

        Returns:
            True when |σ_initial - σ_final| < 1e-6.
        """
        return abs(self.sigma_initial - self.sigma_final) < 1.0e-6

    def nucleation_rate(self, temperature: float) -> float:
        """Bubble nucleation rate per unit volume Γ (m⁻³ s⁻¹).

        In the thin-wall (Coleman-Callan) approximation a bubble of the
        new (σ < ½) phase nucleates inside the old saturated phase with a
        rate:

            Γ = A exp(-S_E / ℏ)

        where S_E is the Euclidean bounce action and A is a pre-exponential
        factor of order K / ℏ (Planck units).  For a quartic potential the
        bounce action is:

            S_E = 27π² σ²_wall ξ⁴ K² / (2 Δε²)

        with Δε = ε_latent and σ_wall the surface tension of the bubble
        wall.  Here we approximate σ_wall ~ K × ξ (energy / area):

            σ_wall ≈ K × ξ_planck
            Δε = latent_heat_density()

        At non-zero temperature T the thermal correction modifies the action:
            S_eff = S_E / (1 + k_B T / Δε)

        Args:
            temperature: Local temperature of the medium (K).  Set to 0 for
                pure quantum tunnelling.

        Returns:
            Nucleation rate in m⁻³ s⁻¹.  Returns 0.0 if Δε ≤ 0 (no energy
            difference to drive the transition).
        """
        delta_epsilon = self.latent_heat_density()
        if delta_epsilon <= 0.0:
            return 0.0

        xi_planck = np.sqrt(hbar * G / c**3)    # Planck length (m)
        sigma_wall = self.K * xi_planck          # bubble wall surface tension (J m⁻²)

        # Bounce action (Coleman 1977, thin-wall).
        # Work entirely in log-space to avoid floating-point overflow at
        # Planck-scale energies where S_bounce / ℏ >> 1.
        #
        #   log(S_bounce) = log(27 π²) + 4 log(σ_wall) - log(2) - 3 log(Δε)
        #   log(S_eff)    = log(S_bounce) - log(1 + k_B T / Δε)
        log_S_bounce = (
            np.log(27.0 * np.pi**2)
            + 4.0 * np.log(sigma_wall)
            - np.log(2.0)
            - 3.0 * np.log(delta_epsilon)
        )

        # Thermal correction: reduces effective action
        thermal_correction = 1.0 + k_B * temperature / delta_epsilon
        log_S_eff = log_S_bounce - np.log(thermal_correction)

        # Pre-exponential factor ~ K / ℏ (natural Planck-scale rate)
        prefactor = self.K / hbar

        # Exponent = S_eff / ℏ.  Guard against overflow at Planck scales.
        log_exponent = log_S_eff - np.log(hbar)
        exponent = min(log_exponent, 700.0)   # exp(700) ~ 10^304, near float max

        return prefactor * np.exp(-exponent)


# ---------------------------------------------------------------------------
# Module-level functions
# ---------------------------------------------------------------------------


def cmb_temperature_from_latent_heat(
    latent_heat_J_per_m3: float,
    expansion_factor: float,
) -> float:
    """Predict the observed CMB temperature from the phase-transition latent heat.

    Per §18.44, the phase transition releases energy density ε_latent at the
    moment of de-saturation.  That radiation subsequently redshifts as the
    universe expands.  Radiation energy density scales as:

        ε_radiation ∝ a⁻⁴  (photon number dilutes as a⁻³, each photon
                              loses energy ∝ a⁻¹ due to redshift)

    equivalently T ∝ a⁻¹.  Setting the Stefan-Boltzmann relation at the
    transition epoch:

        ε_latent = a_rad × T_transition⁴

    and redshifting:

        T_CMB = T_transition / expansion_factor

    Args:
        latent_heat_J_per_m3: Energy density released at the phase transition
            (J m⁻³).
        expansion_factor: Ratio a_today / a_transition of the cosmological
            scale factor.  Approximately 10²⁸ for a Planck-scale transition
            redshifted to the observed CMB temperature.

    Returns:
        Predicted CMB temperature today in Kelvin.
    """
    # Invert a_rad T⁴ = ε_latent to get T at transition epoch
    T_transition = (latent_heat_J_per_m3 / a_rad) ** 0.25

    # Redshift: T scales as 1/a
    T_cmb = T_transition / expansion_factor
    return T_cmb


def phase_boundary_velocity(K: float, rho: float) -> float:
    """Speed of the phase boundary (bubble wall) as a substrate sound speed.

    Per §18.44, the de-saturation phase boundary propagates at the substrate's
    natural wave speed, which from §4 is:

        c_s = sqrt(K / ρ)

    For the physical substrate this equals c (the speed of light), since
    c² = K / ρ is one of the fundamental derived relations (spec §4).
    Passing different K and ρ values allows testing whether a candidate set
    of substrate parameters is self-consistent with c.

    Args:
        K: Substrate stiffness (Pa).
        rho: Substrate mass density (kg m⁻³).

    Returns:
        Phase boundary velocity in m s⁻¹.  Equals c when K / ρ = c².
    """
    return float(np.sqrt(K / rho))


def tensor_to_scalar_ratio_from_transition(epsilon_slow_roll: float) -> float:
    """Predict the tensor-to-scalar ratio r for inflationary observables.

    The CMB B-mode polarisation encodes primordial tensor modes generated
    during the de-saturation transition.  In slow-roll inflation the
    consistency relation is:

        r = 16 ε

    where ε is the slow-roll parameter quantifying how rapidly the inflaton
    (here: the substrate strain field σ) changes during the transition.

    Per §18.44 the CMB B-modes originate from the tensor spectrum of the
    first-order de-saturation phase transition.  For a slow-rolling
    transition r is small (~10⁻³ per spec §18.56.8); for a rapid transition
    the tensor spectrum is larger.

    Using the standard inflationary consistency relation (inherited by our
    model since the same Lagrangian §18.11 underlies both the kink dynamics
    and the scalar-driven expansion):

        r = 16 ε_slow_roll

    Args:
        epsilon_slow_roll: Slow-roll parameter ε = -dH/dt / H² (dimensionless).
            For a de-saturation transition driven by the saturation potential
            §18.39, ε ~ (K σ_max² / V_total) where V_total is the total
            vacuum energy.  Values near 0 correspond to very gradual
            transitions; ε ~ 0.01 gives r ~ 0.16 (near current Planck/BICEP
            upper bound of r < 0.036).

    Returns:
        Tensor-to-scalar ratio r (dimensionless).  Currently constrained
        observationally to r < 0.036 (BICEP/Keck 2021).
    """
    return 16.0 * epsilon_slow_roll
