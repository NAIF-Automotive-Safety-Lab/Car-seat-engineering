"""Black hole thermodynamics in the stiff-medium (saturated substrate) framework.

Black holes are regions where the substrate has reached its elastic limit:
σ = ½ (universal saturation, per §18.39). There is no singularity — inside the
horizon the medium is uniformly at σ = ½, analogous to a saturated elastic solid.

Key identifications (spec §§18.39, 18.42, 18.51, 18.55):
  - Black hole interior    → saturated substrate region (σ = ½ everywhere)
  - Event horizon          → σ = ½ surface boundary
  - Hawking radiation      → thermal fluctuations leaking across saturated boundary
  - Bekenstein-Hawking S   → counting Möbius half-flux configurations on horizon
  - Info paradox resolution→ information stored on saturated boundary (holographic)

All quantities in SI units throughout. Constants are pulled from scipy when
available, otherwise from CODATA 2018 values embedded here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# SI constants (CODATA 2018 / IAU 2015)
# ---------------------------------------------------------------------------

#: Speed of light in vacuum (m/s)
C: Final[float] = 2.99792458e8

#: Newtonian gravitational constant (m³ kg⁻¹ s⁻²)
G: Final[float] = 6.67430e-11

#: Reduced Planck constant ℏ (J·s)
HBAR: Final[float] = 1.054571817e-34

#: Boltzmann constant k_B (J/K)
K_B: Final[float] = 1.380649e-23

#: Solar mass (kg)  — IAU 2015
M_SUN: Final[float] = 1.98892e30

#: Planck length ℓ_P = √(ℏG/c³) (m)
PLANCK_LENGTH: Final[float] = math.sqrt(HBAR * G / C**3)

#: Universal saturation value σ at every horizon (spec §18.39)
SIGMA_AT_HORIZON: Final[float] = 0.5


# ---------------------------------------------------------------------------
# Dataclass
# ---------------------------------------------------------------------------


@dataclass
class BlackHole:
    """A Schwarzschild, Kerr, or Reissner-Nordström black hole.

    In the stiff-medium framework the interior is a uniformly saturated
    substrate region (σ = ½) — no singularity, no information loss.

    Attributes:
        mass_kg: Gravitational mass in kilograms.
        spin_parameter: Dimensionless spin a/M (0 ≤ spin_parameter ≤ 1).
            For spin_parameter = 0 the hole is Schwarzschild.
        charge_C: Electric charge in Coulombs (0 for uncharged).

    Note:
        Hawking temperature and evaporation formulae use the Schwarzschild
        (spin = 0, charge = 0) result. For Kerr / RN holes these differ;
        the Schwarzschild expression is returned when spin_parameter != 0 or
        charge_C != 0 with a comment in the docstrings below.

    Spec refs: §§18.39, 18.42, 18.51, 18.55.
    """

    mass_kg: float
    spin_parameter: float = 0.0   # a/M, dimensionless; 0 = Schwarzschild
    charge_C: float = 0.0         # Coulombs; 0 = uncharged

    # -----------------------------------------------------------------------
    # Geometric / thermodynamic properties
    # -----------------------------------------------------------------------

    @property
    def schwarzschild_radius(self) -> float:
        """Schwarzschild radius r_s = 2GM/c² (metres).

        For a Schwarzschild BH this is the horizon radius.  For Kerr or
        Reissner-Nordström holes the actual horizon is at a different radius,
        but r_s is still a useful scale.

        In the substrate model: r_s is the locus where σ(r) = GM/(r c²) = ½,
        i.e. the saturated boundary (spec §18.39).
        """
        return 2.0 * G * self.mass_kg / C**2

    @property
    def horizon_area(self) -> float:
        """Horizon surface area A = 4π r_s² (m²).

        Uses the Schwarzschild horizon radius.  For a rotating (Kerr) hole the
        area is smaller by a factor depending on spin_parameter, but the
        Schwarzschild formula is used here as the leading-order result.

        Spec §18.51: entropy is proportional to this area via the Möbius
        half-flux count.
        """
        rs = self.schwarzschild_radius
        return 4.0 * math.pi * rs**2

    @property
    def bekenstein_hawking_entropy(self) -> float:
        """Bekenstein-Hawking entropy S_BH = A k_B c³ / (4 G ℏ) (J/K).

        Equivalently S_BH = A / (4 ℓ_P²) in units of k_B.

        In the substrate framework (spec §18.51) this counts the number of
        distinguishable Möbius half-flux configurations that can tile the
        σ = ½ horizon surface, each carrying one bit of information.
        """
        return (self.horizon_area * K_B * C**3) / (4.0 * G * HBAR)

    @property
    def hawking_temperature(self) -> float:
        """Hawking temperature T_H = ℏc³ / (8π G M k_B) (Kelvin).

        Schwarzschild formula.  Physical interpretation in the substrate model
        (spec §18.39): thermal fluctuations of the medium at the saturated
        boundary leak outward as radiation.  The surface gravity κ = c⁴/(4GM)
        sets the characteristic frequency of these boundary oscillations.
        """
        return HBAR * C**3 / (8.0 * math.pi * G * self.mass_kg * K_B)

    @property
    def hawking_evaporation_time(self) -> float:
        """Evaporation timescale t_evap = 5120 π G² M³ / (ℏ c⁴) (seconds).

        Time for a Schwarzschild black hole to radiate away its entire mass
        via Hawking emission, assuming no accretion.

        Spec §18.55: eventual evaporation is consistent with the saturated
        boundary gradually losing mass through thermal leakage.
        """
        return (5120.0 * math.pi * G**2 * self.mass_kg**3) / (HBAR * C**4)

    @property
    def sigma_at_horizon(self) -> float:
        """Substrate strain σ at the horizon surface.

        Universally equal to ½ for any black hole mass, spin, or charge
        (spec §18.39).  This is the defining condition of the horizon in the
        stiff-medium framework:

            σ(r_s) = GM / (r_s c²) = GM / (2GM/c² · c²) = ½

        Returns:
            0.5 — the universal substrate saturation value.
        """
        return SIGMA_AT_HORIZON

    # -----------------------------------------------------------------------
    # Methods
    # -----------------------------------------------------------------------

    def mass_loss_rate(self) -> float:
        """Hawking mass-loss rate dM/dt (kg/s) — negative for evaporating BH.

        From differentiating t_evap ∝ M³ with respect to M:

            dM/dt = -ℏ c⁴ / (15360 π G² M²)

        The negative sign indicates mass loss.  For stellar-mass BHs this rate
        is vanishingly small; it becomes significant only near the final
        evaporation when M → 0.

        Returns:
            Mass loss rate in kg/s (always ≤ 0).
        """
        return -(HBAR * C**4) / (15360.0 * math.pi * G**2 * self.mass_kg**2)

    def evaporate(self, dt: float) -> "BlackHole":
        """Step BH evaporation forward by a time interval dt (seconds).

        Applies the Hawking mass-loss rate as a first-order Euler step.
        Returns a new BlackHole with updated mass; does not mutate self.

        Args:
            dt: Time step in seconds.  Must be positive.

        Returns:
            New BlackHole with mass_kg = max(0, self.mass_kg + dM/dt * dt).
            Mass is clamped to zero; spin and charge are preserved.

        Note:
            For astrophysical BHs (M ~ M_sun) the mass-loss rate is
            ~10⁻⁴⁵ kg/s, far below float64 resolution relative to the mass.
            Observable evolution requires dt >> 10⁷⁵ s.  Use the analytic
            result M(t) = M₀ (1 - t/t_evap)^(1/3) for long-timescale studies.

        Raises:
            ValueError: If dt is not positive.
        """
        if dt <= 0:
            raise ValueError(f"Time step dt must be positive; got {dt}")
        new_mass = self.mass_kg + self.mass_loss_rate() * dt
        new_mass = max(new_mass, 0.0)
        return BlackHole(
            mass_kg=new_mass,
            spin_parameter=self.spin_parameter,
            charge_C=self.charge_C,
        )

    def is_evaporated(self, threshold_kg: float = 1e10) -> bool:
        """Check whether the black hole has fully evaporated.

        A BH is considered evaporated when its mass falls below threshold_kg.
        The default threshold (10¹⁰ kg) is well above Planck mass (~2 × 10⁻⁸ kg)
        but tiny relative to any astrophysical BH.  In the substrate model,
        the saturated region vanishes when too little mass remains to maintain
        a stable σ = ½ boundary.

        Args:
            threshold_kg: Mass below which the BH is considered evaporated.

        Returns:
            True if mass_kg <= threshold_kg.
        """
        return self.mass_kg <= threshold_kg

    def qnm_frequency(self) -> complex:
        """Dominant quasinormal mode (QNM) frequency ω (rad/s) for ringdown.

        Uses the Leaver / Nollert fit for the fundamental l=2 mode of a
        Schwarzschild black hole:

            ω M = 0.3737 - 0.0890i   (in units where G = c = 1)

        Converting to SI: ω = (0.3737 - 0.0890j) × c³ / (G M)

        This governs the gravitational-wave ringdown after a BH merger; the
        real part gives the oscillation frequency, the imaginary part gives
        the damping.  The physical frequency in Hz is Re(ω) / (2π).

        Spec §18.51: the QNM is set by the geometry of the saturated boundary;
        the ringdown does not require a singularity, only the horizon structure.

        Returns:
            Complex angular frequency ω (rad/s).  Re > 0; Im < 0 (decaying).
        """
        omega_natural = complex(0.3737, -0.0890)        # in G=c=1 units
        scale = C**3 / (G * self.mass_kg)               # rad/s per unit ω
        return omega_natural * scale

    def qnm_damping_time(self) -> float:
        """QNM e-folding damping time τ = 1 / |Im(ω)| (seconds).

        The ringdown amplitude decays as exp(-|Im(ω)| t).  A larger black hole
        rings for longer.

        Returns:
            Damping time in seconds.
        """
        omega = self.qnm_frequency()
        return 1.0 / abs(omega.imag)

    def __repr__(self) -> str:
        return (
            f"BlackHole(mass={self.mass_kg:.4e} kg, "
            f"r_s={self.schwarzschild_radius:.4e} m, "
            f"T_H={self.hawking_temperature:.4e} K)"
        )


# ---------------------------------------------------------------------------
# Module-level functions
# ---------------------------------------------------------------------------


def formation_threshold_density(target_mass: float) -> float:
    """Minimum mean density (kg/m³) inside the Schwarzschild radius for BH formation.

    For a uniform sphere of mass M and radius r_s = 2GM/c², the mean density is:

        ρ_threshold = M / (4π/3 × r_s³)
                    = 3 c⁶ / (32π G³ M²)

    This shows the famous inverse-square scaling: more massive black holes form
    at lower densities because their Schwarzschild radius grows linearly with M
    while volume grows as M³.

    In the stiff-medium model (spec §18.39) this is the density at which the
    substrate strain at the boundary of the region reaches σ = ½.  The table in
    §18.39 is reproduced by evaluating this function at Earth, Sun, Sgr A*, M87
    masses.

    Args:
        target_mass: Mass of the object in kg.

    Returns:
        Critical density in kg/m³.
    """
    rs = 2.0 * G * target_mass / C**2
    volume = (4.0 / 3.0) * math.pi * rs**3
    return target_mass / volume


def holographic_information_capacity(area_m2: float) -> float:
    """Maximum information (bits) storable on a surface of area area_m2.

    From the Bekenstein bound / holographic principle, one bit per Planck area:

        N_bits = A / (4 ℓ_P²)

    In the substrate model (spec §18.51) each Planck-area cell on the σ = ½
    boundary can host one Möbius half-flux configuration (flipped or unflipped),
    giving one bit of information.  This is the mechanism by which infalling
    matter's information is preserved on the horizon rather than lost at a
    singularity — there is no singularity; the information sits on the
    saturated boundary surface.

    Args:
        area_m2: Surface area in m².

    Returns:
        Maximum information capacity in bits (float; may be very large).
    """
    return area_m2 / (4.0 * PLANCK_LENGTH**2)


def merger_radiated_energy(m1_kg: float, m2_kg: float,
                            efficiency: float = 0.05) -> float:
    """Gravitational-wave energy radiated during a BH merger (Joules).

    Numerical relativity simulations of equal-mass mergers show ~5 % of the
    total system mass radiated as gravitational waves for comparable-mass
    mergers.  GW150914 radiated ~3 M_☉ ≈ 5 % of the 65 M_☉ initial mass.

    In the substrate model the radiated energy is the kinetic energy of the
    medium's bulk elastic oscillations excited when two saturated regions
    collide and merge into a single larger saturated region.

    Args:
        m1_kg: Mass of the first black hole in kg.
        m2_kg: Mass of the second black hole in kg.
        efficiency: Fraction of total mass-energy radiated (default 0.05 = 5 %).

    Returns:
        Radiated energy in Joules (E = efficiency × (M1 + M2) c²).
    """
    total_mass = m1_kg + m2_kg
    return efficiency * total_mass * C**2


def info_paradox_resolution() -> str:
    """Explanation of how information is preserved on the saturated boundary.

    Returns:
        Multi-line string describing the substrate-mechanical resolution of
        the black hole information paradox (spec §§18.39, 18.51, 18.55).
    """
    return (
        "BLACK HOLE INFORMATION PARADOX — RESOLUTION IN THE SUBSTRATE MODEL\n"
        "=====================================================================\n"
        "\n"
        "Standard GR: a singularity at r = 0 destroys information.  Infalling\n"
        "  matter reaches the singularity and vanishes; Hawking radiation is\n"
        "  thermal and carries no information about what fell in.  The paradox\n"
        "  is that unitarity of quantum mechanics requires information to be\n"
        "  preserved, but GR says it is destroyed.\n"
        "\n"
        "Stiff-medium model (spec §§18.39, 18.51, 18.55):\n"
        "\n"
        "  1. NO SINGULARITY.  The substrate has a finite elastic limit:\n"
        "     σ_max = ½.  As matter falls in, σ rises toward ½ but can never\n"
        "     exceed it.  The medium saturates uniformly; there is no r = 0\n"
        "     singularity, only a σ = ½ interior.\n"
        "\n"
        "  2. INFORMATION ON THE BOUNDARY.  Infalling matter deposits its\n"
        "     degrees of freedom as Möbius half-flux configurations on the\n"
        "     σ = ½ boundary surface (the horizon shell).  Each Planck-area\n"
        "     cell stores one bit.  Total capacity = A / (4 ℓ_P²) bits\n"
        "     (the Bekenstein-Hawking entropy, spec §18.51).\n"
        "\n"
        "  3. HAWKING RADIATION CARRIES INFORMATION.  Because the boundary is\n"
        "     a structured surface (not a singularity), thermal fluctuations\n"
        "     that escape as Hawking radiation are correlated with the boundary\n"
        "     state.  The radiation is not strictly thermal; it is thermal-plus-\n"
        "     correlations.  The correlations carry the information out slowly\n"
        "     over the evaporation timescale.\n"
        "\n"
        "  4. UNITARITY PRESERVED.  No information is destroyed because the\n"
        "     'singularity' never forms.  The process is: matter in → boundary\n"
        "     state encoded → Hawking radiation out (correlated with boundary)\n"
        "     → BH evaporates → information fully returned to the exterior.\n"
        "     This is analogous to burning a book: the information is present\n"
        "     in the ash and smoke, just highly scrambled.\n"
        "\n"
        "  5. HOLOGRAPHIC IN THE STRICT SENSE.  The interior (3D volume) is\n"
        "     described completely by the boundary (2D surface) because the\n"
        "     interior is featureless saturated medium.  The holographic\n"
        "     principle emerges mechanically, without invoking AdS/CFT or\n"
        "     entanglement entropy arguments.\n"
        "\n"
        "Spec commitment: §18.51 stakes the model on this resolution.\n"
        "Falsification: if future GW observations show evidence of singular\n"
        "  merger behavior inconsistent with saturated-boundary ringdown,\n"
        "  this mechanism is ruled out.\n"
    )
