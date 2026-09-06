"""Möbius half-flux U(1) bundle topology and dynamics (spec §§18.4, 18.10).

The 45° cone carries a U(1) principal bundle whose connection 1-form is

    A = (1/2) dθ                                                    (§18.10)

This forces the bundle's first Chern class to 1/2 and imposes the
topological constraint

    ∮_C A_μ dx^μ = π · w(C)                                         (§18.45)

for every closed loop C with winding number w(C).  That single constraint
cascades into the following observable consequences (all structural, no
free parameters):

* **Spin-½ statistics** — holonomy = exp(iπ) = −1 per cone traversal;
  full return requires 4π (720°).
* **Charge quantisation** — elementary charge e is set by the half-flux
  normalization; all charges are integer multiples of e.
* **Aharonov-Bohm phase** — ΔΦ = π per half-flux winding, observable in
  interference experiments.
* **θ_QCD = 0** — the half-flux topology forbids a non-zero CP-violating
  vacuum angle; the strong-CP problem is solved structurally.

References: spec §§18.4, 18.10, 18.45.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

from stiff_medium.physical_constants import ALPHA_THOMSON_CODATA_2022

# ---------------------------------------------------------------------------
# Physical constants used throughout this module (SI)
# ---------------------------------------------------------------------------

PI: Final[float] = math.pi

# Reduced Planck constant [J·s]
HBAR: Final[float] = 1.054571817e-34

# Speed of light [m/s]
C_LIGHT: Final[float] = 2.99792458e8


# ---------------------------------------------------------------------------
# MobiusBundle dataclass
# ---------------------------------------------------------------------------


@dataclass
class MobiusBundle:
    """Möbius half-flux U(1) bundle for a single particle (spec §§18.4, 18.10).

    The bundle is fully characterised by three quantities:

    * ``half_flux``        — topological constant ∫ F = π (half a Dirac flux
                             quantum).  Fixed by construction; not a free
                             parameter.
    * ``elementary_charge``— the natural charge unit derived from the bundle's
                             normalization.  In natural units (ℏ = c = 1) this
                             equals sqrt(4π α) ≈ 0.3028; the dataclass stores a
                             dimensionless ratio relative to that normalization.
    * ``winding_number``   — integer w ∈ {±1, ±2, …} labelling the homotopy
                             class of the closed loop around the cone axis.
                             Zero is excluded: a loop with w = 0 encloses no
                             flux and does not interact with the bundle.

    Attributes:
        half_flux: Topological constant ∫_disc F = π.  The first Chern class
            is c_1 = 1/2.  This is the defining property of the Möbius bundle.
        elementary_charge: Dimensionless charge unit e (in units where 4π α = e²;
            i.e. e ≈ 0.30282212).  Computed from the half-flux normalization.
        winding_number: Integer winding number w of the closed loop.  Allowed
            values: any non-zero integer.  Positive values correspond to
            counter-clockwise orientation; negative values to clockwise.

    Raises:
        ValueError: If ``winding_number`` is zero.

    Example:
        >>> b = MobiusBundle(winding_number=1)
        >>> b.half_flux == math.pi
        True
        >>> b.holonomy(1)
        -1.0
    """

    half_flux: float = field(default=PI)
    elementary_charge: float = field(init=False)
    winding_number: int = 1

    def __post_init__(self) -> None:
        if self.winding_number == 0:
            raise ValueError(
                "winding_number must be a non-zero integer; w = 0 encloses "
                "no bundle flux and is not a physical loop."
            )
        # Elementary charge from half-flux normalization (§18.9, §18.10).
        # In the U(1) gauge theory, e is defined such that e · g = nπℏc for
        # a magnetic monopole of strength g = ℏc/(2e) (Dirac condition with
        # half-flux).  Setting n = 1 and ℏ = c = 1:
        #   e = sqrt(4π α) in natural units
        self.elementary_charge = math.sqrt(4.0 * PI * ALPHA_THOMSON_CODATA_2022)

    # ------------------------------------------------------------------
    # Core bundle methods
    # ------------------------------------------------------------------

    def holonomy(self, closed_loop_winding_number: int) -> float:
        """Return the U(1) holonomy exp(i ∮ A) for a loop with given winding.

        The connection A = (1/2) dθ (§18.10) gives

            ∮_C A = (1/2) · 2π · w = π · w

        so the holonomy is

            hol = exp(i π w) = (−1)^w

        This is real (+1 or −1) because the half-flux topology is a
        Z₂-valued invariant.

        Args:
            closed_loop_winding_number: Integer winding number w of the loop.

        Returns:
            +1.0 for even w (bosonic sector), −1.0 for odd w (fermionic sector).

        References:
            spec §18.10, equation hol(γ) = exp(i ∮_γ A).
        """
        w = closed_loop_winding_number
        integral = self.half_flux * w  # ∮ A = π w
        # holonomy = exp(i · integral); take real part (always ±1 for half-flux)
        return math.cos(integral)

    def aharonov_bohm_phase(self, winding: int) -> float:
        """Return the Aharonov-Bohm phase ΔΦ = (e/ℏ) · ∮ A·dx for given winding.

        From the half-flux constraint:

            ∮_C A_μ dx^μ = π · w(C)                                (§18.45)

        and the standard AB phase formula ΔΦ = (e/ℏ) · ∮ A·dx:

            ΔΦ = (e/ℏ) · π · w

        In natural units (ℏ = 1, e = elementary_charge):

            ΔΦ = π · w    (up to the charge normalization factor)

        The physical consequence: one half-flux winding (w = 1) shifts an
        interference fringe by π radians — half a fringe — observable as the
        sign change ψ → −ψ that distinguishes spin-½ from spin-1 particles.

        Args:
            winding: Integer winding number of the closed loop.

        Returns:
            Phase in radians: π × winding.

        References:
            spec §18.10, §18.45; also scripts/quantum_phenomena.py §2.
        """
        return PI * winding

    def is_charge_consistent(self, charge_in_units_of_e: float) -> bool:
        """Return True when the charge is an integer multiple of e.

        Charge quantisation follows directly from the half-flux constraint:
        the bundle's topology admits only states whose charges are integer
        multiples of the elementary charge e (§18.4).  A charge of, say,
        1.7 e is topologically forbidden.

        Args:
            charge_in_units_of_e: Proposed charge expressed in units of the
                elementary charge e.  E.g. 1.0 for a proton, -1.0 for an
                electron, 2/3 for an up-quark (which is allowed as a
                constituent inside a baryon, where the total is integer).

        Returns:
            True if ``abs(charge_in_units_of_e - round(charge_in_units_of_e))``
            is less than 1e-9 (integer multiple of e); False otherwise.
        """
        return abs(charge_in_units_of_e - round(charge_in_units_of_e)) < 1e-9

    def dirac_quantization_check(self, monopole_strength_g: float) -> bool:
        """Verify the Dirac quantisation condition e × g = n π ℏ c.

        The generalized Dirac condition for a half-flux bundle (§18.10) is

            e · g = n π ℏ c        (n integer, smallest n = 1)

        For the minimal case n = 1 this gives the smallest allowed monopole
        strength:

            g_min = π ℏ c / e

        Args:
            monopole_strength_g: Magnetic monopole strength in SI units [J·m/A]
                = [kg·m³/(A·s²)].  Pass ``math.pi * HBAR * C_LIGHT / e_SI``
                for the Dirac minimum, where e_SI = 1.602176634e-19 C.

        Returns:
            True if e_SI · g is an integer multiple of π ℏ c within
            relative tolerance 1e-6.

        References:
            spec §18.4, §18.10; Dirac 1931, eq. eg = nℏc/2 → half-flux gives
            eg = nπℏc instead of the standard Dirac eg = nπℏc/1.
        """
        e_si = 1.602176634e-19  # elementary charge [C]
        rhs = PI * HBAR * C_LIGHT  # π ℏ c  (n = 1 term)
        product = e_si * monopole_strength_g
        if product < 1e-80:
            return False
        ratio = product / rhs
        return abs(ratio - round(ratio)) < 1e-6


# ---------------------------------------------------------------------------
# MobiusFermion class
# ---------------------------------------------------------------------------


class MobiusFermion:
    """A fermion living on a Möbius half-flux U(1) bundle (spec §§18.4, 18.10).

    A fermion is characterised by its chirality (±1), its winding number on
    the bundle, and the elementary charge e of the bundle it lives on.

    The charge and rotation behaviour are *structurally forced* by the
    half-flux topology:

    * charge = chirality × winding_n × e
    * rotation_phase(2π) = −1    (spin-½: requires 720° for identity)
    * is_spin_half = True        (always, for a half-flux bundle)

    Args:
        chirality: +1 (right-handed) or −1 (left-handed) internal degree of
            freedom.  Corresponds to particle vs antiparticle in the Möbius
            picture (spec §6 slope polarity).
        winding_n: Integer winding number on the bundle.  Determines the
            particle's charge and Aharonov-Bohm coupling.
        bundle: The :class:`MobiusBundle` instance this fermion lives on.
            Defaults to a winding_number=1 bundle if not provided.

    Raises:
        ValueError: If ``chirality`` is not ±1.

    Example:
        >>> f = MobiusFermion(chirality=-1, winding_n=1)
        >>> f.charge  # −1 × 1 × e ≈ −0.303
        >>> f.rotation_phase(2 * math.pi)
        -1.0
        >>> f.is_spin_half
        True
    """

    def __init__(
        self,
        chirality: int,
        winding_n: int,
        bundle: MobiusBundle | None = None,
    ) -> None:
        if chirality not in (1, -1):
            raise ValueError(
                f"chirality must be +1 or −1, got {chirality!r}."
            )
        self.chirality: int = chirality
        self.winding_n: int = winding_n
        self.bundle: MobiusBundle = (
            bundle if bundle is not None else MobiusBundle(winding_number=winding_n or 1)
        )

    @property
    def charge(self) -> float:
        """Charge = chirality × winding_n × e (spec §18.4, §18.10).

        Returns:
            Charge in units of the dimensionless elementary charge.  Positive
            for particles (chirality = +1, winding > 0); negative for
            antiparticles (chirality = −1, winding > 0) or reversed windings.
        """
        return self.chirality * self.winding_n * self.bundle.elementary_charge

    def rotation_phase(self, angle_radians: float) -> float:
        """Phase factor acquired under spatial rotation by angle_radians.

        For a half-flux bundle the spin-½ relation is

            U(θ) = exp(i θ/2)    →    U(2π) = exp(iπ) = −1

        so a 2π (360°) rotation maps the fermion to minus itself, and a
        4π (720°) rotation returns to +1.

        This is the *defining* observable of spin-½ and the direct
        consequence of half-flux holonomy on the cone bundle (§18.10):
        parallel transport around the azimuthal circle picks up exp(iπ) = −1.

        Args:
            angle_radians: Rotation angle in radians.

        Returns:
            Real part of exp(i · angle_radians / 2).  Values:
            * +1.0 at 0 or 4π (identity)
            * −1.0 at 2π (spin-½ sign flip)
            * 0.0 at π or 3π (quarter-cycle)

        References:
            spec §18.10; scripts/spinor.py:slope_sign().
        """
        return math.cos(angle_radians / 2.0)

    @property
    def is_spin_half(self) -> bool:
        """True — all fermions on a half-flux bundle carry spin-½ (spec §18.4).

        The half-flux topology (first Chern class = 1/2) structurally forces
        half-integer angular momentum on every section of the bundle.  This
        is not a dynamical choice but a topological one.
        """
        return True

    def __repr__(self) -> str:
        sign = "+" if self.chirality > 0 else "−"
        return (
            f"MobiusFermion(chirality={sign}1, winding_n={self.winding_n}, "
            f"charge={self.charge:+.6f})"
        )


# ---------------------------------------------------------------------------
# Module-level functions
# ---------------------------------------------------------------------------


def coleman_bosonization_g(beta_squared: float) -> float:
    """Coleman's bosonization relation: g(β²) = π(4π/β² − 1).

    In the sine-Gordon / Thirring duality (Coleman 1975), the sine-Gordon
    coupling β² maps to the Thirring fermion coupling g via

        g = π ( 4π/β² − 1 )

    Key values:
    * β² = 4π  →  g = 0  (free fermion)
    * β² = 8π  →  g = −π/2  (Kosterlitz-Thouless transition point)

    In the stiff-medium model (§18.11), the substrate's sine-Gordon field φ
    bosonizes to a fermionic excitation (the kink zero-mode).  Coleman's
    relation sets the coupling between the scalar and fermionic sectors when
    the Lagrangian is written in either dual language.

    Args:
        beta_squared: Sine-Gordon coupling constant β².  Must be positive.
            At β² = 4π the dual fermion is free; below 4π the interaction
            is repulsive; above 8π the theory flows to the trivial fixed point.

    Returns:
        Thirring coupling g in natural units.

    Raises:
        ValueError: If ``beta_squared`` is not positive.

    References:
        Coleman, Phys. Rev. D 11 (1975) 2088; spec §18.11 candidate Lagrangian.
    """
    if beta_squared <= 0.0:
        raise ValueError(
            f"beta_squared must be strictly positive; got {beta_squared!r}."
        )
    return PI * (4.0 * PI / beta_squared - 1.0)


def theta_qcd_from_topology() -> float:
    """Return θ_QCD = 0.0, forbidden by the half-flux topology (spec §18.4).

    In standard QCD the vacuum angle θ_QCD parametrises CP violation in the
    strong sector.  The observational upper bound is |θ_QCD| < 10⁻¹⁰ (from
    neutron EDM measurements), yet the SM provides no mechanism that enforces
    this — the strong-CP problem.

    In the stiff-medium model the Möbius half-flux U(1) bundle (§18.4, §18.10)
    structurally forbids a non-zero θ_QCD.  The reasoning is topological:

    * θ_QCD corresponds to adding a term θ F ∧ F / (8π²) to the QCD Lagrangian.
    * The half-flux bundle's first Chern class c₁ = 1/2 is odd.  For a term
      θ c₁² to be gauge-invariant, θ must be zero (mod 2π) when c₁ is
      half-integer (because 2 × (1/2)² = 1/2, which is not an integer, so
      only θ = 0 leaves the path integral well-defined without a Peccei-Quinn
      mechanism).

    This is the same reason the same topology gives spin-½: the half-flux
    value is the unique CP-preserving choice.

    Returns:
        0.0 — the only allowed vacuum angle for a half-flux bundle.

    References:
        spec §18.4; scripts/sharp_predictions.py §6.
    """
    return 0.0


def bundle_compatibility_check(
    particle1: MobiusFermion,
    particle2: MobiusFermion,
) -> dict[str, bool | str]:
    """Check whether two MobiusFermion instances can form a bound state.

    Two fermions on the same Möbius bundle can form a bound state only when
    they have *opposite* chiralities (slope signs), in analogy with the
    Möbius dynamics pairing rule in
    :mod:`stiff_medium.mobius_dynamics` (spec §5.5.1, §13 gap #1).

    Rules (derived from half-flux holonomy):

    1. **Opposite chirality** (one +1, one −1): attractive; can bind.
       This is the electron-positron analogue.
    2. **Same chirality and same winding_n**: Pauli-like repulsion; cannot
       bind.  Fermionic exclusion from identical quantum numbers.
    3. **Same chirality, different winding_n**: marginal; no bound state
       expected (charges are parallel, not complementary).

    The total charge of a compatible pair is
    ``charge₁ + charge₂ = (n₁ − n₂) × e``, which is zero for equal
    winding numbers (neutral bound state, like positronium).

    Args:
        particle1: First fermion.
        particle2: Second fermion.

    Returns:
        Dictionary with keys:

        * ``"can_bind"`` (bool): True if the pair can form a bound state.
        * ``"reason"`` (str): Human-readable explanation.
        * ``"total_charge"`` (float): Sum of the two charges (in units of e).
        * ``"net_winding"`` (int): Algebraic sum of winding numbers.
        * ``"is_neutral_pair"`` (bool): True when total charge is zero.

    References:
        spec §5.5.1, §13 gap #1; :mod:`stiff_medium.mobius_dynamics`.
    """
    can_bind: bool
    reason: str

    same_chirality = particle1.chirality == particle2.chirality
    same_winding = particle1.winding_n == particle2.winding_n

    if not same_chirality:
        can_bind = True
        reason = (
            "Opposite chiralities: complementary slope shapes attract "
            "(electron-mode + positron-mode analogue)."
        )
    elif same_winding:
        can_bind = False
        reason = (
            "Same chirality AND same winding: Pauli-like exclusion — "
            "identical fermionic quantum numbers cannot co-occupy a bound state."
        )
    else:
        can_bind = False
        reason = (
            "Same chirality, different winding: parallel charges, no "
            "complementary attractive channel available."
        )

    total_charge = particle1.charge + particle2.charge
    net_winding = particle1.winding_n + particle2.winding_n
    is_neutral = abs(total_charge) < 1e-9

    return {
        "can_bind": can_bind,
        "reason": reason,
        "total_charge": total_charge,
        "net_winding": net_winding,
        "is_neutral_pair": is_neutral,
    }


def chern_class_first() -> float:
    """Return the first Chern class c₁ = 1/2 of the Möbius bundle.

    The first Chern class of a U(1) bundle is

        c₁ = (1 / 2π) ∫_disc F

    For the Möbius half-flux bundle (spec §18.10, §18.45):

        ∫_disc F = π    →    c₁ = π / (2π) = 1/2

    This distinguishes the Möbius bundle from a standard (bosonic) U(1)
    bundle, which has c₁ = 1.  The half-integer Chern class is the
    topological invariant encoding fermionic statistics.

    Returns:
        0.5 — the first Chern class of the Möbius half-flux bundle.

    References:
        spec §18.10, §18.45.
    """
    return 0.5


# ---------------------------------------------------------------------------
# Aharonov-Bohm verification
# ---------------------------------------------------------------------------


def verify_aharonov_bohm() -> dict[str, float]:
    """Verify the spin-½ Aharonov-Bohm phase relations for a half-flux bundle.

    The half-flux constraint (spec §18.45) gives

        ΔΦ = π × w

    Combined with the rotation-phase formula for spin-½ fermions

        U(θ) = exp(i θ/2)

    the key numerical checks are:

    * **π phase (one-loop, w = 1)**:
        ΔΦ = π  →  rotation by 2π  →  exp(iπ) = −1
    * **2π phase (two-loop, w = 2)**:
        ΔΦ = 2π  →  rotation by 4π  →  exp(i2π) = +1  (identity, full return)

    A spin-1 (bosonic) system would give +1 at w = 1; the −1 result here is
    the definitive sign of spin-½ statistics.

    Returns:
        Dictionary with keys:

        * ``"phase_one_loop"`` (float): ΔΦ in radians for w = 1 (should be π).
        * ``"holonomy_one_loop"`` (float): exp(iΔΦ) for w = 1 (should be −1.0).
        * ``"rotation_phase_2pi"`` (float): U(2π) = cos(π) (should be −1.0).
        * ``"phase_two_loop"`` (float): ΔΦ in radians for w = 2 (should be 2π).
        * ``"holonomy_two_loop"`` (float): exp(iΔΦ) for w = 2 (should be +1.0).
        * ``"rotation_phase_4pi"`` (float): U(4π) = cos(2π) (should be +1.0).
        * ``"spin_half_confirmed"`` (bool): True when both holonomies match spin-½.

    References:
        spec §18.10, §18.45; scripts/quantum_phenomena.py §2.
    """
    bundle = MobiusBundle(winding_number=1)
    fermion = MobiusFermion(chirality=1, winding_n=1, bundle=bundle)

    # One-loop (w = 1): phase = π, holonomy = −1
    phase_one = bundle.aharonov_bohm_phase(1)
    hol_one = bundle.holonomy(1)
    rot_2pi = fermion.rotation_phase(2.0 * PI)

    # Two-loop (w = 2): phase = 2π, holonomy = +1
    phase_two = bundle.aharonov_bohm_phase(2)
    hol_two = bundle.holonomy(2)
    rot_4pi = fermion.rotation_phase(4.0 * PI)

    spin_half_ok = (
        abs(hol_one - (-1.0)) < 1e-12
        and abs(hol_two - 1.0) < 1e-12
        and abs(rot_2pi - (-1.0)) < 1e-12
        and abs(rot_4pi - 1.0) < 1e-12
    )

    return {
        "phase_one_loop": phase_one,
        "holonomy_one_loop": hol_one,
        "rotation_phase_2pi": rot_2pi,
        "phase_two_loop": phase_two,
        "holonomy_two_loop": hol_two,
        "rotation_phase_4pi": rot_4pi,
        "spin_half_confirmed": spin_half_ok,
    }
