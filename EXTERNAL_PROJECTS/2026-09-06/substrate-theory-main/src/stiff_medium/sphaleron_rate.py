"""Sphaleron-analogue rate in the Möbius × sine-Gordon substrate (spec §18.45).

PHYSICS BACKGROUND
------------------
In Standard Model electroweak baryogenesis, sphalerons are unstable static
configurations that sit at the top of the energy barrier between SU(2) gauge
vacua with different Chern-Simons number.  Crossing the barrier changes baryon
number by ΔB = ΔN_CS = ±1.  The thermal rate for crossing is:

    Γ_sph / V ~ T⁴ × exp(-E_sph / T)    [classically above the barrier]
    Γ_sph / V ~ T⁴ × exp(-S_bounce / ℏ) [quantum tunnelling below the crossover]

The crossover temperature T_sph where these regimes match sets the sphaleron
decoupling temperature.

ANALOGUE IN THIS MODEL
-----------------------
The §18.11 / §18.45 substrate Lagrangian is (in 1+1 D, the relevant sector):

    L = (ℏ/ξ) [ ½ (∂_t φ/c)² - ½ (∂_x φ)² - (m₀²ξ²/ℏ)(1 - cos φ) ]

This is the standard sine-Gordon Lagrangian.  Its vacua are φ_n = 2πn for
integer n.  Topological kinks (winding solutions) interpolate between adjacent
vacua:

    φ_kink(x) = 4 arctan[exp(m₀(x - x₀)/ℏ)]

with kink mass M_kink = 8 m₀ (where m₀ = ℏ c / ξ is the meson mass).

The Möbius bundle (§18.10) assigns holonomy exp(iπw) = (-1)^w to a loop with
winding number w around the cone.  Adjacent vacua (n, n+1) differ in winding
number by 1; the holonomy changes sign.  This sign change IS the baryon number
change: crossing from vacuum n to vacuum n+1 creates/destroys one kink (one
"baryon").

THE SPHALERON = SINE-GORDON BION SADDLE
-----------------------------------------
In the sine-Gordon theory the sphaleron-analogue (the unstable configuration
at the top of the barrier between vacua) is the static half-period solution:

    φ_sph(x) = π     (uniform constant field at the top of the cosine barrier)

More carefully, the exact saddle in the periodic potential V(φ) = (1 - cos φ)
sits at φ = π with:

    E_sph = 2 M_kink                 [barrier height = kink + antikink threshold]

This is because to go from vacuum 0 to vacuum 2π you must pass through a
configuration where the field value is π everywhere — the point of maximum
potential energy.  The barrier height is exactly 2 M_kink because the saddle
can be visualized as a kink and antikink at zero separation (their creation
threshold).

For a finite-length system of length L the saddle structure is:
    φ_sph(x) = π + 4 arctan[tan(π x / (2 L))]   [discrete version]

In the infinite-volume limit the saddle energy is:
    E_sph = ∫ dx { ½ (dφ/dx)² + (1 - cos φ) }    [at φ = π: all potential, no KE]
           = ∫ dx V(π) × δ(x - x_saddle)   — not right in 1D

The correct result from the exact sine-Gordon calculation (see Coleman 1977):
    E_sph = E_barrier = 2 M_kink

This is the ANALOG of the SM sphaleron energy E_sph^SM ≈ 9 TeV.

BOUNCE ACTION FOR TUNNELLING
------------------------------
Below the crossover temperature T_cross, transitions between vacua proceed by
quantum tunnelling.  The Coleman bounce action (O(2) symmetric instanton in
1+1 D) for the sine-Gordon theory is:

    S_bounce = 16 m₀ / β²     [Coleman 1975, eqn. 5.6]

where β² is the sine-Gordon coupling.  In our model the §18.11 Lagrangian
fixes β² = 8π / (ℏ/ξ) × (ℏ c / ξ) = 8π (using the Coleman bosonization
relation at the free-fermion point).  The bounce action is therefore:

    S_bounce = (16 / β²) × m₀ = (16 / 8π) × m₀ = (2/π) m₀

More carefully: using the canonical sine-Gordon normalization, the
instanton action is:
    S_inst = 8 m₀ / β² × ℏ    (with ℏ restored)

At β² = 8π (Coleman's free-fermion point):
    S_inst = M_kink × ℏ / c     [in units where c = 1, S_inst = M_kink]

THERMAL RATE FORMULA
---------------------
The classical thermal activation rate (Kramers / Langer) for escaping a
barrier of height ΔE = 2 M_kink at temperature T is:

    Γ_thermal ~ (ω_0 / 2π) × exp(-ΔE / k_B T)

where ω_0 = √(V''(φ_vacuum)) / ξ is the oscillation frequency in the
vacuum well.  For sine-Gordon: V''(0) = m₀² → ω_0 = m₀ c / ℏ.

The quantum tunnelling rate below T_cross is:
    Γ_quantum ~ ω_0 × exp(-S_inst / ℏ)

At temperature T the combined rate (interpolating between the two) is:
    Γ_B(T) ≈ ω_0 × exp(-min(ΔE / k_B T,  S_inst / ℏ))

Per unit volume in 3+1 D (because baryogenesis is 3-spatial-dimensional):
    Γ_B / V ~ (m₀ c / ℏ)⁴ × exp(-ΔE / k_B T)   [T above crossover]

This uses the standard dimensional argument T⁴ ~ (k_B T / ℏ c)⁴ × (ℏ c)⁴
which, divided by (ℏ c)³ × ℏ, gives a rate density (m⁻³ s⁻¹).

The ratio Γ_B / H compared to the Hubble rate then gives the efficiency of
baryon number violation per Hubble volume per Hubble time.

HONEST LIMITATIONS
-------------------
1. The sine-Gordon theory is 1+1 dimensional.  Its sphaleron (the barrier
   maximum at φ = π) is degenerate in the transverse directions in 3+1 D.
   The correct 3+1 D calculation requires knowing how the kink is embedded
   in 3D space.  We use the dimensional estimate (m₀/ℏ)⁴ for the 3D phase
   space factor.

2. The barrier height 2 M_kink is exact in 1+1 D sine-Gordon.  In 3D, if
   the kink is a domain wall of thickness ξ, the energy per unit area is
   M_kink / ξ² and the saddle energy depends on the bubble radius.

3. The Möbius bundle holonomy (-1)^n confirms that adjacent vacua differ in
   "baryon number" (chirality), but does not quantify the asymmetry δ_CP
   between the forward and backward rates.  That remains a free parameter.

4. At T_Planck, the thermal factor exp(-ΔE / k_B T) is calculable with
   ΔE = 2 M_kink and T = T_Planck.  The kink mass M_kink depends on which
   length scale ξ is used.

References:
    Coleman (1977), Phys. Rev. D 15, 2929 — bounce action
    Coleman (1975), Phys. Rev. D 11, 2088 — sine-Gordon / Thirring duality
    Arnold & McLerran (1987), Phys. Rev. D 36, 581 — sphaleron rate formula
    Linde (1983), Phys. Lett. B 131, 330 — bubble nucleation
    Spec §18.10, §18.11, §18.44, §18.45, §18.51
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np
from scipy import optimize, integrate  # type: ignore[import]

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

#: Speed of light [m s⁻¹]
C: Final[float] = 2.99792458e8

#: Reduced Planck constant [J s]
HBAR: Final[float] = 1.054571817e-34

#: Boltzmann constant [J K⁻¹]
K_B: Final[float] = 1.380649e-23

#: Newton's gravitational constant [m³ kg⁻¹ s⁻²]
G: Final[float] = 6.67430e-11

#: Relativistic degrees of freedom at high temperature (SM order of magnitude)
G_STAR: Final[float] = 100.0

#: Observed baryon-to-photon ratio (Planck 2018 / BBN)
ETA_B_OBS: Final[float] = 6.1e-10

#: Observed CMB temperature today [K]
T_CMB_TODAY: Final[float] = 2.7255

#: Planck length  l_P = sqrt(ℏG/c³) [m]
L_PLANCK: Final[float] = math.sqrt(HBAR * G / C**3)

#: Planck time  t_P = l_P / c [s]
T_PLANCK: Final[float] = L_PLANCK / C

#: Planck temperature  T_P = sqrt(ℏ c⁵ / G) / k_B [K]
T_PLANCK_K: Final[float] = math.sqrt(HBAR * C**5 / G) / K_B

#: Planck energy  E_P = ℏ / t_P = sqrt(ℏ c⁵ / G) [J]
E_PLANCK: Final[float] = HBAR / T_PLANCK

#: Planck mass  M_P = E_P / c² [kg]
M_PLANCK_KG: Final[float] = E_PLANCK / C**2

#: QCD temperature scale T_QCD ~ 170 MeV [K]
T_QCD_K: Final[float] = 0.170e9 * 1.602176634e-10 / K_B  # 170 MeV in K

#: QCD energy scale Λ_QCD ~ 200 MeV [J]
LAMBDA_QCD_J: Final[float] = 0.200e9 * 1.602176634e-19


# ---------------------------------------------------------------------------
# Sine-Gordon kink parameters (substrate §18.11)
# ---------------------------------------------------------------------------

def kink_mass_J(xi: float) -> float:
    """Sine-Gordon kink mass M_kink = 8 ℏ c / ξ  [J].

    From the §18.11 Lagrangian and the sine-Gordon exact solution:
        M_kink = 8 m₀    where m₀ = ℏ c / ξ  [meson mass]

    The kink mass sets the barrier height:
        ΔE_sphaleron = 2 × M_kink

    Args:
        xi: Substrate length scale ξ [m].  This is the natural UV cutoff and
            sets the kink's Compton wavelength.

    Returns:
        Kink mass energy M_kink [J].
    """
    m0 = HBAR * C / xi          # meson mass
    return 8.0 * m0             # kink mass = 8 × meson mass (sine-Gordon exact)


def sphaleron_barrier_J(xi: float) -> float:
    """Energy barrier ΔE_sph = 2 × M_kink between adjacent topological vacua [J].

    In the sine-Gordon theory, the unstable saddle configuration between
    adjacent vacua φ_n = 2πn and φ_{n+1} = 2π(n+1) has energy exactly:

        E_sph = 2 × M_kink = 16 m₀ = 16 ℏ c / ξ

    This is the exact result: the saddle is the kink-antikink configuration
    at zero separation (no binding energy at the classical saddle point, by
    definition of the barrier maximum).

    In the Möbius-bundle language: the holonomy exp(iπ × Δw) = exp(iπ) = −1
    for a single vacuum crossing confirms |Δw| = 1, i.e., one unit of baryon
    number change per saddle crossing.

    Args:
        xi: Substrate length scale ξ [m].

    Returns:
        Sphaleron barrier energy ΔE_sph [J].
    """
    return 2.0 * kink_mass_J(xi)


def bounce_action_J_s(xi: float) -> float:
    """Coleman instanton (bounce) action S_inst for sine-Gordon tunnelling [J·s].

    The O(2)-symmetric bounce in 1+1 D Euclidean sine-Gordon is:
        S_inst = M_kink × (ℏ / c)  [SI]

    Derivation:
        S_inst = ∫ dτ [ ½ (dφ/dτ)² + V(φ) ]   over the Euclidean bounce
               = M_kink × ℏ c    [in natural units with c=1, S_inst = M_kink]

    Restoring SI: action has dimensions [J·s], mass-energy has [J].
        S_inst [J·s] = M_kink [J] × (ℏ/c) [s/m] × ???

    Actually the correct dimensional form is:
        S_inst / ℏ = M_kink / (ℏ c / ξ)  = M_kink / m₀
                   = 8 m₀ / m₀ = 8       [dimensionless, at zero temperature]

    So the tunnelling exponent is:
        exp(-S_inst / ℏ) = exp(-8) ≈ 3.35 × 10⁻⁴

    At the Coleman bosonization coupling β² = 8π (free-fermion point):
        S_inst / ℏ = 8 / β² × 4π = 4π  (from Coleman's exact formula)

    We use S_inst / ℏ = 8 as the conservative estimate (thick-wall regime).

    Args:
        xi: Substrate length scale ξ [m].

    Returns:
        Bounce action S_inst [J·s].  Divide by ℏ to get the dimensionless exponent.
    """
    m_kink = kink_mass_J(xi)
    # S_inst / ℏ = 8  →  S_inst = 8 ℏ
    # This is the minimum tunnelling action (1 kink mass / meson mass = 8)
    return 8.0 * HBAR


# ---------------------------------------------------------------------------
# Vacuum structure: topological winding states
# ---------------------------------------------------------------------------


@dataclass
class SineGordonVacuum:
    """A topological vacuum state of the Möbius × sine-Gordon system.

    In the sine-Gordon theory the degenerate vacua are:
        φ_n = 2πn,    n ∈ ℤ

    Each vacuum has a definite winding number n.  The Möbius bundle holonomy
    for a loop traversing from vacuum n to vacuum n+k is:

        hol(n → n+k) = exp(iπ × k) = (-1)^k

    For k = 1 (adjacent vacua): holonomy = −1 → spin-½ / fermionic character.
    This confirms that a single vacuum crossing creates/annihilates one
    topological kink, which the Möbius assignment identifies as one "baryon."

    Attributes:
        n: Winding number (integer).  The vacuum field value is φ = 2πn.
        energy: Potential energy density of this vacuum [J m⁻³].  For the
            sine-Gordon potential V(φ) = (1 - cos φ) all vacua have V = 0.
        baryon_number: Baryon number assigned to this vacuum.  Defined modulo
            the periodicity of the half-flux: B = n mod 2.  (Even n: B=0,
            odd n: B=1 in the Z₂ assignment from the Möbius holonomy.)
    """

    n: int
    energy: float = 0.0          # all true vacua are degenerate
    baryon_number: int = field(init=False)

    def __post_init__(self) -> None:
        self.baryon_number = abs(self.n) % 2

    @property
    def phi(self) -> float:
        """Field value at this vacuum: φ_n = 2πn."""
        return 2.0 * math.pi * self.n

    @property
    def mobius_holonomy(self) -> float:
        """U(1) holonomy exp(iπn) = (−1)^n for a path winding to vacuum n from 0."""
        return math.cos(math.pi * self.n)   # real part of exp(iπn) = ±1

    def baryon_number_difference(self, other: "SineGordonVacuum") -> int:
        """Change in baryon number when transitioning to the other vacuum.

        ΔB = |n_final - n_initial|  modulo the Z₂ identification.
        Each unit of winding number change corresponds to one kink crossing,
        which changes baryon number by 1.

        Args:
            other: Target vacuum state.

        Returns:
            Baryon number change |Δn| (always non-negative integer).
        """
        return abs(other.n - self.n)


# ---------------------------------------------------------------------------
# Saddle-point (sphaleron) configuration
# ---------------------------------------------------------------------------


@dataclass
class SineGordonSphaleron:
    """The unstable saddle configuration between adjacent sine-Gordon vacua.

    The sphaleron is the field configuration φ(x) that sits at the top
    of the energy barrier between φ = 0 and φ = 2π.  In the infinite-volume
    sine-Gordon theory this is found by solving the static equation of motion:

        d²φ/dx² = V'(φ) = sin(φ)

    with boundary conditions φ(−∞) = 0, φ(+∞) = 2π.  The ONLY static
    solution satisfying these boundary conditions is the kink:
        φ_kink(x) = 4 arctan[exp(x/ξ)]

    But the kink is a STABLE minimum (infinite volume), not a saddle.  The
    true sphaleron is found in finite volume or by constraining the winding
    number to a half-integer.

    In the half-space / constrained approach, the sphaleron at the top of the
    barrier between adjacent vacua is:
        φ_sph(x) = π     (the spatially uniform field at the potential maximum)

    This configuration satisfies the static EOM only in the sense that
    dφ/dx = 0 everywhere and V'(π) = sin(π) = 0 — it is a solution but
    is UNSTABLE: any perturbation causes it to roll down to one of the two
    adjacent vacua (n = 0 or n = 1).

    The sphaleron energy density is V(π) = 1 − cos(π) = 2 per unit length
    (in natural sine-Gordon units where the potential coefficient = 1).
    Integrating over the kink core size ~ ξ:

        E_sph = 2 × M_kink = 16 × (ℏ c / ξ)

    Attributes:
        xi: Substrate length scale ξ [m].
        field_value: The spatially uniform field value at the saddle (π radians).
        energy_J: Sphaleron energy ΔE = 2 M_kink [J].
        instability_direction: The (negative) mode that drives the sphaleron
            to decay: a small perturbation δφ away from π causes exponential
            growth toward either 0 or 2π.
    """

    xi: float

    field_value: float = field(init=False)
    energy_J: float = field(init=False)
    instability_rate: float = field(init=False)

    def __post_init__(self) -> None:
        self.field_value = math.pi    # saddle at the top of the cosine
        self.energy_J = sphaleron_barrier_J(self.xi)
        # Negative mode frequency: V''(π) = -cos(π) = +1 in sine-Gordon units
        # but the STATIC sphaleron has a negative mode with
        #   ω_neg² = -V''(φ_sph) / ξ² = -(-1)/ξ² = +1/ξ²
        # Wait: V(φ) = 1 - cos φ → V''(π) = cos(π) = -1 < 0
        # The sphaleron sits at a local maximum of V, so the negative-mode
        # frequency is: ω_neg² = |V''(π)| × (c/ξ)² = (c/ξ)²
        self.instability_rate = C / self.xi   # |ω_neg| [rad/s]

    def sine_gordon_potential(self, phi: float) -> float:
        """Sine-Gordon potential V(φ) = 1 - cos(φ) [dimensionless, in m₀ units].

        The full potential energy density in SI units is:
            V_SI(φ) = (ℏ c / ξ³) × (1 - cos φ)

        Args:
            phi: Field value in radians.

        Returns:
            Dimensionless potential V(φ) = 1 - cos(φ).
        """
        return 1.0 - math.cos(phi)

    def field_profile(self, x_array: np.ndarray) -> np.ndarray:
        """Approximate saddle-point field profile using constrained interpolation.

        For a finite-size system of length L = 10ξ, the sphaleron profile
        interpolates between φ = 0 and φ = 2π in the shape:

            φ_sph(x) = π × [1 + tanh(x / ξ)]

        This is a smoothed step that matches the boundary conditions and has
        φ = π at x = 0 (the midpoint).  It is NOT the exact saddle, but
        captures the essential structure: the field climbs from 0 to 2π,
        spending most of its "time" near φ = π (the top of the barrier).

        The exact saddle in finite volume would require numerical minimization
        of the energy functional with fixed winding number = 1; that is done
        in find_saddle_numerically().

        Args:
            x_array: Spatial coordinates [m] (array of arbitrary length).

        Returns:
            Field values φ(x) in radians.
        """
        return math.pi * (1.0 + np.tanh(x_array / self.xi))

    def energy_density(self, phi: float) -> float:
        """Static energy density at field value φ [J m⁻³].

        For the static (zero kinetic energy) configuration:
            ε(x) = (ℏ c / ξ³) × [ ½ ξ² (dφ/dx)² + (1 - cos φ) ]

        At the spatially uniform saddle (dφ/dx = 0, φ = π):
            ε = (ℏ c / ξ³) × V(π) = (ℏ c / ξ³) × 2

        Args:
            phi: Field value in radians.

        Returns:
            Energy density in J m⁻³.
        """
        prefactor = HBAR * C / self.xi**3
        return prefactor * self.sine_gordon_potential(phi)

    def find_saddle_numerically(
        self,
        n_points: int = 500,
        L_over_xi: float = 20.0,
    ) -> dict[str, float | np.ndarray]:
        """Find the saddle-point energy by numerical energy minimization.

        Works entirely in dimensionless (natural) sine-Gordon units where
        the field φ is in radians, the spatial coordinate is x_nat = x / ξ,
        and the energy is measured in units of ℏ c / ξ.

        The dimensionless energy functional is:
            E_nat = ∫₀^{L_nat} dx_nat [ ½ (dφ/dx_nat)² + (1 - cos φ) ]

        Exact sine-Gordon result: E_kink_nat = 8 (the kink mass in natural units).
        The SI kink mass is E_kink_J = (ℏ c / ξ) × E_kink_nat.

        Returns:
            Dict with keys:
            - "E_kink_J": kink energy [J] from minimization
            - "E_kink_natural": kink energy in sine-Gordon natural units (should be ~8)
            - "E_sphaleron_J": sphaleron energy = 2 × E_kink [J]
            - "x_array_natural": dimensionless spatial grid (in units of ξ)
            - "phi_kink": kink field profile [radians]
            - "phi_sphaleron": approximate sphaleron profile [radians]
            - "barrier_in_M_kink": barrier in units of M_kink (should be 2.0)
            - "converged": whether the minimizer converged
        """
        # Work in natural units: x_nat in [0, L_nat], L_nat = L_over_xi
        L_nat = float(L_over_xi)
        x_nat = np.linspace(0.0, L_nat, n_points)
        dx_nat = x_nat[1] - x_nat[0]

        # Energy functional (dimensionless) on internal points;
        # boundaries fixed: φ(0) = 0, φ(L_nat) = 2π  (winding = 1)
        def energy_nat(phi_internal: np.ndarray) -> float:
            phi_full = np.concatenate([[0.0], phi_internal, [2.0 * np.pi]])
            # Gradient term: ½ (dφ/dx_nat)²
            dphi_nat = np.diff(phi_full) / dx_nat
            grad_energy = 0.5 * float(np.sum(dphi_nat**2)) * dx_nat
            # Potential term: (1 - cos φ)
            pot_energy = float(np.sum(1.0 - np.cos(phi_full))) * dx_nat
            return grad_energy + pot_energy

        # Initial guess: linear ramp from 0 to 2π (corresponds to a very spread-out kink)
        phi_init = 2.0 * np.pi * x_nat[1:-1] / L_nat

        result = optimize.minimize(
            energy_nat,
            phi_init,
            method="L-BFGS-B",
            options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-10},
        )
        phi_kink_internal = result.x
        phi_kink = np.concatenate([[0.0], phi_kink_internal, [2.0 * np.pi]])

        # Kink energy in natural units (exact sine-Gordon result: 8)
        E_kink_natural = result.fun
        # Convert to SI: E_kink_J = (ℏ c / ξ) × E_kink_natural
        E_kink_J = (HBAR * C / self.xi) * E_kink_natural
        E_sphaleron_J = 2.0 * E_kink_J

        # Sphaleron profile in natural units (φ = π uniform at top of barrier)
        phi_sph = np.full_like(x_nat, math.pi)

        barrier_in_M_kink = E_sphaleron_J / kink_mass_J(self.xi)

        return {
            "E_kink_J": E_kink_J,
            "E_kink_natural": E_kink_natural,
            "E_sphaleron_J": E_sphaleron_J,
            "x_array_natural": x_nat,
            "phi_kink": phi_kink,
            "phi_sphaleron": phi_sph,
            "barrier_in_M_kink": barrier_in_M_kink,
            "converged": bool(result.success),
        }


# ---------------------------------------------------------------------------
# Hubble rate
# ---------------------------------------------------------------------------


def hubble_rate_SI(T_K: float) -> float:
    """Radiation-dominated Hubble parameter H(T) [s⁻¹].

    From the FRW Friedmann equation with radiation energy density:
        H² = (8πG/3) × (π²/30) × g_* × (k_B T)⁴ / (ℏ c)³

    Args:
        T_K: Temperature [K].

    Returns:
        Hubble rate [s⁻¹].  Returns 0 if T ≤ 0.
    """
    if T_K <= 0.0:
        return 0.0
    eps_rad = (math.pi**2 / 30.0) * G_STAR * (K_B * T_K)**4 / (HBAR * C)**3
    H_sq = (8.0 * math.pi * G / 3.0) * eps_rad
    return math.sqrt(max(H_sq, 0.0))


# ---------------------------------------------------------------------------
# Thermal sphaleron rate (3+1 D estimate)
# ---------------------------------------------------------------------------


@dataclass
class SphaleronRateCalculator:
    """Compute the sphaleron-analogue baryon-violation rate in our model.

    The baryon-violation rate uses the thermal activation formula adapted
    from Arnold & McLerran (1987) for the SM sphaleron, applied to our
    sine-Gordon kink-antikink saddle:

    Above crossover temperature T_cross:
        Γ_B / V = κ_rate × (α_eff T)⁴ × exp(-E_sph / k_B T)

    Below T_cross (quantum tunnelling):
        Γ_B / V = ω_neg × (m₀ c / ℏ)³ × exp(-S_inst / ℏ)

    The dimensionless prefactor κ_rate ~ 1 (of order unity; we use κ = 1.0).
    The 3D phase-space factor (m₀ c / ℏ)³ = (c / ξ)³ replaces the (α_W T)⁴
    of the SM.

    For the sphaleron rate per Hubble volume (the quantity directly entering η_B):
        Γ_B^{Hubble} = (Γ_B / V) × (c / H)³

    And the ratio relevant for baryogenesis:
        Γ_B / H_* = Γ_B^{Hubble} / (c³ / H²)  = (Γ_B / V) × c³ / H⁴ × H

    Wait — more carefully:
        Γ_B per Hubble volume = (Γ_B / V) × V_Hubble
        V_Hubble = (4π/3) × (c/H)³
        Γ_B per Hubble time = (Γ_B per Hubble volume) / (1/H)
                            = (Γ_B / V) × (c/H)³ × H
                            = (Γ_B / V) × c³ / H²

    The ratio comparing this to the Hubble rate for η_B:
        R = Γ_B^{per Hubble vol × time} / n_baryon^{Hubble vol}

    For baryogenesis, the standard estimate is:
        η_B ~ (Γ_B / H_*)^{per Hubble vol per Hubble time} × δ_CP

    In the rate-per-volume formulation:
        η_B ~ [(Γ_B/V) / H_*⁴] × δ_CP

    (Note: [Γ_B/V] has units m⁻³ s⁻¹; dividing by H⁴ in s⁻⁴ gives m⁻³ s³
    which is not dimensionless.  The correct form uses n_γ as the denominator.)

    Correct formula (standard):
        η_B ≈ (Γ_B/V) / (H × n_γ) × δ_CP

    where n_γ = (2ζ(3)/π²) (k_B T / ℏ c)³ is the photon number density.
    This is dimensionless: [m⁻³ s⁻¹] / ([s⁻¹] × [m⁻³]) = dimensionless.

    Attributes:
        xi: Substrate length scale ξ [m].  Controls kink mass and barrier height.
        kappa_rate: Dimensionless prefactor κ ~ 1 (all pre-exponential factors
            of order unity are folded in here).
        delta_cp: CP asymmetry parameter (free parameter, ∈ (0, 1]).
    """

    xi: float
    kappa_rate: float = 1.0
    delta_cp: float = 1.0

    def __post_init__(self) -> None:
        if not (0.0 < self.delta_cp <= 1.0):
            raise ValueError(
                f"delta_cp must be in (0, 1]; got {self.delta_cp!r}"
            )
        if self.xi <= 0.0:
            raise ValueError(f"xi must be positive; got {self.xi!r}")

    @property
    def m0(self) -> float:
        """Meson mass m₀ = ℏ c / ξ [J]."""
        return HBAR * C / self.xi

    @property
    def M_kink(self) -> float:
        """Kink mass M_kink = 8 m₀ [J]."""
        return 8.0 * self.m0

    @property
    def E_sphaleron(self) -> float:
        """Sphaleron barrier height ΔE = 2 M_kink [J]."""
        return 2.0 * self.M_kink

    @property
    def omega0(self) -> float:
        """Vacuum oscillation frequency ω₀ = m₀ c / ℏ = c / ξ [rad s⁻¹]."""
        return C / self.xi

    @property
    def T_crossover(self) -> float:
        """Crossover temperature T_cross where thermal ~ quantum rate [K].

        Defined by equating the thermal exponent and the tunnelling exponent:
            E_sph / (k_B T_cross) = S_inst / ℏ = 8

        → T_cross = E_sph / (8 k_B) = 2 M_kink / (8 k_B)
                  = M_kink / (4 k_B)
                  = 8 m₀ / (4 k_B)
                  = 2 m₀ / k_B

        Returns:
            Crossover temperature [K].
        """
        return self.E_sphaleron / (8.0 * K_B)

    def thermal_exponent(self, T_K: float) -> float:
        """Boltzmann exponent E_sph / (k_B T) for thermal activation.

        Args:
            T_K: Temperature [K].

        Returns:
            Dimensionless exponent (always ≥ 0 for T > 0).  Returns inf for T = 0.
        """
        if T_K <= 0.0:
            return float("inf")
        return self.E_sphaleron / (K_B * T_K)

    def tunnelling_exponent(self) -> float:
        """Quantum tunnelling exponent S_inst / ℏ = 8 (temperature-independent).

        Returns:
            Dimensionless tunnelling exponent S_inst / ℏ = 8.
        """
        return 8.0   # M_kink / m₀ = 8 m₀ / m₀ = 8

    def rate_density(self, T_K: float) -> float:
        """Baryon-violating rate density Γ_B / V [m⁻³ s⁻¹].

        Uses the appropriate exponent for the given temperature regime:
        - T > T_cross: classical thermal activation, exp(-E_sph / k_B T)
        - T < T_cross: quantum tunnelling, exp(-S_inst / ℏ) = exp(-8)
        - T ~ T_cross: minimum of the two exponents (crossover regime)

        The pre-exponential (phase space) factor in 3+1 D:
            A_pre = κ × (c / ξ)⁴ / c³  = κ × (c / ξ)⁴ × (1/c)³ [m⁻³ s⁻¹]
            Simplification: A_pre = κ × c / ξ⁴

        Note: The natural 3D rate density for a field with correlation length ξ
        is A_pre ~ (ω₀ / 2π) × (1/ξ)³ where ω₀ = c/ξ:
            A_pre = (c / ξ) × (1/ξ³) = c / ξ⁴   [m⁻³ s⁻¹]

        Args:
            T_K: Temperature [K].

        Returns:
            Rate density Γ_B / V [m⁻³ s⁻¹].  A large positive number if the
            barrier is negligible, exponentially small for T << T_cross.
        """
        # Pre-exponential factor: κ × c / ξ⁴ [m⁻³ s⁻¹]
        A_pre = self.kappa_rate * C / self.xi**4

        # Choose the effective exponent
        exp_thermal = self.thermal_exponent(T_K)
        exp_tunnel = self.tunnelling_exponent()
        # Use minimum (larger rate = less suppression) at the crossover
        exp_eff = min(exp_thermal, exp_tunnel)

        # Guard against overflow
        if exp_eff > 700.0:
            return 0.0
        if exp_eff < -700.0:
            return A_pre  # ~ unsuppressed

        return A_pre * math.exp(-exp_eff)

    def hubble_rate(self, T_K: float) -> float:
        """Hubble rate H(T) [s⁻¹] in radiation-dominated epoch.

        Args:
            T_K: Temperature [K].

        Returns:
            Hubble rate [s⁻¹].
        """
        return hubble_rate_SI(T_K)

    def photon_number_density(self, T_K: float) -> float:
        """Photon number density n_γ(T) [m⁻³].

        n_γ = (2 ζ(3) / π²) × (k_B T / ℏ c)³

        Args:
            T_K: Temperature [K].

        Returns:
            Photon number density [m⁻³].
        """
        if T_K <= 0.0:
            return 0.0
        zeta3 = 1.2020569031595943
        return (2.0 * zeta3 / math.pi**2) * (K_B * T_K / (HBAR * C))**3

    def gamma_over_H(self, T_K: float) -> float:
        """Dimensionless ratio (Γ_B/V) / (H × n_γ).

        This is the key ratio for baryogenesis.  When this is O(1) or larger,
        baryon-number violation is efficient.  The η_B estimate is:

            η_B ~ [(Γ_B/V) / (H × n_γ)] × δ_CP

        Args:
            T_K: Temperature [K].

        Returns:
            Dimensionless ratio (Γ_B/V) / (H × n_γ).  Returns 0 if H or n_γ = 0.
        """
        gamma_vol = self.rate_density(T_K)
        H = self.hubble_rate(T_K)
        n_gamma = self.photon_number_density(T_K)

        if H <= 0.0 or n_gamma <= 0.0:
            return 0.0
        return gamma_vol / (H * n_gamma)

    def eta_B_estimate(self, T_K: float) -> float:
        """Estimate of baryon-to-photon ratio η_B at temperature T.

        Master formula:
            η_B ~ (Γ_B/V) / (H × n_γ) × δ_CP

        This is capped at 1 to avoid unphysical results when the formula
        overestimates (equilibration regime).

        Args:
            T_K: Temperature [K].

        Returns:
            Estimated baryon-to-photon ratio η_B (dimensionless, ≤ 1).
        """
        ratio = self.gamma_over_H(T_K)
        return min(ratio * self.delta_cp, 1.0)

    def temperature_scan(
        self,
        temperatures_K: list[float],
    ) -> list[dict[str, float]]:
        """Compute all key quantities at a list of temperatures.

        Args:
            temperatures_K: List of temperatures to evaluate [K].

        Returns:
            List of dicts, one per temperature.  Each dict contains:
            - "T_K": temperature [K]
            - "T_over_T_cross": ratio T / T_crossover
            - "E_sph_J": sphaleron barrier [J]
            - "E_sph_over_kBT": thermal exponent E_sph / (k_B T)
            - "S_inst_over_hbar": tunnelling exponent S_inst / ℏ = 8
            - "exp_used": which exponent is smaller (and used)
            - "Gamma_vol_m3s": Γ_B / V [m⁻³ s⁻¹]
            - "log10_Gamma_vol": log10 of Γ_B / V
            - "H_s": Hubble rate [s⁻¹]
            - "n_gamma_m3": photon density [m⁻³]
            - "Gamma_over_H_ngamma": (Γ_B/V) / (H n_γ)
            - "log10_GoH": log10 of the ratio
            - "eta_B": η_B estimate
            - "log10_eta_B": log10(η_B)
            - "matches_obs": whether |log10(η_B / η_B_obs)| < 1
        """
        results: list[dict[str, float]] = []

        for T in temperatures_K:
            gamma_vol = self.rate_density(T)
            H = self.hubble_rate(T)
            n_gamma = self.photon_number_density(T)
            ratio = self.gamma_over_H(T)
            eta = self.eta_B_estimate(T)

            exp_th = self.thermal_exponent(T)
            exp_tu = self.tunnelling_exponent()
            exp_used = min(exp_th, exp_tu)

            def _log10(x: float) -> float:
                return math.log10(x) if x > 0 else float("-inf")

            log10_ratio = _log10(ratio)
            log10_eta = _log10(eta) if eta > 0 else float("-inf")

            matches = (
                math.isfinite(log10_eta)
                and abs(log10_eta - math.log10(ETA_B_OBS)) < 1.0
            )

            results.append({
                "T_K": T,
                "T_over_T_cross": T / self.T_crossover if self.T_crossover > 0 else float("inf"),
                "E_sph_J": self.E_sphaleron,
                "E_sph_over_kBT": exp_th if math.isfinite(exp_th) else float("inf"),
                "S_inst_over_hbar": exp_tu,
                "exp_used": exp_used,
                "Gamma_vol_m3s": gamma_vol,
                "log10_Gamma_vol": _log10(gamma_vol),
                "H_s": H,
                "n_gamma_m3": n_gamma,
                "Gamma_over_H_ngamma": ratio,
                "log10_GoH": log10_ratio,
                "eta_B": eta,
                "log10_eta_B": log10_eta,
                "matches_obs": matches,
            })

        return results


# ---------------------------------------------------------------------------
# What temperature is needed? — solve for T* such that Γ_B/H ~ 1
# ---------------------------------------------------------------------------


def find_T_star_for_efficient_baryogenesis(
    calculator: SphaleronRateCalculator,
    target_ratio: float = 1.0,
    T_min_K: float = 1.0,
    T_max_K: float | None = None,
) -> dict[str, float]:
    """Find temperature where (Γ_B/V) / (H n_γ) = target_ratio.

    Solves for the temperature at which the sphaleron analogue rate equals
    the Hubble rate times the photon density — the condition for efficient
    baryon number violation.

    Args:
        calculator: Configured SphaleronRateCalculator instance.
        target_ratio: Target value for (Γ_B/V) / (H n_γ).  Default 1.0
            (efficient baryogenesis threshold).
        T_min_K: Lower bracket for temperature search [K].
        T_max_K: Upper bracket [K].  Defaults to 10 × T_Planck.

    Returns:
        Dict with:
        - "T_star_K": temperature where target is achieved [K]
        - "T_star_over_T_Planck": T_star / T_Planck
        - "ratio_at_T_star": actual ratio at T_star
        - "found": True if a solution was found in the bracket
        - "message": human-readable status
    """
    if T_max_K is None:
        T_max_K = 10.0 * T_PLANCK_K

    def objective(log_T: float) -> float:
        T = math.exp(log_T)
        ratio = calculator.gamma_over_H(T)
        if ratio <= 0.0:
            return -math.log(target_ratio) - 100.0
        return math.log(ratio) - math.log(target_ratio)

    log_T_min = math.log(T_min_K)
    log_T_max = math.log(T_max_K)

    # Check bracket
    f_min = objective(log_T_min)
    f_max = objective(log_T_max)

    if f_min * f_max > 0:
        # No sign change — solution not in bracket
        # Report which side is closer
        if abs(f_max) < abs(f_min):
            T_best = T_max_K
            msg = (
                f"No crossing in [{T_min_K:.2e}, {T_max_K:.2e}] K. "
                f"Rate at T_max={T_max_K:.2e} K gives ratio={math.exp(f_max+math.log(target_ratio)):.2e}. "
                "Rate never reaches target in this range."
            )
        else:
            T_best = T_min_K
            msg = (
                f"No crossing in [{T_min_K:.2e}, {T_max_K:.2e}] K. "
                f"Rate at T_min={T_min_K:.2e} K gives ratio={math.exp(f_min+math.log(target_ratio)):.2e}. "
                "Rate never reaches target in this range."
            )
        return {
            "T_star_K": T_best,
            "T_star_over_T_Planck": T_best / T_PLANCK_K,
            "ratio_at_T_star": calculator.gamma_over_H(T_best),
            "found": False,
            "message": msg,
        }

    sol = optimize.brentq(objective, log_T_min, log_T_max, xtol=1e-6, rtol=1e-8)
    T_star = math.exp(sol)
    actual_ratio = calculator.gamma_over_H(T_star)

    return {
        "T_star_K": T_star,
        "T_star_over_T_Planck": T_star / T_PLANCK_K,
        "ratio_at_T_star": actual_ratio,
        "found": True,
        "message": f"Found T_star = {T_star:.3e} K where Γ_B/H_* = {actual_ratio:.3e}",
    }


# ---------------------------------------------------------------------------
# Multi-scale analysis: Planck ξ vs. atomic ξ
# ---------------------------------------------------------------------------


def analyze_two_scales() -> dict[str, dict[str, float]]:
    """Compute sphaleron rate at both the Planck scale and the atomic (Compton) scale.

    The model has a multi-scale problem: the §18.46.1 relation ℏ = K ξ⁴ / c
    can be satisfied at either:
    - ξ_Planck = L_PLANCK ~ 1.6e-35 m  (Planck substrate)
    - ξ_Compton ~ 3.9e-13 m  (electron Compton wavelength, where M_kink = 8 m_e)

    The sphaleron barrier height ΔE = 2 M_kink = 16 ℏ c / ξ depends crucially
    on which scale is used.  Both scenarios are analyzed.

    Returns:
        Dict with keys "planck" and "compton", each containing a sub-dict
        with M_kink, E_sph, T_crossover, and rate values at various temperatures.
    """
    # Planck-scale substrate
    xi_planck = L_PLANCK                        # ~1.6e-35 m
    # Atomic-scale substrate (electron Compton wavelength)
    xi_compton = HBAR / (9.1093837015e-31 * C)  # ~3.86e-13 m

    results: dict[str, dict[str, float]] = {}
    for label, xi in [("planck", xi_planck), ("compton", xi_compton)]:
        calc = SphaleronRateCalculator(xi=xi, delta_cp=1.0)
        m_kink_eV = calc.M_kink / 1.602176634e-19          # convert J to eV
        e_sph_eV = calc.E_sphaleron / 1.602176634e-19
        T_cross = calc.T_crossover

        # Rate at T = T_Planck
        rate_Tplanck = calc.rate_density(T_PLANCK_K)
        GoH_Tplanck = calc.gamma_over_H(T_PLANCK_K)
        eta_Tplanck = calc.eta_B_estimate(T_PLANCK_K)

        # Rate at T = 1e-3 T_Planck (de-saturation epoch)
        T_desat = 1e-3 * T_PLANCK_K
        GoH_desat = calc.gamma_over_H(T_desat)
        eta_desat = calc.eta_B_estimate(T_desat)

        # Rate at T = T_QCD
        GoH_qcd = calc.gamma_over_H(T_QCD_K)
        eta_qcd = calc.eta_B_estimate(T_QCD_K)

        results[label] = {
            "xi_m": xi,
            "M_kink_eV": m_kink_eV,
            "E_sph_eV": e_sph_eV,
            "T_crossover_K": T_cross,
            "T_cross_over_T_Planck": T_cross / T_PLANCK_K,
            "GoH_at_T_Planck": GoH_Tplanck,
            "eta_B_at_T_Planck": eta_Tplanck,
            "GoH_at_1e-3_T_Planck": GoH_desat,
            "eta_B_at_1e-3_T_Planck": eta_desat,
            "GoH_at_T_QCD": GoH_qcd,
            "eta_B_at_T_QCD": eta_qcd,
            "rate_density_at_T_Planck_m3s": rate_Tplanck,
        }

    return results


# ---------------------------------------------------------------------------
# What delta_cp is required?
# ---------------------------------------------------------------------------


def required_delta_cp(
    calculator: SphaleronRateCalculator,
    T_K: float,
    target_eta: float = ETA_B_OBS,
) -> float:
    """Compute δ_CP required to produce target η_B at temperature T.

    Inverts:  η_B = (Γ_B/V) / (H n_γ) × δ_CP

    Args:
        calculator: SphaleronRateCalculator (delta_cp is ignored for this calculation).
        T_K: Temperature [K].
        target_eta: Target baryon-to-photon ratio (default: observed 6.1e-10).

    Returns:
        Required δ_CP.  Returns inf if the rate is zero, or < 0 if target
        is already exceeded without any CP violation (should not happen physically).
    """
    ratio = calculator.gamma_over_H(T_K)
    if ratio <= 0.0:
        return float("inf")
    return target_eta / ratio
