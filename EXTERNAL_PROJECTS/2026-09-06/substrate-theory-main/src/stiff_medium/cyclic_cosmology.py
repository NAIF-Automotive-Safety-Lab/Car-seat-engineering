"""Cyclic cosmology — FRW dynamics and cyclic universe (§§18.40-18.43).

The substrate-mechanical cosmology of the stiff-medium theory:

- Universe begins saturated (σ = 0.5), the Big Bang state identical to a
  black hole interior (§18.42).
- De-saturation drives structure formation through the FRW Friedmann equations.
- The universe asymptotically reaches one of two end-states:
    End-state A: all matter accreted into one universe-spanning saturated
                 region, which IS the Big Bang state → cycle restarts.
    End-state B: Hawking evaporation returns all mass to radiation; universe
                 enters heat death at baseline strain σ₀; quantum nucleation
                 of a new saturated region triggers the next cycle.
- The substrate (K, ρ, ξ, Möbius half-flux) persists eternally; only the
  matter pattern within it resets between cycles.

See spec §§18.40-18.43 for the full theoretical development.
See `saturation.py` for the σ_max = ½ elastic-limit mechanics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Literal

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

c: float = 2.998e8       # speed of light, m s⁻¹
hbar: float = 1.055e-34  # reduced Planck constant, J s
G: float = 6.674e-11     # Newton's gravitational constant, m³ kg⁻¹ s⁻²
k_B: float = 1.381e-23   # Boltzmann constant, J K⁻¹

# Derived helpers
_YEAR_S: float = 3.15576e7   # seconds per Julian year
_MPC_M: float = 3.0857e22    # metres per megaparsec

# ---------------------------------------------------------------------------
# Universe phase enum
# ---------------------------------------------------------------------------


class UniversePhase(Enum):
    """Qualitative epoch of the current cycle.

    Thresholds follow the standard density-parameter dominance scheme
    applied to the Planck-2018 cosmology:

    - ``SATURATED``       : σ ≥ σ_max (= 0.5) — the Big Bang / inflationary era.
    - ``MATTER_DOMINATED``  : Ω_m(a) > Ω_Λ(a) and a < a_eq(m,r) — standard structure.
    - ``LAMBDA_DOMINATED``  : Ω_Λ(a) > Ω_m(a) — dark-energy-driven acceleration.
    - ``HEAT_DEATH``        : energy density below the baseline σ₀ threshold.
    """

    SATURATED = auto()
    MATTER_DOMINATED = auto()
    LAMBDA_DOMINATED = auto()
    HEAT_DEATH = auto()


# ---------------------------------------------------------------------------
# Module-level cosmological functions
# ---------------------------------------------------------------------------


def hubble_radius(H: float) -> float:
    """Return the Hubble radius c / H.

    The Hubble radius is the distance at which the recession velocity equals
    the speed of light. In the substrate-mechanical picture it sets the largest
    scale over which the medium's elastic response can be coherent at epoch H.

    Args:
        H: Hubble parameter in s⁻¹.

    Returns:
        Hubble radius in metres.  Returns +inf if H == 0.
    """
    if H == 0.0:
        return math.inf
    return c / H


def critical_density(H: float, G_newton: float = G) -> float:
    """Return the critical density 3H² / (8πG).

    The critical density is the energy density that makes the universe exactly
    spatially flat. In our model it corresponds to the substrate strain energy
    density required for marginal FRW closure.

    Args:
        H:        Hubble parameter in s⁻¹.
        G_newton: Newton's constant in m³ kg⁻¹ s⁻².  Defaults to SI value.

    Returns:
        Critical density in kg m⁻³.
    """
    return 3.0 * H**2 / (8.0 * math.pi * G_newton)


def cycle_period_estimate(endpoint: Literal["A", "B"] = "B") -> float:
    """Order-of-magnitude cycle period in years.

    End-state A (saturation): matter merges faster than Hawking evaporation.
    Upper bound ~ 10⁴⁰ years (dynamical collapse of last BHs).

    End-state B (heat-death): dominated by Hawking evaporation of the most
    massive BH in the observable universe (≈ 10⁵³ kg horizon mass) via
    t_evap = 5120π G² M³ / (ℏ c⁴), giving ≈ 10¹⁰⁰ years.

    Args:
        endpoint: ``"A"`` for saturation path, ``"B"`` for heat-death path.

    Returns:
        Estimated cycle duration in years (order-of-magnitude).

    References:
        §18.43 — Cyclic cosmology, specific ratio of cycle timescales.
    """
    if endpoint == "A":
        # Dynamical collapse timescale — BH mergers, gravitational waves
        return 1e40
    # End-state B: Hawking evaporation of universe-spanning BH (M ≈ 10⁵³ kg)
    M_horizon = 1e53  # kg — approximate mass of observable universe
    t_evap_s = 5120.0 * math.pi * G**2 * M_horizon**3 / (hbar * c**4)
    return t_evap_s / _YEAR_S


def hawking_evaporation_time(M_kg: float) -> float:
    """Hawking evaporation time for a Schwarzschild black hole of mass M.

    Formula: t_evap = 5120 π G² M³ / (ℏ c⁴)

    Args:
        M_kg: Black hole mass in kilograms.

    Returns:
        Evaporation time in seconds.
    """
    return 5120.0 * math.pi * G**2 * M_kg**3 / (hbar * c**4)


# ---------------------------------------------------------------------------
# FRWUniverse dataclass
# ---------------------------------------------------------------------------

#: H₀ = 67.4 km s⁻¹ Mpc⁻¹ in SI units (Planck 2018, Table 2 column "TT,TE,EE+lowE+lensing")
_H0_SI: float = 67.4e3 / _MPC_M   # s⁻¹ ≈ 2.184e-18 s⁻¹

#: Baseline vacuum strain (§18.38: σ₀ ~ 5e-62)
_SIGMA_0: float = 5e-62

#: Maximum substrate strain before saturation (§18.39)
_SIGMA_MAX: float = 0.5

#: Scale factor at matter-Λ equality (Planck 2018 parameters)
_A_EQUALITY_M_LAMBDA: float = 0.754   # a when Ω_m(a) = Ω_Λ(a), approx


@dataclass
class FRWUniverse:
    """A Friedmann–Robertson–Walker (FRW) universe on the stiff-medium substrate.

    Represents the state of a single cosmic cycle.  The universe begins at
    saturation (σ = 0.5) and evolves via Friedmann dynamics.  Cycle endpoints
    trigger a restart through ``CyclicCosmology``.

    **Scale factor convention:**  a = 1 corresponds to the current epoch (z = 0,
    "today").  The post-CMB de-saturation era started at a_CMB ≈ 1/1100 ≈ 9.1e-4.
    When initialising a fresh cycle from the Big Bang state, set ``a`` to
    ``a_CMB`` (≈ 9.1e-4) and ``sigma`` to 0.5; the subsequent Euler steps
    then expand the universe towards a = 1 today.

    For convenience, the default ``sigma`` at ``a = 1`` (today) is computed
    from the de-saturation formula σ(a) = σ_max × (a_CMB/a)³, which gives
    σ ≈ 3.8×10⁻¹⁰ — physically the residual substrate strain from today's
    matter density (§18.38, §18.40).

    All energy fractions are evaluated at a = 1 (today), following the
    standard Planck-2018 parametrisation.

    Args:
        sigma:            Substrate strain at current epoch.  For Big Bang
                          initial state: σ = 0.5.  For z = 0 (today): σ ≈ 3.8e-13.
                          Defaults to σ_today computed from the de-saturation
                          formula at a = 1.
        a:                Dimensionless scale factor (a = 1 at z = 0 today;
                          a_CMB ≈ 9.1e-4 at the de-saturation transition).
        omega_matter:     Matter density fraction Ω_m at a = 1.  Default 0.315
                          (Planck 2018).
        omega_lambda:     Dark-energy fraction Ω_Λ at a = 1.  Default 0.685
                          (Planck 2018).
        omega_radiation:  Radiation density fraction Ω_r at a = 1.  Default
                          9.4e-5 (photons + neutrinos, Planck 2018).
        H_0:              Hubble parameter today in s⁻¹.  Default corresponds
                          to 67.4 km s⁻¹ Mpc⁻¹ (Planck 2018).

    References:
        §18.40 — Saturated initial condition and de-saturation.
        §18.42 — Big Bang ≡ black hole interior (σ = ½).
        §18.43 — Two end-states and cycle restart.
    """

    sigma: float = field(default=_SIGMA_MAX * (1.0 / 1100.0) ** 3)
    a: float = 1.0
    omega_matter: float = 0.315
    omega_lambda: float = 0.685
    omega_radiation: float = 9.4e-5
    H_0: float = _H0_SI

    # Internal clock: cumulative proper time elapsed since cycle start (seconds)
    _age_s: float = field(default=0.0, init=False, repr=False)

    # Cache the last computed H for convenience
    _H_last: float = field(default=_H0_SI, init=False, repr=False)

    # ------------------------------------------------------------------
    # Post-init validation
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        if self.a <= 0.0:
            raise ValueError(f"Scale factor must be positive; got a = {self.a}")
        if not (0.0 <= self.sigma <= _SIGMA_MAX + 1e-12):
            raise ValueError(
                f"Substrate strain σ must be in [0, {_SIGMA_MAX}]; got σ = {self.sigma}"
            )
        if self.H_0 <= 0.0:
            raise ValueError(f"H_0 must be positive; got H_0 = {self.H_0}")
        self._H_last = self.friedmann_equation(self.a)

    # ------------------------------------------------------------------
    # Friedmann equation
    # ------------------------------------------------------------------

    def friedmann_equation(self, a: float) -> float:
        """Return the Hubble parameter H(a) via the Friedmann equation.

        The first Friedmann equation (flat universe, k = 0):

            H²(a) = (8πG/3) × Σ εᵢ(a)

        where the energy densities scale as:
            - Matter:     εₘ(a) ∝ a⁻³  →  ρ_m,0 × (a₀/a)³
            - Radiation:  εᵣ(a) ∝ a⁻⁴  →  ρ_r,0 × (a₀/a)⁴
            - Λ:          ε_Λ(a) = const (w = -1 from substrate strain, §18.38)

        In terms of the dimensionless density parameters at a = 1:

            H²(a) = H₀² × [Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_Λ]

        This is the standard ΛCDM Friedmann equation, recovered here from the
        substrate-mechanical framework (§18.40: saturated-medium equation of
        state drives the expansion; §18.38: Ω_Λ arises from vacuum substrate
        strain σ₀).

        Args:
            a: Scale factor at which to evaluate H.

        Returns:
            Hubble parameter H(a) in s⁻¹.  Always ≥ 0.

        References:
            §18.40 — H equation from saturation physics.
            §18.38 — Λ from vacuum strain σ₀.
        """
        if a <= 0.0:
            raise ValueError(f"Scale factor must be positive; got a = {a}")
        H_squared = self.H_0**2 * (
            self.omega_matter * a**-3
            + self.omega_radiation * a**-4
            + self.omega_lambda
        )
        return math.sqrt(max(H_squared, 0.0))

    # ------------------------------------------------------------------
    # Time evolution
    # ------------------------------------------------------------------

    def evolve(self, dt: float) -> None:
        """Step the universe forward by proper time dt.

        Uses a first-order Euler step on the FRW scale factor:

            da/dt = a × H(a)

        and updates the substrate strain via the de-saturation formula:

            σ(a) = σ_max × (a_CMB / a)³

        where a_CMB = 1/1100 ≈ 9.1×10⁻⁴ is the scale factor at the CMB
        de-saturation phase transition (§18.44), and the exponent 3 reflects
        the matter-density scaling (strain energy tracks Ω_m).  At a = a_CMB
        the formula recovers σ = σ_max = 0.5 (saturated).  At a = 1 (today)
        it gives σ ≈ 3.8×10⁻¹³, consistent with the sub-cosmological strain
        expected from today's matter density.

        The strain floor is the baseline vacuum σ₀ ~ 5×10⁻⁶² (§18.38);
        the strain ceiling is σ_max = 0.5 (saturation, §18.39).

        Args:
            dt: Proper time increment in seconds.  Must be positive.

        Raises:
            ValueError: If dt ≤ 0.

        References:
            §18.39 — σ_max = ½ elastic limit.
            §18.40 — De-saturation as the universe expands.
            §18.44 — CMB as de-saturation phase transition (a_CMB = 1/1100).
        """
        if dt <= 0.0:
            raise ValueError(f"Time step dt must be positive; got dt = {dt}")

        H = self.friedmann_equation(self.a)
        self._H_last = H

        # Euler step on scale factor: a′ = a + a H dt
        da = self.a * H * dt
        self.a = max(self.a + da, 1e-300)   # guard against underflow

        # Substrate strain: σ(a) = σ_max × (a_CMB / a)³
        # a_CMB = 1/1100 (de-saturation transition, §18.44)
        a_cmb = 1.0 / 1100.0
        sigma_raw = _SIGMA_MAX * (a_cmb / self.a) ** 3
        self.sigma = max(min(sigma_raw, _SIGMA_MAX), _SIGMA_0)

        self._age_s += dt

    # ------------------------------------------------------------------
    # Phase classification
    # ------------------------------------------------------------------

    def current_phase(self) -> UniversePhase:
        """Return the qualitative cosmological phase of the current epoch.

        Classification rule (by dominant energy component at scale factor a):

        - ``SATURATED``        : σ ≥ σ_max (Big Bang / inflationary era, §18.40).
        - ``HEAT_DEATH``       : H(a) / H₀ < 10⁻³⁰ (expansion has essentially frozen).
        - ``LAMBDA_DOMINATED`` : Ω_Λ(a) > Ω_m(a)  (dark-energy era, a > a_eq).
        - ``MATTER_DOMINATED`` : Ω_m(a) ≥ Ω_Λ(a)  (structure-formation era).

        Returns:
            ``UniversePhase`` enum member.

        References:
            §18.40 — Saturated phase and de-saturation.
            §18.43 — Heat death (End-state B).
        """
        if self.sigma >= _SIGMA_MAX - 1e-12:
            return UniversePhase.SATURATED

        H_current = self.friedmann_equation(self.a)
        # Heat death: Hubble rate has dropped to negligible levels
        if H_current < self.H_0 * 1e-30:
            return UniversePhase.HEAT_DEATH

        # Compare matter vs Λ energy densities at current a
        rho_m = self.omega_matter * self.a**-3
        rho_lambda = self.omega_lambda
        if rho_lambda > rho_m:
            return UniversePhase.LAMBDA_DOMINATED
        return UniversePhase.MATTER_DOMINATED

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    def proper_age(self) -> float:
        """Return the elapsed proper time since the current cycle began, in years.

        Returns:
            Cycle age in years.
        """
        return self._age_s / _YEAR_S

    def is_endpoint_A(self) -> bool:
        """Check whether the universe has reached End-state A (full re-saturation).

        End-state A: all matter has accreted into one universe-spanning
        saturated region — indistinguishable from the Big Bang initial state
        (§18.42, §18.43).

        The practical criterion here is σ ≥ σ_max, which is triggered when
        the scale factor has collapsed back to a ≲ 1 and strain reaches the
        elastic limit.

        Returns:
            ``True`` if σ ≥ σ_max (re-saturation reached).

        References:
            §18.43 — End-state A restart.
        """
        return self.sigma >= _SIGMA_MAX - 1e-12

    def is_endpoint_B(self) -> bool:
        """Check whether the universe has reached End-state B (heat death).

        End-state B: all BHs have Hawking-evaporated to radiation, and the
        substrate has returned to baseline strain σ₀ (§18.43). The universe
        is in "heat death" — an essentially empty, maximally-expanded medium.

        The practical criterion: the universe phase is HEAT_DEATH.

        Returns:
            ``True`` if the universe has reached heat death.

        References:
            §18.43 — End-state B restart via quantum nucleation.
        """
        return self.current_phase() == UniversePhase.HEAT_DEATH

    # ------------------------------------------------------------------
    # String representation
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        H = self.friedmann_equation(self.a)
        return (
            f"FRWUniverse(σ={self.sigma:.3e}, a={self.a:.4f}, "
            f"H={H:.3e} s⁻¹, age={self.proper_age():.3e} yr, "
            f"phase={self.current_phase().name})"
        )


# ---------------------------------------------------------------------------
# CyclicCosmology manager
# ---------------------------------------------------------------------------


@dataclass
class CyclicCosmology:
    """Manager for a sequence of cosmic cycles on the eternal substrate.

    The substrate itself (parameters K, ρ, ξ, Möbius half-flux) persists
    across all cycles.  Each call to ``restart_cycle`` or ``simulate_n_cycles``
    initialises a fresh ``FRWUniverse`` pattern in the same medium.

    Args:
        universe: The initial (or current) ``FRWUniverse`` instance.

    Attributes:
        cycle_count:    Number of completed cycles so far (not counting the
                        current running cycle).
        cycle_durations: List of durations (in years) for each completed cycle.

    References:
        §18.43 — Substrate persists eternally; matter pattern resets.
    """

    universe: FRWUniverse = field(default_factory=FRWUniverse)
    cycle_count: int = field(default=0, init=False)
    cycle_durations: list[float] = field(default_factory=list, init=False)

    # ------------------------------------------------------------------
    # Cycle management
    # ------------------------------------------------------------------

    def restart_cycle(
        self,
        restart_method: Literal["A", "B"],
    ) -> None:
        """Record the current cycle's completion and initialise a new one.

        Both restart methods produce the same substrate state (saturated
        medium, σ = ½) — the distinction is the physical path that led there:

        - Method ``"A"`` (re-saturation): all matter merged into one
          universe-spanning BH, which IS the Big Bang state (§18.42).
          The new cycle begins immediately (w = -1 equation of state drives
          exponential expansion of the saturated medium, §18.40).

        - Method ``"B"`` (quantum nucleation): after heat death, quantum
          fluctuations of the baseline σ₀ medium nucleate a high-strain
          region that reaches saturation, seeding a new Big Bang (§18.43).
          The waiting time for nucleation vastly exceeds all prior timescales
          (but is finite), so the new cycle starts in the same saturated
          initial state.

        After a restart the internal ``FRWUniverse`` is replaced by a fresh
        instance (σ = 0.5, a = 1, age = 0).

        Args:
            restart_method: ``"A"`` for saturation re-collapse, ``"B"`` for
                            Hawking-evaporation + quantum-nucleation path.

        Raises:
            ValueError: If ``restart_method`` is not ``"A"`` or ``"B"``.

        References:
            §18.42 — BB = BH interior; saturation IS the initial state.
            §18.43 — Both end-states restart via the same saturated medium.
        """
        if restart_method not in ("A", "B"):
            raise ValueError(
                f"restart_method must be 'A' or 'B'; got {restart_method!r}"
            )

        # Record completed cycle duration
        completed_age_yr = self.universe.proper_age()
        self.cycle_durations.append(completed_age_yr)
        self.cycle_count += 1

        # Initialise fresh cycle — same substrate parameters, new matter pattern.
        # The cycle starts at the de-saturation transition (a = a_CMB = 1/1100,
        # σ = 0.5), identical to the CMB epoch of our current universe (§18.44).
        self.universe = FRWUniverse(
            sigma=_SIGMA_MAX,
            a=1.0 / 1100.0,
            omega_matter=self.universe.omega_matter,
            omega_lambda=self.universe.omega_lambda,
            omega_radiation=self.universe.omega_radiation,
            H_0=self.universe.H_0,
        )

    def simulate_n_cycles(
        self,
        n: int,
        dt_s: float = 1e15,
        max_steps_per_cycle: int = 10_000,
    ) -> None:
        """Run up to n complete cosmic cycles from the current state.

        Each cycle is evolved with time steps of ``dt_s`` seconds until either
        End-state A or End-state B is detected, or ``max_steps_per_cycle`` is
        exhausted.  When a cycle ends, ``restart_cycle`` is called
        automatically.

        Args:
            n:                    Number of cycles to complete.
            dt_s:                 Euler time step in seconds.  Default 10¹⁵ s
                                  (≈ 31.7 million years) — coarse but adequate
                                  for order-of-magnitude cycle mapping.
            max_steps_per_cycle:  Safety cap on steps per cycle to prevent
                                  infinite loops in pathological parameter regimes.

        Notes:
            This is an order-of-magnitude simulation.  FRW evolution from
            today's epoch (a = 1) to heat death spans many Hubble times; the
            default step size is intentionally large to make the demonstration
            tractable.  Production use should reduce dt_s and increase
            max_steps_per_cycle accordingly.

        References:
            §18.43 — Cyclic cosmology with substrate-eternal background.
        """
        for _ in range(n):
            endpoint_reached = False
            for _step in range(max_steps_per_cycle):
                self.universe.evolve(dt_s)
                if self.universe.is_endpoint_A():
                    self.restart_cycle("A")
                    endpoint_reached = True
                    break
                if self.universe.is_endpoint_B():
                    self.restart_cycle("B")
                    endpoint_reached = True
                    break
            if not endpoint_reached:
                # Step cap exhausted: record the partial cycle as a lower bound
                # on the cycle duration and restart for the next iteration.
                completed_age_yr = self.universe.proper_age()
                self.cycle_durations.append(completed_age_yr)
                self.cycle_count += 1
                self.universe = FRWUniverse(
                    sigma=_SIGMA_MAX,
                    a=1.0 / 1100.0,
                    omega_matter=self.universe.omega_matter,
                    omega_lambda=self.universe.omega_lambda,
                    omega_radiation=self.universe.omega_radiation,
                    H_0=self.universe.H_0,
                )
