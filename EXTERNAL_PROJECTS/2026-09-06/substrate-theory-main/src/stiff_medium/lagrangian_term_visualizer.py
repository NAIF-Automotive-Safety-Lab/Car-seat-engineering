"""Lagrangian-term visualizer — paired GEOMETRY + SIM + VIZ for the substrate
Lagrangian field equations themselves.

THE SUBSTRATE LAGRANGIAN
------------------------
The full substrate Lagrangian density (§18.11), written for a scalar field
u(x, t) with units chosen so that c² = K/ρ when only kinetic + gradient act:

    L = ½ρ(∂_t u)²  −  ½K|∇u|²  −  V(u)  −  γ u (∂_t u)
        ───┬───      ───┬───     ──┬──     ──────┬──────
       kinetic       gradient   potential       drag

with the canonical sine-Gordon-like potential

    V(u) = (K/ξ²) (1 − cos u)

The Euler-Lagrange equation derived from L is

    ρ ∂_t² u  −  K ∇² u  +  (K/ξ²) sin u  +  γ ∂_t u  =  0

EACH TERM IS RESPONSIBLE FOR A DIFFERENT PIECE OF DYNAMICS:

  ½ρ(∂_t u)²    — TIME-INERTIA. Without it the field has no momentum and
                   collapses to V(u) minimum instantaneously. It is what
                   makes ω² = c² k² + m² instead of just ω² = m².

  ½K|∇u|²       — SPATIAL STIFFNESS. Without it neighbouring sites are
                   uncoupled; energy never propagates so there are no waves
                   at all (each lattice site is an independent oscillator
                   in its V(u)).

  V(u)          — CONFINING POTENTIAL. Without it the field has no
                   restoring force at u≠0, so excitations *only* propagate
                   (free wave) and never bind. With V there is a vacuum
                   and the linearised mass m² = V''(0) = K/ξ².

  γ u ∂_t u     — DRAG. Without it dynamics are exactly conservative
                   (energy stays put forever). With γ>0 the field
                   monotonically dissipates to the V minimum, and γ is
                   what gives moving excitations an effective inertial
                   mass m_eff = γ/κ in the cone-bouncing protocol.

WHAT THIS MODULE DOES
---------------------
1. ``LagrangianTermGeometry``    — symbolic identification of the four
   Lagrangian terms together with each one's physical meaning, sign in
   the equation of motion, and which Euler-Lagrange contribution it
   produces.
2. ``TermContributionSimulator`` — 1D substrate evolution under each
   subset of (kinetic, gradient, potential, drag), so you can SEE which
   piece of dynamics is lost when a given term is omitted.

The four canonical sub-Lagrangians produced are:

    KINETIC + GRADIENT      → free wave (massless, undamped)
    KINETIC + POTENTIAL     → bound oscillator (no spatial coupling)
    KINETIC + DRAG          → exponentially decaying particle
    KINETIC + GRADIENT + POTENTIAL + DRAG  → full particle dynamics

VISUALS PRODUCED BY ``render_lagrangian_terms``
-----------------------------------------------
32_lagrangian_terms.png       — 4 panels: free wave (kinetic+gradient),
                                bound state (kinetic+potential), exponential
                                decay (kinetic+drag), full particle.
33_term_contributions.png     — Energy-flow diagram with each Lagrangian
                                term labeled and its physical role spelled
                                out alongside.

REFERENCES
----------
src/stiff_medium/lattice_substrate_2d.py   -- 2D solver of the same L.
src/stiff_medium/substrate_field_solver.py -- 1D static-kink validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Symbolic descriptors of the four Lagrangian terms
# ---------------------------------------------------------------------------

TERM_NAMES: Tuple[str, str, str, str] = ("kinetic", "gradient", "potential", "drag")
"""Canonical names of the four Lagrangian terms, in the standard order."""


# ---------------------------------------------------------------------------
# 1. Lagrangian-term geometry (symbolic identification)
# ---------------------------------------------------------------------------

@dataclass
class LagrangianTermGeometry:
    """Identify each term in L and its physical meaning.

    The dataclass itself does not carry numerical state for the field;
    it is the *symbolic* description used by the simulator and the
    renderer to label what each piece of L = ½ρ(∂_t u)² − ½K|∇u|² −
    V(u) − γ u(∂_t u) actually does.

    Attributes
    ----------
    rho : float
        Substrate mass density. Sets the time-inertia coefficient.
    K : float
        Stiffness. Sets the spatial-coupling coefficient and (with ξ)
        the linearised mass m² = K/ξ².
    xi : float
        Coherence length. Sets the well-width of V(u) and hence the
        linear mass scale.
    gamma : float
        Drag coefficient. γ > 0 dissipates energy monotonically to the
        V-minimum; γ = 0 gives exactly conservative dynamics.
    """

    rho: float = 1.0
    K: float = 1.0
    xi: float = 1.0
    gamma: float = 0.1

    # ------------------------------------------------------------------ Lagrangian density at a point

    def lagrangian_density(self, u: float, du_dt: float, du_dx: float) -> float:
        """Evaluate L at a single (u, ∂_t u, ∂_x u).

        L = ½ρ(∂_t u)² − ½K(∂_x u)² − V(u) − γ u(∂_t u)
        """
        kin = 0.5 * self.rho * du_dt * du_dt
        grad = 0.5 * self.K * du_dx * du_dx
        pot = (self.K / self.xi ** 2) * (1.0 - np.cos(u))
        drag = self.gamma * u * du_dt
        return float(kin - grad - pot - drag)

    # ------------------------------------------------------------------ term descriptors

    def term_descriptors(self) -> Dict[str, Dict[str, str]]:
        """Return symbolic descriptor of each term in L.

        Each descriptor carries:
          - ``symbol``     -- the LaTeX-ish term as it appears in L.
          - ``meaning``    -- one-line physical role.
          - ``el_term``    -- contribution to the Euler-Lagrange equation
                              ρ∂_t²u − K∇²u + V'(u) + γ∂_t u = 0.
          - ``omitted``    -- what dynamics is lost if this term is
                              removed from L.
        """
        return {
            "kinetic": {
                "symbol": "½ρ(∂_t u)²",
                "meaning": "time-inertia: gives the field momentum",
                "el_term": "ρ ∂_t² u",
                "omitted": "no propagation — only quasi-static collapse to V minimum",
            },
            "gradient": {
                "symbol": "−½K|∇u|²",
                "meaning": "spatial stiffness: couples neighbouring sites",
                "el_term": "−K ∇² u",
                "omitted": "no waves — every site is an independent oscillator",
            },
            "potential": {
                "symbol": "−V(u)  with  V = (K/ξ²)(1 − cos u)",
                "meaning": "confining well: gives a vacuum and a linear mass m²=K/ξ²",
                "el_term": "+V'(u) = +(K/ξ²) sin u",
                "omitted": "free wave only — no bound state, no rest energy",
            },
            "drag": {
                "symbol": "−γ u (∂_t u)",
                "meaning": "dissipation: gives moving excitations effective mass γ/κ",
                "el_term": "+γ ∂_t u",
                "omitted": "exactly conservative — no decay, no inertial mass scale",
            },
        }

    # ------------------------------------------------------------------ Euler-Lagrange contributions

    def euler_lagrange_contributions(self) -> Dict[str, str]:
        """Map term name → its Euler-Lagrange contribution (text form)."""
        return {k: v["el_term"] for k, v in self.term_descriptors().items()}

    def linear_mass_squared(self) -> float:
        """m² = V''(0) = K/ξ². The linearised dispersion ω² = c²k² + m²."""
        return float(self.K / self.xi ** 2)

    def wave_speed(self) -> float:
        """c = √(K/ρ). Phase / group speed of the K + ρ kinetic-gradient pair."""
        return float(np.sqrt(self.K / self.rho))

    def damping_rate(self) -> float:
        """κ = γ / ρ. Exponential decay rate of a uniform mode under pure drag."""
        return float(self.gamma / self.rho)

    # ------------------------------------------------------------------ summary

    def summary(self) -> Dict[str, object]:
        return {
            "rho": self.rho,
            "K": self.K,
            "xi": self.xi,
            "gamma": self.gamma,
            "wave_speed": self.wave_speed(),
            "linear_mass_sq": self.linear_mass_squared(),
            "damping_rate": self.damping_rate(),
            "terms": list(TERM_NAMES),
            "descriptors": self.term_descriptors(),
        }


# ---------------------------------------------------------------------------
# 2. Term-contribution simulator (1D substrate evolution per term-subset)
# ---------------------------------------------------------------------------

@dataclass
class TermContributionSimulator:
    """1D substrate evolution under selectable Lagrangian-term subsets.

    The full Euler-Lagrange equation derived from
    L = ½ρ(∂_t u)² − ½K|∇u|² − V(u) − γ u(∂_t u) is

        ρ ∂_t² u = K ∇² u  −  V'(u)  −  γ ∂_t u

    Each term is multiplied by a 0/1 flag so we can turn it on or off:

        ρ ∂_t² u = K λ_grad ∇² u  −  λ_pot V'(u)  −  γ λ_drag ∂_t u

    The kinetic flag λ_kin scales the LHS (effectively setting an
    over-damped quasi-static limit when λ_kin = 0; we then drop the
    second-derivative term and integrate the first-order equation
    γ ∂_t u = K ∇² u − V'(u) instead). All four single-term and
    multi-term combinations are supported.

    Attributes
    ----------
    geometry : LagrangianTermGeometry
        Coefficients (ρ, K, ξ, γ) used to assemble the equation of motion.
    Nx : int
        Number of spatial grid points.
    L : float
        Total length of the 1D box.
    dt : float
        Time-step. Must satisfy a CFL-like bound when the gradient term
        is active: dt < dx / c with c = √(K/ρ).
    bc : str
        Boundary condition. ``'periodic'`` or ``'neumann'``.
    """

    geometry: LagrangianTermGeometry = field(default_factory=LagrangianTermGeometry)
    Nx: int = 256
    L: float = 40.0
    dt: float = 0.02
    bc: str = "periodic"

    # State
    x: NDArray[np.float64] = field(init=False)
    dx: float = field(init=False)
    u: NDArray[np.float64] = field(init=False)
    v: NDArray[np.float64] = field(init=False)
    t: float = field(init=False, default=0.0)

    def __post_init__(self) -> None:
        if self.bc not in ("periodic", "neumann"):
            raise ValueError(f"Unknown BC: {self.bc!r}")
        self.dx = self.L / self.Nx
        self.x = (np.arange(self.Nx) - self.Nx / 2) * self.dx
        self.u = np.zeros(self.Nx, dtype=np.float64)
        self.v = np.zeros(self.Nx, dtype=np.float64)
        self.t = 0.0

    # ------------------------------------------------------------------ operators

    def laplacian(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        if self.bc == "periodic":
            up = np.roll(u, -1)
            um = np.roll(u, +1)
        else:  # neumann
            up = np.empty_like(u)
            um = np.empty_like(u)
            up[:-1] = u[1:]; up[-1] = u[-1]
            um[1:] = u[:-1]; um[0] = u[0]
        return (up + um - 2.0 * u) / (self.dx ** 2)

    def potential_force(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        """V'(u) for V = (K/ξ²)(1 − cos u)."""
        return (self.geometry.K / self.geometry.xi ** 2) * np.sin(u)

    # ------------------------------------------------------------------ initial conditions

    def reset(self) -> None:
        self.u.fill(0.0)
        self.v.fill(0.0)
        self.t = 0.0

    def gaussian_pulse(self, x0: float = 0.0, amplitude: float = 0.3,
                       width: float = 2.0, velocity: float = 0.0) -> None:
        """Initial small-amplitude Gaussian. ``velocity`` boosts to a
        right-moving pulse via ∂_t u = -v ∂_x u (the d'Alembert form)."""
        self.reset()
        self.u = amplitude * np.exp(-((self.x - x0) ** 2) / (2.0 * width ** 2))
        if velocity != 0.0:
            # ∂_t u = -v ∂_x u for a right-moving pulse f(x - v t)
            du_dx = np.gradient(self.u, self.dx)
            self.v = -velocity * du_dx

    def uniform_displacement(self, amplitude: float = 0.3) -> None:
        """Spatially uniform u₀ — useful for testing kinetic+potential
        bound oscillation (no gradient acts because ∇²u = 0)."""
        self.reset()
        self.u.fill(amplitude)

    # ------------------------------------------------------------------ evolution

    def step(self,
             include_kinetic: bool = True,
             include_gradient: bool = True,
             include_potential: bool = True,
             include_drag: bool = True) -> None:
        """Advance one time-step under the selected term subset.

        Two integrators are used:

        * ``include_kinetic = True``: drag-aware velocity-Verlet on
          ρ ∂_t² u = K λ_grad ∇² u − λ_pot V'(u) − γ λ_drag ∂_t u.

        * ``include_kinetic = False``: drop the second-derivative
          term and integrate the first-order over-damped form
          γ ∂_t u = K λ_grad ∇² u − λ_pot V'(u). Requires γ>0 (since
          otherwise there is no ∂_t u in the equation at all).
        """
        K, rho, gamma = self.geometry.K, self.geometry.rho, self.geometry.gamma
        dt = self.dt

        lam_grad = 1.0 if include_gradient else 0.0
        lam_pot = 1.0 if include_potential else 0.0
        lam_drag = 1.0 if include_drag else 0.0

        if include_kinetic:
            # Conservative force per unit mass
            def accel(u_arg: NDArray[np.float64]) -> NDArray[np.float64]:
                f = K * lam_grad * self.laplacian(u_arg)
                f = f - lam_pot * self.potential_force(u_arg)
                return f / rho

            drag_half = float(np.exp(-gamma * lam_drag * dt / (2.0 * rho)))
            a_n = accel(self.u)
            v_half = (self.v + 0.5 * dt * a_n) * drag_half
            self.u = self.u + dt * v_half
            a_np1 = accel(self.u)
            self.v = v_half * drag_half + 0.5 * dt * a_np1
        else:
            # Over-damped first-order equation: γ ∂_t u = K ∇²u − V'(u)
            if gamma <= 0.0 or lam_drag == 0.0:
                raise ValueError(
                    "Without the kinetic term, the drag term must be active "
                    "and gamma > 0 (otherwise there is no time evolution)."
                )
            f = K * lam_grad * self.laplacian(self.u)
            f = f - lam_pot * self.potential_force(self.u)
            self.u = self.u + (dt / gamma) * f
            self.v = f / gamma  # bookkeeping: instantaneous ∂_t u

        self.t += dt

    def evolve(self, steps: int,
               include_kinetic: bool = True,
               include_gradient: bool = True,
               include_potential: bool = True,
               include_drag: bool = True) -> None:
        for _ in range(int(steps)):
            self.step(
                include_kinetic=include_kinetic,
                include_gradient=include_gradient,
                include_potential=include_potential,
                include_drag=include_drag,
            )

    # ------------------------------------------------------------------ energies

    def kinetic_energy(self) -> float:
        return float(0.5 * self.geometry.rho * np.sum(self.v ** 2) * self.dx)

    def gradient_energy(self) -> float:
        # |∂_x u|² via central differences
        if self.bc == "periodic":
            ux = (np.roll(self.u, -1) - np.roll(self.u, +1)) / (2.0 * self.dx)
        else:
            ux = np.zeros_like(self.u)
            ux[1:-1] = (self.u[2:] - self.u[:-2]) / (2.0 * self.dx)
            ux[0] = (self.u[1] - self.u[0]) / self.dx
            ux[-1] = (self.u[-1] - self.u[-2]) / self.dx
        return float(0.5 * self.geometry.K * np.sum(ux ** 2) * self.dx)

    def potential_energy(self) -> float:
        return float((self.geometry.K / self.geometry.xi ** 2)
                     * np.sum(1.0 - np.cos(self.u)) * self.dx)

    def total_energy(self) -> float:
        return self.kinetic_energy() + self.gradient_energy() + self.potential_energy()

    # ------------------------------------------------------------------ scenarios

    def run_free_wave(self, steps: int = 200) -> Dict[str, NDArray[np.float64]]:
        """KINETIC + GRADIENT only → free (massless) wave that propagates.

        Returns dict with ``x``, ``u_initial``, ``u_final``, and ``times``,
        ``snapshots`` arrays for plotting.
        """
        self.gaussian_pulse(x0=0.0, amplitude=0.3, width=2.0, velocity=0.0)
        u0 = self.u.copy()
        snapshots: List[NDArray[np.float64]] = [u0]
        times: List[float] = [0.0]
        save_every = max(1, steps // 6)
        for k in range(steps):
            self.step(include_kinetic=True, include_gradient=True,
                      include_potential=False, include_drag=False)
            if (k + 1) % save_every == 0:
                snapshots.append(self.u.copy())
                times.append(self.t)
        return {
            "x": self.x.copy(),
            "u_initial": u0,
            "u_final": self.u.copy(),
            "times": np.array(times, dtype=float),
            "snapshots": np.array(snapshots, dtype=float),
        }

    def run_bound_oscillation(self, steps: int = 600,
                              amplitude: float = 0.3) -> Dict[str, NDArray[np.float64]]:
        """KINETIC + POTENTIAL only (uniform mode) → harmonic bound state.

        With ∇²u = 0 (uniform initial condition) and no drag, the field
        obeys ρ ü = −V'(u) ≈ −(K/ξ²) u for small u, which is a simple
        harmonic oscillator with ω² = K/(ρ ξ²).
        """
        self.uniform_displacement(amplitude=amplitude)
        ts = []
        us = []
        for k in range(steps):
            self.step(include_kinetic=True, include_gradient=False,
                      include_potential=True, include_drag=False)
            ts.append(self.t)
            us.append(self.u.mean())
        return {
            "t": np.array(ts, dtype=float),
            "u_mean": np.array(us, dtype=float),
            "omega_pred": float(np.sqrt(self.geometry.K
                                        / (self.geometry.rho * self.geometry.xi ** 2))),
            "amplitude": amplitude,
        }

    def run_drag_decay(self, steps: int = 400,
                       amplitude: float = 0.3) -> Dict[str, NDArray[np.float64]]:
        """KINETIC + DRAG only (uniform mode) → exponential velocity decay.

        With ∇²u = 0, no V, and γ>0 the equation is ρ ü = −γ u̇, so
        u̇(t) = u̇(0) exp(−γ t / ρ). We seed v ≠ 0 (kick) so velocity decays
        from a non-zero initial value.
        """
        self.reset()
        self.v.fill(amplitude)  # uniform initial velocity
        ts = []
        vs = []
        for k in range(steps):
            self.step(include_kinetic=True, include_gradient=False,
                      include_potential=False, include_drag=True)
            ts.append(self.t)
            vs.append(self.v.mean())
        return {
            "t": np.array(ts, dtype=float),
            "v_mean": np.array(vs, dtype=float),
            "decay_rate_pred": float(self.geometry.gamma / self.geometry.rho),
            "v_initial": amplitude,
        }

    def run_full_particle(self, steps: int = 600,
                          v_kick: float = 0.4) -> Dict[str, NDArray[np.float64]]:
        """All four terms → kink + antikink pair: bound, propagating, dissipating.

        Initial condition is a kink + antikink pair so the topological
        charges cancel and the field is periodic-compatible. The kink is
        boosted with ``v_kick`` so drag has nonzero ∂_t u to dissipate;
        otherwise a stationary IC would have nothing for γ to act on.
        """
        # Kink-antikink pair: kink at x=-L/4, antikink at x=+L/4.
        # Periodic-friendly because total topological charge = 0.
        self.reset()
        x0_K = -self.L / 4.0
        x0_AK = +self.L / 4.0
        xi = self.geometry.xi
        kink = 4.0 * np.arctan(np.exp((self.x - x0_K) / xi))
        antikink = -4.0 * np.arctan(np.exp((self.x - x0_AK) / xi))
        self.u = kink + antikink
        self.u = self.u - self.u[0]  # offset so u(boundary) = 0
        # Boost: ∂_t u = -v_kink ∂_x u_kink (right-moving kink, left-moving AK)
        # so the pair sums approach each other.
        if v_kick != 0.0:
            sech_K = 1.0 / np.cosh((self.x - x0_K) / xi)
            sech_AK = 1.0 / np.cosh((self.x - x0_AK) / xi)
            du_dx_K = (2.0 / xi) * sech_K
            du_dx_AK = -(2.0 / xi) * sech_AK
            # Kink moves right at +v_kick; antikink moves left at -v_kick.
            self.v = -v_kick * du_dx_K - (-v_kick) * du_dx_AK
        u0 = self.u.copy()

        snapshots: List[NDArray[np.float64]] = [u0.copy()]
        times: List[float] = [0.0]
        energies: List[Tuple[float, float, float]] = [
            (self.kinetic_energy(), self.gradient_energy(), self.potential_energy())
        ]
        save_every = max(1, steps // 6)
        for k in range(steps):
            self.step(include_kinetic=True, include_gradient=True,
                      include_potential=True, include_drag=True)
            if (k + 1) % save_every == 0:
                snapshots.append(self.u.copy())
                times.append(self.t)
                energies.append(
                    (self.kinetic_energy(), self.gradient_energy(), self.potential_energy())
                )
        return {
            "x": self.x.copy(),
            "u_initial": u0,
            "u_final": self.u.copy(),
            "times": np.array(times, dtype=float),
            "snapshots": np.array(snapshots, dtype=float),
            "energies": np.array(energies, dtype=float),
        }


# ---------------------------------------------------------------------------
# 3. Visualization helpers
# ---------------------------------------------------------------------------

def draw_term_panels(geometry: Optional[LagrangianTermGeometry] = None,
                     fig=None) -> "matplotlib.figure.Figure":
    """Produce the 4-panel figure for 32_lagrangian_terms.png.

    Panel layout:

      (0,0) KINETIC + GRADIENT  → free wave
      (0,1) KINETIC + POTENTIAL → bound oscillation
      (1,0) KINETIC + DRAG      → exponential decay
      (1,1) FULL Lagrangian     → particle dynamics
    """
    import matplotlib.pyplot as plt

    if geometry is None:
        geometry = LagrangianTermGeometry()

    sim = TermContributionSimulator(geometry=geometry, Nx=256, L=40.0, dt=0.02)
    free = sim.run_free_wave(steps=200)

    sim = TermContributionSimulator(geometry=geometry, Nx=256, L=40.0, dt=0.02)
    bound = sim.run_bound_oscillation(steps=600)

    sim = TermContributionSimulator(geometry=geometry, Nx=256, L=40.0, dt=0.02)
    decay = sim.run_drag_decay(steps=400)

    sim = TermContributionSimulator(geometry=geometry, Nx=512, L=40.0, dt=0.02)
    full = sim.run_full_particle(steps=600)

    if fig is None:
        fig = plt.figure(figsize=(13, 10))
    axes = fig.subplots(2, 2)

    # --- (0,0) free wave -------------------------------------------------
    ax = axes[0, 0]
    snaps = free["snapshots"]
    times = free["times"]
    cmap = plt.cm.viridis(np.linspace(0.0, 0.95, len(snaps)))
    for snap, t, c in zip(snaps, times, cmap):
        ax.plot(free["x"], snap, color=c, alpha=0.85, label=f"t={t:.1f}")
    ax.set_title("KINETIC + GRADIENT only → free wave\n"
                 "L = ½ρ(∂_t u)² − ½K|∇u|²")
    ax.set_xlabel("x")
    ax.set_ylabel("u(x, t)")
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    ax.grid(True, alpha=0.3)

    # --- (0,1) bound oscillation ----------------------------------------
    ax = axes[0, 1]
    ax.plot(bound["t"], bound["u_mean"], "b-", lw=1.2, label="u(t)")
    omega = bound["omega_pred"]
    period = 2.0 * np.pi / omega
    ax.axvline(period, color="red", linestyle="--", alpha=0.7,
               label=f"T=2π/ω={period:.2f}")
    ax.set_title("KINETIC + POTENTIAL only → bound oscillation\n"
                 "L = ½ρ(∂_t u)² − V(u)  →  ω² = K/(ρ ξ²)")
    ax.set_xlabel("t")
    ax.set_ylabel("u(t)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- (1,0) exponential decay ----------------------------------------
    ax = axes[1, 0]
    decay_rate = decay["decay_rate_pred"]
    v0 = decay["v_initial"]
    v_pred = v0 * np.exp(-decay_rate * decay["t"])
    ax.plot(decay["t"], decay["v_mean"], "b-", lw=1.4, label="∂_t u (numerical)")
    ax.plot(decay["t"], v_pred, "r--", lw=1.2,
            label=f"v₀ exp(−γt/ρ), κ={decay_rate:.3f}")
    ax.set_title("KINETIC + DRAG only → exponential decay\n"
                 "L = ½ρ(∂_t u)² − γ u(∂_t u)  →  v̇ = −(γ/ρ)v")
    ax.set_xlabel("t")
    ax.set_ylabel("⟨∂_t u⟩")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- (1,1) full Lagrangian ------------------------------------------
    ax = axes[1, 1]
    snaps = full["snapshots"]
    times = full["times"]
    cmap = plt.cm.plasma(np.linspace(0.0, 0.95, len(snaps)))
    for snap, t, c in zip(snaps, times, cmap):
        ax.plot(full["x"], snap, color=c, alpha=0.85, label=f"t={t:.1f}")
    ax.set_title("FULL Lagrangian → particle (kink) dynamics\n"
                 "L = ½ρ(∂_t u)² − ½K|∇u|² − V(u) − γu(∂_t u)")
    ax.set_xlabel("x")
    ax.set_ylabel("u(x, t)")
    ax.legend(fontsize=7, ncol=2, loc="lower right")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Each Lagrangian term turned on separately — "
                 "what dynamics survives, what is lost",
                 fontsize=12, y=0.995)
    fig.tight_layout()
    return fig


def draw_term_contributions(geometry: Optional[LagrangianTermGeometry] = None,
                            fig=None) -> "matplotlib.figure.Figure":
    """Produce the energy-flow diagram for 33_term_contributions.png.

    Layout: a single axes that maps each Lagrangian term to its physical
    role, the ODE contribution, and what is lost when omitted. A small
    bar chart on the right shows the energy partition under the full L
    after a fixed number of evolution steps from a kink-like IC.
    """
    import matplotlib.pyplot as plt

    if geometry is None:
        geometry = LagrangianTermGeometry()

    descriptors = geometry.term_descriptors()

    # Run a short full-Lagrangian evolution to get the energy partition.
    sim = TermContributionSimulator(geometry=geometry, Nx=512, L=40.0, dt=0.02)
    full = sim.run_full_particle(steps=400)
    e_kin, e_grad, e_pot = full["energies"][-1]
    energy_partition = {
        "kinetic": float(e_kin),
        "gradient": float(e_grad),
        "potential": float(e_pot),
    }

    if fig is None:
        fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1.0])
    ax_text = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])

    # ---- Left panel: term-by-term annotations --------------------------
    ax_text.axis("off")

    # Title row
    ax_text.text(0.02, 0.97,
                 "L = ½ρ(∂_t u)²  −  ½K|∇u|²  −  V(u)  −  γ u(∂_t u)",
                 fontsize=14, weight="bold", transform=ax_text.transAxes,
                 family="monospace")

    # Sub-title
    ax_text.text(0.02, 0.92,
                 "Euler-Lagrange:  ρ ∂_t² u  −  K ∇² u  +  V'(u)  +  γ ∂_t u  =  0",
                 fontsize=10, style="italic", transform=ax_text.transAxes,
                 family="monospace")

    # Each term
    palette = {
        "kinetic": "#1f77b4",     # blue
        "gradient": "#2ca02c",    # green
        "potential": "#d62728",   # red
        "drag": "#ff7f0e",        # orange
    }
    y0 = 0.82
    dy = 0.18
    for i, key in enumerate(TERM_NAMES):
        d = descriptors[key]
        y = y0 - i * dy
        col = palette[key]
        # box border
        ax_text.add_patch(plt.Rectangle((0.01, y - 0.16), 0.97, 0.165,
                                         transform=ax_text.transAxes,
                                         facecolor=col, alpha=0.10,
                                         edgecolor=col, linewidth=1.5,
                                         clip_on=False))
        ax_text.text(0.025, y - 0.025, key.upper(),
                     fontsize=12, weight="bold", color=col,
                     transform=ax_text.transAxes)
        ax_text.text(0.135, y - 0.025, d["symbol"],
                     fontsize=11, family="monospace",
                     transform=ax_text.transAxes)
        ax_text.text(0.025, y - 0.060,
                     f"meaning  : {d['meaning']}",
                     fontsize=9, transform=ax_text.transAxes)
        ax_text.text(0.025, y - 0.090,
                     f"E-L term : {d['el_term']}",
                     fontsize=9, family="monospace",
                     transform=ax_text.transAxes)
        ax_text.text(0.025, y - 0.120,
                     f"omitted  : {d['omitted']}",
                     fontsize=9, color="#555555",
                     transform=ax_text.transAxes)

    # ---- Right panel: energy partition under full L -------------------
    keys = ["kinetic", "gradient", "potential"]
    vals = [energy_partition[k] for k in keys]
    colors = [palette[k] for k in keys]
    ax_bar.bar(keys, vals, color=colors, edgecolor="black", linewidth=0.8)
    ax_bar.set_title(f"Energy partition under full L\n"
                     f"after t={full['times'][-1]:.1f}",
                     fontsize=10)
    ax_bar.set_ylabel("Energy [substrate units]")
    ax_bar.grid(True, alpha=0.3, axis="y")
    for k, v in zip(keys, vals):
        ax_bar.text(k, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    ax_bar.text(0.02, 0.95,
                f"γ ≠ 0 dissipates total E\nfrom {sum(full['energies'][0]):.2f}"
                f" → {sum(full['energies'][-1]):.2f}",
                transform=ax_bar.transAxes, fontsize=8,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    fig.suptitle("Substrate Lagrangian — term-by-term contribution to dynamics",
                 fontsize=13, y=0.995)
    fig.tight_layout()
    return fig


def make_lagrangian_term_figures(
    geometry: Optional[LagrangianTermGeometry] = None,
) -> Tuple["matplotlib.figure.Figure", "matplotlib.figure.Figure"]:
    """Build both Lagrangian-term figures (32 and 33).

    Returns
    -------
    (fig_terms, fig_contrib) : tuple of matplotlib figures.
    """
    if geometry is None:
        geometry = LagrangianTermGeometry()
    fig_terms = draw_term_panels(geometry=geometry)
    fig_contrib = draw_term_contributions(geometry=geometry)
    return fig_terms, fig_contrib
