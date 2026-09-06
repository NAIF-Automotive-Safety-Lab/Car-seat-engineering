"""Kink-antikink scattering simulator for the 1D sine-Gordon substrate.

Builds on `substrate_field_solver.py` (static-kink validation) by adding
*time-dependent* dynamics: the full 1+1D sine-Gordon equation

    ∂_t² φ − ∂_x² φ + sin(φ) = 0       (natural units c = ξ = 1)

This module supports three canonical scattering experiments:

1. **Kink-kink (KK) elastic collision**.  Two kinks of the same topological
   charge approach head-on, repel, and pass through (the sine-Gordon
   equation is integrable: solitons emerge with the *same* shape, only
   acquiring a phase shift).  The analytic phase shift for two equal-speed
   kinks colliding head-on is

       δ_KK(v) = 2 ln(γ_v / v_rel)

   where γ_v = 1/√(1−v²) and v_rel ≈ 2 v / (1 + v²) is the relative
   velocity in either kink's rest frame.  More precisely (Faddeev-Korepin):
   the spatial shift each soliton picks up is Δx = 2 ξ ln(1/v_rel) / γ_rel.

2. **Kink-antikink (KA) collision**.  Same charge difference, opposite
   topological sign.  At high velocity they pass through (with phase
   shift); at low velocity they form a bound state — a *breather*.

3. **Breather**.  Bound kink-antikink solution

       φ_B(x, t) = 4 arctan( (η / ω) sin(ω t) / cosh(η x) ),  η = √(1−ω²)

   oscillates at frequency ω with period T = 2π/ω, exactly.

The conserved topological charge

       Q = (1/2π) ∫ ∂_x φ dx

must be preserved to numerical precision under sine-Gordon evolution; the
time-stepper used here (symmetric leap-frog / Stormer-Verlet) is symplectic
in the limit dt → 0 and conserves Q to round-off in practice.

REFERENCES
----------
- Spec §18.11 (Lagrangian), §18.48.3 (kink-antikink potential).
- Faddeev & Korepin, Phys. Rep. 42 (1978) on sine-Gordon scattering.
- Rajaraman, "Solitons and Instantons" (1982), §4.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Analytic profiles
# ---------------------------------------------------------------------------

def lorentz_kink(x: NDArray[np.float64], t: float, x0: float, v: float,
                 xi: float = 1.0, sign: int = +1) -> NDArray[np.float64]:
    """Boosted single (anti)kink at time t.

    sign=+1 → kink   φ = 4 arctan(exp( γ (x−x0−vt)/ξ ))   winding 0 → 2π
    sign=−1 → antikink φ = 4 arctan(exp(−γ (x−x0−vt)/ξ ))  winding 2π → 0

    Args:
        x:    Spatial grid.
        t:    Time.
        x0:   Position at t=0.
        v:    Velocity (|v|<1).
        xi:   Coherence length.
        sign: +1 kink, −1 antikink.
    """
    gamma = 1.0 / np.sqrt(max(1.0 - v * v, 1e-30))
    arg = sign * gamma * (x - x0 - v * t) / xi
    # Clip to avoid overflow in exp; arctan(exp(±large)) saturates to 0 or π/2.
    arg = np.clip(arg, -50.0, 50.0)
    return 4.0 * np.arctan(np.exp(arg))


def lorentz_kink_dot(x: NDArray[np.float64], t: float, x0: float, v: float,
                     xi: float = 1.0, sign: int = +1) -> NDArray[np.float64]:
    """Time derivative ∂_t φ of a boosted (anti)kink.

    d/dt [4 arctan(exp(s γ (x − x0 − v t)/ξ))]
        = 4 · (s γ (−v)/ξ) · sech(s γ (x − x0 − v t)/ξ)
        = −2 s γ v / ξ · sech(arg)         (using d/du arctan(e^u) = ½ sech u)

    Note 4 · ½ · (−sγv/ξ) = −2 s γ v /ξ.  Sign factor s² = 1 cancels in
    sech (even).
    """
    gamma = 1.0 / np.sqrt(max(1.0 - v * v, 1e-30))
    arg = sign * gamma * (x - x0 - v * t) / xi
    arg = np.clip(arg, -50.0, 50.0)
    # sech(u) = 1/cosh(u)
    return -2.0 * sign * gamma * v / xi / np.cosh(arg)


def breather_profile(x: NDArray[np.float64], t: float, omega: float,
                     x0: float = 0.0) -> NDArray[np.float64]:
    """Standing breather: φ_B(x,t) = 4 arctan((η/ω) sin(ω t)/cosh(η (x−x0))).

    Valid for 0 < ω < 1, η = √(1 − ω²).  Period T = 2π/ω.
    """
    if not (0.0 < omega < 1.0):
        raise ValueError("breather requires 0 < omega < 1")
    eta = np.sqrt(1.0 - omega * omega)
    num = (eta / omega) * np.sin(omega * t)
    den = np.cosh(eta * (x - x0))
    return 4.0 * np.arctan(num / den)


def breather_profile_dot(x: NDArray[np.float64], t: float, omega: float,
                         x0: float = 0.0) -> NDArray[np.float64]:
    """∂_t φ for the breather.

    Let f(t) = (η/ω) sin(ω t), g(x) = 1/cosh(η(x−x0)), so φ = 4 arctan(f g).
    ∂_t φ = 4 · g f'/(1 + (fg)²),  f'(t) = η cos(ω t).
    """
    if not (0.0 < omega < 1.0):
        raise ValueError("breather requires 0 < omega < 1")
    eta = np.sqrt(1.0 - omega * omega)
    f = (eta / omega) * np.sin(omega * t)
    fdot = eta * np.cos(omega * t)
    g = 1.0 / np.cosh(eta * (x - x0))
    return 4.0 * g * fdot / (1.0 + (f * g) ** 2)


# ---------------------------------------------------------------------------
# Discretized operators
# ---------------------------------------------------------------------------

def laplacian_1d(phi: NDArray[np.float64], dx: float,
                 bc: str = "open") -> NDArray[np.float64]:
    """Second-order centred Laplacian on a 1D grid.

    bc:
        'open'     — Neumann-style: ghost = boundary (∂_x = 0 at edges).
        'periodic' — wrap-around.
        'dirichlet'— ghost = 0 (field pinned to 0 outside grid).
    """
    lap = np.empty_like(phi)
    lap[1:-1] = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / (dx * dx)
    if bc == "periodic":
        lap[0] = (phi[1] - 2.0 * phi[0] + phi[-1]) / (dx * dx)
        lap[-1] = (phi[0] - 2.0 * phi[-1] + phi[-2]) / (dx * dx)
    elif bc == "open":
        # Neumann: ghost cell equals interior neighbour.
        lap[0] = (phi[1] - 2.0 * phi[0] + phi[1]) / (dx * dx)
        lap[-1] = (phi[-2] - 2.0 * phi[-1] + phi[-2]) / (dx * dx)
    elif bc == "dirichlet":
        lap[0] = (phi[1] - 2.0 * phi[0] + 0.0) / (dx * dx)
        lap[-1] = (0.0 - 2.0 * phi[-1] + phi[-2]) / (dx * dx)
    else:
        raise ValueError(f"unknown bc: {bc}")
    return lap


def topological_charge(phi: NDArray[np.float64], dx: float) -> float:
    """Q = (1/2π) ∫ ∂_x φ dx = (φ(+L) − φ(−L)) / (2π).

    Equivalent to counting net 2π-windings.  For a single static kink
    Q = +1; antikink Q = −1; breather Q = 0; vacuum Q = 0.
    """
    return float((phi[-1] - phi[0]) / (2.0 * np.pi))


def total_energy(phi: NDArray[np.float64], pi: NDArray[np.float64],
                 dx: float) -> float:
    """H = ∫ [½ π² + ½ (∂_x φ)² + (1 − cos φ)] dx, where π = ∂_t φ."""
    grad = np.empty_like(phi)
    grad[1:-1] = (phi[2:] - phi[:-2]) / (2.0 * dx)
    grad[0] = (phi[1] - phi[0]) / dx
    grad[-1] = (phi[-1] - phi[-2]) / dx
    e = 0.5 * pi * pi + 0.5 * grad * grad + (1.0 - np.cos(phi))
    return float(np.trapezoid(e, dx=dx))


# ---------------------------------------------------------------------------
# KinkScattering simulator
# ---------------------------------------------------------------------------

@dataclass
class ScatteringResult:
    """Diagnostic record from a scattering run."""

    t: NDArray[np.float64]
    Q: NDArray[np.float64]            # topological charge vs time
    E: NDArray[np.float64]            # total energy vs time
    centres: NDArray[np.float64]      # kink centre(s) vs time, may be (T,2)
    phi_history: Optional[NDArray[np.float64]] = None  # (T_save, Nx)
    snapshot_t: Optional[NDArray[np.float64]] = None
    extra: dict = field(default_factory=dict)


class KinkScattering:
    """Numerical sine-Gordon scattering on a 1D substrate lattice.

    Time-stepping:  Stormer-Verlet leap-frog (symplectic, 2nd-order).

        π_{n+1/2} = π_n  + (dt/2) (∂_x² φ_n  − sin φ_n)
        φ_{n+1}  = φ_n  + dt π_{n+1/2}
        π_{n+1}  = π_{n+1/2} + (dt/2) (∂_x² φ_{n+1} − sin φ_{n+1})

    Stability requires the CFL condition c · dt / dx < 1 (here c = 1).
    The default dt = 0.4 dx is well inside this bound and preserves Q to
    machine precision in practice.

    Args:
        Nx:   Number of lattice points.
        L:    Half-domain extent (grid spans [−L, +L]).
        xi:   Coherence length (1 in natural units).
        bc:   'open', 'periodic', or 'dirichlet'.
        dt:   Time-step (auto-chosen if None).
        cfl:  CFL safety factor when dt is auto-chosen.
    """

    def __init__(self, Nx: int = 1024, L: float = 80.0, xi: float = 1.0,
                 bc: str = "open", dt: Optional[float] = None,
                 cfl: float = 0.4):
        self.Nx = int(Nx)
        self.L = float(L)
        self.xi = float(xi)
        self.bc = bc
        self.x = np.linspace(-L, L, Nx)
        self.dx = float(self.x[1] - self.x[0])
        self.dt = float(dt) if dt is not None else cfl * self.dx
        self.cfl = cfl
        self.phi: NDArray[np.float64] = np.zeros(Nx)
        self.pi: NDArray[np.float64] = np.zeros(Nx)
        self.t = 0.0

    # ------------------------------------------------------------------
    # Initialisers
    # ------------------------------------------------------------------

    def init_two_kinks(self, x_K: float, v_K: float,
                       x_AK: float, v_AK: float,
                       sign_K: int = +1, sign_AK: int = -1) -> None:
        """Initialise field as a superposition of two boosted (anti)kinks.

        For the kink-antikink (sign_K = +1, sign_AK = −1) ansatz we use
        φ = φ_K + φ_AK − 2π    so endpoints both → 0  (Q = 0).
        For the kink-kink ansatz φ = φ_K + φ_K2 → 0 → 4π (Q = 2).
        """
        phi_1 = lorentz_kink(self.x, 0.0, x_K, v_K, self.xi, sign_K)
        phi_2 = lorentz_kink(self.x, 0.0, x_AK, v_AK, self.xi, sign_AK)
        pi_1 = lorentz_kink_dot(self.x, 0.0, x_K, v_K, self.xi, sign_K)
        pi_2 = lorentz_kink_dot(self.x, 0.0, x_AK, v_AK, self.xi, sign_AK)
        if sign_K == +1 and sign_AK == -1:
            # Kink (0→2π) + antikink (2π→0)−2π gives 0→0 with bump ~ 2π.
            self.phi = phi_1 + phi_2 - 2.0 * np.pi
        elif sign_K == -1 and sign_AK == +1:
            # Antikink (2π→0) + kink (0→2π) − 0  → bump dipping below 0.
            # Adjust: subtract 0 (already 2π → 2π asymptotic? No: AK → 0 left,
            # kink → 2π right, sum = 2π everywhere with a dip in middle).
            self.phi = phi_1 + phi_2
        else:
            # Same-sign: kink + kink → 0→4π (Q=2), or AK+AK → 4π→0 (Q=−2).
            self.phi = phi_1 + phi_2
        self.pi = pi_1 + pi_2
        self.t = 0.0

    def init_breather(self, omega: float, x0: float = 0.0) -> None:
        """Initialise the analytic breather solution at t = 0."""
        self.phi = breather_profile(self.x, 0.0, omega, x0)
        self.pi = breather_profile_dot(self.x, 0.0, omega, x0)
        self.t = 0.0

    def breather_initial(self, omega: float, x0: float = 0.0) -> None:
        """Alias matching the spec request name."""
        self.init_breather(omega, x0)

    # ------------------------------------------------------------------
    # Conserved diagnostics
    # ------------------------------------------------------------------

    def topological_charge(self) -> float:
        return topological_charge(self.phi, self.dx)

    def total_energy(self) -> float:
        return total_energy(self.phi, self.pi, self.dx)

    def find_kink_centres(self, n_expected: int = 2) -> NDArray[np.float64]:
        """Locate kink centres by finding zero crossings of (φ − π) (mod 2π).

        For a single kink centred at x0, φ(x0) = π.  We threshold ∂_x φ to
        identify regions of rapid winding, then take the weighted-mean x of
        ∂_x φ within each region.  Returns up to n_expected centre positions.
        """
        gradphi = np.gradient(self.phi, self.dx)
        weight = np.abs(gradphi)
        if weight.sum() < 1e-6:
            return np.array([])
        # Split into n_expected lobes by clustering on gradient-density.
        thresh = 0.1 * weight.max()
        mask = weight > thresh
        if not mask.any():
            return np.array([])
        # Group contiguous-True segments
        idx = np.where(mask)[0]
        groups = np.split(idx, np.where(np.diff(idx) > 3)[0] + 1)
        centres = []
        for g in groups:
            if g.size < 2:
                continue
            w = weight[g]
            c = float((self.x[g] * w).sum() / w.sum())
            centres.append(c)
        return np.asarray(centres)

    # ------------------------------------------------------------------
    # Time-stepping
    # ------------------------------------------------------------------

    def _accel(self, phi: NDArray[np.float64]) -> NDArray[np.float64]:
        return laplacian_1d(phi, self.dx, self.bc) - np.sin(phi)

    def step(self, n: int = 1) -> None:
        """Advance n leap-frog steps in place."""
        dt = self.dt
        a = self._accel(self.phi)
        for _ in range(n):
            # half kick
            self.pi = self.pi + 0.5 * dt * a
            # drift
            self.phi = self.phi + dt * self.pi
            # second half kick uses *new* acceleration
            a = self._accel(self.phi)
            self.pi = self.pi + 0.5 * dt * a
            self.t += dt

    def run(self, t_end: float, n_snapshots: int = 0,
            track_centres: bool = True,
            n_expected_kinks: int = 2) -> ScatteringResult:
        """Evolve until t = t_end, recording diagnostics.

        Args:
            t_end:           Final simulation time.
            n_snapshots:     Number of φ snapshots to retain (0 = none).
            track_centres:   Record kink centre positions each diagnostic step.
            n_expected_kinks: Hint for centre-finding.

        Returns:
            ScatteringResult.
        """
        n_total = int(np.ceil(t_end / self.dt))
        # Diagnostic every ~200 steps (or every step if total is small).
        diag_every = max(1, n_total // 400)
        n_diag = n_total // diag_every + 1

        t_arr = np.zeros(n_diag)
        Q_arr = np.zeros(n_diag)
        E_arr = np.zeros(n_diag)
        centres_arr = np.full((n_diag, max(1, n_expected_kinks)), np.nan)
        snap_phi: list = []
        snap_t: list = []
        snap_every = (n_total // max(1, n_snapshots)) if n_snapshots > 0 else None

        for i_step in range(n_total + 1):
            if i_step % diag_every == 0:
                k = i_step // diag_every
                if k >= n_diag:
                    break
                t_arr[k] = self.t
                Q_arr[k] = self.topological_charge()
                E_arr[k] = self.total_energy()
                if track_centres:
                    cs = self.find_kink_centres(n_expected_kinks)
                    for j in range(min(len(cs), centres_arr.shape[1])):
                        centres_arr[k, j] = cs[j]
            if snap_every is not None and i_step % snap_every == 0:
                snap_phi.append(self.phi.copy())
                snap_t.append(self.t)
            if i_step < n_total:
                self.step(1)

        return ScatteringResult(
            t=t_arr, Q=Q_arr, E=E_arr,
            centres=centres_arr,
            phi_history=np.array(snap_phi) if snap_phi else None,
            snapshot_t=np.array(snap_t) if snap_t else None,
        )

    # ------------------------------------------------------------------
    # High-level scattering experiments
    # ------------------------------------------------------------------

    def kk_collision(self, v: float, separation: float = 30.0,
                     t_end: Optional[float] = None) -> ScatteringResult:
        """Head-on kink-kink (same charge) collision at speeds ±v.

        Both kinks have topological charge +1, so total Q = +2 throughout.
        Returns the scattering record (centres before/after give phase shift
        via `phase_shift`).
        """
        x_K = -separation / 2.0
        x_K2 = +separation / 2.0
        self.init_two_kinks(x_K, +v, x_K2, -v, sign_K=+1, sign_AK=+1)
        # Default run time: enough to collide and separate again
        if t_end is None:
            t_end = 1.5 * separation / max(abs(v), 0.05)
        return self.run(t_end, n_snapshots=20, n_expected_kinks=2)

    def ka_collision(self, v: float, separation: float = 30.0,
                     t_end: Optional[float] = None) -> ScatteringResult:
        """Head-on kink-antikink collision at speeds ±v.

        Outcome depends on v and on the kinetic-vs-rest-mass ratio:
            v low  → bound state (breather) forms.
            v high → solitons pass through (with phase shift).
        Q = 0 throughout.
        """
        self.init_two_kinks(-separation / 2.0, +v, +separation / 2.0, -v,
                            sign_K=+1, sign_AK=-1)
        if t_end is None:
            t_end = 2.0 * separation / max(abs(v), 0.05)
        return self.run(t_end, n_snapshots=30, n_expected_kinks=2)

    def phase_shift(self, v: float, separation: float = 40.0,
                    channel: str = "KK") -> dict:
        """Measure the phase shift δ(v) of head-on (anti)kink collision.

        Method: track the trailing kink's position before and after
        collision; subtract the free trajectory (linear in t) to obtain the
        spatial shift Δx, then convert to a dimensionless phase
            δ = Δx / ξ.

        Analytic comparison (sine-Gordon, two equal-speed kinks):
            δ_KK^analytic(v) = 2 ln( (1 + v) / (1 − v) ) / 2          (kink-kink)
                            ≈ 2 ln(γ_v / v_rel) at non-relativistic v
        We return both the numeric Δx and the analytic prediction.
        """
        if channel not in ("KK", "KA"):
            raise ValueError("channel must be 'KK' or 'KA'")
        # Run a fresh sim
        sign1, sign2 = (+1, +1) if channel == "KK" else (+1, -1)
        self.init_two_kinks(-separation / 2.0, +v, +separation / 2.0, -v,
                            sign_K=sign1, sign_AK=sign2)
        t_end = 2.0 * separation / max(abs(v), 0.05)
        res = self.run(t_end, n_snapshots=0, n_expected_kinks=2)

        # Identify "before" and "after" windows by time, excluding the
        # collision interval where centres merge / NaN out.
        t = res.t
        c = res.centres
        n = len(t)
        # Approx collision time: kinks moving at ±v, separation/2 each.
        t_collide = (separation / 2.0) / max(abs(v), 1e-3)
        before = (t > 0.05 * t_collide) & (t < 0.6 * t_collide)
        after = (t > 1.6 * t_collide) & (t < 2.5 * t_collide)

        # Track leftmost (smallest x) and rightmost (largest x) kink centres
        # in each window via linear fit to the centres trajectory.
        def fit_v(times, xs):
            mask = np.isfinite(xs)
            if mask.sum() < 3:
                return np.nan, np.nan
            slope, intercept = np.polyfit(times[mask], xs[mask], 1)
            return float(slope), float(intercept)

        # Use leftmost-kink trajectory.
        left_before = np.nanmin(c, axis=1)
        v_lb, x0_lb = fit_v(t[before], left_before[before])
        v_la, x0_la = fit_v(t[after], left_before[after])

        # In KK: leftmost kink moves +v before, then after collision is at
        # the position of the kink that emerged on the LEFT (i.e. the
        # original right kink that bounced).  Its free trajectory (no
        # interaction) is the original right kink at +sep/2 moving with -v:
        #     x_free(t) = +sep/2 − v · t.
        # The shift is Δx = x_measured(t) − x_free(t), evaluated mid-after.
        # In KA: solitons pass through; leftmost kink before is the same
        # one that becomes leftmost after, but DELAYED by time spent in
        # the interaction region.  Free traj: x_free(t) = -sep/2 + v · t
        # ... but then it would be at +x post-collision, contradiction.
        # Better: in KA, sign flips; the leftmost late-time kink IS the
        # original antikink that bounced past (kink continues to right).
        # So the leftmost late-time trajectory free-version starts at +sep/2
        # with v=-v: x_free(t) = +sep/2 − v · t.  Same formula as KK.
        t_eval = float(np.mean(t[after]))
        x_meas = v_la * t_eval + x0_la
        x_free = (separation / 2.0) - v * t_eval
        delta_x = x_meas - x_free

        # Analytic
        gamma = 1.0 / np.sqrt(max(1.0 - v * v, 1e-12))
        # Sine-Gordon two-kink lab-frame spatial shift (Faddeev-Korepin /
        # Rajaraman §4): for two equal-speed kinks colliding head-on,
        #     Δx_KK = (2 ξ / γ) · ln(1/v)
        # i.e. each kink is delayed (or, in the bounce picture, displaced)
        # by this much relative to its free trajectory.  This matches the
        # leading non-relativistic asymptote δ_KK ≈ 2 ln(γ/v_rel).
        # KA channel: same magnitude, sign reflects attraction vs repulsion.
        if abs(v) > 0:
            delta_x_analytic_KK = (2.0 / gamma) * np.log(1.0 / abs(v)) * self.xi
        else:
            delta_x_analytic_KK = 0.0
        delta_x_analytic_KA = -delta_x_analytic_KK

        return {
            "v": v,
            "channel": channel,
            "delta_x_numeric": float(delta_x),
            "delta_x_analytic": float(
                delta_x_analytic_KK if channel == "KK" else delta_x_analytic_KA
            ),
            "delta_phase": float(delta_x / self.xi),
            "v_left_after": v_la,
            "x_meas_after": x_meas,
            "x_free_after": x_free,
        }


# ---------------------------------------------------------------------------
# Convenience: characterise KA collision outcome
# ---------------------------------------------------------------------------

def classify_ka_outcome(res: ScatteringResult,
                        kink_mass: float = 8.0) -> str:
    """Heuristic classification of a KA collision: 'breather', 'pass', 'annihilate'.

    Uses centres trajectory and total energy:
      - If post-collision energy ≪ 2 m_K, kinks annihilated into radiation.
      - If centres remain bounded near the origin (oscillating), it's a breather.
      - Otherwise, pass-through.
    """
    n = len(res.t)
    c = res.centres[3 * n // 4:]
    finite = c[np.isfinite(c)]
    if finite.size == 0:
        return "annihilate"
    # Robust localization metric: 90th percentile of |centre| should be
    # small (a few ξ) for a bound state.  Late-time radiation can produce
    # spurious large-|x| centres but they are sparse — use a percentile.
    p90 = float(np.nanpercentile(np.abs(finite), 90))
    if p90 < 6.0:
        return "breather"
    return "pass"
