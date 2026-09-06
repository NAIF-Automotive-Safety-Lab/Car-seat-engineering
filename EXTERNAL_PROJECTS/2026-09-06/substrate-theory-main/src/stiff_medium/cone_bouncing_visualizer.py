"""cone_bouncing_visualizer.py
=================================

Paired GEOMETRY + SIMULATION module that **visualizes** the cone-bouncing
mechanism for rest-mass generation.

Background
----------
A bound substrate excitation cannot exceed the saturation cap σ ≤ 1/2.  The
locus σ = ±1/2 forms a pair of "cone walls" that any traveling envelope
must reflect off.  The internal oscillation frequency between successive
reflections is the canonical cone-bouncing frequency

    ω_b² = (c/ξ)² + (γ/2ρ)²,        c ≡ √(K/ρ).

This frequency sets the rest mass via the Planck–Compton relation

    m c² = ℏ ω_b.

This module provides

    * ``ConeBouncingGeometry`` — analytic cone half-angle, period, mass
      ladder for a swept drag coefficient.
    * ``ConeBouncingSimulator`` — small numerical model of a substrate
      envelope rattling between σ = ±1/2 walls under linear restoring
      force and (optional) drag, used to *measure* the period and check
      it against the analytic ω_b.

The simulator is deliberately a 1-D toy: it integrates a damped harmonic
oscillator whose amplitude is reflected (sign-flipped on velocity) when
it exceeds the cap A_max ∝ φ_max.  This is the simplest faithful
realisation of "cone-bouncing".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from src.stiff_medium.cone_bouncing_protocol import (
    C, HBAR, M_E, MEV, XI_E,
    omega_b_canonical, mass_from_omega,
    electron_anchor_primitives,
)


# ---------------------------------------------------------------------------
# Anchor primitives for charged-lepton scan (Solution A: ξ = λ_C(electron))
# ---------------------------------------------------------------------------

# Particle masses [kg]
M_MU  = 1.883531627e-28      # muon
M_TAU = 3.16754e-27          # tau

# Compton wavelengths (ξ = ℏ / (m c)); all use Solution A primitives.
XI_MU  = HBAR / (M_MU  * C)
XI_TAU = HBAR / (M_TAU * C)


def _saturation_to_half_angle(sigma_max: float) -> float:
    """Cone half-angle θ_c from saturation cap σ_max ≤ 1/2.

    The substrate cannot tilt more than σ_max in strain; the corresponding
    geometric tilt of the local light-cone is

        θ_c = arctan(σ_max / (1 - σ_max))      [tilt from un-saturated axis]

    For σ_max = 1/2:  θ_c = arctan(1) = π/4 = 45°.

    Returns the half-angle in radians.  The full opening angle of the
    cone is 2 θ_c.
    """
    if not 0.0 < sigma_max < 1.0:
        raise ValueError("sigma_max must lie strictly in (0, 1).")
    return float(np.arctan(sigma_max / (1.0 - sigma_max)))


# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------

@dataclass
class ConeBouncingGeometry:
    """Analytic geometry of the cone-bouncing mass mechanism.

    Attributes
    ----------
    K, rho, xi, gamma : float
        Substrate primitives (SI).  Defaults: electron anchor primitives.
    sigma_max : float
        Saturation cap in strain units (default 1/2 per MODEL.md §18.39).

    Notes
    -----
    The "cone" here is the locus σ(x) = ±σ_max in (x, σ) phase space.
    A bound state's strain envelope u(t)/φ_max bounces between these
    two walls at frequency ω_b.  All geometric quantities derive from
    the canonical formula in cone_bouncing_protocol.
    """

    K:         float = field(default_factory=lambda: electron_anchor_primitives()["K"])
    rho:       float = field(default_factory=lambda: electron_anchor_primitives()["rho"])
    xi:        float = field(default_factory=lambda: electron_anchor_primitives()["xi"])
    gamma:     float = 0.0
    sigma_max: float = 0.5

    # ------------------------------------------------------------------
    # Frequency / period / mass
    # ------------------------------------------------------------------

    @property
    def omega_b(self) -> float:
        """Canonical cone-bouncing angular frequency [rad/s]."""
        return omega_b_canonical(self.K, self.rho, self.xi, self.gamma)

    @property
    def period(self) -> float:
        """Bounce period T = 2π / ω_b [s]."""
        omega = self.omega_b
        if omega <= 0:
            return float("inf")
        return float(2.0 * np.pi / omega)

    @property
    def mass_kg(self) -> float:
        """Rest mass m = ℏ ω_b / c² [kg]."""
        return mass_from_omega(self.omega_b)

    @property
    def mass_mev(self) -> float:
        """Rest mass m c² [MeV]."""
        return float(self.mass_kg * C * C / MEV)

    # ------------------------------------------------------------------
    # Cone wall geometry
    # ------------------------------------------------------------------

    @property
    def cone_half_angle(self) -> float:
        """Cone half-angle θ_c [rad] from saturation cap σ_max."""
        return _saturation_to_half_angle(self.sigma_max)

    @property
    def cone_half_angle_deg(self) -> float:
        """Cone half-angle in degrees (for plot labels)."""
        return float(self.cone_half_angle * 180.0 / np.pi)

    def cone_walls(self, t: np.ndarray, c_eff: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Return the upper / lower cone walls at the saturation cap.

        For visualisation we plot the cap level u = ±σ_max · ξ as a
        function of t.  A truly travelling cone wall would be
        u_wall = ±tan(θ_c) · (c_eff t), but for the *bound-state*
        envelope picture a horizontal cap is the relevant locus.

        Parameters
        ----------
        t : ndarray
            Time samples.
        c_eff : float, optional
            Substrate wave speed.  Used only for the slanted-wall mode
            (returned by ``slanted_walls``); ignored here.

        Returns
        -------
        upper, lower : ndarrays of shape ``t.shape``
            Cap heights ±σ_max · ξ.
        """
        cap = self.sigma_max * self.xi
        return np.full_like(t, cap, dtype=float), np.full_like(t, -cap, dtype=float)

    def slanted_walls(self, t: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Return the slanted light-cone walls u = ±tan(θ_c) · c_eff · t.

        These lines bound the *causal* support of any disturbance and are
        used in the geometry panel of the visualisation.
        """
        c_eff = float(np.sqrt(self.K / self.rho))
        slope = np.tan(self.cone_half_angle) * c_eff
        return slope * t, -slope * t

    # ------------------------------------------------------------------
    # Bouncing trajectory (analytic)
    # ------------------------------------------------------------------

    def trajectory(
        self,
        t: np.ndarray,
        amplitude_frac: float = 0.95,
    ) -> np.ndarray:
        """Analytic strain-envelope trajectory u(t).

        The bound mode oscillates as a damped sinusoid between the
        ±σ_max·ξ cap; we expose the *un-clipped* sinusoid here and
        let the Simulator reflect it numerically.

        Parameters
        ----------
        t : ndarray
            Time samples [s].
        amplitude_frac : float
            Amplitude as a fraction of the cap (0 < frac < 1).  The
            envelope is ``frac · σ_max · ξ``.

        Returns
        -------
        u : ndarray
            Strain trajectory [m].
        """
        amp = amplitude_frac * self.sigma_max * self.xi
        Gamma = self.gamma / (2.0 * self.rho)
        return amp * np.exp(-Gamma * t) * np.cos(self.omega_b * t)

    # ------------------------------------------------------------------
    # Mass ladder over γ
    # ------------------------------------------------------------------

    def mass_vs_gamma(self, gammas: np.ndarray) -> Dict[str, np.ndarray]:
        """Sweep γ → compute (ω_b, m) for each value.

        Returns a dict with keys ``gamma`` (input), ``omega_b``, ``mass_kg``,
        ``mass_mev``.
        """
        gammas = np.asarray(gammas, dtype=float)
        omegas = np.array([
            omega_b_canonical(self.K, self.rho, self.xi, float(g))
            for g in gammas
        ])
        masses = np.array([mass_from_omega(o) for o in omegas])
        return {
            "gamma":    gammas,
            "omega_b":  omegas,
            "mass_kg":  masses,
            "mass_mev": masses * C * C / MEV,
        }


# ---------------------------------------------------------------------------
# 2. Simulator
# ---------------------------------------------------------------------------

@dataclass
class ConeBouncingSimulator:
    """Time-evolve a substrate envelope inside the saturation cone.

    The model is a 1-D damped harmonic oscillator whose amplitude is
    *clipped* and *reflected* at u = ±A_max:

        ρ ü + γ u̇ + K_eff u = 0,                 (between walls)
        u → sign(u) · A_max  and  u̇ → -u̇        (on impact)

    where K_eff is chosen so that the **un-clipped** small-oscillation
    frequency equals the canonical ω_b of the geometry.  Reflection off
    the saturation cap corresponds to one half-cycle of cone bouncing.

    Attributes
    ----------
    geometry : ConeBouncingGeometry
        Source of K, ρ, ξ, γ, σ_max.
    amplitude_frac : float
        Initial amplitude as fraction of A_max = σ_max·ξ (default 0.95).
    n_steps : int
        Number of integration steps (default 4000).
    dt_per_period : int
        Steps per analytic period (default 200) — sets dt = T / 200.
    """

    geometry:        ConeBouncingGeometry = field(default_factory=ConeBouncingGeometry)
    amplitude_frac:  float = 0.95
    n_steps:         int = 4000
    dt_per_period:   int = 200

    # ------------------------------------------------------------------
    # Internal: solve trajectory
    # ------------------------------------------------------------------

    def _omega_target(self) -> float:
        """The frequency the simulator is *tuned* to reproduce."""
        return self.geometry.omega_b

    def _dt(self) -> float:
        """Time step Δt = T / dt_per_period."""
        return self.geometry.period / float(self.dt_per_period)

    def run(self) -> Dict[str, np.ndarray]:
        """Integrate the envelope dynamics.

        Returns
        -------
        dict with keys
            ``t``      : ndarray of times [s],
            ``u``      : ndarray of strains [m],
            ``v``      : ndarray of strain rates,
            ``energy`` : ndarray of instantaneous energies (½ρv² + ½K_eff u²),
            ``hits``   : list of indices where a wall reflection occurred,
            ``A_max``  : float saturation amplitude [m].
        """
        omega = self._omega_target()
        if omega <= 0:
            raise ValueError("ω_b must be positive — increase ξ or γ.")

        # Pick K_eff such that √(K_eff / ρ) = ω.  This is the bound-mode
        # spring constant for an oscillator of inertia ρ (per unit volume).
        K_eff = self.geometry.rho * omega * omega

        Gamma = self.geometry.gamma / (2.0 * self.geometry.rho)  # damping rate
        A_max = self.geometry.sigma_max * self.geometry.xi
        dt    = self._dt()

        u = np.empty(self.n_steps, dtype=float)
        v = np.empty(self.n_steps, dtype=float)
        E = np.empty(self.n_steps, dtype=float)
        hits: List[int] = []

        # Initial conditions — start at amplitude_frac·A_max with a small
        # outward velocity so wall reflections happen (and to seed energy
        # above the static cap by an infinitesimal amount when frac→1).
        u[0] = self.amplitude_frac * A_max
        # Velocity that pushes amplitude just barely over the cap so the
        # first oscillation actually reflects rather than falling short
        # by a Verlet round-off margin.
        v[0] = 0.05 * A_max * omega
        E[0] = 0.5 * self.geometry.rho * v[0]**2 + 0.5 * K_eff * u[0]**2

        # Velocity-Verlet for ρ ü + γ u̇ + K_eff u = 0  →
        #     ü = -(K_eff/ρ) u  - 2Γ u̇
        for i in range(self.n_steps - 1):
            a = -(K_eff / self.geometry.rho) * u[i] - 2.0 * Gamma * v[i]
            v_half = v[i] + 0.5 * a * dt
            u_new = u[i] + v_half * dt

            # Wall reflection on cap exceedance
            if u_new > A_max:
                u_new = A_max - (u_new - A_max)
                v_half = -abs(v_half)
                hits.append(i + 1)
            elif u_new < -A_max:
                u_new = -A_max + (-A_max - u_new)
                v_half = abs(v_half)
                hits.append(i + 1)

            a_new = -(K_eff / self.geometry.rho) * u_new - 2.0 * Gamma * v_half
            v_new = v_half + 0.5 * a_new * dt

            u[i + 1] = u_new
            v[i + 1] = v_new
            E[i + 1] = 0.5 * self.geometry.rho * v_new**2 + 0.5 * K_eff * u_new**2

        t = np.arange(self.n_steps) * dt
        return {
            "t":      t,
            "u":      u,
            "v":      v,
            "energy": E,
            "hits":   np.asarray(hits, dtype=int),
            "A_max":  A_max,
        }

    # ------------------------------------------------------------------
    # Period extraction (zero-crossing method)
    # ------------------------------------------------------------------

    def measured_period(self, result: Optional[Dict[str, np.ndarray]] = None) -> float:
        """Estimate the oscillation period from successive zero crossings.

        Returns the average period [s].  If fewer than two zero crossings
        are found, returns ``float('nan')``.
        """
        if result is None:
            result = self.run()
        u = result["u"]
        t = result["t"]

        # Sign changes (zero crossings)
        signs = np.sign(u)
        zc_idx = np.where((signs[:-1] != signs[1:]) & (signs[:-1] != 0))[0]
        if len(zc_idx) < 3:
            return float("nan")

        # Linear interpolation for sub-step crossings (improves accuracy)
        zc_times = []
        for i in zc_idx:
            t0, t1 = t[i], t[i + 1]
            u0, u1 = u[i], u[i + 1]
            if u1 != u0:
                t_cross = t0 - u0 * (t1 - t0) / (u1 - u0)
            else:
                t_cross = t0
            zc_times.append(t_cross)
        zc_times = np.asarray(zc_times)

        # Period = 2 × mean spacing between consecutive zero crossings
        spacings = np.diff(zc_times)
        return float(2.0 * np.mean(spacings))

    # ------------------------------------------------------------------
    # Animation snapshots — return frames for plotting
    # ------------------------------------------------------------------

    def animate_bouncing(
        self,
        n_snapshots: int = 6,
    ) -> List[Dict[str, np.ndarray]]:
        """Produce ``n_snapshots`` time-evolution snapshots of the envelope.

        Each snapshot dict has keys ``t``, ``u``, ``v``, ``energy``, plus
        the index ``snapshot_index``.  Snapshots are evenly spaced across
        the full integration window.
        """
        full = self.run()
        N    = self.n_steps
        idx  = np.linspace(0, N - 1, n_snapshots, dtype=int)
        out: List[Dict[str, np.ndarray]] = []
        for k, i in enumerate(idx):
            out.append({
                "snapshot_index": k,
                "t":      float(full["t"][i]),
                "u":      float(full["u"][i]),
                "v":      float(full["v"][i]),
                "energy": float(full["energy"][i]),
            })
        return out

    # ------------------------------------------------------------------
    # Drag scan: vary γ across lepton anchors
    # ------------------------------------------------------------------

    def drag_scan(
        self,
        gammas: Optional[np.ndarray] = None,
    ) -> Dict[str, np.ndarray]:
        """Sweep drag coefficient over a 5-point ladder anchored to leptons.

        If ``gammas`` is not provided, a 5-point logarithmic sweep from
        γ = 0 (electron, drag-free Compton anchor) up to a value that
        produces ω_b ≈ 100 × m_e c² / ℏ is used.  We then report the
        canonical ω_b and the inferred mass.

        Returns dict with arrays ``gamma``, ``omega_b``, ``mass_kg``,
        ``mass_mev``, plus ``period_canonical`` and
        ``period_measured`` from the simulator.
        """
        if gammas is None:
            # 5-point ladder, geometric, plus the γ=0 electron point.
            # Stay UNDERDAMPED so the time-domain response actually
            # oscillates and the simulator can measure a period; in the
            # overdamped (Γ ≥ ω₀) regime ω_b from the canonical formula
            # is still the resonance-peak frequency but the impulse
            # response decays monotonically (no zero crossings).  We
            # cap Γ at 0.5 × ω₀ to stay safely underdamped.
            base = self.geometry
            omega0 = float(np.sqrt(base.K / base.rho) / base.xi)
            gamma_max = 0.5 * omega0 * 2.0 * base.rho
            gammas = np.concatenate(([0.0], np.geomspace(gamma_max / 1e3,
                                                          gamma_max, 4)))

        gammas = np.asarray(gammas, dtype=float)
        omegas = np.empty_like(gammas)
        masses_kg = np.empty_like(gammas)
        period_can = np.empty_like(gammas)
        period_meas = np.empty_like(gammas)

        for i, g in enumerate(gammas):
            geom = ConeBouncingGeometry(
                K=self.geometry.K,
                rho=self.geometry.rho,
                xi=self.geometry.xi,
                gamma=float(g),
                sigma_max=self.geometry.sigma_max,
            )
            omegas[i] = geom.omega_b
            masses_kg[i] = geom.mass_kg
            period_can[i] = geom.period
            sim = ConeBouncingSimulator(
                geometry=geom,
                amplitude_frac=self.amplitude_frac,
                n_steps=self.n_steps,
                dt_per_period=self.dt_per_period,
            )
            period_meas[i] = sim.measured_period()

        return {
            "gamma":            gammas,
            "omega_b":          omegas,
            "mass_kg":          masses_kg,
            "mass_mev":         masses_kg * C * C / MEV,
            "period_canonical": period_can,
            "period_measured":  period_meas,
        }


# ---------------------------------------------------------------------------
# 3. Visualisation helpers
# ---------------------------------------------------------------------------

def make_bouncing_figure(
    geometry: Optional[ConeBouncingGeometry] = None,
    simulator: Optional[ConeBouncingSimulator] = None,
):
    """Top panel: cone walls + envelope.  Bottom panel: u(t).

    Returns the matplotlib Figure (not yet saved).
    """
    import matplotlib.pyplot as plt

    if geometry is None:
        geometry = ConeBouncingGeometry()
    if simulator is None:
        simulator = ConeBouncingSimulator(geometry=geometry)

    result = simulator.run()
    t = result["t"]
    u = result["u"]
    A_max = result["A_max"]

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(11, 8),
                                          gridspec_kw={"height_ratios": [1.2, 1.0]})

    # ----- Top: cone walls and bouncing trajectory -------------------
    # Show the slanted causal cone (light-cone tilted by saturation),
    # plus horizontal saturation caps.
    t_geom = np.linspace(t.min(), t.max(), 200)
    upper_slant, lower_slant = geometry.slanted_walls(t_geom - t_geom.mean())
    upper_cap,   lower_cap   = geometry.cone_walls(t_geom)

    ax_top.fill_between(t_geom, upper_cap, lower_cap, color="#fff7e6",
                        alpha=0.45, label=f"interior |σ| < {geometry.sigma_max}")
    ax_top.plot(t_geom, upper_cap, color="orange", linestyle="--", linewidth=2,
                label=f"σ = +{geometry.sigma_max} (cap)")
    ax_top.plot(t_geom, lower_cap, color="orange", linestyle="--", linewidth=2,
                label=f"σ = −{geometry.sigma_max} (cap)")
    ax_top.plot(t_geom, upper_slant, color="grey", linestyle=":", linewidth=1,
                alpha=0.7, label="tilted causal cone")
    ax_top.plot(t_geom, lower_slant, color="grey", linestyle=":", linewidth=1,
                alpha=0.7)
    ax_top.plot(t, u, "b-", linewidth=1.6, label="bouncing envelope u(t)")

    # Mark wall reflections
    hits = result["hits"]
    if len(hits) > 0:
        ax_top.scatter(t[hits], u[hits], c="red", s=20, zorder=5,
                       label=f"reflections (n={len(hits)})")

    ax_top.set_xlabel("Time [s]")
    ax_top.set_ylabel("Strain envelope u [m]")
    ax_top.set_title(
        f"Cone-bouncing geometry (cone half-angle {geometry.cone_half_angle_deg:.2f}°)\n"
        f"K={geometry.K:.2e} J/m³, ρ={geometry.rho:.2e} kg/m³, "
        f"ξ={geometry.xi:.2e} m, γ={geometry.gamma:.2e} kg/(m³·s)"
    )
    ax_top.set_ylim(-1.4 * A_max, 1.4 * A_max)
    ax_top.legend(fontsize=8, loc="upper right", ncol=2)
    ax_top.grid(True, alpha=0.3)

    # ----- Bottom: oscillation u(t) with measured/analytic period ----
    period_an = geometry.period
    period_me = simulator.measured_period(result)
    ax_bot.plot(t, u, "b-", linewidth=1.4)
    ax_bot.axhline(0.0, color="black", linewidth=0.7, alpha=0.6)
    ax_bot.axhline(A_max, color="orange", linestyle="--", linewidth=1, alpha=0.6)
    ax_bot.axhline(-A_max, color="orange", linestyle="--", linewidth=1, alpha=0.6)

    # Mark one analytic period for visual reference
    if np.isfinite(period_an) and period_an < (t.max() - t.min()):
        t0 = t[0]
        ax_bot.axvspan(t0, t0 + period_an, color="green", alpha=0.10,
                       label=f"T_analytic = 2π/ω_b = {period_an:.3e} s")
        if np.isfinite(period_me):
            ax_bot.axvspan(t0, t0 + period_me, color="red", alpha=0.05,
                           label=f"T_measured = {period_me:.3e} s")

    ax_bot.set_xlabel("Time [s]")
    ax_bot.set_ylabel("u(t) [m]")
    rel_err = (abs(period_me - period_an) / period_an) if (np.isfinite(period_me) and period_an > 0) else float("nan")
    ax_bot.set_title(
        f"Oscillation u(t):  m c² = ℏ ω_b = {geometry.mass_mev:.3e} MeV  "
        f"(period match: Δ={rel_err*100:.2f}%)"
    )
    ax_bot.legend(fontsize=8, loc="upper right")
    ax_bot.grid(True, alpha=0.3)

    fig.tight_layout()
    return fig


def make_drag_scan_figure(
    simulator: Optional[ConeBouncingSimulator] = None,
    gammas: Optional[np.ndarray] = None,
):
    """Drag-scan figure: m vs γ on log-log scale.

    Top panel: ω_b(γ) and m(γ) on log-log; bottom panel: predicted vs
    measured period agreement.
    Returns the matplotlib Figure.
    """
    import matplotlib.pyplot as plt

    if simulator is None:
        simulator = ConeBouncingSimulator()

    scan = simulator.drag_scan(gammas=gammas)
    g    = scan["gamma"]
    omg  = scan["omega_b"]
    m_kg = scan["mass_kg"]
    p_an = scan["period_canonical"]
    p_me = scan["period_measured"]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8),
                              gridspec_kw={"height_ratios": [1.0, 0.9]})

    # ----- Top: m vs γ (log-log) ------------------------------------
    ax = axes[0]
    # Replace γ=0 with a tiny positive value for log-axis plotting
    g_plot = g.copy()
    g_plot[g_plot <= 0] = max(g[g > 0].min() / 10.0, 1e-30) if np.any(g > 0) else 1e-30

    ax.loglog(g_plot, m_kg, "bo-", linewidth=2, markersize=8,
              label="m = ℏω_b/c²  (canonical)")
    # Reference horizontal lines for m_e, m_μ, m_τ
    ax.axhline(M_E,  color="green",  linestyle="--", linewidth=1, alpha=0.6,
               label=f"m_e  = {M_E:.2e} kg")
    ax.axhline(M_MU, color="orange", linestyle="--", linewidth=1, alpha=0.6,
               label=f"m_μ  = {M_MU:.2e} kg")
    ax.axhline(M_TAU, color="red",   linestyle="--", linewidth=1, alpha=0.6,
               label=f"m_τ  = {M_TAU:.2e} kg")
    ax.set_xlabel("Drag coefficient γ [kg/(m³·s)]")
    ax.set_ylabel("Rest mass m [kg]")
    ax.set_title("Cone-bouncing mass ladder: m(γ) at fixed (K, ρ, ξ = ξ_e)")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3, which="both")

    # ----- Bottom: period match ------------------------------------
    ax = axes[1]
    rel_err = np.where(p_an > 0, np.abs(p_me - p_an) / p_an, np.nan)
    ax.semilogx(g_plot, rel_err * 100.0, "rs-", linewidth=2, markersize=8,
                label="|T_measured − T_analytic| / T_analytic")
    ax.axhline(5.0, color="black", linestyle="--", linewidth=1, alpha=0.6,
               label="5% acceptance band")
    ax.set_xlabel("Drag coefficient γ [kg/(m³·s)]")
    ax.set_ylabel("Relative period error [%]")
    ax.set_title("Simulator vs canonical: period agreement")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Convenience: reproduce a 5-point lepton-anchor scan
# ---------------------------------------------------------------------------

def lepton_anchor_drag_scan() -> Dict[str, np.ndarray]:
    """Return drag values + frequencies / masses for the 5-point ladder
    anchored at (electron, muon, tau)-like Compton scales.

    Strategy: at fixed (K, ρ) = electron-anchor primitives, vary ξ across
    {ξ_e, ξ_μ, ξ_τ} and γ = 0 — these give the three SM lepton masses
    automatically by Compton baseline.  Two extra γ-driven points fill
    the ladder.  This is a *visualisation* aid — not a derivation.
    """
    p = electron_anchor_primitives()
    K, rho = p["K"], p["rho"]
    rows = []
    for label, xi in [("e (ξ_e, γ=0)", XI_E),
                      ("μ (ξ_μ, γ=0)", XI_MU),
                      ("τ (ξ_τ, γ=0)", XI_TAU)]:
        omega = omega_b_canonical(K, rho, xi, 0.0)
        rows.append((label, 0.0, xi, omega, mass_from_omega(omega)))
    # Two drag-driven points at ξ_e
    omega0 = float(np.sqrt(K / rho) / XI_E)
    for label, gfac in [("e + drag×10", 10.0), ("e + drag×100", 100.0)]:
        gamma = gfac * omega0 * 2.0 * rho
        omega = omega_b_canonical(K, rho, XI_E, gamma)
        rows.append((label, gamma, XI_E, omega, mass_from_omega(omega)))

    return {
        "label":   np.array([r[0] for r in rows]),
        "gamma":   np.array([r[1] for r in rows]),
        "xi":      np.array([r[2] for r in rows]),
        "omega_b": np.array([r[3] for r in rows]),
        "mass_kg": np.array([r[4] for r in rows]),
    }


if __name__ == "__main__":
    geom = ConeBouncingGeometry()
    sim  = ConeBouncingSimulator(geometry=geom)
    print(f"Cone half-angle: {geom.cone_half_angle_deg:.2f}°")
    print(f"ω_b (electron):  {geom.omega_b:.4e} rad/s")
    print(f"Period:          {geom.period:.4e} s")
    print(f"Mass:            {geom.mass_mev:.4e} MeV")
    res = sim.run()
    print(f"Measured period: {sim.measured_period(res):.4e} s")
