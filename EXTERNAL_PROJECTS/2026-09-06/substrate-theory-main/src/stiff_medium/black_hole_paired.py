"""Paired GEOMETRY + SIM + VIZ for black holes through substrate σ → 1/2 saturation.

This module provides three coupled classes for the B3 / stiff-medium picture
of black holes as substrate-saturation regions (σ_max = 1/2):

  1. ``BlackHoleGeometry`` — static geometry: Schwarzschild metric near the
     horizon and the substrate strain profile σ(r) = GM/(r c²) with σ = 1/2
     at r = r_s = 2GM/c².

  2. ``BlackHoleSimulator`` — dynamics: time-evolves a photon trajectory
     through the substrate gradient with c_local(r) = c · √(1 − 2σ(r)),
     demonstrating horizon formation and capture; emits Hawking temperature
     T_H = ℏ c³ / (8π G M k_B) and Schwarzschild ringdown frequency
     f ≈ 0.0594 c³ / (G M_f) for the (2,2,0) QNM.

  3. ``BlackHoleVisualizer`` — renders the substrate σ(r) profile with
     photon trajectories near r_s (panel 34) and the Hawking blackbody
     spectrum + ringdown waveform (panel 35).

In the substrate model the horizon is *defined* by σ = 1/2 (the universal
saturation cap). There is no singularity: the interior is uniformly
saturated. Hawking radiation is thermal noise leaking past the saturated
boundary; ringdown is the (2,2,0) quasinormal oscillation of the σ = 1/2
shell. Both are computed from the standard Schwarzschild expressions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018 / IAU 2015)
# ---------------------------------------------------------------------------

C: Final[float] = 2.99792458e8           # speed of light, m/s
G: Final[float] = 6.67430e-11            # gravitational constant
HBAR: Final[float] = 1.054571817e-34     # reduced Planck constant
K_B: Final[float] = 1.380649e-23         # Boltzmann constant
M_SUN: Final[float] = 1.98892e30         # solar mass, kg
PLANCK_LENGTH: Final[float] = math.sqrt(HBAR * G / C**3)

#: Universal substrate saturation cap (B3 §18.39): horizon ↔ σ = 1/2.
SIGMA_HORIZON: Final[float] = 0.5


# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------


@dataclass
class BlackHoleGeometry:
    """Schwarzschild metric and substrate strain profile near the horizon.

    The geometry is parameterised by mass ``M`` (kg) only (Schwarzschild,
    spin-0). The substrate strain σ(r) = GM/(r c²) connects the metric
    to the saturation field; σ = 1/2 at r = r_s defines the horizon.
    """

    M: float                                       # gravitational mass, kg
    sigma_cap: float = SIGMA_HORIZON               # universal cap = 1/2

    # ------------------------------------------------------------------
    # Length / area scales
    # ------------------------------------------------------------------
    @property
    def r_s(self) -> float:
        """Schwarzschild radius r_s = 2GM/c² (m)."""
        return 2.0 * G * self.M / (C * C)

    @property
    def r_g(self) -> float:
        """Gravitational radius GM/c² = r_s/2 (m)."""
        return G * self.M / (C * C)

    @property
    def r_photon(self) -> float:
        """Photon-sphere radius r = 1.5 r_s (Schwarzschild)."""
        return 1.5 * self.r_s

    @property
    def r_isco(self) -> float:
        """Innermost stable circular orbit r_ISCO = 3 r_s (Schwarzschild)."""
        return 3.0 * self.r_s

    @property
    def horizon_area(self) -> float:
        """Horizon surface area A = 4π r_s² (m²)."""
        return 4.0 * math.pi * self.r_s ** 2

    # ------------------------------------------------------------------
    # σ(r) — substrate strain profile
    # ------------------------------------------------------------------
    def sigma(self, r: float | np.ndarray) -> float | np.ndarray:
        """Substrate strain σ(r) = GM/(r c²), capped at σ_cap inside r_s.

        At r = r_s this gives exactly σ = 1/2 (universal saturation cap).
        Outside the horizon σ falls as 1/r; inside, the medium is uniformly
        saturated so we clamp at σ_cap.
        """
        r_arr = np.asarray(r, dtype=float)
        if np.any(r_arr <= 0):
            raise ValueError("r must be positive")
        sig = self.r_g / r_arr
        sig = np.minimum(sig, self.sigma_cap)
        return sig if r_arr.ndim else float(sig)

    # ------------------------------------------------------------------
    # c_local(r) — wave-speed map from saturation
    # ------------------------------------------------------------------
    def c_local(self, r: float | np.ndarray) -> float | np.ndarray:
        """Effective local light speed c_local(r) = c · √(1 − (σ/σ_cap)²).

        Vanishes at the horizon (σ = σ_cap) and equals c far away.
        """
        sig = self.sigma(r)
        ratio_sq = (np.asarray(sig) / self.sigma_cap) ** 2
        ratio_sq = np.clip(ratio_sq, 0.0, 1.0)
        cl = C * np.sqrt(1.0 - ratio_sq)
        return cl if np.asarray(r).ndim else float(cl)

    # ------------------------------------------------------------------
    # Schwarzschild metric components (static, spherical)
    # ------------------------------------------------------------------
    def metric_g_tt(self, r: float | np.ndarray) -> float | np.ndarray:
        """Time-time component g_tt = −(1 − r_s/r) c²."""
        r_arr = np.asarray(r, dtype=float)
        return -(1.0 - self.r_s / r_arr) * C * C

    def metric_g_rr(self, r: float | np.ndarray) -> float | np.ndarray:
        """Radial-radial component g_rr = 1 / (1 − r_s/r)  (singular at r_s)."""
        r_arr = np.asarray(r, dtype=float)
        denom = 1.0 - self.r_s / r_arr
        # Avoid hard 1/0 at the horizon — clamp tiny values
        denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
        return 1.0 / denom

    def gravitational_redshift(self, r_emit: float, r_obs: float) -> float:
        """Frequency ratio f_obs / f_emit between two static observers.

            f_obs / f_emit = √((1 − r_s/r_emit) / (1 − r_s/r_obs))

        Diverges to 0 as r_emit → r_s (infinite redshift at the horizon).
        """
        if r_emit <= self.r_s:
            return 0.0
        return math.sqrt((1.0 - self.r_s / r_emit) /
                         (1.0 - self.r_s / r_obs))


# ---------------------------------------------------------------------------
# 2. Simulator
# ---------------------------------------------------------------------------


@dataclass
class PhotonTrajectory:
    """Result of a photon time-evolution through the substrate gradient."""

    t: np.ndarray       # times (s)
    r: np.ndarray       # radii (m)
    captured: bool      # True if the photon crossed the horizon


class BlackHoleSimulator:
    """Time-evolution of photon paths through the substrate strain field.

    The simulator wraps a ``BlackHoleGeometry`` and provides:

      * ``evolve_radial_photon`` — trace an inward-falling photon through
        c_local(r), demonstrating horizon capture (r → r_s, c_local → 0).
      * ``hawking_temperature`` — T_H = ℏ c³ / (8π G M k_B) (Schwarzschild).
      * ``hawking_spectrum`` — Planck blackbody at T_H, returns frequency
        and spectral radiance arrays.
      * ``ringdown_frequency`` / ``ringdown_waveform`` — (2,2,0) QNM:
        ω r_g/c = 0.3737 − 0.0890 i  (Leaver 1985), giving
        f ≈ 0.0594 c³/(GM) and τ ≈ 11.24 GM/c³ for a Schwarzschild remnant.
    """

    def __init__(self, geometry: BlackHoleGeometry) -> None:
        self.geom = geometry

    # ------------------------------------------------------------------
    # Radial photon trajectory through c_local(r)
    # ------------------------------------------------------------------
    def evolve_radial_photon(self,
                             r0: float,
                             t_max: float | None = None,
                             n_steps: int = 4000,
                             stop_factor: float = 1.001
                             ) -> PhotonTrajectory:
        """Evolve an inward-falling radial photon through the substrate gradient.

        Uses dr/dt = −c_local(r). The photon decelerates as σ(r) → 1/2,
        asymptotically approaching r_s (the σ = 1/2 horizon) but never
        crossing it in coordinate time — the universal Schwarzschild
        result. We declare the photon "captured" once r ≤ stop_factor · r_s
        (default ~0.1% above horizon).

        Parameters
        ----------
        r0 : starting radius (m), must satisfy r0 > r_s.
        t_max : evolution time (s); defaults to 50 r_s/c.
        n_steps : number of time steps (uniform grid).
        stop_factor : capture detection threshold in units of r_s.

        Returns a ``PhotonTrajectory`` containing t, r arrays and a captured flag.
        """
        rs = self.geom.r_s
        if r0 <= rs:
            raise ValueError(f"r0={r0:.3e} must exceed r_s={rs:.3e}")
        if t_max is None:
            t_max = 50.0 * rs / C
        dt = t_max / n_steps
        r = float(r0)
        ts: list[float] = [0.0]
        rs_track: list[float] = [r]
        captured = False
        for k in range(1, n_steps + 1):
            cl = self.geom.c_local(r)
            r -= float(cl) * dt
            t = k * dt
            if r <= stop_factor * rs:
                ts.append(t)
                rs_track.append(max(r, stop_factor * rs))
                captured = True
                break
            ts.append(t)
            rs_track.append(r)
        return PhotonTrajectory(t=np.asarray(ts),
                                r=np.asarray(rs_track),
                                captured=captured)

    # ------------------------------------------------------------------
    # Hawking temperature and spectrum
    # ------------------------------------------------------------------
    def hawking_temperature(self) -> float:
        """T_H = ℏ c³ / (8π G M k_B)  (K).

        For M = M_sun this gives T_H ≈ 6.17 × 10⁻⁸ K.
        """
        return HBAR * C**3 / (8.0 * math.pi * G * self.geom.M * K_B)

    def hawking_spectrum(self,
                         f_min: float | None = None,
                         f_max: float | None = None,
                         n: int = 400) -> tuple[np.ndarray, np.ndarray]:
        """Planck spectral radiance B_ν(T_H) over a frequency window.

        B_ν(T) = (2 h ν³ / c²) · 1 / (exp(hν/kT) − 1)

        For a Schwarzschild BH the Hawking emission is approximately
        thermal (ignoring grey-body factors). We return (ν, B_ν) on a
        log-spaced grid centred on the spectral peak ν_peak ≈ 5.88·k_B T/h.
        """
        T = self.hawking_temperature()
        h = 2.0 * math.pi * HBAR
        nu_peak = 5.879e10 * T              # Wien displacement (Hz/K)
        if f_min is None:
            f_min = 1e-3 * nu_peak
        if f_max is None:
            f_max = 1e3 * nu_peak
        nu = np.logspace(math.log10(f_min), math.log10(f_max), n)
        x = h * nu / (K_B * T)
        # Avoid overflow in exp for very large x
        x_clip = np.minimum(x, 700.0)
        B_nu = (2.0 * h * nu**3 / C**2) / (np.expm1(x_clip))
        return nu, B_nu

    # ------------------------------------------------------------------
    # Ringdown — (2,2,0) Schwarzschild QNM
    # ------------------------------------------------------------------
    def ringdown_frequency(self) -> float:
        """Real ringdown frequency f (Hz) of the (2,2,0) Schwarzschild QNM.

        ω · GM/c³ = 0.3737 − 0.0890 i  →  f = 0.3737 c³/(2π G M)
                                               ≈ 0.0594 c³/(G M)

        For ~62 M_sun (GW150914 remnant) this gives f ≈ 250 Hz.
        For a single 10 M_sun BH this gives f ≈ 1.2 kHz.
        """
        omega_re = 0.3737 * C**3 / (G * self.geom.M)
        return omega_re / (2.0 * math.pi)

    def ringdown_damping_time(self) -> float:
        """Damping time τ = 1/|Im ω| ≈ 11.24 GM/c³  (s).

        For ~62 M_sun this is τ ≈ 3.4 ms.
        """
        omega_im = 0.0890 * C**3 / (G * self.geom.M)
        return 1.0 / omega_im

    def ringdown_quality(self) -> float:
        """Quality factor Q = Re ω / (2 Im ω) ≈ 2.10 (universal for 2,2,0)."""
        return 0.5 * 0.3737 / 0.0890

    def ringdown_waveform(self,
                          n_cycles: float = 8.0,
                          n_samples: int = 2000
                          ) -> tuple[np.ndarray, np.ndarray]:
        """Time-domain ringdown h(t) = exp(−t/τ) · cos(2π f t).

        Returns (t, h). Length covers ``n_cycles`` periods.
        """
        f = self.ringdown_frequency()
        tau = self.ringdown_damping_time()
        T_total = n_cycles / f
        t = np.linspace(0.0, T_total, n_samples)
        h = np.exp(-t / tau) * np.cos(2.0 * math.pi * f * t)
        return t, h


# ---------------------------------------------------------------------------
# 3. Visualizer
# ---------------------------------------------------------------------------


class BlackHoleVisualizer:
    """Render-side helpers for the paired geometry + simulator.

    Used by ``scripts/render_all_visuals.render_black_hole`` to produce
    the two PNG panels.
    """

    def __init__(self, sim: BlackHoleSimulator) -> None:
        self.sim = sim
        self.geom = sim.geom

    def horizon_panel(self, ax_sigma, ax_traj) -> None:
        """Plot σ(r) profile and a few inward photon trajectories."""
        rs = self.geom.r_s
        # σ(r) over [0.3 r_s, 6 r_s]
        r_grid = np.linspace(0.3 * rs, 6.0 * rs, 600)
        sig = self.geom.sigma(r_grid)
        ax_sigma.plot(r_grid / rs, sig, "b-", linewidth=2.2,
                      label=r"$\sigma(r) = GM/(r c^2)$")
        ax_sigma.axhline(SIGMA_HORIZON, color="red", linestyle="--",
                         linewidth=1.4,
                         label=r"$\sigma_{\max} = 1/2$ (saturation cap)")
        ax_sigma.axvline(1.0, color="black", linestyle=":", linewidth=1.2,
                         label=r"$r_s$ (horizon)")
        ax_sigma.axvline(1.5, color="green", linestyle=":", linewidth=1.0,
                         label=r"$1.5\,r_s$ (photon sphere)")
        ax_sigma.set_xlabel(r"radius $r / r_s$")
        ax_sigma.set_ylabel(r"substrate strain $\sigma$")
        ax_sigma.set_title("Substrate strain near the horizon")
        ax_sigma.set_ylim(0.0, 0.6)
        ax_sigma.set_xlim(r_grid.min() / rs, r_grid.max() / rs)
        ax_sigma.grid(True, alpha=0.3)
        ax_sigma.legend(loc="upper right", fontsize=9)

        # A few photon trajectories
        starts_rs = [5.0, 3.0, 2.0, 1.4, 1.1]
        colors = ["#1b3b6f", "#2a6f97", "#4cae4c", "#e07a1f", "#c0392b"]
        for r0_factor, color in zip(starts_rs, colors):
            r0 = r0_factor * rs
            traj = self.sim.evolve_radial_photon(r0=r0, n_steps=4000)
            t_us = traj.t * 1e6 if rs < 1e5 else traj.t * 1e3
            ax_traj.plot(t_us, traj.r / rs, color=color, linewidth=1.6,
                         label=f"r₀ = {r0_factor:.1f} r_s "
                               f"({'captured' if traj.captured else 'escaping'})")
        ax_traj.axhline(1.0, color="red", linestyle="--", linewidth=1.4,
                        label="horizon r_s")
        ax_traj.axhline(1.5, color="green", linestyle=":", linewidth=1.0,
                        label="photon sphere 1.5 r_s")
        time_unit = r"time [µs]" if rs < 1e5 else r"time [ms]"
        ax_traj.set_xlabel(time_unit)
        ax_traj.set_ylabel(r"radius $r / r_s$")
        ax_traj.set_title(r"Inward photon trajectories: $dr/dt = -c_{\rm local}(r)$")
        ax_traj.grid(True, alpha=0.3)
        ax_traj.legend(loc="upper right", fontsize=8)

    def hawking_ringdown_panel(self, ax_spec, ax_wave) -> None:
        """Plot Hawking blackbody spectrum (left) and ringdown waveform (right)."""
        T_H = self.sim.hawking_temperature()
        nu, B_nu = self.sim.hawking_spectrum(n=500)
        ax_spec.loglog(nu, B_nu, "b-", linewidth=2.0)
        ax_spec.set_xlabel(r"frequency $\nu$ [Hz]")
        ax_spec.set_ylabel(r"$B_\nu(T_H)$ [W·sr$^{-1}$·m$^{-2}$·Hz$^{-1}$]")
        ax_spec.set_title(
            f"Hawking spectrum  M = {self.geom.M / M_SUN:.1f} M$_\\odot$\n"
            f"$T_H = {T_H:.3e}$ K"
        )
        ax_spec.grid(True, which="both", alpha=0.3)
        # Mark Wien peak
        nu_peak = 5.879e10 * T_H
        ax_spec.axvline(nu_peak, color="red", linestyle="--", linewidth=1.0,
                        label=f"Wien peak {nu_peak:.2e} Hz")
        ax_spec.legend(loc="lower center", fontsize=9)

        # Ringdown waveform
        f_qnm = self.sim.ringdown_frequency()
        tau = self.sim.ringdown_damping_time()
        t, h = self.sim.ringdown_waveform(n_cycles=8.0)
        scale = 1e3 if t.max() < 1.0 else 1.0
        unit = "ms" if scale == 1e3 else "s"
        ax_wave.plot(t * scale, h, "k-", linewidth=1.4)
        envelope = np.exp(-t / tau)
        ax_wave.plot(t * scale, envelope, "r--", linewidth=1.0,
                     label=r"envelope $e^{-t/\tau}$")
        ax_wave.plot(t * scale, -envelope, "r--", linewidth=1.0)
        ax_wave.axhline(0.0, color="gray", linewidth=0.6)
        ax_wave.set_xlabel(f"time [{unit}]")
        ax_wave.set_ylabel(r"strain $h(t)$ (arb. units)")
        ax_wave.set_title(
            f"Ringdown (2,2,0) QNM   $f = {f_qnm:.1f}$ Hz, "
            f"$\\tau = {tau * scale:.2f}$ {unit},  Q = {self.sim.ringdown_quality():.2f}"
        )
        ax_wave.grid(True, alpha=0.3)
        ax_wave.legend(loc="upper right", fontsize=9)


# ---------------------------------------------------------------------------
# Module-level convenience
# ---------------------------------------------------------------------------


def default_solar_bh() -> BlackHoleSimulator:
    """Convenience: solar-mass Schwarzschild BH simulator."""
    return BlackHoleSimulator(BlackHoleGeometry(M=M_SUN))


def default_stellar_remnant() -> BlackHoleSimulator:
    """Convenience: ~48 M_sun stellar BH (Schwarzschild ringdown ≈ 250 Hz).

    Note: GW150914's 62 M_sun remnant rang at ~250 Hz because of its
    Kerr spin (a ≈ 0.69). For a non-spinning Schwarzschild BH the same
    250 Hz frequency corresponds to M ≈ 48 M_sun via f = 0.0594 c³/(GM).
    """
    return BlackHoleSimulator(BlackHoleGeometry(M=48.0 * M_SUN))


__all__ = [
    "BlackHoleGeometry",
    "BlackHoleSimulator",
    "BlackHoleVisualizer",
    "PhotonTrajectory",
    "default_solar_bh",
    "default_stellar_remnant",
    "C", "G", "HBAR", "K_B", "M_SUN", "PLANCK_LENGTH", "SIGMA_HORIZON",
]


if __name__ == "__main__":  # pragma: no cover
    sim_sun = default_solar_bh()
    sim_sr = default_stellar_remnant()
    print(f"Solar-mass BH:")
    print(f"  r_s = {sim_sun.geom.r_s:.3e} m")
    print(f"  σ(r_s) = {sim_sun.geom.sigma(sim_sun.geom.r_s):.4f}")
    print(f"  T_H = {sim_sun.hawking_temperature():.3e} K")
    print(f"  ringdown f = {sim_sun.ringdown_frequency():.1f} Hz")
    print(f"Stellar remnant (62 M_sun):")
    print(f"  ringdown f = {sim_sr.ringdown_frequency():.1f} Hz")
    print(f"  τ = {sim_sr.ringdown_damping_time() * 1e3:.2f} ms")
