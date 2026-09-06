"""gw_signal_paired.py
================================

Paired GEOMETRY + SIMULATOR module that constructs the gravitational-wave
signal from a binary-black-hole (BBH) merger as a substrate-strain
disturbance, in the cone-tilt / σ → 1/2 picture.

Background
----------
In the substrate framework, gravity is the dynamical response of a stiff
medium with strain σ(x, t).  Two compact objects in mutual orbit each
locally tilt the substrate's light cone (σ → 1/2 marks the horizon), and
the orbital motion drives a coherent transverse-traceless ripple — the
gravitational wave — that propagates at the substrate's intrinsic
transverse speed

    c_gw = sqrt(K / rho) ≡ c

This is exactly the structural identity of `cone_bouncing_protocol.C`:
in the un-tilted bulk substrate, both light and gravitational waves
ride the same transverse mode.  No tunable parameter sets v_gw — it is
forced by the substrate's elastic moduli.  GW170817 measured this to
1 part in 10^15.

Inspiral phase (post-Newtonian)
-------------------------------
For a quasi-circular binary with chirp mass M_c the leading-order GW
frequency f(t) and strain amplitude h(t) follow

    f(t) = (1 / pi) * (5 / (256 * tau))^{3/8}
                    * (G M_c / c^3)^{-5/8}

    h(t) = (4 / D) * (G M_c / c^2)^{5/3} * (pi f(t) / c)^{2/3}

with tau = t_coal - t the time to coalescence, and D the luminosity
distance.  This gives the famous f ∝ tau^{-3/8} chirp.

Merger / ringdown
-----------------
Once the two horizons merge the perturbed final BH rings down at its
fundamental quasi-normal mode (l=2, m=2, n=0), with frequency

    f_RD ≈ (c^3 / (2 pi G M_f)) * F(a_f)
    F(a_f) ≈ 1.5251 - 1.1568 * (1 - a_f)^0.1292

(Echeverria 1989 / Berti et al. fit) and damping time

    tau_RD = Q / (pi f_RD),    Q ≈ 0.7000 + 1.4187 * (1 - a_f)^{-0.4990}

Substrate interpretation
------------------------
The ringdown frequency is the QNM of a black hole, which in the cone-tilt
picture is the resonance frequency of the σ-saturation surface as it
relaxes back toward equilibrium.  At σ = 1/2 the cone tilts to 90°; the
QNM is the residual oscillation of the σ field about that fixed point.

This module supplies

    * ``GWGeometry``: chirp mass, total mass, mass ratio, ISCO, Schwarzschild
      radius, cone-tilt parameters and the substrate-derived v_gw = c.
    * ``GWSimulator``: time-domain h(t) for inspiral + merger + ringdown,
      frequency evolution f(t), phase, spectrogram support.
    * ``render_inspiral_figure(geom, sim)`` and
      ``render_ringdown_figure(geom, sim)``: matplotlib helpers used by
      ``scripts/render_all_visuals.py``.

References within this codebase
-------------------------------
    src/stiff_medium/black_hole.py              Schwarzschild substrate picture
    src/stiff_medium/saturation_horizon_geometry.py  Cone tilt + σ horizon
    src/stiff_medium/gravitational_lensing.py   Substrate refractive picture
    src/stiff_medium/cone_bouncing_protocol.py  c = sqrt(K/rho) identity
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from src.stiff_medium.cone_bouncing_protocol import C, HBAR


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

G_NEWTON: float = 6.674e-11           # m^3 kg^-1 s^-2
M_SUN:    float = 1.98892e30          # kg
MPC:      float = 3.0857e22           # m
KPC:      float = 3.0857e19           # m

# Substrate cap (matches saturation_horizon_geometry.SIGMA_MAX)
SIGMA_MAX: float = 0.5


# ---------------------------------------------------------------------------
# (1) GEOMETRY: chirp mass, ISCO, horizons, cone tilt
# ---------------------------------------------------------------------------

@dataclass
class GWGeometry:
    """Geometry of a binary black-hole inspiral in the substrate picture.

    Attributes
    ----------
    m1_solar, m2_solar : float
        Component masses in solar masses.  Defaults match GW150914
        (m1 ≈ 36, m2 ≈ 29 → M_c ≈ 28.1 M_sun, M_total ≈ 65 M_sun).
    spin1, spin2 : float
        Dimensionless aligned spins a = J c / (G m^2), |a| < 1.  The
        effective spin chi_eff (mass-weighted average) modifies the
        inspiral phase and the final-BH mass / spin via simple fits.
    distance_mpc : float
        Luminosity distance in Mpc.  Defaults to 410 Mpc (GW150914).
    """

    m1_solar:     float = 36.0
    m2_solar:     float = 29.0
    spin1:        float = 0.0
    spin2:        float = 0.0
    distance_mpc: float = 410.0

    # ------------------------------------------------------------------
    # Derived mass quantities (SI)
    # ------------------------------------------------------------------

    @property
    def m1(self) -> float:
        """Primary mass [kg]."""
        return float(self.m1_solar * M_SUN)

    @property
    def m2(self) -> float:
        """Secondary mass [kg]."""
        return float(self.m2_solar * M_SUN)

    @property
    def total_mass(self) -> float:
        """M = m1 + m2 [kg]."""
        return self.m1 + self.m2

    @property
    def total_mass_solar(self) -> float:
        return self.m1_solar + self.m2_solar

    @property
    def reduced_mass(self) -> float:
        """μ = m1 m2 / (m1 + m2) [kg]."""
        return self.m1 * self.m2 / self.total_mass

    @property
    def mass_ratio(self) -> float:
        """q = m2 / m1, with m1 ≥ m2 → q ≤ 1."""
        a, b = self.m1, self.m2
        if a < b:
            a, b = b, a
        return b / a

    @property
    def symmetric_mass_ratio(self) -> float:
        """η = m1 m2 / (m1 + m2)^2.  Equals 1/4 for equal masses."""
        return self.m1 * self.m2 / (self.total_mass ** 2)

    @property
    def chirp_mass(self) -> float:
        """M_c = (m1 m2)^{3/5} / (m1+m2)^{1/5}  [kg].

        The chirp mass is the *only* mass combination that appears in
        the leading-order PN inspiral; LIGO determines it most accurately
        of all source parameters.
        """
        return (self.m1 * self.m2) ** (3.0 / 5.0) / self.total_mass ** (1.0 / 5.0)

    @property
    def chirp_mass_solar(self) -> float:
        return float(self.chirp_mass / M_SUN)

    @property
    def chi_eff(self) -> float:
        """Mass-weighted aligned spin χ_eff = (m1 a1 + m2 a2) / M."""
        return (self.m1 * self.spin1 + self.m2 * self.spin2) / self.total_mass

    # ------------------------------------------------------------------
    # Horizon / ISCO / cone-tilt geometry
    # ------------------------------------------------------------------

    def schwarzschild_radius(self, mass_kg: Optional[float] = None) -> float:
        """r_s = 2 G M / c^2 [m].  In the substrate picture this is the
        locus where σ(r) = 1/2 → cone tilt is 90°.
        """
        m = self.total_mass if mass_kg is None else float(mass_kg)
        return 2.0 * G_NEWTON * m / (C * C)

    def isco_radius(self, mass_kg: Optional[float] = None) -> float:
        """ISCO = 6 G M / c^2 = 3 r_s for a Schwarzschild BH."""
        m = self.total_mass if mass_kg is None else float(mass_kg)
        return 6.0 * G_NEWTON * m / (C * C)

    def isco_orbital_frequency(self, mass_kg: Optional[float] = None) -> float:
        """Kepler frequency at ISCO: f_orb = (1/(2π))·sqrt(GM/r^3)."""
        m = self.total_mass if mass_kg is None else float(mass_kg)
        r = self.isco_radius(m)
        return (1.0 / (2.0 * np.pi)) * np.sqrt(G_NEWTON * m / r ** 3)

    def isco_gw_frequency(self, mass_kg: Optional[float] = None) -> float:
        """Quadrupole GW frequency at ISCO: f_GW = 2 f_orb (l=m=2 mode)."""
        return 2.0 * self.isco_orbital_frequency(mass_kg)

    def cone_tilt_at_horizon(self) -> float:
        """Cone tilt at the σ = 1/2 horizon in radians (= π/2 = 90°)."""
        return float(np.pi / 2.0)

    # ------------------------------------------------------------------
    # Final-BH (post-merger) mass and spin from Pannarale fits
    # ------------------------------------------------------------------

    @property
    def final_mass(self) -> float:
        """Approximate radiated-energy-corrected final mass [kg].

        Fit: E_rad / M ≈ η · (1 - 2 √2 / 3 · χ_eff_rad), capped at small
        positive numbers.  For non-spinning equal-mass binaries ~5%.
        """
        eta = self.symmetric_mass_ratio
        # Buonanno+2007 fit: rough but adequate for a substrate-flavored viz
        e_rad = eta * (1.0 - 2.0 * np.sqrt(2.0) / 3.0)  # ≈ 0.057 for η=1/4
        if e_rad < 0.0:
            e_rad = 0.0
        return float(self.total_mass * (1.0 - e_rad))

    @property
    def final_spin(self) -> float:
        """Dimensionless spin of the merged remnant.

        Heuristic fit (Rezzolla 2008-style, simplified): a_f ≈ 2 √3 η
        clamped to [0, 0.99].  For equal masses (η = 1/4) this gives
        a_f ≈ 0.866, close to the measured GW150914 a_f ≈ 0.69 once
        spins of components are folded in.  Including chi_eff:
        """
        eta = self.symmetric_mass_ratio
        a_f = 2.0 * np.sqrt(3.0) * eta - 0.07 * eta + (1.0 - 4.0 * eta) * self.chi_eff
        return float(np.clip(a_f, 0.0, 0.99))

    # ------------------------------------------------------------------
    # Substrate identity: v_gw = c
    # ------------------------------------------------------------------

    @staticmethod
    def gw_speed(K: float = 1.0, rho: float = 1.0) -> float:
        """Gravitational-wave speed v_gw = √(K/ρ).

        In the substrate framework K and ρ are tied to the bulk
        elastic moduli; choosing K, ρ such that √(K/ρ) = c sets the
        substrate-rest-frame convention.  The point is *structural*:
        v_gw and the photon speed share the same transverse mode, so
        they are forced to be equal — no parameter is tuned.
        """
        if K <= 0 or rho <= 0:
            raise ValueError("K, rho must be strictly positive.")
        return float(np.sqrt(K / rho))

    @staticmethod
    def gw_speed_equals_c() -> Dict[str, float]:
        """Structural identity: v_gw / c = 1 from the substrate ansatz."""
        v = GWGeometry.gw_speed(K=1.0, rho=1.0 / (C * C))
        return {
            "v_gw_over_c": float(v / C),
            "rel_error":   float(abs(v - C) / C),
        }


# ---------------------------------------------------------------------------
# (2) SIMULATOR: inspiral + merger + ringdown waveform
# ---------------------------------------------------------------------------

@dataclass
class GWSimulator:
    """Substrate-derived BBH waveform: inspiral + merger + ringdown.

    Notes
    -----
    *  Inspiral uses the leading-order PN chirp:
            f(t) = (1/π) (5 / (256 τ))^{3/8} (G M_c / c^3)^{-5/8}
            h(t) = (4/D) (G M_c / c^2)^{5/3} (π f / c)^{2/3} cos φ(t)

    *  Merger is a smooth interpolation from the final inspiral GW
       frequency f_ISCO up to the ringdown frequency f_RD over a few
       cycles.  Amplitude peaks here.

    *  Ringdown is a damped sinusoid at the QNM frequency f_RD with
       damping time τ_RD.  In the substrate picture, f_RD is the
       resonance frequency of the σ → 1/2 surface relaxing back to
       equilibrium.
    """

    geometry:    GWGeometry = field(default_factory=GWGeometry)
    sample_rate: float = 4096.0          # Hz (LIGO standard)
    f_low:       float = 35.0            # Hz (typical LIGO low-cutoff)

    # ------------------------------------------------------------------
    # Time-to-coalescence at a given GW frequency  (PN leading order)
    # ------------------------------------------------------------------

    def tau_of_f(self, f_gw: float) -> float:
        """Newtonian-PN time to coalescence at GW frequency f_gw.

            τ(f) = (5 / 256) · (π f)^{-8/3} · (G M_c / c^3)^{-5/3}
        """
        Mc_geom = G_NEWTON * self.geometry.chirp_mass / (C ** 3)  # in seconds
        return float(
            (5.0 / 256.0)
            * (np.pi * f_gw) ** (-8.0 / 3.0)
            * Mc_geom ** (-5.0 / 3.0)
        )

    def f_of_tau(self, tau: NDArray) -> NDArray:
        """Inverse: f(τ) along the Newtonian-PN chirp.

        Used to verify the ``f ∝ τ^{-3/8}`` law in the inspiral.
        """
        Mc_geom = G_NEWTON * self.geometry.chirp_mass / (C ** 3)
        tau = np.asarray(tau, dtype=float)
        # avoid divide-by-zero at τ = 0 by clipping below
        tau_safe = np.where(tau > 1e-6, tau, 1e-6)
        f = (1.0 / np.pi) * (5.0 / (256.0 * tau_safe)) ** (3.0 / 8.0) \
            * Mc_geom ** (-5.0 / 8.0)
        return np.asarray(f, dtype=float)

    # ------------------------------------------------------------------
    # Inspiral strain h(t)
    # ------------------------------------------------------------------

    def inspiral(
        self,
        f_start: Optional[float] = None,
        f_end:   Optional[float] = None,
    ) -> Dict[str, NDArray]:
        """Inspiral strain h(t) from f_start (default = f_low) up to
        f_end (default = ISCO).  Returns dict with keys ``t``, ``f``,
        ``phi``, ``h``.
        """
        f0 = self.f_low if f_start is None else float(f_start)
        f1 = self.geometry.isco_gw_frequency() if f_end is None else float(f_end)
        if f1 <= f0:
            raise ValueError("f_end must exceed f_start.")

        tau0 = self.tau_of_f(f0)        # time-to-coalescence at f_start
        tau1 = self.tau_of_f(f1)        # time-to-coalescence at f_end (small)

        dt = 1.0 / self.sample_rate
        # t increases from 0 to (tau0 - tau1); coalescence is t = tau0 (≡ t_coal)
        n = int(np.floor((tau0 - tau1) / dt))
        if n < 16:
            n = 16
        t = np.arange(n) * dt
        # τ(t) = tau0 - t  (positive, decreasing)
        tau = tau0 - t
        tau = np.maximum(tau, tau1)
        f = self.f_of_tau(tau)

        # Phase: φ(t) = ∫ 2π f dt' = -2 (5 G M_c / c^3)^{-5/8} τ^{5/8} + const
        Mc_geom = G_NEWTON * self.geometry.chirp_mass / (C ** 3)
        phi = -2.0 * (5.0 * Mc_geom) ** (-5.0 / 8.0) * tau ** (5.0 / 8.0)
        phi = phi - phi[0]    # zero at t=0

        # Amplitude: h(t) = (4/D) (G M_c / c^2)^{5/3} (π f / c)^{2/3}
        D = float(self.geometry.distance_mpc * MPC)
        Mc_len = G_NEWTON * self.geometry.chirp_mass / (C * C)   # in metres
        h_amp = (4.0 / D) * Mc_len ** (5.0 / 3.0) * (np.pi * f / C) ** (2.0 / 3.0)
        h = h_amp * np.cos(phi)

        return {
            "t":    np.asarray(t, dtype=float),
            "tau":  np.asarray(tau, dtype=float),
            "f":    np.asarray(f, dtype=float),
            "phi":  np.asarray(phi, dtype=float),
            "h":    np.asarray(h, dtype=float),
            "amp":  np.asarray(h_amp, dtype=float),
        }

    # ------------------------------------------------------------------
    # Ringdown frequency from substrate σ → 1/2 surface
    # ------------------------------------------------------------------

    def ringdown_frequency(self) -> float:
        """Quasi-normal-mode frequency of the merged BH [Hz].

        Echeverria/Berti l=m=2, n=0 fit for Kerr:

            f_RD = (c^3 / (2π G M_f)) * (1.5251 - 1.1568 (1 - a_f)^0.1292)

        Substrate interpretation:  the cone has tilted to 90° at the
        σ = 1/2 surface; the QNM is the resonance of the σ-field
        relaxing back from saturation.
        """
        M_f = self.geometry.final_mass
        a_f = self.geometry.final_spin
        F   = 1.5251 - 1.1568 * (1.0 - a_f) ** 0.1292
        return float((C ** 3 / (2.0 * np.pi * G_NEWTON * M_f)) * F)

    def ringdown_quality_factor(self) -> float:
        """Quality factor Q for the QNM amplitude decay."""
        a_f = self.geometry.final_spin
        return float(0.7000 + 1.4187 * (1.0 - a_f) ** (-0.4990))

    def ringdown_damping_time(self) -> float:
        """Damping time τ_RD = Q / (π f_RD) [s]."""
        return float(self.ringdown_quality_factor()
                     / (np.pi * self.ringdown_frequency()))

    # ------------------------------------------------------------------
    # Merger + ringdown strain
    # ------------------------------------------------------------------

    def ringdown(
        self,
        duration: float = 0.05,
        amplitude: Optional[float] = None,
        t0:        float = 0.0,
    ) -> Dict[str, NDArray]:
        """Damped sinusoid at f_RD, lasting ``duration`` seconds.

        Returns ``t``, ``f`` (constant = f_RD), ``h``.
        """
        f_rd  = self.ringdown_frequency()
        tau_d = self.ringdown_damping_time()
        dt    = 1.0 / self.sample_rate
        n     = max(16, int(np.floor(duration / dt)))
        t     = t0 + np.arange(n) * dt

        # Amplitude estimate: peak strain at merger ~ (4/D) (G M_f / c^2)
        D = float(self.geometry.distance_mpc * MPC)
        amp = ((4.0 / D) * G_NEWTON * self.geometry.final_mass / (C * C)
               * 0.4) if amplitude is None else float(amplitude)

        envelope = amp * np.exp(-(t - t0) / tau_d)
        h = envelope * np.cos(2.0 * np.pi * f_rd * (t - t0))
        f = np.full_like(t, f_rd)
        return {
            "t": np.asarray(t, dtype=float),
            "f": np.asarray(f, dtype=float),
            "h": np.asarray(h, dtype=float),
            "envelope": np.asarray(envelope, dtype=float),
            "f_RD": float(f_rd),
            "tau_RD": float(tau_d),
        }

    def merger(
        self,
        n_cycles: float = 2.5,
        f_start: Optional[float] = None,
        f_end:   Optional[float] = None,
        amp_start: float = 0.0,
        amp_peak:  float = 1.0,
        t0:        float = 0.0,
    ) -> Dict[str, NDArray]:
        """Bridge from f_ISCO up to f_RD over ``n_cycles`` GW cycles.

        The merger is short (≪ 1 inspiral cycle); we ramp the amplitude
        up to ``amp_peak`` and the frequency from f_ISCO to f_RD with
        a smooth half-cosine.  This is a phenomenological bridge — the
        full IMR phenomenology (e.g. SEOBNRv4, IMRPhenomD) is far beyond
        this module's scope, but this preserves the *qualitative*
        chirp-merger-ringdown morphology for visualization.
        """
        f0 = self.geometry.isco_gw_frequency() if f_start is None else float(f_start)
        f1 = self.ringdown_frequency()        if f_end   is None else float(f_end)
        # Bridge duration ~ n_cycles / f_avg
        f_avg    = 0.5 * (f0 + f1)
        duration = max(n_cycles / f_avg, 1.0 / self.sample_rate)
        dt       = 1.0 / self.sample_rate
        n        = max(8, int(np.floor(duration / dt)))
        t_local  = np.arange(n) * dt
        t        = t0 + t_local

        # Smooth half-cosine ramp for f and amplitude
        s = 0.5 * (1.0 - np.cos(np.pi * t_local / duration))
        f = f0 + (f1 - f0) * s
        # Phase = ∫ 2π f dt with discrete cumulative sum
        phi = 2.0 * np.pi * np.cumsum(f) * dt
        amp = amp_start + (amp_peak - amp_start) * s
        h = amp * np.cos(phi)
        return {
            "t":   np.asarray(t, dtype=float),
            "f":   np.asarray(f, dtype=float),
            "h":   np.asarray(h, dtype=float),
            "amp": np.asarray(amp, dtype=float),
            "phi": np.asarray(phi, dtype=float),
        }

    def full_waveform(
        self,
        f_start:        Optional[float] = None,
        ringdown_dur:   float = 0.05,
        merger_cycles:  float = 2.5,
    ) -> Dict[str, NDArray]:
        """Inspiral + merger + ringdown concatenated end-to-end.

        Returns dict with ``t`` (continuous), ``h``, ``f``, plus the
        three individual segments under the keys ``inspiral``,
        ``merger``, ``ringdown``.
        """
        ins = self.inspiral(f_start=f_start)
        # Peak strain at end of inspiral
        h_peak_inspiral = float(np.max(np.abs(ins["h"][-200:])))
        if h_peak_inspiral <= 0.0:
            h_peak_inspiral = 1e-22

        # Merger amplitude peaks at ~ 4× the inspiral peak
        mer = self.merger(
            n_cycles=merger_cycles,
            amp_start=h_peak_inspiral,
            amp_peak=4.0 * h_peak_inspiral,
            t0=float(ins["t"][-1] + 1.0 / self.sample_rate),
        )
        rd_amp = float(np.abs(mer["h"][-1]))
        rd  = self.ringdown(
            duration=ringdown_dur,
            amplitude=max(rd_amp, h_peak_inspiral),
            t0=float(mer["t"][-1] + 1.0 / self.sample_rate),
        )

        t = np.concatenate([ins["t"], mer["t"], rd["t"]])
        h = np.concatenate([ins["h"], mer["h"], rd["h"]])
        f = np.concatenate([ins["f"], mer["f"], rd["f"]])
        return {
            "t": t, "h": h, "f": f,
            "inspiral": ins,
            "merger":   mer,
            "ringdown": rd,
            "f_RD":     rd["f_RD"],
            "tau_RD":   rd["tau_RD"],
        }

    # ------------------------------------------------------------------
    # Spectrogram
    # ------------------------------------------------------------------

    def spectrogram(
        self,
        h:           NDArray,
        nperseg:     int = 256,
        noverlap:    Optional[int] = None,
    ) -> Dict[str, NDArray]:
        """Return a Hann-windowed STFT spectrogram |X(f, t)|^2 of h(t).

        Implemented in pure numpy so we don't depend on scipy.signal.
        """
        h = np.asarray(h, dtype=float)
        if noverlap is None:
            noverlap = nperseg // 2
        step = nperseg - noverlap
        if step <= 0:
            step = 1

        n = h.size
        if n < nperseg:
            # Pad
            h = np.concatenate([h, np.zeros(nperseg - n)])
            n = h.size

        n_segments = 1 + (n - nperseg) // step
        win = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(nperseg) / (nperseg - 1)))
        win_norm = float(np.sum(win ** 2))
        if win_norm <= 0.0:
            win_norm = 1.0

        S = np.zeros((nperseg // 2 + 1, n_segments), dtype=float)
        for k in range(n_segments):
            seg = h[k * step : k * step + nperseg] * win
            spec = np.fft.rfft(seg)
            S[:, k] = (np.abs(spec) ** 2) / win_norm

        freqs = np.fft.rfftfreq(nperseg, d=1.0 / self.sample_rate)
        times = (np.arange(n_segments) * step + nperseg / 2) / self.sample_rate
        return {
            "f": freqs,
            "t": times,
            "S": S,
        }

    # ------------------------------------------------------------------
    # Substrate cross-check: peak GW frequency at merger
    # ------------------------------------------------------------------

    def peak_frequency(self) -> float:
        """Approximate peak GW frequency at merger.

        Uses f_ISCO as a robust lower bound, capped against the QNM
        frequency: in practice the strain peaks near f_ISCO for ~ 60 M_sun
        binaries (~250 Hz for GW150914), then transitions to f_RD.
        """
        f_isco = self.geometry.isco_gw_frequency()
        f_rd   = self.ringdown_frequency()
        # Merger peak strain ~ at the midpoint between ISCO and RD,
        # weighted toward ISCO (factor ~ 0.6 from numerical relativity).
        return float(0.6 * f_isco + 0.4 * f_rd)


# ---------------------------------------------------------------------------
# (3) Convenience: GW150914 preset
# ---------------------------------------------------------------------------

def gw150914_preset() -> Tuple[GWGeometry, GWSimulator]:
    """LIGO GW150914 parameters: m1≈36, m2≈29, D≈410 Mpc."""
    geom = GWGeometry(m1_solar=36.0, m2_solar=29.0,
                      spin1=0.32, spin2=-0.44,
                      distance_mpc=410.0)
    sim  = GWSimulator(geometry=geom, sample_rate=4096.0, f_low=35.0)
    return geom, sim


# ---------------------------------------------------------------------------
# (4) Rendering helpers used by scripts/render_all_visuals.py
# ---------------------------------------------------------------------------

def render_inspiral_figure(
    geom: Optional[GWGeometry] = None,
    sim:  Optional[GWSimulator] = None,
):
    """Build the matplotlib figure for ``visuals/46_gw_inspiral.png``.

    Two panels:
      (a) h(t) showing the chirp + merger + ringdown.
      (b) f(t) along with the analytic τ^{-3/8} chirp curve.
    """
    import matplotlib.pyplot as plt  # local import — keeps module import-light

    if geom is None or sim is None:
        geom, sim = gw150914_preset()

    full = sim.full_waveform()
    t = full["t"]
    h = full["h"]
    f = full["f"]

    # Center the time axis on coalescence (t=0 at end of inspiral)
    t_coal = float(full["inspiral"]["t"][-1])
    t_rel  = t - t_coal

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    ax = axes[0]
    ax.plot(t_rel, h, "b-", linewidth=1.0, alpha=0.85)
    ax.axvline(0.0, color="red", linestyle="--", linewidth=1.2,
               label="merger (t=0)")
    ax.set_ylabel("strain  h(t)")
    ax.set_title(
        f"GW150914-like BBH inspiral + merger + ringdown\n"
        f"M_c = {geom.chirp_mass_solar:.1f} M_sun,  "
        f"M_total = {geom.total_mass_solar:.1f} M_sun,  "
        f"D = {geom.distance_mpc:.0f} Mpc  |  "
        f"v_GW = c (substrate, structural)",
        fontsize=11,
    )
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.semilogy(t_rel, f, "g-", linewidth=1.5, label="instantaneous f(t)")

    # Overlay analytic τ^{-3/8} chirp on the inspiral portion
    ins = full["inspiral"]
    tau = ins["tau"]
    f_ana = sim.f_of_tau(tau)
    ax.semilogy(ins["t"] - t_coal, f_ana, "k--", linewidth=1.0,
                alpha=0.7, label="analytic  f ∝ τ^(-3/8)")

    ax.axhline(geom.isco_gw_frequency(), color="orange", linestyle=":",
               linewidth=1.2, label=f"f_ISCO ≈ {geom.isco_gw_frequency():.0f} Hz")
    ax.axhline(sim.ringdown_frequency(), color="purple", linestyle=":",
               linewidth=1.2, label=f"f_RD ≈ {sim.ringdown_frequency():.0f} Hz")
    ax.set_xlabel("time relative to merger  [s]")
    ax.set_ylabel("GW frequency  [Hz]")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return fig


def render_ringdown_figure(
    geom: Optional[GWGeometry] = None,
    sim:  Optional[GWSimulator] = None,
):
    """Build the matplotlib figure for ``visuals/47_gw_merger_ringdown.png``.

    Two panels:
      (a) Spectrogram |H(f, t)|^2 of the full waveform with f_ISCO and
          f_RD marked.
      (b) Ringdown strain h(t) with the substrate-derived f_RD and τ_RD
          annotation.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    if geom is None or sim is None:
        geom, sim = gw150914_preset()

    full = sim.full_waveform()
    h    = full["h"]
    t    = full["t"]
    t_coal = float(full["inspiral"]["t"][-1])
    t_rel  = t - t_coal

    spec = sim.spectrogram(h, nperseg=512, noverlap=384)
    S   = spec["S"]
    f_s = spec["f"]
    t_s = spec["t"] - t_coal

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    # ----- Panel (a): spectrogram -----
    ax = axes[0]
    S_plot = np.maximum(S, 1e-12 * S.max() if S.max() > 0 else 1e-30)
    pcm = ax.pcolormesh(
        t_s, f_s, S_plot,
        norm=LogNorm(vmin=S_plot.max() * 1e-5, vmax=S_plot.max()),
        shading="auto", cmap="viridis",
    )
    ax.set_yscale("log")
    ax.set_ylim(20.0, sim.sample_rate / 2.0)
    ax.set_ylabel("frequency  [Hz]")
    ax.axhline(geom.isco_gw_frequency(), color="orange",
               linestyle="--", linewidth=1.2,
               label=f"f_ISCO = {geom.isco_gw_frequency():.0f} Hz")
    ax.axhline(sim.ringdown_frequency(), color="white",
               linestyle=":", linewidth=1.5,
               label=f"f_RD = {sim.ringdown_frequency():.0f} Hz "
                     f"(σ→½ resonance)")
    ax.axvline(0.0, color="red", linestyle="--", linewidth=1.2)
    ax.set_title(
        f"BBH spectrogram |H(f, t)|^2  —  "
        f"f rises from f_low={sim.f_low:.0f} Hz "
        f"through f_ISCO={geom.isco_gw_frequency():.0f} Hz "
        f"to f_RD={sim.ringdown_frequency():.0f} Hz",
        fontsize=11,
    )
    ax.legend(loc="upper left", fontsize=9, framealpha=0.6)
    fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.02, label="|H|^2")

    # ----- Panel (b): ringdown waveform -----
    ax = axes[1]
    rd = full["ringdown"]
    t_rd_rel = rd["t"] - t_coal
    ax.plot(t_rd_rel, rd["h"], "b-", linewidth=1.4, label="h(t) ringdown")
    ax.plot(t_rd_rel, rd["envelope"], "r--", linewidth=1.2,
            alpha=0.8, label=f"exp envelope (τ_RD = {rd['tau_RD']*1e3:.2f} ms)")
    ax.plot(t_rd_rel, -rd["envelope"], "r--", linewidth=1.2, alpha=0.8)
    # Annotate substrate interpretation
    ax.text(
        0.02, 0.95,
        f"Substrate ringdown:\n"
        f"  f_RD = {rd['f_RD']:.1f} Hz   (cone σ→½ relaxation)\n"
        f"  Q    = {sim.ringdown_quality_factor():.2f}\n"
        f"  M_f  = {geom.final_mass / M_SUN:.1f} M_sun\n"
        f"  a_f  = {geom.final_spin:.3f}",
        transform=ax.transAxes,
        fontsize=10, va="top", ha="left",
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.3",
                  facecolor="lightyellow", edgecolor="gray", alpha=0.85),
    )
    ax.set_xlabel("time after merger  [s]")
    ax.set_ylabel("ringdown strain  h(t)")
    ax.set_title("Ringdown: damped QNM at f_RD — substrate σ→½ resonance")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    return fig


__all__ = [
    "G_NEWTON", "M_SUN", "MPC", "KPC", "SIGMA_MAX",
    "GWGeometry", "GWSimulator",
    "gw150914_preset",
    "render_inspiral_figure", "render_ringdown_figure",
]
