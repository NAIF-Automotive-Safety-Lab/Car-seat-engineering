"""Substrate frequency sonification — convert ω_b cone-bouncing to audio.

The B3 / Stiff-Medium framework gives every massive particle a rest-frame
cone-bouncing frequency

    ω_b² = (c/ξ)² + (γ / 2ρ)²              (cone_bouncing_protocol)

with rest mass m c² = ℏ ω_b. For the electron, ω_b ≈ 7.76e20 rad/s, far
above audible (20 Hz - 20 kHz). This module rescales the substrate
frequency tower into the audible range while preserving mass *ratios*, so
the lepton/quark mass spectrum becomes a chord that can be heard.

Anchor convention
-----------------
The electron is mapped to A3 = 220 Hz, then every other particle's tone is

    f_audio = f_e_audio * (m_particle / m_e)

so the *ratio* of audible frequencies equals the ratio of rest masses. The
muon, for example, sits at 220 * 207 ≈ 45.5 kHz — above the audible
ceiling, so we additionally fold by halving (octave-down) until the tone
falls in 20 Hz - 20 kHz. This preserves *pitch class* (mass mod log₂)
while keeping everything audible. The user can hear the substrate-derived
mass ratios as musical intervals.

Antimatter convention
---------------------
B3 §18.49: antimatter is the mirror cone-bouncing solution (hill ↔ trough,
same ω_b, same |m|). Stereo channels split matter (left) from antimatter
(right) — both are the same pure tone, illustrating m_p̄ = m_p (BASE
16-ppt CPT result).

Drag envelope
-------------
The drag γ enters ω_b via a frequency-shift but ALSO sets the dissipation
timescale τ_decay = 2ρ/γ. We use this as the audio envelope decay
constant (ADSR-like): the ratio γ / ω_b sets how quickly the tone rings
down, so heavier particles (larger γ relative to ξ) ring down faster.

Outputs
-------
* visuals/audio/lepton_tower.wav     — e/μ/τ chord
* visuals/audio/quark_tower.wav      — u,c,t and d,s,b two chords
* visuals/audio/cone_bouncing.wav    — frequency sweep + drag envelope
* visuals/audio/cosmology_evolution.wav — σ from 1/2 → 0 sonified
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from scipy.io import wavfile


# ---------------------------------------------------------------------------
# Physical constants (SI) — copied here to keep this module self-contained
# ---------------------------------------------------------------------------

HBAR: float = 1.054571817e-34       # J·s
C_LIGHT: float = 2.99792458e8       # m/s
EV_J: float = 1.602176634e-19       # J/eV
MEV_J: float = 1.0e6 * EV_J         # J/MeV

# PDG masses (MeV/c²) — same values used in drag_mass_generator.py
PDG_MASSES_MEV: Dict[str, float] = {
    "e":   0.51099895,
    "mu":  105.6583755,
    "tau": 1776.86,
    "u":   2.16,
    "d":   4.67,
    "s":   93.4,
    "c":   1270.0,
    "b":   4180.0,
    "t":   172570.0,
}


# ---------------------------------------------------------------------------
# Audio constants
# ---------------------------------------------------------------------------

DEFAULT_SAMPLE_RATE: int = 44100
ELECTRON_AUDIO_HZ: float = 220.0      # A3 — anchor for the audible mapping
AUDIBLE_LOW: float = 20.0
AUDIBLE_HIGH: float = 20_000.0
INT16_MAX: int = 32_767


# ---------------------------------------------------------------------------
# ω_b helpers
# ---------------------------------------------------------------------------

def compton_xi(mass_mev: float) -> float:
    """Compton wavelength ξ = ℏ / (m c)  in metres."""
    m_kg = (mass_mev * MEV_J) / C_LIGHT ** 2
    return HBAR / (m_kg * C_LIGHT)


def omega_b_compton(mass_mev: float) -> float:
    """Bare (γ → 0) cone-bouncing frequency ω₀ = c / ξ in rad/s.

    Equivalent to mass·c²/ℏ. We use the Compton form to keep the link to
    cone-bouncing explicit.
    """
    return C_LIGHT / compton_xi(mass_mev)


def fold_to_audible(freq_hz: float,
                    low: float = AUDIBLE_LOW,
                    high: float = AUDIBLE_HIGH) -> float:
    """Octave-fold a frequency into [low, high].

    Halve until below `high`; double until above `low`. Preserves pitch
    class (i.e. log₂ mod 1) so that octave-equivalent ratios remain
    audible. Returns the folded frequency.
    """
    if freq_hz <= 0.0 or not math.isfinite(freq_hz):
        return low
    f = freq_hz
    # avoid pathological loops on degenerate input
    for _ in range(200):
        if f > high:
            f *= 0.5
            continue
        if f < low:
            f *= 2.0
            continue
        break
    return f


def mass_to_audio_hz(mass_mev: float,
                     anchor_mev: float = PDG_MASSES_MEV["e"],
                     anchor_hz: float = ELECTRON_AUDIO_HZ,
                     fold: bool = True) -> float:
    """Rescale a particle mass to an audible frequency.

    The *ratio* anchor_hz / anchor_mev is preserved: every other particle's
    audio frequency equals anchor_hz * (mass / anchor_mev). When `fold`
    is True we additionally octave-fold into the audible window.
    """
    f_raw = anchor_hz * (mass_mev / anchor_mev)
    return fold_to_audible(f_raw) if fold else f_raw


# ---------------------------------------------------------------------------
# Tone generation
# ---------------------------------------------------------------------------

def tone_from_omega(omega: float,
                    duration: float,
                    sample_rate: int = DEFAULT_SAMPLE_RATE,
                    amplitude: float = 0.5,
                    decay: Optional[float] = None,
                    attack: float = 0.01,
                    release: float = 0.05) -> np.ndarray:
    """Generate a pure sine tone at angular frequency ω (rad/s).

    Parameters
    ----------
    omega : float
        Angular frequency (rad/s). Will be reduced to a Hz value via
        f = omega / (2π). The caller is responsible for ensuring this
        falls in the audible window — use `mass_to_audio_hz` first if
        you start from a substrate mass.
    duration : float
        Tone length in seconds.
    sample_rate : int
        Sample rate in Hz. 44100 by default (CD quality).
    amplitude : float
        Linear amplitude in [0, 1]. Scales the carrier before envelope.
    decay : Optional[float]
        Exponential decay τ in seconds. If None, no decay (only attack/
        release). When given, samples carry an extra exp(-t/decay)
        envelope — useful for the drag-driven cone-bouncing analogy.
    attack, release : float
        Linear ramp times (seconds) to suppress click artefacts at the
        start and end of the tone.

    Returns
    -------
    np.ndarray of float32 in [-1, 1], length int(duration * sample_rate).
    """
    if duration <= 0.0:
        return np.zeros(0, dtype=np.float32)
    if sample_rate <= 0:
        raise ValueError("sample_rate must be positive")

    n = int(duration * sample_rate)
    t = np.arange(n, dtype=np.float64) / float(sample_rate)
    # Pure carrier at f = ω/(2π)
    f_hz = float(omega) / (2.0 * math.pi)
    wave = np.sin(2.0 * math.pi * f_hz * t)

    # ADSR-lite: linear attack, sustain, linear release, optional decay
    env = np.ones_like(t)
    n_attack = max(1, int(attack * sample_rate))
    n_release = max(1, int(release * sample_rate))
    if n_attack < n:
        env[:n_attack] = np.linspace(0.0, 1.0, n_attack)
    if n_release < n:
        env[-n_release:] = np.linspace(1.0, 0.0, n_release)
    if decay is not None and decay > 0.0:
        env *= np.exp(-t / decay)

    out = (amplitude * wave * env).astype(np.float32)
    # Defensive clip to [-1, 1] in case the user passed amplitude > 1
    np.clip(out, -1.0, 1.0, out=out)
    return out


def mix(*signals: np.ndarray) -> np.ndarray:
    """Sum-and-normalize an arbitrary number of mono float32 signals.

    Output is float32 in [-1, 1]. Length equals the longest input;
    shorter signals are zero-padded.
    """
    if not signals:
        return np.zeros(0, dtype=np.float32)
    n = max(len(s) for s in signals)
    acc = np.zeros(n, dtype=np.float64)
    for s in signals:
        if len(s) == 0:
            continue
        acc[: len(s)] += s
    peak = float(np.max(np.abs(acc))) or 1.0
    return (acc / peak).astype(np.float32)


def to_stereo(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Stack two mono signals into an (N, 2) interleaved stereo array.

    Pads the shorter channel with zeros. Each channel is independently
    clipped to [-1, 1].
    """
    n = max(len(left), len(right))
    out = np.zeros((n, 2), dtype=np.float32)
    if len(left) > 0:
        out[: len(left), 0] = left
    if len(right) > 0:
        out[: len(right), 1] = right
    np.clip(out, -1.0, 1.0, out=out)
    return out


def write_wav(path: str, signal: np.ndarray,
              sample_rate: int = DEFAULT_SAMPLE_RATE) -> None:
    """Write a float32 (mono or (N,2) stereo) signal as a 16-bit PCM WAV.

    Creates parent directories if they don't already exist.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Convert to int16 PCM
    sig = np.asarray(signal, dtype=np.float32)
    np.clip(sig, -1.0, 1.0, out=sig)
    pcm = (sig * INT16_MAX).astype(np.int16)
    wavfile.write(str(p), sample_rate, pcm)


# ---------------------------------------------------------------------------
# Output paths (resolved relative to the project root: parent of `src/`)
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    """Resolve the project root: parent of src/stiff_medium/."""
    return Path(__file__).resolve().parents[2]


def _audio_dir() -> Path:
    d = _project_root() / "visuals" / "audio"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Tower compositions
# ---------------------------------------------------------------------------

@dataclass
class TowerEntry:
    """A single particle in a sonified mass tower."""
    name: str
    mass_mev: float
    audio_hz: float
    omega_b_rad_s: float


def _build_tower(names: Iterable[str]) -> List[TowerEntry]:
    """Build TowerEntry list with audio frequencies anchored to electron."""
    entries: List[TowerEntry] = []
    for n in names:
        m = PDG_MASSES_MEV[n]
        entries.append(TowerEntry(
            name=n,
            mass_mev=m,
            audio_hz=mass_to_audio_hz(m),
            omega_b_rad_s=omega_b_compton(m),
        ))
    return entries


def _chord_signal(entries: List[TowerEntry],
                  duration: float,
                  sample_rate: int) -> np.ndarray:
    """Sum the entries' tones into a chord. Decay ~ inversely with mass
    so heavier particles ring down faster (drag analogy)."""
    # Reference decay: electron rings 80% of duration; heavier particles
    # decay proportionally to 1/mass (mass ∝ γ at fixed ξ ratio).
    sigs: List[np.ndarray] = []
    m_ref = PDG_MASSES_MEV["e"]
    for e in entries:
        # decay scale: never < 5% of duration, never > 100%
        ratio = m_ref / e.mass_mev
        decay = max(0.05 * duration, min(0.8 * duration, 0.6 * duration * ratio ** 0.25))
        omega_audio = 2.0 * math.pi * e.audio_hz
        s = tone_from_omega(
            omega=omega_audio,
            duration=duration,
            sample_rate=sample_rate,
            amplitude=0.5,
            decay=decay,
        )
        sigs.append(s)
    return mix(*sigs)


def lepton_tower_audio(duration: float = 4.0,
                       sample_rate: int = DEFAULT_SAMPLE_RATE,
                       output_path: Optional[str] = None) -> Tuple[str, List[TowerEntry]]:
    """Render the e/μ/τ lepton chord with substrate-derived mass ratios.

    Each lepton's audible frequency = 220 Hz × (m_lep / m_e), octave-folded
    into the audible band so that the *ratio* (and therefore the pitch
    interval) of the substrate mass tower is preserved.

    Stereo: left = matter, right = antimatter. Per B3 §18.49 antimatter
    is the same |m| (CPT), so the right channel is identical — illustrating
    the BASE 16-ppt m_p̄ = m_p result.

    Returns (path, entries).
    """
    entries = _build_tower(["e", "mu", "tau"])
    mono = _chord_signal(entries, duration, sample_rate)
    stereo = to_stereo(mono, mono.copy())
    out = output_path or str(_audio_dir() / "lepton_tower.wav")
    write_wav(out, stereo, sample_rate=sample_rate)
    return out, entries


def quark_tower_audio(duration: float = 5.0,
                      sample_rate: int = DEFAULT_SAMPLE_RATE,
                      output_path: Optional[str] = None) -> Tuple[str, List[TowerEntry]]:
    """Render the (u,c,t) and (d,s,b) quark chords as a stereo pair.

    Up-type tower (left channel) and down-type tower (right channel) are
    rendered as separate chords stacked in stereo — letting the listener
    hear the I3 = +1/2 vs I3 = -1/2 split from quark_mass_spectrum.py.

    Returns (path, entries) where entries lists all six quarks.
    """
    up_entries = _build_tower(["u", "c", "t"])
    down_entries = _build_tower(["d", "s", "b"])
    up_mono = _chord_signal(up_entries, duration, sample_rate)
    down_mono = _chord_signal(down_entries, duration, sample_rate)
    stereo = to_stereo(up_mono, down_mono)
    out = output_path or str(_audio_dir() / "quark_tower.wav")
    write_wav(out, stereo, sample_rate=sample_rate)
    return out, up_entries + down_entries


def cone_bouncing_audio(duration: float = 6.0,
                        sample_rate: int = DEFAULT_SAMPLE_RATE,
                        output_path: Optional[str] = None,
                        gamma_over_omega_max: float = 1.5) -> Tuple[str, dict]:
    """Sweep the dimensionless drag γ/ω₀ and sonify the resulting ω_b shift.

    Per cone_bouncing_protocol.py:
        ω_b² = ω₀² + (γ / 2ρ)² = ω₀² (1 + (γ̃ / 2)²)        (γ̃ = γ/ω₀)

    We sweep γ̃ from 0 → gamma_over_omega_max linearly. The carrier is
    anchored to the electron audible frequency (220 Hz) for γ=0; as drag
    grows the audible frequency rises by √(1 + (γ̃/2)²). The drag
    itself drives the amplitude envelope (exp decay with τ ∝ 1/γ),
    illustrating that drag is *both* the mass-generator and the dissipator.

    Returns (path, info_dict).
    """
    n = int(duration * sample_rate)
    t = np.arange(n, dtype=np.float64) / float(sample_rate)

    # Linear sweep of dimensionless drag
    gamma_tilde = np.linspace(0.0, gamma_over_omega_max, n)
    f_e = ELECTRON_AUDIO_HZ
    f_t = f_e * np.sqrt(1.0 + (gamma_tilde / 2.0) ** 2)

    # Phase = ∫ 2π f(t) dt  via cumulative trapezoid for chirp synthesis
    dt = 1.0 / float(sample_rate)
    phase = 2.0 * math.pi * np.cumsum(f_t) * dt
    carrier = np.sin(phase)

    # Drag-driven envelope: starts loud (γ small, slow decay), gets
    # progressively damped (large γ → fast decay). Use amplitude
    # ∝ 1/(1 + γ̃) to model dissipation eating the bound state.
    env_amp = 1.0 / (1.0 + gamma_tilde)
    # Plus a final fade
    fade = np.ones_like(t)
    fade_n = max(1, int(0.05 * n))
    fade[-fade_n:] = np.linspace(1.0, 0.0, fade_n)

    mono = (0.6 * carrier * env_amp * fade).astype(np.float32)

    # Stereo: L = matter, R = antimatter (mirror solution, same ω_b)
    stereo = to_stereo(mono, mono.copy())
    out = output_path or str(_audio_dir() / "cone_bouncing.wav")
    write_wav(out, stereo, sample_rate=sample_rate)

    info = {
        "gamma_tilde_start": float(gamma_tilde[0]),
        "gamma_tilde_end": float(gamma_tilde[-1]),
        "f_start_hz": float(f_t[0]),
        "f_end_hz": float(f_t[-1]),
        "frequency_shift_ratio": float(f_t[-1] / f_t[0]),
    }
    return out, info


def cosmology_evolution_audio(duration: float = 8.0,
                              sample_rate: int = DEFAULT_SAMPLE_RATE,
                              output_path: Optional[str] = None) -> Tuple[str, dict]:
    """Sonify σ from 1/2 (Big Bang saturated state) → 0 (modern vacuum).

    Per saturation.py and desaturation.py:
      - σ_max = 1/2 is the substrate elastic limit (§18.39).
      - At cosmic origin σ ≈ 1/2 (everywhere saturated, no clocks).
      - The CMB de-saturation phase transition releases latent energy and
        the universe-wide σ decays toward 0 as a⁻³ on average.

    Mapping: σ ↔ how heavily the substrate is "pinched" so substrate
    fluctuations are detuned by σ. We render an organ-like drone whose
    pitch slides from a low rumble (σ = 1/2 ⇒ saturated, dense) up to
    the lepton anchor 220 Hz (σ → 0 ⇒ ordinary vacuum, low strain).

    Stereo: L = matter perspective, R = antimatter perspective. They
    coincide (CPT) but with a slight 0.5% detune to make the shared
    history audible as a beating chorus.

    Returns (path, info_dict).
    """
    n = int(duration * sample_rate)
    t = np.arange(n, dtype=np.float64) / float(sample_rate)

    # σ(t) decays exponentially from 0.5 → 0 with a half-life set so
    # the curve has interesting motion across the full duration.
    tau = duration / 4.0
    sigma = 0.5 * np.exp(-t / tau)

    # Audible mapping: at σ=0.5 the substrate is so saturated that the
    # bound-state frequency is suppressed (low rumble); at σ=0 the
    # ordinary cone-bouncing frequency lights up.
    # Use f(σ) = f_low + (f_high - f_low) * (1 - σ/σ_max) to slide.
    f_low = 40.0
    f_high = ELECTRON_AUDIO_HZ
    f_t = f_low + (f_high - f_low) * (1.0 - sigma / 0.5)

    dt = 1.0 / float(sample_rate)
    phase = 2.0 * math.pi * np.cumsum(f_t) * dt
    carrier = np.sin(phase)

    # Add a 5th harmonic at low amplitude to give the drone body
    overtone = 0.25 * np.sin(2.0 * phase)
    body = carrier + overtone

    # Amplitude swells as σ decreases (vacuum becomes "alive")
    amp = 0.3 + 0.5 * (1.0 - sigma / 0.5)
    # Final fade
    fade = np.ones_like(t)
    fade_n = max(1, int(0.08 * n))
    fade[-fade_n:] = np.linspace(1.0, 0.0, fade_n)
    n_attack = max(1, int(0.05 * n))
    fade[:n_attack] *= np.linspace(0.0, 1.0, n_attack)

    mono = (0.5 * body * amp * fade).astype(np.float32)

    # Antimatter slight detune for chorus effect
    f_t_anti = f_t * 1.005
    phase_anti = 2.0 * math.pi * np.cumsum(f_t_anti) * dt
    body_anti = np.sin(phase_anti) + 0.25 * np.sin(2.0 * phase_anti)
    mono_anti = (0.5 * body_anti * amp * fade).astype(np.float32)

    stereo = to_stereo(mono, mono_anti)
    out = output_path or str(_audio_dir() / "cosmology_evolution.wav")
    write_wav(out, stereo, sample_rate=sample_rate)

    info = {
        "sigma_start": float(sigma[0]),
        "sigma_end": float(sigma[-1]),
        "f_start_hz": float(f_t[0]),
        "f_end_hz": float(f_t[-1]),
    }
    return out, info


# ---------------------------------------------------------------------------
# CLI entry-point — render all four files at once
# ---------------------------------------------------------------------------

def render_all(sample_rate: int = DEFAULT_SAMPLE_RATE) -> Dict[str, str]:
    """Render all four sonification files. Returns map name → path."""
    paths: Dict[str, str] = {}
    paths["lepton"], _ = lepton_tower_audio(sample_rate=sample_rate)
    paths["quark"], _ = quark_tower_audio(sample_rate=sample_rate)
    paths["cone_bouncing"], _ = cone_bouncing_audio(sample_rate=sample_rate)
    paths["cosmology"], _ = cosmology_evolution_audio(sample_rate=sample_rate)
    return paths


if __name__ == "__main__":  # pragma: no cover
    out = render_all()
    for k, v in out.items():
        print(f"{k:>14s}  ->  {v}")
