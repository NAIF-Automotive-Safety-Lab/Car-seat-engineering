"""Tests for src.stiff_medium.substrate_audio.

Covers
------
* `tone_from_omega` produces correct length and frequency content.
* WAV files round-trip with the right sample rate, bit depth and channel
  count (mono for tones, stereo for chord/cosmology renders).
* Mass-tower frequencies satisfy the substrate ratio relation
  f_audio(particle) / f_audio(electron) == m_particle / m_electron
  (modulo octave-folding for particles whose raw audible frequency
  exceeds 20 kHz).
* The drag-driven cone-bouncing render shows the expected √(1 + (γ̃/2)²)
  upshift between start and end of the sweep.
* The cosmology render encodes σ: 1/2 → 0 with a corresponding pitch
  glide upward.

We DO NOT assert on PCM byte-level content; only on the physics-meaningful
frequency content (estimated via FFT) and metadata of each rendered file.
"""

from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.io import wavfile

from src.stiff_medium.substrate_audio import (
    AUDIBLE_HIGH,
    AUDIBLE_LOW,
    DEFAULT_SAMPLE_RATE,
    ELECTRON_AUDIO_HZ,
    INT16_MAX,
    PDG_MASSES_MEV,
    compton_xi,
    cone_bouncing_audio,
    cosmology_evolution_audio,
    fold_to_audible,
    lepton_tower_audio,
    mass_to_audio_hz,
    mix,
    omega_b_compton,
    quark_tower_audio,
    render_all,
    to_stereo,
    tone_from_omega,
    write_wav,
)


# ---------------------------------------------------------------------------
# Fixtures: render into a tmp dir so tests don't litter visuals/audio
# (the production functions also accept output_path overrides).
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_audio_dir(tmp_path):
    return tmp_path


# ---------------------------------------------------------------------------
# Mapping math
# ---------------------------------------------------------------------------

def test_compton_xi_inverts_mass():
    m = PDG_MASSES_MEV["e"]
    xi = compton_xi(m)
    # ω₀ = c / ξ should give m c²/ℏ (in rad/s)
    omega = omega_b_compton(m)
    # consistency: ω = c/ξ
    assert math.isclose(omega * xi, 2.99792458e8, rel_tol=1e-12)


def test_fold_to_audible_keeps_pitch_class():
    # An out-of-band frequency should fold into [20, 20000] preserving
    # the log2-mod-1 part (octave equivalence).
    f_in = 50_000.0
    f_out = fold_to_audible(f_in)
    assert AUDIBLE_LOW <= f_out <= AUDIBLE_HIGH
    # log2 difference should be very close to an integer
    log2_diff = math.log2(f_in / f_out)
    assert math.isclose(log2_diff, round(log2_diff), abs_tol=1e-9)


def test_fold_handles_subaudible():
    f_out = fold_to_audible(5.0)
    assert AUDIBLE_LOW <= f_out <= AUDIBLE_HIGH


def test_mass_to_audio_anchor_is_220():
    f = mass_to_audio_hz(PDG_MASSES_MEV["e"])
    assert math.isclose(f, ELECTRON_AUDIO_HZ, rel_tol=1e-12)


def test_mass_to_audio_unfolded_preserves_ratio():
    """Without octave-folding, audio Hz / electron Hz == mass ratio exactly."""
    m_mu = PDG_MASSES_MEV["mu"]
    m_e = PDG_MASSES_MEV["e"]
    f_mu = mass_to_audio_hz(m_mu, fold=False)
    f_e = mass_to_audio_hz(m_e, fold=False)
    assert math.isclose(f_mu / f_e, m_mu / m_e, rel_tol=1e-12)


def test_mass_to_audio_folded_preserves_pitch_class():
    """With folding, ratios coincide modulo a factor of 2^k."""
    m_mu = PDG_MASSES_MEV["mu"]
    m_e = PDG_MASSES_MEV["e"]
    f_mu = mass_to_audio_hz(m_mu, fold=True)
    f_e = mass_to_audio_hz(m_e, fold=True)
    raw_ratio = m_mu / m_e
    audible_ratio = f_mu / f_e
    # log2 difference must be an integer (octave-equivalence)
    log2_diff = math.log2(raw_ratio / audible_ratio)
    assert math.isclose(log2_diff, round(log2_diff), abs_tol=1e-9)
    # And f_mu sits in audible range
    assert AUDIBLE_LOW <= f_mu <= AUDIBLE_HIGH


# ---------------------------------------------------------------------------
# tone_from_omega
# ---------------------------------------------------------------------------

def test_tone_length_and_dtype():
    sr = 44100
    dur = 0.25
    s = tone_from_omega(omega=2 * math.pi * 440.0, duration=dur, sample_rate=sr)
    assert s.dtype == np.float32
    assert len(s) == int(dur * sr)
    assert np.all(np.isfinite(s))
    assert np.max(np.abs(s)) <= 1.0


def test_tone_frequency_via_fft():
    """The peak of the magnitude spectrum should sit at the requested f."""
    sr = 44100
    dur = 1.0
    f = 440.0
    s = tone_from_omega(omega=2 * math.pi * f, duration=dur, sample_rate=sr,
                        amplitude=0.9, attack=0.0, release=0.0)
    spec = np.abs(np.fft.rfft(s.astype(np.float64)))
    freqs = np.fft.rfftfreq(len(s), d=1.0 / sr)
    peak_f = freqs[int(np.argmax(spec))]
    # 1 Hz of FFT bin width with dur=1 → tolerance ~1 Hz
    assert abs(peak_f - f) <= 1.5


def test_tone_zero_duration():
    s = tone_from_omega(omega=2 * math.pi * 440.0, duration=0.0)
    assert s.size == 0


def test_tone_decay_envelope_strictly_decaying():
    """With non-zero decay the envelope of |s| must decline from start to end."""
    sr = 8000  # short for speed
    dur = 1.0
    s = tone_from_omega(omega=2 * math.pi * 200.0, duration=dur,
                        sample_rate=sr, amplitude=1.0, decay=0.2,
                        attack=0.0, release=0.0)
    n = len(s)
    early = float(np.max(np.abs(s[: n // 8])))
    late = float(np.max(np.abs(s[-n // 8:])))
    assert late < early


# ---------------------------------------------------------------------------
# write_wav round-trip
# ---------------------------------------------------------------------------

def test_write_wav_mono(tmp_path):
    sr = 22050
    s = tone_from_omega(omega=2 * math.pi * 440.0, duration=0.1, sample_rate=sr)
    out = tmp_path / "tone.wav"
    write_wav(str(out), s, sample_rate=sr)
    assert out.exists()
    sr2, data = wavfile.read(str(out))
    assert sr2 == sr
    assert data.dtype == np.int16
    assert data.ndim == 1
    assert len(data) == len(s)


def test_write_wav_stereo(tmp_path):
    sr = 22050
    left = tone_from_omega(omega=2 * math.pi * 440.0, duration=0.1, sample_rate=sr)
    right = tone_from_omega(omega=2 * math.pi * 660.0, duration=0.1, sample_rate=sr)
    stereo = to_stereo(left, right)
    out = tmp_path / "stereo.wav"
    write_wav(str(out), stereo, sample_rate=sr)
    sr2, data = wavfile.read(str(out))
    assert sr2 == sr
    assert data.shape[1] == 2


def test_mix_normalizes():
    a = np.ones(10, dtype=np.float32) * 0.7
    b = np.ones(10, dtype=np.float32) * 0.7
    m = mix(a, b)
    assert math.isclose(float(np.max(np.abs(m))), 1.0, rel_tol=1e-6)


# ---------------------------------------------------------------------------
# Lepton tower
# ---------------------------------------------------------------------------

def test_lepton_tower_creates_stereo_wav(tmp_path):
    out = tmp_path / "lepton.wav"
    path, entries = lepton_tower_audio(duration=0.5, sample_rate=22050,
                                       output_path=str(out))
    assert path == str(out)
    assert out.exists()
    sr, data = wavfile.read(path)
    assert sr == 22050
    assert data.shape[1] == 2  # stereo: matter/antimatter
    # All three leptons present
    names = [e.name for e in entries]
    assert names == ["e", "mu", "tau"]
    # Ratios preserved up to octave folding
    raw_ratio = entries[1].mass_mev / entries[0].mass_mev
    audible_ratio = entries[1].audio_hz / entries[0].audio_hz
    log2_diff = math.log2(raw_ratio / audible_ratio)
    assert math.isclose(log2_diff, round(log2_diff), abs_tol=1e-9)


def test_lepton_tower_frequency_content(tmp_path):
    """FFT of the rendered chord should show energy near each lepton's audible f.

    Use a sample rate high enough that *every* tower entry sits well below
    Nyquist (the tau, folded into the audible window, lands near 13.7 kHz
    so we need sr >= ~30 kHz to leave headroom).
    """
    out = tmp_path / "lepton.wav"
    sr_test = 44100
    path, entries = lepton_tower_audio(duration=2.0, sample_rate=sr_test,
                                       output_path=str(out))
    sr, data = wavfile.read(path)
    assert sr == sr_test
    mono = data[:, 0].astype(np.float64) / INT16_MAX
    spec = np.abs(np.fft.rfft(mono))
    freqs = np.fft.rfftfreq(len(mono), d=1.0 / sr)
    nyquist = sr / 2.0
    for e in entries:
        # Skip if a tower entry happens to land above Nyquist for this sr
        if e.audio_hz >= nyquist:
            continue
        bin_target = int(round(e.audio_hz * len(mono) / sr))
        # local-max search ±50 bins, clamped to spectrum bounds
        lo = max(0, bin_target - 50)
        hi = min(len(spec), bin_target + 51)
        assert hi > lo, f"empty spectrum slice for {e.name}"
        peak = lo + int(np.argmax(spec[lo:hi]))
        peak_f = freqs[peak]
        assert abs(peak_f - e.audio_hz) <= 5.0, \
            f"lepton {e.name}: target {e.audio_hz} Hz, got {peak_f} Hz"


# ---------------------------------------------------------------------------
# Quark tower
# ---------------------------------------------------------------------------

def test_quark_tower_creates_stereo_wav(tmp_path):
    out = tmp_path / "quark.wav"
    path, entries = quark_tower_audio(duration=0.5, sample_rate=22050,
                                      output_path=str(out))
    assert out.exists()
    sr, data = wavfile.read(path)
    assert sr == 22050
    assert data.shape[1] == 2
    names = [e.name for e in entries]
    assert names == ["u", "c", "t", "d", "s", "b"]
    # Each entry's audible frequency must be within audible range
    for e in entries:
        assert AUDIBLE_LOW <= e.audio_hz <= AUDIBLE_HIGH


def test_quark_audible_ratios_match_mass_ratios_mod_octave():
    _, entries = quark_tower_audio(duration=0.2, sample_rate=22050,
                                   output_path=os.devnull
                                   if False else str(Path(tempfile.gettempdir()) / "tmp_q.wav"))
    by = {e.name: e for e in entries}
    for a, b in [("c", "u"), ("t", "c"), ("s", "d"), ("b", "s")]:
        raw = by[a].mass_mev / by[b].mass_mev
        aud = by[a].audio_hz / by[b].audio_hz
        log2_diff = math.log2(raw / aud)
        assert math.isclose(log2_diff, round(log2_diff), abs_tol=1e-9)


# ---------------------------------------------------------------------------
# Cone-bouncing sweep
# ---------------------------------------------------------------------------

def test_cone_bouncing_creates_stereo(tmp_path):
    out = tmp_path / "cb.wav"
    path, info = cone_bouncing_audio(duration=1.0, sample_rate=22050,
                                     output_path=str(out))
    assert out.exists()
    sr, data = wavfile.read(path)
    assert sr == 22050
    assert data.shape[1] == 2


def test_cone_bouncing_drag_upshifts_frequency(tmp_path):
    """End frequency = √(1+(γ̃_max/2)²) × start frequency."""
    out = tmp_path / "cb.wav"
    gamma_max = 1.5
    path, info = cone_bouncing_audio(duration=1.0, sample_rate=22050,
                                     gamma_over_omega_max=gamma_max,
                                     output_path=str(out))
    expected = math.sqrt(1.0 + (gamma_max / 2.0) ** 2)
    assert math.isclose(info["frequency_shift_ratio"], expected, rel_tol=1e-9)
    assert math.isclose(info["f_start_hz"], ELECTRON_AUDIO_HZ, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Cosmology
# ---------------------------------------------------------------------------

def test_cosmology_creates_stereo(tmp_path):
    out = tmp_path / "cos.wav"
    path, info = cosmology_evolution_audio(duration=1.0, sample_rate=22050,
                                           output_path=str(out))
    assert out.exists()
    sr, data = wavfile.read(path)
    assert sr == 22050
    assert data.shape[1] == 2


def test_cosmology_sigma_decays_and_pitch_rises(tmp_path):
    out = tmp_path / "cos.wav"
    _, info = cosmology_evolution_audio(duration=1.0, sample_rate=22050,
                                        output_path=str(out))
    # σ goes 1/2 → 0
    assert math.isclose(info["sigma_start"], 0.5, rel_tol=1e-9)
    assert info["sigma_end"] < info["sigma_start"]
    assert info["sigma_end"] >= 0.0
    # pitch rises (saturated → relaxed substrate sounds higher)
    assert info["f_end_hz"] > info["f_start_hz"]


# ---------------------------------------------------------------------------
# render_all
# ---------------------------------------------------------------------------

def test_render_all_produces_four_files(tmp_path, monkeypatch):
    """Patch the audio dir to a tmp path and confirm all four files appear."""
    from src.stiff_medium import substrate_audio as sa

    monkeypatch.setattr(sa, "_audio_dir", lambda: tmp_path)
    paths = sa.render_all(sample_rate=22050)
    assert set(paths.keys()) == {"lepton", "quark", "cone_bouncing", "cosmology"}
    for p in paths.values():
        assert Path(p).exists()
        sr, data = wavfile.read(p)
        assert sr == 22050
        assert data.dtype == np.int16
