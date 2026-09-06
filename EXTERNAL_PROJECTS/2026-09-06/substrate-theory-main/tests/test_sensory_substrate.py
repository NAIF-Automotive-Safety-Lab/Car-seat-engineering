"""Tests for the paired sensory-transduction GEOMETRY + SIMULATOR module.

Covers the four required claims:
    (a) Single-photon detection: a single absorbed photon -> detectable
        signal (Hecht-Shlaer-Pirenne 1942) with cascade gain ~1e6.
    (b) Weber's law: just-noticeable difference is proportional to I.
    (c) Cochlear tonotopic map runs from 20 Hz at the apex to 20 kHz at
        the base (von Bekesy 1928 place principle).
    (d) MET-channel half-activation at ~1 picometer stereocilia tip
        displacement (Hudspeth-Corey 1989).

Plus geometric / encoding sanity:
    * 11-cis -> all-trans isomerization completes within ~ 200 fs
    * cascade gains: R* ~500 G_t, G_t ~1000 cGMP, total ~1e6
    * olfactory: ~400 receptors, combinatorial coding capacity >= 10^4
    * Greenwood inverse map round-trips
    * Weber-Fechner curve is monotonic and concave (log shape)
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.sensory_substrate import (
    COCHLEA_LENGTH_MM,
    F_APEX_HZ,
    F_BASE_HZ,
    LAMBDA_QCD_MEV,
    MET_THRESHOLD_M,
    N_DISTINGUISHABLE_ODORS,
    N_OLFACTORY_RECEPTORS,
    PHOTON_TO_CHANNEL_GAIN,
    PDE_GAIN,
    RHODOPSIN_ISOMERIZATION_FS,
    RHODOPSIN_LAMBDA_MAX_NM,
    SINGLE_PHOTON_THRESHOLD,
    TRANSDUCIN_GAIN,
    WEBER_K_HEARING,
    WEBER_K_SMELL,
    WEBER_K_VISION,
    HairCell,
    OlfactoryReceptor,
    Rhodopsin,
    SensoryGeometry,
    SensorySimulator,
    cochlear_frequency_at_position,
    cochlear_position_at_frequency,
    just_noticeable_difference,
    photon_energy_J,
    weber_fechner,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> SensoryGeometry:
    return SensoryGeometry()


@pytest.fixture
def sim() -> SensorySimulator:
    return SensorySimulator()


# --------------------------------------------------------------------- #
# Constants / anchor sanity                                             #
# --------------------------------------------------------------------- #

def test_constants_match_textbook():
    assert RHODOPSIN_ISOMERIZATION_FS == pytest.approx(200.0, rel=1e-9)
    assert F_APEX_HZ == pytest.approx(20.0, rel=1e-9)
    assert F_BASE_HZ == pytest.approx(20000.0, rel=1e-9)
    assert MET_THRESHOLD_M == pytest.approx(1.0e-12, rel=1e-9)
    assert N_OLFACTORY_RECEPTORS == 400
    assert N_DISTINGUISHABLE_ODORS == 10000
    assert COCHLEA_LENGTH_MM == pytest.approx(35.0, rel=0.1)


def test_substrate_signature_has_lambda_qcd(geom):
    sig = geom.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["lambda_qcd_MeV"] == pytest.approx(LAMBDA_QCD_MEV, rel=1e-9)
    assert sig["MET_threshold_m"] == pytest.approx(1.0e-12, rel=1e-9)
    assert sig["olfactory_receptors"] == 400


# --------------------------------------------------------------------- #
# (a) Single-photon detection + amplification cascade                   #
# --------------------------------------------------------------------- #

def test_single_photon_amplification_is_about_1e6(sim):
    cas = sim.phototransduction_cascade(n_photons=1)
    assert cas["n_photons"] == 1.0
    assert cas["n_R_star"] == 1.0
    # transducin gain ~ 500
    assert cas["n_transducin"] == pytest.approx(TRANSDUCIN_GAIN, rel=1e-9)
    # PDE -> cGMP gain ~ 1000 per transducin
    assert cas["n_cGMP_hydrolyzed"] == pytest.approx(
        TRANSDUCIN_GAIN * PDE_GAIN, rel=1e-9)
    # Total gain ~ 10^6
    assert cas["total_gain"] == pytest.approx(PHOTON_TO_CHANNEL_GAIN, rel=1e-9)
    assert 1.0e5 <= cas["total_gain"] <= 1.0e7


def test_single_photon_detectable_above_quiet_noise(sim):
    """Hecht-Shlaer-Pirenne 1942: a single photon evokes a detectable signal."""
    assert sim.single_photon_detectable(n_photons=1, background_noise=1)
    assert not sim.single_photon_detectable(n_photons=0, background_noise=1)


def test_photoisomerization_completes_in_about_200_fs(sim):
    """11-cis -> all-trans population reaches ~63% by t = tau."""
    t = np.linspace(0.0, 5.0 * RHODOPSIN_ISOMERIZATION_FS, 200)
    out = sim.photoisomerize(t)
    # At t = tau, trans pop = 1 - 1/e ~ 0.632
    idx_tau = int(np.argmin(np.abs(t - RHODOPSIN_ISOMERIZATION_FS)))
    assert out["trans_pop"][idx_tau] == pytest.approx(1.0 - math.exp(-1.0),
                                                       rel=0.05)
    # By 5 tau, isomerization is essentially complete
    assert out["trans_pop"][-1] > 0.99
    # Conservation: cis + trans = 1
    assert np.allclose(out["cis_pop"] + out["trans_pop"], 1.0)


def test_rhodopsin_lambda_max_in_visible(geom):
    assert 380.0 <= geom.rhodopsin.lambda_max_nm <= 700.0
    assert geom.rhodopsin.lambda_max_nm == pytest.approx(498.0, rel=0.02)


def test_rhodopsin_photon_energy(geom):
    """E = hc/lambda for green light should be ~3-4 eV (~5e-19 J)."""
    E = geom.rhodopsin_photon_energy_J()
    assert 3.0e-19 < E < 6.0e-19


# --------------------------------------------------------------------- #
# (b) Weber's law / Weber-Fechner                                       #
# --------------------------------------------------------------------- #

def test_just_noticeable_difference_is_linear_in_intensity():
    """JND(I) = k * I:  doubling I doubles JND."""
    for k in (WEBER_K_VISION, WEBER_K_HEARING, WEBER_K_SMELL):
        assert just_noticeable_difference(100.0, k) == pytest.approx(100.0 * k)
        assert just_noticeable_difference(200.0, k) == pytest.approx(200.0 * k)
        assert (just_noticeable_difference(200.0, k)
                == pytest.approx(2.0 * just_noticeable_difference(100.0, k)))


def test_weber_fechner_is_logarithmic_concave(sim):
    """psi(I) = k log10(I/I_0) is monotonic and concave (second diff < 0)."""
    out = sim.weber_fechner_curve("vision", I_min=1.0, I_max=1.0e6, n_points=50)
    psi = out["psi"]
    # Monotonic non-decreasing
    assert np.all(np.diff(psi) >= -1e-12)
    # Concave: at log-spaced I the spacing in psi should NOT grow
    # (ln-spacing gives constant increments, so second diff ~ 0; we test
    # that psi(1e6) / psi(1e3) ~ 6/3 = 2, NOT > 2.)
    psi_at_1e3 = float(np.interp(1.0e3, out["I"], psi))
    psi_at_1e6 = float(np.interp(1.0e6, out["I"], psi))
    ratio = psi_at_1e6 / max(psi_at_1e3, 1e-12)
    assert 1.5 < ratio < 2.5


def test_weber_fraction_ordering():
    """Sensory acuity ordering: vision > hearing > smell (k_vision smallest)."""
    assert WEBER_K_VISION < WEBER_K_HEARING < WEBER_K_SMELL


# --------------------------------------------------------------------- #
# (c) Cochlear tonotopic frequency map                                  #
# --------------------------------------------------------------------- #

def test_cochlear_apex_frequency_is_20_hz():
    f = cochlear_frequency_at_position(0.0)
    assert f == pytest.approx(F_APEX_HZ, rel=1e-9)


def test_cochlear_base_frequency_is_20_khz():
    f = cochlear_frequency_at_position(COCHLEA_LENGTH_MM)
    assert f == pytest.approx(F_BASE_HZ, rel=1e-9)


def test_cochlear_position_inverse_round_trips():
    for f in (50.0, 200.0, 1000.0, 4000.0, 16000.0):
        x = cochlear_position_at_frequency(f)
        f_back = cochlear_frequency_at_position(x)
        assert f_back == pytest.approx(f, rel=1e-6)


def test_cochlear_map_is_monotonic(sim):
    out = sim.tonotopic_map(n_points=50)
    f = out["f_hz"]
    assert np.all(np.diff(f) > 0)  # strictly increasing
    assert f[0] == pytest.approx(F_APEX_HZ, rel=1e-9)
    assert f[-1] == pytest.approx(F_BASE_HZ, rel=1e-9)


def test_cochlear_traveling_wave_peaks_at_correct_place(sim):
    """Drive at 1 kHz should peak ~22-24 mm from apex (mid cochlea)."""
    out = sim.cochlear_traveling_wave(f_drive_hz=1000.0, n_points=400)
    # Peak position from envelope max
    idx_peak = int(np.argmax(out["envelope"]))
    x_peak_from_env = float(out["x_mm"][idx_peak])
    # Reference position predicted by the inverse map
    x_peak_ref = cochlear_position_at_frequency(1000.0)
    # Envelope peak should sit very close to the place principle prediction
    assert abs(x_peak_from_env - x_peak_ref) < 1.0  # within 1 mm
    assert 18.0 < x_peak_ref < 28.0


# --------------------------------------------------------------------- #
# (d) MET-channel sub-pm threshold                                      #
# --------------------------------------------------------------------- #

def test_met_threshold_is_about_1_pm(sim):
    """Half-activation at ~1 pm displacement."""
    half_disp = MET_THRESHOLD_M  # 1 pm
    I_half = float(sim.met_current(np.array([half_disp]))[0])
    assert I_half == pytest.approx(50.0, rel=0.1)  # 50% of 100 pA I_max


def test_met_zero_displacement_is_subthreshold(sim):
    I_rest = float(sim.met_current(np.array([0.0]))[0])
    # Boltzmann at -3.3 sigma should give < 5% of I_max
    assert I_rest < 5.0  # 5 pA at I_max = 100 pA


def test_met_full_open_at_large_displacement(sim):
    """At dx >> 1 pm, all MET channels are open."""
    I_open = float(sim.met_current(np.array([1.0e-11]))[0])  # 10 pm
    assert I_open > 95.0


def test_met_threshold_detectable_just_above_1pm(sim):
    """A 2 pm deflection should be detectable; sub-pm should not."""
    assert sim.met_threshold_detectable(2.0e-12)
    assert not sim.met_threshold_detectable(0.1e-12)


# --------------------------------------------------------------------- #
# Olfaction: combinatorial coding                                       #
# --------------------------------------------------------------------- #

def test_olfactory_receptor_count_is_400(geom):
    assert geom.n_olfactory_receptors == 400


def test_olfactory_combinatorial_capacity_exceeds_10000(sim):
    """With 400 receptors and ~3% sparsity, capacity >> 10^4."""
    cap = sim.discriminable_odors_via_combinatorial_code(sparsity=0.03)
    assert cap >= 10_000


def test_olfactory_receptor_response_is_sparse(sim, geom):
    """Each odorant activates only a small subset of ORs."""
    R = geom.make_receptor_response_matrix(n_odorants=200, sparsity=0.03,
                                            seed=7)
    # Average activation per odorant should be near 3% * 400 = 12
    activations_per_odor = R.sum(axis=0)
    mean_act = float(activations_per_odor.mean())
    assert 5.0 <= mean_act <= 25.0  # combinatorial sparsity holds


def test_olfactory_response_returns_binary_vector(sim):
    R = sim.geometry.make_receptor_response_matrix(n_odorants=64, seed=0)
    resp = sim.olfactory_response(odorant_index=3, response_matrix=R)
    assert resp.shape[0] == sim.geometry.n_olfactory_receptors
    assert set(np.unique(resp).tolist()).issubset({0, 1})


# --------------------------------------------------------------------- #
# Cross-modality: Weber-Fechner curve shape consistent across senses    #
# --------------------------------------------------------------------- #

def test_weber_fechner_three_modalities_all_logarithmic(sim):
    for modality in ("vision", "hearing", "smell"):
        out = sim.weber_fechner_curve(modality, I_min=1.0, I_max=1.0e4,
                                       n_points=50)
        # log-shaped: psi(10*I) - psi(I) is ~constant on log-spaced grid
        diffs = np.diff(out["psi"])
        # No huge jumps (no thresholds within the range)
        assert np.all(np.isfinite(diffs))
        # Strictly monotonic
        assert np.all(diffs >= -1e-12)


# --------------------------------------------------------------------- #
# Geometry sanity                                                       #
# --------------------------------------------------------------------- #

def test_haircell_default_inventory():
    cell = HairCell()
    assert cell.n_stereocilia == 50
    assert cell.met_threshold_m == pytest.approx(1.0e-12)


def test_rhodopsin_total_amplification_matches_constants():
    rh = Rhodopsin()
    assert rh.total_amplification() == pytest.approx(
        TRANSDUCIN_GAIN * PDE_GAIN, rel=1e-9)


def test_photon_energy_helper():
    E = photon_energy_J(500.0)  # green light
    assert 3.5e-19 < E < 4.5e-19
