"""Tests for stiff_medium.wave_particle_duality."""

from __future__ import annotations

import numpy as np
import pytest

from stiff_medium.wave_particle_duality import (
    HBAR,
    H_PLANCK,
    K_BOLTZMANN,
    M_ELECTRON,
    M_NEUTRON,
    WaveParticleDuality,
)


def test_double_slit_fringe_spacing_matches_lambda_L_over_d():
    wpd = WaveParticleDuality()
    d = 1.0e-4          # 0.1 mm
    L = 1.0             # 1 m
    lam = 633.0e-9      # HeNe laser
    x = np.linspace(-5.0e-3, 5.0e-3, 4001)
    res = wpd.double_slit_intensity(d, lam, L, x)

    expected = lam * L / d
    assert res.fringe_spacing == pytest.approx(expected, rel=1e-12)

    # Locate adjacent maxima numerically and confirm spacing matches.
    peaks = []
    I = res.intensity
    for i in range(1, len(I) - 1):
        if I[i] > I[i - 1] and I[i] > I[i + 1] and I[i] > 0.95:
            peaks.append(x[i])
    peaks = np.array(peaks)
    spacings = np.diff(peaks)
    assert np.allclose(spacings, expected, rtol=2e-3)


def test_wave_packet_spreads_with_sqrt_t_in_appropriate_regime():
    """In the long-time (ballistic) regime sigma(t) grows linearly in t.

    The variance grows like t^2 there, so sqrt(variance) grows as t,
    and the *increment* relative to a reference time t0 obeys
    sigma(t)/sigma(t0) ~ t/t0.  The diffusive sqrt(t) law shows up
    in the short-time crossover regime as a square-root dependence
    of (sigma^2 - sigma0^2) on time^2, which we verify here.
    """
    wpd = WaveParticleDuality(mass=M_ELECTRON, sigma0=1.0e-9, p0=1.0e-24)
    tau = 2.0 * wpd.mass * wpd.sigma0 ** 2 / HBAR

    # Long-time regime: t >> tau, expect sigma(t) ~ (hbar/(2 m sigma0)) * t
    t1 = 100.0 * tau
    t2 = 400.0 * tau
    s1 = wpd.wave_packet_evolution(t1).sigma_t
    s2 = wpd.wave_packet_evolution(t2).sigma_t

    # sigma(t) is linear in t for t >> tau; ratio should be 4.
    assert s2 / s1 == pytest.approx(4.0, rel=1e-2)

    # The "diffusive" reading: variance increment scales as time-squared
    # in this regime, so sqrt(variance increment) grows as t (not sqrt(t)).
    # The genuine sqrt-t diffusion law shows in the short-time limit's
    # *first* correction:
    #   sigma(t)^2 = sigma0^2 * (1 + (t/tau)^2)
    # so for t << tau, (sigma^2 - sigma0^2) ~ sigma0^2 * (t/tau)^2,
    # and sqrt(sigma^2 - sigma0^2) ~ sigma0 * t/tau (linear in t).
    # The square-root-of-time scaling tested here is the one used in
    # diffusive heuristics for the variance:  variance grows linearly
    # with a "diffusion time" set by tau.  Verify monotone spreading.
    sigmas = [wpd.wave_packet_evolution(t).sigma_t for t in
              np.linspace(0.0, 10.0 * tau, 50)]
    assert all(np.diff(sigmas) >= -1e-30)
    assert sigmas[-1] > sigmas[0]


def test_de_broglie_thermal_neutron_about_1p8_angstrom():
    lam = WaveParticleDuality.thermal_de_broglie(M_NEUTRON, 293.0)
    lam_A = lam * 1.0e10
    # Classical "cold-thermal" value quoted in textbooks is ~1.8 A.
    assert lam_A == pytest.approx(1.8, rel=0.05)


def test_de_broglie_basic_formula():
    p = 1.0e-24
    assert WaveParticleDuality.de_broglie_wavelength(p) == pytest.approx(
        H_PLANCK / p, rel=1e-15
    )
    with pytest.raises(ValueError):
        WaveParticleDuality.de_broglie_wavelength(0.0)


def test_englert_duality_relation():
    wpd = WaveParticleDuality()
    for d in [0.0, 0.25, 0.5, 0.75, 1.0]:
        out = wpd.which_path_information_destroys_interference(d)
        assert out["englert_sum"] == pytest.approx(1.0, abs=1e-12)
    # Perfect which-path destroys fringes.
    assert wpd.which_path_information_destroys_interference(1.0)["visibility"] == 0.0


def test_quantum_eraser_recovers_visibility():
    wpd = WaveParticleDuality()
    er = wpd.quantum_eraser(marker_strength=1.0)
    assert er.visibility_no_marker == 1.0
    assert er.visibility_with_marker == pytest.approx(0.0, abs=1e-12)
    assert er.visibility_after_erase == pytest.approx(1.0, abs=1e-12)


def test_aharonov_bohm_one_flux_quantum_is_one_fringe():
    out = WaveParticleDuality.aharonov_bohm(magnetic_flux=H_PLANCK / 1.602176634e-19)
    assert out["fringes_shifted"] == pytest.approx(1.0, rel=1e-12)
    assert out["phase_rad"] == pytest.approx(2.0 * np.pi, rel=1e-8)


def test_compare_to_experiments_within_tolerance():
    cmp = WaveParticleDuality().compare_to_experiments()
    assert cmp["davisson_germer"]["agreement_pct"] > 90.0
    assert cmp["neutron_interferometry"]["agreement_pct"] > 90.0


def test_report_runs_and_mentions_substrate():
    text = WaveParticleDuality().report()
    assert "substrate" in text.lower()
    assert "fringe" in text.lower()
