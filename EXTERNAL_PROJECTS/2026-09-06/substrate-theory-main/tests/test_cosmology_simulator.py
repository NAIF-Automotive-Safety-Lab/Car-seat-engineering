"""Tests for the substrate cosmology simulator."""

import numpy as np
import pytest

from src.stiff_medium.cosmology_simulator import (
    H0_OBS,
    OMEGA_LAMBDA_OBS,
    RHO_LAMBDA_OBS,
    SIGMA8_OBS,
    SUM_MNU_OBS_MEV,
    SubstrateCosmologySimulator,
    SubstrateParams,
)


@pytest.fixture(scope="module")
def sim():
    s = SubstrateCosmologySimulator()
    s.run()
    return s


def test_runs_and_produces_history(sim):
    h = sim.history
    for k in ("t", "sigma", "a", "H", "rho_total", "T", "sum_mnu_ev"):
        assert k in h, f"missing key {k}"
        assert len(h[k]) == sim.p.n_steps


def test_sigma_never_exceeds_half(sim):
    sigma = sim.history["sigma"]
    assert np.max(sigma) <= 0.5 + 1e-12, f"sigma exceeded 0.5: max={np.max(sigma)}"
    assert sigma[0] == pytest.approx(0.5, abs=1e-12)


def test_scale_factor_monotone(sim):
    a = sim.history["a"]
    assert np.all(np.diff(a) >= -1e-15), "scale factor must be monotone"
    assert a[-1] > a[0]


def test_H0_within_5pct(sim):
    h0 = sim.predict_H0()
    assert abs(h0 - H0_OBS) / H0_OBS < 0.05, f"H0={h0}"


def test_sum_mnu_within_5pct(sim):
    m = sim.predict_sum_mnu_ev()
    assert abs(m - SUM_MNU_OBS_MEV) / SUM_MNU_OBS_MEV < 0.05, f"sum_mnu={m}"


def test_sigma8_within_5pct(sim):
    s8 = sim.predict_sigma8()
    assert abs(s8 - SIGMA8_OBS) / SIGMA8_OBS < 0.05, f"sigma8={s8}"


def test_omega_lambda_within_5pct(sim):
    ol = sim.predict_omega_lambda()
    assert abs(ol - OMEGA_LAMBDA_OBS) / OMEGA_LAMBDA_OBS < 0.05, f"OmegaL={ol}"


def test_rho_lambda_close(sim):
    r = sim.predict_rho_lambda()
    assert abs(r - RHO_LAMBDA_OBS) / RHO_LAMBDA_OBS < 0.05, f"rhoL={r}"


def test_temperature_profile_monotone(sim):
    prof = sim.temperature_profile()
    z = prof["z"]
    T = prof["T"]
    # T decreases with cosmic time => increases with z
    assert T[0] >= T[-1] - 1e-9


def test_density_perturbation_reports_gap(sim):
    rep = sim.density_perturbation_amplitude()
    assert "gap_factor" in rep
    assert "status" in rep
    assert rep["gap_factor"] > 1.0
    assert rep["log10_gap"] > 5.0  # gap is large; B3 reports ~10^13


def test_summary_keys():
    s = SubstrateCosmologySimulator().summary()
    for top in ("predicted", "observed", "rel_err"):
        assert top in s
    for k in ("H0_km_s_Mpc", "sum_mnu_eV", "sigma_8", "Omega_Lambda", "rho_Lambda_kg_m3"):
        assert k in s["predicted"]
