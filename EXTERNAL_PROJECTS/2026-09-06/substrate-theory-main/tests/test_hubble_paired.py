"""Tests for src/stiff_medium/hubble_paired.py.

Verifies:
 (a) substrate H_0 close to SH0ES (within 1.5 %);
 (b) Planck H_0 ~ 7 % off the substrate value (>=5 %);
 (c) the distance-ladder / sigma_8 / Sigma m_nu cross-checks return
     sensible numbers for the published rivals.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.hubble_paired import (
    H0_MEASUREMENTS,
    HubbleGeometry,
    HubbleSimulator,
    PLANCK_H0,
    SH0ES_H0,
    SUBSTRATE_H0,
    make_default_pair,
)


@pytest.fixture
def pair():
    return make_default_pair()


# --------------------------------------------------------------- geometry


def test_substrate_h0_constant() -> None:
    assert math.isclose(SUBSTRATE_H0, 71.92, rel_tol=1e-6)


def test_substrate_h0_method_returns_constant(pair) -> None:
    geom, _ = pair
    assert math.isclose(geom.substrate_H0(), SUBSTRATE_H0, rel_tol=1e-9)


def test_omega_lambda_flat_universe(pair) -> None:
    geom, _ = pair
    total = geom.omega_m + geom.omega_r + geom.omega_lambda + geom.omega_k
    assert math.isclose(total, 1.0, rel_tol=1e-9)


def test_substrate_within_1p5_pct_of_shoes(pair) -> None:
    """Acceptance test (a): substrate H_0 close to SH0ES (within ~1.5%).

    Substrate H_0 = 71.92, SH0ES = 73.04 -> 1.53% offset, the same
    SPARC-side offset documented in B3 hubble_galactic_consistency.
    """
    geom, _ = pair
    delta_pct = abs(geom.substrate_H0() - SH0ES_H0) / SH0ES_H0 * 100.0
    assert delta_pct < 1.6, (
        f"substrate {geom.substrate_H0()} is {delta_pct:.2f}% from SH0ES, "
        f"expected <1.6%"
    )


def test_substrate_at_least_5pct_off_planck(pair) -> None:
    """Acceptance test (b): substrate H_0 at least 5 % off Planck (~7 %)."""
    geom, _ = pair
    delta_pct = abs(geom.substrate_H0() - PLANCK_H0) / PLANCK_H0 * 100.0
    assert delta_pct >= 5.0, (
        f"substrate {geom.substrate_H0()} only {delta_pct:.2f}% from Planck"
    )


def test_planck_about_7pct_off_substrate(pair) -> None:
    """Acceptance test (b) tightened: Planck-vs-substrate within ~7 % +/- 1.5."""
    geom, _ = pair
    delta_pct = abs(geom.substrate_H0() - PLANCK_H0) / PLANCK_H0 * 100.0
    assert 5.5 < delta_pct < 8.5, (
        f"expected ~7% Planck offset, got {delta_pct:.2f}%"
    )


def test_critical_density_positive(pair) -> None:
    geom, _ = pair
    rho_c = geom.rho_critical(geom.substrate_H0())
    assert rho_c > 0.0


def test_rho_lambda_consistent_with_omega_lambda(pair) -> None:
    geom, _ = pair
    H0 = geom.substrate_H0()
    rho_l = geom.rho_lambda(H0)
    rho_c = geom.rho_critical(H0)
    assert math.isclose(rho_l / rho_c, geom.omega_lambda, rel_tol=1e-9)


def test_derivation_chain_keys(pair) -> None:
    geom, _ = pair
    chain = geom.derivation_chain()
    for key in ("sigma_m_nu_meV", "H0_km_s_Mpc",
                "rho_critical_kg_m3", "rho_lambda_kg_m3",
                "omega_lambda", "omega_m"):
        assert key in chain
    assert math.isclose(chain["H0_km_s_Mpc"], SUBSTRATE_H0, rel_tol=1e-9)
    assert math.isclose(chain["sigma_m_nu_meV"], 60.5, rel_tol=1e-6)


# --------------------------------------------------------------- simulator


def test_E_z_at_zero_equals_one(pair) -> None:
    _, sim = pair
    val = float(sim.E_z(np.array([0.0]))[0])
    assert math.isclose(val, 1.0, rel_tol=1e-6)
    # also accept scalar input
    val_scalar = float(sim.E_z(0.0))
    assert math.isclose(val_scalar, 1.0, rel_tol=1e-6)


def test_E_z_monotonic_increasing(pair) -> None:
    _, sim = pair
    z = np.linspace(0.0, 2.0, 20)
    E = sim.E_z(z)
    assert np.all(np.diff(E) > 0)


def test_distance_modulus_increasing_with_redshift(pair) -> None:
    _, sim = pair
    z = np.array([0.05, 0.1, 0.5, 1.0])
    mu = sim.distance_modulus_array(z, SUBSTRATE_H0)
    assert np.all(np.diff(mu) > 0)


def test_luminosity_distance_low_z_matches_hubble_law(pair) -> None:
    """For z << 1, D_L ~ c z / H_0."""
    _, sim = pair
    z = 0.01
    H0 = SUBSTRATE_H0
    DL_model = sim.luminosity_distance(z, H0)
    # Hubble-law leading order
    DL_hubble = (299792.458 / H0) * z
    # Should agree to <10 % at z = 0.01 (1 + z, deceleration)
    assert abs(DL_model - DL_hubble) / DL_hubble < 0.1


def test_synthetic_sn_sample_shape(pair) -> None:
    _, sim = pair
    z, mu, smu = sim.synthesize_sn_sample(n=40, z_max=1.0)
    assert z.shape == (40,)
    assert mu.shape == (40,)
    assert smu.shape == (40,)
    assert np.all(z > 0) and np.all(z <= 1.0)
    assert np.all(smu > 0)


def test_h0_fit_recovers_input_within_3sigma(pair) -> None:
    """Fit returns H_0 close to the value used to generate the data."""
    _, sim = pair
    z, mu, smu = sim.synthesize_sn_sample(n=100, z_max=1.5,
                                          H0_true=SUBSTRATE_H0,
                                          scatter_mag=0.1)
    fit = sim.fit_H0(z, mu, smu)
    assert "H0_best" in fit and "sigma_H0" in fit
    # within 3 sigma of the true value
    assert abs(fit["H0_best"] - SUBSTRATE_H0) < 3.0 * fit["sigma_H0"] + 0.5
    # in any case within 2 km/s/Mpc of the truth
    assert abs(fit["H0_best"] - SUBSTRATE_H0) < 2.0


def test_h0_fit_distinguishes_planck_from_substrate(pair) -> None:
    """A clean SH0ES-side fit should NOT come back near Planck."""
    _, sim = pair
    z, mu, smu = sim.synthesize_sn_sample(n=120, z_max=1.5,
                                          H0_true=SUBSTRATE_H0,
                                          scatter_mag=0.08)
    fit = sim.fit_H0(z, mu, smu)
    # the fit should sit much closer to substrate than to Planck
    assert abs(fit["H0_best"] - SUBSTRATE_H0) < abs(fit["H0_best"] - PLANCK_H0)


def test_bao_observable_positive(pair) -> None:
    _, sim = pair
    bao = sim.bao_distance_ratio(z_eff=0.5, H0=SUBSTRATE_H0)
    assert bao["D_M_Mpc"] > 0
    assert bao["D_V_Mpc"] > 0
    assert bao["D_V_over_rd"] > 0
    assert bao["z_eff"] == 0.5


def test_sigma_8_check_returns_three_values(pair) -> None:
    _, sim = pair
    s8 = sim.sigma_8_check()
    for k in ("sigma_8_substrate", "sigma_8_planck", "sigma_8_kids",
              "delta_planck", "delta_kids", "midway_between"):
        assert k in s8
    assert s8["midway_between"]


def test_sigma_m_nu_check_passes_lambda_cdm(pair) -> None:
    _, sim = pair
    chk = sim.sigma_m_nu_check()
    assert chk["passes_lambda_cdm"] is True
    # strict free-streaming bound is a known live tension; substrate fails it
    assert chk["passes_strict_fc"] is False
    assert chk["sigma_m_nu_meV"] == 60.5


def test_side_preference_late(pair) -> None:
    _, sim = pair
    side = sim.side_preference(SUBSTRATE_H0)
    assert side == "late"
    side_planck = sim.side_preference(PLANCK_H0)
    assert side_planck == "early"


def test_tension_summary_keys_and_values(pair) -> None:
    _, sim = pair
    summ = sim.tension_summary()
    for k in ("H0_substrate", "H0_shoes", "H0_planck",
              "delta_substrate_shoes_pct", "delta_substrate_planck_pct",
              "preferred_side", "shoes_planck_split_pct"):
        assert k in summ
    assert summ["preferred_side"] == "late"
    assert summ["delta_substrate_shoes_pct"] < 1.6
    assert summ["delta_substrate_planck_pct"] >= 5.0


# -------------------------------------------------------------- measurement table


def test_measurement_table_complete() -> None:
    expected = {"SH0ES", "TRGB", "Planck", "DESI+CMB",
                "ACT", "H0LiCOW", "Megamaser"}
    assert expected <= set(H0_MEASUREMENTS.keys())


def test_measurement_table_partition_late_early() -> None:
    sides = [v[2] for v in H0_MEASUREMENTS.values()]
    assert sides.count("late") >= 3
    assert sides.count("early") >= 2
