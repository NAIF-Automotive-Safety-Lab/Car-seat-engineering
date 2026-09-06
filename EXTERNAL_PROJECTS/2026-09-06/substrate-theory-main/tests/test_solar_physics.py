"""Tests for stiff_medium.solar_physics."""

import math

import pytest

from stiff_medium.solar_physics import (
    FLUX_8B_OBS,
    FLUX_PP_OBS,
    L_SUN_OBS,
    T_EFF_OBS,
    SolarPhysics,
)


@pytest.fixture
def sun():
    return SolarPhysics()


def test_luminosity_within_1pct(sun):
    L = sun.solar_luminosity()
    assert math.isclose(L, L_SUN_OBS, rel_tol=0.01), (
        f"L_sun {L:.3e} vs obs {L_SUN_OBS:.3e}"
    )


def test_effective_temperature_within_1pct(sun):
    T = sun.effective_temperature()
    assert math.isclose(T, T_EFF_OBS, rel_tol=0.01), (
        f"T_eff {T:.2f} vs obs {T_EFF_OBS:.2f}"
    )


def test_pp_neutrino_flux_within_30pct(sun):
    flux = sun.solar_neutrino_flux_pp()
    rel = abs(flux - FLUX_PP_OBS) / FLUX_PP_OBS
    assert rel < 0.30, f"pp flux {flux:.3e} vs obs {FLUX_PP_OBS:.3e}, rel={rel:.2%}"


def test_8B_neutrino_flux_matches_reference(sun):
    flux = sun.solar_neutrino_flux_8B()
    assert math.isclose(flux, FLUX_8B_OBS, rel_tol=1e-9)


def test_core_temperature_order_of_magnitude(sun):
    T = sun.core_temperature()
    # SSM core ~1.57e7 K; allow factor-of-2 since this is an analytic estimate
    assert 7e6 < T < 3e7, f"T_core estimate {T:.2e} K outside reasonable band"


def test_pp_chain_rate_consistent_with_luminosity(sun):
    pp = sun.pp_chain_rate()
    Q_J = pp["Q_deposited_MeV"] * 1.602176634e-13
    L_recovered = pp["cycles_per_second"] * Q_J
    assert math.isclose(L_recovered, sun.L_obs, rel_tol=1e-6)


def test_helioseismology_5min_band(sun):
    h = sun.helioseismology_oscillations()
    # 5-min p-modes: period should be 4-7 min
    assert 4.0 < h["period_minutes"] < 7.0
    assert 100 < h["large_separation_uHz"] < 200


def test_solar_wind_speeds_in_band(sun):
    w = sun.solar_wind_parameters()
    assert 300 < w["v_slow_km_s"] < 500
    assert 600 < w["v_fast_km_s"] < 900
    assert w["mass_loss_rate_Msun_per_yr"] < 1e-12


def test_coronal_heating_ratio(sun):
    c = sun.coronal_heating_problem()
    # corona ~250-350x hotter than photosphere
    assert 100 < c["T_ratio"] < 1000
    assert c["SM_open_problem"] is True


def test_compare_dict_has_all_fields(sun):
    cmp = sun.compare_to_observations()
    expected = {
        "luminosity_W",
        "T_eff_K",
        "T_core_K",
        "nu_pp_flux_cm2_s",
        "nu_8B_flux_cm2_s",
    }
    assert set(cmp.keys()) == expected
    for d in cmp.values():
        assert {"calc", "obs", "pct_diff"} <= set(d.keys())


def test_report_returns_string(sun):
    sun.core_temperature()
    r = sun.report()
    assert isinstance(r, str)
    assert "B3 SOLAR PHYSICS REPORT" in r
    assert "p-p chain" in r
    assert "Coronal heating" in r
