"""Tests for stiff_medium.black_hole_v2 — substrate BH physics."""

from __future__ import annotations

import math

import pytest

from stiff_medium.black_hole_v2 import (
    BlackHoleV2,
    M_SUN,
    M_M87,
    M_SGR_A,
    SIGMA_HORIZON,
    default_bh,
)


@pytest.fixture
def bh() -> BlackHoleV2:
    return default_bh()


# ---------------------------------------------------------------------------
# (a) Schwarzschild radius for the Sun
# ---------------------------------------------------------------------------


def test_schwarzschild_radius_sun(bh: BlackHoleV2) -> None:
    rs = bh.schwarzschild_radius(M_SUN)
    # Textbook value is 2.953 km
    assert rs == pytest.approx(2953.0, rel=2e-3)


# ---------------------------------------------------------------------------
# (b) Hawking temperature for solar-mass BH ~6e-8 K
# ---------------------------------------------------------------------------


def test_hawking_temperature_solar_mass(bh: BlackHoleV2) -> None:
    T = bh.hawking_temperature(M_SUN)
    # Canonical value ≈ 6.17e-8 K
    assert T == pytest.approx(6.17e-8, rel=5e-3)


# ---------------------------------------------------------------------------
# (c) Ringdown ~250 Hz for ~60 M_sun BH-BH merger
# ---------------------------------------------------------------------------


def test_ringdown_60_solar_masses(bh: BlackHoleV2) -> None:
    # ~30 + 30 -> ~60 M_sun system
    rd = bh.merger_GW(30.0 * M_SUN, 30.0 * M_SUN, radiation_efficiency=0.05)
    # Schwarzschild (2,2,0) QNM at ~57 M_sun remnant ≈ 210 Hz; spinning
    # remnants (e.g. GW150914) ringdown nearer 250 Hz. Allow 25% band.
    assert rd.f_ringdown_Hz == pytest.approx(250.0, rel=0.25)
    # Damping τ a few ms (GW150914 ≈ 4 ms)
    assert 1.0e-3 < rd.tau_damping_s < 1.0e-2
    assert rd.Q_quality > 1.0


# ---------------------------------------------------------------------------
# Additional structural checks
# ---------------------------------------------------------------------------


def test_horizon_substrate_sigma_one_half(bh: BlackHoleV2) -> None:
    h = bh.horizon_substrate(M_SUN)
    assert h.sigma_value == SIGMA_HORIZON
    assert h.r_h_m == pytest.approx(bh.schwarzschild_radius(M_SUN))
    assert h.area_m2 == pytest.approx(4 * math.pi * h.r_h_m ** 2)


def test_cone_tilt_at_horizon_is_90deg(bh: BlackHoleV2) -> None:
    rs = bh.schwarzschild_radius(M_SUN)
    assert bh.cone_tilt_at_radius(rs, M_SUN) == pytest.approx(math.pi / 2.0)
    # Inside is also fully tilted (saturated)
    assert bh.cone_tilt_at_radius(0.5 * rs, M_SUN) == pytest.approx(math.pi / 2.0)


def test_cone_tilt_far_field_is_small(bh: BlackHoleV2) -> None:
    rs = bh.schwarzschild_radius(M_SUN)
    tilt = bh.cone_tilt_at_radius(1e6 * rs, M_SUN)
    assert 0.0 < tilt < 1e-2


def test_cone_tilt_monotone_in_radius(bh: BlackHoleV2) -> None:
    rs = bh.schwarzschild_radius(M_SUN)
    tilts = [bh.cone_tilt_at_radius(r, M_SUN)
             for r in (1.01 * rs, 2 * rs, 10 * rs, 1000 * rs)]
    assert tilts == sorted(tilts, reverse=True)


def test_kerr_horizons_obey_cosmic_censorship(bh: BlackHoleV2) -> None:
    res = bh.kerr_metric(M_SUN, a=0.7)
    assert res["r_plus_m"] > res["r_minus_m"] > 0
    assert res["r_plus_m"] < res["r_s_m"]
    with pytest.raises(ValueError):
        bh.kerr_metric(M_SUN, a=1.2)


def test_kerr_extremal_collapses_horizons(bh: BlackHoleV2) -> None:
    res = bh.kerr_metric(M_SUN, a=1.0)
    assert res["r_plus_m"] == pytest.approx(res["r_minus_m"])
    assert res["r_plus_m"] == pytest.approx(G_times_M_over_c2 := res["r_g_m"])


def test_bekenstein_hawking_entropy_scales_as_M_squared(bh: BlackHoleV2) -> None:
    S1 = bh.bekenstein_hawking_entropy(M_SUN)
    S2 = bh.bekenstein_hawking_entropy(2.0 * M_SUN)
    assert S2 / S1 == pytest.approx(4.0, rel=1e-10)
    # Solar-mass BH entropy ~1.05e77 (in units of k_B)
    assert S1 == pytest.approx(1.05e77, rel=5e-2)


def test_eht_M87_ring_diameter_microarcsec(bh: BlackHoleV2) -> None:
    img = bh.eht_image_M87()
    # EHT measured shadow diameter ≈ 42 µas (2019)
    assert img.ring_diameter_uas == pytest.approx(40.0, rel=0.20)
    assert img.r_photon_sphere_m == pytest.approx(1.5 * img.r_s_m)
    assert img.sigma_at_ring == SIGMA_HORIZON


def test_sgr_a_radius_about_12_million_meters(bh: BlackHoleV2) -> None:
    rs = bh.schwarzschild_radius(M_SGR_A)
    # 4.3e6 M_sun → r_s ≈ 1.27e10 m
    assert rs == pytest.approx(1.27e10, rel=5e-2)


def test_information_paradox_substrate_resolution(bh: BlackHoleV2) -> None:
    info = bh.information_paradox_substrate()
    assert info["singularity"] is False
    assert "Planck area" in info["encoding_density"]
    assert info["compatible_with_ER_EPR"] is True


def test_report_contains_all_sections(bh: BlackHoleV2) -> None:
    rep = bh.report()
    assert rep["sigma_cap"] == 0.5
    assert len(rep["masses_table"]) == 3
    assert 150.0 < rep["GW150914_ringdown_Hz"] < 350.0
    assert "info_paradox" in rep


def test_invalid_inputs_raise(bh: BlackHoleV2) -> None:
    with pytest.raises(ValueError):
        bh.schwarzschild_radius(-1.0)
    with pytest.raises(ValueError):
        bh.hawking_temperature(0.0)
    with pytest.raises(ValueError):
        bh.cone_tilt_at_radius(-1.0, M_SUN)
    with pytest.raises(ValueError):
        BlackHoleV2(sigma_cap=0.9)
