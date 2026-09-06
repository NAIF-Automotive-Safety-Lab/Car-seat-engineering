"""Tests for cube_dm_cross_section module."""

import math

import pytest

from src.stiff_medium.cube_dm_cross_section import (
    CubeDMCrossSection,
    LIMIT_LZ_2024_CM2,
    LIMIT_PANDAX_2024_CM2,
    LIMIT_XENONNT_PROJ_CM2,
)


@pytest.fixture
def cube():
    return CubeDMCrossSection()


def test_default_mass(cube):
    assert cube.m_dm_gev == pytest.approx(27.5)
    assert cube.m_dm_kg > 0.0


def test_sigma_si_extremely_small(cube):
    """Per-nucleon SI cross-section must be << 1e-80 cm^2 (gravity-only)."""
    s = cube.sigma_si_nucleon_cm2()
    assert s > 0.0
    assert s < 1.0e-80, f"sigma_SI = {s:.3e} cm^2 too large for gravity-only"


def test_sigma_sd_smaller_than_si(cube):
    """SD piece is suppressed by (v/c)^2."""
    si = cube.sigma_si_nucleon_cm2()
    sd = cube.sigma_sd_nucleon_cm2()
    assert sd < si
    ratio = sd / si
    assert ratio == pytest.approx((cube.v_rel / 299792458.0) ** 2, rel=1e-9)


def test_sigma_below_all_current_limits(cube):
    s = cube.sigma_si_nucleon_cm2()
    assert s < LIMIT_LZ_2024_CM2
    assert s < LIMIT_PANDAX_2024_CM2
    assert s < LIMIT_XENONNT_PROJ_CM2
    # By many orders of magnitude
    assert s / LIMIT_LZ_2024_CM2 < 1.0e-30


def test_event_rate_at_LZ_is_zero(cube):
    """LZ active mass ~5.5 t, full ~7 t. Expect <<1 event in any human
    timescale for gravity-only DM."""
    rate = cube.event_rate(detector_mass_kg=7.0e3, exposure_yr=10.0)
    assert rate >= 0.0
    assert rate < 1.0e-20


def test_discovery_horizon_impossibly_large(cube):
    """Detector mass needed for 1 event/yr should dwarf planetary scales."""
    mass_kg = cube.discovery_horizon(target_events=1.0, exposure_yr=1.0)
    assert mass_kg > 1.0e30  # >> Earth (6e24 kg)


def test_verify_null_consistency(cube):
    out = cube.verify_null_consistency()
    assert out["consistent_with_nulls"] is True
    assert out["below_LZ_by_factor"] > 1.0e30
    assert out["below_nu_floor_by"] > 1.0e30


def test_compare_to_wimp_scenarios(cube):
    out = cube.compare_to_wimp_scenarios()
    assert "substrate_cube_DM" in out
    assert out["substrate_cube_DM"]["detectable_now"] is False
    assert out["SUSY_neutralino_typical"]["ratio_to_substrate"] > 1.0e30


def test_coherent_A2_enhancement(cube):
    """Nuclear cross-section enhanced over per-nucleon by ~A^2 * (mu_nuc/mu_N)^2."""
    per_nuc = cube.sigma_si_nucleon_cm2()
    nuc = cube.sigma_si_nucleus_cm2()
    enhancement = nuc / per_nuc
    expected = cube.target_A**2 * (cube.mu_nucleus / cube.mu_nucleon) ** 2
    assert enhancement == pytest.approx(expected, rel=1e-9)
    assert enhancement > 1.0e4  # A=131 gives >> 1


def test_event_rate_scales_linearly(cube):
    r1 = cube.event_rate(detector_mass_kg=1.0e3, exposure_yr=1.0)
    r2 = cube.event_rate(detector_mass_kg=2.0e3, exposure_yr=1.0)
    r3 = cube.event_rate(detector_mass_kg=1.0e3, exposure_yr=2.0)
    if r1 > 0.0:
        assert r2 / r1 == pytest.approx(2.0, rel=1e-6)
        assert r3 / r1 == pytest.approx(2.0, rel=1e-6)


def test_summary_runs(cube):
    s = cube.summary()
    assert "Cube DM" in s
    assert "sigma_SI" in s


def test_alternative_target():
    """Try a lighter target (Ar-40)."""
    c = CubeDMCrossSection(target_A=40, target_amu=39.948)
    assert c.sigma_si_nucleus_cm2() > 0.0
    rate = c.event_rate(detector_mass_kg=1.0e3, exposure_yr=1.0)
    assert math.isfinite(rate)
