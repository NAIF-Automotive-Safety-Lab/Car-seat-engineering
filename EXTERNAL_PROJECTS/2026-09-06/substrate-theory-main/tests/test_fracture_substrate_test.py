"""Tests for the substrate-cap fracture-mechanics module.

Required claims:
    (a) All four r_p formulas (substrate, Irwin plane stress, Irwin plane
        strain, Dugdale) return strictly positive values for every reference
        material.
    (b) substrate_rp(K_I, sigma_y) is mathematically derivable from the
        substrate cap sigma <= 1/2 plus K_I = sigma_inf*sqrt(pi*a) and
        reduces to (1/(2*pi)) * (K_I/sigma_y)**2.
    (c) r_p scales as (K_I/sigma_y)**2 -- doubling K_I quadruples r_p,
        doubling sigma_y quarters r_p.
    (d) For at least one canonical material (Steel 4340), the substrate
        prediction is within an order of magnitude of Irwin plane-stress
        (in fact: equal to it, by derivation).
    (e) compare_predictions and run_test populate the expected schema and
        identify Irwin plane-stress as the closest standard match.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.fracture_substrate_test import (
    MATERIALS,
    compare_predictions,
    dugdale_rp,
    irwin_plane_strain_rp,
    irwin_plane_stress_rp,
    run_test,
    substrate_rp,
)


# --------------------------------------------------------------------------- #
# Claim (a): all four formulas positive on every reference material           #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", list(MATERIALS.keys()))
def test_all_predictions_positive(name: str) -> None:
    props = MATERIALS[name]
    K = props["K_IC"]
    sy = props["sigma_y"]

    assert substrate_rp(K, sy) > 0.0
    assert irwin_plane_stress_rp(K, sy) > 0.0
    assert irwin_plane_strain_rp(K, sy) > 0.0
    assert dugdale_rp(K, sy) > 0.0


def test_invalid_yield_raises() -> None:
    with pytest.raises(ValueError):
        substrate_rp(50.0, 0.0)
    with pytest.raises(ValueError):
        irwin_plane_stress_rp(50.0, -1.0)
    with pytest.raises(ValueError):
        irwin_plane_strain_rp(50.0, 0.0)
    with pytest.raises(ValueError):
        dugdale_rp(50.0, 0.0)


# --------------------------------------------------------------------------- #
# Claim (b): substrate prediction is mathematically derivable from sigma<=1/2 #
# --------------------------------------------------------------------------- #


def test_substrate_rp_matches_derivation() -> None:
    """Substrate r_p must equal (1/(2*pi)) * (K_I/sigma_y)**2.

    The substrate cap sigma <= 1/2 closes the ODE for the cohesive-zone
    boundary at exactly the same prefactor as the Irwin plane-stress
    yielding cutoff. Verify identity to 1e-12 across all materials.
    """

    for name, props in MATERIALS.items():
        K = props["K_IC"]
        sy = props["sigma_y"]
        expected = (K * K) / (2.0 * math.pi * sy * sy)
        assert math.isclose(
            substrate_rp(K, sy), expected, rel_tol=1.0e-12
        ), f"substrate derivation broken for {name}"


def test_substrate_rp_equals_irwin_plane_stress() -> None:
    """The cap-derivation gives the same prefactor 1/(2*pi)."""

    for name, props in MATERIALS.items():
        K = props["K_IC"]
        sy = props["sigma_y"]
        assert math.isclose(
            substrate_rp(K, sy),
            irwin_plane_stress_rp(K, sy),
            rel_tol=1.0e-12,
        ), f"substrate vs Irwin PS broken for {name}"


def test_irwin_plane_strain_factor_three_smaller() -> None:
    """Irwin plane strain is 1/3 of substrate (1/(6*pi) vs 1/(2*pi))."""

    for name, props in MATERIALS.items():
        K = props["K_IC"]
        sy = props["sigma_y"]
        ratio = substrate_rp(K, sy) / irwin_plane_strain_rp(K, sy)
        assert math.isclose(ratio, 3.0, rel_tol=1.0e-12), (
            f"plane-strain ratio broken for {name}: ratio={ratio}"
        )


def test_dugdale_factor() -> None:
    """Dugdale is pi**2/4 times the substrate (i.e. larger by ~2.467)."""

    for name, props in MATERIALS.items():
        K = props["K_IC"]
        sy = props["sigma_y"]
        ratio = dugdale_rp(K, sy) / substrate_rp(K, sy)
        expected_ratio = (math.pi ** 2) / 4.0
        assert math.isclose(ratio, expected_ratio, rel_tol=1.0e-12)


# --------------------------------------------------------------------------- #
# Claim (c): r_p scales as (K_I/sigma_y)**2                                   #
# --------------------------------------------------------------------------- #


def test_substrate_rp_quadratic_in_K() -> None:
    """Doubling K_I must quadruple r_p (sigma_y fixed)."""

    sy = 500.0
    base = substrate_rp(40.0, sy)
    doubled = substrate_rp(80.0, sy)
    assert math.isclose(doubled / base, 4.0, rel_tol=1.0e-12)


def test_substrate_rp_inverse_quadratic_in_sigma_y() -> None:
    """Doubling sigma_y must quarter r_p (K_I fixed)."""

    K = 50.0
    base = substrate_rp(K, 400.0)
    doubled = substrate_rp(K, 800.0)
    assert math.isclose(doubled / base, 0.25, rel_tol=1.0e-12)


@pytest.mark.parametrize(
    "rp_func",
    [substrate_rp, irwin_plane_stress_rp, irwin_plane_strain_rp, dugdale_rp],
)
def test_all_formulas_quadratic_scaling(rp_func) -> None:
    """All four formulas must share (K_I/sigma_y)**2 scaling."""

    base = rp_func(50.0, 500.0)
    K_doubled = rp_func(100.0, 500.0)
    sy_doubled = rp_func(50.0, 1000.0)
    assert math.isclose(K_doubled / base, 4.0, rel_tol=1.0e-12)
    assert math.isclose(sy_doubled / base, 0.25, rel_tol=1.0e-12)


# --------------------------------------------------------------------------- #
# Claim (d): order-of-magnitude check on a canonical material                 #
# --------------------------------------------------------------------------- #


def test_steel_4340_within_order_of_magnitude_of_irwin_plane_stress() -> None:
    """For Steel 4340, substrate r_p must be within 10x of Irwin plane stress.

    By derivation they are equal, so this is a much stronger test than the
    minimum requirement.
    """

    K = MATERIALS["Steel 4340"]["K_IC"]
    sy = MATERIALS["Steel 4340"]["sigma_y"]
    sub = substrate_rp(K, sy)
    ips = irwin_plane_stress_rp(K, sy)
    ratio = sub / ips
    assert 0.1 < ratio < 10.0, f"out-of-order-of-magnitude: {ratio}"


def test_steel_4340_value_in_expected_range() -> None:
    """Sanity: Steel 4340 plastic zone should be sub-millimetre."""

    K = MATERIALS["Steel 4340"]["K_IC"]
    sy = MATERIALS["Steel 4340"]["sigma_y"]
    rp = substrate_rp(K, sy)
    # K=50 MPa.m^0.5, sigma_y=860 MPa => r_p ~ 0.54 mm.
    assert 0.0001 < rp < 0.005


# --------------------------------------------------------------------------- #
# Claim (e): compare_predictions / run_test schemas and conclusions           #
# --------------------------------------------------------------------------- #


def test_compare_predictions_schema() -> None:
    res = compare_predictions()
    assert set(res.keys()) == set(MATERIALS.keys())
    expected_keys = {
        "K_IC_MPa_sqrt_m",
        "sigma_y_MPa",
        "E_GPa",
        "substrate_m",
        "irwin_ps_m",
        "irwin_pe_m",
        "dugdale_m",
    }
    for row in res.values():
        assert set(row.keys()) == expected_keys
        for v in row.values():
            assert v > 0.0


def test_run_test_identifies_irwin_plane_stress_as_closest() -> None:
    res = run_test()
    assert res["closest_match"] == "irwin_plane_stress"
    assert math.isclose(
        res["ratios"]["substrate_over_irwin_ps"], 1.0, rel_tol=1.0e-12
    )
    assert math.isclose(
        res["ratios"]["substrate_over_irwin_pe"], 3.0, rel_tol=1.0e-12
    )
    expected_dugdale = 1.0 / ((math.pi ** 2) / 4.0)
    assert math.isclose(
        res["ratios"]["substrate_over_dugdale"],
        expected_dugdale,
        rel_tol=1.0e-12,
    )


def test_run_test_correlations_all_unity() -> None:
    """All four formulas are linear in (K_I/sigma_y)**2 so r=+1 across set."""

    res = run_test()
    for k, v in res["correlations"].items():
        assert math.isclose(v, 1.0, rel_tol=1.0e-9), f"{k}: r={v}"


def test_run_test_no_significant_diff_vs_irwin_ps() -> None:
    """Substrate matches Irwin plane stress to machine precision."""

    res = run_test()
    assert res["significant_diffs"] == []
