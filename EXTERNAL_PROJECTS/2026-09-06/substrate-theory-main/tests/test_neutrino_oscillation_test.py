"""Tests for stiff_medium.neutrino_oscillation_test (substrate vs NuFIT 5.2)."""

import math
import os

import pytest
from scipy import constants

from stiff_medium.neutrino_oscillation_test import (
    NUFIT_5_2,
    NeutrinoOscillationTest,
    Comparison,
    render_visual,
)


@pytest.fixture(scope="module")
def test_obj():
    return NeutrinoOscillationTest()


# ---------------------------------------------------------------------------
# Substrate-prediction algebra
# ---------------------------------------------------------------------------

def test_alpha_used_is_codata(test_obj):
    """The α used must equal scipy.constants.alpha (CODATA 2018)."""
    assert test_obj.alpha == constants.alpha


def test_sin2_theta12_equals_42_alpha(test_obj):
    """Substrate formula: sin²θ_12 = 42 α."""
    assert math.isclose(test_obj.sin2_theta12,
                        42.0 * constants.alpha,
                        rel_tol=1e-12)


def test_sin2_theta13_equals_3_alpha(test_obj):
    """Substrate formula: sin²θ_13 = 3 α."""
    assert math.isclose(test_obj.sin2_theta13,
                        3.0 * constants.alpha,
                        rel_tol=1e-12)


def test_sin2_theta23_equals_half_plus_2pi_alpha(test_obj):
    """Substrate formula: sin²θ_23 = ½ + 2π α."""
    expected = 0.5 + 2.0 * math.pi * constants.alpha
    assert math.isclose(test_obj.sin2_theta23, expected, rel_tol=1e-12)


def test_delta_CP_substrate_ansatz_is_3pi_over_4(test_obj):
    """Substrate suggests δ_CP = 3π/4."""
    assert math.isclose(test_obj.delta_CP_rad, 3.0 * math.pi / 4.0,
                        rel_tol=1e-12)
    assert math.isclose(test_obj.delta_CP_deg, 135.0, rel_tol=1e-12)


def test_substrate_angles_in_physical_range(test_obj):
    """All three angles should land in [0, π/2]."""
    for theta in (test_obj.theta12_rad,
                  test_obj.theta13_rad,
                  test_obj.theta23_rad):
        assert 0.0 < theta < math.pi / 2.0


# ---------------------------------------------------------------------------
# Comparison records
# ---------------------------------------------------------------------------

def test_comparisons_returns_five_rows(test_obj):
    cs = test_obj.comparisons()
    assert len(cs) == 5
    names = [c.name for c in cs]
    assert "sin²θ_12" in names
    assert "sin²θ_13" in names
    assert "sin²θ_23" in names
    assert "Δm²_21" in names
    assert "|Δm²_31|" in names


def test_comparison_dataclass_residual_and_sigma(test_obj):
    cs = test_obj.comparisons()
    for c in cs:
        assert isinstance(c, Comparison)
        # Residual
        assert math.isclose(c.residual, c.pred - c.obs, rel_tol=1e-12)
        # σ-distance
        if c.sigma_obs > 0:
            expected = (c.pred - c.obs) / c.sigma_obs
            assert math.isclose(c.sigma_distance, expected, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Quality of substrate predictions vs NuFIT 5.2
# ---------------------------------------------------------------------------

def test_sin2_theta12_within_1sigma_of_NuFIT_52(test_obj):
    """Banner test: sin²θ_12 = 42α should match NuFIT 5.2 within 1σ."""
    cs = test_obj.comparisons()
    c12 = next(c for c in cs if c.name == "sin²θ_12")
    obs, sig = NUFIT_5_2["sin2_theta12"]
    assert (obs - sig) <= c12.pred <= (obs + sig), (
        f"sin²θ_12 substrate = {c12.pred:.5f} not in NuFIT 1σ "
        f"window [{obs - sig:.5f}, {obs + sig:.5f}]")


def test_sin2_theta13_within_1sigma_of_NuFIT_52(test_obj):
    """Banner test: sin²θ_13 = 3α should match NuFIT 5.2 within 1σ."""
    cs = test_obj.comparisons()
    c13 = next(c for c in cs if c.name == "sin²θ_13")
    obs, sig = NUFIT_5_2["sin2_theta13"]
    assert (obs - sig) <= c13.pred <= (obs + sig), (
        f"sin²θ_13 substrate = {c13.pred:.5f} not in NuFIT 1σ "
        f"window [{obs - sig:.5f}, {obs + sig:.5f}]")


def test_sin2_theta23_within_2sigma_of_NuFIT_52(test_obj):
    """sin²θ_23 = ½ + 2πα is the weakest of the three; within 2σ."""
    cs = test_obj.comparisons()
    c23 = next(c for c in cs if c.name == "sin²θ_23")
    assert abs(c23.sigma_distance) < 2.0, (
        f"sin²θ_23 substrate {c23.pred:.5f} more than 2σ from "
        f"NuFIT 5.2 {c23.obs:.5f} ± {c23.sigma_obs:.5f}")


def test_delta_CP_within_3sigma_of_NuFIT_52(test_obj):
    """δ_CP substrate ansatz 135° vs NuFIT 5.2 197° ± 27°: ~2.3σ off."""
    d = test_obj.delta_CP_comparison()
    assert abs(d.sigma_distance) < 3.0, (
        f"δ_CP substrate {d.pred:.1f}° more than 3σ from "
        f"NuFIT 5.2 {d.obs:.1f}° ± {d.sigma_obs:.1f}°")


def test_worst_sigma_distance_under_two(test_obj):
    """Aggregate honesty score: worst (excl. δ_CP) under 2σ."""
    rep = test_obj.report(file=open(os.devnull, "w"))
    assert rep["worst_sigma"] < 2.0


# ---------------------------------------------------------------------------
# Visual
# ---------------------------------------------------------------------------

def test_render_visual_writes_file(tmp_path):
    out = tmp_path / "neutrino_oscillation_test.png"
    path = render_visual(str(out))
    assert os.path.exists(path)
    assert os.path.getsize(path) > 5000  # non-trivial image


def test_main_visual_present_after_run():
    """The canonical visuals/135_neutrino_oscillation.png should exist
    once the module's main() has been invoked at least once."""
    expected = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..",
        "visuals", "135_neutrino_oscillation.png"))
    if not os.path.exists(expected):
        # Generate it now so tests are self-contained
        render_visual(expected)
    assert os.path.exists(expected)
    assert os.path.getsize(expected) > 5000
