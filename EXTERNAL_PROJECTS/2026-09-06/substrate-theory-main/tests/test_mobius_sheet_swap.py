"""Tests for src/stiff_medium/mobius_sheet_swap.py.

Covers the five required structural invariants:

(a) the Möbius bundle has exactly **2 sheets**;
(b) the sheet-swap involution τ satisfies **τ² = id**;
(c) **charged** input states map to **distinct antiparticles**;
(d) **neutral** input states map to **themselves** (Majorana ν = ν̄);
(e) the saturation cap **σ = 1/2** is a **fixed point** of τ.

Plus the half-integer spin from the double-cover holonomy
(``(−1)^n``) and the cap-invariance of simulator trajectories.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.mobius_sheet_swap import (
    SHEET_A,
    SHEET_B,
    SIGMA_MAX,
    Excitation,
    MobiusSheetGeometry,
    SheetSwapSimulator,
    default_charged,
    default_neutral,
    summary,
)


# ---------------------------------------------------------------------------
# (a) two sheets
# ---------------------------------------------------------------------------

def test_geometry_has_exactly_two_sheets():
    geom = MobiusSheetGeometry()
    assert geom.n_sheets() == 2
    assert geom.sheets() == (+1, -1)
    assert SHEET_A == +1
    assert SHEET_B == -1


# ---------------------------------------------------------------------------
# (b) involution: τ² = id
# ---------------------------------------------------------------------------

def test_involution_swap_flips_sheet():
    """τ(θ, +1) = (θ, −1) and τ(θ, −1) = (θ, +1)."""
    theta_out, s_out = MobiusSheetGeometry.involution(1.234, +1)
    assert math.isclose(theta_out, 1.234)
    assert s_out == -1
    theta_out, s_out = MobiusSheetGeometry.involution(2.5, -1)
    assert math.isclose(theta_out, 2.5)
    assert s_out == +1


def test_involution_squared_is_identity_explicit():
    """Spot-check τ² = id at a handful of (θ, s)."""
    grid_theta = [0.0, 0.5, 1.0, math.pi, 4.0, 6.2]
    for s in (+1, -1):
        for theta in grid_theta:
            t2, s2 = MobiusSheetGeometry.involution_squared(theta, s)
            assert math.isclose(t2, theta)
            assert s2 == s


def test_involution_squared_is_identity_dense_grid():
    """Dense regression: τ² = id on a 64-point base grid."""
    assert MobiusSheetGeometry.is_involution(n_samples=64) is True


def test_involution_rejects_non_sheet_label():
    with pytest.raises(ValueError):
        MobiusSheetGeometry.involution(0.0, 0)
    with pytest.raises(ValueError):
        MobiusSheetGeometry.involution(0.0, 2)


# ---------------------------------------------------------------------------
# (c) charged states swap to distinct antiparticle
# ---------------------------------------------------------------------------

def test_charged_swap_produces_distinct_antiparticle():
    sim = SheetSwapSimulator()
    e_minus = Excitation(charge=-1.0, sheet=SHEET_A, amplitude=0.3,
                         label="electron")
    e_plus = sim.apply_sheet_swap(e_minus)
    # Distinct: charge AND sheet flipped.
    assert e_plus.charge == +1.0
    assert e_plus.sheet == SHEET_B
    # Distinctness is the structural test: not equal to input.
    assert e_plus != e_minus
    # Conventional name flip.
    assert e_plus.label == "anti-electron"


def test_charged_swap_is_involutive_on_state():
    """τ²(electron) == electron."""
    sim = SheetSwapSimulator()
    e_minus = Excitation(charge=-1.0, sheet=SHEET_A, amplitude=0.3,
                         label="electron")
    twice = sim.apply_sheet_swap(sim.apply_sheet_swap(e_minus))
    assert twice.charge == e_minus.charge
    assert twice.sheet == e_minus.sheet
    assert math.isclose(twice.amplitude, e_minus.amplitude)


def test_evolve_one_loop_charged_endpoint_is_antiparticle():
    sim = SheetSwapSimulator(n_steps=64)
    e_minus = Excitation(charge=-1.0, sheet=SHEET_A, amplitude=0.4,
                         label="electron")
    out = sim.evolve_one_loop(e_minus)
    assert out["is_majorana"] is False
    end = out["end"]
    assert end.charge == +1.0
    assert end.sheet == SHEET_B


# ---------------------------------------------------------------------------
# (d) neutral states map to themselves (Majorana)
# ---------------------------------------------------------------------------

def test_neutral_swap_maps_to_self_majorana():
    sim = SheetSwapSimulator()
    nu = Excitation(charge=0.0, sheet=SHEET_A, amplitude=0.3,
                    label="neutrino")
    nu_bar = sim.apply_sheet_swap(nu)
    # Identical: charge unchanged, sheet unchanged, amplitude unchanged.
    assert nu_bar.charge == nu.charge
    assert nu_bar.sheet == nu.sheet
    assert math.isclose(nu_bar.amplitude, nu.amplitude)
    assert nu_bar.label == nu.label
    # Hash/equality identity (frozen dataclass): same state.
    assert nu_bar == nu


def test_neutral_loop_flagged_as_majorana():
    sim = SheetSwapSimulator(n_steps=64)
    nu = default_neutral()
    out = sim.evolve_one_loop(nu)
    assert out["is_majorana"] is True
    assert out["end"] == out["start"]


def test_neutral_swap_works_on_both_sheets():
    sim = SheetSwapSimulator()
    nu_a = Excitation(charge=0.0, sheet=SHEET_A, amplitude=0.4)
    nu_b = Excitation(charge=0.0, sheet=SHEET_B, amplitude=0.4)
    assert sim.apply_sheet_swap(nu_a) == nu_a
    assert sim.apply_sheet_swap(nu_b) == nu_b


# ---------------------------------------------------------------------------
# (e) σ = 1/2 fixed point and cap enforcement
# ---------------------------------------------------------------------------

def test_saturation_cap_value():
    assert SIGMA_MAX == 0.5


def test_saturation_cap_is_fixed_point_flag():
    assert MobiusSheetGeometry.saturation_cap_is_fixed_point() is True


def test_saturation_cap_rejects_above_threshold_state():
    """Trying to construct an excitation above σ = 1/2 must raise."""
    with pytest.raises(ValueError):
        Excitation(charge=-1.0, sheet=SHEET_A, amplitude=0.6)


def test_saturation_cap_holds_along_loop():
    sim = SheetSwapSimulator(n_steps=128)
    out = sim.evolve_one_loop(default_charged())
    assert SheetSwapSimulator.cap_is_invariant(out["amplitudes"]) is True
    assert float(np.max(out["amplitudes"])) <= SIGMA_MAX + 1e-12


def test_saturation_cap_caps_an_input_at_threshold():
    """An excitation built at exactly σ = 1/2 should still loop cleanly."""
    sim = SheetSwapSimulator(n_steps=64)
    edge = Excitation(charge=-1.0, sheet=SHEET_A, amplitude=SIGMA_MAX)
    out = sim.evolve_one_loop(edge)
    assert float(np.max(out["amplitudes"])) <= SIGMA_MAX + 1e-12


# ---------------------------------------------------------------------------
# (f) double cover → spin-½
# ---------------------------------------------------------------------------

def test_loop_holonomy_alternates_with_loop_count():
    geom = MobiusSheetGeometry()
    assert geom.loop_holonomy(0) == +1.0
    assert geom.loop_holonomy(1) == -1.0
    assert geom.loop_holonomy(2) == +1.0
    assert geom.loop_holonomy(3) == -1.0


def test_double_cover_predicate():
    assert MobiusSheetGeometry.is_double_cover() is True


def test_spin_is_one_half():
    assert MobiusSheetGeometry.spin_value() == 0.5


# ---------------------------------------------------------------------------
# Geometry surface for the visualiser
# ---------------------------------------------------------------------------

def test_surface_xyz_shapes_are_consistent():
    geom = MobiusSheetGeometry(n_theta=80, n_v=12)
    X, Y, Z = geom.surface_xyz()
    assert X.shape == (12, 80)
    assert Y.shape == (12, 80)
    assert Z.shape == (12, 80)


def test_sheet_color_field_signs():
    geom = MobiusSheetGeometry(n_theta=10, n_v=4)
    s = geom.sheet_color_field()
    assert s.shape == (4, 10)
    # Top half-rows must be positive (sheet A), bottom negative (sheet B).
    assert np.all(s[-1] > 0)
    assert np.all(s[0] < 0)


# ---------------------------------------------------------------------------
# End-to-end summary
# ---------------------------------------------------------------------------

def test_summary_reports_all_invariants():
    s = summary()
    assert s["n_sheets"] == 2
    assert s["involution_squared_is_identity"] is True
    assert s["charged_swap_distinct"] is True
    assert s["neutral_swap_identical"] is True
    assert s["saturation_cap_value"] == 0.5
    assert s["saturation_cap_fixed_point"] is True
    assert s["double_cover_holonomy_one_loop"] == -1.0
    assert s["double_cover_holonomy_two_loops"] == +1.0
    assert s["spin_value"] == 0.5


def test_default_factories():
    c = default_charged()
    n = default_neutral()
    assert c.charge == -1.0
    assert n.charge == 0.0
    assert n.is_neutral
    assert not c.is_neutral
