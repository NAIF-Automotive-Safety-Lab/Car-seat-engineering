"""Tests for ``src.stiff_medium.substrate_polywell``.

Covers:
    * Q_3 cube cell forces 6 face / coil positions; |O_h| = 48
    * Wiffle Ball radius derives from k_rank=5 substrate buffer
    * σ ≤ 1/2 cap drives the density limit; scales as B^2 / V_well
    * Substrate exponents (a, b) = (4, 3) match Bussard's empirical fit
    * predict_wb6_neutron_rate returns a value within an order of
      magnitude of the reported 1e9 neutrons/s for default parameters
    * extrapolate_to_breakeven gives sensible engineering scales
      (r ~ 1 m at 5 T, B ~ a few T at r = 1.5 m)
    * fusion cross section reduces to standard astrophysical S-factor
      form when substrate correction is set to 1
    * optimal injection energies for DD/DT match standard Bosch-Hale
      peaks; pB11 substrate prediction is at LOWER energy than the
      bare cross-section peak
    * 4-coil and 8-coil arrangements are NOT predicted by Q_3 (the
      6-coil layout is forced)
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium import b3_constants as B3
from src.stiff_medium.substrate_polywell import (
    MU_0,
    SIGMA_CAP,
    SUBSTRATE_EXPONENT_B,
    SUBSTRATE_EXPONENT_R,
    PolywellGeometry,
    PolywellSimulator,
)


# ---------------------------------------------------------------------------
# Q_3 geometry: 6 face positions
# ---------------------------------------------------------------------------

def test_q3_has_six_face_coils():
    """The cube cell Q_3 has exactly 6 faces → 6 coil positions."""
    g = PolywellGeometry(r_coil=0.15)
    assert g.n_faces == 6
    coils = g.coil_positions()
    assert coils.shape == (6, 3)


def test_q3_has_twelve_edges():
    """Q_3 has 12 edges (used in face-pair coupling counting)."""
    g = PolywellGeometry(r_coil=0.15)
    assert g.n_edges == 12


def test_q3_has_eight_vertices():
    """Q_3 has 8 vertices (used in O_h symmetry orbit count)."""
    g = PolywellGeometry(r_coil=0.15)
    assert g.n_vertices == 8


def test_q3_automorphism_group_O_h():
    """Cube symmetry group |O_h| = 48 (full octahedral group with reflections)."""
    g = PolywellGeometry(r_coil=0.15)
    assert g.automorphism_order == 48


def test_q3_face_positions_on_axes():
    """Face centres sit at ±r along each Cartesian axis (substrate alignment)."""
    g = PolywellGeometry(r_coil=0.15)
    coils = g.coil_positions()
    # Every coil must lie on exactly one axis (two coords = 0)
    for c in coils:
        zeros = np.sum(np.isclose(c, 0.0))
        assert zeros == 2, f"coil {c} not on a Cartesian axis"
    # Norms all equal r_coil
    norms = np.linalg.norm(coils, axis=1)
    assert np.allclose(norms, 0.15, rtol=1e-12)


def test_q3_coil_normals_unit():
    """Coil normals (face normals) are unit vectors."""
    g = PolywellGeometry(r_coil=0.15)
    n = g.coil_normals()
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0, rtol=1e-12)


def test_polywell_geometry_substrate_forced_flag():
    """Geometry advertises that 6-coil layout is substrate-forced."""
    g = PolywellGeometry(r_coil=0.15)
    assert g.forced_by_substrate is True


def test_4coil_and_8coil_NOT_q3():
    """Substrate Q_3 does NOT admit a 4-coil or 8-coil arrangement.

    Tests that ``PolywellGeometry`` always uses 6 (and the other
    Q_3 invariants); a 4-coil device would correspond to a tetrahedral
    K_4 cell, an 8-coil device to a vertex-attached layout — both
    different cell topologies, not Q_3.
    """
    g = PolywellGeometry(r_coil=0.15)
    assert g.n_faces != 4
    assert g.n_faces != 8


# ---------------------------------------------------------------------------
# Wiffle Ball geometry from substrate σ ≤ 1/2 cap
# ---------------------------------------------------------------------------

def test_wiffle_ball_radius_from_k_rank():
    """r_wb / r_coil = (1 - 1/k_rank) = 4/5 = 0.8 at canonical K_rank=5."""
    g = PolywellGeometry(r_coil=0.15)
    assert math.isclose(
        g.wiffle_ball_radius() / g.r_coil, 0.8, rel_tol=1e-12
    )
    # Confirm the buffer comes from K_rank
    assert math.isclose(
        g.wiffle_ball_radius() / g.r_coil,
        1.0 - 1.0 / B3.K_rank,
        rel_tol=1e-12,
    )


def test_confinement_volume_sphere():
    """Volume = (4/3) π r_wb^3 with r_wb = 0.8 r_coil."""
    g = PolywellGeometry(r_coil=0.15)
    expected = (4.0 / 3.0) * math.pi * (0.8 * 0.15) ** 3
    assert math.isclose(g.confinement_volume(), expected, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Wiffle Ball density limit from σ ≤ 1/2
# ---------------------------------------------------------------------------

def test_wiffle_ball_density_limit_formula():
    """n_max = SIGMA_CAP * B^2 / (2 mu_0 * e * V)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    B = 0.10
    n_max = sim.wiffle_ball_density_limit(B)
    # Closed form
    p_mag = B * B / (2.0 * MU_0)
    E_e = 1.602_176_634e-19 * 12_000.0
    expected = SIGMA_CAP * p_mag / E_e
    assert math.isclose(n_max, expected, rel_tol=1e-12)


def test_wiffle_ball_density_scales_as_Bsq():
    """Doubling B quadruples n_max (B^2 scaling)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    n1 = sim.wiffle_ball_density_limit(0.10)
    n2 = sim.wiffle_ball_density_limit(0.20)
    assert math.isclose(n2 / n1, 4.0, rel_tol=1e-12)


def test_wiffle_ball_density_inverse_voltage():
    """Doubling V_well halves n_max (1/V scaling)."""
    g = PolywellGeometry(r_coil=0.15)
    sim_a = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    sim_b = PolywellSimulator(geometry=g, voltage_V=24_000.0)
    n_a = sim_a.wiffle_ball_density_limit(0.10)
    n_b = sim_b.wiffle_ball_density_limit(0.10)
    assert math.isclose(n_b / n_a, 0.5, rel_tol=1e-12)


def test_density_limit_rejects_zero_field():
    """B = 0 → no confinement → ValueError."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    with pytest.raises(ValueError):
        sim.wiffle_ball_density_limit(0.0)


# ---------------------------------------------------------------------------
# Substrate exponents reduce to Bussard's B^4 r^3
# ---------------------------------------------------------------------------

def test_substrate_exponents_are_4_and_3():
    """The substrate-forced power scaling exponents."""
    a, b = PolywellSimulator.substrate_exponents()
    assert (a, b) == (4, 3)
    assert SUBSTRATE_EXPONENT_B == 4
    assert SUBSTRATE_EXPONENT_R == 3


def test_substrate_matches_bussard_empirical_exponents():
    """Substrate prediction and Bussard's empirical fit agree on (a, b)."""
    assert (
        PolywellSimulator.substrate_exponents()
        == PolywellSimulator.bussard_empirical_exponents()
    )


def test_power_scaling_doubling_B():
    """Doubling B multiplies P by 16 (B^4)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    p1 = sim.power_scaling(B=0.10, r=0.15)
    p2 = sim.power_scaling(B=0.20, r=0.15)
    assert math.isclose(p2 / p1, 16.0, rel_tol=1e-12)


def test_power_scaling_doubling_r():
    """Doubling r multiplies P by 8 (r^3)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    p1 = sim.power_scaling(B=0.10, r=0.15)
    p2 = sim.power_scaling(B=0.10, r=0.30)
    assert math.isclose(p2 / p1, 8.0, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# WB-6 prediction within an order of magnitude
# ---------------------------------------------------------------------------

def test_wb6_neutron_rate_within_order_of_magnitude():
    """Predicted WB-6 neutron rate is within 10× of reported 1e9 /s."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    res = sim.predict_wb6_neutron_rate()
    rate_subs = res["rate_per_s_substrate"]
    assert 1.0e8 < rate_subs < 1.0e10, (
        f"WB-6 prediction {rate_subs:.3e} not within 10x of 1e9"
    )


def test_wb6_neutron_rate_substrate_correction_small():
    """K_4 face-pair correction is < 5 % at thermal D-D energies."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    res = sim.predict_wb6_neutron_rate()
    f_K4 = res["ratio_substrate_to_bare"]
    assert 1.0 <= f_K4 < 1.05, (
        f"K_4 correction {f_K4} should be small (< 5%)"
    )


def test_compare_to_empirical_bundle():
    """compare_to_empirical returns the expected schema and matches scaling."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    cmp = sim.compare_to_empirical()
    # Geometry: 6 coils both ways
    assert cmp["geometry_match"] is True
    assert cmp["geometry_n_coils_substrate"] == 6
    # Scaling: substrate matches empirical
    assert cmp["scaling_match"] is True
    # Within order of magnitude flag
    assert cmp["wb6_within_one_order_of_magnitude"] is True


# ---------------------------------------------------------------------------
# Extrapolation to breakeven gives sensible engineering scales
# ---------------------------------------------------------------------------

def test_breakeven_radius_at_5T_is_meter_scale():
    """At 5 T, r for 1 MW D-D fusion is roughly meter-scale (0.1–10 m)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    ext = sim.extrapolate_to_breakeven()
    r = ext["r_for_breakeven_m_at_5T"]
    assert 0.1 < r < 10.0, f"r = {r:.3f} m not meter-scale"


def test_breakeven_field_at_1p5m_is_few_tesla():
    """At r = 1.5 m, B for 1 MW D-D fusion is a few T (1–20 T)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    ext = sim.extrapolate_to_breakeven()
    B = ext["B_for_breakeven_T_at_1p5m"]
    assert 1.0 < B < 20.0, f"B = {B:.3f} T not few-Tesla scale"


def test_breakeven_min_confinement_field_O_10T():
    """Confinement at n=1e21 m^-3, V=100 kV needs O(10 T) field."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    ext = sim.extrapolate_to_breakeven()
    B_min = ext["B_min_for_confinement_T"]
    assert 3.0 < B_min < 30.0, f"B_min = {B_min:.3f} T not O(10 T)"


def test_breakeven_verdicts_present():
    """Verdict strings are populated."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g, voltage_V=12_000.0)
    ext = sim.extrapolate_to_breakeven()
    assert ext["verdict_5T_design"] in ("feasible", "marginal", "infeasible")
    assert ext["verdict_1p5m_design"] in ("feasible", "marginal", "infeasible")


# ---------------------------------------------------------------------------
# Cross-section sanity
# ---------------------------------------------------------------------------

def test_xsec_DD_falls_below_threshold():
    """At 1 keV the D-D cross section is exponentially small (Gamow tunnel)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    sigma_low = sim.fusion_cross_section_substrate(E_keV=1.0, fuel="DD")
    sigma_mid = sim.fusion_cross_section_substrate(E_keV=20.0, fuel="DD")
    assert sigma_low < sigma_mid
    # Rough magnitude: σ at ~10 keV is ~1e-3 b for D-D
    assert 1e-7 < sigma_mid < 1.0


def test_xsec_DT_larger_than_DD():
    """D-T cross section is much larger than D-D at fusion-relevant T."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    sigma_DT = sim.fusion_cross_section_substrate(E_keV=15.0, fuel="DT")
    sigma_DD = sim.fusion_cross_section_substrate(E_keV=15.0, fuel="DD")
    assert sigma_DT > 10.0 * sigma_DD


def test_xsec_pB11_smaller_than_DT():
    """p-B11 cross section is much smaller than D-T at fixed E."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    sigma_pB = sim.fusion_cross_section_substrate(E_keV=200.0, fuel="pB11")
    sigma_DT = sim.fusion_cross_section_substrate(E_keV=200.0, fuel="DT")
    assert sigma_pB < sigma_DT


def test_xsec_zero_energy_returns_zero():
    """E_keV = 0 returns σ = 0 (no kinetic energy)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    assert sim.fusion_cross_section_substrate(0.0, "DD") == 0.0


# ---------------------------------------------------------------------------
# Optimal injection energies
# ---------------------------------------------------------------------------

def test_optimal_DD_around_15keV():
    """D-D resonance ≈ 15 keV (matches Bosch-Hale optimum)."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    E = sim.optimal_injection_energy("DD")
    assert 10_000 < E < 20_000


def test_optimal_DT_around_13keV():
    """D-T resonance ≈ 13.5 keV."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    E = sim.optimal_injection_energy("DT")
    assert 10_000 < E < 20_000


def test_optimal_pB11_substrate_below_bare_peak():
    """Substrate predicts pB11 optimum BELOW the standard ~580 keV peak."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    E = sim.optimal_injection_energy("pB11")
    # Substrate prediction: ~50 keV (cone-bouncing K_4 face-pair resonance)
    assert E < 100_000  # well below the bare peak at 580 keV
    assert E > 10_000   # above the Coulomb-barrier suppression region


def test_unknown_fuel_raises():
    """Unknown fuel raises ValueError."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    with pytest.raises(ValueError):
        sim.optimal_injection_energy("DHe3")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        sim.fusion_cross_section_substrate(10.0, "DHe3")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Bundled all_predictions schema
# ---------------------------------------------------------------------------

def test_all_predictions_schema():
    """all_predictions returns the expected top-level keys."""
    g = PolywellGeometry(r_coil=0.15)
    sim = PolywellSimulator(geometry=g)
    out = sim.all_predictions()
    expected = {
        "geometry",
        "scaling_exponents",
        "WB6_compare",
        "breakeven_extrapolation_DD_1MW",
        "optimal_energies_eV",
        "derivability_tags",
    }
    assert expected.issubset(out.keys())
    # Derivability tags refer to A/B/C categories
    tags = out["derivability_tags"]
    assert any("A" in v for v in tags.values())
    assert any("B" in v for v in tags.values())
    assert any("C" in v for v in tags.values())
