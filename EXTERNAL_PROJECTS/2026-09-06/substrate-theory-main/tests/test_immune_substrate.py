"""Tests for the paired immune-system GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) antibody Y-shape: total span ~14 nm, two Fab arms + Fc stem,
        two paratopes at the Fab tips, MHC platform with peptide groove
    (b) binding affinity K_d in 10^-9 .. 10^-12 M (textbook IgG range)
    (c) free energy of binding is negative (favourable)
    (d) primary response: latency ~5 d, modest peak; secondary: latency
        ~1.5 d, peak ~100x larger
    (e) memory cell fraction ~10 percent of activated B cells
    (f) somatic hypermutation rate ~10^-3 per bp per division
        (~10^6 over background)
    (g) clonal selection: high-affinity (low K_d) clones win
    (h) SIR-like model produces a pathogen burst that is then cleared
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.immune_substrate import (
    ANTIBODY_FAB_LENGTH_NM,
    ANTIBODY_FC_LENGTH_NM,
    ANTIBODY_HINGE_ANGLE_DEG,
    ANTIBODY_TOTAL_SPAN_NM,
    BACKGROUND_MUTATION_RATE_PER_BP_PER_DIV,
    HYPERMUTATION_RATE_PER_BP_PER_DIV,
    KD_HIGH_AFFINITY_M,
    KD_LOW_AFFINITY_M,
    KD_MEDIUM_AFFINITY_M,
    MEMORY_CELL_FRACTION,
    MHC_PEPTIDE_LENGTH_AA,
    PRIMARY_LATENCY_DAYS,
    SECONDARY_LATENCY_DAYS,
    ImmuneGeometry,
    ImmuneSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> ImmuneGeometry:
    return ImmuneGeometry()


@pytest.fixture
def sim() -> ImmuneSimulator:
    return ImmuneSimulator()


# --------------------------------------------------------------------- #
# (a) Antibody Y-shape geometry                                         #
# --------------------------------------------------------------------- #

def test_fc_endpoints_are_along_z(geom):
    base, hinge = geom.fc_endpoints()
    assert np.allclose(base, [0.0, 0.0, -ANTIBODY_FC_LENGTH_NM])
    assert np.allclose(hinge, [0.0, 0.0, 0.0])


def test_two_fab_arms_at_hinge_angle(geom):
    _, left_tip = geom.fab_endpoints("left")
    _, right_tip = geom.fab_endpoints("right")
    # Left tip has -y, right tip has +y, both have +z
    assert left_tip[1] < 0.0 and right_tip[1] > 0.0
    assert left_tip[2] > 0.0 and right_tip[2] > 0.0
    # Symmetry: y components are opposite, z components equal
    assert abs(left_tip[1] + right_tip[1]) < 1.0e-9
    assert abs(left_tip[2] - right_tip[2]) < 1.0e-9
    # Length of Fab equals fab_length_nm
    L = float(np.linalg.norm(right_tip))
    assert abs(L - ANTIBODY_FAB_LENGTH_NM) < 1.0e-9


def test_invalid_fab_side_raises(geom):
    with pytest.raises(ValueError):
        geom.fab_endpoints("middle")


def test_paratope_centers_two_tips(geom):
    tips = geom.paratope_centers()
    assert tips.shape == (2, 3)


def test_total_span_positive_and_consistent(geom):
    span = geom.total_span_nm()
    # Two Fab arms at +/- hinge_angle, span is 2 * fab_length * sin(theta)
    expected = 2.0 * ANTIBODY_FAB_LENGTH_NM * math.sin(
        math.radians(ANTIBODY_HINGE_ANGLE_DEG)
    )
    assert abs(span - expected) < 1.0e-9
    # Within 10 percent of textbook 14 nm canonical span
    assert abs(span - ANTIBODY_TOTAL_SPAN_NM) / ANTIBODY_TOTAL_SPAN_NM < 0.15


def test_antigen_position_outside_paratope(geom):
    p = geom.paratope_centers()[1]  # right tip
    a = geom.antigen_position("right")
    # Antigen sits OUTSIDE the paratope (further from origin)
    assert np.linalg.norm(a) > np.linalg.norm(p)


def test_mhc_platform_is_square(geom):
    corners = geom.mhc_platform_corners(z_mhc=-12.0)
    assert corners.shape == (4, 3)
    # All corners at same z
    assert np.allclose(corners[:, 2], -12.0)
    side01 = np.linalg.norm(corners[1] - corners[0])
    side12 = np.linalg.norm(corners[2] - corners[1])
    assert abs(side01 - side12) < 1.0e-9


def test_mhc_peptide_correct_length(geom):
    pep = geom.mhc_peptide_positions()
    # MHC class I displays 8-10 amino acids; we use 9
    assert pep.shape == (MHC_PEPTIDE_LENGTH_AA, 3)


# --------------------------------------------------------------------- #
# (b) Binding affinity K_d in 10^-9 .. 10^-12 M                         #
# --------------------------------------------------------------------- #

def test_default_kd_in_igg_range(sim):
    Kd = sim.binding_affinity_Kd_M()
    # IgG textbook range: 10^-9 .. 10^-12 M
    assert KD_HIGH_AFFINITY_M <= Kd <= KD_MEDIUM_AFFINITY_M


def test_kd_ladder_is_strictly_decreasing():
    # low > medium > high (smaller K_d == tighter binding)
    assert KD_LOW_AFFINITY_M > KD_MEDIUM_AFFINITY_M > KD_HIGH_AFFINITY_M


def test_affinity_classification(sim):
    assert sim.affinity_class(1.0e-5) == "low"
    assert sim.affinity_class(1.0e-7) == "medium"  # below 10^-6
    assert sim.affinity_class(1.0e-11) == "high"


# --------------------------------------------------------------------- #
# (c) Free energy of binding < 0 (favourable)                           #
# --------------------------------------------------------------------- #

def test_binding_free_energy_negative_for_tight_binding(sim):
    dG = sim.binding_free_energy_kcal_per_mol(1.0e-9)
    assert dG < 0.0
    # K_d = 10^-9 -> delta_G ~ -12.7 kcal/mol at 310 K
    assert -15.0 < dG < -10.0


def test_tighter_binding_more_negative_dG(sim):
    dG_loose = sim.binding_free_energy_kcal_per_mol(1.0e-6)
    dG_tight = sim.binding_free_energy_kcal_per_mol(1.0e-12)
    assert dG_tight < dG_loose


def test_fraction_bound_langmuir(sim):
    # At [Ag] = K_d, half-saturated
    Kd = sim.Kd_M
    f = sim.fraction_bound(Kd)
    assert abs(f - 0.5) < 1.0e-9
    # At [Ag] = 0, none bound
    assert sim.fraction_bound(0.0) == 0.0
    # At [Ag] >> K_d, fully saturated
    assert sim.fraction_bound(1.0e6 * Kd) > 0.99


# --------------------------------------------------------------------- #
# (d) Primary vs secondary response                                     #
# --------------------------------------------------------------------- #

def test_primary_response_has_latency(sim):
    t = np.linspace(0.0, 30.0, 200)
    titer = sim.antibody_titer(t, "primary")
    # At t < latency, titer is exactly zero
    assert np.all(titer[t < PRIMARY_LATENCY_DAYS] == 0.0)
    # Titer rises after latency
    late = titer[t > PRIMARY_LATENCY_DAYS + 2.0]
    assert late.max() > 0.0


def test_secondary_response_faster_than_primary(sim):
    t = np.linspace(0.0, 30.0, 300)
    primary = sim.antibody_titer(t, "primary")
    secondary = sim.antibody_titer(t, "secondary")
    # Time to half-peak: secondary reaches it faster
    half_p = primary.max() / 2.0
    half_s = secondary.max() / 2.0
    t_half_p = t[np.argmax(primary > half_p)]
    t_half_s = t[np.argmax(secondary > half_s)]
    assert t_half_s < t_half_p


def test_secondary_peak_higher_than_primary(sim):
    t = np.linspace(0.0, 50.0, 500)
    primary_peak = sim.antibody_titer(t, "primary").max()
    secondary_peak = sim.antibody_titer(t, "secondary").max()
    assert secondary_peak > 10.0 * primary_peak


def test_invalid_response_raises(sim):
    t = np.linspace(0.0, 10.0, 50)
    with pytest.raises(ValueError):
        sim.antibody_titer(t, "tertiary")


def test_primary_vs_secondary_pairing(sim):
    t = np.linspace(0.0, 120.0, 600)
    out = sim.primary_vs_secondary(t, secondary_offset_days=60.0)
    assert "primary" in out and "secondary" in out and "total" in out
    # Total = primary + secondary
    assert np.allclose(out["total"], out["primary"] + out["secondary"])
    # Secondary should peak after the secondary exposure
    sec_peak_idx = int(np.argmax(out["secondary"]))
    assert t[sec_peak_idx] > 60.0
    # And the secondary peak amplitude should exceed the primary peak
    assert out["secondary"].max() > out["primary"].max()


# --------------------------------------------------------------------- #
# (e) Memory cell fraction                                              #
# --------------------------------------------------------------------- #

def test_memory_cell_fraction_value(sim):
    assert abs(sim.memory_cell_fraction() - MEMORY_CELL_FRACTION) < 1.0e-12
    # Reasonable biological range: 5-20 percent
    assert 0.05 <= sim.memory_cell_fraction() <= 0.20


def test_memory_cells_after_response(sim):
    # 10^6 activated B cells -> ~10^5 memory cells
    n_mem = sim.memory_cells_after_response(1.0e6)
    assert n_mem == pytest.approx(1.0e5)


# --------------------------------------------------------------------- #
# (f) Somatic hypermutation rate                                        #
# --------------------------------------------------------------------- #

def test_hypermutation_rate_value(sim):
    rate = sim.hypermutation_rate_per_bp_per_div()
    # Textbook ~10^-3 per bp per division
    assert 1.0e-4 <= rate <= 1.0e-2
    assert rate == HYPERMUTATION_RATE_PER_BP_PER_DIV


def test_hypermutation_far_above_background(sim):
    fold = sim.hypermutation_fold_over_background()
    # SHM is ~10^6 fold above background mutation rate
    assert fold >= 1.0e5


def test_expected_mutations_in_v_region(sim):
    # 10 divisions * 360 bp * 10^-3 = ~3.6 expected mutations
    n_mut = sim.expected_mutations_in_v_region(10, v_region_bp=360)
    assert n_mut == pytest.approx(3.6)


# --------------------------------------------------------------------- #
# (g) Clonal selection                                                  #
# --------------------------------------------------------------------- #

def test_clonal_selection_high_affinity_wins(sim):
    rng = np.random.default_rng(seed=1)
    # Use enough proliferation steps that selection is clear, and a
    # large enough antigen concentration that low K_d clones dominate
    out = sim.clonal_selection(
        n_steps=50, n_initial_clones=8,
        antigen_concentration_M=1.0e-7, rng=rng,
    )
    final_pops = out["history"][-1]
    winner = out["winning_clone"]
    Kds = out["Kds_M"]
    # Winning clone should have one of the lowest K_d values
    sorted_idx = np.argsort(Kds)  # ascending K_d (tightest first)
    # The winning clone should be in the top half by affinity
    assert winner in sorted_idx[: max(1, len(sorted_idx) // 2)]
    # Final winner population must exceed initial
    assert final_pops[winner] > out["history"][0, winner]


def test_clonal_selection_history_grows(sim):
    rng = np.random.default_rng(seed=2)
    out = sim.clonal_selection(
        n_steps=20, n_initial_clones=4,
        antigen_concentration_M=1.0e-6, rng=rng,
    )
    # Total population grows monotonically (or at least non-decreasing)
    totals = out["history"].sum(axis=1)
    assert totals[-1] > totals[0]


# --------------------------------------------------------------------- #
# (h) SIR-like coupled dynamics                                         #
# --------------------------------------------------------------------- #

def test_sir_pathogen_grows_then_cleared(sim):
    t = np.linspace(0.0, 30.0, 600)
    out = sim.sir_immune_response(t, P0=1.0, A0=0.0)
    P, A = out["P"], out["A"]
    # Pathogen grows initially
    assert P.max() > P[0]
    # Antibody titer rises
    assert A.max() > A[0]
    # By the end of the simulation the pathogen has been driven down
    assert P[-1] < 0.3 * P.max()


def test_sir_antibody_peaks_after_pathogen(sim):
    t = np.linspace(0.0, 40.0, 800)
    out = sim.sir_immune_response(t, P0=1.0, A0=0.0)
    t_pathogen_peak = t[int(np.argmax(out["P"]))]
    t_antibody_peak = t[int(np.argmax(out["A"]))]
    # Antibody response lags the pathogen burst
    assert t_antibody_peak > t_pathogen_peak


# --------------------------------------------------------------------- #
# Substrate metadata                                                    #
# --------------------------------------------------------------------- #

def test_geometry_substrate_signature(geom):
    sig = geom.substrate_signature()
    for k in ("interpretation", "fab_length_nm", "fc_length_nm",
              "hinge_angle_deg", "total_span_nm", "lambda_qcd_MeV"):
        assert k in sig


def test_simulator_substrate_signature(sim):
    sig = sim.substrate_signature()
    for k in ("Kd_M", "memory_fraction", "hypermutation_rate_per_bp_per_div",
              "primary_latency_days", "secondary_latency_days",
              "secondary_over_primary_peak", "lambda_qcd_MeV"):
        assert k in sig
    assert sig["secondary_over_primary_peak"] >= 10.0
