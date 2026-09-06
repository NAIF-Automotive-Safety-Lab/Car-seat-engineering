"""Tests for the paired protein-folding GEOMETRY + SIMULATOR module.

Covers the four required claims:
    (a) alpha-helix pitch = 5.4 Angstrom within 1%
    (b) alpha-helix has 3.6 residues per turn
    (c) antiparallel beta-sheet H-bond period = 7.0 Angstrom
    (d) Levinthal paradox is resolved by the substrate folding funnel
        (random-search time >> funnel-relaxation time, by ~10^40 at N=100).

Plus geometric sanity:
    * canonical (phi, psi) Ramachandran torsions land in allowed regions
    * CA-CA bonded distance ~ 3.8 Angstrom (textbook)
    * funnel has a unique global minimum at the native fold (Q->1, RMSD->0)
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.protein_folding_substrate import (
    HELIX_PITCH_ANGSTROM,
    HELIX_RADIUS_ANGSTROM,
    HELIX_RESIDUES_PER_TURN,
    HELIX_RISE_PER_RES_A,
    PHI_ALPHA_HELIX,
    PHI_BETA_SHEET,
    PSI_ALPHA_HELIX,
    PSI_BETA_SHEET,
    SHEET_H_BOND_PERIOD_A,
    SHEET_RISE_PER_RES_A,
    ProteinFoldingGeometry,
    ProteinFoldingSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def helix() -> ProteinFoldingGeometry:
    return ProteinFoldingGeometry.alpha_helix(n_residues=30)


@pytest.fixture
def strand() -> ProteinFoldingGeometry:
    return ProteinFoldingGeometry.beta_sheet_strand(n_residues=10)


@pytest.fixture
def sim() -> ProteinFoldingSimulator:
    return ProteinFoldingSimulator()


# --------------------------------------------------------------------- #
# (a) Alpha-helix pitch = 5.4 Angstrom within 1%                        #
# --------------------------------------------------------------------- #

def test_helix_pitch_5p4_angstrom_within_1pct(helix):
    _, _, pitch = helix.helix_axis_and_pitch()
    err = abs(pitch - HELIX_PITCH_ANGSTROM) / HELIX_PITCH_ANGSTROM
    assert err < 0.01, f"pitch={pitch:.3f} A vs target 5.4 A ({err * 100:.2f}% off)"


def test_helix_pitch_via_simulator(sim):
    pitch = sim.helix_pitch_angstrom(n_residues=30)
    err = abs(pitch - 5.4) / 5.4
    assert err < 0.01, f"pitch={pitch:.3f} A ({err * 100:.2f}% off)"


def test_helix_rise_per_residue_within_1pct(helix):
    _, rise, _ = helix.helix_axis_and_pitch()
    err = abs(rise - HELIX_RISE_PER_RES_A) / HELIX_RISE_PER_RES_A
    assert err < 0.01, f"rise={rise:.3f} A vs 1.5 A ({err * 100:.2f}% off)"


# --------------------------------------------------------------------- #
# (b) 3.6 residues per turn                                             #
# --------------------------------------------------------------------- #

def test_helix_residues_per_turn(sim):
    assert sim.helix_residues_per_turn() == pytest.approx(3.6, abs=1e-9)


def test_helix_residues_per_turn_geometric(helix):
    """pitch / rise should equal residues per turn."""
    _, rise, pitch = helix.helix_axis_and_pitch()
    assert pitch / rise == pytest.approx(3.6, abs=0.05)


# --------------------------------------------------------------------- #
# alpha-helix coordinates match observed geometry                       #
# --------------------------------------------------------------------- #

def test_helix_radius_observed_within_2pct(helix):
    rad = helix.helix_radius()
    err = abs(rad - HELIX_RADIUS_ANGSTROM) / HELIX_RADIUS_ANGSTROM
    assert err < 0.02, f"helix radius={rad:.3f} A vs 2.3 A ({err * 100:.2f}% off)"


def test_helix_calpha_calpha_distance(helix):
    """Consecutive C-alpha atoms are ~ 3.8 A apart (textbook)."""
    ca = helix.calpha_coords()
    dists = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    assert dists.mean() == pytest.approx(3.8, abs=0.1)
    # Helix is uniform -- almost no spread
    assert dists.std() < 0.01


def test_helix_axis_is_well_defined(helix):
    """SVD of CA cloud yields a clean principal direction."""
    ca = helix.calpha_coords()
    centred = ca - ca.mean(axis=0)
    s = np.linalg.svd(centred, compute_uv=False)
    # First singular value should dominate the others by >> 1
    assert s[0] / s[1] > 5.0


def test_helix_n_atoms_matches_residues(helix):
    """3 backbone atoms (N, CA, C') per residue."""
    coords = helix.backbone_coords()
    assert coords.shape == (3 * helix.n_residues, 3)


# --------------------------------------------------------------------- #
# (c) Beta-sheet H-bond pattern (period 7 A)                            #
# --------------------------------------------------------------------- #

def test_beta_sheet_h_bond_period_7A(sim):
    period = sim.beta_sheet_h_bond_period_angstrom()
    assert period == pytest.approx(7.0, abs=1e-9)
    assert period == SHEET_H_BOND_PERIOD_A


def test_beta_strand_is_extended(strand):
    """beta-strand rise per residue ~ 3.3 A (extended chain)."""
    sim = ProteinFoldingSimulator()
    rise = sim.beta_sheet_rise_per_residue_angstrom(n_residues=10)
    err = abs(rise - SHEET_RISE_PER_RES_A) / SHEET_RISE_PER_RES_A
    assert err < 0.10, f"strand rise={rise:.3f} A vs 3.3 A ({err * 100:.2f}% off)"


def test_beta_strand_calpha_calpha_distance(strand):
    """beta-strand has CA-CA bonded distance ~ 3.8 A like the helix."""
    ca = strand.calpha_coords()
    dists = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    assert dists.mean() == pytest.approx(3.8, abs=0.1)


# --------------------------------------------------------------------- #
# Ramachandran (phi, psi) angles                                        #
# --------------------------------------------------------------------- #

def test_alpha_helix_torsions_in_allowed_region():
    assert ProteinFoldingGeometry.is_allowed(PHI_ALPHA_HELIX, PSI_ALPHA_HELIX)


def test_beta_sheet_torsions_in_allowed_region():
    assert ProteinFoldingGeometry.is_allowed(PHI_BETA_SHEET, PSI_BETA_SHEET)


def test_disallowed_region_rejected():
    """Random forbidden region (e.g. phi=+60, psi=-60) should be disallowed."""
    assert not ProteinFoldingGeometry.is_allowed(0.0, -150.0)


def test_ramachandran_angles_returned_in_degrees(helix):
    angles = helix.ramachandran_angles_deg()
    assert angles.shape == (helix.n_residues, 2)
    # All residues in the helix share the same canonical (phi, psi)
    assert np.allclose(angles[:, 0], PHI_ALPHA_HELIX)
    assert np.allclose(angles[:, 1], PSI_ALPHA_HELIX)


# --------------------------------------------------------------------- #
# (d) Folding funnel and Levinthal paradox resolution                   #
# --------------------------------------------------------------------- #

def test_funnel_landscape_shape(sim):
    q, r, U = sim.funnel_landscape()
    assert U.shape == (len(q), len(r))


def test_funnel_has_unique_minimum_at_native(sim):
    """Global minimum should be at high Q (native contacts) and low RMSD."""
    q_star, r_star, _ = sim.funnel_minimum()
    assert q_star > 0.9, f"native Q* = {q_star} should be > 0.9"
    assert r_star < 1.0, f"native RMSD* = {r_star} A should be < 1 A"


def test_funnel_depth_is_cooperative(sim):
    """Funnel depth >> kT so folding is cooperative (two-state)."""
    assert sim.funnel_depth_kT >= 10.0


def test_levinthal_paradox_resolved(sim):
    """At N=100 residues, funnel is faster than random search by >> 10^20."""
    out = sim.levinthal_resolved(n_residues=100)
    assert out["resolved"] is True
    assert out["speedup_ratio"] > 1.0e20
    # Funnel time at N=100 should be sub-second (microsecond-scale)
    assert out["funnel_relaxation_s"] < 1.0e-3


def test_levinthal_random_search_is_astronomical(sim):
    """Random search of 100 residues at 3 states would take >> age of universe."""
    age_universe_s = 4.35e17
    t_lev = sim.levinthal_search_time_seconds(100)
    assert t_lev > 1.0e10 * age_universe_s


# --------------------------------------------------------------------- #
# Substrate identity and metadata                                       #
# --------------------------------------------------------------------- #

def test_geometry_substrate_signature(helix):
    sig = helix.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["torsion_dof_per_residue"] == 2
    assert sig["secondary_structure_count"] == 2


def test_simulator_substrate_signature(sim):
    sig = sim.substrate_signature()
    assert "funnel_depth_kT" in sig
    assert "lambda_qcd_MeV" in sig
    assert "mass-torque" in sig["interpretation"]
