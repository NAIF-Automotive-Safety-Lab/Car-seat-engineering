"""Tests for stiff_medium.molecular_bond_substrate (paired GEOM+SIM+VIZ).

Verifies:
  (a) H_2 bond length = 0.74 Å, dissociation energy 432 kJ/mol (5%)
  (b) H_2O bond angle = 104.5° (1%)
  (c) CH_4 all six H–C–H angles = arccos(-1/3) = 109.47° (substrate sp^3)
  (d) NH_3 H–N–H angle compressed to ~107° (VSEPR)
  (e) CO_2 linear (180°) with symmetric stretch ν_1 ≈ 1330 cm⁻¹
  (f) Tetrahedral angle = K_4 deuteron face-pair angle
  (g) Common-bonds catalogue contains ≥ 20 entries with sane values
  (h) Visualizer / benzene renderer don't crash on a matplotlib axis
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from stiff_medium.molecular_bond_substrate import (
    ATOM_STYLE,
    COMMON_BONDS,
    MoleculeGeometry,
    MolecularBondSimulator,
    MoleculeVisualizer,
    REFERENCE_ENERGIES,
    REFERENCE_GEOMETRY,
    TETRAHEDRAL_ANGLE_DEG,
    all_small_molecules,
    benzene_positions_ang,
    render_benzene,
)


# ---------------------------------------------------------------------------
# (a) H_2 bond length and dissociation energy
# ---------------------------------------------------------------------------

def test_h2_bond_length_074() -> None:
    g = MoleculeGeometry("H2")
    assert g.bond_length == pytest.approx(0.7414, rel=1e-3)
    # Distance between the two H atoms
    d = float(np.linalg.norm(g.positions[1] - g.positions[0]))
    assert d == pytest.approx(0.7414, rel=1e-3)


def test_h2_dissociation_energy_5pct() -> None:
    sim = MolecularBondSimulator(MoleculeGeometry("H2"))
    BE = sim.bond_energy_kJmol()
    # Reference: H_2 BE = 432 kJ/mol (slight variation 432–436 in tables).
    assert abs(BE - 432.0) / 432.0 < 0.05, f"H2 BE = {BE} kJ/mol"


# ---------------------------------------------------------------------------
# (b) H_2O bond angle
# ---------------------------------------------------------------------------

def test_h2o_angle_1045() -> None:
    g = MoleculeGeometry("H2O")
    angle = g.all_bond_angles()[0]
    assert abs(angle - 104.5) < 0.01 * 104.5
    # Bond length matches NIST O–H = 0.9584 Å
    assert g.bond_length == pytest.approx(0.9584, rel=1e-3)


# ---------------------------------------------------------------------------
# (c) CH_4 — every H–C–H angle is the K_4 tetrahedral angle
# ---------------------------------------------------------------------------

def test_ch4_all_six_angles_tetrahedral() -> None:
    g = MoleculeGeometry("CH4")
    angles = g.all_bond_angles()
    assert len(angles) == 6                                  # C(4,2) = 6 pairs
    for a in angles:
        assert a == pytest.approx(TETRAHEDRAL_ANGLE_DEG, abs=1e-3), a
    # The angle is exactly arccos(-1/3) = 109.47°, the SAME as the K_4
    # deuteron face-pair tetrahedral angle.
    assert TETRAHEDRAL_ANGLE_DEG == pytest.approx(109.4712, abs=1e-3)


def test_ch4_all_four_bond_lengths_equal() -> None:
    g = MoleculeGeometry("CH4")
    lens = g.all_bond_lengths()
    assert len(lens) == 4
    assert max(lens) - min(lens) < 1e-6                       # all equal
    assert lens[0] == pytest.approx(g.bond_length, rel=1e-6)


# ---------------------------------------------------------------------------
# (d) NH_3 — VSEPR-compressed angle
# ---------------------------------------------------------------------------

def test_nh3_compressed_angle() -> None:
    g = MoleculeGeometry("NH3")
    angles = g.all_bond_angles()
    assert len(angles) == 3                                   # C(3,2) = 3 pairs
    for a in angles:
        assert abs(a - 106.7) < 1.0
    # Compressed below tetrahedral
    assert angles[0] < TETRAHEDRAL_ANGLE_DEG


# ---------------------------------------------------------------------------
# (e) CO_2 — linear with symmetric stretch ν_1 ≈ 1333 cm⁻¹
# ---------------------------------------------------------------------------

def test_co2_linear_geometry() -> None:
    g = MoleculeGeometry("CO2")
    angle = g.all_bond_angles()[0]
    assert angle == pytest.approx(180.0, abs=1e-6)
    lens = g.all_bond_lengths()
    assert lens[0] == pytest.approx(lens[1], rel=1e-9)        # symmetric
    assert lens[0] == pytest.approx(1.16, rel=1e-3)


def test_co2_symmetric_stretch_1330_cm() -> None:
    sim = MolecularBondSimulator(MoleculeGeometry("CO2"))
    nu1 = sim.co2_modes_cm()["nu1_sym_stretch"]
    # NIST: ν_1 = 1333 cm^-1  (test allows ±5 cm^-1)
    assert abs(nu1 - 1333.0) < 5.0
    # Asymmetric stretch is higher than symmetric
    nu3 = sim.co2_modes_cm()["nu3_asym_stretch"]
    assert nu3 > nu1
    assert abs(nu3 - 2349.0) < 5.0


# ---------------------------------------------------------------------------
# (f) sp^3 hybridization angle ≡ K_4 deuteron face-pair angle
# ---------------------------------------------------------------------------

def test_sp3_equals_arccos_minus_one_third() -> None:
    """Substrate-derived: 109.47° = arccos(-1/3) is the K_4 tetrahedral angle."""
    a = MoleculeGeometry.sp3_tetrahedral_angle()
    assert a == pytest.approx(math.degrees(math.acos(-1.0 / 3.0)))
    assert a == pytest.approx(109.47122, abs=1e-4)
    # Identity with the deuteron face-pair vertex-centre-vertex angle:
    from stiff_medium.k4_face_pair_geometry import TETRAHEDRAL_ANGLE_DEG as K4_TET
    assert a == pytest.approx(K4_TET, abs=1e-6)


# ---------------------------------------------------------------------------
# (g) Common-bonds catalogue
# ---------------------------------------------------------------------------

def test_common_bonds_catalogue_size() -> None:
    assert len(COMMON_BONDS) >= 20
    for name, length, energy in COMMON_BONDS:
        assert 0.5 < length < 2.5,    (name, length)          # reasonable Å
        assert 100.0 < energy < 1100.0, (name, energy)        # reasonable kJ/mol


def test_common_bonds_correlate_short_strong() -> None:
    """Stronger bonds should be on average shorter (Badger-rule trend)."""
    lengths  = np.array([b[1] for b in COMMON_BONDS])
    energies = np.array([b[2] for b in COMMON_BONDS])
    # Pearson r should be negative (anti-correlation)
    r = np.corrcoef(lengths, energies)[0, 1]
    assert r < -0.4, f"Pearson correlation = {r:.3f}"


# ---------------------------------------------------------------------------
# (h) Vibrational frequencies of H_2O and H_2 are accurate
# ---------------------------------------------------------------------------

def test_h2o_stretch_and_bend_match_reference() -> None:
    sim = MolecularBondSimulator(MoleculeGeometry("H2O"))
    nu_stretch = sim.stretch_frequency_cm()
    nu_bend = sim.bend_frequency_cm()
    assert nu_stretch == pytest.approx(3657.0, rel=1e-3)
    assert nu_bend    == pytest.approx(1595.0, rel=1e-3)


def test_h2_stretch_4401() -> None:
    sim = MolecularBondSimulator(MoleculeGeometry("H2"))
    nu_stretch = sim.stretch_frequency_cm()
    assert nu_stretch == pytest.approx(4401.0, rel=1e-3)


def test_h2_zero_point_energy_eV() -> None:
    """ZPE of H_2 stretch is well-known: ½ ℏω ≈ 0.273 eV."""
    sim = MolecularBondSimulator(MoleculeGeometry("H2"))
    zpe = sim.zero_point_energy_eV()
    assert abs(zpe - 0.273) < 0.005


# ---------------------------------------------------------------------------
# (i) Morse potential goes to −D_e at r_e and 0 at large r
# ---------------------------------------------------------------------------

def test_morse_minimum_at_re() -> None:
    sim = MolecularBondSimulator(MoleculeGeometry("H2"))
    r_e = sim.geometry.bond_length
    V_at_re = float(sim.morse_potential(np.array([r_e]))[0])
    assert V_at_re == pytest.approx(-sim.bond_energy_kJmol(), rel=1e-9)
    # Far away → 0
    V_far = float(sim.morse_potential(np.array([20.0]))[0])
    assert abs(V_far) < 0.5                                   # kJ/mol


# ---------------------------------------------------------------------------
# (j) Builders / convenience
# ---------------------------------------------------------------------------

def test_all_small_molecules_dict_is_complete() -> None:
    d = all_small_molecules()
    assert set(d) == {"H2", "H2O", "CH4", "NH3", "CO2"}
    for name, geom in d.items():
        assert isinstance(geom, MoleculeGeometry)
        assert geom.species == name


def test_unknown_species_raises() -> None:
    with pytest.raises(ValueError):
        MoleculeGeometry("HCN")


# ---------------------------------------------------------------------------
# (k) Visualizer does not crash on a real matplotlib 3D axis
# ---------------------------------------------------------------------------

def test_visualizer_renders_each_molecule() -> None:
    fig = plt.figure(figsize=(4, 4))
    for spec in REFERENCE_GEOMETRY:
        fig.clear()
        ax = fig.add_subplot(111, projection="3d")
        viz = MoleculeVisualizer(MoleculeGeometry(spec))
        viz.render_ball_and_stick(ax)
    plt.close(fig)


def test_benzene_geometry_and_render() -> None:
    pos, atoms = benzene_positions_ang()
    assert len(atoms) == 12                                   # 6 C + 6 H
    assert atoms.count("C") == 6 and atoms.count("H") == 6
    # All carbons on a hexagon of radius 1.397 Å
    rs = np.linalg.norm(pos[:6], axis=1)
    assert np.all(np.abs(rs - 1.397) < 1e-6)
    # Render without crashing
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection="3d")
    render_benzene(ax)
    plt.close(fig)


def test_atom_style_has_required_elements() -> None:
    for el in ("H", "C", "N", "O"):
        assert el in ATOM_STYLE
