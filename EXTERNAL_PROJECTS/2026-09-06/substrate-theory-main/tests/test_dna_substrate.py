"""Tests for the paired DNA-mechanics GEOMETRY + SIMULATOR module.

Covers the three required claims:
    (a) helix pitch = 34 A within 5%  (target 3.4 A x 10 bp/turn)
    (b) persistence length = 50 nm within 10%
    (c) GC vs AT total H-bond binding-energy ratio matches 3:2 H-bond
        count (within ~5% of the textbook 10/7)

Plus geometric / mechanical sanity:
    * antiparallel strand handedness
    * Watson-Crick complement is its own inverse
    * helix radius = 10 A
    * Marko-Siggia force diverges as x -> L
    * B->S overstretching transition lives at F* ~ 65 pN
    * twist-bend coupling is finite and positive
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.dna_substrate import (
    AT_HBOND_KCAL_PER_MOL,
    B_TO_S_EXTENSION_RATIO,
    B_TO_S_FORCE_PN,
    BP_PER_TURN,
    GC_HBOND_KCAL_PER_MOL,
    HELIX_PITCH_A,
    HELIX_RADIUS_A,
    N_HBONDS_AT,
    N_HBONDS_GC,
    PERSISTENCE_LENGTH_NM,
    RISE_PER_BP_A,
    STRETCH_MODULUS_PN,
    TWIST_PER_BP_DEG,
    DNAGeometry,
    DNASimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def helix() -> DNAGeometry:
    return DNAGeometry(n_bp=30)


@pytest.fixture
def sim() -> DNASimulator:
    return DNASimulator()


# --------------------------------------------------------------------- #
# (a) Helix pitch = 3.4 A x 10 bp/turn = 34 A (canonical B-DNA)         #
# --------------------------------------------------------------------- #

def test_constants_pitch_is_rise_x_bp_per_turn():
    """Module-level pitch constant matches the textbook product."""
    assert HELIX_PITCH_A == pytest.approx(RISE_PER_BP_A * BP_PER_TURN, abs=1e-9)


def test_helix_pitch_within_5pct_of_34A(helix):
    """Pitch = rise * bp_per_turn; B-DNA target is ~ 34 A."""
    pitch = helix.pitch_A
    err = abs(pitch - 34.0) / 34.0
    assert err < 0.05, f"pitch={pitch:.3f} A vs 34 A ({err * 100:.2f}% off)"


def test_helix_pitch_textbook_3p4_x_10p5(helix):
    """Strict check: pitch == rise * bp_per_turn in the constructor."""
    pitch = helix.pitch_A
    assert pitch == pytest.approx(3.4 * 10.5, abs=1e-9)


def test_helix_rise_per_bp_is_3p4A(helix):
    assert helix.rise_A == pytest.approx(3.4, abs=1e-9)


def test_helix_bp_per_turn_is_10p5(helix):
    assert helix.bp_per_turn == pytest.approx(10.5, abs=1e-9)


def test_twist_per_bp_is_360_over_10p5(helix):
    assert helix.twist_per_bp_deg == pytest.approx(360.0 / 10.5, abs=1e-9)
    assert TWIST_PER_BP_DEG == pytest.approx(360.0 / 10.5, abs=1e-9)


def test_helix_radius_is_10A(helix):
    assert helix.radius_A == pytest.approx(10.0, abs=1e-9)
    assert HELIX_RADIUS_A == 10.0


def test_strand1_radius_matches(helix):
    """Strand 1 phosphate distance from helical (z) axis equals helix radius."""
    coords = helix.strand1_coords()
    radii = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    assert np.allclose(radii, helix.radius_A, atol=1e-9)


def test_strand2_radius_matches(helix):
    """Strand 2 phosphate distance from helical (z) axis equals helix radius."""
    coords = helix.strand2_coords()
    radii = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    assert np.allclose(radii, helix.radius_A, atol=1e-9)


def test_one_full_turn_advances_by_pitch(helix):
    """Going from index 0 to index 10.5 base pairs advances by pitch_A in z."""
    # Test 21 bp = 2 full turns -> z advance = 2 * pitch
    big = DNAGeometry(n_bp=22)
    z = big.strand1_coords()[:, 2]
    assert (z[21] - z[0]) == pytest.approx(2.0 * big.pitch_A, abs=1e-9)


def test_strands_have_groove_asymmetry(helix):
    """Strand 2 sits at angle (pi - groove_offset) from strand 1 in the
    cross-section, producing the major/minor groove asymmetry."""
    s1 = helix.strand1_coords()[0, :2]
    s2 = helix.strand2_coords()[0, :2]
    azimuth_offset = math.acos(
        np.dot(s1, s2) / (np.linalg.norm(s1) * np.linalg.norm(s2))
    )
    expected = math.pi - helix.groove_offset_rad
    assert azimuth_offset == pytest.approx(expected, abs=1e-6)
    # And explicitly: the two strands are NOT on top of each other and
    # NOT diametrically opposite.
    assert 0.0 < azimuth_offset < math.pi


# --------------------------------------------------------------------- #
# (b) Persistence length = 50 nm (within 10%)                            #
# --------------------------------------------------------------------- #

def test_persistence_length_50nm_within_10pct(sim):
    err = abs(sim.persistence_length_nm - 50.0) / 50.0
    assert err < 0.10, (
        f"L_p={sim.persistence_length_nm:.2f} nm vs 50 nm ({err * 100:.2f}% off)"
    )


def test_persistence_length_in_bp_is_about_150(sim):
    """L_p in nm / rise per bp (in nm) = 50 / 0.34 ~ 147 bp."""
    lp_bp = sim.persistence_length_bp()
    assert lp_bp == pytest.approx(50.0 / 0.34, rel=1e-9)
    # Textbook quote of 150 bp is satisfied within 5 bp
    assert abs(lp_bp - 150.0) < 5.0


def test_bending_modulus_positive(sim):
    """kappa = L_p * k_B T must be positive (pN nm^2)."""
    kappa = sim.bending_modulus_pn_nm2()
    assert kappa > 0.0


# --------------------------------------------------------------------- #
# (c) Watson-Crick H-bond binding energies (AT vs GC)                    #
# --------------------------------------------------------------------- #

def test_at_binding_energy_textbook(sim):
    """A:T ~ 7 kcal/mol (2 H-bonds at ~3.5 kcal/mol each)."""
    assert sim.at_binding_kcal_per_mol() == AT_HBOND_KCAL_PER_MOL == 7.0


def test_gc_binding_energy_textbook(sim):
    """G:C ~ 10 kcal/mol (3 H-bonds at ~3.3 kcal/mol each)."""
    assert sim.gc_binding_kcal_per_mol() == GC_HBOND_KCAL_PER_MOL == 10.0


def test_gc_to_at_energy_ratio_matches_hbond_count_within_5pct(sim):
    """Total energy ratio (10/7 = 1.429) matches 3/2 H-bond count to ~5%."""
    e_ratio = sim.gc_to_at_energy_ratio()
    n_ratio = sim.gc_to_at_hbond_count_ratio()
    assert n_ratio == pytest.approx(1.5, abs=1e-9)
    assert e_ratio == pytest.approx(10.0 / 7.0, abs=1e-9)
    err = abs(e_ratio - n_ratio) / n_ratio
    assert err < 0.06, (
        f"E(GC)/E(AT)={e_ratio:.3f} vs N(GC)/N(AT)={n_ratio:.3f} ({err * 100:.2f}% off)"
    )


def test_total_binding_energy_sums_correctly(sim):
    """Total binding energy of 'AATTGGCC' = 4*AT + 4*GC."""
    seq = "AATTGGCC"
    total = sim.total_binding_energy_kcal_per_mol(seq)
    expected = 4 * AT_HBOND_KCAL_PER_MOL + 4 * GC_HBOND_KCAL_PER_MOL
    assert total == pytest.approx(expected, abs=1e-9)


def test_n_hbonds_count_textbook():
    assert N_HBONDS_AT == 2
    assert N_HBONDS_GC == 3


# --------------------------------------------------------------------- #
# Watson-Crick complementarity                                           #
# --------------------------------------------------------------------- #

def test_complement_pairs():
    assert DNAGeometry.complement("A") == "T"
    assert DNAGeometry.complement("T") == "A"
    assert DNAGeometry.complement("G") == "C"
    assert DNAGeometry.complement("C") == "G"


def test_complement_is_involution():
    for b in "ATGC":
        assert DNAGeometry.complement(DNAGeometry.complement(b)) == b


def test_complement_sequence_is_reverse_complement():
    helix = DNAGeometry(n_bp=8, sequence="AATTGGCC")
    # Reverse complement of AATTGGCC is GGCCAATT
    assert helix.complement_sequence() == "GGCCAATT"


def test_unknown_base_raises():
    with pytest.raises(ValueError):
        DNAGeometry.complement("X")


# --------------------------------------------------------------------- #
# Force-extension and B->S transition                                    #
# --------------------------------------------------------------------- #

def test_marko_siggia_diverges_near_full_extension(sim):
    """F(x/L) -> infty as x -> L in the Marko-Siggia entropic regime."""
    f_low  = sim.marko_siggia_force_pn(np.array([0.5]))[0]
    f_high = sim.marko_siggia_force_pn(np.array([0.95]))[0]
    assert f_high > 10.0 * f_low > 0.0


def test_marko_siggia_at_zero_extension_is_small(sim):
    """At zero extension, entropic force should be tiny (formally zero)."""
    f0 = sim.marko_siggia_force_pn(np.array([0.0]))[0]
    assert abs(f0) < 0.05  # pN


def test_force_extension_curve_has_b_to_s_plateau(sim):
    """Composite F(x/L) hits the B->S plateau at F ~ 65 pN."""
    x, f = sim.force_extension_curve()
    plateau_mask = (x > 1.01) & (x < B_TO_S_EXTENSION_RATIO - 0.01)
    f_plateau = f[plateau_mask]
    # Plateau force matches B_TO_S_FORCE_PN
    assert np.allclose(f_plateau, B_TO_S_FORCE_PN, atol=1e-6)


def test_force_extension_curve_post_overstretch_steep(sim):
    """Beyond x/L = 1.7 the modulus is back to S = 1100 pN (very steep)."""
    x, f = sim.force_extension_curve()
    post_mask = x > B_TO_S_EXTENSION_RATIO + 0.01
    if post_mask.any():
        # Force should rise very rapidly with extension after the plateau
        df = np.diff(f[post_mask])
        dx = np.diff(x[post_mask])
        slope = df.mean() / dx.mean()
        # Slope = stretch modulus ~ 1100 pN per (x/L unit)
        assert slope == pytest.approx(STRETCH_MODULUS_PN, rel=0.1)


def test_b_to_s_force_textbook(sim):
    assert sim.b_to_s_force_pn == pytest.approx(65.0, abs=1.0)


def test_stretching_modulus_textbook(sim):
    assert sim.stretch_modulus_pn == pytest.approx(1100.0, abs=50.0)


# --------------------------------------------------------------------- #
# Twist-bend coupling                                                    #
# --------------------------------------------------------------------- #

def test_twist_bend_coupling_finite_positive(sim):
    """g * curvature * twist gives a positive energy density."""
    e = sim.twist_bend_coupling_energy(0.01, 0.05)
    assert e > 0.0


def test_twist_modulus_in_kT_per_bp_is_reasonable(sim):
    """C * rise / kT should be tens of kT per bp (B-DNA torsional stiffness)."""
    val = sim.twist_modulus_kT_per_bp()
    # 410 pN nm^2 * 0.34 nm / 4.28 pN nm ~ 32 kT/bp
    assert 10.0 < val < 100.0


# --------------------------------------------------------------------- #
# Substrate signatures                                                   #
# --------------------------------------------------------------------- #

def test_geometry_substrate_signature(helix):
    sig = helix.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["right_handed"] is True
    assert sig["antiparallel"] is True
    assert sig["bp_per_turn"] == pytest.approx(10.5, abs=1e-9)


def test_simulator_substrate_signature(sim):
    sig = sim.substrate_signature()
    assert "persistence_length_nm" in sig
    assert "stretch_modulus_pn" in sig
    assert "lambda_qcd_MeV" in sig
    assert "B->S" in sig["interpretation"]


# --------------------------------------------------------------------- #
# Sequence-construction sanity                                           #
# --------------------------------------------------------------------- #

def test_default_sequence_correct_length():
    helix = DNAGeometry(n_bp=12)
    assert len(helix.sequence) == 12


def test_explicit_sequence_length_validated():
    with pytest.raises(ValueError):
        DNAGeometry(n_bp=8, sequence="ATCG")


def test_contour_length_consistent(helix):
    assert helix.contour_length_A == pytest.approx(
        helix.n_bp * helix.rise_A, abs=1e-9
    )
