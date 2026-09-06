"""Tests for the paired ribosome translation GEOMETRY + SIMULATOR module.

Covers the three required claims:
    (a) codon-anticodon Watson-Crick pairing is correctly enforced
    (b) translation rate ~ 20 aa/s in bacteria, ~ 5 aa/s in eukaryotes
    (c) per-codon error rate ~ 1e-4 (kinetic proofreading)

Plus geometric and biochemical sanity:
    * GTP budget = 2 per elongation cycle
    * inter-site spacing = 20 A (A/P, P/E)
    * mRNA codon pitch = 10.2 A (3 nt at 3.4 A)
    * codon counts match (mRNA length / 3)
    * tRNA L-shape arm = 76 A
    * elongation rate distribution mean recovers textbook value
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.ribosome_substrate import (
    CODON_LENGTH_NT,
    CODON_PITCH_A,
    ERROR_RATE_PER_CODON,
    GTP_PER_ELONGATION_CYCLE,
    INITIAL_SELECTION_ERR,
    INTER_SITE_SPACING_A,
    NT_RISE_A,
    PROOFREADING_FACTOR,
    RATE_BACTERIA_AA_PER_S,
    RATE_EUKARYOTE_AA_PER_S,
    SUBUNIT_30S_DIAMETER_A,
    SUBUNIT_30S_HEIGHT_A,
    SUBUNIT_50S_DIAMETER_A,
    SUBUNIT_50S_HEIGHT_A,
    TRNA_LSHAPE_ARM_A,
    WC_BASE_PAIR_DISTANCE_A,
    RibosomeGeometry,
    RibosomeSimulator,
    codon_anticodon_pair,
    watson_crick_complement,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom_bact() -> RibosomeGeometry:
    return RibosomeGeometry(
        organism_class="bacteria",
        mrna_sequence="AUGGCUAAACGUUAA",   # Met-Ala-Lys-Arg-stop (5 codons)
    )


@pytest.fixture
def geom_euk() -> RibosomeGeometry:
    return RibosomeGeometry(
        organism_class="eukaryote",
        mrna_sequence="AUGCCGGAAUGGUAA",   # Met-Pro-Glu-Trp-stop (5 codons)
    )


@pytest.fixture
def sim_bact() -> RibosomeSimulator:
    return RibosomeSimulator(organism_class="bacteria")


@pytest.fixture
def sim_euk() -> RibosomeSimulator:
    return RibosomeSimulator(organism_class="eukaryote")


# --------------------------------------------------------------------- #
# (a) Codon-anticodon Watson-Crick pairing                              #
# --------------------------------------------------------------------- #

def test_watson_crick_complement_basic():
    """Antiparallel reverse-complement: 5'-AUG-3' -> 3'-UAC-5' = "CAU"."""
    assert watson_crick_complement("AUG") == "CAU"
    assert watson_crick_complement("GCU") == "AGC"
    assert watson_crick_complement("AAA") == "UUU"
    assert watson_crick_complement("CCCGGG") == "CCCGGG"  # palindrome


def test_watson_crick_complement_dna_synonym():
    """T should be accepted as a synonym for U (DNA-style input)."""
    assert watson_crick_complement("ATG") == "CAU"


def test_codon_anticodon_pair_start_codon():
    """Start codon AUG (5'->3') pairs with anticodon UAC (read 3'->5').

    Position-by-position antiparallel pairing:
        A (5')-U (3')   U-A   G-C
    """
    assert codon_anticodon_pair("AUG") == "UAC"


def test_codon_anticodon_pair_lengths():
    """Codon must be 3 nt; longer or shorter raises."""
    with pytest.raises(ValueError):
        codon_anticodon_pair("AU")
    with pytest.raises(ValueError):
        codon_anticodon_pair("AUGC")


def test_simulator_recognises_cognate_pairing(sim_bact):
    """A perfect codon-anticodon match has 3 Watson-Crick pairs."""
    codon = "AUG"
    anticodon = codon_anticodon_pair(codon)  # "CAU"
    assert sim_bact.watson_crick_pairs(codon, anticodon) == 3
    assert sim_bact.is_cognate(codon, anticodon) is True


def test_simulator_rejects_near_cognate(sim_bact):
    """A single mismatch -> only 2 of 3 paired -> not cognate.

    Cognate anticodon for AUG is UAC.  Mutating the last position to A
    (UAA) breaks the G-C pair while preserving A-U and U-A.
    """
    codon = "AUG"
    near = "UAA"  # last position mismatched (should be C for cognate)
    assert sim_bact.watson_crick_pairs(codon, near) == 2
    assert sim_bact.is_cognate(codon, near) is False


def test_simulator_rejects_total_mismatch(sim_bact):
    """All three positions wrong -> 0 paired."""
    codon = "AAA"
    anti  = "AAA"  # would need UUU to pair
    assert sim_bact.watson_crick_pairs(codon, anti) == 0


def test_geometry_a_site_anticodon_is_watson_crick(geom_bact):
    """The geometry-derived A-site anticodon must Watson-Crick pair its codon."""
    sim = RibosomeSimulator()
    assert sim.is_cognate(geom_bact.a_site_codon(), geom_bact.a_site_anticodon())


# --------------------------------------------------------------------- #
# (b) Translation rate: bacteria ~ 20 aa/s, eukaryotes ~ 5 aa/s         #
# --------------------------------------------------------------------- #

def test_bacterial_elongation_rate(sim_bact):
    assert sim_bact.elongation_rate_aa_per_second() == pytest.approx(
        RATE_BACTERIA_AA_PER_S, abs=1e-9,
    )
    assert 15.0 <= sim_bact.elongation_rate_aa_per_second() <= 25.0


def test_eukaryotic_elongation_rate(sim_euk):
    assert sim_euk.elongation_rate_aa_per_second() == pytest.approx(
        RATE_EUKARYOTE_AA_PER_S, abs=1e-9,
    )
    assert 3.0 <= sim_euk.elongation_rate_aa_per_second() <= 8.0


def test_bacterial_dwell_time_per_codon(sim_bact):
    """20 aa/s -> 50 ms per codon."""
    assert sim_bact.time_per_codon_seconds() == pytest.approx(0.05, abs=1e-6)


def test_eukaryotic_dwell_time_per_codon(sim_euk):
    """5 aa/s -> 200 ms per codon."""
    assert sim_euk.time_per_codon_seconds() == pytest.approx(0.20, abs=1e-6)


def test_translation_rate_distribution_mean(sim_bact):
    """Sampled distribution mean recovers the textbook arithmetic mean."""
    samples = sim_bact.elongation_rate_distribution(n_samples=20000, seed=42)
    assert samples.mean() == pytest.approx(RATE_BACTERIA_AA_PER_S, rel=0.05)
    assert samples.min() > 0.0


def test_bacteria_eukaryote_speed_ratio(sim_bact, sim_euk):
    """Bacteria translate ~ 4x faster than eukaryotes."""
    ratio = sim_bact.elongation_rate_aa_per_second() / sim_euk.elongation_rate_aa_per_second()
    assert ratio == pytest.approx(4.0, abs=0.1)


def test_translation_time_for_300aa_protein(sim_bact):
    """300 residues at 20 aa/s = 15 s in bacteria."""
    assert sim_bact.time_to_translate_protein_seconds(300) == pytest.approx(
        15.0, abs=1e-9,
    )


# --------------------------------------------------------------------- #
# (c) Error rate ~ 1e-4 from kinetic proofreading                       #
# --------------------------------------------------------------------- #

def test_per_codon_error_rate_is_1e4(sim_bact):
    assert sim_bact.per_codon_error_rate() == pytest.approx(1e-4, rel=1e-9)


def test_proofreading_two_gates_multiply(sim_bact):
    """Combined error = initial_selection * proofreading ~ 1e-4."""
    g1 = sim_bact.initial_selection_error()
    g2 = sim_bact.proofreading_factor()
    combined = sim_bact.proofreading_combined_error()
    assert combined == pytest.approx(g1 * g2, rel=1e-9)
    # Each gate ~ 1e-2, combined ~ 1e-4
    assert g1 == pytest.approx(1e-2, rel=1e-9)
    assert g2 == pytest.approx(1e-2, rel=1e-9)
    assert combined == pytest.approx(1e-4, rel=1e-9)


def test_per_codon_error_matches_proofreading_combined(sim_bact):
    """The combined two-gate error rate should equal the per-codon error rate."""
    assert sim_bact.per_codon_error_rate() == pytest.approx(
        sim_bact.proofreading_combined_error(), rel=1e-9,
    )


def test_errors_per_300aa_protein(sim_bact):
    """At 1e-4 per codon, a 300-residue protein has ~ 0.03 errors on average."""
    assert sim_bact.errors_per_protein(300) == pytest.approx(0.03, abs=1e-9)


def test_perfect_protein_fraction(sim_bact):
    """P(perfect) for a 300-residue protein = (1 - 1e-4)^300 ~ 0.97."""
    p_perfect = sim_bact.fraction_perfect_proteins(300)
    assert p_perfect == pytest.approx(0.97, abs=0.01)


def test_perfect_protein_fraction_decays_with_length(sim_bact):
    """Longer proteins have a lower per-copy perfect fraction."""
    p100  = sim_bact.fraction_perfect_proteins(100)
    p1000 = sim_bact.fraction_perfect_proteins(1000)
    assert p100 > p1000
    # 1000 residues at 1e-4 -> exp(-0.1) ~ 0.90
    assert p1000 == pytest.approx(math.exp(-0.1), abs=0.01)


# --------------------------------------------------------------------- #
# GTP hydrolysis budget                                                  #
# --------------------------------------------------------------------- #

def test_gtp_per_elongation_cycle(sim_bact):
    """Two GTPs per cycle: EF-Tu (delivery + proofreading) + EF-G (translocation)."""
    assert sim_bact.gtp_per_cycle() == GTP_PER_ELONGATION_CYCLE
    assert sim_bact.gtp_per_cycle() == 2


def test_gtp_per_protein_scales_linearly(sim_bact):
    """A 300-residue protein consumes 600 GTP molecules."""
    assert sim_bact.gtp_per_protein(300) == 600


def test_energy_per_peptide_bond_in_kT(sim_bact):
    """2 GTP * 7.3 kcal/mol / kT(310 K) ~ 23 kT per peptide bond."""
    e_kT = sim_bact.energy_per_peptide_bond_kT()
    assert 20.0 <= e_kT <= 26.0


# --------------------------------------------------------------------- #
# Geometry: subunits, mRNA track, A/P/E sites                           #
# --------------------------------------------------------------------- #

def test_bacterial_subunit_dimensions(geom_bact):
    """30S = 220 A diameter, 100 A height; 50S = 250 A, 150 A."""
    d_small, h_small = geom_bact.small_subunit_dimensions_angstrom()
    d_large, h_large = geom_bact.large_subunit_dimensions_angstrom()
    assert d_small == pytest.approx(SUBUNIT_30S_DIAMETER_A, rel=1e-9)
    assert h_small == pytest.approx(SUBUNIT_30S_HEIGHT_A, rel=1e-9)
    assert d_large == pytest.approx(SUBUNIT_50S_DIAMETER_A, rel=1e-9)
    assert h_large == pytest.approx(SUBUNIT_50S_HEIGHT_A, rel=1e-9)


def test_eukaryotic_subunits_are_larger(geom_euk):
    """40S/60S are 7-10% larger than 30S/50S."""
    d_small, h_small = geom_euk.small_subunit_dimensions_angstrom()
    d_large, h_large = geom_euk.large_subunit_dimensions_angstrom()
    assert d_small > SUBUNIT_30S_DIAMETER_A
    assert h_small > SUBUNIT_30S_HEIGHT_A
    assert d_large > SUBUNIT_50S_DIAMETER_A
    assert h_large > SUBUNIT_50S_HEIGHT_A


def test_codon_count_matches_mrna_length(geom_bact):
    """5-codon mRNA -> 5 codon positions."""
    assert geom_bact.codon_count() == 5
    assert geom_bact.codon_positions().shape == (5, 3)


def test_codon_pitch_is_10_2_angstrom():
    """3 nt at 3.4 A/nt = 10.2 A per codon."""
    assert CODON_PITCH_A == pytest.approx(NT_RISE_A * CODON_LENGTH_NT, rel=1e-9)
    assert CODON_PITCH_A == pytest.approx(10.2, rel=1e-9)


def test_codon_positions_are_evenly_spaced(geom_bact):
    """Consecutive codons are ~ 10.2 A apart along the mRNA track."""
    coords = geom_bact.codon_positions()
    spacings = np.linalg.norm(np.diff(coords, axis=0), axis=1)
    assert np.allclose(spacings, CODON_PITCH_A, atol=1e-6)


def test_a_p_e_site_spacings(geom_bact):
    """A-P and P-E are each 20 A centre-to-centre."""
    sites = geom_bact.site_centres()
    ap = float(np.linalg.norm(sites["A"] - sites["P"]))
    pe = float(np.linalg.norm(sites["P"] - sites["E"]))
    assert ap == pytest.approx(INTER_SITE_SPACING_A, rel=1e-9)
    assert pe == pytest.approx(INTER_SITE_SPACING_A, rel=1e-9)


def test_a_site_codon_is_second_codon(geom_bact):
    """At initiation, AUG is at P; A site holds the second codon."""
    assert geom_bact.a_site_codon_index() == 1
    assert geom_bact.a_site_codon() == geom_bact.codon_at(1)


def test_codon_at_indexing(geom_bact):
    """5-codon mRNA: indices 0..4 valid, 5 raises."""
    for i in range(5):
        assert len(geom_bact.codon_at(i)) == CODON_LENGTH_NT
    with pytest.raises(IndexError):
        geom_bact.codon_at(5)


def test_trna_lshape_arm_length():
    """tRNA L-shape arm acceptor-to-anticodon is ~ 76 A."""
    assert TRNA_LSHAPE_ARM_A == pytest.approx(76.0, rel=1e-9)


def test_watson_crick_base_pair_distance():
    """Standard N..N donor-acceptor distance in WC pairs ~ 2.85 A."""
    assert WC_BASE_PAIR_DISTANCE_A == pytest.approx(2.85, rel=1e-9)


def test_subunits_are_stacked(geom_bact):
    """Large subunit centre is above small subunit centre (in z)."""
    z_small = geom_bact.small_subunit_centre()[2]
    z_large = geom_bact.large_subunit_centre()[2]
    assert z_large > z_small


# --------------------------------------------------------------------- #
# Substrate identity and metadata                                       #
# --------------------------------------------------------------------- #

def test_geometry_substrate_signature(geom_bact):
    sig = geom_bact.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["subunit_count"] == 2
    assert sig["trna_sites"] == 3
    assert sig["gtp_per_cycle"] == 2


def test_simulator_substrate_signature(sim_bact):
    sig = sim_bact.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert "elongation_rate_aa_per_s" in sig
    assert "per_codon_error_rate" in sig
    assert sig["gtp_per_cycle"] == 2
    assert "mass-torque" in sig["interpretation"]
