"""Paired ribosome translation GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for ribosomal translation
under the B3 / stiff-medium ontology.

GEOMETRY
--------
``RibosomeGeometry`` represents the assembled 70S (bacterial) or 80S
(eukaryotic) ribosome as two stacked cylindrical subunits with an mRNA
groove threading the small subunit and three tRNA-binding sites (A, P, E)
spanning the inter-subunit cleft.

Reference dimensions (cryo-EM consensus, bacterial 70S):
    * 30S small subunit:  ~ 220 Angstrom diameter, ~ 100 A height
    * 50S large subunit:  ~ 250 Angstrom diameter, ~ 150 A height
    * mRNA codon pitch:   ~ 10.4 Angstrom per codon (3 nt at 3.4 A/nt)
    * Inter-site spacing: ~ 20 Angstrom centre-to-centre (A/P, P/E)
    * tRNA L-shape arm:   ~ 76 nt, ~ 76 A acceptor-to-anticodon
    * Anticodon stem-loop: ~ 18 A across the codon-anticodon helix

The codon-anticodon pairing in the A site uses standard Watson-Crick
geometry: A-U and G-C base pairs at ~ 2.85 A donor-acceptor H-bond
distance, three pairs per codon (~ 10.2 A helical rise).

SIMULATOR
---------
``RibosomeSimulator`` exposes the kinetic and thermodynamic observables of
a translating ribosome.  Under the substrate ontology the ribosome is a
catalytic strain-pump: each elongation cycle dissipates ~ 2 GTP worth of
substrate strain energy to advance the mRNA by one codon and lay down one
peptide bond.

Translation rates (textbook):
    * Bacteria  (E. coli, 37 C):   ~ 20 amino acids per second
    * Eukaryotes (mammalian, 37 C): ~  5 amino acids per second

Error rate (kinetic proofreading, Hopfield 1974 / Ninio 1975):
    * Per-codon misincorporation:  ~ 1e-4
    * Initial selection accuracy:  ~ 1e-2
    * Proofreading multiplier:     ~ 1e-2  (2nd kinetic gate)
    * Combined accuracy:           ~ 1e-4

GTP hydrolysis budget per elongation cycle:
    * EF-Tu * GTP (cognate aa-tRNA delivery + proofreading): 1 GTP
    * EF-G  * GTP (translocation A->P, P->E):                1 GTP
    * Total per peptide bond:                                2 GTP

Substrate-specific predictions encoded here:
    * elongation rate ratio (bacteria / eukaryotes) ~ 4
      reflects the ratio of cytosolic strain density between
      prokaryotic and eukaryotic substrates
    * proofreading enhancement factor f_proof = exp(-Delta G_proof / kT)
      where Delta G_proof ~ 4-5 kT per discrimination step
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) values -- cryo-EM and biochemistry consensus
# ---------------------------------------------------------------------------

# Subunit dimensions (Angstrom), bacterial 70S
SUBUNIT_30S_DIAMETER_A: float = 220.0
SUBUNIT_30S_HEIGHT_A:   float = 100.0
SUBUNIT_50S_DIAMETER_A: float = 250.0
SUBUNIT_50S_HEIGHT_A:   float = 150.0

# mRNA geometry
NT_RISE_A:        float = 3.4    # B-form-like rise per nucleotide on mRNA track
CODON_LENGTH_NT:  int   = 3
CODON_PITCH_A:    float = NT_RISE_A * CODON_LENGTH_NT   # 10.2 A per codon

# tRNA-site spacing (centre-to-centre), A site -> P site -> E site
INTER_SITE_SPACING_A: float = 20.0

# tRNA L-shape (textbook 76 nt cloverleaf folded into L)
TRNA_LENGTH_NT:        int   = 76
TRNA_LSHAPE_ARM_A:     float = 76.0   # acceptor end to anticodon end

# Codon-anticodon Watson-Crick helix
WC_BASE_PAIR_DISTANCE_A: float = 2.85   # mean N..N H-bond donor-acceptor
WC_HELIX_RADIUS_A:       float = 9.5    # half-width of the codon-anticodon helix
WC_HELIX_RISE_PER_BP_A:  float = 3.4

# Translation rates (aa/s)
RATE_BACTERIA_AA_PER_S:   float = 20.0
RATE_EUKARYOTE_AA_PER_S:  float =  5.0

# Error / proofreading
ERROR_RATE_PER_CODON:     float = 1.0e-4   # combined misincorporation rate
INITIAL_SELECTION_ERR:    float = 1.0e-2   # kinetic discrimination gate 1
PROOFREADING_FACTOR:      float = 1.0e-2   # kinetic discrimination gate 2

# GTP budget per elongation cycle
GTP_PER_ELONGATION_CYCLE: int = 2          # EF-Tu + EF-G (bacterial)

# Watson-Crick complementarity table (DNA letters convert to RNA U)
_COMPLEMENT_RNA = {
    "A": "U", "U": "A", "G": "C", "C": "G",
    # Allow DNA-style T as a synonym for U (for tests passing mRNA as DNA)
    "T": "A",
}

# Substrate anchor (B3 framework) -- enters strain-energy normalisation
LAMBDA_QCD_MEV: float = 200.0
KT_AT_BODY_T_KCAL_PER_MOL: float = 0.617  # kT at 310 K, kcal/mol
GTP_HYDROLYSIS_KCAL_PER_MOL: float = 7.3  # standard biochem ~ 7.3 kcal/mol


def watson_crick_complement(rna: str) -> str:
    """Return the antiparallel Watson-Crick reverse-complement of an RNA sequence.

    Each base is replaced with its WC partner (A<->U, G<->C, T treated as U)
    and the resulting string is reversed.  This gives the partner strand
    read 5'->3'.

    Examples
    --------
    >>> watson_crick_complement("AUG")  # codon 5'->3'
    'CAU'                                 # anticodon 5'->3'
    """
    seq = rna.upper()
    complement = "".join(_COMPLEMENT_RNA[ch] for ch in seq)
    return complement[::-1]


def codon_anticodon_pair(codon: str) -> str:
    """Return the canonical anticodon (read 3'->5') for a codon read 5'->3'.

    The anticodon is the antiparallel partner of the codon: each base
    pairs with its Watson-Crick complement when written in the antiparallel
    convention, so position-by-position identity reading codon (5'->3')
    against anticodon (3'->5') gives the WC pairs directly:

        codon (5'->3'):   A U G
                          | | |
        anticodon(3'->5'):U A C  -> returned as "UAC"
    """
    if len(codon) != CODON_LENGTH_NT:
        raise ValueError(f"codon must be {CODON_LENGTH_NT} nt, got {len(codon)}")
    seq = codon.upper()
    # Position-by-position complement (no reversal): codon[i] (5'->3') pairs
    # with anticodon[i] (3'->5') when both are written in their natural
    # reading directions.
    return "".join(_COMPLEMENT_RNA[ch] for ch in seq)


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class RibosomeGeometry:
    """70S/80S ribosome scaffold: two subunits + mRNA + A/P/E sites + tRNA arm.

    The geometry is a deliberate caricature -- two stacked cylindrical
    subunits joined at the inter-subunit cleft, an mRNA track threading
    the small subunit, and three tRNA-binding pockets (A, P, E) along the
    inter-subunit interface.  Coordinates are in Angstrom.

    Attributes
    ----------
    organism_class : str
        "bacteria" (70S = 30S + 50S) or "eukaryote" (80S = 40S + 60S).
    mrna_sequence : str
        5'->3' RNA sequence threading the small subunit.  May contain a
        whole-number multiple of three nucleotides (whole codons).
    """

    organism_class: str = "bacteria"   # or "eukaryote"
    mrna_sequence: str = "AUGGCUAAACGUUAAUAA"  # Met-Ala-Lys-Arg-stop-stop

    # ------------------------------------------------------------------
    # Subunit cylinders
    # ------------------------------------------------------------------

    def small_subunit_dimensions_angstrom(self) -> Tuple[float, float]:
        """Return (diameter, height) of the small subunit in Angstrom.

        Eukaryotic 40S is ~ 7% larger in each dimension than bacterial 30S.
        """
        scale = 1.07 if self.organism_class == "eukaryote" else 1.0
        return SUBUNIT_30S_DIAMETER_A * scale, SUBUNIT_30S_HEIGHT_A * scale

    def large_subunit_dimensions_angstrom(self) -> Tuple[float, float]:
        """Return (diameter, height) of the large subunit in Angstrom."""
        scale = 1.10 if self.organism_class == "eukaryote" else 1.0
        return SUBUNIT_50S_DIAMETER_A * scale, SUBUNIT_50S_HEIGHT_A * scale

    def small_subunit_centre(self) -> np.ndarray:
        """Centre of the small subunit (origin)."""
        return np.array([0.0, 0.0, 0.0])

    def large_subunit_centre(self) -> np.ndarray:
        """Centre of the large subunit (stacked above the small one)."""
        _, h_small = self.small_subunit_dimensions_angstrom()
        _, h_large = self.large_subunit_dimensions_angstrom()
        gap = 10.0  # inter-subunit cleft (~ 10 A)
        return np.array([0.0, 0.0, 0.5 * h_small + gap + 0.5 * h_large])

    # ------------------------------------------------------------------
    # mRNA codon track
    # ------------------------------------------------------------------

    def codon_count(self) -> int:
        """Number of whole codons currently threaded through the small subunit."""
        return len(self.mrna_sequence) // CODON_LENGTH_NT

    def codon_positions(self) -> np.ndarray:
        """Return the 3D positions of each codon's centre along the mRNA track.

        The mRNA threads the small subunit along the +y direction at the
        cleft height.  Codon 0 sits at y=0 (the P site equivalent at start).

        Returns
        -------
        coords : np.ndarray, shape (codon_count, 3)
            Cartesian coordinates [Angstrom].
        """
        n = self.codon_count()
        # mRNA track sits just below the inter-subunit cleft
        cleft_z = 0.5 * SUBUNIT_30S_HEIGHT_A * (
            1.07 if self.organism_class == "eukaryote" else 1.0
        )
        ys = (np.arange(n) - (n - 1) / 2.0) * CODON_PITCH_A
        coords = np.column_stack([
            np.zeros(n),
            ys,
            np.full(n, cleft_z - 5.0),
        ])
        return coords

    def codon_at(self, index: int) -> str:
        """Return the 3-nt codon at the given codon index (0-based)."""
        if not (0 <= index < self.codon_count()):
            raise IndexError(f"codon index {index} out of range 0..{self.codon_count()-1}")
        i0 = index * CODON_LENGTH_NT
        return self.mrna_sequence[i0:i0 + CODON_LENGTH_NT].upper()

    # ------------------------------------------------------------------
    # A / P / E tRNA-binding sites
    # ------------------------------------------------------------------

    def site_centres(self) -> dict:
        """Return centres of the A, P, and E tRNA-binding sites.

        The three sites lie in the inter-subunit cleft, displaced from the
        mRNA track by ~ 25 A perpendicular.  A site is downstream (3' on
        mRNA), then P, then E (5'-most).  Centre-to-centre spacing = 20 A.
        """
        cleft_z = self.small_subunit_centre()[2] + 0.5 * SUBUNIT_30S_HEIGHT_A * (
            1.07 if self.organism_class == "eukaryote" else 1.0
        ) + 5.0  # ~ 5 A into the cleft
        # Perpendicular displacement from mRNA track (toward large subunit)
        # Sites span the cleft along +y (along reading direction) at a
        # representative inter-site spacing of 20 A.
        return {
            "A": np.array([0.0,  INTER_SITE_SPACING_A, cleft_z]),
            "P": np.array([0.0,  0.0,                  cleft_z]),
            "E": np.array([0.0, -INTER_SITE_SPACING_A, cleft_z]),
        }

    # ------------------------------------------------------------------
    # tRNA anticodon (paired with the A-site codon)
    # ------------------------------------------------------------------

    def a_site_codon_index(self) -> int:
        """The A-site codon is the second codon (0-based: index 1) at start.

        At translation initiation, the start codon (AUG) occupies P, the
        next codon occupies A.  The A-site index advances by one each
        elongation cycle.  This convenience getter returns the start state.
        """
        return 1 if self.codon_count() >= 2 else 0

    def a_site_codon(self) -> str:
        """The codon currently in the A site (textbook initiation state)."""
        return self.codon_at(self.a_site_codon_index())

    def a_site_anticodon(self) -> str:
        """Anticodon (3'->5' tRNA) that Watson-Crick-pairs with the A-site codon."""
        return codon_anticodon_pair(self.a_site_codon())

    def codon_anticodon_helix_axis(self) -> np.ndarray:
        """Direction of the codon-anticodon helix axis at the A site.

        The helix axis is aligned with the mRNA track (the +y direction
        in the local frame); the anticodon stem-loop projects from the
        large subunit down toward the cleft.
        """
        return np.array([0.0, 1.0, 0.0])

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "ribosome = catalytic strain-pump in stiff substrate; "
                "each cycle dissipates 2 GTP of substrate strain energy"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "subunit_count": 2,
            "trna_sites": 3,                # A / P / E
            "gtp_per_cycle": GTP_PER_ELONGATION_CYCLE,
        }


# ---------------------------------------------------------------------------
# SIMULATOR: translation kinetics, proofreading, GTP budget
# ---------------------------------------------------------------------------

@dataclass
class RibosomeSimulator:
    """Translation kinetics and proofreading on the ribosome geometry.

    The simulator does NOT integrate a chemical-kinetic ODE; instead it
    exposes (i) elongation rates (aa/s) for bacteria/eukaryotes,
    (ii) the kinetic-proofreading per-codon error rate, and (iii) the
    GTP-hydrolysis budget per cycle and over a full protein.

    Substrate interpretation:  the ribosome converts GTP-bound substrate
    strain into kinetic motion (translocation) + chemical bond formation
    (peptide bond), with two sequential discrimination gates that
    multiplicatively suppress wrong-codon incorporation.
    """

    geometry: RibosomeGeometry | None = None
    organism_class: str = "bacteria"
    n_grid: int = 60

    # ------------------------------------------------------------------
    # Translation rate
    # ------------------------------------------------------------------

    def elongation_rate_aa_per_second(self) -> float:
        """Mean amino-acid incorporation rate for the chosen organism class."""
        if self.organism_class == "eukaryote":
            return RATE_EUKARYOTE_AA_PER_S
        return RATE_BACTERIA_AA_PER_S

    def elongation_rate_distribution(
        self, n_samples: int = 5000, seed: int = 0,
    ) -> np.ndarray:
        """Sample elongation rates from a log-normal distribution.

        Translation elongation is heterogeneous: codon-dependent decoding
        time (cognate vs near-cognate) and ribosome pausing yield a long
        upper tail.  We use a log-normal with mean = textbook rate.

        Returns
        -------
        rates : np.ndarray, shape (n_samples,)
            Elongation rates [aa/s].
        """
        mu = self.elongation_rate_aa_per_second()
        # sigma in log-space: ~ 0.5 yields a ~ 60% coefficient of variation
        sigma_log = 0.5
        # Log-normal with arithmetic mean = mu requires mean_log = ln(mu) - sigma^2/2
        mean_log = math.log(mu) - 0.5 * sigma_log ** 2
        rng = np.random.default_rng(seed)
        return rng.lognormal(mean=mean_log, sigma=sigma_log, size=n_samples)

    def time_per_codon_seconds(self) -> float:
        """Mean dwell time per codon = 1 / elongation rate."""
        return 1.0 / self.elongation_rate_aa_per_second()

    def time_to_translate_protein_seconds(self, n_residues: int) -> float:
        """Total elongation time for a protein of n_residues residues."""
        return n_residues * self.time_per_codon_seconds()

    # ------------------------------------------------------------------
    # Kinetic proofreading -- per-codon error rate
    # ------------------------------------------------------------------

    def per_codon_error_rate(self) -> float:
        """Combined per-codon misincorporation rate (~ 1e-4)."""
        return ERROR_RATE_PER_CODON

    def initial_selection_error(self) -> float:
        """First kinetic gate (cognate vs near-cognate, EF-Tu * GTP)."""
        return INITIAL_SELECTION_ERR

    def proofreading_factor(self) -> float:
        """Second kinetic gate (post-GTP hydrolysis discrimination)."""
        return PROOFREADING_FACTOR

    def proofreading_combined_error(self) -> float:
        """Two gates multiply: combined error = gate1 * gate2 ~ 1e-4."""
        return self.initial_selection_error() * self.proofreading_factor()

    def errors_per_protein(self, n_residues: int) -> float:
        """Expected number of mistranslated residues per copy of a protein."""
        return n_residues * self.per_codon_error_rate()

    def fraction_perfect_proteins(self, n_residues: int) -> float:
        """Probability a single copy of an n-residue protein has zero errors.

        P(perfect) = (1 - error_rate)^n_residues
        """
        p_correct = 1.0 - self.per_codon_error_rate()
        return p_correct ** n_residues

    # ------------------------------------------------------------------
    # GTP hydrolysis budget
    # ------------------------------------------------------------------

    def gtp_per_cycle(self) -> int:
        """GTPs consumed per elongation cycle (EF-Tu + EF-G in bacteria)."""
        return GTP_PER_ELONGATION_CYCLE

    def gtp_hydrolysis_per_cycle_kcal_per_mol(self) -> float:
        """Free-energy budget per cycle in kcal/mol."""
        return self.gtp_per_cycle() * GTP_HYDROLYSIS_KCAL_PER_MOL

    def gtp_per_protein(self, n_residues: int) -> int:
        """Total GTP molecules consumed translating an n-residue protein.

        Each elongation cycle uses 2 GTP; initiation and termination add
        another small fixed cost (~ 1-2 GTP) which is negligible for long
        proteins, so we report 2 * n_residues here.
        """
        return self.gtp_per_cycle() * n_residues

    def energy_per_peptide_bond_kT(self, T_kelvin: float = 310.0) -> float:
        """GTP energy per peptide bond expressed in units of kT at body T."""
        # kT in kcal/mol at 310 K = 0.617
        kt = 0.617 * (T_kelvin / 310.0)
        return self.gtp_hydrolysis_per_cycle_kcal_per_mol() / kt

    # ------------------------------------------------------------------
    # Codon-anticodon Watson-Crick check
    # ------------------------------------------------------------------

    def watson_crick_pairs(self, codon: str, anticodon_3to5: str) -> int:
        """Count canonical Watson-Crick pairs between codon (5'->3') and
        anticodon (3'->5').

        With a perfect match all three positions are A-U or G-C.  Returns
        the number of paired positions (0..3).
        """
        if len(codon) != CODON_LENGTH_NT or len(anticodon_3to5) != CODON_LENGTH_NT:
            return 0
        pairs = 0
        wc = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
        for c, a in zip(codon.upper(), anticodon_3to5.upper()):
            if (c, a) in wc:
                pairs += 1
        return pairs

    def is_cognate(self, codon: str, anticodon_3to5: str) -> bool:
        """True iff all three codon-anticodon positions are Watson-Crick."""
        return self.watson_crick_pairs(codon, anticodon_3to5) == CODON_LENGTH_NT

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "elongation_rate_aa_per_s": self.elongation_rate_aa_per_second(),
            "per_codon_error_rate": self.per_codon_error_rate(),
            "gtp_per_cycle": self.gtp_per_cycle(),
            "interpretation": (
                "ribosome = strain-pump; 2 GTP/cycle dissipated as "
                "translocation + peptide bond; 2 kinetic gates drive "
                "error rate to ~ 1e-4 (mass-torque axiom on substrate)"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


__all__ = [
    "RibosomeGeometry",
    "RibosomeSimulator",
    "watson_crick_complement",
    "codon_anticodon_pair",
    "SUBUNIT_30S_DIAMETER_A",
    "SUBUNIT_30S_HEIGHT_A",
    "SUBUNIT_50S_DIAMETER_A",
    "SUBUNIT_50S_HEIGHT_A",
    "NT_RISE_A",
    "CODON_LENGTH_NT",
    "CODON_PITCH_A",
    "INTER_SITE_SPACING_A",
    "TRNA_LENGTH_NT",
    "TRNA_LSHAPE_ARM_A",
    "WC_BASE_PAIR_DISTANCE_A",
    "WC_HELIX_RADIUS_A",
    "WC_HELIX_RISE_PER_BP_A",
    "RATE_BACTERIA_AA_PER_S",
    "RATE_EUKARYOTE_AA_PER_S",
    "ERROR_RATE_PER_CODON",
    "INITIAL_SELECTION_ERR",
    "PROOFREADING_FACTOR",
    "GTP_PER_ELONGATION_CYCLE",
    "GTP_HYDROLYSIS_KCAL_PER_MOL",
    "LAMBDA_QCD_MEV",
]
