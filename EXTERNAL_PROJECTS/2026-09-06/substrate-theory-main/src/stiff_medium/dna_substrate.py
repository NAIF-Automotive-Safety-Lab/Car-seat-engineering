"""Paired DNA mechanics GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for DNA mechanics
under the B3 / stiff-medium ontology.

GEOMETRY
--------
``DNAGeometry`` represents the canonical B-form DNA double helix as a
right-handed, antiparallel pair of phosphate-sugar backbones threaded
through a stack of Watson-Crick base pairs.  The textbook B-DNA
parameters are wired in as constants:

    * rise per base pair:      3.4 A
    * base pairs per turn:     10.5  (so pitch = 35.7 A)
    * helix radius:            10.0 A
    * helix groove asymmetry:  major (~22 A wide) / minor (~12 A wide)
    * twist per base pair:     360 / 10.5 = 34.286 deg

Each base pair is realised as two complementary nucleotides (A:T or G:C)
sitting on opposite strands at angular offset 180 deg minus a small
groove offset.  The two backbones are antiparallel: strand 1 runs 5'->3'
in +z, strand 2 runs 3'->5' in +z.

SIMULATOR
---------
``DNASimulator`` exposes the polymer-mechanics observables of B-DNA:

    * Watson-Crick H-bond energies:
        - A:T pair has 2 H-bonds, ~ 7 kcal/mol total
        - G:C pair has 3 H-bonds, ~ 10 kcal/mol total
        => GC/AT binding ratio ~ 1.43 (cf. ~3/2 from H-bond counting)

    * Worm-like-chain elasticity:
        - persistence length L_p ~ 50 nm = 150 bp
        - stretching modulus  S  ~ 1100 pN
        - twist modulus       C  ~ 410 pN nm^2
        - twist-bend coupling g  ~ 30 pN nm   (Marko-Siggia 1995)

    * Force-extension F(x):
        - low force  (< 6 pN):    Marko-Siggia interpolation
                                  F = (k_B T / L_p) [0.25 (1 - x/L)^-2 - 0.25 + x/L]
        - mid force  (6-65 pN):   linear stretching with modulus S
        - B->S transition at F* ~ 65 pN: extension jumps from 1.0 L_0 to 1.7 L_0
        - high force (> 65 pN):   overstretched S-DNA plateau

    * Persistence length:
        L_p = (kappa / k_B T) where kappa is the bending modulus.
        For B-DNA, kappa ~ 50 nm * k_B T at 310 K  =>  L_p = 50 nm.

Substrate-specific predictions encoded here:
    * DNA double helix = paired stiff-medium strain pattern with
      antiparallel chirality, locked together by Watson-Crick H-bond
      cohesion (mass-torque axiom: the H-bond cohesion is the cumulative
      torque holding the pair against thermal bending).
    * The B->S overstretching transition is a substrate-elasticity
      phase transition (the helix releases stored twist into stretch).
    * Persistence length 50 nm = 150 bp is a derived ratio of bending
      modulus to thermal energy; it falls out of the WLC at body T.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) B-DNA values
# ---------------------------------------------------------------------------

# B-DNA helical parameters (Watson-Crick 1953; later Drew/Dickerson refinement)
RISE_PER_BP_A:        float = 3.4    # A per base pair (axial advance)
BP_PER_TURN:          float = 10.5   # base pairs per helical turn
HELIX_PITCH_A:        float = RISE_PER_BP_A * BP_PER_TURN  # 35.7 A
HELIX_RADIUS_A:       float = 10.0   # A (phosphate backbone radius)
TWIST_PER_BP_DEG:     float = 360.0 / BP_PER_TURN          # ~34.286 deg
MAJOR_GROOVE_WIDTH_A: float = 22.0
MINOR_GROOVE_WIDTH_A: float = 12.0

# Watson-Crick base-pair H-bond energies (kcal/mol, total per base pair)
# A:T has 2 H-bonds, G:C has 3 H-bonds; literature consensus values.
AT_HBOND_KCAL_PER_MOL: float = 7.0   # 2 H-bonds (~ 3.5 kcal/mol each)
GC_HBOND_KCAL_PER_MOL: float = 10.0  # 3 H-bonds (~ 3.3 kcal/mol each)
N_HBONDS_AT:           int   = 2
N_HBONDS_GC:           int   = 3

# Worm-like-chain elasticity parameters at body temperature
PERSISTENCE_LENGTH_NM:    float = 50.0    # nm  (~150 bp, textbook B-DNA)
STRETCH_MODULUS_PN:       float = 1100.0  # pN  (Smith et al. 1996)
TWIST_MODULUS_PN_NM2:     float = 410.0   # pN nm^2 (Bryant et al. 2003)
TWIST_BEND_COUPLING_PN_NM: float = 30.0   # pN nm  (Marko 1997)

# Overstretching (B->S) transition
B_TO_S_FORCE_PN:        float = 65.0   # pN  (Cluzel 1996, Smith 1996)
B_TO_S_EXTENSION_RATIO: float = 1.7    # contour length ratio after transition

# Thermodynamic anchor: kT at body temperature
TEMPERATURE_K:          float = 310.0
KB_J_PER_K:             float = 1.380649e-23
KB_T_BODY_PN_NM:        float = KB_J_PER_K * TEMPERATURE_K * 1.0e21  # pN nm
KT_AT_BODY_T_KCAL_PER_MOL: float = 0.617  # kT at 310 K in kcal/mol

# Substrate anchor (B3 framework)
LAMBDA_QCD_MEV: float = 200.0


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class DNAGeometry:
    """Right-handed, antiparallel B-DNA double-helix geometry.

    The helix is generated parametrically from the canonical B-form
    constants: at base-pair index i (0 .. n_bp - 1),
        z_i      = i * rise
        theta_i  = i * twist
    Strand 1 phosphate sits at angle theta_i on a cylinder of radius R;
    strand 2 phosphate sits at theta_i + pi - groove_offset (antiparallel
    + groove asymmetry).

    Attributes
    ----------
    n_bp : int
        Number of base pairs in the modelled helix.
    sequence : str | None
        ASCII sequence of strand 1 (5'->3'); strand 2 is its Watson-Crick
        complement.  If None, a default repeating "ATCG" sequence is used.
    """

    n_bp: int = 30
    sequence: str | None = None
    groove_offset_rad: float = math.radians(150.0)  # angular offset to strand 2
    rise_A: float = RISE_PER_BP_A
    bp_per_turn: float = BP_PER_TURN
    radius_A: float = HELIX_RADIUS_A

    def __post_init__(self) -> None:
        if self.sequence is None:
            self.sequence = ("ATCG" * (self.n_bp // 4 + 1))[: self.n_bp]
        elif len(self.sequence) != self.n_bp:
            raise ValueError(
                f"sequence length {len(self.sequence)} != n_bp {self.n_bp}"
            )

    # ------------------------------------------------------------------
    # Helical parameters
    # ------------------------------------------------------------------

    @property
    def twist_per_bp_deg(self) -> float:
        return 360.0 / self.bp_per_turn

    @property
    def pitch_A(self) -> float:
        """Helical pitch = rise per turn (axial distance per 360 deg)."""
        return self.rise_A * self.bp_per_turn

    @property
    def contour_length_A(self) -> float:
        """Total contour length along the helix axis."""
        return self.n_bp * self.rise_A

    # ------------------------------------------------------------------
    # Base-pair complementarity (Watson-Crick)
    # ------------------------------------------------------------------

    @staticmethod
    def complement(base: str) -> str:
        """Return Watson-Crick complement of A/T/G/C (case insensitive)."""
        wc = {"A": "T", "T": "A", "G": "C", "C": "G",
              "a": "t", "t": "a", "g": "c", "c": "g"}
        if base not in wc:
            raise ValueError(f"unknown base {base!r}")
        return wc[base]

    def complement_sequence(self) -> str:
        """Return the antiparallel complement strand 2 (5'->3')."""
        # strand 2 is reverse complement of strand 1
        return "".join(self.complement(b) for b in reversed(self.sequence))

    @staticmethod
    def is_at_pair(base: str) -> bool:
        return base.upper() in ("A", "T")

    @staticmethod
    def is_gc_pair(base: str) -> bool:
        return base.upper() in ("G", "C")

    # ------------------------------------------------------------------
    # Backbone phosphate coordinates
    # ------------------------------------------------------------------

    def strand1_coords(self) -> np.ndarray:
        """Phosphate coordinates of strand 1 (5' -> 3' along +z)."""
        i = np.arange(self.n_bp)
        theta = np.radians(i * self.twist_per_bp_deg)
        z = i * self.rise_A
        x = self.radius_A * np.cos(theta)
        y = self.radius_A * np.sin(theta)
        return np.column_stack([x, y, z])

    def strand2_coords(self) -> np.ndarray:
        """Phosphate coordinates of strand 2 (3' -> 5' along +z, antiparallel).

        Strand 2 sits on the other side of the helix axis at angle
        theta + pi - groove_offset, producing the major/minor groove
        asymmetry of B-DNA.
        """
        i = np.arange(self.n_bp)
        theta = np.radians(i * self.twist_per_bp_deg) + math.pi - self.groove_offset_rad
        z = i * self.rise_A
        x = self.radius_A * np.cos(theta)
        y = self.radius_A * np.sin(theta)
        return np.column_stack([x, y, z])

    # ------------------------------------------------------------------
    # Base-pair (axis) coordinates
    # ------------------------------------------------------------------

    def base_pair_axis_coords(self) -> np.ndarray:
        """Helical-axis position of each base pair (centre of the rung)."""
        i = np.arange(self.n_bp)
        z = i * self.rise_A
        return np.column_stack([np.zeros_like(z), np.zeros_like(z), z])

    def base_pair_endpoints(self) -> Tuple[np.ndarray, np.ndarray]:
        """Endpoints of each base-pair rung on strand 1 and strand 2."""
        return self.strand1_coords(), self.strand2_coords()

    # ------------------------------------------------------------------
    # Groove geometry
    # ------------------------------------------------------------------

    @property
    def major_groove_width_A(self) -> float:
        return MAJOR_GROOVE_WIDTH_A

    @property
    def minor_groove_width_A(self) -> float:
        return MINOR_GROOVE_WIDTH_A

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "DNA double helix = antiparallel paired stiff-medium "
                "strain pattern, locked by Watson-Crick H-bond cohesion"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "rise_per_bp_A": self.rise_A,
            "bp_per_turn": self.bp_per_turn,
            "helix_radius_A": self.radius_A,
            "twist_per_bp_deg": self.twist_per_bp_deg,
            "right_handed": True,
            "antiparallel": True,
        }


# ---------------------------------------------------------------------------
# SIMULATOR: WLC mechanics, base-pair energetics, force-extension
# ---------------------------------------------------------------------------

@dataclass
class DNASimulator:
    """B-DNA polymer mechanics: H-bond energetics and WLC elasticity.

    Substrate interpretation
    ------------------------
    * Watson-Crick H-bonds = directional substrate-strain cohesion at
      paired sites; the AT vs GC ratio is purely geometric (2 vs 3
      bond sites per pair).
    * Persistence length, stretch modulus, and twist modulus are the
      three independent elastic constants of the helix as a substrate
      filament.
    * The B->S transition is a substrate-elasticity phase change: the
      helix loses its twist register and stretches by ~ 70 %.
    """

    geometry: DNAGeometry | None = None
    persistence_length_nm: float = PERSISTENCE_LENGTH_NM
    stretch_modulus_pn:    float = STRETCH_MODULUS_PN
    twist_modulus_pn_nm2:  float = TWIST_MODULUS_PN_NM2
    twist_bend_coupling_pn_nm: float = TWIST_BEND_COUPLING_PN_NM
    b_to_s_force_pn:       float = B_TO_S_FORCE_PN
    b_to_s_extension_ratio: float = B_TO_S_EXTENSION_RATIO
    temperature_K:         float = TEMPERATURE_K

    # ------------------------------------------------------------------
    # Watson-Crick base-pair energetics
    # ------------------------------------------------------------------

    @staticmethod
    def at_binding_kcal_per_mol() -> float:
        """Total A:T binding energy (2 H-bonds)."""
        return AT_HBOND_KCAL_PER_MOL

    @staticmethod
    def gc_binding_kcal_per_mol() -> float:
        """Total G:C binding energy (3 H-bonds)."""
        return GC_HBOND_KCAL_PER_MOL

    @staticmethod
    def base_pair_energy_kcal_per_mol(base: str) -> float:
        """Return Watson-Crick total H-bond energy of a base pair."""
        b = base.upper()
        if b in ("A", "T"):
            return AT_HBOND_KCAL_PER_MOL
        if b in ("G", "C"):
            return GC_HBOND_KCAL_PER_MOL
        raise ValueError(f"unknown base {base!r}")

    @classmethod
    def gc_to_at_energy_ratio(cls) -> float:
        """G:C / A:T total H-bond energy ratio (target ~ 1.43)."""
        return cls.gc_binding_kcal_per_mol() / cls.at_binding_kcal_per_mol()

    @classmethod
    def gc_to_at_hbond_count_ratio(cls) -> float:
        """G:C / A:T H-bond count ratio = 3/2 = 1.5."""
        return N_HBONDS_GC / N_HBONDS_AT

    def total_binding_energy_kcal_per_mol(self, sequence: str) -> float:
        """Sum of Watson-Crick binding energies along a strand."""
        return sum(self.base_pair_energy_kcal_per_mol(b) for b in sequence)

    # ------------------------------------------------------------------
    # WLC persistence length (in bp and nm)
    # ------------------------------------------------------------------

    def persistence_length_bp(self) -> float:
        """Persistence length expressed in base pairs (~ 150 bp)."""
        rise_nm = RISE_PER_BP_A * 0.1
        return self.persistence_length_nm / rise_nm

    def bending_modulus_pn_nm2(self) -> float:
        """Bending modulus kappa = L_p * k_B T (pN nm^2)."""
        return self.persistence_length_nm * KB_T_BODY_PN_NM

    # ------------------------------------------------------------------
    # Marko-Siggia worm-like-chain force-extension
    # ------------------------------------------------------------------

    def marko_siggia_force_pn(self, x_over_L: np.ndarray) -> np.ndarray:
        """Marko-Siggia interpolation (Macromolecules 1995).

            F * L_p / (k_B T) = 0.25 (1 - x/L)^-2 - 0.25 + x/L

        Valid for entropic regime x/L < ~0.97; diverges as x -> L.
        """
        x = np.clip(np.asarray(x_over_L, dtype=float), 0.0, 0.999)
        f_dimensionless = 0.25 / (1.0 - x) ** 2 - 0.25 + x
        # F = (k_B T / L_p) * f_dimensionless ; convert L_p nm -> nm to get pN
        f_pn = (KB_T_BODY_PN_NM / self.persistence_length_nm) * f_dimensionless
        return f_pn

    def stretching_force_pn(self, x_over_L: np.ndarray) -> np.ndarray:
        """Linear elastic stretching beyond the entropic regime.

        F_lin(x/L) = S * (x/L - 1)   for x/L > 1
        """
        x = np.asarray(x_over_L, dtype=float)
        return self.stretch_modulus_pn * (x - 1.0)

    def force_extension_curve(
        self, x_over_L: np.ndarray | None = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Composite force-extension F(x/L) including B->S transition.

        Three regimes are stitched together:
            (i)   entropic Marko-Siggia for x/L in [0, ~0.97]
            (ii)  linear elastic stretch (modulus S) for x/L in [1.0, 1.05]
            (iii) B->S overstretching plateau at F* until x/L = 1.7
            (iv)  S-DNA elastic regime beyond x/L = 1.7
        """
        if x_over_L is None:
            x_over_L = np.linspace(0.001, 1.85, 400)
        x = np.asarray(x_over_L, dtype=float)

        f = np.zeros_like(x)

        # Regime I: WLC entropic (x/L <= 1)
        wlc_mask = x <= 1.0
        f[wlc_mask] = self.marko_siggia_force_pn(x[wlc_mask])

        # Regime II/III: B->S plateau on [1.0, b_to_s_extension_ratio]
        plateau_mask = (x > 1.0) & (x <= self.b_to_s_extension_ratio)
        # Smooth ramp to the plateau over a narrow band
        ramp_end = 1.005
        ramp_mask = plateau_mask & (x < ramp_end)
        # WLC extrapolated value at x=1 should match plateau force smoothly
        f_at_1 = float(self.marko_siggia_force_pn(np.array([0.999]))[0])
        # Linear interpolate from F_at_1 to b_to_s_force_pn over [1.0, ramp_end]
        f[ramp_mask] = (
            f_at_1
            + (self.b_to_s_force_pn - f_at_1)
            * (x[ramp_mask] - 1.0) / (ramp_end - 1.0)
        )
        post_ramp_mask = plateau_mask & (x >= ramp_end)
        f[post_ramp_mask] = self.b_to_s_force_pn

        # Regime IV: post-overstretch S-DNA, linear with modulus S
        post_mask = x > self.b_to_s_extension_ratio
        f[post_mask] = (
            self.b_to_s_force_pn
            + self.stretch_modulus_pn * (x[post_mask] - self.b_to_s_extension_ratio)
        )

        return x, f

    # ------------------------------------------------------------------
    # Twist-bend coupling
    # ------------------------------------------------------------------

    def twist_bend_coupling_energy(
        self, bend_curvature_per_nm: float, twist_per_nm: float,
    ) -> float:
        """Marko-Siggia 1994 twist-bend coupling cross term, units pN nm.

            E_g = g * (bend_curvature) * (excess_twist)

        With g = TWIST_BEND_COUPLING_PN_NM ~ 30 pN nm.
        """
        return self.twist_bend_coupling_pn_nm * bend_curvature_per_nm * twist_per_nm

    def twist_modulus_kT_per_bp(self) -> float:
        """Twist modulus expressed in kT per bp (dimensionless check)."""
        rise_nm = RISE_PER_BP_A * 0.1
        return self.twist_modulus_pn_nm2 * rise_nm / KB_T_BODY_PN_NM

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "persistence_length_nm": self.persistence_length_nm,
            "persistence_length_bp": self.persistence_length_bp(),
            "stretch_modulus_pn":    self.stretch_modulus_pn,
            "twist_modulus_pn_nm2":  self.twist_modulus_pn_nm2,
            "twist_bend_coupling_pn_nm": self.twist_bend_coupling_pn_nm,
            "b_to_s_force_pn":       self.b_to_s_force_pn,
            "b_to_s_extension_ratio": self.b_to_s_extension_ratio,
            "lambda_qcd_MeV":        LAMBDA_QCD_MEV,
            "interpretation": (
                "DNA mechanics = elasticity of paired antiparallel "
                "stiff-medium filament; B->S = substrate phase change"
            ),
        }


__all__ = [
    "DNAGeometry",
    "DNASimulator",
    "RISE_PER_BP_A",
    "BP_PER_TURN",
    "HELIX_PITCH_A",
    "HELIX_RADIUS_A",
    "TWIST_PER_BP_DEG",
    "MAJOR_GROOVE_WIDTH_A",
    "MINOR_GROOVE_WIDTH_A",
    "AT_HBOND_KCAL_PER_MOL",
    "GC_HBOND_KCAL_PER_MOL",
    "N_HBONDS_AT",
    "N_HBONDS_GC",
    "PERSISTENCE_LENGTH_NM",
    "STRETCH_MODULUS_PN",
    "TWIST_MODULUS_PN_NM2",
    "TWIST_BEND_COUPLING_PN_NM",
    "B_TO_S_FORCE_PN",
    "B_TO_S_EXTENSION_RATIO",
    "TEMPERATURE_K",
    "KB_T_BODY_PN_NM",
    "KT_AT_BODY_T_KCAL_PER_MOL",
    "LAMBDA_QCD_MEV",
]
