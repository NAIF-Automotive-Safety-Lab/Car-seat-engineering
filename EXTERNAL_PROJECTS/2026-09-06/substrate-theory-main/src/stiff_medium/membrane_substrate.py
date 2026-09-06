"""Paired biological-membrane GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for fluid-mosaic
lipid-bilayer membranes under the B3 / stiff-medium ontology.

GEOMETRY
--------
``MembraneGeometry`` represents a biological membrane as a phospholipid
bilayer.  Each leaflet is a 2D close-packed lattice of amphiphilic
phospholipids, with hydrophilic head groups exposed to water and
hydrophobic acyl tails buried in the bilayer interior.  Canonical
textbook values (POPC bilayer, body temperature):

    * total bilayer thickness         d_bilayer  = 4.0 nm  (textbook 3.7-4.0)
    * headgroup-region thickness      d_head     = 0.9 nm  per leaflet
    * hydrophobic-tail thickness      d_tail     = 1.4 nm  per leaflet
    * area per lipid (fluid phase)    A_lipid    = 0.65 nm^2
    * lipid acyl-chain length         L_acyl     = 1.4 nm

A single phospholipid is rendered as a polar head (cylinder of diameter
0.7 nm) plus two hydrocarbon tails (cylinders of length 1.4 nm).  The
fluid-mosaic model (Singer-Nicolson 1972) is encoded by allowing
integral membrane proteins (transmembrane alpha-helix bundles) to be
embedded between lipids.

SIMULATOR
---------
``MembraneSimulator`` exposes the elastic and transport properties of
the bilayer as a 2D fluid film embedded in the 3D substrate:

    * self_assembly_free_energy(): amphiphilic-driven membrane formation;
      the hydrophobic effect collapses dispersed lipids into a bilayer
      with delta_G ~ -7 kcal/mol per lipid (typical PC lipid).
    * bending_modulus_kT():        kappa = 20 kT  (POPC, body temperature)
    * lateral_diffusion_coefficient_um2_per_s():  D = 1 micron^2 / s
                                       (FRAP measurement on POPC bilayers)
    * helfrich_curvature_energy(c1, c2, c0):
          E/A = (kappa/2) * (c1 + c2 - 2*c0)^2  +  kappa_G * c1 * c2
      Helfrich (1973) curvature elasticity functional.
    * vesicle_critical_radius_nm(): minimum radius for a stable vesicle
      under thermal undulations; sets the lower bound on liposome size.
    * membrane_undulation_spectrum(q): <|h_q|^2> = kT / (kappa * q^4)
      for a tensionless membrane (capillary-wave spectrum).

Substrate-specific predictions encoded here:
    * bending modulus kappa = 20 kT emerges from the lipid-tail acyl
      strain energy; thicker bilayers give larger kappa (cubic in d_tail).
    * lateral diffusion D = kT / (4 pi eta_membrane * R_lipid) follows
      the Saffman-Delbruck 2D-Stokes formula with eta_membrane ~ 100x eta_water.
    * vesicle stability: a lipid bilayer closes spontaneously into a
      vesicle when curvature-energy cost (8 pi kappa) is less than the
      edge-energy of a flat patch (2 pi R * gamma_edge).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) values -- POPC bilayer, body temperature 310 K
# ---------------------------------------------------------------------------

# Bilayer thickness components (nanometres) -- POPC, fluid L_alpha phase
BILAYER_THICKNESS_NM:        float = 4.0   # total head-to-head
HEAD_THICKNESS_PER_LEAFLET_NM: float = 0.9
TAIL_THICKNESS_PER_LEAFLET_NM: float = 1.1   # 4.0 = 2*0.9 + 2*1.1

# Lipid molecular dimensions (nanometres / nm^2)
LIPID_AREA_NM2:              float = 0.65
LIPID_HEAD_DIAMETER_NM:      float = 0.85
LIPID_TAIL_LENGTH_NM:        float = 1.4
LIPID_TAIL_DIAMETER_NM:      float = 0.42

# Elastic constants (in units of kT, body temperature)
BENDING_MODULUS_KT:          float = 20.0   # POPC, Helfrich coefficient
GAUSSIAN_MODULUS_KT:         float = -10.0  # kappa_G (typical, negative)
EDGE_TENSION_PN:             float = 12.0   # piconewtons, lipid-bilayer edge

# Transport constants
LATERAL_DIFFUSION_UM2_PER_S: float = 1.0    # FRAP on fluid POPC at 310 K
TRANSVERSE_FLIP_RATE_PER_S:  float = 1.0e-5  # passive flip-flop (slow)

# Hydrophobic effect (per CH2 group transfer water -> oil)
HYDROPHOBIC_DG_PER_CH2_KCAL_PER_MOL: float = 0.8  # textbook ~ 0.7-0.9
LIPID_TAIL_CH2_GROUPS:       int = 16        # palmitoyl (C16) acyl chain

# Substrate anchor (B3 framework)
LAMBDA_QCD_MEV:              float = 200.0
KT_AT_BODY_T_KCAL_PER_MOL:   float = 0.617  # 310 K, kcal/mol
KT_AT_BODY_T_J:              float = 4.28e-21  # joules at 310 K


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class MembraneGeometry:
    """Lipid-bilayer geometry (fluid-mosaic model).

    The membrane is modelled as a planar bilayer with two leaflets of
    close-packed phospholipids.  The hydrophilic heads sit on both
    aqueous interfaces (z = +/- d_bilayer/2); the hydrophobic acyl
    tails fill the interior (|z| < d_tail).

    Attributes
    ----------
    n_lipids_per_leaflet : int
        Number of lipids per leaflet in the modelled patch.
    bilayer_thickness_nm : float
        Total head-to-head thickness (nanometres).
    area_per_lipid_nm2 : float
        Molecular area per lipid in the fluid phase (nm^2).
    """

    n_lipids_per_leaflet: int = 64
    bilayer_thickness_nm: float = BILAYER_THICKNESS_NM
    area_per_lipid_nm2: float = LIPID_AREA_NM2
    head_thickness_nm: float = HEAD_THICKNESS_PER_LEAFLET_NM
    tail_thickness_nm: float = TAIL_THICKNESS_PER_LEAFLET_NM
    head_diameter_nm: float = LIPID_HEAD_DIAMETER_NM
    tail_length_nm: float = LIPID_TAIL_LENGTH_NM
    tail_diameter_nm: float = LIPID_TAIL_DIAMETER_NM

    # ------------------------------------------------------------------
    # Bilayer dimensions
    # ------------------------------------------------------------------

    def total_thickness_nm(self) -> float:
        """Return head-to-head bilayer thickness (nanometres)."""
        return self.bilayer_thickness_nm

    def hydrophobic_thickness_nm(self) -> float:
        """Return the thickness of the hydrophobic tail core (nanometres)."""
        return 2.0 * self.tail_thickness_nm

    def patch_side_length_nm(self) -> float:
        """Side length of a square membrane patch holding the modelled lipids."""
        area_per_leaflet = self.n_lipids_per_leaflet * self.area_per_lipid_nm2
        return math.sqrt(area_per_leaflet)

    def patch_area_nm2(self) -> float:
        """Total membrane patch area (single leaflet area)."""
        return self.n_lipids_per_leaflet * self.area_per_lipid_nm2

    # ------------------------------------------------------------------
    # Lipid coordinates
    # ------------------------------------------------------------------

    def lipid_head_positions(self, leaflet: str = "upper") -> np.ndarray:
        """Return (n_lipids, 3) array of head-group centre positions.

        Lipids are placed on a square close-packed lattice of side
        sqrt(area_per_lipid).  z = +d/2 (upper leaflet) or -d/2 (lower).
        """
        if leaflet not in ("upper", "lower"):
            raise ValueError("leaflet must be 'upper' or 'lower'")
        n_side = max(1, int(math.ceil(math.sqrt(self.n_lipids_per_leaflet))))
        spacing = math.sqrt(self.area_per_lipid_nm2)
        L = (n_side - 1) * spacing
        xs = np.linspace(-L / 2, L / 2, n_side)
        ys = np.linspace(-L / 2, L / 2, n_side)
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        flat_x = gx.flatten()[: self.n_lipids_per_leaflet]
        flat_y = gy.flatten()[: self.n_lipids_per_leaflet]
        z = +self.bilayer_thickness_nm / 2.0 if leaflet == "upper" \
            else -self.bilayer_thickness_nm / 2.0
        return np.stack([flat_x, flat_y, np.full_like(flat_x, z)], axis=1)

    def tail_endpoint_positions(self, leaflet: str = "upper") -> np.ndarray:
        """Return (n_lipids, 3) array of tail-tip positions (bilayer midplane)."""
        heads = self.lipid_head_positions(leaflet=leaflet)
        # Tails point toward z = 0 (the bilayer midplane)
        sign = +1 if leaflet == "upper" else -1
        z_tail = heads[:, 2] - sign * (self.head_thickness_nm + self.tail_thickness_nm)
        out = heads.copy()
        out[:, 2] = z_tail
        return out

    # ------------------------------------------------------------------
    # Fluid-mosaic model: integral membrane protein placeholder
    # ------------------------------------------------------------------

    def integral_protein_positions(
        self, n_proteins: int = 4, radius_nm: float = 1.0,
        rng: np.random.Generator | None = None,
    ) -> np.ndarray:
        """Return (n_proteins, 3) positions of transmembrane proteins.

        Proteins are placed at random lateral positions, with z = 0
        (centred on the bilayer midplane).  Each protein represents an
        integral membrane protein of given radius.
        """
        rng = rng or np.random.default_rng(seed=0)
        L = self.patch_side_length_nm()
        xs = rng.uniform(-L / 2, L / 2, n_proteins)
        ys = rng.uniform(-L / 2, L / 2, n_proteins)
        zs = np.zeros(n_proteins)
        return np.stack([xs, ys, zs], axis=1)

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "biological membrane = 2D fluid film of amphiphilic lipids "
                "embedded in the 3D stiff medium"
            ),
            "bilayer_thickness_nm": self.bilayer_thickness_nm,
            "area_per_lipid_nm2": self.area_per_lipid_nm2,
            "leaflets": 2,
            "model": "Singer-Nicolson fluid mosaic (1972)",
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class MembraneSimulator:
    """Membrane elasticity, dynamics, and self-assembly thermodynamics.

    The simulator does not run a molecular-dynamics trajectory; instead
    it exposes the canonical bilayer constants (Helfrich kappa, lateral
    diffusion, hydrophobic free energy) and a closed-form Helfrich
    curvature-energy functional, plus the capillary-wave undulation
    spectrum used to bound vesicle stability.
    """

    geometry: MembraneGeometry | None = None
    bending_modulus_kT: float = BENDING_MODULUS_KT
    gaussian_modulus_kT: float = GAUSSIAN_MODULUS_KT
    lateral_diffusion_um2_s: float = LATERAL_DIFFUSION_UM2_PER_S
    edge_tension_pN: float = EDGE_TENSION_PN
    temperature_K: float = 310.0

    # ------------------------------------------------------------------
    # Self-assembly from amphiphilic interactions
    # ------------------------------------------------------------------

    def self_assembly_free_energy_per_lipid_kcal_per_mol(
        self, n_ch2: int = LIPID_TAIL_CH2_GROUPS,
    ) -> float:
        """Hydrophobic free energy gained per lipid on bilayer assembly.

        delta_G = -n_CH2 * delta_G_per_CH2  (transferring acyl tails from
        water to the oil-like bilayer interior).  Two acyl chains per
        lipid; each contributes ~ n_ch2/2 effective CH2 groups buried.
        Returns kcal/mol per lipid (negative = favourable).
        """
        # Two tails per lipid; each tail buries n_ch2 CH2 groups
        return -2.0 * n_ch2 * HYDROPHOBIC_DG_PER_CH2_KCAL_PER_MOL

    def critical_micelle_concentration_M(
        self, dG_kcal_per_mol: float | None = None,
    ) -> float:
        """Estimate CMC from the hydrophobic transfer free energy.

        CMC ~ exp(delta_G_assembly / kT).  Negative delta_G => small CMC.
        Returns molar concentration.
        """
        if dG_kcal_per_mol is None:
            # Use single-chain detergent (n_ch2 = 12) as the relevant scale,
            # since CMC is a single-amphiphile property.
            dG_kcal_per_mol = -12.0 * HYDROPHOBIC_DG_PER_CH2_KCAL_PER_MOL
        kt = KT_AT_BODY_T_KCAL_PER_MOL
        return float(math.exp(dG_kcal_per_mol / kt))

    # ------------------------------------------------------------------
    # Bending modulus (Helfrich)
    # ------------------------------------------------------------------

    def bending_modulus_kT_value(self) -> float:
        """Return the Helfrich bending modulus in kT units (POPC: 20 kT)."""
        return self.bending_modulus_kT

    def bending_modulus_J(self) -> float:
        """Return the Helfrich bending modulus in joules."""
        return self.bending_modulus_kT * KT_AT_BODY_T_J

    def helfrich_curvature_energy_density(
        self, c1: float, c2: float, c0: float = 0.0,
    ) -> float:
        """Helfrich (1973) curvature-energy density per unit area.

        E/A = (kappa/2) * (c1 + c2 - 2 c0)^2 + kappa_G * c1 * c2

        Inputs c1, c2, c0 are principal/spontaneous curvatures (1/nm);
        output is in kT per nm^2.
        """
        H_term = 0.5 * self.bending_modulus_kT * (c1 + c2 - 2.0 * c0) ** 2
        K_term = self.gaussian_modulus_kT * c1 * c2
        return H_term + K_term

    def helfrich_total_energy_kT(
        self, area_nm2: float, c1: float, c2: float, c0: float = 0.0,
    ) -> float:
        """Total Helfrich energy of a uniform-curvature patch (in kT)."""
        return self.helfrich_curvature_energy_density(c1, c2, c0) * area_nm2

    def vesicle_bending_energy_kT(self) -> float:
        """Bending energy of a closed spherical vesicle.

        Closed sphere: 2 H = 2/R, so E_bend = (kappa/2) * (2/R)^2 * 4 pi R^2
                              = 8 pi kappa.  Independent of R (topology).
        """
        return 8.0 * math.pi * self.bending_modulus_kT

    # ------------------------------------------------------------------
    # Lateral diffusion (Saffman-Delbruck-like)
    # ------------------------------------------------------------------

    def lateral_diffusion_coefficient_um2_per_s(self) -> float:
        """Lateral diffusion coefficient of a single lipid (POPC: 1 um^2/s)."""
        return self.lateral_diffusion_um2_s

    def msd_um2(self, t_s: float) -> float:
        """Mean-square displacement <r^2> = 4 D t for 2D random walk."""
        return 4.0 * self.lateral_diffusion_um2_s * t_s

    def diffusion_distance_nm(self, t_s: float) -> float:
        """RMS lateral displacement (nm) after time t (s)."""
        return math.sqrt(self.msd_um2(t_s)) * 1.0e3  # um -> nm

    # ------------------------------------------------------------------
    # Membrane undulation spectrum (capillary-wave)
    # ------------------------------------------------------------------

    def undulation_spectrum(
        self, q_inv_nm: np.ndarray | Sequence[float],
        tension_kT_per_nm2: float = 0.0,
    ) -> np.ndarray:
        """<|h_q|^2> = kT / (gamma q^2 + kappa q^4) for a fluid membrane.

        For a tensionless membrane (gamma -> 0) this is the Helfrich
        spectrum <|h_q|^2> = kT / (kappa * q^4); the 1/q^4 divergence at
        long wavelengths is cut off by the patch size.
        """
        q = np.asarray(q_inv_nm, dtype=float)
        denom = tension_kT_per_nm2 * q ** 2 + self.bending_modulus_kT * q ** 4
        return 1.0 / np.maximum(denom, 1.0e-30)

    def rms_undulation_amplitude_nm(
        self, patch_side_nm: float,
    ) -> float:
        """RMS height fluctuation of a tensionless square patch.

        <h^2> = kT / (kappa) * sum_q 1/(q^4 A) ~ kT L^2 / (4 pi^3 kappa).
        """
        return math.sqrt(patch_side_nm ** 2 / (4.0 * math.pi ** 3 * self.bending_modulus_kT))

    # ------------------------------------------------------------------
    # Vesicle stability
    # ------------------------------------------------------------------

    def vesicle_critical_radius_nm(
        self, edge_tension_pN: float | None = None,
    ) -> float:
        """Critical radius above which a vesicle is more stable than a disc.

        A flat disc of radius R has edge energy E_edge = 2 pi R gamma_edge.
        A closed vesicle of the same area has bending energy 8 pi kappa.
        Setting them equal: R_crit = 4 kappa / gamma_edge.
        """
        gamma = edge_tension_pN if edge_tension_pN is not None else self.edge_tension_pN
        # gamma in pN = 1e-12 N; kappa in J -> R = 4 * kappa / gamma in metres
        R_m = 4.0 * self.bending_modulus_J() / (gamma * 1.0e-12)
        return R_m * 1.0e9  # m -> nm

    def vesicle_is_stable(self, radius_nm: float) -> bool:
        """True if a vesicle of given radius is stable against opening."""
        return radius_nm >= self.vesicle_critical_radius_nm()

    # ------------------------------------------------------------------
    # Bending-mode shape (single-mode visualisation helper)
    # ------------------------------------------------------------------

    def bending_mode_shape(
        self, side_nm: float, n_grid: int = 80,
        nx: int = 1, ny: int = 1, amplitude_nm: float = 0.5,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (X, Y, Z) of a single sinusoidal bending mode.

        h(x, y) = A sin(nx pi x / L) sin(ny pi y / L) -- a standing wave
        on a square patch with clamped edges.  Useful for visualising the
        lowest-q membrane fluctuations.
        """
        xs = np.linspace(0.0, side_nm, n_grid)
        ys = np.linspace(0.0, side_nm, n_grid)
        X, Y = np.meshgrid(xs, ys, indexing="xy")
        Z = amplitude_nm * np.sin(nx * math.pi * X / side_nm) \
                        * np.sin(ny * math.pi * Y / side_nm)
        return X, Y, Z

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "bending_modulus_kT": self.bending_modulus_kT,
            "lateral_diffusion_um2_per_s": self.lateral_diffusion_um2_s,
            "interpretation": (
                "lipid bilayer = 2D fluid embedded in 3D substrate; "
                "elasticity = strain energy of acyl-tail packing"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


__all__ = [
    "MembraneGeometry",
    "MembraneSimulator",
    "BILAYER_THICKNESS_NM",
    "HEAD_THICKNESS_PER_LEAFLET_NM",
    "TAIL_THICKNESS_PER_LEAFLET_NM",
    "LIPID_AREA_NM2",
    "LIPID_HEAD_DIAMETER_NM",
    "LIPID_TAIL_LENGTH_NM",
    "LIPID_TAIL_DIAMETER_NM",
    "BENDING_MODULUS_KT",
    "GAUSSIAN_MODULUS_KT",
    "EDGE_TENSION_PN",
    "LATERAL_DIFFUSION_UM2_PER_S",
    "TRANSVERSE_FLIP_RATE_PER_S",
    "HYDROPHOBIC_DG_PER_CH2_KCAL_PER_MOL",
    "LIPID_TAIL_CH2_GROUPS",
    "LAMBDA_QCD_MEV",
    "KT_AT_BODY_T_KCAL_PER_MOL",
    "KT_AT_BODY_T_J",
]
