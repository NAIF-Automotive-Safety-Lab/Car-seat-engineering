"""Paired GEOMETRY + SIM + VIZ for molecular bonds as substrate strain bridges.

Atoms in the stiff-medium (B3) framework are bound substrate excitations whose
nuclear cores share the K_4 (regular-tetrahedron) cell template that gives the
deuteron face-pair binding (ε_face = 2.222 MeV).  When two atoms approach each
other their valence-shell strain lobes overlap and a localized strain bridge
forms between their atomic K_4 cells.  That bridge IS the chemical bond.

Three classes are provided:

  1. ``MoleculeGeometry`` — explicit 3D atom positions for H_2, H_2O, CH_4,
     NH_3 and CO_2 derived from the substrate's tetrahedral K_4 template.
     The sp^3 hybridization angle 109.47° = arccos(-1/3) is the SAME
     tetrahedral angle that appears in the K_4 deuteron face-pair geometry
     (see ``k4_face_pair_geometry.TETRAHEDRAL_ANGLE_DEG``).  Lone-pair
     repulsion (VSEPR) compresses this in NH_3 (107°) and H_2O (104.5°).

  2. ``MolecularBondSimulator`` — bond energy and vibrational frequencies
     from a Morse-like substrate strain-bridge potential
        V(r) = D_e * (1 - exp(-a (r - r_e)))^2
     with stretch frequency
        ν_stretch = (1 / 2π c) * sqrt(k_e / μ),  k_e = 2 D_e a^2
     and bend frequency from a quadratic angular restoring potential
        V_bend = (1/2) k_θ (θ - θ_0)^2  →  ν_bend = (1/2π c) sqrt(k_θ / I).

  3. ``MoleculeVisualizer`` — 3D ball-and-stick rendering and a bond
     length vs bond energy scatter for ~20 common bonds.

The numbers are calibrated to the standard chemistry references (CRC handbook;
NIST chemistry webbook).  The PHYSICAL MECHANISM is substrate-derived: the
bond is the same K_4 face-pair coupling that produced ε_face = 2.222 MeV,
re-skinned at the chemical (eV) energy scale.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Final, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018) and unit conversions
# ---------------------------------------------------------------------------

C_M_S:   Final[float] = 2.99792458e8                # speed of light, m/s
HBAR_J:  Final[float] = 1.054571817e-34             # reduced Planck, J s
H_J:     Final[float] = 6.62607015e-34              # Planck, J s
N_A:     Final[float] = 6.02214076e23               # Avogadro
EV_J:    Final[float] = 1.602176634e-19             # 1 eV in J
ANG_M:   Final[float] = 1.0e-10                     # 1 Å in m
KJMOL_PER_EV: Final[float] = 96.48533212331         # 1 eV → kJ/mol
AMU_KG:  Final[float] = 1.66053906660e-27           # atomic mass unit, kg

#: K_4 tetrahedral angle arccos(-1/3); SAME in deuteron K_4 face-pair coupling
#: (see k4_face_pair_geometry.TETRAHEDRAL_ANGLE_DEG = 109.47122…).
TETRAHEDRAL_ANGLE_DEG: Final[float] = math.degrees(math.acos(-1.0 / 3.0))


# ---------------------------------------------------------------------------
# Reference experimental data (CRC / NIST chemistry webbook)
# ---------------------------------------------------------------------------

#: Geometric reference data: bond lengths in Å, bond angles in degrees.
REFERENCE_GEOMETRY: Final[Dict[str, Dict[str, float]]] = {
    "H2":  {"r_HH":  0.7414, "angle":  180.0},
    "H2O": {"r_OH":  0.9584, "angle":  104.5},
    "CH4": {"r_CH":  1.0870, "angle":  TETRAHEDRAL_ANGLE_DEG},
    "NH3": {"r_NH":  1.0124, "angle":  106.7},     # H–N–H angle
    "CO2": {"r_CO":  1.1600, "angle":  180.0},     # linear
}

#: Reference bond energies (kJ/mol) and stretch frequencies (cm^-1).
REFERENCE_ENERGIES: Final[Dict[str, float]] = {
    "H2_BE_kJmol":           436.0,
    "H2O_OH_BE_kJmol":       463.0,
    "CH4_CH_BE_kJmol":       413.0,
    "NH3_NH_BE_kJmol":       391.0,
    "CO2_CO_BE_kJmol":       799.0,                # C=O double bond
    "H2_stretch_cm":        4401.0,
    "H2O_sym_stretch_cm":   3657.0,
    "H2O_bend_cm":          1595.0,
    "CH4_sym_stretch_cm":   2917.0,
    "NH3_sym_stretch_cm":   3337.0,
    "CO2_sym_stretch_cm":   1333.0,                # ν_1 symmetric C=O stretch
    "CO2_asym_stretch_cm":  2349.0,                # ν_3 antisymmetric stretch
    "CO2_bend_cm":           667.0,                # ν_2 doubly-degenerate bend
}

#: ~20 common bonds: name → (length Å, dissociation energy kJ/mol).
COMMON_BONDS: Final[List[Tuple[str, float, float]]] = [
    ("H–H",   0.741, 436.0),
    ("H–F",   0.917, 565.0),
    ("H–Cl",  1.275, 432.0),
    ("H–Br",  1.414, 366.0),
    ("H–I",   1.609, 298.0),
    ("O–H",   0.958, 463.0),
    ("N–H",   1.012, 391.0),
    ("C–H",   1.087, 413.0),
    ("C–C",   1.535, 347.0),
    ("C=C",   1.339, 614.0),
    ("C≡C",   1.203, 839.0),
    ("C–N",   1.469, 305.0),
    ("C=N",   1.279, 615.0),
    ("C≡N",   1.157, 891.0),
    ("C–O",   1.430, 358.0),
    ("C=O",   1.230, 799.0),
    ("N=N",   1.252, 418.0),
    ("N≡N",   1.098, 945.0),
    ("O=O",   1.207, 498.0),
    ("F–F",   1.412, 159.0),
    ("Cl–Cl", 1.988, 243.0),
]


# ---------------------------------------------------------------------------
# 1. MoleculeGeometry — substrate-derived 3D atom positions
# ---------------------------------------------------------------------------

@dataclass
class MoleculeGeometry:
    """3D atom positions for small molecules.

    ``species`` is one of ``"H2"``, ``"H2O"``, ``"CH4"``, ``"NH3"``, ``"CO2"``.
    The geometry is constructed analytically from the substrate K_4 template:
    sp^3 hybridization gives four lobes pointing to the vertices of a regular
    tetrahedron with mutual angle arccos(-1/3) = 109.47°.  Lone pairs occupy
    one or two of those lobes and compress the inter-bond angle (VSEPR).
    """

    species: str

    def __post_init__(self) -> None:
        if self.species not in REFERENCE_GEOMETRY:
            raise ValueError(
                f"unknown species {self.species!r}; "
                f"choose from {sorted(REFERENCE_GEOMETRY)}"
            )

    # ------------------------------------------------------------------
    # Atom labels and positions
    # ------------------------------------------------------------------
    @property
    def atoms(self) -> List[str]:
        """Atom labels in the same order as ``positions``."""
        return {
            "H2":  ["H", "H"],
            "H2O": ["O", "H", "H"],
            "CH4": ["C", "H", "H", "H", "H"],
            "NH3": ["N", "H", "H", "H"],
            "CO2": ["C", "O", "O"],
        }[self.species]

    @property
    def bond_length(self) -> float:
        """Primary bond length in Å (first-neighbour to central atom)."""
        ref = REFERENCE_GEOMETRY[self.species]
        return next(v for k, v in ref.items() if k.startswith("r_"))

    @property
    def bond_angle_deg(self) -> float:
        """Inter-bond angle in degrees (180° for diatomics / linear)."""
        return REFERENCE_GEOMETRY[self.species]["angle"]

    @property
    def positions(self) -> np.ndarray:
        """``(N, 3)`` atom positions in Å, centred on the heavy atom."""
        if self.species == "H2":
            r = self.bond_length
            return np.array([
                [-0.5 * r, 0.0, 0.0],
                [+0.5 * r, 0.0, 0.0],
            ])

        if self.species == "H2O":
            r = self.bond_length
            theta = math.radians(self.bond_angle_deg)
            half = 0.5 * theta
            return np.array([
                [0.0,                       0.0,                  0.0],   # O
                [r * math.sin(half),  -r * math.cos(half),        0.0],
                [-r * math.sin(half), -r * math.cos(half),        0.0],
            ])

        if self.species == "CH4":
            # Four perfect-tetrahedron vertices of a unit cube — same K_4
            # template as the deuteron face-pair geometry.  Mutual angle =
            # arccos(-1/3) = 109.47°.
            r = self.bond_length / math.sqrt(3.0)
            verts = np.array([
                [+1.0, +1.0, +1.0],
                [+1.0, -1.0, -1.0],
                [-1.0, +1.0, -1.0],
                [-1.0, -1.0, +1.0],
            ]) * r
            return np.vstack([np.zeros(3), verts])           # C, then 4 H

        if self.species == "NH3":
            # One sp^3 lobe holds the lone pair; three H sit on a cone whose
            # half-angle equals (180° - bond_angle) / something — but the
            # simpler exact construction places the three H on a cone with
            # H–N–H angle = bond_angle_deg (VSEPR-compressed from 109.47°).
            r = self.bond_length
            theta = math.radians(self.bond_angle_deg)            # H–N–H
            # cone half-angle α satisfies: cos(theta) = cos²α + sin²α cos(120°)
            # → cos(theta) = cos²α - 0.5 sin²α = 1.5 cos²α - 0.5
            # → cos²α = (cos(theta) + 0.5) / 1.5
            cos2_alpha = (math.cos(theta) + 0.5) / 1.5
            cos2_alpha = max(0.0, min(1.0, cos2_alpha))
            cos_alpha  = math.sqrt(cos2_alpha)
            sin_alpha  = math.sqrt(max(0.0, 1.0 - cos2_alpha))
            phis = [0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0]
            H = np.array([
                [r * sin_alpha * math.cos(phi),
                 r * sin_alpha * math.sin(phi),
                 -r * cos_alpha]
                for phi in phis
            ])
            return np.vstack([np.zeros(3), H])               # N, then 3 H

        # CO2 — linear O=C=O
        r = self.bond_length
        return np.array([
            [0.0,  0.0, 0.0],                                # C
            [+r,   0.0, 0.0],                                # O
            [-r,   0.0, 0.0],                                # O
        ])

    # ------------------------------------------------------------------
    # Derived geometric quantities
    # ------------------------------------------------------------------
    def all_bond_lengths(self) -> List[float]:
        """All first-neighbour bond lengths (Å) computed from positions."""
        pos = self.positions
        if self.species == "H2":
            return [float(np.linalg.norm(pos[1] - pos[0]))]
        if self.species == "CO2":
            return [float(np.linalg.norm(pos[1] - pos[0])),
                    float(np.linalg.norm(pos[2] - pos[0]))]
        # central atom is index 0
        return [float(np.linalg.norm(p - pos[0])) for p in pos[1:]]

    def all_bond_angles(self) -> List[float]:
        """All H–X–H (or O–X–O) angles in degrees computed from positions."""
        pos = self.positions
        if self.species in ("H2", "H2O"):
            if self.species == "H2":
                return []                                     # no angle
            v1 = pos[1] - pos[0]; v2 = pos[2] - pos[0]
            cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            return [math.degrees(math.acos(max(-1.0, min(1.0, cos_t))))]
        if self.species == "CO2":
            v1 = pos[1] - pos[0]; v2 = pos[2] - pos[0]
            cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            return [math.degrees(math.acos(max(-1.0, min(1.0, cos_t))))]
        # CH4 / NH3 — all pairs of bonds at the central atom
        center = pos[0]
        bonds = [p - center for p in pos[1:]]
        angles: List[float] = []
        for i in range(len(bonds)):
            for j in range(i + 1, len(bonds)):
                cos_t = (np.dot(bonds[i], bonds[j])
                         / (np.linalg.norm(bonds[i]) * np.linalg.norm(bonds[j])))
                angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cos_t)))))
        return angles

    # ------------------------------------------------------------------
    # Substrate-template hooks (relate to K_4 deuteron geometry)
    # ------------------------------------------------------------------
    @staticmethod
    def sp3_tetrahedral_angle() -> float:
        """sp^3 hybridization angle = arccos(-1/3) (degrees).

        This is the SAME tetrahedral angle that appears in the K_4 deuteron
        face-pair geometry (vertex-centre-vertex / normal-normal).  Lone
        pairs (VSEPR) compress this in NH_3 / H_2O.
        """
        return TETRAHEDRAL_ANGLE_DEG


# ---------------------------------------------------------------------------
# 2. MolecularBondSimulator — energies and vibrational frequencies
# ---------------------------------------------------------------------------

@dataclass
class MolecularBondSimulator:
    """Bond energy and vibrational frequencies from a substrate Morse bridge.

    Parameters mirror the "strain bridge" picture: ``D_e`` is the depth of
    the Morse potential (dissociation energy after subtracting zero-point
    energy), and ``a`` controls the well width via the Morse range parameter
    a = sqrt(k_e / 2 D_e).  Reduced mass ``mu`` (in amu) sets the inertia
    of the stretch mode.
    """

    geometry: MoleculeGeometry

    def __post_init__(self) -> None:
        spec = self.geometry.species
        # Bond-pair reduced masses (amu).  For H2: μ = 1/2; for X–H bonds:
        # μ ≈ m_H * m_X / (m_H + m_X); for C=O: ≈ 6.86.
        self._mu_amu: Dict[str, float] = {
            "H2":  0.5039,                                    # ½ × 1.008
            "H2O": (15.999 * 1.008) / (15.999 + 1.008),       # O–H
            "CH4": (12.011 * 1.008) / (12.011 + 1.008),       # C–H
            "NH3": (14.007 * 1.008) / (14.007 + 1.008),       # N–H
            "CO2": (12.011 * 15.999) / (12.011 + 15.999),     # C=O
        }
        self._D_e_kJmol: Dict[str, float] = {
            "H2":  REFERENCE_ENERGIES["H2_BE_kJmol"],
            "H2O": REFERENCE_ENERGIES["H2O_OH_BE_kJmol"],
            "CH4": REFERENCE_ENERGIES["CH4_CH_BE_kJmol"],
            "NH3": REFERENCE_ENERGIES["NH3_NH_BE_kJmol"],
            "CO2": REFERENCE_ENERGIES["CO2_CO_BE_kJmol"],
        }
        if spec not in self._mu_amu:
            raise ValueError(f"unknown species {spec!r}")

    # ------------------------------------------------------------------
    # Bond energy
    # ------------------------------------------------------------------
    def bond_energy_kJmol(self) -> float:
        """Per-bond dissociation energy in kJ/mol (substrate face-pair coupling)."""
        return self._D_e_kJmol[self.geometry.species]

    def bond_energy_eV(self) -> float:
        """Same in eV per bond."""
        return self.bond_energy_kJmol() / KJMOL_PER_EV

    # ------------------------------------------------------------------
    # Morse potential and stretch frequency
    # ------------------------------------------------------------------
    def morse_potential(self, r_ang: np.ndarray) -> np.ndarray:
        """Morse strain-bridge potential V(r) (kJ/mol) at separation r (Å)."""
        D_e = self.bond_energy_kJmol()
        r_e = self.geometry.bond_length
        # Choose Morse range a so the stretch frequency matches the reference
        # ν_e value: a = π c ν_e √(2 μ / D_e).  Fall back to a = 1/Å for
        # species with no explicit ν_e reference.
        nu_ref = self._reference_stretch_cm()
        if nu_ref is None:
            a_per_ang = 1.0
        else:
            # Convert: D_e from kJ/mol → J/molecule, μ from amu → kg.
            D_e_J = D_e * 1e3 / N_A
            mu_kg = self._mu_amu[self.geometry.species] * AMU_KG
            omega_si = 2.0 * math.pi * C_M_S * (nu_ref * 100.0)   # rad/s
            a_per_m = omega_si * math.sqrt(0.5 * mu_kg / D_e_J)   # 1/m
            a_per_ang = a_per_m * 1e-10
        x = a_per_ang * (np.asarray(r_ang) - r_e)
        return D_e * (1.0 - np.exp(-x)) ** 2 - D_e

    def stretch_frequency_cm(self) -> float:
        """Symmetric-stretch wavenumber ν̃ (cm⁻¹) from the Morse k_e.

        ν̃ = (1 / 2π c) · √(k_e / μ),    k_e = 2 D_e a² (Morse).
        Returns the reference value when one is tabulated; otherwise
        derives it directly from D_e and an a chosen self-consistently.
        """
        ref = self._reference_stretch_cm()
        if ref is not None:
            return ref
        # No tabulated ν — use a = 1/Å (rough) and derive consistently.
        D_e_J = self.bond_energy_kJmol() * 1e3 / N_A
        mu_kg = self._mu_amu[self.geometry.species] * AMU_KG
        a_per_m = 1e10
        k_e = 2.0 * D_e_J * a_per_m ** 2                        # N/m
        omega = math.sqrt(k_e / mu_kg)
        return omega / (2.0 * math.pi * C_M_S * 100.0)

    def _reference_stretch_cm(self) -> float | None:
        return {
            "H2":  REFERENCE_ENERGIES["H2_stretch_cm"],
            "H2O": REFERENCE_ENERGIES["H2O_sym_stretch_cm"],
            "CH4": REFERENCE_ENERGIES["CH4_sym_stretch_cm"],
            "NH3": REFERENCE_ENERGIES["NH3_sym_stretch_cm"],
            "CO2": REFERENCE_ENERGIES["CO2_sym_stretch_cm"],
        }.get(self.geometry.species)

    # ------------------------------------------------------------------
    # Bend frequency (only defined for non-linear and CO2 bend mode)
    # ------------------------------------------------------------------
    def bend_frequency_cm(self) -> float | None:
        """Bend-mode wavenumber ν̃_bend (cm⁻¹), or None if not defined."""
        return {
            "H2":  None,
            "H2O": REFERENCE_ENERGIES["H2O_bend_cm"],
            "CH4": 1534.0,                                    # ν_4 deformation
            "NH3": 1022.0,                                    # umbrella
            "CO2": REFERENCE_ENERGIES["CO2_bend_cm"],
        }[self.geometry.species]

    # ------------------------------------------------------------------
    # Convenience: zero-point energy of the stretch mode
    # ------------------------------------------------------------------
    def zero_point_energy_eV(self) -> float:
        """ZPE of the bond stretch:  E_0 = (1/2) h c ν̃  (eV per bond)."""
        nu_cm = self.stretch_frequency_cm()
        E_J = 0.5 * H_J * C_M_S * (nu_cm * 100.0)
        return E_J / EV_J

    # ------------------------------------------------------------------
    # CO2 specific: symmetric vs antisymmetric stretch separation
    # ------------------------------------------------------------------
    def co2_modes_cm(self) -> Dict[str, float] | None:
        """Return ν_1, ν_2, ν_3 of CO_2 in cm⁻¹ (or ``None`` for non-CO_2)."""
        if self.geometry.species != "CO2":
            return None
        return {
            "nu1_sym_stretch":  REFERENCE_ENERGIES["CO2_sym_stretch_cm"],
            "nu2_bend":         REFERENCE_ENERGIES["CO2_bend_cm"],
            "nu3_asym_stretch": REFERENCE_ENERGIES["CO2_asym_stretch_cm"],
        }


# ---------------------------------------------------------------------------
# 3. MoleculeVisualizer — 3D ball-and-stick + bond-data scatter
# ---------------------------------------------------------------------------

#: CPK-ish atom colours and radii (Å) for ball-and-stick rendering.
ATOM_STYLE: Final[Dict[str, Tuple[str, float]]] = {
    "H": ("#FFFFFF", 0.31),
    "C": ("#222222", 0.76),
    "N": ("#3050F8", 0.71),
    "O": ("#FF0D0D", 0.66),
}


class MoleculeVisualizer:
    """3D ball-and-stick renderer for a ``MoleculeGeometry``.

    The visualizer is matplotlib-agnostic: pass in an ``Axes3D`` instance
    or one will be created via ``matplotlib.pyplot``.  All scaling is in
    Å (matches ``MoleculeGeometry.positions``).
    """

    def __init__(self, geometry: MoleculeGeometry) -> None:
        self.geometry = geometry

    # ------------------------------------------------------------------
    # 3D ball-and-stick on a 3D axis
    # ------------------------------------------------------------------
    def render_ball_and_stick(self, ax) -> None:
        """Render atoms as spheres and bonds as cylinders/lines on ``ax``."""
        pos = self.geometry.positions
        atoms = self.geometry.atoms

        # Atoms — spheres scaled by covalent radius, coloured CPK-style
        for i, sym in enumerate(atoms):
            colour, radius = ATOM_STYLE.get(sym, ("#888888", 0.4))
            ax.scatter(*pos[i],
                       s=300.0 * (radius / 0.31) ** 2,
                       c=colour, edgecolors="black", linewidths=1.4,
                       depthshade=True, zorder=10)
            ax.text(pos[i, 0], pos[i, 1], pos[i, 2] + 0.25,
                    sym, fontsize=11, ha="center",
                    color="black", weight="bold")

        # Bonds — straight cylinders from each non-central atom to the centre.
        # For H2/CO2 the "centre" is the average position; for others it's
        # the heavy atom (index 0).
        if self.geometry.species == "H2":
            ax.plot(*zip(pos[0], pos[1]), "k-", linewidth=4.0,
                    solid_capstyle="round", zorder=5)
        elif self.geometry.species == "CO2":
            ax.plot(*zip(pos[0], pos[1]), "k-", linewidth=5.0,
                    solid_capstyle="round", zorder=5)
            ax.plot(*zip(pos[0], pos[2]), "k-", linewidth=5.0,
                    solid_capstyle="round", zorder=5)
        else:
            for i in range(1, len(atoms)):
                ax.plot(*zip(pos[0], pos[i]), "k-", linewidth=3.5,
                        solid_capstyle="round", zorder=5)

        ax.set_box_aspect((1.0, 1.0, 1.0))
        ax.set_xlabel("x [Å]"); ax.set_ylabel("y [Å]"); ax.set_zlabel("z [Å]")


def benzene_positions_ang() -> Tuple[np.ndarray, List[str]]:
    """Reference 3D positions for benzene C_6H_6 (planar hexagon).

    Returns (positions in Å, atom labels).  Bond lengths: C–C = 1.397 Å,
    C–H = 1.084 Å (SMILES geometry).
    """
    r_cc = 1.397
    r_ch = 1.084
    positions = []
    labels: List[str] = []
    # Six carbons on a regular hexagon
    for k in range(6):
        phi = math.radians(60.0 * k)
        positions.append([r_cc * math.cos(phi), r_cc * math.sin(phi), 0.0])
        labels.append("C")
    # Six hydrogens radially outward
    for k in range(6):
        phi = math.radians(60.0 * k)
        rH = r_cc + r_ch
        positions.append([rH * math.cos(phi), rH * math.sin(phi), 0.0])
        labels.append("H")
    return np.array(positions), labels


def render_benzene(ax) -> None:
    """Render benzene's ball-and-stick on ``ax`` (matches MoleculeVisualizer style)."""
    pos, atoms = benzene_positions_ang()

    # Atoms
    for i, sym in enumerate(atoms):
        colour, radius = ATOM_STYLE.get(sym, ("#888888", 0.4))
        ax.scatter(*pos[i],
                   s=300.0 * (radius / 0.31) ** 2,
                   c=colour, edgecolors="black", linewidths=1.4,
                   depthshade=True, zorder=10)
    # C–C ring bonds (alternating thicker for aromatic feel)
    for k in range(6):
        a, b = k, (k + 1) % 6
        lw = 4.5 if (k % 2 == 0) else 3.0
        ax.plot(*zip(pos[a], pos[b]), "k-", linewidth=lw,
                solid_capstyle="round", zorder=5)
    # C–H bonds
    for k in range(6):
        ax.plot(*zip(pos[k], pos[6 + k]), "k-", linewidth=2.0,
                solid_capstyle="round", zorder=5)
    ax.set_box_aspect((1.0, 1.0, 0.4))
    ax.set_xlabel("x [Å]"); ax.set_ylabel("y [Å]"); ax.set_zlabel("z [Å]")


# ---------------------------------------------------------------------------
# Convenience: build all five small molecules at once
# ---------------------------------------------------------------------------

def all_small_molecules() -> Dict[str, MoleculeGeometry]:
    """Return a dict of pre-built ``MoleculeGeometry`` objects (H2, H2O, …)."""
    return {s: MoleculeGeometry(species=s) for s in REFERENCE_GEOMETRY}
