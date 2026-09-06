"""Madelung-constant test of the substrate K_4 close-packed lattice ansatz.

The substrate framework asserts that ionic crystals are particular tilings
of the K_4 close-packed cell, and that their Madelung constants — the
dimensionless lattice sums

    M = - V(test ion) * r_NN / Z(test ion)        (per-ion convention)

are *purely geometric* outputs of the lattice arrangement plus Coulomb
1/r summation.  No free parameters are introduced; every prediction is
the direct lattice sum for the published geometry of each ionic crystal
(Ewald-summed for full convergence).

WHAT THIS MODULE DOES
---------------------
For each of eight ionic-crystal structure types — NaCl, CsCl, ZnS
(sphalerite), ZnS (wurtzite), CaF_2 (fluorite), TiO_2 (rutile), Cu_2O
(cuprite), and the simple-cubic CsCl 'beta' alternative — it

  1. defines the conventional unit cell (lattice vectors + ion basis);
  2. computes the Madelung constant from a Ewald lattice sum (real-space
     erfc-screened plus reciprocal-space Gaussian-screened parts) — the
     standard exact method for periodic charge distributions;
  3. compares each prediction against the published literature value
     supplied by the task spec; and
  4. reports the relative error per crystal type, plus a verdict on which
     structures the substrate-direct-sum lattice ansatz reproduces.

  An Evjen direct-sum mode (`use_evjen=True`) is also exposed.  Direct
  cubic-shell summation converges only for charge-balanced cell tilings
  (NaCl rocksalt); for CsCl, ZnS, CaF_2, TiO_2 and Cu_2O the conventional-
  cell tiling is conditionally convergent and direct sums diverge or
  converge to the wrong limit.  Ewald removes that pathology and is
  therefore the default method here.

CONVENTIONS
-----------
"per-ion" Madelung (the Sherman / Tosi convention reported in textbooks
for NaCl, CsCl, ZnS):  potential at one ion site, normalised so that
    V(test ion) = -M_per_ion * Z_self * e^2 / (4 pi eps0 r_NN).

"per-formula-unit" Madelung (the convention used in modern lattice-energy
references for asymmetric salts CaF_2, TiO_2, Cu_2O):
    U_latt = -M_per_FU * e^2 / (4 pi eps0 r_NN)            per FU,
    M_per_FU = -(1/2) sum_{i in FU} Z_i V_i r_NN / e^2*(4 pi eps0).

For NaCl, CsCl, ZnS (where |Z+| = |Z-|), the two conventions coincide
(M_per_ion = M_per_FU).  For CaF_2 (Z+=2, Z-=-1, two F per Ca), they
differ — published M_per_FU(CaF_2) = 5.0388 corresponds to
M_per_anion(CaF_2) = 2.51939, etc.

The substrate ansatz makes a *single* prediction per crystal: the Coulomb
lattice sum on the published geometry.  We report it under both
conventions so the user can compare to whichever literature value they
have at hand.

REFERENCES
----------
- Sherman, Chem. Rev. 1932 — original Madelung tabulation for NaCl,
  CsCl, ZnS, CaF_2, TiO_2.
- Tosi, in *Solid State Physics* vol. 16 (1964) — Ewald-summation
  high-precision values.
- Evjen, Phys. Rev. 39, 675 (1932) — direct cubic-shell weighting trick.
- Cu_2O Madelung: O'Keeffe & Bovin, J. Chem. Phys. 73, 4922 (1980)
  reports per-formula-unit ~ 4.4421.
- Rutile Madelung per formula unit ~ 4.816 with c/a = 0.6442 and
  internal parameter u = 0.3056 (Baur 1956; Hund 1948).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np


# --------------------------------------------------------------------------- #
# Published reference values (from spec + literature cross-check)             #
# --------------------------------------------------------------------------- #

PUBLISHED: Dict[str, Dict[str, float]] = {
    "NaCl": {
        "M_per_ion": 1.747565,
        "M_per_FU":  1.747565,
        "convention": "per_ion",
    },
    "CsCl": {
        "M_per_ion": 1.762675,
        "M_per_FU":  1.762675,
        "convention": "per_ion",
    },
    "ZnS_sphalerite": {
        "M_per_ion": 1.6381,
        "M_per_FU":  1.6381,
        "convention": "per_ion",
    },
    "ZnS_wurtzite": {
        "M_per_ion": 1.6413,
        "M_per_FU":  1.6413,
        "convention": "per_ion",
    },
    "CaF2_fluorite": {
        "M_per_ion": 2.51939,
        "M_per_FU":  5.03878,
        "convention": "per_FU",
    },
    "TiO2_rutile": {
        "M_per_ion": 2.408,
        "M_per_FU":  4.816,
        "convention": "per_FU",
    },
    "Cu2O_cuprite": {
        "M_per_ion": 2.2210,
        "M_per_FU":  4.4421,
        "convention": "per_FU",
    },
    "beta_CsCl_SC": {
        # 'beta-CsCl' = simple-cubic NaCl-arrangement (Cs+ at corners,
        # Cl- on alternating cubic vertices).  Equivalent to NaCl on SC.
        "M_per_ion": 1.747565,
        "M_per_FU":  1.747565,
        "convention": "per_ion",
    },
}


# --------------------------------------------------------------------------- #
# Crystal geometries: returns (lattice vectors a, list of (frac_pos, charge)) #
# --------------------------------------------------------------------------- #


def _nacl_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """NaCl rocksalt: 4 Na+ at FCC + 4 Cl- shifted by (1/2,1/2,1/2)."""
    a = np.eye(3) * 1.0
    basis: List[Tuple[np.ndarray, float]] = []
    fcc_offsets = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    ]
    for o in fcc_offsets:
        basis.append((np.array(o), +1.0))
    cl_offsets = [
        (0.5, 0.5, 0.5),
        (0.0, 0.0, 0.5),
        (0.0, 0.5, 0.0),
        (0.5, 0.0, 0.0),
    ]
    for o in cl_offsets:
        basis.append((np.array(o), -1.0))
    return a, basis


def _cscl_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """CsCl: simple cubic Cs+ at corners + Cl- at body centre."""
    a = np.eye(3) * 1.0
    basis = [
        (np.array([0.0, 0.0, 0.0]), +1.0),
        (np.array([0.5, 0.5, 0.5]), -1.0),
    ]
    return a, basis


def _beta_cscl_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """beta-CsCl: the high-temperature (T > 469 deg C) NaCl-rocksalt
    polymorph of CsCl, where Cs+ and Cl- swap from CN=8 (CsCl phase) to
    CN=6 (rocksalt).  Madelung constant is then the rocksalt value
    1.747565 — published as a literature alternative for CsCl.
    """
    a = np.eye(3) * 1.0
    fcc_offsets = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    ]
    basis: List[Tuple[np.ndarray, float]] = []
    for o in fcc_offsets:
        basis.append((np.array(o), +1.0))           # Cs+
    cl_offsets = [
        (0.5, 0.5, 0.5),
        (0.0, 0.0, 0.5),
        (0.0, 0.5, 0.0),
        (0.5, 0.0, 0.0),
    ]
    for o in cl_offsets:
        basis.append((np.array(o), -1.0))           # Cl-
    return a, basis


def _zns_sphalerite_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """Zincblende ZnS: FCC S^2- + FCC Zn^2+ offset by (1/4,1/4,1/4)."""
    a = np.eye(3) * 1.0
    fcc_offsets = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    ]
    basis: List[Tuple[np.ndarray, float]] = []
    for o in fcc_offsets:
        # Use unit charges so M_per_ion compares directly to 1.6381 (ZnS).
        basis.append((np.array(o), -1.0))           # S
    for o in fcc_offsets:
        basis.append((np.array([o[0] + 0.25,
                                o[1] + 0.25,
                                o[2] + 0.25]), +1.0))   # Zn
    return a, basis


def _zns_wurtzite_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """Wurtzite ZnS: hexagonal P63mc, ideal c/a = sqrt(8/3), u = 3/8."""
    c_over_a = math.sqrt(8.0 / 3.0)
    u = 0.375
    a_vecs = np.array([
        [1.0,           0.0,                0.0],
        [-0.5,          math.sqrt(3.0) / 2.0, 0.0],
        [0.0,           0.0,                c_over_a],
    ])
    basis = [
        (np.array([1.0 / 3.0, 2.0 / 3.0, 0.0]),         -1.0),
        (np.array([2.0 / 3.0, 1.0 / 3.0, 0.5]),         -1.0),
        (np.array([1.0 / 3.0, 2.0 / 3.0, u]),           +1.0),
        (np.array([2.0 / 3.0, 1.0 / 3.0, 0.5 + u]),     +1.0),
    ]
    return a_vecs, basis


def _caf2_fluorite_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """Fluorite CaF_2: Ca^2+ on FCC, F^- at all 8 tetrahedral sites."""
    a = np.eye(3) * 1.0
    fcc_offsets = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    ]
    basis: List[Tuple[np.ndarray, float]] = []
    for o in fcc_offsets:
        basis.append((np.array(o), +2.0))
    f_offsets = [
        (0.25, 0.25, 0.25),
        (0.25, 0.25, 0.75),
        (0.25, 0.75, 0.25),
        (0.75, 0.25, 0.25),
        (0.25, 0.75, 0.75),
        (0.75, 0.25, 0.75),
        (0.75, 0.75, 0.25),
        (0.75, 0.75, 0.75),
    ]
    for o in f_offsets:
        basis.append((np.array(o), -1.0))
    return a, basis


def _tio2_rutile_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """Rutile TiO_2: tetragonal P4_2/mnm, c/a = 0.6442, u = 0.3056.

    2 Ti^4+ + 4 O^2- per cell.
    """
    c_over_a = 0.6442
    u = 0.3056
    a_vecs = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, c_over_a],
    ])
    basis = [
        (np.array([0.0, 0.0, 0.0]),                     +4.0),
        (np.array([0.5, 0.5, 0.5]),                     +4.0),
        (np.array([u,   u,   0.0]) % 1.0,               -2.0),
        (np.array([1.0 - u, 1.0 - u, 0.0]) % 1.0,       -2.0),
        (np.array([0.5 + u, 0.5 - u, 0.5]) % 1.0,       -2.0),
        (np.array([0.5 - u, 0.5 + u, 0.5]) % 1.0,       -2.0),
    ]
    return a_vecs, basis


def _cu2o_cuprite_cell() -> Tuple[np.ndarray, List[Tuple[np.ndarray, float]]]:
    """Cuprite Cu_2O: cubic Pn-3m.

    Two formula units (4 Cu+ + 2 O^2-) per unit cell.
    O at (0,0,0) and (1/2,1/2,1/2).
    Cu at (1/4,1/4,1/4), (1/4,3/4,3/4), (3/4,1/4,3/4), (3/4,3/4,1/4).
    """
    a = np.eye(3) * 1.0
    basis = [
        (np.array([0.0, 0.0, 0.0]),         -2.0),
        (np.array([0.5, 0.5, 0.5]),         -2.0),
        (np.array([0.25, 0.25, 0.25]),      +1.0),
        (np.array([0.25, 0.75, 0.75]),      +1.0),
        (np.array([0.75, 0.25, 0.75]),      +1.0),
        (np.array([0.75, 0.75, 0.25]),      +1.0),
    ]
    return a, basis


CRYSTAL_CONSTRUCTORS = {
    "NaCl":             _nacl_cell,
    "CsCl":             _cscl_cell,
    "beta_CsCl_SC":     _beta_cscl_cell,
    "ZnS_sphalerite":   _zns_sphalerite_cell,
    "ZnS_wurtzite":     _zns_wurtzite_cell,
    "CaF2_fluorite":    _caf2_fluorite_cell,
    "TiO2_rutile":      _tio2_rutile_cell,
    "Cu2O_cuprite":     _cu2o_cuprite_cell,
}


# --------------------------------------------------------------------------- #
# Lattice geometry helpers                                                    #
# --------------------------------------------------------------------------- #


def _cell_volume(a_vecs: np.ndarray) -> float:
    return float(abs(np.linalg.det(a_vecs)))


def _reciprocal_vectors(a_vecs: np.ndarray) -> np.ndarray:
    """2 pi (a_j x a_k) / V cyclic.  Returns 3x3 array, rows = b_i."""
    V = _cell_volume(a_vecs)
    a1, a2, a3 = a_vecs[0], a_vecs[1], a_vecs[2]
    b1 = 2.0 * math.pi * np.cross(a2, a3) / V
    b2 = 2.0 * math.pi * np.cross(a3, a1) / V
    b3 = 2.0 * math.pi * np.cross(a1, a2) / V
    return np.vstack([b1, b2, b3])


# --------------------------------------------------------------------------- #
# Direct (Evjen-weighted) sum — works only for NaCl-style                     #
# --------------------------------------------------------------------------- #


def _evjen_weight(ijk: Sequence[int], N: int) -> float:
    """Evjen boundary weight: 1/2 per face, 1/4 per edge, 1/8 per vertex."""
    w = 1.0
    for c in ijk:
        if abs(c) == N:
            w *= 0.5
    return w


def _site_potential_direct(
    a_vecs: np.ndarray,
    basis: List[Tuple[np.ndarray, float]],
    site_index: int,
    N_shells: int,
    use_evjen: bool = True,
) -> float:
    """Direct sum sum_{i!=self} q_i / r_i, cubic-shell truncated.

    Use only for charge-balanced cubic-cell tilings (NaCl); diverges for
    CsCl, ZnS, CaF_2, TiO_2, Cu_2O.
    """
    self_pos_frac, _ = basis[site_index]
    self_pos_cart = self_pos_frac @ a_vecs
    V = 0.0
    for i in range(-N_shells, N_shells + 1):
        for j in range(-N_shells, N_shells + 1):
            for k in range(-N_shells, N_shells + 1):
                R = (i * a_vecs[0]) + (j * a_vecs[1]) + (k * a_vecs[2])
                w = _evjen_weight((i, j, k), N_shells) if use_evjen else 1.0
                for n_b, (b_pos_frac, b_q) in enumerate(basis):
                    if i == 0 and j == 0 and k == 0 and n_b == site_index:
                        continue
                    b_pos_cart = b_pos_frac @ a_vecs
                    d = float(np.linalg.norm((b_pos_cart + R) - self_pos_cart))
                    if d > 0.0:
                        V += w * b_q / d
    return V


# --------------------------------------------------------------------------- #
# Ewald summation — exact for any periodic charge distribution                #
# --------------------------------------------------------------------------- #


def _site_potential_ewald(
    a_vecs: np.ndarray,
    basis: List[Tuple[np.ndarray, float]],
    site_index: int,
    eta: float | None = None,
    N_real: int = 6,
    N_recip: int = 6,
) -> float:
    """Coulomb potential at one basis site, Ewald-summed.

    Returns V_self in cell-frame units (a_vecs in length L => V in 1/L).
    Both real-space (erfc-screened) and reciprocal-space (Gaussian-screened)
    parts are summed; the self-energy correction is subtracted.

    Choosing eta ~ 2/(V_cell)^(1/3) gives balanced convergence.
    """
    a_vecs = np.asarray(a_vecs, dtype=float)
    V_cell = _cell_volume(a_vecs)
    if eta is None:
        eta = 2.0 / (V_cell ** (1.0 / 3.0))

    self_frac, self_q = basis[site_index]
    self_cart = self_frac @ a_vecs

    # ---------------- real-space sum ---------------- #
    V_R = 0.0
    for i in range(-N_real, N_real + 1):
        for j in range(-N_real, N_real + 1):
            for k in range(-N_real, N_real + 1):
                R = (i * a_vecs[0]) + (j * a_vecs[1]) + (k * a_vecs[2])
                for n_b, (b_frac, b_q) in enumerate(basis):
                    if i == 0 and j == 0 and k == 0 and n_b == site_index:
                        continue
                    b_cart = b_frac @ a_vecs
                    d = float(np.linalg.norm((b_cart + R) - self_cart))
                    if d > 0.0:
                        V_R += b_q * math.erfc(eta * d) / d

    # ---------------- reciprocal-space sum ---------------- #
    b_vecs = _reciprocal_vectors(a_vecs)
    V_K = 0.0
    for h in range(-N_recip, N_recip + 1):
        for kk in range(-N_recip, N_recip + 1):
            for ll in range(-N_recip, N_recip + 1):
                if h == 0 and kk == 0 and ll == 0:
                    continue
                G = (h * b_vecs[0]) + (kk * b_vecs[1]) + (ll * b_vecs[2])
                G2 = float(np.dot(G, G))
                # Structure factor S(G) including all basis ions (the
                # self-site contributes cos(0) = 1; the self-Gaussian
                # piece is subtracted via the V_self correction below).
                S_re = 0.0
                for b_frac, b_q in basis:
                    b_cart = b_frac @ a_vecs
                    phase = float(np.dot(G, b_cart - self_cart))
                    S_re += b_q * math.cos(phase)
                V_K += (4.0 * math.pi / V_cell) * (
                    math.exp(-G2 / (4.0 * eta * eta)) / G2
                ) * S_re

    # ---------------- self-energy correction ---------------- #
    # Removes the self-Gaussian's contribution that the long-range sum
    # double-counts at zero separation:
    V_self_corr = -2.0 * eta / math.sqrt(math.pi) * self_q

    # Background-charge correction is zero by construction (cells are
    # already charge-neutral; the G=0 term is therefore omitted above).

    return V_R + V_K + V_self_corr


# --------------------------------------------------------------------------- #
# Nearest-neighbour distance                                                  #
# --------------------------------------------------------------------------- #


def _nearest_neighbour_distance(
    a_vecs: np.ndarray,
    basis: List[Tuple[np.ndarray, float]],
    N_shells: int = 2,
) -> float:
    """Find the global shortest unlike-charge distance in the periodic
    crystal."""
    best = math.inf
    for s, (s_frac, s_q) in enumerate(basis):
        s_cart = s_frac @ a_vecs
        for i in range(-N_shells, N_shells + 1):
            for j in range(-N_shells, N_shells + 1):
                for k in range(-N_shells, N_shells + 1):
                    R = (i * a_vecs[0]) + (j * a_vecs[1]) + (k * a_vecs[2])
                    for n_b, (b_frac, b_q) in enumerate(basis):
                        if i == 0 and j == 0 and k == 0 and n_b == s:
                            continue
                        if b_q * s_q >= 0.0:
                            continue
                        b_cart = b_frac @ a_vecs
                        d = float(np.linalg.norm((b_cart + R) - s_cart))
                        if d < best:
                            best = d
    return best


# --------------------------------------------------------------------------- #
# Madelung prediction                                                          #
# --------------------------------------------------------------------------- #


@dataclass
class MadelungResult:
    """One-crystal Madelung-prediction result, both conventions reported."""
    crystal: str
    M_per_FU_predicted: float       # formal charges, lattice-energy convention
    M_per_ion_predicted: float      # M_per_FU / n_anions_per_FU
    M_per_FU_unit_charges: float    # unit-charge geometric convention
    M_per_FU_published: float
    M_per_ion_published: float
    rel_err_per_FU: float
    rel_err_per_ion: float
    rel_err_unit_charges: float     # error of unit-charge prediction vs pub
    convention: str
    rel_err_published: float
    method: str

    def as_dict(self) -> Dict[str, float]:
        return {
            "crystal":              self.crystal,
            "M_per_FU_predicted":   self.M_per_FU_predicted,
            "M_per_ion_predicted":  self.M_per_ion_predicted,
            "M_per_FU_unit_charges":self.M_per_FU_unit_charges,
            "M_per_FU_published":   self.M_per_FU_published,
            "M_per_ion_published":  self.M_per_ion_published,
            "rel_err_per_FU":       self.rel_err_per_FU,
            "rel_err_per_ion":      self.rel_err_per_ion,
            "rel_err_unit_charges": self.rel_err_unit_charges,
            "convention":           self.convention,
            "rel_err_published":    self.rel_err_published,
            "method":               self.method,
        }


def _formula_units_per_cell(crystal: str) -> int:
    table = {
        "NaCl":             4,
        "CsCl":             1,
        "beta_CsCl_SC":     4,    # NaCl-type rocksalt polymorph
        "ZnS_sphalerite":   4,
        "ZnS_wurtzite":     2,
        "CaF2_fluorite":    4,
        "TiO2_rutile":      2,
        "Cu2O_cuprite":     2,
    }
    return table[crystal]


def _ions_per_FU_for_per_ion_convention(crystal: str) -> int:
    """Number of like-charge ions per formula unit used as the divisor in
    converting M_per_FU to M_per_ion.

    For 1:1 salts (NaCl, CsCl, ZnS) it is 1.  For CaF_2 and TiO_2 (2 anions
    per FU) and Cu_2O (2 cations per FU), the per-ion convention divides
    M_per_FU by 2.
    """
    table = {
        "NaCl":             1,
        "CsCl":             1,
        "beta_CsCl_SC":     1,
        "ZnS_sphalerite":   1,
        "ZnS_wurtzite":     1,
        "CaF2_fluorite":    2,    # 2 F per FU
        "TiO2_rutile":      2,    # 2 O per FU
        "Cu2O_cuprite":     2,    # 2 Cu per FU
    }
    return table[crystal]


def _pick_anion_site(basis: List[Tuple[np.ndarray, float]]) -> int:
    return int(min(range(len(basis)), key=lambda s: basis[s][1]))


def madelung_predict(
    crystal: str,
    method: str = "ewald",
    N_real: int = 6,
    N_recip: int = 6,
    N_direct: int = 8,
    eta: float | None = None,
) -> MadelungResult:
    """Substrate-direct-sum Madelung prediction for one crystal type.

    Procedure
    ---------
    1. Build the conventional unit cell (lattice vectors + ion fractional
       positions + charges).
    2. Find the nearest-neighbour distance r_NN (shortest unlike-charge
       separation in the periodic lattice).
    3. Compute V_i = potential at each basis site by Ewald summation
       (default; `method='ewald'`) or by Evjen-weighted direct cubic-shell
       summation (`method='evjen'`, only convergent for NaCl).
    4. Per-ion Madelung (textbook NaCl convention):
         M_per_ion = -V_anion * r_NN / Z_anion.
    5. Per-formula-unit Madelung (asymmetric salts):
         pair_energy_per_cell = (1/2) sum_basis q_i V_i,
         M_per_FU = -pair_energy_per_cell * r_NN / n_FU.

    Returns both predictions plus the published reference and relative
    errors.
    """
    if crystal not in CRYSTAL_CONSTRUCTORS:
        raise KeyError(
            f"unknown crystal '{crystal}'; supported: "
            f"{sorted(CRYSTAL_CONSTRUCTORS)}"
        )
    a_vecs, basis = CRYSTAL_CONSTRUCTORS[crystal]()

    r_NN = _nearest_neighbour_distance(a_vecs, basis, N_shells=2)

    site_potentials: List[float] = []
    for s in range(len(basis)):
        if method == "ewald":
            V = _site_potential_ewald(
                a_vecs, basis, s, eta=eta,
                N_real=N_real, N_recip=N_recip,
            )
        elif method == "evjen":
            V = _site_potential_direct(
                a_vecs, basis, s, N_shells=N_direct, use_evjen=True,
            )
        else:
            raise ValueError(f"unknown method '{method}'")
        site_potentials.append(V)

    n_FU = _formula_units_per_cell(crystal)
    pair_energy = 0.5 * sum(b[1] * V for b, V in zip(basis, site_potentials))
    M_FU = -pair_energy * r_NN / n_FU

    # Per-ion Madelung constant: literature divides M_FU by the number of
    # ions per FU of the more-abundant species (anions for CaF_2, TiO_2;
    # cations for Cu_2O).  For 1:1 salts (NaCl, CsCl, ZnS), the divisor
    # is 1 and M_ion = M_FU.
    n_ions = _ions_per_FU_for_per_ion_convention(crystal)
    M_ion = M_FU / n_ions

    # Also compute the unit-charge geometric Madelung constant — same
    # lattice topology, but every ion replaced by +/-1.  For 1:1 salts
    # (NaCl, CsCl, ZnS) this equals M_FU exactly; for asymmetric salts
    # (CaF_2, TiO_2, Cu_2O) it differs by Z+ * |Z-| factors and is the
    # convention used by some literature tabulations (notably rutile).
    basis_unit = [(p, math.copysign(1.0, q)) for p, q in basis]
    site_pot_unit: List[float] = []
    for s in range(len(basis_unit)):
        if method == "ewald":
            V_unit = _site_potential_ewald(
                a_vecs, basis_unit, s, eta=eta,
                N_real=N_real, N_recip=N_recip,
            )
        else:
            V_unit = _site_potential_direct(
                a_vecs, basis_unit, s, N_shells=N_direct, use_evjen=True,
            )
        site_pot_unit.append(V_unit)
    pair_energy_unit = 0.5 * sum(
        b[1] * V for b, V in zip(basis_unit, site_pot_unit)
    )
    # For unit-charge convention some references use n_FU as divisor
    # (CaF_2, Cu_2O); others use n_cations (TiO_2 rutile).  Report the
    # n_FU-divided value as the canonical unit-charge Madelung.
    M_FU_unit = -pair_energy_unit * r_NN / n_FU

    pub = PUBLISHED[crystal]
    M_FU_pub = pub["M_per_FU"]
    M_ion_pub = pub["M_per_ion"]
    convention = pub["convention"]
    err_FU = abs(M_FU - M_FU_pub) / M_FU_pub
    err_ion = abs(M_ion - M_ion_pub) / M_ion_pub
    err_unit = abs(M_FU_unit - M_FU_pub) / M_FU_pub
    err_conv = err_ion if convention == "per_ion" else err_FU

    return MadelungResult(
        crystal=crystal,
        M_per_FU_predicted=M_FU,
        M_per_ion_predicted=M_ion,
        M_per_FU_unit_charges=M_FU_unit,
        M_per_FU_published=M_FU_pub,
        M_per_ion_published=M_ion_pub,
        rel_err_per_FU=err_FU,
        rel_err_per_ion=err_ion,
        rel_err_unit_charges=err_unit,
        convention=convention,
        rel_err_published=err_conv,
        method=method,
    )


# --------------------------------------------------------------------------- #
# Run all predictions                                                          #
# --------------------------------------------------------------------------- #


def run_test(
    crystals: Sequence[str] | None = None,
    method: str = "ewald",
    N_real: int = 6,
    N_recip: int = 6,
    N_direct: int = 8,
    eta: float | None = None,
) -> Dict[str, object]:
    """Compute and compare Madelung predictions for every supported crystal.

    Returns a dict with:
        rows:      per-crystal MadelungResult dicts
        summary:   {crystal: rel_err in spec convention}
        verdict:   {best, worst, n_within_X, ...}
    """
    if crystals is None:
        crystals = list(CRYSTAL_CONSTRUCTORS.keys())

    rows: Dict[str, Dict[str, float]] = {}
    summary: Dict[str, float] = {}
    for c in crystals:
        r = madelung_predict(
            c, method=method,
            N_real=N_real, N_recip=N_recip, N_direct=N_direct, eta=eta,
        )
        rows[c] = r.as_dict()
        summary[c] = r.rel_err_published

    if summary:
        best = min(summary, key=summary.get)
        worst = max(summary, key=summary.get)
    else:
        best = worst = ""
    verdict = {
        "best":             best,
        "worst":            worst,
        "n_within_0p1pct":  sum(1 for v in summary.values() if v < 1e-3),
        "n_within_1pct":    sum(1 for v in summary.values() if v < 1e-2),
        "n_within_5pct":    sum(1 for v in summary.values() if v < 5e-2),
        "n_total":          len(summary),
    }

    return {
        "rows":     rows,
        "summary":  summary,
        "verdict":  verdict,
        "method":   method,
    }


# --------------------------------------------------------------------------- #
# CLI entrypoint                                                              #
# --------------------------------------------------------------------------- #


def main() -> None:
    print("Substrate K_4 close-packed Madelung-constant test (Ewald summation)")
    print("=" * 78)
    print()
    print(
        "Coulomb lattice sum on the published geometry of each ionic crystal.\n"
        "Both per-ion and per-formula-unit conventions reported; relative\n"
        "error is in the convention the literature uses for that salt."
    )
    print()
    res = run_test(method="ewald", N_real=6, N_recip=6)
    rows = res["rows"]
    print(
        f"{'Crystal':<18}{'M_pred (FU)':>13}{'M_pub (FU)':>13}"
        f"{'M_pred (ion)':>14}{'M_pub (ion)':>13}{'M_unit':>10}"
        f"{'rel err':>10}"
    )
    for name, r in rows.items():
        print(
            f"{name:<18}"
            f"{r['M_per_FU_predicted']:>13.6f}"
            f"{r['M_per_FU_published']:>13.6f}"
            f"{r['M_per_ion_predicted']:>14.6f}"
            f"{r['M_per_ion_published']:>13.6f}"
            f"{r['M_per_FU_unit_charges']:>10.4f}"
            f"{r['rel_err_published'] * 100.0:>9.4f}%"
        )
    print()
    v = res["verdict"]
    print("Verdict:")
    print(f"  best  match: {v['best']:<18} "
          f"err = {res['summary'][v['best']] * 100.0:.4f}%")
    print(f"  worst match: {v['worst']:<18} "
          f"err = {res['summary'][v['worst']] * 100.0:.4f}%")
    print(f"  within 0.1%: {v['n_within_0p1pct']} / {v['n_total']}")
    print(f"  within 1%:   {v['n_within_1pct']} / {v['n_total']}")
    print(f"  within 5%:   {v['n_within_5pct']} / {v['n_total']}")


if __name__ == "__main__":
    main()
