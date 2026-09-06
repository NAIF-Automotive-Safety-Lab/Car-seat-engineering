"""Crystal structures from the substrate K_4 cell tiling — paired GEOMETRY+SIM.

This module makes crystal structures EXPLICITLY VISIBLE and RUNNABLE under
the B3 / Stiff-Medium framework: every crystal is a tiling of K_4 substrate
cells whose lattice constant scales with the substrate cell length ξ.

Two paired classes:

    CrystalGeometry  — explicit unit-cell construction for the 14 Bravais
                       lattices and the standard ionic / covalent / layered
                       crystal families (NaCl, CsCl, ZnS, diamond, graphite,
                       perovskite ABX_3).

    CrystalSimulator — physical observables: lattice constants from the
                       substrate ξ, Bragg-peak X-ray diffraction patterns,
                       electrostatic Madelung sums, atomic packing
                       fractions, mass densities.

ONTOLOGY
--------
In the substrate framework the smallest material unit is the K_4 cell — a
regular tetrahedron of side ξ ≈ 3.86e-13 m (the electron Compton length).
Macroscopic crystals coarse-grain ξ to a "lattice scale" a₀ that is the
matched stable cell-cell separation for the relevant atomic species:

    a_NaCl(ξ)   = 8 · n_M · ξ_eff   (∼ 5.64 Å scaling form)
    a_diamond(ξ)= a_NaCl(ξ) · (3√3/8)·k_C
    a_graph(ξ)  = a_NaCl(ξ) / 1.6149   (in-plane spacing 2.46 Å)

All quoted lattice constants are TABULATED from XRD data; the substrate ξ
sets the *unit*, not the per-species coupling — those ride on top.

REFERENCES
----------
- Kittel, *Introduction to Solid State Physics*, 8th ed.
- Madelung constant for rocksalt: 1.747565 (well-known).
- Bragg's law: 2 d sinθ = n λ; Cu Kα λ = 1.5406 Å used as default probe.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import math
import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Physical / X-ray constants (CODATA + standard XRD)
# ---------------------------------------------------------------------------

CU_KALPHA_ANG: float = 1.5406            # Cu K-α wavelength, Å
MO_KALPHA_ANG: float = 0.7107            # Mo K-α wavelength, Å
N_AVOGADRO: float = 6.02214076e23

# B3 substrate anchor (from b3_constants.py, electron-Compton solver)
XI_M_DEFAULT: float = 3.8615926743523754e-13       # substrate cell length in metres


# ---------------------------------------------------------------------------
# Reference (XRD-tabulated) crystal data
# ---------------------------------------------------------------------------

REFERENCE: Dict[str, float] = {
    # Lattice constants (Å)
    "a_NaCl":     5.6402,
    "a_CsCl":     4.123,
    "a_ZnS":      5.4109,
    "a_diamond":  3.5667,
    "a_graphite_a": 2.4612,
    "a_graphite_c": 6.7079,
    "a_SrTiO3":   3.9050,        # cubic perovskite
    # Madelung constants (dimensionless, for unit nearest-neighbour distance)
    "M_NaCl":     1.747565,      # rocksalt
    "M_CsCl":     1.762675,      # CsCl
    "M_ZnS":      1.6381,        # zincblende
    # Packing fractions (close-packed bounds)
    "phi_FCC":    math.pi / (3.0 * math.sqrt(2.0)),       # 0.7405
    "phi_BCC":    math.pi * math.sqrt(3.0) / 8.0,         # 0.6802
    "phi_HCP":    math.pi / (3.0 * math.sqrt(2.0)),       # 0.7405
    "phi_SC":     math.pi / 6.0,                          # 0.5236
    "phi_diamond":math.pi * math.sqrt(3.0) / 16.0,        # 0.3401
    # Mass densities (g/cm^3) — for sanity checking
    "rho_NaCl":   2.165,
    "rho_diamond":3.515,
    "rho_graphite": 2.267,
}


# ---------------------------------------------------------------------------
# 1. CrystalGeometry — explicit lattice constructions
# ---------------------------------------------------------------------------

@dataclass
class CrystalGeometry:
    """Explicit construction of crystal unit cells from K_4 substrate tiling.

    The 14 Bravais lattices are exposed via :meth:`bravais_basis`;
    composite ionic / covalent / layered structures via dedicated
    constructors that return `(positions, species, cell_vectors)` tuples
    suitable for plotting and structure-factor evaluation.

    Coordinates are in fractional units of the cubic edge `a` unless
    indicated as Cartesian (Å).
    """

    # --- Bravais lattice catalogue -----------------------------------------
    @staticmethod
    def bravais_basis(name: str, a: float = 1.0,
                      b: float | None = None,
                      c: float | None = None,
                      alpha: float = 90.0,
                      beta: float = 90.0,
                      gamma: float = 90.0) -> Tuple[NDArray, NDArray]:
        """Return `(lattice_vectors, motif)` for one of the 14 Bravais lattices.

        ``name`` is one of:
            'cP', 'cI', 'cF',                # cubic: P, I (BCC), F (FCC)
            'tP', 'tI',                      # tetragonal P, I
            'oP', 'oI', 'oF', 'oC',          # orthorhombic P, I, F, C
            'mP', 'mC',                      # monoclinic P, C
            'aP',                            # triclinic
            'hP', 'hR'                       # hexagonal P, rhombohedral
        Returns
        -------
        lattice_vectors : (3, 3) array     Conventional cell row-vectors (Å).
        motif : (n, 3) array               Fractional motif positions.
        """
        b = a if b is None else b
        c = a if c is None else c
        ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
        sg = math.sin(math.radians(gamma))
        # General triclinic vectors (Niggli-style)
        a1 = np.array([a, 0.0, 0.0])
        a2 = np.array([b * cg, b * sg, 0.0])
        cx = c * cb
        cy = c * (ca - cb * cg) / max(sg, 1e-12)
        cz = c * math.sqrt(max(0.0, 1.0 - cb**2 - ((ca - cb*cg)/max(sg,1e-12))**2))
        a3 = np.array([cx, cy, cz])
        L = np.vstack([a1, a2, a3])

        if name == "cP" or name == "tP" or name == "oP" or name == "mP" or name == "aP" or name == "hP":
            motif = np.array([[0.0, 0.0, 0.0]])
        elif name == "cI" or name == "tI" or name == "oI":
            motif = np.array([[0.0, 0.0, 0.0],
                              [0.5, 0.5, 0.5]])
        elif name == "cF" or name == "oF":
            motif = np.array([[0.0, 0.0, 0.0],
                              [0.5, 0.5, 0.0],
                              [0.5, 0.0, 0.5],
                              [0.0, 0.5, 0.5]])
        elif name == "oC" or name == "mC":
            motif = np.array([[0.0, 0.0, 0.0],
                              [0.5, 0.5, 0.0]])
        elif name == "hR":
            motif = np.array([[0.0, 0.0, 0.0],
                              [1.0/3.0, 2.0/3.0, 1.0/3.0],
                              [2.0/3.0, 1.0/3.0, 2.0/3.0]])
        else:
            raise ValueError(f"unknown Bravais label '{name}'")
        return L, motif

    # --- 14 lattices (canonical names for catalogue tests) -----------------
    BRAVAIS_NAMES: List[str] = field(
        default_factory=lambda: [
            "cP", "cI", "cF",
            "tP", "tI",
            "oP", "oI", "oF", "oC",
            "mP", "mC",
            "aP",
            "hP", "hR",
        ]
    )

    # --- NaCl (rocksalt, FCC) ----------------------------------------------
    @staticmethod
    def nacl(a: float = REFERENCE["a_NaCl"]) -> Dict[str, NDArray]:
        """Rocksalt structure: two interpenetrating FCC sublattices.

        Na+ at FCC corners + face centres; Cl- shifted by (1/2, 0, 0).
        """
        fcc = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ])
        na = fcc * a
        cl = (fcc + np.array([0.5, 0.0, 0.0])) * a
        cl = np.mod(cl, a)
        return {
            "Na": na,
            "Cl": cl,
            "cell": np.eye(3) * a,
            "a": a,
            "z": 4,                     # formula units per cell
        }

    # --- CsCl (primitive cubic, two sublattices) ---------------------------
    @staticmethod
    def cscl(a: float = REFERENCE["a_CsCl"]) -> Dict[str, NDArray]:
        return {
            "Cs": np.array([[0.0, 0.0, 0.0]]) * a,
            "Cl": np.array([[0.5, 0.5, 0.5]]) * a,
            "cell": np.eye(3) * a,
            "a": a,
            "z": 1,
        }

    # --- ZnS (zincblende, FCC + tet shift) ---------------------------------
    @staticmethod
    def zns(a: float = REFERENCE["a_ZnS"]) -> Dict[str, NDArray]:
        fcc = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ])
        zn = fcc * a
        s = (fcc + np.array([0.25, 0.25, 0.25])) * a
        return {
            "Zn": zn,
            "S": s,
            "cell": np.eye(3) * a,
            "a": a,
            "z": 4,
        }

    # --- Diamond (FCC + tet shift, like ZnS but identical species) ---------
    @staticmethod
    def diamond(a: float = REFERENCE["a_diamond"]) -> Dict[str, NDArray]:
        fcc = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ])
        c1 = fcc * a
        c2 = (fcc + np.array([0.25, 0.25, 0.25])) * a
        positions = np.vstack([c1, c2])
        return {
            "C": positions,
            "cell": np.eye(3) * a,
            "a": a,
            "z": 8,
        }

    # --- Graphite (hexagonal, AB-stacked) ----------------------------------
    @staticmethod
    def graphite(a: float = REFERENCE["a_graphite_a"],
                 c: float = REFERENCE["a_graphite_c"]) -> Dict[str, NDArray]:
        """AB-stacked hexagonal graphite, 4 atoms / cell."""
        # Hexagonal cell vectors
        a1 = np.array([a, 0.0, 0.0])
        a2 = np.array([-0.5 * a, 0.5 * math.sqrt(3.0) * a, 0.0])
        a3 = np.array([0.0, 0.0, c])
        cell = np.vstack([a1, a2, a3])
        # A layer: positions (0, 0, 0) and (1/3, 2/3, 0)
        # B layer: shifted by c/2 and (2/3, 1/3, 1/2)
        frac = np.array([
            [0.0,     0.0,     0.0],
            [1.0/3.0, 2.0/3.0, 0.0],
            [0.0,     0.0,     0.5],
            [2.0/3.0, 1.0/3.0, 0.5],
        ])
        cart = frac @ cell
        return {
            "C": cart,
            "cell": cell,
            "a": a,
            "c": c,
            "z": 4,
        }

    # --- Perovskite ABX_3 (e.g. SrTiO_3 cubic) -----------------------------
    @staticmethod
    def perovskite(a: float = REFERENCE["a_SrTiO3"],
                   A: str = "Sr", B: str = "Ti", X: str = "O"
                   ) -> Dict[str, NDArray]:
        """Cubic perovskite ABX_3: A at corners, B at body centre, X on faces."""
        return {
            A: np.array([[0.0, 0.0, 0.0]]) * a,
            B: np.array([[0.5, 0.5, 0.5]]) * a,
            X: np.array([
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
            ]) * a,
            "cell": np.eye(3) * a,
            "a": a,
            "z": 1,
        }


# ---------------------------------------------------------------------------
# 2. CrystalSimulator — physical observables
# ---------------------------------------------------------------------------

@dataclass
class CrystalSimulator:
    """Physical observables for crystal structures built from K_4 substrate cells."""

    xi_m: float = XI_M_DEFAULT
    results: Dict[str, float] = field(default_factory=dict)

    # -----------------------------------------------------------------------
    # Lattice constant scaling from substrate ξ
    # -----------------------------------------------------------------------
    def lattice_constant_from_xi(self, crystal: str) -> float:
        """Return the lattice constant for a tabulated crystal in metres,
        treating the substrate ξ as the underlying cell length unit.

        Macroscopic lattice constants are tabulated XRD measurements; the
        substrate-side claim is that they all scale linearly with ξ
        (changing ξ rescales every length identically).  The dimensionless
        ratios a_crystal / a_NaCl are *species-determined* and are pulled
        from XRD tables; this method returns the absolute length as
        ``ratio · a_NaCl_ref`` and registers a ``(xi/xi_default)``
        sensitivity factor.
        """
        ratios = {
            "NaCl":     1.0,
            "CsCl":     REFERENCE["a_CsCl"]    / REFERENCE["a_NaCl"],
            "ZnS":      REFERENCE["a_ZnS"]     / REFERENCE["a_NaCl"],
            "diamond":  REFERENCE["a_diamond"] / REFERENCE["a_NaCl"],
            "graphite": REFERENCE["a_graphite_a"] / REFERENCE["a_NaCl"],
            "SrTiO3":   REFERENCE["a_SrTiO3"]  / REFERENCE["a_NaCl"],
        }
        if crystal not in ratios:
            raise ValueError(f"unknown crystal '{crystal}'")
        a_ref_m = REFERENCE["a_NaCl"] * 1e-10                  # Å → m
        scale = self.xi_m / XI_M_DEFAULT
        a_m = ratios[crystal] * a_ref_m * scale
        self.results[f"a_{crystal}_m"] = a_m
        return a_m

    # -----------------------------------------------------------------------
    # Bragg X-ray diffraction
    # -----------------------------------------------------------------------
    @staticmethod
    def bragg_two_theta_deg(d_spacing_ang: float,
                            wavelength_ang: float = CU_KALPHA_ANG,
                            n: int = 1) -> float:
        """Bragg's law: 2θ for the n-th order reflection from a plane of
        spacing ``d_spacing_ang`` (Å) at wavelength ``wavelength_ang`` (Å).

        Returns 2θ in degrees, or ``nan`` when sin θ > 1 (peak forbidden).
        """
        sin_theta = n * wavelength_ang / (2.0 * d_spacing_ang)
        if sin_theta > 1.0:
            return float("nan")
        return 2.0 * math.degrees(math.asin(sin_theta))

    @staticmethod
    def cubic_d_spacing(a_ang: float, h: int, k: int, l: int) -> float:
        """d_{hkl} for a cubic lattice."""
        s2 = h * h + k * k + l * l
        if s2 == 0:
            return float("inf")
        return a_ang / math.sqrt(s2)

    @staticmethod
    def fcc_allowed(h: int, k: int, l: int) -> bool:
        """Selection rule: FCC reflections require h, k, l all-even or all-odd."""
        all_even = (h % 2 == 0 and k % 2 == 0 and l % 2 == 0)
        all_odd = (h % 2 == 1 and k % 2 == 1 and l % 2 == 1)
        return all_even or all_odd

    @staticmethod
    def diamond_allowed(h: int, k: int, l: int) -> bool:
        """Diamond selection rule: FCC plus h+k+l ≠ 4n+2 (when all even)."""
        if not CrystalSimulator.fcc_allowed(h, k, l):
            return False
        if (h % 2 == 0):
            # all even: h+k+l must be 4n
            return ((h + k + l) % 4) == 0
        # all odd is always allowed
        return True

    def xrd_peaks(self,
                  crystal: str,
                  wavelength_ang: float = CU_KALPHA_ANG,
                  hkl_max: int = 5) -> List[Dict[str, float]]:
        """Enumerate Bragg peaks for a tabulated crystal.

        Returns a list of {hkl, d, 2theta, intensity} dicts sorted by 2θ.
        Intensities are crude squared-magnitude structure factors with
        atomic-form-factor approximated by atomic number Z.
        """
        peaks: List[Dict[str, float]] = []

        if crystal == "NaCl":
            geom = CrystalGeometry.nacl()
            atoms = {"Na": 11, "Cl": 17}
            allowed = CrystalSimulator.fcc_allowed
        elif crystal == "diamond":
            geom = CrystalGeometry.diamond()
            atoms = {"C": 6}
            allowed = CrystalSimulator.diamond_allowed
        elif crystal == "ZnS":
            geom = CrystalGeometry.zns()
            atoms = {"Zn": 30, "S": 16}
            allowed = CrystalSimulator.fcc_allowed
        elif crystal == "CsCl":
            geom = CrystalGeometry.cscl()
            atoms = {"Cs": 55, "Cl": 17}
            allowed = lambda h, k, l: True            # primitive cubic
        else:
            raise ValueError(f"xrd_peaks: unsupported crystal '{crystal}'")

        a_ang = geom["a"]
        seen = {}
        for h in range(0, hkl_max + 1):
            for k in range(0, hkl_max + 1):
                for l in range(0, hkl_max + 1):
                    if (h, k, l) == (0, 0, 0):
                        continue
                    if not allowed(h, k, l):
                        continue
                    s2 = h * h + k * k + l * l
                    if s2 in seen:
                        continue
                    d = CrystalSimulator.cubic_d_spacing(a_ang, h, k, l)
                    two_theta = CrystalSimulator.bragg_two_theta_deg(
                        d, wavelength_ang
                    )
                    if math.isnan(two_theta):
                        continue
                    # crude structure factor F = sum_j f_j exp(2pi i hkl . r_j)
                    F_re, F_im = 0.0, 0.0
                    for sp, Z in atoms.items():
                        for r in geom[sp]:
                            phase = 2.0 * math.pi * (
                                h * r[0] + k * r[1] + l * r[2]
                            ) / a_ang
                            F_re += Z * math.cos(phase)
                            F_im += Z * math.sin(phase)
                    intensity = F_re * F_re + F_im * F_im
                    if intensity < 1e-3:
                        continue
                    seen[s2] = True
                    peaks.append({
                        "hkl": (h, k, l),
                        "s2": s2,
                        "d_ang": d,
                        "two_theta_deg": two_theta,
                        "intensity": intensity,
                    })
        peaks.sort(key=lambda p: p["two_theta_deg"])
        # Normalize intensities to peak max = 100
        if peaks:
            I_max = max(p["intensity"] for p in peaks)
            for p in peaks:
                p["intensity"] = 100.0 * p["intensity"] / I_max
        self.results[f"xrd_{crystal}_n_peaks"] = float(len(peaks))
        return peaks

    # -----------------------------------------------------------------------
    # Madelung constant — convergent direct sum for NaCl rocksalt
    # -----------------------------------------------------------------------
    @staticmethod
    def madelung_nacl(N_max: int = 60) -> float:
        """Madelung constant for rocksalt by Evjen-style cubic shells.

        Uses Evjen's trick of weighting boundary ions on each cubic shell
        by 1/2 (face), 1/4 (edge), 1/8 (corner), which gives rapid
        convergence to the literature value 1.747565 within ~1%.

        ``N_max`` sets the cubic radius (number of unit cells from origin).
        """
        M = 0.0
        for i in range(-N_max, N_max + 1):
            for j in range(-N_max, N_max + 1):
                for k in range(-N_max, N_max + 1):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    r2 = i * i + j * j + k * k
                    sign = (-1.0) ** (i + j + k)
                    # Evjen weights for shell membership at boundary
                    w = 1.0
                    if abs(i) == N_max:
                        w *= 0.5
                    if abs(j) == N_max:
                        w *= 0.5
                    if abs(k) == N_max:
                        w *= 0.5
                    M += w * sign / math.sqrt(r2)
        return -M

    @staticmethod
    def madelung_cscl(N_real: int = 8, N_recip: int = 8,
                      eta: float | None = None) -> float:
        """Madelung constant for CsCl by Ewald summation.

        Test site: Cs+ at (0,0,0); Cl- at body centre (½,½,½).
        Ewald split parameter ``eta`` defaults to 2/√π for fast
        convergence; ``N_real`` and ``N_recip`` are real- and
        reciprocal-space lattice cutoffs.

        Returns the dimensionless Madelung constant α (per ion) defined by

            V(Cs⁺) = -α · e / (4πε₀ r_NN),   r_NN(CsCl) = (√3/2)·a.

        Converges to literature 1.762675 (Tosi 1964) for N_real=N_recip=8.
        Direct lattice sums (Evjen) converge poorly for CsCl because
        the natural cubelet decomposition does not align with the
        body-centre / corner sublattice structure; Ewald is reliable.
        """
        if eta is None:
            eta = 2.0 / math.sqrt(math.pi)
        # Lattice vectors (unit cubic, a=1)
        a_vec = [
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ]
        V_cell = 1.0
        b_vec = [v * 2.0 * math.pi for v in a_vec]   # reciprocal vectors
        ions = [
            (np.array([0.0, 0.0, 0.0]), +1.0),
            (np.array([0.5, 0.5, 0.5]), -1.0),
        ]

        # Real-space sum: erfc-screened
        V_R = 0.0
        for i in range(-N_real, N_real + 1):
            for j in range(-N_real, N_real + 1):
                for k in range(-N_real, N_real + 1):
                    R = i * a_vec[0] + j * a_vec[1] + k * a_vec[2]
                    for r_ion, q in ions:
                        # Skip self
                        if (i == 0 and j == 0 and k == 0
                                and np.allclose(r_ion, 0.0)):
                            continue
                        d = float(np.linalg.norm(R + r_ion))
                        if d > 0.0:
                            V_R += q * math.erfc(eta * d) / d

        # Reciprocal-space sum
        V_K = 0.0
        for h in range(-N_recip, N_recip + 1):
            for kk in range(-N_recip, N_recip + 1):
                for ll in range(-N_recip, N_recip + 1):
                    if h == 0 and kk == 0 and ll == 0:
                        continue
                    G = h * b_vec[0] + kk * b_vec[1] + ll * b_vec[2]
                    G2 = float(np.dot(G, G))
                    S_re = sum(q * math.cos(float(np.dot(G, r_ion)))
                               for r_ion, q in ions)
                    V_K += (4.0 * math.pi / V_cell) * (
                        math.exp(-G2 / (4.0 * eta * eta)) / G2
                    ) * S_re

        # Self-energy correction at test site (charge +1)
        V_self = -2.0 * eta / math.sqrt(math.pi)

        V_total = V_R + V_K + V_self
        return -V_total * (math.sqrt(3.0) / 2.0)

    # -----------------------------------------------------------------------
    # Atomic packing fractions
    # -----------------------------------------------------------------------
    @staticmethod
    def packing_fraction_fcc() -> float:
        """FCC atomic packing fraction: π / (3√2)."""
        return math.pi / (3.0 * math.sqrt(2.0))

    @staticmethod
    def packing_fraction_bcc() -> float:
        """BCC atomic packing fraction: π√3 / 8."""
        return math.pi * math.sqrt(3.0) / 8.0

    @staticmethod
    def packing_fraction_sc() -> float:
        """Simple cubic atomic packing fraction: π / 6."""
        return math.pi / 6.0

    @staticmethod
    def packing_fraction_diamond() -> float:
        """Diamond cubic packing fraction: π√3 / 16."""
        return math.pi * math.sqrt(3.0) / 16.0

    @staticmethod
    def packing_fraction_hcp() -> float:
        """HCP packing fraction: π / (3√2) — same as FCC (close-packed)."""
        return math.pi / (3.0 * math.sqrt(2.0))

    # -----------------------------------------------------------------------
    # Mass density (g/cm^3)
    # -----------------------------------------------------------------------
    @staticmethod
    def mass_density(a_ang: float, z_per_cell: int,
                     formula_mass_amu: float) -> float:
        """ρ = z·M / (N_A · V_cell), V from cubic edge a in Å (returns g/cm³)."""
        V_cm3 = (a_ang * 1e-8) ** 3
        return z_per_cell * formula_mass_amu / (N_AVOGADRO * V_cm3)

    def density_NaCl(self) -> float:
        rho = self.mass_density(REFERENCE["a_NaCl"], 4, 22.99 + 35.45)
        self.results["rho_NaCl_gcm3"] = rho
        return rho

    def density_diamond(self) -> float:
        rho = self.mass_density(REFERENCE["a_diamond"], 8, 12.011)
        self.results["rho_diamond_gcm3"] = rho
        return rho

    # -----------------------------------------------------------------------
    # Convenience: full report for a named crystal
    # -----------------------------------------------------------------------
    def report(self, crystal: str) -> Dict[str, float]:
        out: Dict[str, float] = {}
        out["a_metres"] = self.lattice_constant_from_xi(crystal)
        out["a_ang"] = out["a_metres"] * 1e10
        if crystal == "NaCl":
            out["madelung"] = CrystalSimulator.madelung_nacl(N_max=12)
            out["packing"] = CrystalSimulator.packing_fraction_fcc()
            out["density_gcm3"] = self.density_NaCl()
        elif crystal == "diamond":
            out["packing"] = CrystalSimulator.packing_fraction_diamond()
            out["density_gcm3"] = self.density_diamond()
        elif crystal == "CsCl":
            out["madelung"] = CrystalSimulator.madelung_cscl(
                N_real=8, N_recip=8
            )
            out["packing"] = CrystalSimulator.packing_fraction_bcc()
        return out


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    sim = CrystalSimulator()
    print("# Crystal substrate sanity checks")
    print(f"FCC packing  : {sim.packing_fraction_fcc():.6f}  "
          f"(ref π/(3√2) = {REFERENCE['phi_FCC']:.6f})")
    print(f"BCC packing  : {sim.packing_fraction_bcc():.6f}")
    print(f"SC  packing  : {sim.packing_fraction_sc():.6f}")
    print(f"diamond pack : {sim.packing_fraction_diamond():.6f}")
    M = sim.madelung_nacl(N_max=15)
    print(f"NaCl Madelung (Evjen N=15): {M:.6f}  ref {REFERENCE['M_NaCl']:.6f}")
    M_cs = sim.madelung_cscl(N_real=8, N_recip=8)
    print(f"CsCl Madelung (Ewald)   : {M_cs:.6f}  ref {REFERENCE['M_CsCl']:.6f}")
    peaks_nacl = sim.xrd_peaks("NaCl")
    print(f"\nNaCl XRD peaks (Cu Kα): first 5")
    for p in peaks_nacl[:5]:
        print(f"  {p['hkl']}  2θ={p['two_theta_deg']:6.2f}°  "
              f"d={p['d_ang']:.3f} Å  I={p['intensity']:5.1f}")
    peaks_d = sim.xrd_peaks("diamond")
    print(f"\nDiamond XRD peaks: first 5")
    for p in peaks_d[:5]:
        print(f"  {p['hkl']}  2θ={p['two_theta_deg']:6.2f}°  "
              f"d={p['d_ang']:.3f} Å  I={p['intensity']:5.1f}")
    print(f"\nNaCl mass density: {sim.density_NaCl():.4f} g/cm³  "
          f"(ref {REFERENCE['rho_NaCl']})")
    print(f"Diamond density  : {sim.density_diamond():.4f} g/cm³  "
          f"(ref {REFERENCE['rho_diamond']})")
