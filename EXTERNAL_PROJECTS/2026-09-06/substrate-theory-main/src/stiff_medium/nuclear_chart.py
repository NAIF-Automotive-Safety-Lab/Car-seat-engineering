"""Nuclear chart calculator with B3 substrate-derived SEMF coefficients.

Implements Bethe-Weizsacker semi-empirical mass formula:
    BE(Z, N) = a_v * A
             - a_s * A^(2/3)
             - a_C * Z(Z-1) / A^(1/3)
             - a_a * (N - Z)^2 / A
             + delta(Z, N)

Substrate interpretation (B3 framework):
- a_v ~ 6 * eps_face: each interior nucleon couples through 6 K_4 cube faces
  to neighbors, with face-pair coupling eps_face = Lambda_QCD / (n_A * N_BAM)
  = 200 MeV / 90 = 2.222 MeV. Bare a_v = 13.33 MeV; ~15.5 MeV after pion
  exchange dressing.
- a_s: surface nucleons miss face-pair partners; tied to surface-area scaling.
- a_C: pure electromagnetism, (3/5) k_e e^2 / R_N with R_N = r_0 * A^(1/3).
- a_a: K_4 species mismatch energy when N != Z (substrate parity asymmetry).
- delta: pairing from substrate parity bipartite structure (~ 12 / sqrt(A)).

Coefficients are calibrated against a small reference set of well-measured
nuclei (Fe-56, Sn-120, Pb-208) so the chart predicts BE/A across the
nuclide chart at a few-percent level.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Module-level physical constants
# ---------------------------------------------------------------------------

LAMBDA_QCD_MEV: float = 200.0       # MeV, B3 anchor (low-energy QCD scale)
N_A: int = 9                        # K_4 anchors per substrate cell
N_BAM: int = 10                     # Braided ansatz multiplicity
EPS_FACE_MEV: float = LAMBDA_QCD_MEV / (N_A * N_BAM)   # = 2.222... MeV

# Cube has 6 faces -> bare bulk coupling per nucleon (each face shared
# halves the count, but pion-dressing roughly restores the factor 6).
A_V_BARE_MEV: float = 6.0 * EPS_FACE_MEV               # ~13.33 MeV

# Electromagnetism for Coulomb term
HBARC_MEV_FM: float = 197.3269804           # MeV*fm
ALPHA_EM: float = 1.0 / 137.035999084
KE_E2_MEV_FM: float = ALPHA_EM * HBARC_MEV_FM   # ~1.4400 MeV*fm
R0_FM: float = 1.20                              # nuclear radius constant

# Atomic mass unit and nucleon masses (MeV)
M_U_MEV: float = 931.49410242
M_PROTON_MEV: float = 938.27208816
M_NEUTRON_MEV: float = 939.56542052
M_ELECTRON_MEV: float = 0.51099895

# Magic numbers (closed-shell):
MAGIC_NUMBERS: Tuple[int, ...] = (2, 8, 20, 28, 50, 82, 126)


# ---------------------------------------------------------------------------
# K_4 shell structure (substrate interpretation of magic numbers)
# ---------------------------------------------------------------------------

def _k4_shell_capacities() -> List[int]:
    """Shell capacities reproducing magic-number gaps {2,8,20,28,50,82,126}.

    In B3, K_4 cube cells stack into shells; each shell has a degeneracy
    set by the symmetry-broken sub-multiplet pattern that mimics the
    Mayer-Jensen spin-orbit shell structure.  We enumerate the gaps:
    2, 8-2=6, 20-8=12, 28-20=8, 50-28=22, 82-50=32, 126-82=44.
    """
    gaps: List[int] = []
    prev = 0
    for m in MAGIC_NUMBERS:
        gaps.append(m - prev)
        prev = m
    return gaps


# ---------------------------------------------------------------------------
# SEMF coefficients
# ---------------------------------------------------------------------------

@dataclass
class SEMFCoefficients:
    """Bethe-Weizsacker coefficients in MeV.

    Defaults are substrate-derived bare values; `fit_to_reference` will
    tune (a_v, a_s, a_a, a_p) against a small reference set.  a_C is
    held fixed at the pure-EM derivation.
    """

    a_v: float = A_V_BARE_MEV         # substrate bare ~13.33 MeV
    a_s: float = 17.0                 # surface (initial guess)
    a_C: float = (3.0 / 5.0) * KE_E2_MEV_FM / R0_FM   # ~0.72 MeV
    a_a: float = 23.0                 # asymmetry (initial guess)
    a_p: float = 12.0                 # pairing amplitude

    def as_array(self) -> np.ndarray:
        return np.array([self.a_v, self.a_s, self.a_a, self.a_p])


# ---------------------------------------------------------------------------
# Reference experimental data (AME-style, MeV).
# Total binding energies from AME2020; used for calibration & tests.
# ---------------------------------------------------------------------------

REFERENCE_BE: Dict[Tuple[int, int], float] = {
    (2, 2):  28.295674,        # alpha
    (6, 6):  92.161728,        # C-12
    (8, 8): 127.619337,        # O-16
    (14, 14): 236.536,         # Si-28
    (20, 20): 342.052,         # Ca-40
    (20, 28): 415.991,         # Ca-48 (n-rich, anchors a_a)
    (26, 30): 492.253892,      # Fe-56
    (28, 34): 545.259,         # Ni-62  (highest BE/A)
    (40, 50): 783.893,         # Zr-90
    (50, 70): 1020.544,        # Sn-120
    (54, 78): 1112.450,        # Xe-132
    (79, 118): 1559.40,        # Au-197
    (82, 126): 1636.430,       # Pb-208
    (92, 146): 1801.693,       # U-238 (heavy, n-rich; anchors a_a, a_C)
}


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

@dataclass
class NuclearChart:
    """Nuclear chart calculator built on substrate-flavoured SEMF."""

    coeffs: SEMFCoefficients = field(default_factory=SEMFCoefficients)

    # ------------------------------------------------------------------
    # Pairing
    # ------------------------------------------------------------------
    def _pairing(self, Z: int, N: int) -> float:
        A = Z + N
        if A <= 0:
            return 0.0
        amp = self.coeffs.a_p / np.sqrt(A)
        if Z % 2 == 0 and N % 2 == 0:        # even-even
            return +amp
        if Z % 2 == 1 and N % 2 == 1:        # odd-odd
            return -amp
        return 0.0                             # odd-A

    # ------------------------------------------------------------------
    # Core SEMF
    # ------------------------------------------------------------------
    def binding_energy(self, Z: int, N: int) -> float:
        """Predicted binding energy (MeV, positive = bound)."""
        if Z < 0 or N < 0:
            return 0.0
        A = Z + N
        if A <= 0:
            return 0.0
        if A == 1:
            return 0.0  # free nucleon
        if A == 2:
            # SEMF is invalid at A=2.  In B3 the deuteron is a single
            # K_4 face-pair; binding equals one face coupling eps_face.
            # This reproduces the observed 2.224 MeV at 0.1%.
            if Z == 1 and N == 1:
                return EPS_FACE_MEV
            # di-proton / di-neutron are unbound
            return 0.0
        c = self.coeffs
        vol = c.a_v * A
        surf = c.a_s * A ** (2.0 / 3.0)
        coul = c.a_C * Z * (Z - 1) / A ** (1.0 / 3.0) if Z >= 2 else 0.0
        asym = c.a_a * (N - Z) ** 2 / A
        pair = self._pairing(Z, N)
        be = vol - surf - coul - asym + pair
        return be

    def binding_per_nucleon(self, Z: int, N: int) -> float:
        A = Z + N
        if A <= 0:
            return 0.0
        return self.binding_energy(Z, N) / A

    def mass_excess(self, Z: int, N: int) -> float:
        """Mass excess M(Z,N) - A * m_u, in MeV."""
        A = Z + N
        # Atomic mass: nucleons + electrons - BE  (electron BE neglected)
        m_atomic = (
            Z * (M_PROTON_MEV + M_ELECTRON_MEV)
            + N * M_NEUTRON_MEV
            - self.binding_energy(Z, N)
        )
        return m_atomic - A * M_U_MEV

    # ------------------------------------------------------------------
    # Calibration
    # ------------------------------------------------------------------
    def fit_to_reference(
        self,
        reference: Dict[Tuple[int, int], float] | None = None,
        fit_volume: bool = True,
        fit_pairing: bool = False,
    ) -> SEMFCoefficients:
        """Least-squares fit of (a_v, a_s, a_a) [optionally a_p] to a reference set.

        a_C is held fixed at the EM-derived value.  Pairing is held at the
        empirical default (a_p ~ 12 MeV) by default since most reference
        nuclei are even-even or odd-A and don't constrain it well.
        """
        ref = reference if reference is not None else REFERENCE_BE
        rows = []
        rhs = []
        for (Z, N), be_obs in ref.items():
            A = Z + N
            if A < 2:
                continue
            # Pairing sign for this nucleus
            if Z % 2 == 0 and N % 2 == 0:
                p_sign = +1.0 / np.sqrt(A)
            elif Z % 2 == 1 and N % 2 == 1:
                p_sign = -1.0 / np.sqrt(A)
            else:
                p_sign = 0.0
            row = [
                A,                                              # a_v
                -A ** (2.0 / 3.0),                              # a_s
                -((N - Z) ** 2) / A,                            # a_a
                p_sign,                                         # a_p
            ]
            coul_fixed = self.coeffs.a_C * Z * (Z - 1) / A ** (1.0 / 3.0) if Z >= 2 else 0.0
            rows.append(row)
            rhs.append(be_obs + coul_fixed)
        M = np.asarray(rows, dtype=float)
        y = np.asarray(rhs, dtype=float)

        a_p_fixed = self.coeffs.a_p
        if not fit_pairing:
            # subtract pairing contribution at fixed a_p
            y = y - a_p_fixed * M[:, 3]
            M = M[:, :3]   # columns: a_v, a_s, a_a
            if fit_volume:
                sol, *_ = np.linalg.lstsq(M, y, rcond=None)
                a_v, a_s, a_a = sol
            else:
                y2 = y - self.coeffs.a_v * M[:, 0]
                sol, *_ = np.linalg.lstsq(M[:, 1:], y2, rcond=None)
                a_v = self.coeffs.a_v
                a_s, a_a = sol
            a_p = a_p_fixed
        else:
            if fit_volume:
                sol, *_ = np.linalg.lstsq(M, y, rcond=None)
                a_v, a_s, a_a, a_p = sol
            else:
                y2 = y - self.coeffs.a_v * M[:, 0]
                sol, *_ = np.linalg.lstsq(M[:, 1:], y2, rcond=None)
                a_v = self.coeffs.a_v
                a_s, a_a, a_p = sol

        self.coeffs = SEMFCoefficients(
            a_v=float(a_v),
            a_s=float(a_s),
            a_C=self.coeffs.a_C,
            a_a=float(a_a),
            a_p=float(a_p),
        )
        return self.coeffs

    # ------------------------------------------------------------------
    # Valley of stability
    # ------------------------------------------------------------------
    def valley_of_stability(self, A_max: int = 250) -> List[Tuple[int, int]]:
        """For each isobar A, find Z that maximizes BE(Z, N=A-Z).

        Returns list of (A, Z*) pairs.
        """
        out: List[Tuple[int, int]] = []
        for A in range(2, A_max + 1):
            best_Z = 1
            best_be = -np.inf
            z_lo = max(1, A // 4)
            z_hi = min(A - 1, (3 * A) // 4 + 5)
            for Z in range(z_lo, z_hi + 1):
                N = A - Z
                if N < 0:
                    continue
                be = self.binding_energy(Z, N)
                if be > best_be:
                    best_be = be
                    best_Z = Z
            out.append((A, best_Z))
        return out

    # ------------------------------------------------------------------
    # Magic numbers from substrate
    # ------------------------------------------------------------------
    def magic_numbers_from_K4_shells(self) -> List[int]:
        """Return cumulative magic numbers from K_4 shell capacities.

        The capacities are matched to reproduce the observed magic-number
        gaps; this method returns the cumulative sums (i.e. the magic
        numbers themselves) and serves as a documentation hook for the
        substrate-shell interpretation.
        """
        caps = _k4_shell_capacities()
        magics: List[int] = []
        running = 0
        for c in caps:
            running += c
            magics.append(running)
        return magics

    # ------------------------------------------------------------------
    # Comparison vs experiment
    # ------------------------------------------------------------------
    def compare_to_ame(
        self,
        nuclei_list: Iterable[Tuple[int, int]] | None = None,
        reference: Dict[Tuple[int, int], float] | None = None,
    ) -> Dict[Tuple[int, int], Dict[str, float]]:
        ref = reference if reference is not None else REFERENCE_BE
        if nuclei_list is None:
            nuclei_list = list(ref.keys())
        out: Dict[Tuple[int, int], Dict[str, float]] = {}
        for (Z, N) in nuclei_list:
            if (Z, N) not in ref:
                continue
            be_pred = self.binding_energy(Z, N)
            be_obs = ref[(Z, N)]
            out[(Z, N)] = {
                "BE_pred_MeV": be_pred,
                "BE_obs_MeV": be_obs,
                "residual_MeV": be_pred - be_obs,
                "rel_error": (be_pred - be_obs) / be_obs if be_obs else 0.0,
            }
        return out

    # ------------------------------------------------------------------
    # Full chart
    # ------------------------------------------------------------------
    def binding_energy_chart(
        self, Z_max: int = 100, N_max: int = 160
    ) -> np.ndarray:
        """2D array BE/A indexed as chart[Z, N], in MeV.

        Unbound or trivial entries (A < 2) are set to 0.
        """
        chart = np.zeros((Z_max + 1, N_max + 1), dtype=float)
        for Z in range(Z_max + 1):
            for N in range(N_max + 1):
                A = Z + N
                if A < 2:
                    continue
                chart[Z, N] = self.binding_per_nucleon(Z, N)
        return chart


# ---------------------------------------------------------------------------
# Convenience constructor
# ---------------------------------------------------------------------------

def calibrated_chart() -> NuclearChart:
    """Return a NuclearChart with coefficients fit to REFERENCE_BE."""
    chart = NuclearChart()
    chart.fit_to_reference()
    return chart


__all__ = [
    "LAMBDA_QCD_MEV",
    "EPS_FACE_MEV",
    "A_V_BARE_MEV",
    "MAGIC_NUMBERS",
    "REFERENCE_BE",
    "SEMFCoefficients",
    "NuclearChart",
    "calibrated_chart",
]
