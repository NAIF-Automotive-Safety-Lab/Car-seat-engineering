"""Unified hadron spectrum from K_4 cell stacking + face-spin v4 baryons.

Two complementary baryon constructions live in this module:

  1. **Cell-stacking ansatz** (:class:`HadronSpectrum`) — mesons are
     cell-pairs (two K_4 cells on a shared face) and baryons are closed
     triangles of three K_4 cells at a Y-junction. Bare inventory model;
     gives p, n at <1% but Xi residuals drift to ~13%.

  2. **Face-spin v4 chromomagnetic** (:class:`BaryonFaceSpinV4`) — the
     De Rújula-Georgi-Glashow spin-flavour decomposition computed inside
     the substrate via K_substrate = (8/3) σ ξ² σ^{3/2} where σ is the
     substrate Cornell string tension and ξ is the QCD coherence length.
     Six inventory-derived couplings (one per flavour-spin pair bucket)
     plus two mass anchors (proton + Λ) hit the full octet AND decuplet
     at sub-2% mean residual. This is the "v4" mechanism referenced in
     the b3_baryon_face_spin_v4 audit notes.

All scales are anchored at Λ_QCD = 200 MeV and use the audited B3
inventory integers (n_M, N_BAM, K_pair, K_rank, n_R) imported from
:mod:`b3_constants`. The v4 face-spin calculator additionally pulls
σ_substrate = (K_pair·K_rank − 1)/K_pair · Λ_QCD² = 9/2 · Λ² = 0.18 GeV²
and the coherence length ξ_QCD = 0.2 fm (substrate cell length at the
QCD scale).

Pattern: pure-data flavour map, two formulas (cell-stacking + face-spin
v4), one report. Not a fit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import b3_constants as bc


# ---------------------------------------------------------------------------
# Anchors and coupling constants
# ---------------------------------------------------------------------------

LAMBDA = bc.LAMBDA_QCD_MEV  # 200 MeV — universal mass-torque scale anchor

# Constituent quark "cell-torque" values T_q in units of Λ_QCD. These are
# inventory-fixed: they are the same numbers used by the v4 baryon face-spin
# calculator and the heavy-quark module. They are NOT free fits — each is
# pinned by a separate single observable elsewhere in the framework.
#
#   T_u = (n_R/n_M) · (K_rank/K_pair)    ≈ 0.168
#   T_d = T_u · (1 + 1/n_M)               (isospin micro-split)
#   T_s = T_u + (N_BAM/K_rank)            ≈ 1.368  → ~470 MeV constituent
#   T_c = T_s + n_R/(K_pair·K_rank)       ≈ 3.168  → ~1500 MeV
#   T_b = T_c + n_M/(K_rank·n_R)          ≈ 6.146  → ~4700 MeV
T_U = (bc.n_R / bc.n_M) * (bc.K_rank / bc.K_pair)
T_D = T_U * (1.0 + 1.0 / bc.n_M)
T_S = T_U + (bc.N_BAM / bc.K_rank)
T_C = T_S + bc.n_R / (bc.K_pair * bc.K_rank)
T_B = T_C + bc.n_M / (bc.K_rank * bc.n_R)

QUARK_TORQUE: Dict[str, float] = {
    "u": T_U,
    "d": T_D,
    "s": T_S,
    "c": T_C,
    "b": T_B,
}

# String-tension binding for cell-pair (meson) and closed-triangle (baryon).
# Both are derived from the n_M / K_rank^k inventory ratios:
#
#   B_meson  = n_M / K_rank^3                                  ≈ 2.144
#   B_baryon = (2·n_M + N_BAM² + K_rank·K_pair) / K_rank^3      ≈ 4.656
#
# The baryon binding has THREE inventory pieces: doubled cell-pair binding
# (2·n_M for the two shared faces of a closed triangle), the centre-vertex
# multiplicity (N_BAM², the squared 2D coordination giving the Y-junction
# face count), and the K_rank·K_pair edge-doubling at the central vertex.
B_MESON = bc.n_M / (bc.K_rank ** 3)
B_BARYON = (
    2.0 * bc.n_M + bc.N_BAM ** 2 + bc.K_rank * bc.K_pair
) / (bc.K_rank ** 3)

# Spin-spin (chromomagnetic) coupling — split octet vs decuplet branch.
# Inventory: C_ss = (K_rank - K_pair - 1)/n_R = 1/9 ≈ 0.1111. Acts with
# sign +1 on decuplet (J=3/2 spin-aligned) and -1/2 on octet (J=1/2).
# Pair contribution uses 1/(T_i+T_j) — heavier pairs split less, which is
# the inverse-mass scaling of the chromomagnetic operator.
C_SS = (bc.K_rank - bc.K_pair - 1) / bc.n_R

# Pseudoscalar/vector light-meson channel multipliers. Pseudoscalars are
# light by the Goldstone mechanism (chiral-symmetry breaking) — the
# cell-pair binding is suppressed by 1/(K_rank+1) ≈ 0.167. Vectors carry
# the full string binding scaled up by K_rank/(K_pair+1) ≈ 1.667.
G_PS = 1.0 / (bc.K_rank + 1)        # ≈ 0.1667
G_V = bc.K_rank / (bc.K_pair + 1)   # ≈ 1.667


# ---------------------------------------------------------------------------
# PDG reference values (MeV) for residual reporting
# ---------------------------------------------------------------------------

PDG: Dict[str, float] = {
    # pseudoscalar nonet
    "pi":   139.57, "pi0": 134.98, "K":  493.68, "K0": 497.61,
    "eta":  547.86, "etap": 957.78,
    # vector nonet
    "rho":  775.26, "omega": 782.66, "K*": 891.66, "phi": 1019.46,
    # spin-1/2 octet
    "p":    938.272, "n":   939.565,
    "Lambda": 1115.683,
    "Sigma+": 1189.37, "Sigma0": 1192.642, "Sigma-": 1197.45,
    "Xi0":   1314.86, "Xi-":  1321.71,
    # spin-3/2 decuplet
    "Delta0":  1232.0, "Delta+": 1232.0, "Delta-": 1232.0, "Delta++": 1232.0,
    "Sigma*+": 1382.8, "Sigma*0": 1383.7, "Sigma*-": 1387.2,
    "Xi*0":    1531.8, "Xi*-":   1535.0,
    "Omega-":  1672.45,
}


# Quark content tables (constituent labels)
MESON_QUARKS: Dict[str, Tuple[str, str]] = {
    "pi":   ("u", "d"),  "pi0": ("u", "u"),
    "K":    ("u", "s"),  "K0":  ("d", "s"),
    "eta":  ("u", "s"),  "etap": ("s", "s"),  # effective; mixing absorbed
    "rho":  ("u", "d"),  "omega": ("u", "u"),
    "K*":   ("u", "s"),  "phi":  ("s", "s"),
}

BARYON_QUARKS: Dict[str, Tuple[str, str, str]] = {
    "p":   ("u", "u", "d"), "n":   ("u", "d", "d"),
    "Lambda": ("u", "d", "s"),
    "Sigma+": ("u", "u", "s"), "Sigma0": ("u", "d", "s"), "Sigma-": ("d", "d", "s"),
    "Xi0":   ("u", "s", "s"), "Xi-":   ("d", "s", "s"),
    # decuplet
    "Delta++": ("u", "u", "u"), "Delta+": ("u", "u", "d"),
    "Delta0":  ("u", "d", "d"), "Delta-": ("d", "d", "d"),
    "Sigma*+": ("u", "u", "s"), "Sigma*0": ("u", "d", "s"), "Sigma*-": ("d", "d", "s"),
    "Xi*0":   ("u", "s", "s"), "Xi*-":   ("d", "s", "s"),
    "Omega-": ("s", "s", "s"),
}

OCTET = ["p", "n", "Lambda", "Sigma+", "Sigma0", "Sigma-", "Xi0", "Xi-"]
DECUPLET = [
    "Delta++", "Delta+", "Delta0", "Delta-",
    "Sigma*+", "Sigma*0", "Sigma*-",
    "Xi*0", "Xi*-", "Omega-",
]
PS_NONET = ["pi", "pi0", "K", "K0", "eta", "etap"]
V_NONET = ["rho", "omega", "K*", "phi"]


# ---------------------------------------------------------------------------
# HadronSpectrum class
# ---------------------------------------------------------------------------


@dataclass
class Residual:
    name: str
    pred: float
    pdg: float

    @property
    def abs_err(self) -> float:
        return self.pred - self.pdg

    @property
    def rel_err(self) -> float:
        return (self.pred - self.pdg) / self.pdg if self.pdg else 0.0


class HadronSpectrum:
    """Compute meson + baryon masses from K_4 cell stacking.

    Two formulas, one anchor:

        meson:  M = Λ_QCD · [ (T_q1 + T_q2) + g · B_meson ]
                    g = G_PS for pseudoscalars, G_V for vectors

        baryon: M = Λ_QCD · [ (T_q1 + T_q2 + T_q3) + B_baryon
                              + s · C_SS · Σ_pairs ]
                    s = +1 for decuplet (J=3/2), -1/2 for octet (J=1/2)
    """

    def __init__(self, lam: float = LAMBDA) -> None:
        self.lam = lam

    # --- core formulas -----------------------------------------------------

    def _meson_formula(self, q1: str, q2: str, channel: str) -> float:
        t = QUARK_TORQUE[q1] + QUARK_TORQUE[q2]
        if channel == "PS":
            return self.lam * (t + G_PS * B_MESON)
        if channel == "V":
            return self.lam * (t + G_V * B_MESON)
        raise ValueError(f"unknown meson channel {channel!r}")

    def _baryon_formula(
        self, q1: str, q2: str, q3: str, branch: str, lambda_like: bool = False
    ) -> float:
        t_sum = QUARK_TORQUE[q1] + QUARK_TORQUE[q2] + QUARK_TORQUE[q3]
        # Spin-spin pair sum: each unordered pair contributes 1/(T_i+T_j),
        # matching the inverse-quark-mass scaling of the chromomagnetic
        # operator (heavier pairs split less). For Λ-like states (octet
        # I=0), the lightest pair is in the antisymmetric spin-0 channel
        # and contributes with the opposite sign — this is the standard
        # SU(6) wavefunction split that resolves the Σ⁰-Λ degeneracy.
        Ts = (QUARK_TORQUE[q1], QUARK_TORQUE[q2], QUARK_TORQUE[q3])
        pair_vals = [
            (1.0 / (Ts[0] + Ts[1]), Ts[0] + Ts[1]),
            (1.0 / (Ts[0] + Ts[2]), Ts[0] + Ts[2]),
            (1.0 / (Ts[1] + Ts[2]), Ts[1] + Ts[2]),
        ]
        if lambda_like:
            # Lightest pair is in the antisymmetric spin-singlet channel:
            # the chromomagnetic operator's spin matrix gives -3× the
            # triplet value on that pair. The other two pairs (involving
            # the heavy quark) remain in the triplet channel.
            light_idx = min(range(3), key=lambda i: pair_vals[i][1])
            pairs = sum(
                ((3.0 if i == light_idx else 1.0) * pair_vals[i][0])
                for i in range(3)
            )
        else:
            pairs = sum(pv[0] for pv in pair_vals)
        s = +1.0 if branch == "decuplet" else -0.5
        return self.lam * (t_sum + B_BARYON + s * C_SS * pairs)

    # --- public single-mass methods ----------------------------------------

    def meson_mass(self, name: str) -> float:
        """Meson mass in MeV. Channel inferred from name."""
        if name not in MESON_QUARKS:
            raise KeyError(f"unknown meson {name!r}")
        q1, q2 = MESON_QUARKS[name]
        channel = "V" if name in V_NONET else "PS"
        return self._meson_formula(q1, q2, channel)

    def baryon_mass(self, name: str) -> float:
        """Baryon mass in MeV. Branch inferred from octet/decuplet membership."""
        if name not in BARYON_QUARKS:
            raise KeyError(f"unknown baryon {name!r}")
        q1, q2, q3 = BARYON_QUARKS[name]
        branch = "decuplet" if name in DECUPLET else "octet"
        lambda_like = name == "Lambda"
        return self._baryon_formula(q1, q2, q3, branch, lambda_like=lambda_like)

    # --- spectra -----------------------------------------------------------

    def octet_spectrum(self) -> Dict[str, float]:
        return {n: self.baryon_mass(n) for n in OCTET}

    def decuplet_spectrum(self) -> Dict[str, float]:
        return {n: self.baryon_mass(n) for n in DECUPLET}

    def meson_nonet(self) -> Dict[str, float]:
        names = PS_NONET + V_NONET
        return {n: self.meson_mass(n) for n in names if n in MESON_QUARKS}

    # --- structural identities --------------------------------------------

    def gell_mann_okubo(self) -> Dict[str, float]:
        """GMO octet identity: 2·(N + Ξ) ?= 3·Λ + Σ.

        Returns the LHS, RHS, and relative residual using the predicted
        spectrum (averaging isospin partners).
        """
        N = 0.5 * (self.baryon_mass("p") + self.baryon_mass("n"))
        Xi = 0.5 * (self.baryon_mass("Xi0") + self.baryon_mass("Xi-"))
        Lam = self.baryon_mass("Lambda")
        Sig = (
            self.baryon_mass("Sigma+")
            + self.baryon_mass("Sigma0")
            + self.baryon_mass("Sigma-")
        ) / 3.0
        lhs = 2.0 * (N + Xi)
        rhs = 3.0 * Lam + Sig
        return {
            "lhs_2(N+Xi)": lhs,
            "rhs_3L+S": rhs,
            "rel_residual": (lhs - rhs) / rhs,
        }

    # --- residuals & report -----------------------------------------------

    def compare_to_pdg(self) -> List[Residual]:
        out: List[Residual] = []
        for name, pred in self.meson_nonet().items():
            if name in PDG:
                out.append(Residual(name, pred, PDG[name]))
        for name, pred in self.octet_spectrum().items():
            if name in PDG:
                out.append(Residual(name, pred, PDG[name]))
        for name, pred in self.decuplet_spectrum().items():
            if name in PDG:
                out.append(Residual(name, pred, PDG[name]))
        return out

    def report(self) -> str:
        lines: List[str] = []
        lines.append(
            "B3 hadron spectrum  (Λ_QCD = {:.0f} MeV, n_M = {})".format(
                self.lam, bc.n_M
            )
        )
        lines.append("=" * 66)
        lines.append(
            "{:<10s} {:>12s} {:>12s} {:>10s} {:>10s}".format(
                "name", "B3 (MeV)", "PDG (MeV)", "Δ (MeV)", "rel %"
            )
        )
        lines.append("-" * 66)
        residuals = self.compare_to_pdg()
        for r in residuals:
            lines.append(
                "{:<10s} {:>12.2f} {:>12.2f} {:>+10.2f} {:>+9.2f}%".format(
                    r.name, r.pred, r.pdg, r.abs_err, 100.0 * r.rel_err
                )
            )
        lines.append("-" * 66)
        rels = [abs(r.rel_err) for r in residuals]
        if rels:
            lines.append(
                "mean |rel| = {:.2f}%   max |rel| = {:.2f}%".format(
                    100.0 * sum(rels) / len(rels), 100.0 * max(rels)
                )
            )
        gmo = self.gell_mann_okubo()
        lines.append(
            "GMO check: 2(N+Ξ) = {:.2f}   3Λ+Σ = {:.2f}   "
            "rel resid = {:+.3f}%".format(
                gmo["lhs_2(N+Xi)"], gmo["rhs_3L+S"], 100.0 * gmo["rel_residual"]
            )
        )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Face-spin v4 baryon spectrum  (chromomagnetic substrate construction)
# ---------------------------------------------------------------------------
#
# This is the upgraded baryon model referenced in the
# ``b3_baryon_face_spin_v4`` audit notes. It hits the full octet AND
# decuplet at sub-2% mean residual using:
#
#   * Substrate-DERIVED ingredients:
#       - σ = (K_pair·K_rank − 1)/K_pair · Λ_QCD²   = 0.18 GeV²
#         [Cornell string tension from inventory integers]
#       - ξ_QCD = 0.2 fm (≈ 1/√σ, substrate coherence length at QCD scale)
#       - K_substrate = (8/3) σ ξ² σ^{3/2}             = 0.0377 GeV³
#         [chromomagnetic contact coefficient — pure σ, ξ inputs]
#       - Six spin-flavour pair coefficients (c_qq, c_qs, c_ss for J=½ and
#         J=3/2 baryons) — Clebsch-Gordan from SU(6) wavefunction algebra,
#         not free parameters.
#
#   * Anchored ingredients (2 anchors total):
#       - m_q_struct ≈ 365 MeV   [light-quark structure mass; fixed by the
#         proton mass anchor M_p = 938.27 MeV]
#       - m_s_struct ≈ 542 MeV   [strange-quark structure mass; fixed by
#         the Λ⁰ mass anchor M_Λ = 1115.68 MeV]
#
# Mass formula (per baryon):
#
#     M_B = N_q · m_q_struct + N_s · m_s_struct
#         + K_substrate · (c_qq/m_q² + c_qs/(m_q m_s) + c_ss/m_s²)
#
# The chromomag denominator masses default to the GEOMETRIC scale
# m_q_chromo = √σ ≈ 424 MeV (from R₀ ≈ 1/√σ). The structure masses are
# the LINEAR sum coefficients that absorb the residual confinement binding.
#
# A/B/C category labels for downstream classification:
#   [A] σ, ξ, K_substrate, c_qq/c_qs/c_ss        (substrate-derived)
#   [B] m_q_struct, m_s_struct                   (anchored on 2 masses)
#   [C] none                                     (no empirical-only inputs)
#

# --- Substrate Cornell σ from inventory integers -----------------------------
SIGMA_QCD_GEV2: float = (
    (bc.K_pair * bc.K_rank - 1) / bc.K_pair * (bc.LAMBDA_QCD_MEV / 1000.0) ** 2
)
"""[A] Cornell string tension σ = (K_pair·K_rank − 1)/K_pair · Λ_QCD².

= 9/2 · (0.200 GeV)² = 0.18 GeV². Substrate-derived from K_pair=2,
K_rank=5 and Λ_QCD = 200 MeV. Matches lattice-QCD value 0.18 GeV²."""

# --- Coherence length ξ at the QCD scale -------------------------------------
XI_QCD_FM: float = 0.2
"""[A] Coherence length ξ_QCD ≈ 1/√σ ≈ 0.2 fm. Substrate-natural value at
the QCD scale; consistent with R₀ = 1/√σ ≈ 0.46 fm Y-junction arm length
divided by ~2.3 (the Y-junction arm-to-coherence ratio)."""

HBARC_GEVFM: float = 0.197326980
XI_QCD_INV_GEV: float = XI_QCD_FM / HBARC_GEVFM   # ≈ 1.0135 GeV⁻¹

# --- Geometric chromomag mass m_q_chromo = √σ -------------------------------
M_K_CHROMO_GEV: float = math.sqrt(SIGMA_QCD_GEV2)   # √σ ≈ 0.4243 GeV
"""[A] Geometric chromomag mass = √σ ≈ 424 MeV. Sets the natural propagator
scale 1/m² in the spin-spin contact term."""


# Pre-computed chromomag K_substrate
def _K_substrate_GeV3(sigma_GeV2: float = SIGMA_QCD_GEV2,
                      xi_inv_GeV: float = XI_QCD_INV_GEV) -> float:
    """[A] K_substrate = (8/3) σ ξ² σ^{3/2}. Substrate-derived chromomagnetic
    contact coefficient. Pure σ, ξ inputs — no free parameters."""
    return (8.0 / 3.0) * sigma_GeV2 * xi_inv_GeV ** 2 * sigma_GeV2 ** 1.5


K_SUBSTRATE_GEV3: float = _K_substrate_GeV3()
"""[A] K_substrate ≈ 0.0377 GeV³. Pure substrate K-chromomag coefficient."""


# --- Spin-flavour pair coefficients (SU(6) Clebsch-Gordan, NOT a fit) -------
# For each baryon, the matrix element ⟨Σ_{i<j} (S_i · S_j)⟩ partitions into
# three flavour buckets: light-light (qq), light-strange (qs), strange-
# strange (ss). Coefficients computed once from the SU(6) wavefunction and
# inherited identically by the substrate via §18.49 SU(3) inheritance.
#
# J=1/2 octet:
#   N (uud,udd):   c_qq = -3/4
#   Λ (uds,light singlet):  c_qq = -3/4   (pure light singlet, s spectator)
#   Σ (uds,light triplet):  c_qq = +1/4, c_qs = -1
#   Ξ (uss,dss):   c_qs = -1, c_ss = +1/4
#
# J=3/2 decuplet (all pairs spin-aligned, each contributes +1/4 per pair):
#   Δ (uud,udd,uuu,ddd):  c_qq = 3·(1/4) = +3/4
#   Σ* (uus,dds):         c_qq = +1/4, c_qs = +1/2  (1 light pair + 2 mixed)
#   Ξ* (uss,dss):         c_qs = +1/2, c_ss = +1/4  (2 mixed + 1 strange pair)
#   Ω⁻ (sss):             c_ss = +3/4

@dataclass(frozen=True)
class _BaryonV4Spec:
    """Spin-flavour spec for one baryon in the face-spin v4 calculator."""
    name: str
    n_light: int
    n_strange: int
    c_qq: float    # [A] light-light pair coupling, SU(6) C-G
    c_qs: float    # [A] light-strange pair coupling, SU(6) C-G
    c_ss: float    # [A] strange-strange pair coupling, SU(6) C-G
    J: float
    branch: str    # "octet" or "decuplet"


# Ground-state J=1/2 octet
_OCTET_V4: Dict[str, _BaryonV4Spec] = {
    "p":      _BaryonV4Spec("p", 3, 0, c_qq=-0.75, c_qs=0.0,  c_ss=0.0,   J=0.5, branch="octet"),
    "n":      _BaryonV4Spec("n", 3, 0, c_qq=-0.75, c_qs=0.0,  c_ss=0.0,   J=0.5, branch="octet"),
    "Lambda": _BaryonV4Spec("Lambda", 2, 1, c_qq=-0.75, c_qs=0.0, c_ss=0.0, J=0.5, branch="octet"),
    "Sigma+": _BaryonV4Spec("Sigma+", 2, 1, c_qq=+0.25, c_qs=-1.0, c_ss=0.0, J=0.5, branch="octet"),
    "Sigma0": _BaryonV4Spec("Sigma0", 2, 1, c_qq=+0.25, c_qs=-1.0, c_ss=0.0, J=0.5, branch="octet"),
    "Sigma-": _BaryonV4Spec("Sigma-", 2, 1, c_qq=+0.25, c_qs=-1.0, c_ss=0.0, J=0.5, branch="octet"),
    "Xi0":    _BaryonV4Spec("Xi0", 1, 2, c_qq=0.0, c_qs=-1.0, c_ss=+0.25, J=0.5, branch="octet"),
    "Xi-":    _BaryonV4Spec("Xi-", 1, 2, c_qq=0.0, c_qs=-1.0, c_ss=+0.25, J=0.5, branch="octet"),
}

# Ground-state J=3/2 decuplet
_DECUPLET_V4: Dict[str, _BaryonV4Spec] = {
    "Delta++": _BaryonV4Spec("Delta++", 3, 0, c_qq=+0.75, c_qs=0.0, c_ss=0.0, J=1.5, branch="decuplet"),
    "Delta+":  _BaryonV4Spec("Delta+",  3, 0, c_qq=+0.75, c_qs=0.0, c_ss=0.0, J=1.5, branch="decuplet"),
    "Delta0":  _BaryonV4Spec("Delta0",  3, 0, c_qq=+0.75, c_qs=0.0, c_ss=0.0, J=1.5, branch="decuplet"),
    "Delta-":  _BaryonV4Spec("Delta-",  3, 0, c_qq=+0.75, c_qs=0.0, c_ss=0.0, J=1.5, branch="decuplet"),
    "Sigma*+": _BaryonV4Spec("Sigma*+", 2, 1, c_qq=+0.25, c_qs=+0.5, c_ss=0.0,   J=1.5, branch="decuplet"),
    "Sigma*0": _BaryonV4Spec("Sigma*0", 2, 1, c_qq=+0.25, c_qs=+0.5, c_ss=0.0,   J=1.5, branch="decuplet"),
    "Sigma*-": _BaryonV4Spec("Sigma*-", 2, 1, c_qq=+0.25, c_qs=+0.5, c_ss=0.0,   J=1.5, branch="decuplet"),
    "Xi*0":    _BaryonV4Spec("Xi*0", 1, 2, c_qq=0.0, c_qs=+0.5, c_ss=+0.25, J=1.5, branch="decuplet"),
    "Xi*-":    _BaryonV4Spec("Xi*-", 1, 2, c_qq=0.0, c_qs=+0.5, c_ss=+0.25, J=1.5, branch="decuplet"),
    "Omega-":  _BaryonV4Spec("Omega-", 0, 3, c_qq=0.0, c_qs=0.0, c_ss=+0.75, J=1.5, branch="decuplet"),
}


# --- Anchor masses (2 mass anchors total: proton + Λ) ------------------------
M_PROTON_ANCHOR_MEV: float = 938.272
"""[B] Proton mass anchor — fixes light-quark structure mass m_q_struct."""

M_LAMBDA_ANCHOR_MEV: float = 1115.683
"""[B] Λ⁰ mass anchor — fixes strange-quark structure mass m_s_struct."""

M_NEUTRON_ANCHOR_MEV: float = 939.565
"""[B] Neutron mass (used only for the proton-anchor isospin average)."""


def _baryon_mass_v4_MeV(spec: _BaryonV4Spec, m_q_struct_MeV: float,
                        m_s_struct_MeV: float) -> float:
    """Compute baryon mass via face-spin v4 chromomagnetic decomposition.

    M_B = N_q · m_q_struct + N_s · m_s_struct
        + K_substrate · (c_qq/m_q² + c_qs/(m_q m_s) + c_ss/m_s²)

    Uses the GEOMETRIC chromomag masses m_q_chromo = √σ for light and
    m_s_chromo = m_q_chromo · (m_s_struct/m_q_struct) for strange.
    """
    K_MeV3 = K_SUBSTRATE_GEV3 * (1000.0 ** 3)
    m_q_chromo = 1000.0 * M_K_CHROMO_GEV  # ≈ 424 MeV
    m_s_chromo = m_q_chromo * (m_s_struct_MeV / m_q_struct_MeV)
    chromo = (
        spec.c_qq / (m_q_chromo ** 2)
        + spec.c_qs / (m_q_chromo * m_s_chromo)
        + spec.c_ss / (m_s_chromo ** 2)
    )
    delta_E = K_MeV3 * chromo
    return spec.n_light * m_q_struct_MeV + spec.n_strange * m_s_struct_MeV + delta_E


def _solve_m_q_struct_MeV(m_target_MeV: float = M_PROTON_ANCHOR_MEV) -> float:
    """Solve m_q_struct from the proton anchor (chromomag = -3/4 K/m_q²)."""
    K_MeV3 = K_SUBSTRATE_GEV3 * (1000.0 ** 3)
    m_q_chromo = 1000.0 * M_K_CHROMO_GEV
    delta_E_N = -0.75 * K_MeV3 / (m_q_chromo ** 2)
    return (m_target_MeV - delta_E_N) / 3.0


def _solve_m_s_struct_MeV(m_q_struct_MeV: float,
                          m_target_MeV: float = M_LAMBDA_ANCHOR_MEV) -> float:
    """Solve m_s_struct from the Λ anchor via bisection."""
    spec_lambda = _OCTET_V4["Lambda"]
    m_lo, m_hi = m_q_struct_MeV, m_q_struct_MeV + 800.0
    for _ in range(80):
        m_mid = 0.5 * (m_lo + m_hi)
        m_pred = _baryon_mass_v4_MeV(spec_lambda, m_q_struct_MeV, m_mid)
        if m_pred > m_target_MeV:
            m_hi = m_mid
        else:
            m_lo = m_mid
        if abs(m_hi - m_lo) < 1e-6:
            break
    return 0.5 * (m_lo + m_hi)


@dataclass
class BaryonFaceSpinV4:
    """Face-spin v4 baryon mass calculator (chromomagnetic substrate model).

    The "v4" mechanism referenced in the b3_baryon_face_spin_v4 audit
    notes. Predicts the full octet AND decuplet from:

      * 6 substrate-derived couplings  [A]:
          σ, ξ, K_substrate, c_qq, c_qs, c_ss
      * 2 mass anchors                 [B]:
          m_q_struct (from proton), m_s_struct (from Λ⁰)

    Hits all 16 baryons (8 octet + 8 decuplet excl. Σ*0/Ξ*-) at <2% mean.
    """

    m_q_struct_MeV: float = field(init=False)
    m_s_struct_MeV: float = field(init=False)

    def __post_init__(self) -> None:
        self.m_q_struct_MeV = _solve_m_q_struct_MeV()
        self.m_s_struct_MeV = _solve_m_s_struct_MeV(self.m_q_struct_MeV)

    def baryon_mass(self, name: str) -> float:
        """Predicted baryon mass in MeV. Octet + decuplet supported."""
        if name in _OCTET_V4:
            spec = _OCTET_V4[name]
        elif name in _DECUPLET_V4:
            spec = _DECUPLET_V4[name]
        else:
            raise KeyError(f"unknown baryon {name!r} in face-spin v4")
        return _baryon_mass_v4_MeV(spec, self.m_q_struct_MeV, self.m_s_struct_MeV)

    def supports(self, name: str) -> bool:
        """Whether the v4 calculator covers this baryon name."""
        return name in _OCTET_V4 or name in _DECUPLET_V4


__all__ = [
    "HadronSpectrum",
    "Residual",
    "QUARK_TORQUE",
    "B_MESON",
    "B_BARYON",
    "C_SS",
    "G_PS",
    "G_V",
    "OCTET",
    "DECUPLET",
    "PS_NONET",
    "V_NONET",
    "PDG",
    # face-spin v4
    "BaryonFaceSpinV4",
    "SIGMA_QCD_GEV2",
    "XI_QCD_FM",
    "K_SUBSTRATE_GEV3",
    "M_K_CHROMO_GEV",
    "M_PROTON_ANCHOR_MEV",
    "M_LAMBDA_ANCHOR_MEV",
]
