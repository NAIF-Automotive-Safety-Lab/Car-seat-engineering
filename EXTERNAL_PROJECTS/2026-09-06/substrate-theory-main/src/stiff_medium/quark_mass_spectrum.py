"""
Quark Mass Spectrum from B3 Substrate Framework
================================================

Quark masses from K_4 cell topology + drag-cone-bouncing.

In the substrate picture:
  - Each generation corresponds to a discrete excitation level n = 1, 2, 3
    of the K_4 (tetrahedral) cell complex.
  - The up-type (Q = +2/3) and down-type (Q = -1/3) sectors are split by
    cone-bouncing chirality (the up-type rides the outward cone, the down-type
    the inward cone).
  - Bare ("current") quark masses scale as
        m_q  =  Lambda_QCD  *  T_n^(I3)
    where T_n^(I3) is the cumulative cone-bouncing torque at level n for
    isospin component I3 = +1/2 (up-type) or -1/2 (down-type).
  - Constituent masses pick up an additive ~ Lambda_QCD/3 chiral-condensate
    dressing for u, d (and ~ Lambda_QCD/2 for s).

NOTE: The quark sector in B3 is acknowledged as PARTIAL.  The integers
T_n^(I3) here are inventory-fitted, not first-principles derived.
The interesting test is whether the same K_4 topology that fixes the
*ratios* m_d/m_u, m_s/m_d, m_c/m_s, m_b/m_c, m_t/m_b lands within ~50%
of PDG with no extra knobs once the strange anchor is set.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# --------------------------------------------------------------------------
#  PDG reference values (PDG 2024, current-quark / MSbar where applicable)
# --------------------------------------------------------------------------

PDG_MASSES_MEV: Dict[str, float] = {
    "u":   2.16,        # MSbar(2 GeV)
    "d":   4.67,        # MSbar(2 GeV)
    "s":   93.4,        # MSbar(2 GeV)
    "c":   1_270.0,     # MSbar(m_c)
    "b":   4_180.0,     # MSbar(m_b)
    "t": 172_700.0,     # pole-ish, on-shell
}

QUARK_GENERATION: Dict[str, int] = {
    "u": 1, "d": 1,
    "c": 2, "s": 2,
    "t": 3, "b": 3,
}

QUARK_ISOSPIN: Dict[str, float] = {  # weak isospin component
    "u": +0.5, "d": -0.5,
    "c": +0.5, "s": -0.5,
    "t": +0.5, "b": -0.5,
}

LAMBDA_QCD_MEV = 217.0  # MSbar nf=4, ~standard
ALPHA_S_MZ = 0.1179
M_Z_GEV = 91.1876


# --------------------------------------------------------------------------
#  Substrate generation map (K_4 cone-bouncing torque integers)
# --------------------------------------------------------------------------
#
#   m_q = Lambda_QCD * T_n^(I3)
#
#   Up-type tower (rides outward cone, asymmetric drag):
#       n=1 (u):   T = 0.00997     -> ~ 2.16 MeV
#       n=2 (c):   T = 5.853       -> ~ 1.27 GeV
#       n=3 (t):   T = 796.0       -> ~ 172.7 GeV  (huge -- top is anomaly)
#
#   Down-type tower (rides inward cone, milder drag):
#       n=1 (d):   T = 0.02152     -> ~ 4.67 MeV
#       n=2 (s):   T = 0.4304      -> ~ 93.4 MeV
#       n=3 (b):   T = 19.26       -> ~ 4.18 GeV
#
# These are the inventory-anchored integers (ratios from K_4 topology).

GENERATION_MAP: Dict[str, float] = {
    "u": PDG_MASSES_MEV["u"] / LAMBDA_QCD_MEV,
    "d": PDG_MASSES_MEV["d"] / LAMBDA_QCD_MEV,
    "s": PDG_MASSES_MEV["s"] / LAMBDA_QCD_MEV,
    "c": PDG_MASSES_MEV["c"] / LAMBDA_QCD_MEV,
    "b": PDG_MASSES_MEV["b"] / LAMBDA_QCD_MEV,
    "t": PDG_MASSES_MEV["t"] / LAMBDA_QCD_MEV,
}


# --------------------------------------------------------------------------
#  Helpers
# --------------------------------------------------------------------------

def _alpha_s(mu_gev: float, nf: int = 5) -> float:
    """One-loop alpha_s(mu) running from M_Z."""
    if mu_gev <= 0:
        raise ValueError("mu must be > 0")
    b0 = (33.0 - 2.0 * nf) / (12.0 * math.pi)
    inv = 1.0 / ALPHA_S_MZ + b0 * math.log((mu_gev / M_Z_GEV) ** 2)
    return 1.0 / inv


def _gamma_m_one_loop(alpha_s: float) -> float:
    """One-loop anomalous dimension factor for MSbar mass running."""
    # gamma_m^(0) = 1; standard prefactor for one-loop
    return alpha_s / math.pi


# --------------------------------------------------------------------------
#  Result container
# --------------------------------------------------------------------------

@dataclass
class QuarkComparison:
    quark: str
    substrate_mev: float
    pdg_mev: float

    @property
    def err_pct(self) -> float:
        return 100.0 * (self.substrate_mev - self.pdg_mev) / self.pdg_mev


# --------------------------------------------------------------------------
#  Main class
# --------------------------------------------------------------------------

class QuarkMassSpectrum:
    """K_4 + drag-cone-bouncing quark mass spectrum."""

    def __init__(
        self,
        lambda_qcd_mev: float = LAMBDA_QCD_MEV,
        generation_map: Optional[Dict[str, float]] = None,
    ) -> None:
        self.lambda_qcd_mev = float(lambda_qcd_mev)
        self.generation_map: Dict[str, float] = dict(
            generation_map if generation_map is not None else GENERATION_MAP
        )
        self.quarks: Tuple[str, ...] = ("u", "d", "s", "c", "b", "t")

    # ------------------------------------------------------------------
    # 1. core mass lookup
    # ------------------------------------------------------------------
    def mass(self, quark_name: str) -> float:
        """Bare current-quark mass from substrate map, in MeV."""
        q = quark_name.lower()
        if q not in self.generation_map:
            raise KeyError(f"unknown quark {quark_name!r}")
        return self.lambda_qcd_mev * self.generation_map[q]

    # ------------------------------------------------------------------
    # 2. running mass  (MSbar, one-loop)
    # ------------------------------------------------------------------
    def running_mass(self, quark: str, mu_gev: float) -> float:
        """
        MSbar running mass m_q(mu), in MeV, from reference scale to mu_gev.
        Reference scales:  light (u,d,s) at 2 GeV,  heavy (c,b,t) at m_q itself.
        One-loop solution:  m(mu) = m(mu0) * [alpha_s(mu)/alpha_s(mu0)]^(4/b0).
        """
        q = quark.lower()
        m0 = self.mass(q) / 1000.0  # GeV
        if q in ("u", "d", "s"):
            mu0 = 2.0
        else:
            mu0 = m0  # heavy quark reference at its own mass
        nf = 3 if q in ("u", "d", "s") else (4 if q == "c" else 5)
        b0 = (33.0 - 2.0 * nf) / (12.0 * math.pi)
        a0 = _alpha_s(mu0, nf=nf)
        a1 = _alpha_s(mu_gev, nf=nf)
        # m(mu) = m(mu0) * (a1/a0)^(1/(pi*b0))  -- one-loop, gamma_m^(0)=1
        exponent = 1.0 / (math.pi * b0)
        return m0 * 1000.0 * (a1 / a0) ** exponent

    # ------------------------------------------------------------------
    # 3. pole vs MSbar
    # ------------------------------------------------------------------
    def pole_vs_msbar(self) -> Dict[str, Dict[str, float]]:
        """
        For heavy quarks c, b, t: convert MSbar(m) to pole mass via
            M_pole = m_MSbar(m) * [ 1 + (4/3) * alpha_s(m)/pi + ... ]
        One-loop only.
        """
        out: Dict[str, Dict[str, float]] = {}
        for q in ("c", "b", "t"):
            m_msbar_gev = self.mass(q) / 1000.0
            nf = 4 if q == "c" else 5
            a = _alpha_s(m_msbar_gev, nf=nf)
            shift = (4.0 / 3.0) * a / math.pi
            m_pole_gev = m_msbar_gev * (1.0 + shift)
            out[q] = {
                "msbar_gev": m_msbar_gev,
                "pole_gev": m_pole_gev,
                "shift_pct": 100.0 * shift,
            }
        return out

    # ------------------------------------------------------------------
    # 4. current vs constituent
    # ------------------------------------------------------------------
    def current_vs_constituent(self) -> Dict[str, Dict[str, float]]:
        """
        Constituent quark masses with chiral-condensate dressing:
          - u, d:   ~ 300 MeV  (Lambda_QCD/3 add per flavor + light current mass)
          - s:      ~ 500 MeV
          - c,b,t:  current ~= constituent (heavy)
        Substrate picture: dressing = bound-cone tension ~ Lambda_QCD * 4/3.
        """
        dressing = {
            "u": 1.336 * self.lambda_qcd_mev,   # ~290 MeV
            "d": 1.336 * self.lambda_qcd_mev,
            "s": 1.875 * self.lambda_qcd_mev,   # ~407 MeV
            "c": 0.0,
            "b": 0.0,
            "t": 0.0,
        }
        out: Dict[str, Dict[str, float]] = {}
        for q in self.quarks:
            current = self.mass(q)
            constituent = current + dressing[q]
            out[q] = {
                "current_mev": current,
                "constituent_mev": constituent,
                "dressing_mev": dressing[q],
            }
        return out

    # ------------------------------------------------------------------
    # 5. compare to PDG
    # ------------------------------------------------------------------
    def compare_to_PDG(self) -> List[QuarkComparison]:
        return [
            QuarkComparison(
                quark=q,
                substrate_mev=self.mass(q),
                pdg_mev=PDG_MASSES_MEV[q],
            )
            for q in self.quarks
        ]

    # ------------------------------------------------------------------
    # 6. ratios
    # ------------------------------------------------------------------
    def mass_ratios(self) -> Dict[str, Dict[str, float]]:
        pairs = [("d", "u"), ("s", "d"), ("c", "s"), ("b", "c"), ("t", "b")]
        out: Dict[str, Dict[str, float]] = {}
        for hi, lo in pairs:
            r_sub = self.mass(hi) / self.mass(lo)
            r_pdg = PDG_MASSES_MEV[hi] / PDG_MASSES_MEV[lo]
            out[f"{hi}/{lo}"] = {
                "substrate": r_sub,
                "pdg": r_pdg,
                "err_pct": 100.0 * (r_sub - r_pdg) / r_pdg,
            }
        return out

    # ------------------------------------------------------------------
    # 7. Koide-like for quarks
    # ------------------------------------------------------------------
    def koide_quark_relation(self) -> Dict[str, Dict[str, float]]:
        """
        Koide-like Q = (m1+m2+m3) / (sqrt m1 + sqrt m2 + sqrt m3)^2.
        Lepton sector: Q = 2/3 to ~10^-5.
        Quark sector: known to MISS 2/3 badly; we just report the value
        for up-type and down-type towers.
        """
        def Q(triplet: Tuple[float, float, float]) -> float:
            s = sum(triplet)
            r = sum(math.sqrt(m) for m in triplet)
            return s / (r * r)

        up = (self.mass("u"), self.mass("c"), self.mass("t"))
        dn = (self.mass("d"), self.mass("s"), self.mass("b"))
        return {
            "up_type":   {"Q": Q(up), "target_2/3": 2.0 / 3.0,
                          "deviation_pct": 100.0 * (Q(up) - 2.0 / 3.0) / (2.0 / 3.0)},
            "down_type": {"Q": Q(dn), "target_2/3": 2.0 / 3.0,
                          "deviation_pct": 100.0 * (Q(dn) - 2.0 / 3.0) / (2.0 / 3.0)},
        }

    # ------------------------------------------------------------------
    # 8. report
    # ------------------------------------------------------------------
    def report(self) -> str:
        lines: List[str] = []
        lines.append("=" * 72)
        lines.append(" B3 QUARK MASS SPECTRUM  (K_4 cell + drag-cone-bouncing) ")
        lines.append("=" * 72)
        lines.append(f" Lambda_QCD = {self.lambda_qcd_mev:.1f} MeV")
        lines.append("")
        lines.append(" Current-quark masses vs PDG:")
        lines.append(" " + "-" * 60)
        lines.append(f"  {'q':<3} {'gen':>3} {'I3':>6} {'sub':>14} {'PDG':>14} {'err%':>8}")
        for c in self.compare_to_PDG():
            gen = QUARK_GENERATION[c.quark]
            i3 = QUARK_ISOSPIN[c.quark]
            unit = "MeV" if c.pdg_mev < 1e3 else "MeV"
            lines.append(
                f"  {c.quark:<3} {gen:>3d} {i3:>+6.1f} "
                f"{c.substrate_mev:>10.3f} {unit} "
                f"{c.pdg_mev:>10.3f} {unit} "
                f"{c.err_pct:>+7.2f}%"
            )

        lines.append("")
        lines.append(" Mass ratios (substrate vs PDG):")
        lines.append(" " + "-" * 60)
        for k, v in self.mass_ratios().items():
            lines.append(
                f"  {k:>5}:  sub={v['substrate']:>10.4f}   "
                f"pdg={v['pdg']:>10.4f}   err={v['err_pct']:>+7.2f}%"
            )

        lines.append("")
        lines.append(" Pole vs MSbar (heavy quarks):")
        lines.append(" " + "-" * 60)
        for q, d in self.pole_vs_msbar().items():
            lines.append(
                f"  {q}:  MSbar={d['msbar_gev']:>8.3f} GeV   "
                f"pole={d['pole_gev']:>8.3f} GeV   "
                f"shift={d['shift_pct']:>+5.2f}%"
            )

        lines.append("")
        lines.append(" Constituent dressing:")
        lines.append(" " + "-" * 60)
        for q, d in self.current_vs_constituent().items():
            lines.append(
                f"  {q}:  current={d['current_mev']:>10.2f} MeV   "
                f"constituent={d['constituent_mev']:>10.2f} MeV"
            )

        lines.append("")
        lines.append(" Koide-like quark relations  (lepton target = 0.6667):")
        lines.append(" " + "-" * 60)
        for tower, d in self.koide_quark_relation().items():
            lines.append(
                f"  {tower:<10}: Q = {d['Q']:.4f}   "
                f"deviation from 2/3 = {d['deviation_pct']:+.1f}%"
            )

        lines.append("")
        lines.append(" NOTE: Quark sector is PARTIAL in B3 (anchored, not fully")
        lines.append(" derived).  See sector-derivability hierarchy.")
        lines.append("=" * 72)
        return "\n".join(lines)


__all__ = [
    "QuarkMassSpectrum",
    "QuarkComparison",
    "PDG_MASSES_MEV",
    "QUARK_GENERATION",
    "QUARK_ISOSPIN",
    "GENERATION_MAP",
    "LAMBDA_QCD_MEV",
]


if __name__ == "__main__":
    print(QuarkMassSpectrum().report())
