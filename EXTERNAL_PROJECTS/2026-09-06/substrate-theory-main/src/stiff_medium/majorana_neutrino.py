"""Majorana neutrino simulator and 0nubb predictor (B3 substrate framework).

The B3 substrate FORCES Majorana neutrinos: neutrality plus the Mobius Z/2
sheet-swap means the two sheets are gauge-identified, so a neutrino is its
own antiparticle.  This module turns substrate-derived inputs (Sigma m_nu,
m_1, PMNS angles from the alpha-anchor) into a falsifiable prediction for
the effective Majorana mass m_bb probed by neutrinoless double beta decay.

Substrate inputs (NH baseline, all from B3 derivations):
    Sigma m_nu = 60.5 meV       (B3 Hubble derivation v2, NH)
    m_1        = 2.26 meV       (lightest neutrino, NH)
    sin^2 t12  = 42 alpha       ~ 0.306
    sin^2 t13  = 3  alpha       ~ 0.0219
    sin^2 t23  = 1/2 + 2 pi alpha ~ 0.546

Mass splittings used to deduce m_2, m_3:
    Dm21^2 = 7.5e-5 eV^2
    Dm31^2 = 2.5e-3 eV^2

The Majorana phases (alpha_1, alpha_2) are physically unconstrained, so we
scan them to bracket m_bb = | sum_i U_ei^2 m_i exp(i alpha_i) |.

Experimental targets:
    KamLAND-Zen 2024     m_bb < 36 meV     (current best)
    LEGEND-200           m_bb < ~50 meV
    LEGEND-1000          m_bb sensitivity ~9-21 meV (~2030)
    nEXO                 m_bb sensitivity ~5-12 meV (~2032)
    CUPID                m_bb sensitivity ~10-20 meV (~2031)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np
from scipy import constants


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ALPHA: float = constants.alpha  # fine-structure constant ~ 1/137.036
EV: float = constants.eV         # 1 eV in Joules (kept for unit interop)

# Default substrate-derived B3 inputs (NH).
DEFAULT_SIGMA_MNU_EV: float = 60.5e-3       # 60.5 meV
DEFAULT_M1_EV: float = 2.26e-3              # 2.26 meV (NH lightest)

# PDG-style oscillation splittings (eV^2).
DM21_SQ_EV2: float = 7.5e-5
DM31_SQ_EV2: float = 2.5e-3

# Substrate PMNS angles (sin^2 of mixing angles).
SIN2_T12_SUBSTRATE: float = 42.0 * ALPHA
SIN2_T13_SUBSTRATE: float = 3.0 * ALPHA
SIN2_T23_SUBSTRATE: float = 0.5 + 2.0 * np.pi * ALPHA


# ---------------------------------------------------------------------------
# Experimental sensitivity table
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Experiment:
    name: str
    m_bb_low_meV: float          # low end of m_bb sensitivity (or upper limit)
    m_bb_high_meV: float
    year: int
    is_limit: bool = False       # True for current upper limits, False for projected sensitivity


EXPERIMENTS: Tuple[Experiment, ...] = (
    Experiment("KamLAND-Zen 2024", 36.0, 36.0, 2024, is_limit=True),
    Experiment("LEGEND-200",        50.0, 50.0, 2026, is_limit=True),
    Experiment("LEGEND-1000",        9.0, 21.0, 2030),
    Experiment("nEXO",               5.0, 12.0, 2032),
    Experiment("CUPID",             10.0, 20.0, 2031),
)


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------

@dataclass
class MajoranaNeutrino:
    """B3 substrate prediction engine for 0nubb effective mass."""

    sigma_mnu_eV: float = DEFAULT_SIGMA_MNU_EV
    m1_eV: float = DEFAULT_M1_EV
    sin2_t12: float = SIN2_T12_SUBSTRATE
    sin2_t13: float = SIN2_T13_SUBSTRATE
    sin2_t23: float = SIN2_T23_SUBSTRATE
    dm21_sq_eV2: float = DM21_SQ_EV2
    dm31_sq_eV2: float = DM31_SQ_EV2
    hierarchy: str = "NH"        # "NH" or "IH"

    # Filled in __post_init__
    masses_eV: np.ndarray = field(init=False)
    U_e: np.ndarray = field(init=False)   # |U_ei|^2 for i=1,2,3

    def __post_init__(self) -> None:
        if self.hierarchy not in ("NH", "IH"):
            raise ValueError(f"hierarchy must be 'NH' or 'IH', got {self.hierarchy!r}")
        self.masses_eV = self._compute_masses()
        self.U_e = self._compute_Ue_squared()

    # ------------------------------------------------------------------
    # Mass spectrum from m_1 (or m_3 for IH) plus oscillation splittings
    # ------------------------------------------------------------------
    def _compute_masses(self) -> np.ndarray:
        if self.hierarchy == "NH":
            m1 = self.m1_eV
            m2 = float(np.sqrt(m1 * m1 + self.dm21_sq_eV2))
            m3 = float(np.sqrt(m1 * m1 + self.dm31_sq_eV2))
        else:  # IH: m_3 is lightest; we reuse m1_eV slot as m_3 input.
            m3 = self.m1_eV
            # IH convention: Dm31^2 < 0; use |Dm31^2| for spacing m_1 from m_3.
            m1 = float(np.sqrt(m3 * m3 + self.dm31_sq_eV2))
            m2 = float(np.sqrt(m1 * m1 + self.dm21_sq_eV2))
        return np.array([m1, m2, m3])

    def _compute_Ue_squared(self) -> np.ndarray:
        """|U_e1|^2, |U_e2|^2, |U_e3|^2 from the standard PMNS parametrisation."""
        c13_sq = 1.0 - self.sin2_t13
        Ue1_sq = c13_sq * (1.0 - self.sin2_t12)
        Ue2_sq = c13_sq * self.sin2_t12
        Ue3_sq = self.sin2_t13
        return np.array([Ue1_sq, Ue2_sq, Ue3_sq])

    # ------------------------------------------------------------------
    # Diagnostic getters
    # ------------------------------------------------------------------
    def sum_masses_eV(self) -> float:
        return float(np.sum(self.masses_eV))

    def sum_masses_meV(self) -> float:
        return self.sum_masses_eV() * 1.0e3

    # ------------------------------------------------------------------
    # m_bb evaluation
    # ------------------------------------------------------------------
    def m_bb(self, majorana_phases: Tuple[float, float] = (0.0, 0.0)) -> float:
        """Effective Majorana mass in eV.

        m_bb = | U_e1^2 m_1 + U_e2^2 m_2 e^{i alpha_1}
                + U_e3^2 m_3 e^{i alpha_2} |
        with the |U_ei|^2 absorbed (CP-violating Dirac phase ignored, since
        only the |U_ei|^2 enter for m_bb).
        """
        a1, a2 = majorana_phases
        amp = (
            self.U_e[0] * self.masses_eV[0]
            + self.U_e[1] * self.masses_eV[1] * np.exp(1j * a1)
            + self.U_e[2] * self.masses_eV[2] * np.exp(1j * a2)
        )
        return float(np.abs(amp))

    def m_bb_meV(self, majorana_phases: Tuple[float, float] = (0.0, 0.0)) -> float:
        return self.m_bb(majorana_phases) * 1.0e3

    # ------------------------------------------------------------------
    # Phase scan
    # ------------------------------------------------------------------
    def m_bb_range(self, n_grid: int = 181) -> Dict[str, float]:
        """Scan Majorana phases; return min, max, and the phases that hit them."""
        phases = np.linspace(0.0, 2.0 * np.pi, n_grid)
        a1g, a2g = np.meshgrid(phases, phases, indexing="ij")
        amp = (
            self.U_e[0] * self.masses_eV[0]
            + self.U_e[1] * self.masses_eV[1] * np.exp(1j * a1g)
            + self.U_e[2] * self.masses_eV[2] * np.exp(1j * a2g)
        )
        m = np.abs(amp)
        i_min = np.unravel_index(np.argmin(m), m.shape)
        i_max = np.unravel_index(np.argmax(m), m.shape)
        return {
            "min_eV": float(m.min()),
            "max_eV": float(m.max()),
            "min_meV": float(m.min()) * 1.0e3,
            "max_meV": float(m.max()) * 1.0e3,
            "alpha1_at_min": float(phases[i_min[0]]),
            "alpha2_at_min": float(phases[i_min[1]]),
            "alpha1_at_max": float(phases[i_max[0]]),
            "alpha2_at_max": float(phases[i_max[1]]),
        }

    # ------------------------------------------------------------------
    # Experiment comparison
    # ------------------------------------------------------------------
    def compare_to_experiments(self) -> Dict[str, Dict[str, object]]:
        """For each experiment, report whether B3 m_bb range is detectable / excluded."""
        rng = self.m_bb_range()
        b3_lo, b3_hi = rng["min_meV"], rng["max_meV"]
        out: Dict[str, Dict[str, object]] = {}
        for exp in EXPERIMENTS:
            exp_lo, exp_hi = exp.m_bb_low_meV, exp.m_bb_high_meV
            if exp.is_limit:
                # Current upper limit: B3 excluded iff its MIN exceeds the limit.
                excluded = b3_lo > exp_hi
                detectable = b3_hi >= exp_hi  # would graze the limit if alpha-phase lucky
                status = "EXCLUDED" if excluded else (
                    "TOUCHES LIMIT" if detectable else "PASSES (allowed)"
                )
            else:
                # Future sensitivity band: discriminating iff B3 range overlaps it
                # OR sits above its low end (so any signal would land in reach).
                covers = (b3_hi >= exp_lo)
                fully_inside = (b3_lo >= exp_lo) and (b3_hi <= exp_hi)
                status = "DECISIVE TEST" if covers else "BELOW REACH"
                if fully_inside:
                    status = "FULLY COVERED (decisive)"
                excluded = False
                detectable = covers
            out[exp.name] = {
                "year": exp.year,
                "exp_low_meV": exp_lo,
                "exp_high_meV": exp_hi,
                "B3_low_meV": b3_lo,
                "B3_high_meV": b3_hi,
                "status": status,
                "is_limit": exp.is_limit,
                "excluded": excluded,
                "detectable": detectable,
            }
        return out

    # ------------------------------------------------------------------
    # Falsification timeline
    # ------------------------------------------------------------------
    def predict_falsification_year(self) -> Dict[str, object]:
        """First future experiment whose sensitivity reaches the B3 m_bb upper end."""
        rng = self.m_bb_range()
        b3_hi = rng["max_meV"]
        candidates = []
        for exp in EXPERIMENTS:
            if exp.is_limit:
                continue
            if b3_hi >= exp.m_bb_low_meV:
                candidates.append(exp)
        if not candidates:
            return {
                "decisive_year": None,
                "experiment": None,
                "B3_max_meV": b3_hi,
                "note": "No planned experiment reaches B3 m_bb upper end.",
            }
        first = min(candidates, key=lambda e: e.year)
        return {
            "decisive_year": first.year,
            "experiment": first.name,
            "B3_max_meV": b3_hi,
            "exp_sensitivity_meV": (first.m_bb_low_meV, first.m_bb_high_meV),
            "note": f"{first.name} (~{first.year}) reaches the B3 upper m_bb.",
        }

    # ------------------------------------------------------------------
    # Inverted-hierarchy comparison
    # ------------------------------------------------------------------
    def inverted_hierarchy_check(self, m3_lightest_eV: float = 1.0e-3) -> Dict[str, float]:
        """Build an IH twin (same PMNS angles, splittings) and report its m_bb range."""
        ih = MajoranaNeutrino(
            sigma_mnu_eV=self.sigma_mnu_eV,
            m1_eV=m3_lightest_eV,        # interpreted as m_3 in IH branch
            sin2_t12=self.sin2_t12,
            sin2_t13=self.sin2_t13,
            sin2_t23=self.sin2_t23,
            dm21_sq_eV2=self.dm21_sq_eV2,
            dm31_sq_eV2=self.dm31_sq_eV2,
            hierarchy="IH",
        )
        rng = ih.m_bb_range()
        return {
            "IH_sum_mnu_meV": ih.sum_masses_meV(),
            "IH_m_bb_min_meV": rng["min_meV"],
            "IH_m_bb_max_meV": rng["max_meV"],
            "NH_m_bb_min_meV": self.m_bb_range()["min_meV"],
            "NH_m_bb_max_meV": self.m_bb_range()["max_meV"],
        }

    # ------------------------------------------------------------------
    # Pretty report
    # ------------------------------------------------------------------
    def report(self) -> str:
        lines = []
        lines.append("=" * 72)
        lines.append("B3 SUBSTRATE MAJORANA-NEUTRINO PREDICTION  (0nubb)")
        lines.append("=" * 72)
        lines.append(f"hierarchy           : {self.hierarchy}")
        lines.append(f"m_1  (lightest)     : {self.masses_eV[0]*1e3:.4f} meV")
        lines.append(f"m_2                 : {self.masses_eV[1]*1e3:.4f} meV")
        lines.append(f"m_3                 : {self.masses_eV[2]*1e3:.4f} meV")
        lines.append(f"Sigma m_nu          : {self.sum_masses_meV():.3f} meV"
                     f"   (substrate target {self.sigma_mnu_eV*1e3:.2f} meV)")
        lines.append("")
        lines.append("PMNS (from substrate alpha-anchor):")
        lines.append(f"  sin^2 t12 = {self.sin2_t12:.5f}    |U_e1|^2 = {self.U_e[0]:.5f}")
        lines.append(f"  sin^2 t13 = {self.sin2_t13:.5f}    |U_e2|^2 = {self.U_e[1]:.5f}")
        lines.append(f"  sin^2 t23 = {self.sin2_t23:.5f}    |U_e3|^2 = {self.U_e[2]:.5f}")
        rng = self.m_bb_range()
        lines.append("")
        lines.append("Effective Majorana mass m_bb (Majorana-phase scan):")
        lines.append(f"  min = {rng['min_meV']:.4f} meV   (alpha1={rng['alpha1_at_min']:.2f},"
                     f" alpha2={rng['alpha2_at_min']:.2f})")
        lines.append(f"  max = {rng['max_meV']:.4f} meV   (alpha1={rng['alpha1_at_max']:.2f},"
                     f" alpha2={rng['alpha2_at_max']:.2f})")
        lines.append("")
        lines.append("Experimental confrontation:")
        lines.append(f"  {'experiment':<22}{'year':>6}{'sens (meV)':>15}{'status':>30}")
        for name, row in self.compare_to_experiments().items():
            sens = (f"{row['exp_low_meV']:.0f}"
                    if row['is_limit']
                    else f"{row['exp_low_meV']:.0f}-{row['exp_high_meV']:.0f}")
            lines.append(f"  {name:<22}{row['year']:>6}{sens:>15}{row['status']:>30}")
        falsif = self.predict_falsification_year()
        lines.append("")
        if falsif["decisive_year"]:
            lines.append(f"Decisive test     : {falsif['experiment']}  (~{falsif['decisive_year']})")
        else:
            lines.append("Decisive test     : none currently planned reaches B3 max")
        ih = self.inverted_hierarchy_check()
        lines.append("")
        lines.append("Hierarchy comparison (B3 always says NH):")
        lines.append(f"  NH m_bb : {ih['NH_m_bb_min_meV']:.3f} - {ih['NH_m_bb_max_meV']:.3f} meV")
        lines.append(f"  IH m_bb : {ih['IH_m_bb_min_meV']:.3f} - {ih['IH_m_bb_max_meV']:.3f} meV"
                     f"   (IH Sigma m_nu = {ih['IH_sum_mnu_meV']:.2f} meV)")
        lines.append("=" * 72)
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Module-level convenience
# ---------------------------------------------------------------------------

def default_substrate_majorana() -> MajoranaNeutrino:
    """Construct the canonical B3 NH substrate Majorana neutrino prediction."""
    return MajoranaNeutrino()


if __name__ == "__main__":  # pragma: no cover
    print(default_substrate_majorana().report())
