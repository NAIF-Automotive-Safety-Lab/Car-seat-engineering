"""
generation_map.py
=================

Substrate generation map for Standard Model fermion mass ratios.

Hypothesis (B3 framework):
    The generation-to-generation mass ratio of charged leptons is
    governed by an exponential of the substrate "mode count" n_M
    divided by a topological denominator built from the pair count
    K_pair and the rank constant K_rank.

    Generation 1 -> 2 (e -> mu):
        ratio_12 = exp(n_M / (K_pair**4 * pi))
                 = exp(268 / (16 * pi))
        -> matches observed m_mu/m_e = 206.768 at ~0.0017 %.

    Generation 2 -> 3 (mu -> tau):
        ratio_23 = exp( n_M / (K_pair**4 * pi)
                        - (n_R - R) / (K_pair * (K_pair + 1)) )
                 = exp(268/(16*pi) - 15/6)
        -> matches observed m_tau/m_mu = 16.817 at ~0.94 %.
        (Canonical n_R-active form, shared with mass_torque_engine.
         Legacy form exp(n_M/(K_rank*(K_rank+1)*pi)) was 2.14 % off.)

    Generation 1 -> 3 (e -> tau):
        ratio_13 = ratio_12 * ratio_23.

This module exposes a small `GenerationMap` class that:
    * holds the substrate constants (n_M, K_pair, K_rank),
    * computes lepton ratios with the formula above,
    * cross-checks observed lepton masses,
    * applies the same form to neutrinos (oscillation data, NH/IH),
    * applies the same form to up-type and down-type quarks,
    * prints an honest scoreboard.

The lepton sector should pass at high precision (sub-percent on 1->2,
~few-percent on 2->3). The quark sector is expected to be 5-15 %
off; it is reported transparently as a partial fit. The neutrino
sector depends strongly on the hierarchy assumption.

References (PDG 2024 central values):
    m_e   = 0.51099895000   MeV
    m_mu  = 105.6583755     MeV
    m_tau = 1776.86         MeV
    m_u   = 2.16            MeV  (MS-bar 2 GeV)
    m_c   = 1273.0          MeV
    m_t   = 172_570.0       MeV  (pole)
    m_d   = 4.70            MeV
    m_s   = 93.5            MeV
    m_b   = 4183.0          MeV

Neutrino oscillation (NuFIT 5.2, NH):
    dm21^2 = 7.41e-5  eV^2
    dm31^2 = 2.511e-3 eV^2
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np

from .b3_constants import (
    n_M as _CANONICAL_N_M,
    K_pair as _CANONICAL_K_PAIR,
    K_rank as _CANONICAL_K_RANK,
    n_R as _CANONICAL_N_R,
    R as _CANONICAL_R,
)


# ---------------------------------------------------------------------------
# Observed fermion masses (MeV unless noted). PDG 2024 central values.
# ---------------------------------------------------------------------------

LEPTON_MASSES_MEV: Dict[str, float] = {
    "e":   0.51099895000,
    "mu":  105.6583755,
    "tau": 1776.86,
}

UP_QUARK_MASSES_MEV: Dict[str, float] = {
    "u": 2.16,
    "c": 1273.0,
    "t": 172_570.0,
}

DOWN_QUARK_MASSES_MEV: Dict[str, float] = {
    "d": 4.70,
    "s": 93.5,
    "b": 4183.0,
}

# Neutrino mass-squared splittings (eV^2), NuFIT 5.2 normal hierarchy.
DM21_SQ_EV2: float = 7.41e-5
DM31_SQ_EV2: float = 2.511e-3


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GenerationMap:
    """
    Substrate generation map.

    Parameters
    ----------
    n_M : int
        Substrate mode count (B3: 268).
    K_pair : int
        Pair count (B3: 2).
    K_rank : int
        Rank constant (B3: 5).
    """

    n_M: int = _CANONICAL_N_M     # 268
    K_pair: int = _CANONICAL_K_PAIR  # 2
    K_rank: int = _CANONICAL_K_RANK  # 5
    n_R: int = _CANONICAL_N_R     # 18
    R: int = _CANONICAL_R         # 3

    # ------------------------------------------------------------------ core

    def _denominator_12(self) -> float:
        """Denominator in the exponent for the 1 -> 2 step."""
        return (self.K_pair ** 4) * np.pi

    def _denominator_23(self) -> float:
        """Denominator in the exponent for the 2 -> 3 step.

        DEPRECATED: legacy form exp(n_M / (K_rank*(K_rank+1)*pi)) gave
        ratio 17.18 (2.14% off PDG). Retained only for back-compat;
        the canonical lepton tower exponent now lives in
        :meth:`_log_ratio_23` and matches mass_torque_engine.
        """
        return self.K_rank * (self.K_rank + 1) * np.pi

    def _log_ratio_23(self) -> float:
        """Log of m_tau / m_mu (canonical, n_R-active form).

        Delegates to :func:`tau_mass_unified.log_m_tau_over_m_mu`, the
        single source of truth shared with ``mass_torque_engine._t_tau``
        and ``integer_rigidity.predict_m_tau_over_m_e``:

            ln(m_tau/m_mu) = n_M / (K_pair^4 * pi)
                             - (n_R - R) / (K_pair*(K_pair+1))

        With (n_M, K_pair, K_rank, n_R, R) = (268, 2, 5, 18, 3) this
        gives 16.97 (0.93% off PDG 16.817), and makes n_R an active
        observable input rather than a derived identity.
        """
        from .tau_mass_unified import log_m_tau_over_m_mu
        return log_m_tau_over_m_mu(
            n_M=self.n_M, K_pair=self.K_pair, n_R=self.n_R, R=self.R,
        )

    def lepton_ratio(self, generation_high: int, generation_low: int) -> float:
        """
        Predict the mass ratio m(generation_high) / m(generation_low).

        Valid pairs: (2,1), (3,2), (3,1) and their inverses.
        """
        if generation_high == generation_low:
            return 1.0
        if generation_high < generation_low:
            return 1.0 / self.lepton_ratio(generation_low, generation_high)

        if (generation_high, generation_low) == (2, 1):
            return float(np.exp(self.n_M / self._denominator_12()))
        if (generation_high, generation_low) == (3, 2):
            return float(np.exp(self._log_ratio_23()))
        if (generation_high, generation_low) == (3, 1):
            return self.lepton_ratio(2, 1) * self.lepton_ratio(3, 2)

        raise ValueError(
            f"Unsupported generation pair: ({generation_high}, {generation_low})"
        )

    # --------------------------------------------------------------- leptons

    def verify_lepton_masses(self) -> Dict[str, Dict[str, float]]:
        """
        Compare predicted lepton mass ratios against observed values.

        Returns a dict keyed by ratio name, each entry containing
        predicted, observed, and percent-error.
        """
        observed_mu_e = LEPTON_MASSES_MEV["mu"] / LEPTON_MASSES_MEV["e"]
        observed_tau_mu = LEPTON_MASSES_MEV["tau"] / LEPTON_MASSES_MEV["mu"]
        observed_tau_e = LEPTON_MASSES_MEV["tau"] / LEPTON_MASSES_MEV["e"]

        pred_mu_e = self.lepton_ratio(2, 1)
        pred_tau_mu = self.lepton_ratio(3, 2)
        pred_tau_e = self.lepton_ratio(3, 1)

        return {
            "mu/e": _row(pred_mu_e, observed_mu_e),
            "tau/mu": _row(pred_tau_mu, observed_tau_mu),
            "tau/e": _row(pred_tau_e, observed_tau_e),
        }

    # ------------------------------------------------------------- neutrinos

    def predict_neutrino_ratios(
        self, hierarchy: str = "NH"
    ) -> Dict[str, Dict[str, float]]:
        """
        Apply the lepton generation map to neutrinos and compare against
        observed mass ratios derived from oscillation splittings.

        For NH with lightest mass m1 small (taken from B3 prediction
        m1 ~ 2.26 meV, but we evaluate on splittings only):
            m2 = sqrt(m1^2 + dm21^2)
            m3 = sqrt(m1^2 + dm31^2)

        Returns dict with predicted vs observed nu2/nu1, nu3/nu1.
        """
        hierarchy = hierarchy.upper()
        if hierarchy not in ("NH", "IH"):
            raise ValueError("hierarchy must be 'NH' or 'IH'")

        # Use B3 lightest-neutrino prediction (eV).
        m1_eV = 2.26e-3 if hierarchy == "NH" else 5.0e-2

        if hierarchy == "NH":
            m2 = np.sqrt(m1_eV ** 2 + DM21_SQ_EV2)
            m3 = np.sqrt(m1_eV ** 2 + DM31_SQ_EV2)
        else:
            # IH: m3 is lightest
            m3 = m1_eV
            m1 = np.sqrt(m3 ** 2 + DM31_SQ_EV2)
            m2 = np.sqrt(m1 ** 2 + DM21_SQ_EV2)
            m1_eV = m1

        observed_2_1 = m2 / m1_eV
        observed_3_1 = m3 / m1_eV
        observed_3_2 = m3 / m2

        pred_2_1 = self.lepton_ratio(2, 1)
        pred_3_2 = self.lepton_ratio(3, 2)
        pred_3_1 = self.lepton_ratio(3, 1)

        return {
            f"nu2/nu1 ({hierarchy})": _row(pred_2_1, observed_2_1),
            f"nu3/nu2 ({hierarchy})": _row(pred_3_2, observed_3_2),
            f"nu3/nu1 ({hierarchy})": _row(pred_3_1, observed_3_1),
        }

    # ---------------------------------------------------------------- quarks

    def predict_quark_ratios(self) -> Dict[str, Dict[str, float]]:
        """
        Apply the lepton generation map verbatim to up-type and down-type
        quarks. Honest expectation: 5-15 % off, sector-dependent.
        """
        out: Dict[str, Dict[str, float]] = {}

        # Up-type
        obs_c_u = UP_QUARK_MASSES_MEV["c"] / UP_QUARK_MASSES_MEV["u"]
        obs_t_c = UP_QUARK_MASSES_MEV["t"] / UP_QUARK_MASSES_MEV["c"]
        obs_t_u = UP_QUARK_MASSES_MEV["t"] / UP_QUARK_MASSES_MEV["u"]

        out["c/u"] = _row(self.lepton_ratio(2, 1), obs_c_u)
        out["t/c"] = _row(self.lepton_ratio(3, 2), obs_t_c)
        out["t/u"] = _row(self.lepton_ratio(3, 1), obs_t_u)

        # Down-type
        obs_s_d = DOWN_QUARK_MASSES_MEV["s"] / DOWN_QUARK_MASSES_MEV["d"]
        obs_b_s = DOWN_QUARK_MASSES_MEV["b"] / DOWN_QUARK_MASSES_MEV["s"]
        obs_b_d = DOWN_QUARK_MASSES_MEV["b"] / DOWN_QUARK_MASSES_MEV["d"]

        out["s/d"] = _row(self.lepton_ratio(2, 1), obs_s_d)
        out["b/s"] = _row(self.lepton_ratio(3, 2), obs_b_s)
        out["b/d"] = _row(self.lepton_ratio(3, 1), obs_b_d)

        return out

    # ----------------------------------------------------------- drag mechanism
    #
    # Underlying mechanism for ALL fermion masses:
    #     m = hbar * omega_b(gamma_gen, particle_type) / c^2
    # where omega_b is the cone-bouncing angular frequency in the stiff
    # medium, set by the substrate drag coefficient gamma at each
    # generation. Mass ratios m_mu/m_e, m_tau/m_e are RELATIONS between
    # drag-derived masses --- the exp(n_M / (16 pi)) (mu/e) and
    # exp(n_M/(K_pair^4 pi) - (n_R-R)/(K_pair(K_pair+1))) (tau/mu) forms
    # ARE the inter-generation drag scale ratios.
    #
    # Anchor: gamma_1 chosen so generation-1 lepton mass = m_e (PDG).
    # Then gamma_2 = gamma_1 * exp(n_M/(K_pair^4 * pi))
    #      gamma_3 = gamma_2 * exp(n_M/(K_pair^4 pi) - (n_R-R)/(K_pair(K_pair+1)))
    #
    # Particle-type prefactor f(particle_type) absorbs sector-dependent
    # coupling (1.0 for charged leptons by definition; quarks/neutrinos
    # carry sector-specific multiplicative factors handled elsewhere).

    def gen_drag_scale(self, generation: int) -> float:
        """
        Substrate drag coefficient gamma_gen at the given generation.

        Returned in MeV-equivalent units: gamma_gen ~ hbar * omega_b / c^2
        for a unit particle-type prefactor. Anchored such that
        gen_drag_scale(1) reproduces m_e exactly.
        """
        if generation not in (1, 2, 3):
            raise ValueError(f"generation must be 1, 2, or 3 (got {generation})")
        gamma_1 = self.electron_mass_anchor()
        if generation == 1:
            return gamma_1
        if generation == 2:
            return gamma_1 * float(np.exp(self.n_M / self._denominator_12()))
        # generation == 3
        gamma_2 = gamma_1 * float(np.exp(self.n_M / self._denominator_12()))
        return gamma_2 * float(np.exp(self._log_ratio_23()))

    def electron_mass_anchor(self) -> float:
        """
        Anchor the generation-1 drag scale so that
        mass_from_drag(1, 'lepton') = m_e (PDG).
        Returns gamma_1 in MeV-equivalent units.
        """
        return LEPTON_MASSES_MEV["e"]

    def mass_from_drag(
        self, generation: int, particle_type: str = "lepton"
    ) -> float:
        """
        Mass from the drag-cone-bouncing mechanism:
            m = hbar * omega_b(gamma_gen, particle_type) / c^2
        Implemented as gamma_gen * f(particle_type) in MeV.

        particle_type: 'lepton' (f=1, anchored), 'up_quark', 'down_quark',
        or 'neutrino'. Non-lepton prefactors are placeholders carrying
        only the sector ratio --- absolute neutrino/quark scales are not
        the subject of this drag formulation here.
        """
        gamma = self.gen_drag_scale(generation)
        prefactors = {
            "lepton": 1.0,
            "up_quark": UP_QUARK_MASSES_MEV["u"] / LEPTON_MASSES_MEV["e"],
            "down_quark": DOWN_QUARK_MASSES_MEV["d"] / LEPTON_MASSES_MEV["e"],
            "neutrino": 1.0e-9,  # rough placeholder, eV-scale
        }
        if particle_type not in prefactors:
            raise ValueError(
                f"unknown particle_type {particle_type!r}; "
                f"valid: {list(prefactors)}"
            )
        return gamma * prefactors[particle_type]

    def predict_lepton_masses_from_drag(self) -> Dict[str, float]:
        """
        Returns the three charged-lepton masses (MeV) derived from the
        drag scaling, anchored at m_e.
        """
        return {
            "e":   self.mass_from_drag(1, "lepton"),
            "mu":  self.mass_from_drag(2, "lepton"),
            "tau": self.mass_from_drag(3, "lepton"),
        }

    def verify_drag_consistency(self) -> Dict[str, Dict[str, float]]:
        """
        Cross-check: drag-derived ratios m_mu/m_e and m_tau/m_e must
        equal the direct exp(n_M/...) lepton_ratio values. If not, the
        drag formulation is internally inconsistent.
        """
        masses = self.predict_lepton_masses_from_drag()
        drag_mu_e = masses["mu"] / masses["e"]
        drag_tau_e = masses["tau"] / masses["e"]
        drag_tau_mu = masses["tau"] / masses["mu"]
        return {
            "mu/e":   _row(drag_mu_e,   self.lepton_ratio(2, 1)),
            "tau/mu": _row(drag_tau_mu, self.lepton_ratio(3, 2)),
            "tau/e":  _row(drag_tau_e,  self.lepton_ratio(3, 1)),
        }

    def report_drag_scales(self) -> str:
        """
        Print the substrate drag coefficient gamma_gen at each generation
        with substrate-ontology interpretation.
        """
        lines = []
        lines.append("=" * 72)
        lines.append("Substrate Drag Scales --- gamma_gen at each generation")
        lines.append("  gamma = hbar * omega_b(cone-bounce) / c^2")
        lines.append(
            "  Anchor: gamma_1 := m_e (PDG) so charged-lepton masses are "
            "drag-derived"
        )
        lines.append("=" * 72)
        labels = {1: "electron", 2: "muon", 3: "tau"}
        for gen in (1, 2, 3):
            gamma = self.gen_drag_scale(gen)
            lines.append(
                f"  gen {gen} ({labels[gen]:<8}): gamma = {gamma:>14.6f} MeV"
            )
        lines.append("")
        lines.append("Inter-generation drag ratios (substrate topology):")
        lines.append(
            f"  gamma_2/gamma_1 = exp(n_M/(K_pair^4 pi)) = "
            f"{self.lepton_ratio(2,1):.6f}"
        )
        lines.append(
            f"  gamma_3/gamma_2 = exp(n_M/(K_pair^4 pi) - (n_R-R)/(K_pair(K_pair+1))) = "
            f"{self.lepton_ratio(3,2):.6f}"
        )
        lines.append("=" * 72)
        text = "\n".join(lines)
        print(text)
        return text

    # ---------------------------------------------------------------- report

    def report(self) -> str:
        """Build and return a comparison-table string for all SM fermions."""
        lines = []
        lines.append("=" * 72)
        lines.append("Substrate Generation Map — predicted vs observed")
        lines.append(
            f"  n_M={self.n_M}, K_pair={self.K_pair}, K_rank={self.K_rank}"
        )
        lines.append(
            f"  exp(n_M/(K_pair^4 * pi)) = {self.lepton_ratio(2,1):.6f}"
        )
        lines.append(
            f"  exp(n_M/(K_pair^4 pi) - (n_R-R)/(K_pair(K_pair+1))) = "
            f"{self.lepton_ratio(3,2):.6f}"
        )
        lines.append("=" * 72)

        sections = [
            ("Charged leptons", self.verify_lepton_masses()),
            ("Neutrinos (NH)", self.predict_neutrino_ratios("NH")),
            ("Neutrinos (IH)", self.predict_neutrino_ratios("IH")),
            ("Quarks", self.predict_quark_ratios()),
        ]

        header = f"{'ratio':<22}{'predicted':>14}{'observed':>16}{'err %':>10}"
        for name, results in sections:
            lines.append("")
            lines.append(f"-- {name} --")
            lines.append(header)
            for label, row in results.items():
                lines.append(
                    f"{label:<22}{row['predicted']:>14.4f}"
                    f"{row['observed']:>16.4f}{row['error_pct']:>10.3f}"
                )

        lines.append("=" * 72)
        text = "\n".join(lines)
        print(text)
        return text


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _row(predicted: float, observed: float) -> Dict[str, float]:
    """Build a {predicted, observed, error_pct} record."""
    error_pct = 100.0 * (predicted - observed) / observed
    return {
        "predicted": float(predicted),
        "observed": float(observed),
        "error_pct": float(error_pct),
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    GenerationMap().report()
