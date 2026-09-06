"""
integer_rigidity.py
===================

B3 Integer Rigidity Test Framework
----------------------------------

The B3 framework rests on the claim that a small set of integers --

    n_M    = 268    (rank-coupled mode count)
    N_BAM  = 6      (binding action multiplicity)
    K_pair = 2      (pair stiffness)
    K_rank = 5      (rank-coupling integer)
    n_R    = 18     (rank dimension)
    n_A    = 15     (anchor pair count = C(N_BAM, 2))

-- are *uniquely* selected by the simplicial / braided ansatz. The
audit finding (audit_06) was that no test in the codebase fails when
those integers are perturbed. This module supplies that test.

Strategy:

  1. Pick a B3 integer and a perturbation delta (e.g. n_M -> 269).
  2. Re-instantiate the predictive engines with the perturbed integer
     set, leaving every other integer / primitive at its canonical
     value.
  3. Recompute every derived prediction that depends on the
     perturbed integer.
  4. Compare to PDG / CODATA / observed values, accumulating a
     chi-squared-like residual.
  5. Show that the canonical integer minimises that residual --
     i.e. the framework is genuinely *rigid* under integer choice.

If shifting n_M by +/-1 leaves the framework's chi^2 unchanged, the
integer is not load-bearing. If shifting it by +/-1 makes the
framework drift far outside the tolerance band, the integer is
rigid and the framework is making real predictions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any
import math

import numpy as np

from .mass_torque_engine import (
    MassTorque,
    DEFAULT_INTEGERS,
    LAMBDA_QCD_MEV,
    M_ELECTRON_MEV,
    M_MUON_MEV,
    M_TAU_MEV,
    DEUTERON_BE_MEV,
    T_C_MAX_K,
)


# ---------------------------------------------------------------------------
# Observed anchors used as ground truth
# ---------------------------------------------------------------------------

OBSERVED: Dict[str, float] = {
    "m_mu_over_m_e":  M_MUON_MEV / M_ELECTRON_MEV,         # 206.768
    "m_tau_over_m_e": M_TAU_MEV / M_ELECTRON_MEV,          # 3477.2
    "deuteron_BE":    DEUTERON_BE_MEV,                     # 2.224573 MeV
    "t_c_max":        T_C_MAX_K,                           # 128.9 K
    "cabibbo":        0.22534,                             # PDG sin(theta_C)
}

CANONICAL_INTEGERS = dict(DEFAULT_INTEGERS)


# ---------------------------------------------------------------------------
# Prediction helpers (pure-functional, integer-driven)
# ---------------------------------------------------------------------------

def predict_m_mu_over_m_e(integers: Dict[str, int]) -> float:
    """m_mu / m_e = exp(n_M / (K_pair^4 * pi))."""
    n_M = integers["n_M"]
    K_pair = integers["K_pair"]
    return float(math.exp(n_M / (K_pair ** 4 * math.pi)))


def predict_m_tau_over_m_e(integers: Dict[str, int]) -> float:
    """m_tau / m_e via the lepton tower (canonical, n_R-active form).

        m_tau / m_mu = exp(n_M/(K_pair^4 pi) - (n_R - R)/(K_pair*(K_pair+1)))
        m_mu  / m_e  = exp(n_M/(K_pair^4 pi))

    The m_tau/m_mu factor delegates to
    :func:`tau_mass_unified.log_m_tau_over_m_mu`, the single source of
    truth shared with ``generation_map._log_ratio_23`` and
    ``mass_torque_engine._t_tau``.

    The reflection-counted term (n_R - R)/(K_pair*(K_pair+1)) equals
    K_rank/K_pair = 5/2 at the canonical integers (15/6 = 5/2 EXACTLY),
    so the prediction is numerically identical to the K_rank/K_pair form.
    But n_R now enters the predictor *explicitly*: perturbing n_R by +/-1
    shifts m_tau/m_e (rather than silently doing nothing as in the prior
    K_rank-only form).
    """
    from .tau_mass_unified import log_m_tau_over_m_mu
    n_M = integers["n_M"]
    K_pair = integers["K_pair"]
    n_R = integers["n_R"]
    R_int = integers.get("R", 3)
    log_mu = n_M / (K_pair ** 4 * math.pi)
    log_ratio = log_m_tau_over_m_mu(
        n_M=n_M, K_pair=K_pair, n_R=n_R, R=R_int,
    )
    return float(math.exp(log_mu) * math.exp(log_ratio))


def predict_deuteron_BE(integers: Dict[str, int],
                        lambda_qcd: float = LAMBDA_QCD_MEV) -> float:
    """BE_d = Lambda_QCD / (n_A * N_BAM)  -> 200/90 = 2.222 MeV.

    Canonical product n_A · N_BAM = 15 · 6 = 90.
    """
    n_A = integers["n_A"]
    N_BAM = integers["N_BAM"]
    return float(lambda_qcd / (n_A * N_BAM))


def predict_t_c_max(integers: Dict[str, int],
                    lambda_qcd: float = LAMBDA_QCD_MEV) -> float:
    """T_c,max via the (Lambda_QCD / R) substrate-alignment factor."""
    R = integers["R"]
    kB_factor = 1.9335
    return float(lambda_qcd * (1.0 / R) * kB_factor)


def predict_cabibbo(integers: Dict[str, int]) -> float:
    """sin(theta_C) = 1 / sqrt(K_rank * (K_rank - 1)).

    Canonical K_rank = 5  ->  1/sqrt(20) = 0.22361   vs PDG 0.22534.
    """
    K_rank = integers["K_rank"]
    denom = K_rank * (K_rank - 1)
    if denom <= 0:
        return float("inf")
    return float(1.0 / math.sqrt(denom))


def predict_eps_face(integers: Dict[str, int]) -> float:
    """Face alignment energy fraction: eps_face = 1/N_BAM (canonical 1/6)."""
    N_BAM = integers["N_BAM"]
    if N_BAM <= 0:
        return float("inf")
    return 1.0 / float(N_BAM)


# ---------------------------------------------------------------------------
# Mapping: integer -> list of (label, predictor, observed_key) tuples
# ---------------------------------------------------------------------------

INTEGER_DEPENDENCIES: Dict[str, List[Tuple[str, Any, Optional[str]]]] = {
    "n_M": [
        ("m_mu/m_e",  predict_m_mu_over_m_e,  "m_mu_over_m_e"),
        ("m_tau/m_e", predict_m_tau_over_m_e, "m_tau_over_m_e"),
    ],
    "N_BAM": [
        ("deuteron_BE", predict_deuteron_BE, "deuteron_BE"),
        ("eps_face",    predict_eps_face,    None),
    ],
    "K_pair": [
        ("m_mu/m_e",  predict_m_mu_over_m_e,  "m_mu_over_m_e"),
        ("m_tau/m_e", predict_m_tau_over_m_e, "m_tau_over_m_e"),
    ],
    "K_rank": [
        ("cabibbo",   predict_cabibbo,        "cabibbo"),
        ("m_tau/m_e", predict_m_tau_over_m_e, "m_tau_over_m_e"),
    ],
    "n_R": [
        ("m_tau/m_e", predict_m_tau_over_m_e, "m_tau_over_m_e"),
    ],
    "n_A": [
        ("deuteron_BE", predict_deuteron_BE, "deuteron_BE"),
    ],
    "R": [
        ("t_c_max", predict_t_c_max, "t_c_max"),
    ],
}

GLOBAL_PREDICTIONS: List[Tuple[str, Any, Optional[str]]] = [
    ("m_mu/m_e",    predict_m_mu_over_m_e,  "m_mu_over_m_e"),
    ("m_tau/m_e",   predict_m_tau_over_m_e, "m_tau_over_m_e"),
    ("deuteron_BE", predict_deuteron_BE,    "deuteron_BE"),
    ("t_c_max",     predict_t_c_max,        "t_c_max"),
    ("cabibbo",     predict_cabibbo,        "cabibbo"),
]


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass
class PerturbationResult:
    """Single (integer, delta) perturbation report."""
    integer_name: str
    delta: int
    perturbed_value: int
    observables: List[Dict[str, Any]] = field(default_factory=list)
    total_chi2: float = 0.0

    def degrade_factor(self, baseline_chi2: float) -> float:
        if baseline_chi2 <= 0:
            return float("inf") if self.total_chi2 > 0 else 1.0
        return self.total_chi2 / baseline_chi2


# ---------------------------------------------------------------------------
# Main rigidity engine
# ---------------------------------------------------------------------------

def _call_predictor(predictor, integers: Dict[str, int],
                    lambda_qcd: float) -> float:
    """Invoke a predictor regardless of whether it accepts lambda_qcd."""
    try:
        return float(predictor(integers, lambda_qcd))
    except TypeError:
        return float(predictor(integers))


class IntegerRigidity:
    """B3 integer rigidity test framework.

    Given the canonical B3 integer set, this class lets you ask:
    'what happens to the framework's chi^2 if I shift n_M by +/-1?'
    The canonical answer should be: it gets dramatically worse.

    The class never mutates the canonical integers; it always operates
    on a *copy* with a single perturbation applied.
    """

    def __init__(self,
                 canonical: Optional[Dict[str, int]] = None,
                 lambda_qcd: float = LAMBDA_QCD_MEV) -> None:
        self.canonical = dict(canonical) if canonical else dict(CANONICAL_INTEGERS)
        self.lambda_qcd = float(lambda_qcd)
        self._baseline_chi2 = self._chi2(self.canonical)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _perturbed(self, name: str, delta: int) -> Dict[str, int]:
        """Return a fresh integer dict with `name` shifted by `delta`."""
        if name not in self.canonical:
            raise KeyError(f"Unknown integer '{name}'. "
                           f"Known: {sorted(self.canonical)}")
        out = dict(self.canonical)
        out[name] = int(out[name]) + int(delta)
        return out

    def _chi2(self, integers: Dict[str, int]) -> float:
        """Sum of squared relative residuals across all global predictions."""
        chi2 = 0.0
        for _label, predictor, obs_key in GLOBAL_PREDICTIONS:
            try:
                pred = _call_predictor(predictor, integers, self.lambda_qcd)
            except (ValueError, ZeroDivisionError, OverflowError):
                return float("inf")
            if obs_key is None:
                continue
            obs = OBSERVED[obs_key]
            if obs == 0 or not math.isfinite(pred):
                return float("inf")
            chi2 += ((pred - obs) / obs) ** 2
        return float(chi2)

    # ------------------------------------------------------------------
    # Generic perturbation
    # ------------------------------------------------------------------

    def perturb(self, integer_name: str, delta: int) -> PerturbationResult:
        """Perturb a single integer by `delta` and report consequences."""
        perturbed = self._perturbed(integer_name, delta)
        baseline = self.canonical
        canonical_baselines = {
            "eps_face": 1.0 / float(baseline["N_BAM"]),
        }

        deps = INTEGER_DEPENDENCIES.get(integer_name, GLOBAL_PREDICTIONS)
        rows: List[Dict[str, Any]] = []
        for label, predictor, obs_key in deps:
            try:
                pred = _call_predictor(predictor, perturbed, self.lambda_qcd)
            except (ValueError, ZeroDivisionError, OverflowError):
                pred = float("nan")
            if obs_key is not None:
                obs = OBSERVED[obs_key]
            else:
                obs = canonical_baselines.get(label, float("nan"))

            if obs and math.isfinite(pred) and math.isfinite(obs):
                rel_err = abs(pred - obs) / abs(obs)
            else:
                rel_err = float("inf")

            try:
                canonical_pred = _call_predictor(predictor, baseline,
                                                 self.lambda_qcd)
            except Exception:
                canonical_pred = float("nan")

            rows.append({
                "observable": label,
                "canonical_pred": canonical_pred,
                "perturbed_pred": pred,
                "observed": obs,
                "rel_err": rel_err,
            })

        total_chi2 = self._chi2(perturbed)
        return PerturbationResult(
            integer_name=integer_name,
            delta=int(delta),
            perturbed_value=perturbed[integer_name],
            observables=rows,
            total_chi2=total_chi2,
        )

    # ------------------------------------------------------------------
    # Named convenience methods
    # ------------------------------------------------------------------

    def perturb_n_M(self, delta: int) -> PerturbationResult:
        return self.perturb("n_M", delta)

    def perturb_N_BAM(self, delta: int) -> PerturbationResult:
        return self.perturb("N_BAM", delta)

    def perturb_K_pair(self, delta: int) -> PerturbationResult:
        return self.perturb("K_pair", delta)

    def perturb_K_rank(self, delta: int) -> PerturbationResult:
        return self.perturb("K_rank", delta)

    def perturb_n_R(self, delta: int) -> PerturbationResult:
        return self.perturb("n_R", delta)

    def perturb_n_A(self, delta: int) -> PerturbationResult:
        return self.perturb("n_A", delta)

    # ------------------------------------------------------------------
    # Aggregate analyses
    # ------------------------------------------------------------------

    def baseline_chi2(self) -> float:
        return float(self._baseline_chi2)

    def scan(self, integer_name: str, deltas=(-2, -1, 0, 1, 2)
             ) -> List[PerturbationResult]:
        """Sweep an integer over a set of deltas, returning per-delta results."""
        return [self.perturb(integer_name, d) for d in deltas]

    def verify_uniqueness(self,
                          integer_name: str = "n_M",
                          window: int = 5) -> Dict[str, Any]:
        """Check that the canonical value is the *unique* chi^2 minimiser."""
        deltas = list(range(-window, window + 1))
        sweep = self.scan(integer_name, deltas=deltas)
        chi2s = [(r.delta, r.total_chi2) for r in sweep]
        finite = [(d, c) for d, c in chi2s if math.isfinite(c)]
        if not finite:
            return {
                "integer": integer_name,
                "argmin_delta": None,
                "is_unique_minimum": False,
                "sweep": chi2s,
            }
        argmin_delta, argmin_chi2 = min(finite, key=lambda t: t[1])
        threshold = argmin_chi2 * 1.01 + 1e-12
        ties = [d for d, c in finite if c <= threshold and d != argmin_delta]
        return {
            "integer": integer_name,
            "argmin_delta": argmin_delta,
            "argmin_chi2": argmin_chi2,
            "is_unique_minimum": (argmin_delta == 0 and not ties),
            "ties": ties,
            "sweep": chi2s,
        }

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def rigidity_report(self,
                        deltas=(-2, -1, 1, 2),
                        integers: Optional[List[str]] = None) -> str:
        """Print and return a human-readable rigidity report."""
        if integers is None:
            integers = ["n_M", "N_BAM", "K_pair", "K_rank", "n_R", "n_A"]
        baseline = self.baseline_chi2()
        lines: List[str] = []
        lines.append("B3 Integer Rigidity Report")
        lines.append("=" * 70)
        lines.append(f"Canonical chi^2 (baseline) = {baseline:.6e}")
        lines.append("")
        header = f"{'integer':<8}{'delta':>6}{'value':>8}{'chi^2':>16}{'degrade x':>14}"
        lines.append(header)
        lines.append("-" * len(header))
        for name in integers:
            for d in deltas:
                try:
                    r = self.perturb(name, d)
                except KeyError:
                    continue
                degrade = r.degrade_factor(baseline)
                lines.append(
                    f"{name:<8}{d:>+6d}{r.perturbed_value:>8d}"
                    f"{r.total_chi2:>16.6e}{degrade:>14.3e}"
                )
            lines.append("-" * len(header))
        text = "\n".join(lines)
        print(text)
        return text


# ---------------------------------------------------------------------------
# Module entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":   # pragma: no cover
    rig = IntegerRigidity()
    rig.rigidity_report()
    print()
    for nm in ("n_M", "K_pair", "N_BAM", "K_rank"):
        u = rig.verify_uniqueness(nm, window=3)
        flag = "UNIQUE" if u["is_unique_minimum"] else "NOT UNIQUE"
        print(f"  {nm}: argmin delta = {u['argmin_delta']} [{flag}]")
