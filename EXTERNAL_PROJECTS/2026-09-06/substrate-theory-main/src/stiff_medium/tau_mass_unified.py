"""
tau_mass_unified.py
===================

Single source of truth for the m_tau/m_mu lepton-tower predictor.

Background
----------
Several modules independently implemented the second-to-third generation
lepton mass ratio:

    * generation_map.GenerationMap._log_ratio_23
    * mass_torque_engine.MassTorque._t_tau
    * integer_rigidity.predict_m_tau_over_m_e

A reconciliation audit (May 2026) confirmed that all three modules
already share the same n_R-active form

        log(m_tau/m_mu) = n_M/(K_pair^4 * pi)
                          - (n_R - R)/(K_pair*(K_pair+1))

at canonical integers (n_M, K_pair, K_rank, n_R, R) = (268, 2, 5, 18, 3),
giving log = 2.83169 -> ratio = 16.9741 (+0.93% above PDG 16.8170).

The supposed "1.18% drift between modules" reported in the audit prompt
was a misreading of older docstring text quoting the legacy form
``exp(n_M/(K_rank*(K_rank+1)*pi))`` (which is +2.14% off and only
appears in deprecated ``_denominator_23`` docstrings, never in active
code paths).

Purpose of this module
----------------------
1. Provide a single canonical function ``log_m_tau_over_m_mu`` that
   every consumer can import, eliminating any future drift risk.
2. Provide a ``reconcile_modules`` helper that asserts all three
   downstream implementations agree to machine precision (1e-12).
3. Document the two candidate NLO corrections that drive the 0.93%
   residual to <0.005%, while being explicit that neither has yet been
   independently derived from substrate first principles.

Substrate first-principles status
---------------------------------
The canonical formula is *forced* from two ingredients:

    * Leading exponent ``n_M / (K_pair^4 * pi)``:
          The same operator that lifts m_e -> m_mu (settled at 0.009%).
          ``K_pair^4 = 16`` is the Mobius double-cover (sheet count K_pair=2)
          to the 4th power = enclosing area scale of the bipartite K_4
          braided-anchor graph.

    * Subtractive rank-3 lock ``K_rank/K_pair`` = ``(n_R-R)/(K_pair*(K_pair+1))``:
          The rank-3 mode lock removes (n_R - R) = 15 reflection orbits
          from the available pool, normalised by the K_pair*(K_pair+1) = 6
          Mobius bundle area (fundamental + boundary). The identity
          K_rank/K_pair = (n_R-R)/(K_pair*(K_pair+1)) = 15/6 = 5/2 makes
          n_R an active observable input rather than a derived identity.

The +0.93% residual is therefore the HONEST precision of the
substrate first-principles lepton-tower formula at the current level
of derivation. Closing it requires either:

    (a) deriving the next-to-leading-order substrate correction from
        the cone-bouncing Lagrangian (in progress, see cone_bouncing_*),
        or
    (b) demonstrating that one of the candidate NLO terms below is
        forced by an independent topological argument.

Until (a) or (b) lands, this module exposes the canonical formula at
0.93% AND flags the two candidate corrections as research candidates,
not adopted predictions.

Candidate NLO corrections (research / not adopted)
--------------------------------------------------
Two single-term corrections drive the residual from 0.93% to <0.005%:

    (A) ``-K_rank/(K_pair * n_M)`` = -5/536 = -0.009328
        => log = 2.82236, pred = 16.8165, err = -0.003%
        Reading: rank-K_rank modes leak with weight 1/n_M per Mobius sheet
        (factor 1/K_pair). This is a *plausible* mass-shell NLO correction
        but has not been derived independently.

    (B) ``-1/(N_BAM * n_R)`` = -1/108 = -0.009259
        => log = 2.82243, pred = 16.8177, err = +0.004%
        Reading: a single rank-3 leakage mode through the
        (N_BAM hex-coord x n_R Mobius-reflection) = 108 channel pool.
        Equivalent to ``-1/(K_pair^2 * R^3)`` since 108 = 4 * 27.
        Also a plausible substrate identity, also not independently derived.

Either correction would close the residual, and there is currently no
substrate principle to choose between them. The honest stance is that
the canonical formula is +0.93% off PDG, full stop.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional

from .b3_constants import (
    n_M as _N_M,
    K_pair as _K_PAIR,
    K_rank as _K_RANK,
    n_R as _N_R,
    R as _R,
    N_BAM as _N_BAM,
)


# ---------------------------------------------------------------------------
# Observed PDG anchors (used only for reporting / verification)
# ---------------------------------------------------------------------------

M_TAU_MEV: float = 1776.86
M_MUON_MEV: float = 105.6583755

OBS_TAU_OVER_MU: float = M_TAU_MEV / M_MUON_MEV   # 16.8170293...
LOG_OBS_TAU_OVER_MU: float = math.log(OBS_TAU_OVER_MU)


# ---------------------------------------------------------------------------
# Canonical formula
# ---------------------------------------------------------------------------

def log_m_tau_over_m_mu(
    n_M: int = _N_M,
    K_pair: int = _K_PAIR,
    n_R: int = _N_R,
    R: int = _R,
) -> float:
    """Canonical n_R-active log(m_tau/m_mu) at substrate first principles.

        log(m_tau/m_mu) = n_M/(K_pair**4 * pi)
                          - (n_R - R)/(K_pair*(K_pair + 1))

    At canonical integers (268, 2, 18, 3) this yields 2.83169
    => m_tau/m_mu = 16.974 (+0.93% above PDG 16.817).

    Notes
    -----
    The reflection-counted second term equals K_rank/K_pair = 5/2
    EXACTLY because (n_R - R)/(K_pair*(K_pair+1)) = 15/6 = 5/2 at
    canonical values. Use the n_R-active form so that a perturbation
    in n_R produces an observable shift (rather than entering only
    via the n_M = K_pair * K_rank**3 + n_R defining identity).

    Parameters
    ----------
    n_M, K_pair, n_R, R : int
        Substrate integers (defaults from b3_constants).

    Returns
    -------
    float
        log(m_tau/m_mu) (canonical; +0.93% above observed).
    """
    return (n_M / (K_pair ** 4 * math.pi)
            - (n_R - R) / (K_pair * (K_pair + 1)))


def m_tau_over_m_mu(
    n_M: int = _N_M,
    K_pair: int = _K_PAIR,
    n_R: int = _N_R,
    R: int = _R,
) -> float:
    """m_tau/m_mu using the canonical substrate formula.

    See :func:`log_m_tau_over_m_mu` for derivation and precision.
    """
    return math.exp(log_m_tau_over_m_mu(n_M=n_M, K_pair=K_pair, n_R=n_R, R=R))


# ---------------------------------------------------------------------------
# Research-candidate NLO corrections (NOT adopted)
# ---------------------------------------------------------------------------

def log_m_tau_over_m_mu_with_nlo_A(
    n_M: int = _N_M,
    K_pair: int = _K_PAIR,
    K_rank: int = _K_RANK,
    n_R: int = _N_R,
    R: int = _R,
) -> float:
    """RESEARCH CANDIDATE: canonical + NLO_A correction.

    Adds ``-K_rank/(K_pair * n_M)`` to the canonical exponent. Drives
    the residual from +0.93% to -0.003% but lacks independent
    substrate-first-principles derivation. Provided for audit /
    sensitivity analysis, NOT for use as a published prediction.
    """
    return log_m_tau_over_m_mu(n_M=n_M, K_pair=K_pair, n_R=n_R, R=R) \
        - K_rank / (K_pair * n_M)


def log_m_tau_over_m_mu_with_nlo_B(
    n_M: int = _N_M,
    K_pair: int = _K_PAIR,
    n_R: int = _N_R,
    R: int = _R,
    N_BAM: int = _N_BAM,
) -> float:
    """RESEARCH CANDIDATE: canonical + NLO_B correction.

    Adds ``-1/(N_BAM * n_R)`` to the canonical exponent. Drives the
    residual from +0.93% to +0.004% but lacks independent
    substrate-first-principles derivation. Numerically equivalent to
    ``-1/(K_pair**2 * R**3)`` since 108 = 4*27. Provided for audit /
    sensitivity analysis, NOT for use as a published prediction.
    """
    return log_m_tau_over_m_mu(n_M=n_M, K_pair=K_pair, n_R=n_R, R=R) \
        - 1.0 / (N_BAM * n_R)


# ---------------------------------------------------------------------------
# Reconciliation helpers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TauReconciliation:
    """Side-by-side comparison of every downstream m_tau/m_mu predictor."""
    canonical: float            # log from this module
    generation_map: float       # log from GenerationMap._log_ratio_23
    mass_torque_engine: float   # log from MassTorque._t_tau (extracted)
    integer_rigidity: float     # log from predict_m_tau_over_m_e (extracted)
    pdg_observed: float         # log of M_TAU/M_MUON
    max_pairwise_diff: float    # max |a - b| across the 4 module logs

    def all_modules_agree(self, tol: float = 1e-12) -> bool:
        """True if all modules agree to within the given tolerance."""
        return self.max_pairwise_diff <= tol

    def err_pct_vs_pdg(self) -> float:
        """Percent error of the canonical prediction vs PDG."""
        pred = math.exp(self.canonical)
        obs = math.exp(self.pdg_observed)
        return 100.0 * (pred - obs) / obs


def reconcile_modules() -> TauReconciliation:
    """Cross-check every downstream m_tau/m_mu predictor.

    Imports each module's predictor, evaluates it at the canonical
    integers, extracts log(m_tau/m_mu) consistently from each, and
    returns a TauReconciliation snapshot. Use ``all_modules_agree``
    on the result to assert cross-module consistency to 1e-12.
    """
    # Canonical (this module)
    log_canonical = log_m_tau_over_m_mu()

    # GenerationMap returns the ratio (not log), and exposes
    # _log_ratio_23 directly.
    from .generation_map import GenerationMap
    gm = GenerationMap()
    log_gm = gm._log_ratio_23()

    # MassTorque._t_tau returns torque, formula, used_prims, used_ints.
    # We reconstruct log(m_tau/m_mu) from m_tau and m_mu values.
    from .mass_torque_engine import MassTorque, M_ELECTRON_MEV
    eng = MassTorque()
    tau_result = eng.compute("tau")
    muon_result = eng.compute("muon")
    log_mt = math.log(tau_result.value_mev / muon_result.value_mev)

    # integer_rigidity exposes m_tau/m_e and m_mu/m_e.
    from .integer_rigidity import (
        predict_m_tau_over_m_e,
        predict_m_mu_over_m_e,
        CANONICAL_INTEGERS,
    )
    tau_e = predict_m_tau_over_m_e(CANONICAL_INTEGERS)
    mu_e = predict_m_mu_over_m_e(CANONICAL_INTEGERS)
    log_ir = math.log(tau_e / mu_e)

    logs = [log_canonical, log_gm, log_mt, log_ir]
    max_diff = max(abs(a - b) for a in logs for b in logs)

    return TauReconciliation(
        canonical=log_canonical,
        generation_map=log_gm,
        mass_torque_engine=log_mt,
        integer_rigidity=log_ir,
        pdg_observed=LOG_OBS_TAU_OVER_MU,
        max_pairwise_diff=max_diff,
    )


# ---------------------------------------------------------------------------
# CLI / quick-look
# ---------------------------------------------------------------------------

def report() -> str:
    """Compact reconciliation + PDG comparison report."""
    rec = reconcile_modules()
    pred = math.exp(rec.canonical)
    obs = math.exp(rec.pdg_observed)
    err = 100.0 * (pred - obs) / obs

    lines = [
        "=" * 72,
        "tau_mass_unified -- canonical m_tau/m_mu reconciliation",
        "=" * 72,
        f"  canonical         log = {rec.canonical:.12f}",
        f"  generation_map    log = {rec.generation_map:.12f}",
        f"  mass_torque       log = {rec.mass_torque_engine:.12f}",
        f"  integer_rigidity  log = {rec.integer_rigidity:.12f}",
        f"  PDG observed      log = {rec.pdg_observed:.12f}",
        f"  max pairwise diff     = {rec.max_pairwise_diff:.3e}",
        "",
        f"  predicted m_tau/m_mu = {pred:.6f}",
        f"  observed  m_tau/m_mu = {obs:.6f}",
        f"  err vs PDG           = {err:+.4f} %",
        "",
        "Research candidate NLO corrections (NOT adopted):",
        f"  NLO_A: -K_rank/(K_pair*n_M) -> err = "
        f"{100*(math.exp(log_m_tau_over_m_mu_with_nlo_A())-obs)/obs:+.4f} %",
        f"  NLO_B: -1/(N_BAM*n_R)       -> err = "
        f"{100*(math.exp(log_m_tau_over_m_mu_with_nlo_B())-obs)/obs:+.4f} %",
        "=" * 72,
    ]
    text = "\n".join(lines)
    print(text)
    return text


if __name__ == "__main__":
    report()


__all__ = [
    "M_TAU_MEV",
    "M_MUON_MEV",
    "OBS_TAU_OVER_MU",
    "LOG_OBS_TAU_OVER_MU",
    "log_m_tau_over_m_mu",
    "m_tau_over_m_mu",
    "log_m_tau_over_m_mu_with_nlo_A",
    "log_m_tau_over_m_mu_with_nlo_B",
    "TauReconciliation",
    "reconcile_modules",
    "report",
]
