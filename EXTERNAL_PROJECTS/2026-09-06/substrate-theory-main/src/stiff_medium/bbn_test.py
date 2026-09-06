"""BBN primordial-abundance comparison test against PDG 2024 / Planck 2018.

Wraps `bbn_substrate.BBNSubstrate` plus the eta-fit formulas exposed by
`bbn_predictions` to produce a single comparison object containing:

  * Y_p (He-4 mass fraction)
  * D/H
  * He-3/H
  * Li-7/H

against the empirical central values + 1-sigma uncertainties from
PDG 2024 BBN review and the Planck 2018 baryon-to-photon ratio
eta = 6.10e-10.

The substrate framework's substantive claim about BBN is twofold:

  1. The standard nuclear network (Y_p, D/H, He-3/H) is INHERITED — the
     substrate predicts the same numbers as SBBN because the relevant
     physics is stable nuclear-mass / weak freeze-out, which the substrate
     reproduces via §18.62-§18.63 (m_p, Δm) and §18.40-§18.44 (FRW H).

  2. For the Li-7 puzzle, the substrate provides a late-time Be-7 destruction
     channel via substrate de-saturation re-thermalization at the observer
     horizon (§18.66 proto-matter window), suppressing the final Li-7
     abundance by a factor consistent with the Spite-plateau observation.

This module does not re-derive numbers; it COMPARES substrate predictions
(piped through bbn_substrate + bbn_predictions eta-fits) to observations
and reports n-sigma deviations + a verdict on the Li-7 puzzle.

References:
  Particle Data Group 2024, Big Bang Nucleosynthesis review (Fields & Olive).
  Planck Collaboration 2018, A&A 641 A6 — eta = 6.10e-10.
  Aver, Olive, Skillman 2015 — Y_p observational central value.
  Cooke, Pettini, Steidel 2018 — D/H from high-z DLAs.
  Sbordone et al. 2010 — Spite-plateau Li-7/H.
  Spec §18.75.4 (BBN in substrate framework), §18.66 (proto-matter).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from .bbn_substrate import BBNSubstrate, ETA_BARYON_PHOTON
from .bbn_predictions import (
    deuterium_abundance,
    helium3_abundance,
    lithium7_abundance,
)
from .substrate_li7_suppression import SUBSTRATE_LI7_SUPPRESSION_FACTOR


# ---------------------------------------------------------------------------
# Observational central values + 1σ uncertainties (PDG 2024 + Planck 2018)
# ---------------------------------------------------------------------------

#: Helium-4 mass fraction (Aver, Olive & Skillman 2015; consistent with EMPRESS 2022).
OBS_Y_P_CENTRAL: float = 0.245
OBS_Y_P_SIGMA: float = 0.003

#: Deuterium-to-hydrogen ratio × 10^5 (Cooke, Pettini, Steidel 2018, high-z DLAs).
OBS_D_H_X1E5_CENTRAL: float = 2.527
OBS_D_H_X1E5_SIGMA: float = 0.030

#: Helium-3 number ratio × 10^5 (interstellar / planetary nebulae).
OBS_HE3_H_X1E5_CENTRAL: float = 1.0
OBS_HE3_H_X1E5_SIGMA: float = 0.1

#: Li-7 number ratio × 10^10 (Spite plateau, metal-poor halo stars).
OBS_LI7_H_X1E10_CENTRAL: float = 1.6
OBS_LI7_H_X1E10_SIGMA: float = 0.3


# ---------------------------------------------------------------------------
# Per-observable result type
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BBNObservableFit:
    """Single primordial-abundance comparison row.

    Attributes:
        name:        Human-readable name (e.g. "Y_p").
        units:       String label of the units (e.g. "mass fraction", "× 1e-5").
        predicted:   Substrate / SBBN prediction in the natural units.
        observed:    Observed central value (same units).
        sigma:       1-σ observational uncertainty.
        rel_error:   |predicted - observed| / observed.
        n_sigma:     |predicted - observed| / sigma.
        within_2s:   True if n_sigma <= 2.
    """
    name: str
    units: str
    predicted: float
    observed: float
    sigma: float
    rel_error: float
    n_sigma: float
    within_2s: bool


@dataclass(frozen=True)
class BBNTestResult:
    """Full BBN observational fit summary.

    Attributes:
        eta:                    Baryon-to-photon ratio used.
        rows:                   {name -> BBNObservableFit}.
        li7_sbbn_prediction:    SBBN ⁷Li/H × 10^10 (no substrate suppression).
        li7_substrate_prediction: substrate ⁷Li/H × 10^10 (with suppression).
        li7_sbbn_n_sigma:       SBBN ⁷Li/H deviation from observation in σ.
        li7_substrate_n_sigma:  substrate ⁷Li/H deviation from observation in σ.
        li7_substrate_beats_sbbn: True iff |substrate-σ| < |SBBN-σ|.
        verdict:                One-line summary.
    """
    eta: float
    rows: Dict[str, BBNObservableFit]
    li7_sbbn_prediction_x1e10: float
    li7_substrate_prediction_x1e10: float
    li7_sbbn_n_sigma: float
    li7_substrate_n_sigma: float
    li7_substrate_beats_sbbn: bool
    verdict: str


# ---------------------------------------------------------------------------
# Comparison helper
# ---------------------------------------------------------------------------

def _compare(name: str,
             units: str,
             predicted: float,
             observed: float,
             sigma: float) -> BBNObservableFit:
    """Build a single BBNObservableFit row."""
    rel = abs(predicted - observed) / abs(observed) if observed != 0.0 else float("inf")
    nsig = abs(predicted - observed) / sigma if sigma > 0.0 else float("inf")
    return BBNObservableFit(
        name=name,
        units=units,
        predicted=predicted,
        observed=observed,
        sigma=sigma,
        rel_error=rel,
        n_sigma=nsig,
        within_2s=(nsig <= 2.0),
    )


# ---------------------------------------------------------------------------
# Top-level test
# ---------------------------------------------------------------------------

def run_bbn_test(eta: float = ETA_BARYON_PHOTON,
                 li7_substrate_suppression: float = SUBSTRATE_LI7_SUPPRESSION_FACTOR) -> BBNTestResult:
    """Run the full substrate-vs-observation BBN comparison.

    Args:
        eta: Baryon-to-photon ratio (default 6.10e-10, Planck 2018).
        li7_substrate_suppression: Substrate Be-7 destruction factor applied
            to the SBBN ⁷Li prediction (§18.66 proto-matter / observer-horizon
            re-thermalization).  Defaults to the substrate-derived value
            ``SUBSTRATE_LI7_SUPPRESSION_FACTOR = n_A / K_rank = 15 / 5 = 3``
            (face-pair adjacency channels per 4-simplex vertex; see
            :mod:`substrate_li7_suppression`).  Brings SBBN's 4.5e-10 down to
            ~1.5e-10, matching the Spite plateau.

    Returns:
        BBNTestResult containing per-observable rows + Li-7-puzzle verdict.
    """
    bbn = BBNSubstrate(eta=eta, li7_substrate_suppression=li7_substrate_suppression)

    yp_pred = bbn.Y_p_He4()                     # 0.247
    dh_pred_x1e5 = deuterium_abundance(eta) * 1.0e5
    he3_pred_x1e5 = helium3_abundance(eta) * 1.0e5
    li7_sbbn_x1e10 = lithium7_abundance(eta) * 1.0e10
    li7_sub_x1e10 = li7_sbbn_x1e10 / li7_substrate_suppression

    rows: Dict[str, BBNObservableFit] = {
        "Y_p": _compare(
            "Y_p (⁴He mass fraction)", "mass fraction",
            yp_pred, OBS_Y_P_CENTRAL, OBS_Y_P_SIGMA,
        ),
        "D/H": _compare(
            "D/H ratio", "× 1e-5",
            dh_pred_x1e5, OBS_D_H_X1E5_CENTRAL, OBS_D_H_X1E5_SIGMA,
        ),
        "He-3/H": _compare(
            "³He/H ratio", "× 1e-5",
            he3_pred_x1e5, OBS_HE3_H_X1E5_CENTRAL, OBS_HE3_H_X1E5_SIGMA,
        ),
        "Li-7/H": _compare(
            "⁷Li/H ratio (substrate-suppressed)", "× 1e-10",
            li7_sub_x1e10, OBS_LI7_H_X1E10_CENTRAL, OBS_LI7_H_X1E10_SIGMA,
        ),
    }

    sbbn_nsig = abs(li7_sbbn_x1e10 - OBS_LI7_H_X1E10_CENTRAL) / OBS_LI7_H_X1E10_SIGMA
    sub_nsig = rows["Li-7/H"].n_sigma
    li7_beats = sub_nsig < sbbn_nsig

    if li7_beats and rows["Li-7/H"].within_2s:
        li7_msg = (f"Li-7 puzzle resolved within 2σ (substrate {sub_nsig:.2f}σ "
                   f"vs SBBN {sbbn_nsig:.1f}σ)")
    elif li7_beats:
        li7_msg = (f"Li-7 puzzle improved (substrate {sub_nsig:.2f}σ "
                   f"vs SBBN {sbbn_nsig:.1f}σ) but not yet within 2σ")
    else:
        li7_msg = (f"substrate suppression does NOT improve Li-7 fit "
                   f"(substrate {sub_nsig:.2f}σ vs SBBN {sbbn_nsig:.1f}σ)")

    other_within_2s = all(rows[k].within_2s for k in ("Y_p", "D/H", "He-3/H"))
    if other_within_2s:
        verdict = f"Y_p/D/H/³He all within 2σ; {li7_msg}"
    else:
        offenders = [k for k in ("Y_p", "D/H", "He-3/H") if not rows[k].within_2s]
        verdict = f"{','.join(offenders)} >2σ; {li7_msg}"

    return BBNTestResult(
        eta=eta,
        rows=rows,
        li7_sbbn_prediction_x1e10=li7_sbbn_x1e10,
        li7_substrate_prediction_x1e10=li7_sub_x1e10,
        li7_sbbn_n_sigma=sbbn_nsig,
        li7_substrate_n_sigma=sub_nsig,
        li7_substrate_beats_sbbn=li7_beats,
        verdict=verdict,
    )


# ---------------------------------------------------------------------------
# eta scan: predicted abundances as a function of baryon-to-photon ratio
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EtaScanResult:
    """Predicted abundances vs. baryon-to-photon ratio η for plot panels.

    Attributes:
        eta_grid:         η values sampled.
        Y_p:              ⁴He mass fraction (constant — Y_p is η-insensitive
                          in the SBBN η range; substrate inherits this).
        D_H_x1e5:         D/H × 10^5 along eta_grid.
        He3_H_x1e5:       ³He/H × 10^5 along eta_grid.
        Li7_H_x1e10_sbbn: SBBN ⁷Li/H × 10^10 along eta_grid (no suppression).
        Li7_H_x1e10_sub:  substrate ⁷Li/H × 10^10 (with suppression).
        eta_planck:       Planck-fixed central η.
    """
    eta_grid: List[float]
    Y_p: List[float]
    D_H_x1e5: List[float]
    He3_H_x1e5: List[float]
    Li7_H_x1e10_sbbn: List[float]
    Li7_H_x1e10_sub: List[float]
    eta_planck: float = ETA_BARYON_PHOTON


def scan_eta(eta_min: float = 3.0e-10,
             eta_max: float = 1.0e-9,
             n_pts: int = 60,
             li7_substrate_suppression: float = SUBSTRATE_LI7_SUPPRESSION_FACTOR) -> EtaScanResult:
    """Sample predictions across η for the abundance-vs-η plot."""
    if eta_max <= eta_min:
        raise ValueError("eta_max must exceed eta_min")
    if n_pts < 2:
        raise ValueError("need >=2 points")

    grid: List[float] = []
    for i in range(n_pts):
        # Logarithmic spacing
        frac = i / (n_pts - 1)
        eta_i = eta_min * (eta_max / eta_min) ** frac
        grid.append(eta_i)

    bbn = BBNSubstrate(li7_substrate_suppression=li7_substrate_suppression)
    yp_const = bbn.Y_p_He4()  # SBBN value, η-insensitive

    yp_list = [yp_const for _ in grid]
    dh_list = [deuterium_abundance(e) * 1.0e5 for e in grid]
    he3_list = [helium3_abundance(e) * 1.0e5 for e in grid]
    li7_sbbn = [lithium7_abundance(e) * 1.0e10 for e in grid]
    li7_sub = [v / li7_substrate_suppression for v in li7_sbbn]

    return EtaScanResult(
        eta_grid=grid,
        Y_p=yp_list,
        D_H_x1e5=dh_list,
        He3_H_x1e5=he3_list,
        Li7_H_x1e10_sbbn=li7_sbbn,
        Li7_H_x1e10_sub=li7_sub,
    )


# ---------------------------------------------------------------------------
# Pretty report
# ---------------------------------------------------------------------------

def report(result: BBNTestResult) -> str:
    """Render a one-page text summary of the BBN test."""
    lines: List[str] = []
    lines.append("=" * 78)
    lines.append("BBN primordial-abundance test: substrate vs PDG 2024 / Planck 2018")
    lines.append("=" * 78)
    lines.append(f"  η = {result.eta:.3e}    (Planck 2018 CMB)")
    lines.append("")
    lines.append(
        f"  {'Observable':<35} {'predicted':>12} {'observed':>12}"
        f" {'σ':>8} {'n-σ':>6}"
    )
    lines.append("  " + "-" * 75)
    for k in ("Y_p", "D/H", "He-3/H", "Li-7/H"):
        r = result.rows[k]
        flag = "✓" if r.within_2s else "✗"
        lines.append(
            f"  {r.name:<35} {r.predicted:>12.4g} {r.observed:>12.4g}"
            f" {r.sigma:>8.3g} {r.n_sigma:>5.2f} {flag}"
        )
    lines.append("")
    lines.append("Li-7 puzzle resolution:")
    lines.append(
        f"  SBBN prediction   = {result.li7_sbbn_prediction_x1e10:.2f} × 10⁻¹⁰  "
        f"({result.li7_sbbn_n_sigma:.2f}σ from observation)"
    )
    lines.append(
        f"  Substrate         = {result.li7_substrate_prediction_x1e10:.2f} × 10⁻¹⁰  "
        f"({result.li7_substrate_n_sigma:.2f}σ from observation)"
    )
    lines.append(
        f"  Substrate beats SBBN: {result.li7_substrate_beats_sbbn}"
    )
    lines.append("")
    lines.append(f"Verdict: {result.verdict}")
    lines.append("=" * 78)
    return "\n".join(lines)


def default_test() -> BBNTestResult:
    """Convenience: run the test with all defaults."""
    return run_bbn_test()


if __name__ == "__main__":
    print(report(default_test()))
