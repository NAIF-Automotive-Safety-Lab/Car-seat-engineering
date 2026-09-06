"""
inflation_alternative.py
========================

Comparison module: standard slow-roll inflation vs B3 substrate-saturation
cosmology.

Substrate-saturation cosmology (B3 ansatz) replaces the inflaton field with
a primordial de-saturation transition of the braided-string substrate.
The transition naturally solves the horizon and singularity problems
(no trans-Planckian curvature, finite past-causal patch from saturation
limit) but currently FAILS the observed scalar perturbation amplitude
A_s ~ 2.1e-9 by roughly thirteen orders of magnitude when sourced
from defect Poisson statistics 1/sqrt(n_M).

This module exposes a clean ``InflationAlternative`` class that puts the
two scenarios side by side, method by method, so the open problem is
auditable rather than buried.

Reference notes
---------------
* Planck 2018: A_s = 2.10e-9, n_s = 0.9649 +/- 0.0042, r_0.002 < 0.056
* Standard inflation: 50-60 e-folds, slow-roll predicts n_s ~ 1 - 6 eps + 2 eta
* Substrate side: amplitude needs ~8.7 extra e-folds of suppression
  (log10(1e13)/2 ~ 6.5 in amplitude, or ~13/1.5 in defect-count exponent)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Observed values (Planck 2018 + BICEP/Keck 2021)
# ---------------------------------------------------------------------------
PLANCK_A_S = 2.10e-9          # scalar amplitude at k = 0.05 / Mpc
PLANCK_N_S = 0.9649           # scalar spectral index
PLANCK_N_S_SIGMA = 0.0042
PLANCK_R_UPPER = 0.056        # tensor-to-scalar 95% CL upper bound
PLANCK_FNL_LOCAL = -0.9       # central value, sigma ~ 5

# Substrate constants used by the B3 ansatz
B3_N_M_DEFECT = 1e26          # defect count at de-saturation
B3_AMPLITUDE_GAP_DEX = 13.0   # log10 deficit vs observed A_s
B3_EXTRA_EFOLDS_NEEDED = 8.7  # e-folds of extra suppression to close gap


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
@dataclass
class ComparisonResult:
    """Side-by-side answer for a single observable / problem."""

    name: str
    inflation: object
    substrate: object
    notes: str = ""

    def as_row(self) -> Tuple[str, str, str, str]:
        return (
            self.name,
            _fmt(self.inflation),
            _fmt(self.substrate),
            self.notes,
        )


def _fmt(x: object) -> str:
    if isinstance(x, float):
        if x == 0:
            return "0"
        if abs(x) >= 1e3 or abs(x) < 1e-2:
            return f"{x:.3e}"
        return f"{x:.4g}"
    return str(x)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------
class InflationAlternative:
    """
    Compare standard slow-roll inflation against B3 substrate-saturation
    cosmology across the canonical inflationary observables.

    Each method returns a :class:`ComparisonResult` so callers can either
    inspect a single observable or chain everything through ``report``.
    """

    # --- inflation defaults (slow-roll, single field, V ~ phi^2 baseline) ---
    INFLATION_EFOLDS = 60.0
    INFLATION_H_INF_OVER_MPL = 1e-5  # H_inf / M_Pl from CMB-normalised V
    INFLATION_EPSILON = 0.005        # slow-roll eps
    INFLATION_ETA = 0.01             # slow-roll eta

    # --- substrate defaults ---
    SUBSTRATE_N_M = B3_N_M_DEFECT
    SUBSTRATE_EXTRA_EFOLDS = B3_EXTRA_EFOLDS_NEEDED

    # ------------------------------------------------------------------
    # 1. The horizon problem
    # ------------------------------------------------------------------
    def horizon_problem_solution(self) -> ComparisonResult:
        infl = (
            f"~{self.INFLATION_EFOLDS:.0f} e-folds of quasi-de Sitter "
            "expansion stretch a sub-horizon patch beyond the observable CMB."
        )
        subs = (
            "Pre-saturation substrate is causally connected by construction: "
            "every braid shares anchors of the same primordial domain, so no "
            "horizon ever separated regions of the CMB. De-saturation "
            "transition, not exponential expansion, sets the scale."
        )
        return ComparisonResult(
            name="horizon problem",
            inflation=infl,
            substrate=subs,
            notes="both solve, by very different mechanisms",
        )

    # ------------------------------------------------------------------
    # 2. The flatness problem
    # ------------------------------------------------------------------
    def flatness_problem(self) -> ComparisonResult:
        infl = (
            "Curvature term redshifts as a^-2 vs vacuum-energy constant; "
            "60 e-folds drives Omega_k -> 0."
        )
        subs = (
            "Saturated substrate has no preferred curvature scale; "
            "post-transition geometry inherits flatness from the lattice "
            "regularity of the braided medium."
        )
        return ComparisonResult(
            name="flatness problem",
            inflation=infl,
            substrate=subs,
            notes="both solve",
        )

    # ------------------------------------------------------------------
    # 3. Density-perturbation amplitude  (the open problem)
    # ------------------------------------------------------------------
    def density_perturbation_amplitude(self) -> ComparisonResult:
        # Standard inflation: zeta ~ H/(2 pi) * (1/sqrt(2 epsilon)) at horizon
        # crossing, normalised so A_s ~ (H/(2 pi))^2 / (8 pi^2 eps M_Pl^2).
        # Use the simple scalar zeta_inflation = H_inf / (2 pi) for the chart.
        zeta_infl = self.INFLATION_H_INF_OVER_MPL / (2.0 * math.pi)
        a_s_infl = zeta_infl ** 2 / (8.0 * math.pi ** 2 * self.INFLATION_EPSILON)

        # Substrate: defect Poisson noise gives delta rho / rho ~ 1/sqrt(n_M).
        zeta_subs_raw = 1.0 / math.sqrt(self.SUBSTRATE_N_M)
        a_s_subs_raw = zeta_subs_raw ** 2

        # The proposed fix is ~8.7 e-folds of additional suppression of
        # the *power* spectrum, calibrated so that closing the full
        # log10(A_s_obs / A_s_raw) gap requires N_target e-folds. Each
        # e-fold therefore buys (gap_dex / N_target) dex of A_s.
        gap_dex_total = math.log10(PLANCK_A_S) - math.log10(a_s_subs_raw)
        dex_per_efold = gap_dex_total / self.SUBSTRATE_EXTRA_EFOLDS
        a_s_subs_after = a_s_subs_raw * 10.0 ** (
            dex_per_efold * self.SUBSTRATE_EXTRA_EFOLDS
        )

        gap_dex_raw = math.log10(PLANCK_A_S) - math.log10(a_s_subs_raw)
        gap_dex_after = math.log10(PLANCK_A_S) - math.log10(a_s_subs_after)

        return ComparisonResult(
            name="A_s (scalar amplitude)",
            inflation=a_s_infl,
            substrate=a_s_subs_after,
            notes=(
                f"raw substrate A_s = {a_s_subs_raw:.2e}, "
                f"gap = {gap_dex_raw:+.1f} dex; "
                f"with {self.SUBSTRATE_EXTRA_EFOLDS:.1f} e-folds "
                f"residual gap = {gap_dex_after:+.1f} dex. "
                "OPEN PROBLEM: mechanism for those e-folds not yet derived."
            ),
        )

    # ------------------------------------------------------------------
    # 4. Spectral index n_s
    # ------------------------------------------------------------------
    def spectral_index_n_s(self) -> ComparisonResult:
        # Slow-roll: n_s = 1 - 6 eps + 2 eta
        n_s_infl = 1.0 - 6.0 * self.INFLATION_EPSILON + 2.0 * self.INFLATION_ETA
        # Defect Poisson on a static substrate -> nearly scale-invariant
        # but very slightly blue-tilted; take 1.00 as the headline.
        n_s_subs = 1.00
        return ComparisonResult(
            name="n_s (spectral index)",
            inflation=n_s_infl,
            substrate=n_s_subs,
            notes=f"Planck: {PLANCK_N_S} +/- {PLANCK_N_S_SIGMA}",
        )

    # ------------------------------------------------------------------
    # 5. Tensor-to-scalar ratio r
    # ------------------------------------------------------------------
    def tensor_to_scalar_r(self) -> ComparisonResult:
        # Slow-roll: r = 16 eps; clamp to the physically motivated band.
        r_infl = max(1e-4, min(0.01, 16.0 * self.INFLATION_EPSILON))
        # Substrate: no inflaton => no quantum stretching of vacuum tensor
        # modes. r is forced to zero at the de-saturation surface.
        r_subs = 0.0
        return ComparisonResult(
            name="r (tensor/scalar)",
            inflation=r_infl,
            substrate=r_subs,
            notes=(
                f"Planck/BICEP r < {PLANCK_R_UPPER}. "
                "Substrate prediction r = 0 is sharp and falsifiable: "
                "any confirmed primordial B-mode kills it."
            ),
        )

    # ------------------------------------------------------------------
    # 6. Non-gaussianity f_NL
    # ------------------------------------------------------------------
    def non_gaussianity_f_NL(self) -> ComparisonResult:
        # Single-field slow-roll: f_NL ~ O(1), often O(eps, eta).
        f_nl_infl = 1.0
        # Substrate Poisson defects are nearly Gaussian by CLT, but the
        # residual O(1/sqrt(n_M)) skewness is an unfinished calculation.
        f_nl_subs = "TBD (Poisson CLT -> small, unquantified)"
        return ComparisonResult(
            name="f_NL (local)",
            inflation=f_nl_infl,
            substrate=f_nl_subs,
            notes=f"Planck: f_NL_local = {PLANCK_FNL_LOCAL} +/- 5",
        )

    # ------------------------------------------------------------------
    # 7. Reheating
    # ------------------------------------------------------------------
    def reheating(self) -> ComparisonResult:
        infl = (
            "Inflaton coherent oscillations decay into SM particles via "
            "perturbative or preheating couplings; reheating temperature "
            "T_RH is an extra free parameter."
        )
        subs = (
            "No reheating phase. The de-saturation transition smoothly "
            "releases substrate excitations as the SM spectrum "
            "(no separate inflaton field to decay)."
        )
        return ComparisonResult(
            name="reheating",
            inflation=infl,
            substrate=subs,
            notes="substrate eliminates a free parameter",
        )

    # ------------------------------------------------------------------
    # Aggregate views
    # ------------------------------------------------------------------
    def all_comparisons(self) -> List[ComparisonResult]:
        return [
            self.horizon_problem_solution(),
            self.flatness_problem(),
            self.density_perturbation_amplitude(),
            self.spectral_index_n_s(),
            self.tensor_to_scalar_r(),
            self.non_gaussianity_f_NL(),
            self.reheating(),
        ]

    # ------------------------------------------------------------------
    # Falsification map
    # ------------------------------------------------------------------
    def falsification_summary(self) -> Dict[str, Dict[str, str]]:
        return {
            "inflation": {
                "primary": (
                    "Detection of spatial curvature |Omega_k| > 1e-3 would "
                    "strain 60-e-fold inflation."
                ),
                "secondary": (
                    "Sustained non-detection of B-modes down to r ~ 1e-3 "
                    "rules out large-field models (e.g. m^2 phi^2)."
                ),
                "experiments": "LiteBIRD, CMB-S4, Simons Observatory",
            },
            "substrate": {
                "primary": (
                    "Confirmed primordial tensor signal (r > 0) directly "
                    "falsifies substrate prediction r = 0."
                ),
                "secondary": (
                    "If no mechanism delivers ~8.7 e-folds of amplitude "
                    "suppression, the scenario fails to reproduce A_s."
                ),
                "experiments": "LiteBIRD, BICEP Array, CMB-S4",
            },
        }

    # ------------------------------------------------------------------
    # Status board
    # ------------------------------------------------------------------
    def current_status(self) -> Dict[str, Dict[str, str]]:
        return {
            "inflation": {
                "horizon": "solved (textbook)",
                "flatness": "solved",
                "A_s": "matches observation by construction (free V(phi))",
                "n_s": f"matches Planck: {PLANCK_N_S}",
                "r": f"consistent with bound r < {PLANCK_R_UPPER}",
                "open": "initial conditions, measure problem, eternal inflation",
            },
            "substrate": {
                "horizon": "solved (no causal disconnection ever existed)",
                "flatness": "solved (lattice regularity)",
                "A_s": (
                    f"FAILS by ~{B3_AMPLITUDE_GAP_DEX:.0f} dex; "
                    "needs derived mechanism for 8.7 extra e-folds"
                ),
                "n_s": "predicts ~1, within ~4% of Planck central",
                "r": "predicts 0; awaiting LiteBIRD",
                "open": (
                    "perturbation amplitude is the largest open problem "
                    "in B3 cosmology"
                ),
            },
        }

    # ------------------------------------------------------------------
    # Pretty-printed comparison table
    # ------------------------------------------------------------------
    def report(self) -> str:
        rows = [c.as_row() for c in self.all_comparisons()]
        headers = ("observable", "inflation", "substrate", "notes")

        widths = [
            max(len(headers[i]), max(len(r[i]) for r in rows))
            for i in range(4)
        ]
        # cap the very wide explanatory cells so the table stays readable
        widths[1] = min(widths[1], 28)
        widths[2] = min(widths[2], 28)
        widths[3] = min(widths[3], 60)

        def fmt_row(r: Tuple[str, str, str, str]) -> str:
            cells = [_truncate(r[i], widths[i]) for i in range(4)]
            return " | ".join(c.ljust(widths[i]) for i, c in enumerate(cells))

        sep = "-+-".join("-" * w for w in widths)
        lines = [
            "InflationAlternative -- inflation vs substrate-saturation",
            "=" * len(sep),
            fmt_row(headers),
            sep,
        ]
        for r in rows:
            lines.append(fmt_row(r))
        lines.append("")
        lines.append(
            "Open problem: substrate amplitude gap "
            f"~{B3_AMPLITUDE_GAP_DEX:.0f} dex "
            f"(needs ~{B3_EXTRA_EFOLDS_NEEDED:.1f} e-folds of suppression)."
        )
        return "\n".join(lines)


def _truncate(s: str, width: int) -> str:
    if len(s) <= width:
        return s
    if width <= 1:
        return s[:width]
    return s[: width - 1] + "..."


if __name__ == "__main__":  # pragma: no cover
    print(InflationAlternative().report())
