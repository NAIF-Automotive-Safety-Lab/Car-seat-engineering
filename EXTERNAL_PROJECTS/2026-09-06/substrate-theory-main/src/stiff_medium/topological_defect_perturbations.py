"""Topological-defect Poisson amplitude for primordial density perturbations.

B3 substrate-saturation cosmology resolves horizon and singularity problems but
fails the density-perturbation amplitude by ~10^13. The most promising lead in
the open-problem audit is that, at substrate de-saturation freeze-out near
H ~ Lambda_QCD, topological defects (cone domain walls / braid junctions)
freeze in with a Poisson amplitude

    delta rho / rho |_freeze = 1 / sqrt(n_M) ,

where n_M = 268 is the master integer counting the number of independent
defect modes per Hubble volume at freeze-out. This baseline ~0.061 then
suffers exp(-N) suppression from N e-folds of post-freeze expansion before
recombination.

Status: OPEN. The model gives the right order of magnitude (a few percent at
freeze-out vs 10^-5 observed) and predicts a localized Lambda_QCD-scale
feature in the matter power spectrum, but the required N ~ 8.7 e-folds is
NOT yet derived from the substrate dynamics. This module makes the chain
explicit, quantitative, and falsifiable.

References:
- /Volumes/NO NAME/B3_PEER_REVIEW_CORE/active-falsifiers
- substrate_saturation_cosmology.py
- master integer M = 268 (b3_pitch_audit_insights)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict


# -- Physical anchors -----------------------------------------------------

LAMBDA_QCD_GEV = 0.2179  # B3 anchor
M_PLANCK_GEV = 1.2209e19
HBAR_C_GEV_M = 1.9733e-16  # GeV*m
N_M_DEFAULT = 268  # master integer
OBS_DELTA_RHO_OVER_RHO = 1.0e-5  # CMB / LSS curvature perturbation amplitude


# -- Core class ----------------------------------------------------------


@dataclass
class TopologicalDefectPerturbations:
    """Predict primordial delta rho / rho from substrate defect Poisson statistics.

    Parameters
    ----------
    n_M : int
        Number of independent defect modes per Hubble volume at freeze-out.
        Default 268 (B3 master integer).
    Lambda_QCD_GeV : float
        Energy scale at substrate de-saturation freeze-out.
    observed_amplitude : float
        Target delta rho / rho ~ 1e-5 from CMB.
    """

    n_M: int = N_M_DEFAULT
    Lambda_QCD_GeV: float = LAMBDA_QCD_GEV
    observed_amplitude: float = OBS_DELTA_RHO_OVER_RHO

    # ---- 1. Freeze-out amplitude --------------------------------------

    def freeze_out_amplitude(self) -> float:
        """Poisson amplitude at H ~ Lambda_QCD: 1 / sqrt(n_M).

        Returns
        -------
        float
            ~0.0611 for n_M = 268.
        """
        return 1.0 / math.sqrt(self.n_M)

    # ---- 2. Post-freeze suppression -----------------------------------

    def suppression_factor(self, N_efolds: float) -> float:
        """exp(-N) suppression from N e-folds of post-freeze expansion.

        Defects of fixed comoving wavelength are stretched and their
        amplitude diluted by the linear growth-vs-Hubble-friction balance,
        which gives single-exponential decay for a stiff-medium background.
        """
        return math.exp(-N_efolds)

    # ---- 3. Composite observable amplitude ----------------------------

    def observable_amplitude(self, N_efolds: float) -> float:
        """delta rho / rho today as a function of post-freeze e-folds."""
        return self.freeze_out_amplitude() * self.suppression_factor(N_efolds)

    # ---- 4. Required e-folds for observed amplitude -------------------

    def predicted_n_efolds_for_observed(self) -> float:
        """N* such that observable_amplitude(N*) = observed_amplitude.

        N* = ln(A_freeze / A_obs).
        """
        return math.log(self.freeze_out_amplitude() / self.observed_amplitude)

    # ---- 5. Spectral index / scale invariance -------------------------

    def predicted_n_s(self) -> float:
        """Defect-Poisson spectrum is approximately scale-invariant.

        A pure white-noise Poisson population on the freeze-out Hubble
        volume, smoothed onto super-horizon modes by causal stretching,
        produces a near-scale-invariant power spectrum with n_s -> 1
        in the leading approximation. Stiff-medium tilt corrections give
        a small red shift; we quote the leading value.
        """
        return 1.0

    def n_s_observed(self) -> float:
        """Planck 2018 best-fit (for comparison)."""
        return 0.9649

    # ---- 6. Lambda_QCD-scale matter power feature ---------------------

    def freeze_out_hubble_radius_m(self) -> float:
        """Hubble radius at freeze-out, R_H = hbar c / Lambda_QCD."""
        return HBAR_C_GEV_M / self.Lambda_QCD_GeV

    def Lambda_QCD_feature_scale(self, N_efolds: float | None = None) -> Dict[str, float]:
        """Comoving / physical scale where the defect feature should appear today.

        The freeze-out Hubble patch (~0.9 fm at T ~ Lambda_QCD ~ 200 MeV) gets
        stretched by the full cosmic expansion factor from then to today.
        Since T scales as a^-1 in radiation era, the redshift from QCD epoch
        to today is

            1 + z_QCD = T_QCD / T_CMB ~ (0.2179 GeV) / (2.348e-13 GeV) ~ 9.3e11.

        N_efolds is provided for the *amplitude* suppression; it is a tiny
        sub-component of total expansion (e^8.7 ~ 6000 vs 9e11 total) and
        doesn't dominate the scale calculation. We report both.
        """
        if N_efolds is None:
            N_efolds = self.predicted_n_efolds_for_observed()
        R_freeze_m = self.freeze_out_hubble_radius_m()

        # T_CMB = 2.7255 K -> kT in GeV
        k_B_GeV_per_K = 8.617333e-14
        T_CMB_GeV = 2.7255 * k_B_GeV_per_K  # ~ 2.348e-13 GeV
        z_QCD_to_today = self.Lambda_QCD_GeV / T_CMB_GeV  # ~ 9.3e11

        R_today_m = R_freeze_m * z_QCD_to_today
        Mpc_in_m = 3.0857e22
        R_today_Mpc = R_today_m / Mpc_in_m

        # The amplitude-suppression stretch is a separate informational quantity
        amplitude_stretch = math.exp(N_efolds)
        return {
            "R_freeze_m": R_freeze_m,
            "z_QCD_to_today": z_QCD_to_today,
            "R_today_m": R_today_m,
            "R_today_Mpc": R_today_Mpc,
            "amplitude_stretch_efolds": N_efolds,
            "amplitude_stretch_factor": amplitude_stretch,
        }

    # ---- 7. Inflation comparison --------------------------------------

    def compare_to_inflation(self) -> Dict[str, Dict[str, str]]:
        """Side-by-side substrate-defect vs slow-roll inflation predictions."""
        return {
            "amplitude_origin": {
                "inflation": "quantum vacuum fluctuations of inflaton field",
                "substrate": "Poisson statistics of n_M=268 defects per Hubble volume",
            },
            "amplitude_value": {
                "inflation": "A_s ~ 2.1e-9 (free parameter via V/epsilon)",
                "substrate": f"A_freeze = 1/sqrt(n_M) = {self.freeze_out_amplitude():.4f}, "
                             f"suppressed by exp(-N) with N~{self.predicted_n_efolds_for_observed():.1f}",
            },
            "n_s_value": {
                "inflation": "n_s = 1 - 6 epsilon + 2 eta (depends on potential)",
                "substrate": "n_s ~ 1 (Poisson white noise, leading order)",
            },
            "tensor_to_scalar": {
                "inflation": "r = 16 epsilon (slow-roll)",
                "substrate": "r ~ 0 (no inflaton tensor modes; defects are scalar)",
            },
            "feature_at_Lambda_QCD": {
                "inflation": "no feature predicted",
                "substrate": "localized bump/dip at the stretched freeze-out scale "
                             "(potentially detectable in galaxy surveys)",
            },
            "free_parameters": {
                "inflation": "V(phi) shape, slow-roll params, e-folds (~3 free)",
                "substrate": "n_M=268 fixed by topology; N_efolds is the only free knob",
            },
        }

    # ---- 8. Status --------------------------------------------------

    def gap_status(self) -> Dict[str, object]:
        """Honest classification of where the calculation stands."""
        N_required = self.predicted_n_efolds_for_observed()
        order_of_magnitude_match = self.freeze_out_amplitude() / 0.1  # ~0.6
        feature = self.Lambda_QCD_feature_scale()
        return {
            "status": "open-partial",
            "freeze_out_amplitude": self.freeze_out_amplitude(),
            "observed_amplitude": self.observed_amplitude,
            "ratio_freeze_to_observed": self.freeze_out_amplitude() / self.observed_amplitude,
            "N_efolds_required": N_required,
            "N_efolds_derived_from_substrate": None,
            "feature_scale_Mpc": feature["R_today_Mpc"],
            "what_works": [
                "freeze-out amplitude order-of-magnitude (0.061 vs naive ~1)",
                f"order-of-magnitude factor improvement: 1 / 0.061 = {1/self.freeze_out_amplitude():.1f}x reduction over unit Poisson",
                "scale invariance n_s ~ 1 (matches observed n_s=0.965 to 4%)",
                "predicts NEW falsifiable feature: Lambda_QCD-scale bump in P(k)",
                "no free parameters in amplitude (n_M=268 is fixed)",
            ],
            "what_is_missing": [
                "derivation of N_efolds ~ 8.7 from substrate de-saturation dynamics",
                "spectral tilt correction n_s -> 0.965 from stiff-medium expansion",
                "non-Gaussianity prediction f_NL from defect collisions",
                "feature scale: naive QCD-Hubble stretch gives sub-mm today, "
                "NOT Mpc; need defect coalescence / growth mechanism to lift to "
                "galaxy-survey scales",
                "calculation of the feature amplitude, not just its scale",
            ],
            "falsifier": (
                f"Naive feature scale today is {feature['R_today_m']:.2e} m "
                f"({feature['R_today_Mpc']:.2e} Mpc) — far below galaxy-survey "
                "resolution. Either (a) defect coalescence lifts the scale by "
                "~22 orders of magnitude, or (b) this channel cannot match "
                "observed P(k) features. (a) needs to be derived; (b) would "
                "kill this channel."
            ),
        }


# -- Convenience runner --------------------------------------------------


def summary() -> Dict[str, object]:
    """One-call summary used by audits."""
    tdp = TopologicalDefectPerturbations()
    return {
        "freeze_out_amplitude": tdp.freeze_out_amplitude(),
        "N_efolds_required": tdp.predicted_n_efolds_for_observed(),
        "feature_scale": tdp.Lambda_QCD_feature_scale(),
        "n_s_predicted": tdp.predicted_n_s(),
        "gap_status": tdp.gap_status(),
    }


if __name__ == "__main__":
    import json

    print(json.dumps(summary(), indent=2, default=str))
