"""
CMB B-mode polarization predictor for inflation tests.

B3 substrate-saturation cosmology REPLACES inflation. Inflation generically
predicts a tensor-to-scalar ratio r > 0 (primordial gravitational waves
imprinted as B-mode polarization in the CMB). Substrate-saturation has no
inflationary epoch and therefore predicts r = 0 at the primordial level.

Only the lensing B-mode signal (gravitational lensing of E-modes by
large-scale structure) survives. This is a generic, well-measured
contamination peaking near multipole l ~ 1000 with amplitude ~ 1e-4 uK^2.

Decisive falsifier: any robust detection of primordial r > 0 (above the
lensing floor) kills B3 substrate cosmology and confirms inflation.

Current state of the art:
    - BICEP/Keck 2021:  r < 0.036  (95% CL)
    - LiteBIRD (~2030): target sensitivity sigma(r) ~ 0.001
    - CMB-S4 (~2030+):  target sensitivity sigma(r) ~ 0.001

Both LiteBIRD and CMB-S4 will reach r ~ 1e-3, three orders of magnitude
below the current bound. Either they detect inflation, or substrate
cosmology survives a clean test.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable


# ---------------------------------------------------------------------------
# Experimental sensitivities (95% CL upper limits or 1-sigma forecasts)
# ---------------------------------------------------------------------------
BICEP_KECK_2021_LIMIT = 0.036       # r < 0.036 (95% CL)
LITEBIRD_TARGET_SIGMA_R = 0.001     # delta-r ~ 1e-3 (forecast)
CMB_S4_TARGET_SIGMA_R = 0.001       # delta-r ~ 1e-3 (forecast)

# Lensing B-mode peak parameters (D_l = l(l+1)C_l/(2pi) in uK^2)
LENSING_PEAK_L = 1000               # multipole of lensing B-mode bump
LENSING_PEAK_AMPLITUDE_UK2 = 1.0e-4 # ~ 1e-4 uK^2 at peak
LENSING_PEAK_WIDTH = 600            # rough Gaussian width in l


@dataclass(frozen=True)
class BModeComparison:
    """Comparison of substrate prediction vs an experiment."""
    experiment: str
    sensitivity_r: float
    substrate_r: float
    consistent: bool
    decisive: bool   # True if experiment's sensitivity is below typical inflation r
    note: str


class CMBBModePolarization:
    """
    CMB B-mode polarization predictor under B3 substrate-saturation cosmology.

    Substrate predicts r = 0 (no inflation -> no primordial gravitational
    wave background). Only lensing B-modes from gravitational lensing of
    E-modes survive. Any robust primordial r > 0 detection falsifies B3.
    """

    def __init__(
        self,
        lensing_peak_l: int = LENSING_PEAK_L,
        lensing_peak_amplitude_uK2: float = LENSING_PEAK_AMPLITUDE_UK2,
        lensing_peak_width: float = LENSING_PEAK_WIDTH,
    ) -> None:
        self.lensing_peak_l = lensing_peak_l
        self.lensing_peak_amplitude_uK2 = lensing_peak_amplitude_uK2
        self.lensing_peak_width = lensing_peak_width

    # -- substrate prediction -----------------------------------------------
    def tensor_to_scalar_ratio_substrate(self) -> float:
        """B3 substrate cosmology: no inflation, no primordial GW background."""
        return 0.0

    # -- lensing B-modes (always present, not from inflation) ---------------
    def lensing_b_mode(self, l: float) -> float:
        """
        Lensing B-mode D_l = l(l+1)C_l/(2pi) [uK^2] from grav. lensing of E-modes.

        Approximated as a Gaussian bump in l centred on lensing_peak_l.
        Real lensing spectrum is broader and asymmetric, but this captures
        the essential observational feature: a single broad bump near
        l ~ 1000 of amplitude ~ 1e-4 uK^2.
        """
        if l <= 0:
            return 0.0
        x = (l - self.lensing_peak_l) / self.lensing_peak_width
        return self.lensing_peak_amplitude_uK2 * math.exp(-0.5 * x * x)

    # -- inflation prediction (for comparison) ------------------------------
    def inflation_b_mode(self, l: float, r: float) -> float:
        """
        Approximate primordial inflation B-mode D_l in uK^2 for tensor ratio r.

        Two characteristic features:
          - reionisation bump near l ~ 5  (peak amplitude scales as r)
          - recombination bump near l ~ 80
        Amplitude normalisation: at r = 0.1 the recombination peak is
        D_l ~ 1e-2 uK^2 (canonical textbook value).
        """
        if r <= 0 or l <= 0:
            return 0.0
        # canonical normalization: r=0.1 -> recomb peak ~ 1e-2 uK^2
        norm = 1.0e-2 * (r / 0.1)
        # recombination bump
        x_rec = (l - 80.0) / 50.0
        recomb = norm * math.exp(-0.5 * x_rec * x_rec)
        # reionisation bump (smaller, at l ~ 5)
        x_rei = (l - 5.0) / 3.0
        reion = 0.3 * norm * math.exp(-0.5 * x_rei * x_rei)
        return recomb + reion

    # -- experimental comparisons -------------------------------------------
    def compare_to_bicep_keck(self) -> BModeComparison:
        """BICEP/Keck 2021: r < 0.036 (95% CL). Substrate r=0 is fully consistent."""
        r_sub = self.tensor_to_scalar_ratio_substrate()
        return BModeComparison(
            experiment="BICEP/Keck 2021",
            sensitivity_r=BICEP_KECK_2021_LIMIT,
            substrate_r=r_sub,
            consistent=(r_sub < BICEP_KECK_2021_LIMIT),
            decisive=False,  # not yet sensitive enough to be decisive
            note="Current upper limit r<0.036 (95% CL); substrate r=0 passes trivially.",
        )

    def compare_to_litebird(self) -> BModeComparison:
        """LiteBIRD (JAXA, ~2030): forecast sigma(r) ~ 1e-3. Decisive."""
        r_sub = self.tensor_to_scalar_ratio_substrate()
        return BModeComparison(
            experiment="LiteBIRD (~2030)",
            sensitivity_r=LITEBIRD_TARGET_SIGMA_R,
            substrate_r=r_sub,
            consistent=True,  # by construction r=0 < any positive sensitivity
            decisive=True,    # 1e-3 is below most slow-roll inflation models
            note="Forecast sigma(r)~1e-3; decisive test of inflation vs substrate.",
        )

    def compare_to_cmb_s4(self) -> BModeComparison:
        """CMB-S4 (ground, ~2030+): forecast sigma(r) ~ 1e-3. Decisive."""
        r_sub = self.tensor_to_scalar_ratio_substrate()
        return BModeComparison(
            experiment="CMB-S4 (~2030+)",
            sensitivity_r=CMB_S4_TARGET_SIGMA_R,
            substrate_r=r_sub,
            consistent=True,
            decisive=True,
            note="Forecast sigma(r)~1e-3; ground-based companion to LiteBIRD.",
        )

    # -- falsifier ----------------------------------------------------------
    def substrate_falsifier(self, r_observed: float, r_uncertainty: float) -> dict:
        """
        Decision rule: substrate cosmology is FALSIFIED if a measurement
        r_observed exceeds zero by more than 3 sigma (i.e. r_observed -
        3 * r_uncertainty > 0). In that case primordial gravitational
        waves from inflation are confirmed and substrate-saturation dies.
        """
        significance = r_observed / r_uncertainty if r_uncertainty > 0 else float('inf')
        falsified = (r_observed - 3.0 * r_uncertainty) > 0.0
        return {
            "r_observed": r_observed,
            "r_uncertainty": r_uncertainty,
            "significance_sigma": significance,
            "substrate_falsified": falsified,
            "verdict": (
                "FALSIFIED — primordial GW detected, inflation confirmed"
                if falsified
                else "CONSISTENT — no significant primordial signal"
            ),
        }

    # -- reporting ----------------------------------------------------------
    def report(self) -> str:
        lines = []
        lines.append("=" * 68)
        lines.append("CMB B-mode polarization — B3 substrate-saturation cosmology")
        lines.append("=" * 68)
        lines.append("")
        lines.append("Substrate prediction:")
        lines.append(f"  tensor-to-scalar ratio r = {self.tensor_to_scalar_ratio_substrate():.4f}")
        lines.append("  (no inflation -> no primordial GW background)")
        lines.append("")
        lines.append("Lensing B-mode (always present, not from inflation):")
        lines.append(f"  peak multipole       l ~ {self.lensing_peak_l}")
        lines.append(f"  peak D_l amplitude   ~ {self.lensing_peak_amplitude_uK2:.2e} uK^2")
        lines.append("")
        lines.append("Experimental comparisons:")
        for comp in (self.compare_to_bicep_keck(),
                     self.compare_to_litebird(),
                     self.compare_to_cmb_s4()):
            lines.append(f"  - {comp.experiment}")
            lines.append(f"      sensitivity: r ~ {comp.sensitivity_r:.4f}")
            lines.append(f"      substrate:   r = {comp.substrate_r:.4f}")
            lines.append(f"      consistent:  {comp.consistent}    decisive: {comp.decisive}")
            lines.append(f"      note: {comp.note}")
        lines.append("")
        lines.append("Falsifier:")
        lines.append("  Any > 3 sigma detection of primordial r > 0 (above the")
        lines.append("  lensing floor) kills B3 substrate cosmology.")
        lines.append("=" * 68)
        return "\n".join(lines)


def main() -> None:
    model = CMBBModePolarization()
    print(model.report())


if __name__ == "__main__":
    main()
