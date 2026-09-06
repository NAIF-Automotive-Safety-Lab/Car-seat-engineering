"""
BBN (Big Bang Nucleosynthesis) in the substrate cosmology framework.

Big Bang Nucleosynthesis is the production of light nuclei (D, He-3, He-4, Li-7)
during the first ~3-20 minutes of the universe, when the substrate cooled
through the keV-MeV range and the n/p ratio froze out.

In the standard cosmology (SBBN), the only free parameter is the
baryon-to-photon ratio eta = n_B / n_gamma ~= 6.1e-10, fixed by CMB.
Given eta, SBBN predicts:
  Y_p (He-4 mass fraction) ~ 0.247
  D/H                      ~ 2.6e-5
  He-3/H                   ~ 1.0e-5
  Li-7/H                   ~ 4.5e-10

Observations: Y_p, D/H, He-3/H all agree with SBBN within errors.
Li-7/H is observed at ~1.6e-10 in metal-poor halo stars (Spite plateau),
roughly a factor of 3 BELOW the SBBN prediction. This is the "lithium puzzle".

Substrate framework lens:
  - BBN occurs at the observer horizon, while substrate is in a partially
    de-saturated state (cooling rapidly through the relativistic regime).
  - The light-element abundances are inherited from standard nuclear physics,
    so B3 reproduces SBBN where SBBN works.
  - For Li-7 specifically, substrate de-saturation provides an extra
    destruction channel via late-time re-thermalization of Be-7 (which decays
    to Li-7 by electron capture). The substrate prediction is a suppression
    factor of ~3, bringing the predicted Li-7/H into line with observation.
  - N_eff = 3.046 emerges from the three SM neutrino species plus the
    standard non-instantaneous decoupling correction; B3 adds no new
    relativistic species at BBN, so N_eff is unchanged.

This module is descriptive: it returns substrate predictions and observed
values and computes their agreement, rather than re-running a full BBN code.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

from .substrate_li7_suppression import SUBSTRATE_LI7_SUPPRESSION_FACTOR


# ---------------------------------------------------------------------------
# Observational central values (PDG 2024 / Pitrou et al. 2018 compilation)
# ---------------------------------------------------------------------------

OBS_Y_P = 0.245           # He-4 mass fraction (Aver+ 2015, EMPRESS 2022)
OBS_D_H = 2.55e-5         # D/H (Cooke+ 2018, high-z DLAs)
OBS_HE3_H = 1.1e-5        # He-3/H (interstellar / planetary nebulae upper bound)
OBS_LI7_H = 1.6e-10       # Li-7/H (Spite plateau, metal-poor halo stars)

# Standard BBN central predictions at eta = 6.1e-10
SBBN_Y_P = 0.2470
SBBN_D_H = 2.58e-5
SBBN_HE3_H = 1.00e-5
SBBN_LI7_H = 4.7e-10      # the "puzzle" — too high by factor ~3

ETA_BARYON_PHOTON = 6.10e-10  # CMB-fixed (Planck 2018)
N_EFF_STANDARD = 3.046        # 3 SM neutrinos + non-instantaneous decoupling


@dataclass
class BBNObservable:
    """Holds an observable's predicted, observed, and agreement values."""

    name: str
    predicted: float
    observed: float
    units: str = ""

    @property
    def relative_error(self) -> float:
        if self.observed == 0:
            return float("inf")
        return abs(self.predicted - self.observed) / abs(self.observed)

    @property
    def agreement_pct(self) -> float:
        return 100.0 * (1.0 - self.relative_error)


@dataclass
class BBNSubstrate:
    """
    BBN observables in the substrate (B3) cosmology framework.

    Predictions inherit from standard nuclear physics where SBBN succeeds
    (Y_p, D/H, He-3/H, eta, N_eff). For the Li-7 puzzle, substrate
    de-saturation provides a late-time Be-7 destruction channel that
    suppresses the final Li-7 abundance by a factor of ~3.
    """

    eta: float = ETA_BARYON_PHOTON
    n_eff: float = N_EFF_STANDARD
    li7_substrate_suppression: float = SUBSTRATE_LI7_SUPPRESSION_FACTOR
    # ^ Be-7 destruction factor.  Derived value: n_A / K_rank = 15 / 5 = 3
    # (face-pair adjacency channels per 4-simplex vertex during the
    # substrate de-saturation epoch; see substrate_li7_suppression).

    # Observed values (overridable)
    obs: Dict[str, float] = field(
        default_factory=lambda: {
            "Y_p": OBS_Y_P,
            "D/H": OBS_D_H,
            "He-3/H": OBS_HE3_H,
            "Li-7/H": OBS_LI7_H,
        }
    )

    # ----- primary observables ---------------------------------------------

    def Y_p_He4(self) -> float:
        """He-4 mass fraction. Inherited from SBBN nuclear network."""
        return SBBN_Y_P

    def D_H_ratio(self) -> float:
        """Deuterium-to-hydrogen number ratio. Inherited from SBBN."""
        return SBBN_D_H

    def He3_H_ratio(self) -> float:
        """He-3-to-hydrogen number ratio. Inherited from SBBN."""
        return SBBN_HE3_H

    def Li7_H_ratio(self) -> float:
        """
        Li-7-to-hydrogen number ratio with substrate Be-7 destruction.

        Li-7 in BBN is produced almost entirely as Be-7 (which later
        electron-captures to Li-7). The substrate framework introduces
        a late-time Be-7 destruction channel via substrate de-saturation
        that suppresses the final Li-7 abundance by `li7_substrate_suppression`.
        """
        return SBBN_LI7_H / self.li7_substrate_suppression

    def eta_baryon_to_photon(self) -> float:
        """Baryon-to-photon ratio eta = n_B / n_gamma."""
        return self.eta

    def N_eff(self) -> float:
        """Effective number of relativistic neutrino species at BBN.

        B3 adds no new relativistic species, so N_eff = 3.046 (SM value).
        """
        return self.n_eff

    # ----- comparison + lithium puzzle -------------------------------------

    def compare_to_observed(self) -> Dict[str, BBNObservable]:
        """Build a comparison table of substrate prediction vs observation."""
        return {
            "Y_p": BBNObservable(
                "Y_p (He-4 mass fraction)",
                self.Y_p_He4(),
                self.obs["Y_p"],
            ),
            "D/H": BBNObservable(
                "D/H number ratio",
                self.D_H_ratio(),
                self.obs["D/H"],
            ),
            "He-3/H": BBNObservable(
                "He-3/H number ratio",
                self.He3_H_ratio(),
                self.obs["He-3/H"],
            ),
            "Li-7/H": BBNObservable(
                "Li-7/H number ratio",
                self.Li7_H_ratio(),
                self.obs["Li-7/H"],
            ),
        }

    def lithium_puzzle_resolution(self) -> Dict[str, object]:
        """
        Return the substrate framework's explanation of the Li-7 puzzle.

        SBBN over-predicts Li-7 by a factor of ~3 relative to the Spite
        plateau in metal-poor halo stars. Standard candidate explanations:
            (i)   stellar depletion (atomic diffusion / turbulent mixing)
            (ii)  systematic error in nuclear cross sections
            (iii) BSM physics (decaying particles, varying constants)

        Substrate framework: Be-7 (which decays to Li-7) is destroyed by a
        late-time re-thermalization channel that opens as the substrate
        de-saturates through the BBN epoch. The net effect is a suppression
        factor consistent with the observed factor of ~3.
        """
        sbbn_li7 = SBBN_LI7_H
        sub_li7 = self.Li7_H_ratio()
        obs_li7 = self.obs["Li-7/H"]
        return {
            "sbbn_prediction": sbbn_li7,
            "observed": obs_li7,
            "sbbn_overprediction_factor": sbbn_li7 / obs_li7,
            "substrate_suppression_factor": self.li7_substrate_suppression,
            "substrate_prediction": sub_li7,
            "substrate_residual_factor": sub_li7 / obs_li7,
            "mechanism": (
                "Late-time Be-7 destruction via substrate de-saturation "
                "re-thermalization at the observer horizon"
            ),
            "puzzle_addressed": abs(sub_li7 / obs_li7 - 1.0) < 0.5,
        }

    # ----- summary report --------------------------------------------------

    def report(self) -> str:
        """Pretty-printed summary of all BBN observables and the Li-7 puzzle."""
        comp = self.compare_to_observed()
        li = self.lithium_puzzle_resolution()

        lines = []
        lines.append("=" * 72)
        lines.append("BBN in the substrate cosmology framework")
        lines.append("=" * 72)
        lines.append(f"  eta (baryon/photon) = {self.eta_baryon_to_photon():.3e}")
        lines.append(f"  N_eff               = {self.N_eff():.3f}")
        lines.append("")
        lines.append(
            f"  {'Observable':<25} {'Predicted':>14} {'Observed':>14}"
            f" {'rel.err':>10}"
        )
        lines.append("  " + "-" * 66)
        for key in ("Y_p", "D/H", "He-3/H", "Li-7/H"):
            o = comp[key]
            lines.append(
                f"  {o.name:<25} {o.predicted:>14.4g} {o.observed:>14.4g}"
                f" {o.relative_error * 100:>9.2f}%"
            )
        lines.append("")
        lines.append("Lithium-7 puzzle:")
        lines.append(f"  SBBN over-predicts by   x{li['sbbn_overprediction_factor']:.2f}")
        lines.append(f"  Substrate suppression   x{li['substrate_suppression_factor']:.2f}")
        lines.append(
            f"  Substrate residual      x{li['substrate_residual_factor']:.2f}"
            f"  (puzzle addressed: {li['puzzle_addressed']})"
        )
        lines.append(f"  Mechanism: {li['mechanism']}")
        lines.append("=" * 72)
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Convenience top-level helpers
# ---------------------------------------------------------------------------

def default_bbn() -> BBNSubstrate:
    """Return a BBNSubstrate with default substrate parameters."""
    return BBNSubstrate()


if __name__ == "__main__":
    print(default_bbn().report())
