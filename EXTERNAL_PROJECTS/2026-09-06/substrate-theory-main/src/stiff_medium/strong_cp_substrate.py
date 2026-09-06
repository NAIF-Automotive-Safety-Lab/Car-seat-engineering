"""B3 strong CP problem — substrate-forced θ_QCD = 0.

The strong CP problem in the SM
-------------------------------
QCD admits a CP-violating term in its Lagrangian,

    L_θ = θ_QCD · (g²/32π²) · G_{μν}^a G̃^{aμν},

where G̃ is the dual gluon field strength. A non-zero θ_QCD generates a
neutron electric dipole moment (nEDM) of order

    d_n ~ θ_QCD · 10⁻¹⁶ e·cm                          (Crewther et al. 1979)

The current experimental bound from PSI (Abel et al., PRL 124, 081803, 2020)
is |d_n| < 1.8 × 10⁻²⁶ e·cm at 90% CL, which forces

    |θ_QCD| < 10⁻¹⁰.

In the Standard Model nothing protects θ_QCD: the chiral anomaly mixes it
with the quark-mass phase arg(det M_q), and there is no symmetry that sets
the sum to zero. The standard fix is the Peccei–Quinn U(1)_PQ axion, which
introduces a new dynamical pseudoscalar field that relaxes θ̄ → 0; this
costs one new particle (the axion) and an extra anomalous symmetry.

The B3 substrate resolution
---------------------------
In the B3 framework the QCD vacuum is realised on the braided-string
substrate, whose topology is a Möbius double cover (K_pair = 2 sheets).
The Möbius bundle has a global ℤ/2 sheet-swap automorphism σ that:

  * is **orientation-reversing** on the substrate,
  * commutes with all gauge action (the gauge bundle is pulled back from
    the base, not the cover),
  * and sends the dual field strength G̃ → −G̃ while leaving G·G invariant.

Hence the topological term G·G̃ is σ-odd. Because σ is a *gauge symmetry of
the substrate* (not a global one we can choose to break), the only
σ-invariant value of its coefficient is zero:

    σ · θ_QCD = −θ_QCD   ⇒   θ_QCD ≡ 0.

This is structurally identical to the way orientability of spacetime forces
the absence of a Pfaffian-like term in standard gauge theory; the
difference is that the B3 substrate is *non*-orientable in a controlled
way, and the resulting ℤ/2 symmetry is exactly what removes θ.

Predictions
-----------
  * θ_QCD = 0 exactly (no fine-tuning, no axion).
  * d_n = 0 from the strong sector at tree level; only the (tiny) CKM
    contribution d_n^CKM ~ 10⁻³² e·cm survives, ~6 orders of magnitude
    below the current PSI bound and any proposed near-future experiment.
  * No axion — the dark-matter axion window (μeV–meV) is *not* required
    for the strong CP solution in B3. Substrate provides its own DM
    candidate (see Σm_ν / saturation modules).

Falsifier
---------
Any positive nEDM measurement above the SM CKM floor (~10⁻³² e·cm) and
below the current bound would falsify the strict B3 prediction.

References
----------
  * Crewther, Di Vecchia, Veneziano, Witten, Phys. Lett. B 88 (1979) 123.
  * Abel et al., Phys. Rev. Lett. 124, 081803 (2020) — PSI nEDM result.
  * Peccei, Quinn, Phys. Rev. Lett. 38 (1977) 1440 — SM axion fix.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from .b3_constants import K_pair


# ---------------------------------------------------------------------------
# Reference numbers (experiment / SM)
# ---------------------------------------------------------------------------

NEDM_BOUND_PSI_2020: float = 1.8e-26     # e·cm, 90% CL upper limit
NEDM_THETA_COEFFICIENT: float = 1.0e-16  # d_n ≃ θ · 10⁻¹⁶ e·cm (Crewther)
NEDM_SM_CKM_FLOOR: float = 1.0e-32       # e·cm, irreducible CKM contribution
SM_THETA_BOUND: float = 1.0e-10          # |θ_QCD| < 10⁻¹⁰ from PSI bound


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TheoryComparison:
    """Side-by-side row for SM vs B3-substrate accounts of strong CP."""

    aspect: str
    standard_model: str
    b3_substrate: str


# ---------------------------------------------------------------------------
# StrongCP class
# ---------------------------------------------------------------------------

class StrongCP:
    """Strong-CP / θ_QCD predictions from the B3 Möbius substrate.

    The class is essentially a wrapper around a single fact (θ ≡ 0 by
    Möbius ℤ/2) plus the experimental and SM reference numbers needed to
    compare against PSI 2020 and to articulate why no axion is required.
    """

    sheet_count: int = K_pair  # = 2 (Möbius double cover)

    # ---- core predictions ------------------------------------------------

    def theta_substrate_prediction(self) -> float:
        """Return predicted θ_QCD on the B3 substrate.

        The Möbius double cover carries a global ℤ/2 sheet-swap symmetry σ
        that flips G·G̃ → −G·G̃. σ is a gauge symmetry of the substrate, so
        its coefficient must be σ-invariant; the only σ-invariant value is
        zero. Hence θ_QCD = 0 *forced*, not fine-tuned.
        """
        # Sanity: the mechanism only works if K_pair = 2 (true Möbius).
        if self.sheet_count != 2:
            raise RuntimeError(
                f"Möbius ℤ/2 mechanism requires K_pair=2, got {self.sheet_count}"
            )
        return 0.0

    def neutron_edm_prediction(self) -> Dict[str, float]:
        """Predicted neutron EDM from the strong sector.

        Strong contribution vanishes (θ = 0). The only surviving piece is
        the CKM contribution, which is the same as in the SM and is ~6
        orders of magnitude below the PSI bound.
        """
        theta = self.theta_substrate_prediction()
        d_n_strong = theta * NEDM_THETA_COEFFICIENT
        d_n_total = d_n_strong + NEDM_SM_CKM_FLOOR
        return {
            "theta_QCD": theta,
            "d_n_strong_e_cm": d_n_strong,
            "d_n_ckm_floor_e_cm": NEDM_SM_CKM_FLOOR,
            "d_n_total_e_cm": d_n_total,
        }

    # ---- experimental comparison ----------------------------------------

    def compare_to_experiments(self) -> Dict[str, object]:
        """Compare substrate prediction to PSI 2020 nEDM bound.

        Returns a dict with the experimental upper bound, the substrate
        prediction, the ratio (margin), and a pass/fail flag.
        """
        pred = self.neutron_edm_prediction()
        d_n_pred = pred["d_n_total_e_cm"]
        bound = NEDM_BOUND_PSI_2020
        margin = bound / d_n_pred if d_n_pred > 0 else float("inf")
        return {
            "experiment": "PSI nEDM (Abel et al., PRL 124, 081803, 2020)",
            "bound_e_cm": bound,
            "bound_CL": "90%",
            "substrate_prediction_e_cm": d_n_pred,
            "margin_below_bound": margin,
            "passes_bound": d_n_pred < bound,
            "implied_theta_bound_SM": SM_THETA_BOUND,
            "substrate_theta": pred["theta_QCD"],
        }

    # ---- alternative-explanation table ----------------------------------

    def axion_alternative(self) -> List[TheoryComparison]:
        """Side-by-side: Peccei–Quinn axion vs B3 substrate ℤ/2.

        Both schemes set θ̄ ≈ 0 dynamically/structurally. The PQ scheme
        adds a new pseudoscalar (the axion); the substrate scheme uses
        the existing Möbius ℤ/2 already present in the framework's
        topology, adding no new particles.
        """
        return [
            TheoryComparison(
                aspect="New particle introduced",
                standard_model="Yes — axion (m_a ~ μeV–meV window)",
                b3_substrate="No — uses existing Möbius ℤ/2 (K_pair=2)",
            ),
            TheoryComparison(
                aspect="New symmetry imposed",
                standard_model="Global U(1)_PQ, anomalous, broken at f_a",
                b3_substrate="Substrate gauge ℤ/2 sheet-swap (σ)",
            ),
            TheoryComparison(
                aspect="θ_QCD value",
                standard_model="Relaxed to 0 by axion VEV",
                b3_substrate="Identically 0 (σ-odd coefficient)",
            ),
            TheoryComparison(
                aspect="Free parameters",
                standard_model="f_a (axion decay constant), m_a",
                b3_substrate="Zero — K_pair=2 is fixed by topology",
            ),
            TheoryComparison(
                aspect="DM candidate",
                standard_model="Axion DM (μeV window, ADMX, etc.)",
                b3_substrate="Independent — substrate saturation sector",
            ),
            TheoryComparison(
                aspect="Quality problem",
                standard_model=(
                    "Yes — gravity is expected to break global U(1)_PQ; "
                    "needs extreme protection"
                ),
                b3_substrate="None — gauge symmetry, not global",
            ),
        ]

    # ---- pretty report --------------------------------------------------

    def report(self) -> str:
        """Human-readable summary suitable for stdout."""
        cmp = self.compare_to_experiments()
        edm = self.neutron_edm_prediction()
        axion = self.axion_alternative()

        lines: List[str] = []
        lines.append("=" * 66)
        lines.append("B3 substrate — strong CP problem")
        lines.append("=" * 66)
        lines.append("")
        lines.append("Mechanism:")
        lines.append("  Möbius double cover (K_pair = 2) carries a global ℤ/2")
        lines.append("  sheet-swap σ that flips G·G̃ → −G·G̃. σ is a substrate")
        lines.append("  gauge symmetry, so the θ-term coefficient must be σ-")
        lines.append("  invariant; the only invariant value is 0.")
        lines.append("")
        lines.append("Predictions:")
        lines.append(f"  θ_QCD                          = {edm['theta_QCD']:.3e}")
        lines.append(f"  d_n (strong, B3)               = {edm['d_n_strong_e_cm']:.3e} e·cm")
        lines.append(f"  d_n (CKM floor, irreducible)   = {edm['d_n_ckm_floor_e_cm']:.3e} e·cm")
        lines.append(f"  d_n (total prediction)         = {edm['d_n_total_e_cm']:.3e} e·cm")
        lines.append("")
        lines.append("Experimental comparison:")
        lines.append(f"  {cmp['experiment']}")
        lines.append(f"  Upper bound (90% CL)           = {cmp['bound_e_cm']:.3e} e·cm")
        lines.append(f"  Substrate prediction           = {cmp['substrate_prediction_e_cm']:.3e} e·cm")
        lines.append(f"  Margin below bound             = {cmp['margin_below_bound']:.3e}×")
        lines.append(f"  Passes bound                   = {cmp['passes_bound']}")
        lines.append("")
        lines.append("SM (Peccei–Quinn axion) vs B3 substrate:")
        width = max(len(row.aspect) for row in axion)
        for row in axion:
            lines.append(f"  {row.aspect.ljust(width)} | SM: {row.standard_model}")
            lines.append(f"  {' '.ljust(width)} | B3: {row.b3_substrate}")
        lines.append("")
        lines.append("Conclusion:")
        lines.append("  No axion required. θ_QCD = 0 is forced by the Möbius")
        lines.append("  ℤ/2 already present in the substrate topology, with")
        lines.append("  zero new particles or free parameters.")
        lines.append("=" * 66)
        return "\n".join(lines)


__all__ = [
    "StrongCP",
    "TheoryComparison",
    "NEDM_BOUND_PSI_2020",
    "NEDM_THETA_COEFFICIENT",
    "NEDM_SM_CKM_FLOOR",
    "SM_THETA_BOUND",
]


if __name__ == "__main__":
    print(StrongCP().report())
