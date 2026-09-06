"""Renormalization-group running of coupling constants — §§18.9, 18.34, 18.46, 18.48.

The fine-structure constant α is not a fixed number; it *runs* with the energy scale Q
owing to quantum vacuum polarization.  At low energy α(0) = 1/137.035999177, while at the
electroweak scale α(M_Z) ≈ 1/127.952.

In the substrate-mechanical theory (§18.34) the effective low-energy theory is identical
to standard QED+EW+QCD in the appropriate limit, so the standard one-loop RG equations
apply unchanged. The still-open part is deriving the *bare* coupling from the Möbius
bundle (§§18.9, 18.10, 18.46): α = e²/(4π K ξ⁴). RG running can say what bare value is
required at a substrate scale, but it does not derive that value by itself.

Standard one-loop QED running (leptons + heavy quarks, perturbative):

    1/α(Q²) = 1/α(0) − (1/(3π)) Σ_f Q_f² N_c^f ln(Q²/m_f²)   [lep + c,b,t]
             − DELTA_INV_ALPHA_UDS_HAD                           [light quarks u,d,s]

where the first sum covers leptons (e, μ, τ) and heavy quarks (c, b, t) above their
respective mass thresholds.  Light quarks u, d, s are non-perturbative at the QCD
confinement scale: their contribution is taken from the measured hadronic vacuum
polarization (dispersive integral over e+e- → hadrons data, PDG 2024).

The constant ``DELTA_INV_ALPHA_UDS_HAD`` = −2.8226 is calibrated so that running from
α(m_e) to M_Z gives exactly 1/α(M_Z) = 127.952.  It is the Δ(1/α) contribution of
light-quark hadronic VP expressed in the additive-shift convention above.

This two-piece structure (perturbative leptons+heavy quarks, dispersive light quarks)
is the standard treatment in precision electroweak physics (see Burkhardt & Pietrzyk
2011; PDG 2024 §10 Electroweak Model review).

For the strong coupling α_s (QCD running) the one-loop beta function is:

    1/α_s(Q²) = 1/α_s(μ²) + b_0 ln(Q²/μ²)
    b_0 = (33 − 2 N_f) / (12π)

Asymptotic freedom: b_0 > 0 for N_f < 16.5, so α_s → 0 as Q → ∞.

References
----------
    spec §18.9  : α = e²/(Kξ⁴) dimensional analysis
    spec §18.34 : structural correspondence to standard QED → running is identical
    spec §18.46 : constants from substrate primitives; α remains an audited gap
    spec §18.48 : precision QED/QCD inheritance; §18.49 SU(3)/asymptotic freedom
    Coleman (1975): Phys. Rev. D 11, 2088 — bosonization / Thirring duality
    Peskin & Schroeder §12.3: one-loop running of α
    Burkhardt & Pietrzyk (2011): hadronic vacuum polarization Δα_had
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

from stiff_medium.mobius_bundle import coleman_bosonization_g
from stiff_medium.physical_constants import (
    ALPHA_THOMSON_CODATA_2022,
    INV_ALPHA_THOMSON_CODATA_2022,
)

# ---------------------------------------------------------------------------
# Physical constants (SI + natural-unit conversions)
# ---------------------------------------------------------------------------

PI: Final[float] = math.pi

#: Fine-structure constant at zero momentum transfer (Thomson limit).
#: NIST/CODATA 2022: alpha^-1 = 137.035999177(21).
ALPHA_0: Final[float] = ALPHA_THOMSON_CODATA_2022
INV_ALPHA_0: Final[float] = INV_ALPHA_THOMSON_CODATA_2022

#: Inverse fine-structure constant at the Z pole [PDG 2024].
INV_ALPHA_MZ_MEASURED: Final[float] = 127.952

#: Strong coupling at M_Z [PDG 2024].
ALPHA_S_MZ: Final[float] = 0.1179

#: Z-boson mass [GeV].
M_Z_GEV: Final[float] = 91.1876

#: Electron mass [GeV] — sets the lowest threshold in the QED running.
M_E_GEV: Final[float] = 0.511e-3

#: Muon mass [GeV].
M_MU_GEV: Final[float] = 0.10566

#: Tau mass [GeV].
M_TAU_GEV: Final[float] = 1.77686

#: Quark masses [GeV] (MS-bar at 2 GeV for light; pole mass for heavy).
M_UP_GEV: Final[float] = 2.2e-3
M_DOWN_GEV: Final[float] = 4.7e-3
M_STRANGE_GEV: Final[float] = 0.095
M_CHARM_GEV: Final[float] = 1.28
M_BOTTOM_GEV: Final[float] = 4.18
M_TOP_GEV: Final[float] = 172.69

#: Light-quark (u, d, s) hadronic vacuum-polarization contribution to Δ(1/α)
#: at scale M_Z, expressed as an additive shift in the inverse coupling.
#:
#: This constant encodes the non-perturbative QCD physics of the light quarks:
#: because u, d, s are confined, their vacuum-polarization contribution cannot
#: be computed from perturbative ln(Q²/m²) loops.  Instead it is extracted from
#: the dispersive integral over e+e- → hadrons cross-section data (PDG 2024,
#: Burkhardt & Pietrzyk 2011: Δα_had^{(5)}(M_Z) = 0.027572, of which the u,d,s
#: piece Δα_uds ≈ 0.02060 after subtracting perturbative c and b contributions).
#:
#: In the additive Δ(1/α) convention used here:
#:     Δ(1/α)_uds = −2.8226
#: This is calibrated so that running from m_e to M_Z gives exactly
#: 1/α(M_Z) = 127.952 (PDG 2024 central value).
DELTA_INV_ALPHA_UDS_HAD: Final[float] = -2.8226

#: Energy threshold above which the light-quark hadronic correction is applied [GeV].
#: Set to 1 GeV — above Λ_QCD where the correction becomes significant.
Q_HAD_THRESHOLD_GEV: Final[float] = 1.0

#: Λ_QCD (non-perturbative QCD scale) in GeV — MS-bar with N_f = 5 [PDG].
LAMBDA_QCD_GEV: Final[float] = 0.217

#: Number of active quark flavours at M_Z.
N_F_AT_MZ: Final[int] = 5

#: Conjectured GUT scale [GeV] (rough expectation from minimal SUSY GUT).
M_GUT_NOMINAL_GEV: Final[float] = 2.0e16


# ---------------------------------------------------------------------------
# Fermion table used in the QED running sum
# ---------------------------------------------------------------------------

#: Perturbative fermion table for the QED running sum.
#:
#: Each entry: (mass_GeV, charge_in_units_of_e, colour_multiplicity, label).
#: Colour = 1 for leptons, 3 for quarks.
#:
#: Light quarks u, d, s are EXCLUDED: their vacuum-polarization contribution is
#: non-perturbative (they live below Λ_QCD) and is instead captured by the
#: constant ``DELTA_INV_ALPHA_UDS_HAD``.  Only leptons and heavy quarks
#: (c, b, t) are summed perturbatively here.
_FERMION_TABLE: Final[list[tuple[float, float, int, str]]] = [
    (M_E_GEV,    1.0, 1, "electron"),
    (M_MU_GEV,   1.0, 1, "muon"),
    (M_TAU_GEV,  1.0, 1, "tau"),
    (M_CHARM_GEV, 2/3, 3, "charm quark"),
    (M_BOTTOM_GEV, 1/3, 3, "bottom quark"),
    (M_TOP_GEV,   2/3, 3, "top quark"),
]


# ---------------------------------------------------------------------------
# RGRunning — electromagnetic coupling
# ---------------------------------------------------------------------------


@dataclass
class RGRunning:
    """One-loop renormalization-group running of the fine-structure constant α.

    The standard one-loop QED result (sum over all fermions above threshold):

        1/α(Q²) = 1/α(μ²) − (1/(3π)) Σ_f Q_f² N_c^f ln(Q²/m_f²)

    This class provides the inverse coupling, the direct coupling, and the
    one-loop beta function at any energy scale Q [GeV].

    The reference point may be any (Q_ref_GeV, alpha_ref) pair; the default
    is the Thomson limit α(0) = 1/137.035999177 taken at an infinitesimally small
    reference scale (effectively Q_ref ≡ m_e so that all thresholds above m_e
    are counted correctly when running upward).

    Args:
        q_ref_gev: Reference energy scale Q_ref [GeV].  Defaults to m_e so
            that the running starts just above the electron threshold.
        alpha_ref: Value of α at the reference scale.  Defaults to
            ``ALPHA_0`` = 1/137.035999177.

    Attributes:
        q_ref_gev: Reference scale [GeV].
        alpha_ref: α at the reference scale.

    References:
        spec §18.9, §18.34, §18.46; Peskin & Schroeder §12.3.

    Example:
        >>> rg = RGRunning()
        >>> abs(rg.alpha_at_scale(M_Z_GEV) - 1/127.952) < 5e-5
        True
    """

    q_ref_gev: float = field(default=M_E_GEV)
    alpha_ref: float = field(default=ALPHA_0)

    def _delta_inv_alpha(self, q_gev: float) -> float:
        """Vacuum-polarization correction Δ(1/α) from q_ref to q_gev.

        Two contributions:

        1. **Perturbative** (leptons + heavy quarks c, b, t): standard one-loop
           sum -(1/(3π)) Σ_f Q_f² N_c ln(Q²/m_f²) integrated from the higher of
           q_ref and m_f up to q_gev.

        2. **Non-perturbative hadronic** (u, d, s quarks): taken as the constant
           ``DELTA_INV_ALPHA_UDS_HAD`` = −2.8226, applied when q_gev exceeds
           ``Q_HAD_THRESHOLD_GEV`` = 1 GeV.  This is calibrated to reproduce
           the PDG value 1/α(M_Z) = 127.952 exactly.

        Args:
            q_gev: Target energy scale [GeV].

        Returns:
            Δ(1/α) = 1/α(q_gev) − 1/α(q_ref) in dimensionless units.
        """
        delta: float = 0.0

        # Perturbative contribution from leptons and heavy quarks.
        #
        # The RG step must be antisymmetric: running from q_ref -> q and then
        # q -> q_ref should recover the starting coupling.  Integrate only the
        # interval where a fermion is active, with the sign set by direction.
        for m_f, q_f, n_c, _ in _FERMION_TABLE:
            coeff = (1.0 / (3.0 * PI)) * q_f**2 * n_c

            if q_gev > self.q_ref_gev:
                if q_gev > m_f:
                    q_lower = max(self.q_ref_gev, m_f)
                    if q_gev > q_lower:
                        delta -= coeff * math.log(q_gev**2 / q_lower**2)
            elif q_gev < self.q_ref_gev:
                if self.q_ref_gev > m_f:
                    q_lower = max(q_gev, m_f)
                    if self.q_ref_gev > q_lower:
                        delta += coeff * math.log(self.q_ref_gev**2 / q_lower**2)

        # Non-perturbative hadronic VP from light quarks u, d, s.
        # Applied as a step that turns on when q_gev crosses Q_HAD_THRESHOLD_GEV
        # and was not already present at the reference scale.
        ref_above_had = self.q_ref_gev >= Q_HAD_THRESHOLD_GEV
        q_above_had = q_gev >= Q_HAD_THRESHOLD_GEV
        if q_above_had and not ref_above_had:
            delta += DELTA_INV_ALPHA_UDS_HAD
        elif ref_above_had and not q_above_had:
            delta -= DELTA_INV_ALPHA_UDS_HAD  # running downward

        return delta

    def one_loop_inverse_alpha(self, q_gev: float) -> float:
        """Return 1/α(Q²) at energy scale Q = q_gev [GeV].

        Uses the one-loop QED running formula:

            1/α(Q²) = 1/α(μ_ref²) − (1/(3π)) Σ_f Q_f² N_c^f ln(Q²/m_f²)

        where the sum includes all fermions active between the reference scale
        and Q.

        Args:
            q_gev: Energy scale Q [GeV].  Must be positive.

        Returns:
            1/α(Q²), a dimensionless positive number.

        Raises:
            ValueError: If ``q_gev`` is not positive.

        References:
            spec §18.34, §18.46; Peskin & Schroeder eq. (12.54).
        """
        if q_gev <= 0.0:
            raise ValueError(f"q_gev must be positive; got {q_gev!r}")
        return 1.0 / self.alpha_ref + self._delta_inv_alpha(q_gev)

    def alpha_at_scale(self, q_gev: float) -> float:
        """Return α(Q²) directly at energy scale Q = q_gev [GeV].

        This is the reciprocal of :meth:`one_loop_inverse_alpha`.

        Args:
            q_gev: Energy scale Q [GeV].  Must be positive.

        Returns:
            α(Q²) = 1 / [1/α(Q²)], dimensionless.

        Raises:
            ValueError: If ``q_gev`` is not positive.

        References:
            spec §18.34, §18.46.
        """
        return 1.0 / self.one_loop_inverse_alpha(q_gev)

    def beta_function(self) -> float:
        """Return the one-loop QED beta function β = dα/d(log μ) at the reference scale.

        In QED the one-loop beta function (evaluated just above the reference scale,
        where all perturbative fermions lighter than q_ref are active) is:

            β(α) = dα/d(ln μ) = (2α²/(3π)) Σ_f Q_f² N_c^f

        The sum covers leptons and heavy quarks (c, b, t) active at q_ref.
        The light-quark (u, d, s) hadronic contribution is non-perturbative and
        does not admit a simple beta-function expression; it is captured by
        ``DELTA_INV_ALPHA_UDS_HAD``.

        This is always positive in QED — the coupling grows with energy (screening).

        Returns:
            β = dα/d(ln μ) at α = alpha_ref [dimensionless per unit of ln(μ)].

        References:
            spec §18.34; Peskin & Schroeder eq. (12.46).
        """
        sum_qf_sq: float = sum(
            q_f**2 * n_c
            for m_f, q_f, n_c, _ in _FERMION_TABLE
            if m_f <= self.q_ref_gev
        )
        return (2.0 * self.alpha_ref**2 / (3.0 * PI)) * sum_qf_sq


# ---------------------------------------------------------------------------
# run_alpha_to_M_Z convenience function
# ---------------------------------------------------------------------------


def run_alpha_to_M_Z() -> dict[str, float]:
    """Run α from the electron mass to M_Z and compare to the measured value.

    The electron threshold is the lowest threshold in QED; starting from
    α(m_e) = 1/137.035999177 and running upward, all heavier fermions contribute
    vacuum-polarization corrections.

    Returns:
        Dictionary with keys:

        * ``"inv_alpha_MZ_predicted"`` — 1/α(M_Z) from one-loop running.
        * ``"inv_alpha_MZ_measured"`` — PDG value 127.952.
        * ``"absolute_error"`` — |computed − measured|.
        * ``"fractional_error_pct"`` — error as a percentage of measured value.
        * ``"passes"`` — True if the fractional error is < 0.1 %.

    References:
        spec §18.34, §18.46; PDG 2024 review of electroweak model.

    Example:
        >>> result = run_alpha_to_M_Z()
        >>> result["passes"]
        True
    """
    rg = RGRunning(q_ref_gev=M_E_GEV, alpha_ref=ALPHA_0)
    inv_alpha_pred = rg.one_loop_inverse_alpha(M_Z_GEV)
    inv_alpha_meas = INV_ALPHA_MZ_MEASURED
    abs_err = abs(inv_alpha_pred - inv_alpha_meas)
    frac_err_pct = abs_err / inv_alpha_meas * 100.0
    return {
        "inv_alpha_MZ_predicted": inv_alpha_pred,
        "inv_alpha_MZ_measured": inv_alpha_meas,
        "absolute_error": abs_err,
        "fractional_error_pct": frac_err_pct,
        "passes": frac_err_pct < 0.1,
    }


# ---------------------------------------------------------------------------
# StrongCoupling — QCD running
# ---------------------------------------------------------------------------


@dataclass
class StrongCoupling:
    """One-loop renormalization-group running of the QCD coupling α_s(Q²).

    The one-loop QCD beta function gives:

        1/α_s(Q²) = 1/α_s(μ²) + b_0 ln(Q²/μ²)

        b_0 = (33 − 2 N_f) / (12π)

    Asymptotic freedom: b_0 > 0 for N_f ≤ 16 → α_s → 0 as Q → ∞.
    Confinement: α_s diverges when Q → Λ_QCD from above.

    The number of active flavours N_f steps up by 1 at each quark mass
    threshold (decoupling approximation).

    Args:
        alpha_s_ref: α_s at the reference scale.  Defaults to the PDG value
            α_s(M_Z) = 0.1179.
        q_ref_gev: Reference scale [GeV].  Defaults to M_Z = 91.1876 GeV.

    Attributes:
        alpha_s_ref: α_s at the reference scale.
        q_ref_gev: Reference scale [GeV].

    References:
        spec §18.49 (§18.49.3 asymptotic freedom); Peskin & Schroeder §18.3.

    Example:
        >>> sc = StrongCoupling()
        >>> sc.alpha_s_at_scale(10.0) > sc.alpha_s_at_scale(100.0)
        True
    """

    alpha_s_ref: float = field(default=ALPHA_S_MZ)
    q_ref_gev: float = field(default=M_Z_GEV)

    @staticmethod
    def _n_f_at_scale(q_gev: float) -> int:
        """Return number of active quark flavours at scale q_gev [GeV].

        Flavour thresholds (MS-bar masses used):
        * u, d, s: always active above 0 (effectively).
        * c: active above M_CHARM_GEV ≈ 1.28 GeV.
        * b: active above M_BOTTOM_GEV ≈ 4.18 GeV.
        * t: active above M_TOP_GEV ≈ 172.69 GeV.

        Args:
            q_gev: Energy scale [GeV].

        Returns:
            Integer N_f in {3, 4, 5, 6}.
        """
        n_f = 3
        if q_gev >= M_CHARM_GEV:
            n_f += 1
        if q_gev >= M_BOTTOM_GEV:
            n_f += 1
        if q_gev >= M_TOP_GEV:
            n_f += 1
        return n_f

    @staticmethod
    def _b0(n_f: int) -> float:
        """One-loop QCD beta-function coefficient b_0 = (33 − 2 N_f) / (12π).

        Args:
            n_f: Number of active quark flavours.

        Returns:
            b_0 (positive for N_f ≤ 16; ensures asymptotic freedom).

        References:
            spec §18.49.3; Peskin & Schroeder eq. (18.15).
        """
        return (33.0 - 2.0 * n_f) / (12.0 * PI)

    def alpha_s_at_scale(self, q_gev: float) -> float | None:
        """Return α_s(Q²) at energy scale Q = q_gev [GeV].

        Performs a step-by-step threshold crossing from the reference scale
        to the target scale, updating N_f at each quark-mass threshold.
        Returns ``None`` when the target scale is below Λ_QCD (coupling
        formally diverges — confinement regime).

        Args:
            q_gev: Target energy scale [GeV].  Must be positive.

        Returns:
            α_s(Q²) as a dimensionless float, or ``None`` if Q < Λ_QCD.

        Raises:
            ValueError: If ``q_gev`` is not positive.

        References:
            spec §18.49.3.
        """
        if q_gev <= 0.0:
            raise ValueError(f"q_gev must be positive; got {q_gev!r}")
        if q_gev < LAMBDA_QCD_GEV:
            return None  # below confinement scale; coupling non-perturbative

        # Build list of threshold crossings between q_ref and q_gev
        thresholds = [M_STRANGE_GEV, M_CHARM_GEV, M_BOTTOM_GEV, M_TOP_GEV]
        q_lo = self.q_ref_gev
        inv_as = 1.0 / self.alpha_s_ref

        if q_gev >= q_lo:
            # running upward
            for q_thr in sorted(thresholds):
                if q_thr <= q_lo or q_thr >= q_gev:
                    continue
                n_f = self._n_f_at_scale(0.5 * (q_lo + q_thr))
                inv_as += self._b0(n_f) * math.log(q_thr**2 / q_lo**2)
                q_lo = q_thr
            n_f = self._n_f_at_scale(0.5 * (q_lo + q_gev))
            inv_as += self._b0(n_f) * math.log(q_gev**2 / q_lo**2)
        else:
            # running downward
            for q_thr in sorted(thresholds, reverse=True):
                if q_thr >= q_lo or q_thr <= q_gev:
                    continue
                n_f = self._n_f_at_scale(0.5 * (q_thr + q_lo))
                inv_as += self._b0(n_f) * math.log(q_thr**2 / q_lo**2)
                q_lo = q_thr
            n_f = self._n_f_at_scale(0.5 * (q_gev + q_lo))
            inv_as += self._b0(n_f) * math.log(q_gev**2 / q_lo**2)

        if inv_as <= 0.0:
            return None  # Landau pole / confinement
        return 1.0 / inv_as

    def lambda_qcd(self) -> float:
        """Return Λ_QCD [GeV] — the non-perturbative QCD scale.

        Λ_QCD is defined as the scale where 1/α_s(Q) = 0, i.e., where the
        one-loop running coupling formally diverges.  From:

            1/α_s(Λ²) = 0 → ln(Λ²/μ²) = −1/(b_0 α_s(μ²))
            Λ = μ × exp(−1/(2 b_0 α_s(μ²)))

        evaluated with N_f = 5 (appropriate at the M_Z reference scale).

        **One-loop vs two-loop:** the one-loop formula with α_s(M_Z) = 0.1179
        gives Λ_QCD ≈ 87 MeV.  The PDG value of 217 MeV (MS-bar, N_f = 5)
        uses the *two-loop* beta-function solution, which includes an additional
        logarithmic term that shifts Λ_QCD upward by a factor of ~2.5.  The
        one-loop result is internally consistent with the one-loop running used
        throughout this module.

        Returns:
            Λ_QCD [GeV] from the one-loop formula.  About 87 MeV at one loop
            vs the PDG two-loop value of 217 MeV.

        References:
            spec §18.49.5; PDG QCD review 2024.
        """
        n_f = self._n_f_at_scale(self.q_ref_gev)
        b0 = self._b0(n_f)
        return self.q_ref_gev * math.exp(-1.0 / (2.0 * b0 * self.alpha_s_ref))

    def asymptotic_freedom_check(self, q_gev: float) -> dict[str, object]:
        """Verify asymptotic freedom: α_s → 0 as Q → ∞.

        Computes α_s at Q = q_gev and confirms it is smaller than the
        reference value α_s(M_Z).

        Args:
            q_gev: Energy scale [GeV].  Should be > M_Z for the check to
                demonstrate asymptotic freedom.

        Returns:
            Dictionary with keys:

            * ``"q_gev"`` (float): the queried scale.
            * ``"alpha_s"`` (float | None): α_s at that scale.
            * ``"alpha_s_ref"`` (float): the reference value α_s(M_Z).
            * ``"is_smaller"`` (bool): True when alpha_s < alpha_s_ref.
            * ``"freedom_confirmed"`` (bool): True when Q > M_Z and α_s < α_s(M_Z).

        References:
            spec §18.49.3.
        """
        a_s = self.alpha_s_at_scale(q_gev)
        is_smaller = (a_s is not None) and (a_s < self.alpha_s_ref)
        return {
            "q_gev": q_gev,
            "alpha_s": a_s,
            "alpha_s_ref": self.alpha_s_ref,
            "is_smaller": is_smaller,
            "freedom_confirmed": (q_gev > self.q_ref_gev) and is_smaller,
        }

    def confinement_scale(self) -> float:
        """Return the scale [GeV] at which 1/α_s = 0 (confinement / Landau pole).

        This is numerically identical to :meth:`lambda_qcd` — the scale where
        the perturbative one-loop running breaks down and the theory enters the
        non-perturbative confinement regime.

        Returns:
            Λ_QCD [GeV] from the one-loop formula.  Compare to PDG 217 MeV.

        References:
            spec §18.49.5.
        """
        return self.lambda_qcd()


# ---------------------------------------------------------------------------
# Module-level utility functions
# ---------------------------------------------------------------------------


def gut_unification_check(
    m_gut_gev: float,
) -> dict[str, float]:
    """Run all three SM gauge couplings to M_GUT and check approximate unification.

    The three SM couplings are:

    * α_1 = (5/3) α / cos²θ_W  — U(1)_Y hypercharge coupling (GUT normalisation).
    * α_2 = α / sin²θ_W        — SU(2)_L weak isospin coupling.
    * α_3 = α_s                 — SU(3)_c strong coupling.

    At low energy sin²θ_W ≈ 0.2312 [PDG].  All three are run up to m_gut_gev
    using one-loop RG.  In the minimal SUSY GUT (MSSM), they meet within ~1%;
    in the non-SUSY SM they miss by ~10–15%.

    Args:
        m_gut_gev: Candidate GUT unification scale [GeV].

    Returns:
        Dictionary with keys:

        * ``"m_gut_gev"`` — the input scale.
        * ``"inv_alpha1"`` — 1/α_1(M_GUT).
        * ``"inv_alpha2"`` — 1/α_2(M_GUT).
        * ``"inv_alpha3"`` — 1/α_s(M_GUT).
        * ``"max_spread"`` — max|1/α_i − 1/α_j| among all pairs.
        * ``"unification_pct"`` — max_spread / mean × 100 (%).
        * ``"unifies_within_10pct"`` — True when unification_pct < 10.

    References:
        spec §18.34 (SM gauge structure inherited); Peskin & Schroeder §20.4.
    """
    sin2_theta_W: float = 0.23122  # PDG 2024 MS-bar at M_Z

    # U(1)_Y hypercharge coupling at M_Z (GUT normalisation: multiply by 5/3)
    alpha_mz = RGRunning(q_ref_gev=M_E_GEV, alpha_ref=ALPHA_0).alpha_at_scale(M_Z_GEV)
    alpha1_mz = (5.0 / 3.0) * alpha_mz / (1.0 - sin2_theta_W)
    alpha2_mz = alpha_mz / sin2_theta_W

    # SU(2) one-loop running: b_0^SU2 = (22 − 4 N_f/3 − N_H/6) / (4π)
    # With N_f = 3 families and N_H = 1 Higgs doublet:  b_0 = (22-4-1/6)/(4π) ≈ 0.528
    # For U(1): b_0^U1 = -(4 N_f×5/3 + N_H×1/6) / (12π) — negative (Landau pole direction)
    # We use simplified approximate formulae consistent with standard SM running.
    # U(1): 1/α_1(Q) = 1/α_1(M_Z) - b1 ln(Q²/M_Z²)
    # SU(2): 1/α_2(Q) = 1/α_2(M_Z) + b2 ln(Q²/M_Z²)
    # SU(3): handled by StrongCoupling

    # Standard SM one-loop b-coefficients (b_i = coefficient in d(1/α_i)/d ln Q²):
    # b_1 = -41/(12π)  (U(1), N_f=3 gen)
    # b_2 = +19/(12π)  (SU(2), N_f=3 gen, 1 Higgs doublet)
    # b_3 is from QCD (above threshold running)

    b1 = -41.0 / (12.0 * PI)   # U(1) — negative: Landau pole (like QED)
    b2 = 19.0 / (12.0 * PI)    # SU(2) — positive: asymptotically free

    ln_ratio = math.log(m_gut_gev**2 / M_Z_GEV**2)

    inv_alpha1_gut = 1.0 / alpha1_mz - b1 * ln_ratio
    inv_alpha2_gut = 1.0 / alpha2_mz + b2 * ln_ratio

    sc = StrongCoupling()
    a_s_gut = sc.alpha_s_at_scale(m_gut_gev)
    inv_alpha3_gut = (1.0 / a_s_gut) if a_s_gut else float("inf")

    values = [inv_alpha1_gut, inv_alpha2_gut, inv_alpha3_gut]
    mean_val = sum(values) / 3.0
    max_spread = max(abs(v - w) for v in values for w in values)
    unif_pct = max_spread / mean_val * 100.0

    return {
        "m_gut_gev": m_gut_gev,
        "inv_alpha1": inv_alpha1_gut,
        "inv_alpha2": inv_alpha2_gut,
        "inv_alpha3": inv_alpha3_gut,
        "max_spread": max_spread,
        "unification_pct": unif_pct,
        "unifies_within_10pct": unif_pct < 10.0,
    }


def electroweak_unification_scale() -> dict[str, float]:
    """Find the scale where the U(1)_Y and SU(2)_L couplings cross.

    At the electroweak scale the two couplings unify into the GSW model.
    Running both from M_Z, the scale Q_EW where 1/α_1(Q) = 1/α_2(Q) is
    identified.  In practice this is very close to M_Z itself because the
    EW mixing angle already encodes the initial ratio.

    The crossing is found by bisection between M_Z and M_GUT.

    Returns:
        Dictionary with keys:

        * ``"Q_EW_GeV"`` — electroweak unification scale [GeV].
        * ``"alpha_at_EW"`` — common coupling value at Q_EW.
        * ``"inv_alpha_at_EW"`` — 1/α at Q_EW.

    References:
        spec §18.34; Peskin & Schroeder §20.2.
    """
    sin2_theta_W: float = 0.23122
    alpha_mz = RGRunning(q_ref_gev=M_E_GEV, alpha_ref=ALPHA_0).alpha_at_scale(M_Z_GEV)
    alpha1_mz = (5.0 / 3.0) * alpha_mz / (1.0 - sin2_theta_W)
    alpha2_mz = alpha_mz / sin2_theta_W

    b1 = -41.0 / (12.0 * PI)
    b2 = 19.0 / (12.0 * PI)

    def diff(log_q_over_mz: float) -> float:
        inv1 = 1.0 / alpha1_mz - b1 * (2.0 * log_q_over_mz)
        inv2 = 1.0 / alpha2_mz + b2 * (2.0 * log_q_over_mz)
        return inv1 - inv2

    # Bisect: diff changes sign between M_Z (diff > 0) and M_GUT (~2e16 GeV, diff < 0)
    lo, hi = 0.0, math.log(M_GUT_NOMINAL_GEV / M_Z_GEV)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if diff(mid) > 0:
            lo = mid
        else:
            hi = mid
    q_ew = M_Z_GEV * math.exp(0.5 * (lo + hi))
    inv1_ew = 1.0 / alpha1_mz - b1 * math.log(q_ew**2 / M_Z_GEV**2)
    alpha_ew = 1.0 / inv1_ew

    return {
        "Q_EW_GeV": q_ew,
        "alpha_at_EW": alpha_ew,
        "inv_alpha_at_EW": inv1_ew,
    }


def running_via_substrate_bundle(beta_squared: float) -> float:
    """Map a trial sine-Gordon beta to a bare alpha via Coleman bosonization.

    In the substrate-mechanical model (spec §18.9, §18.46) the fine-structure
    constant at zero momentum transfer is:

        α = e² / (4π K ξ⁴)

    where e is the Möbius bundle charge determined by the half-flux holonomy.
    The bundle lives on the 45° cone (§18.10) and is the U(1) sector of the
    §18.45 Lagrangian.

    The sine-Gordon / Thirring duality (Coleman 1975, as used in
    :func:`stiff_medium.mobius_bundle.coleman_bosonization_g`) relates the
    sine-Gordon coupling β² to the bare Thirring coupling g.  The bare QED
    coupling at the confinement scale (where the substrate's kink structure
    crystallizes) is extracted from g via:

        e_bare² / (4π) = g / π  →  α_bare = e_bare² / (4π) = g / π²

    (up to O(1) numerical factors depending on the regularization convention).

    This is a diagnostic map, not a first-principles prediction of alpha. A
    beta value can be chosen so the map lands near 1/137, but that choice is a
    fitted constraint until the missing Möbius-bundle effective-action
    calculation fixes beta independently.

    For the free-fermion point β² = 4π: g = 0, α_bare = 0 (decoupled limit).
    For β² slightly below 4π (the physical kink regime): g > 0 and α_bare > 0.

    The value β² ≈ 4π × (1 − α_0/(2π)) ≈ 12.540 gives α_bare ≈ 1/137 at
    leading order because it was chosen to do so.

    Args:
        beta_squared: Sine-Gordon coupling β².  Must be strictly positive.
            The interesting window is 0 < β² ≤ 8π.

    Returns:
        α_bare = g / π² where g = coleman_bosonization_g(beta_squared).
        This is the current diagnostic estimate of the bare fine-structure
        constant before QED running is applied.

    Raises:
        ValueError: If ``beta_squared`` is not positive (propagated from
            :func:`stiff_medium.mobius_bundle.coleman_bosonization_g`).

    References:
        spec §18.9, §18.46; Coleman (1975) Phys. Rev. D 11, 2088;
        :mod:`stiff_medium.mobius_bundle`.
    """
    g = coleman_bosonization_g(beta_squared)
    # g can be negative above the free-fermion point; clamp that branch to zero.
    if g < 0.0:
        return 0.0
    return g / (PI**2)
