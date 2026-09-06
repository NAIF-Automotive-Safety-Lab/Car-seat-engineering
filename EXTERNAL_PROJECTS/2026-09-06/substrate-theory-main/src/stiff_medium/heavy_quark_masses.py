"""Heavy quark masses (charm, bottom, top) under the substrate K(xi) running.

PURPOSE
-------
Test honestly whether the SAME K(xi) running used for leptons and the QCD
string tension also predicts the heavy quark masses, OR whether heavy quarks
live in a regime where the single power-law breaks down.

ANCHORS (same as substrate_rg_running.py and lepton_masses_from_K.py)
------------------------------------------------------------------
    xi_e   = hbar/(m_e c)             = 3.86e-13 m
    K_e    = hbar c / xi_e^4          = 1.42e24 Pa
    xi_QCD = 0.2 fm                   = 2e-16 m
    K_QCD  = sigma_QCD_SI/sigma_lat   ~ 2.87e5 Pa  (so sigma=0.180 GeV^2)

    K(xi) = K_e * (xi_e / xi)^a,  fitted a ~ -5.69

CRITICAL OBSERVATION OF SIGN
----------------------------
The user prompt asserts a "positive a ~ 5.69".  But anchoring at the
electron and the QCD scales gives:

    K_QCD / K_e = 2.87e5 / 1.42e24 ~ 2e-19   (decreases by 19 orders)
    xi_QCD / xi_e ~ 5e-4                       (decreases by 3 orders)

    a = ln(K_QCD/K_e) / ln(xi_e/xi_QCD)
      = ln(2e-19) / ln(2e3)
      = -43.0 / 7.6
      ~ -5.69     (NEGATIVE)

So in the form  K(xi) = K_e (xi_e/xi)^a  the fitted exponent is a ~ -5.69.
Equivalently:    K(xi) = K_e (xi/xi_e)^|a| =  K_e (xi/xi_e)^{+5.69}
       -- K DECREASES toward small xi between the electron and QCD anchors.

This is OPPOSITE to the textbook QCD picture (alpha_s gets weaker at small
xi).  But it is what the two anchor points enforce when sigma_QCD_GeV2=0.180
and sigma_lattice=0.51 are taken as independent inputs.

CONSEQUENCES FOR HEAVY QUARKS
-----------------------------
With a ~ -5.69 (decreasing K toward small xi), at heavy-quark Compton scales
xi_q << xi_QCD, K(xi_q) is even smaller:

    K(xi_c) ~ K_e * (xi_e/xi_c)^{-5.69} ~ 7e4 Pa
    K(xi_b) ~ 1e2 Pa
    K(xi_t) ~ 5e-8 Pa

These give a string tension sigma(xi_q) = sigma_lat * K(xi_q) that goes to
ZERO at small xi, and a Nambu-Goto mass m_NG = sqrt(sigma * hbar c) that is
ORDER OF MAGNITUDE WRONG for heavy quarks (predicts ~ 0.2 GeV, ~ 7 MeV,
~ 0.0002 MeV for c, b, t respectively, vs. observed 1.27 GeV / 4.18 GeV /
173 GeV).

So the SAME single power-law that connects the electron Compton wavelength
to the QCD string tension does NOT continue smoothly to predict heavy quark
masses through Nambu-Goto sqrt(sigma hbar c).

INTERPRETATION
--------------
There are two honest physical conclusions:

1. The Nambu-Goto m = sqrt(sigma hbar c) formula applies only at the
   confinement scale (light hadrons), where the string is the long-distance
   relevant excitation.  Heavy quarks are characterized by their Compton
   wavelength being short compared to confinement; their mass is set by a
   DIFFERENT mechanism (current-quark mass from electroweak Yukawa, not
   string tension).  In a substrate model, this would correspond to a
   regime change at some xi between xi_QCD and xi_c.

2. Alternatively, the K(xi) running may need a regime change (saturation,
   logarithmic running, or sign change) below xi_QCD.  We test the simple
   power law honestly and report where it breaks.

DELIVERABLES
------------
We compute, at xi_q = hbar c / (m_q c^2) for q in {c, b, t}:

    K(xi_q)         [Pa]   from the single power law
    sigma(xi_q)     [GeV^2]
    m_NG(xi_q)      [GeV]  = sqrt(sigma hbar c)
    m_kink(xi_q)    [GeV]  = 8 hbar/(c xi_q) = 8 m_q   (kinematic, K-independent)

We also try alternative regimes for K below xi_QCD:

    (a) Saturation:   K(xi) = K_QCD                  for xi < xi_QCD
    (b) Log running:  K(xi) = K_QCD * (1 + b ln(xi_QCD/xi))   to test slow rise
    (c) Forced match: solve for K(xi_q) such that m_NG = m_q empirically;
                      report the IMPLIED running exponent a_q on each side.

ASSESSMENT
----------
If sigma(xi_c) is within 10x of sigma_QCD: the running extends to charm.
If sigma(xi_c) is more than 1000x off: regime change is needed.

References
----------
- src/stiff_medium/substrate_rg_running.py        (K(xi) anchored fit)
- src/stiff_medium/lepton_masses_from_K.py        (same anchors, leptons)
- src/stiff_medium/regge_spectrum.py              (sigma_QCD = 0.180 GeV^2)
- spec sections 18.6, 18.32, 18.46, 18.49.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final


# ---------------------------------------------------------------------------
# Physical constants (SI; PDG / CODATA 2024)
# ---------------------------------------------------------------------------

C_SI: Final[float] = 2.99792458e8           # m/s
HBAR_SI: Final[float] = 1.054571817e-34     # J s
M_E_SI: Final[float] = 9.1093837015e-31     # kg
EV_TO_J: Final[float] = 1.602176634e-19     # J
MEV_TO_J: Final[float] = EV_TO_J * 1e6
GEV_TO_J: Final[float] = EV_TO_J * 1e9
FM_TO_M: Final[float] = 1e-15

# Electron and QCD anchors (consistent with substrate_rg_running.py)
XI_ELECTRON: Final[float] = HBAR_SI / (M_E_SI * C_SI)            # 3.862e-13 m
K_ELECTRON: Final[float] = HBAR_SI * C_SI / XI_ELECTRON ** 4     # 1.42e24 Pa

XI_QCD: Final[float] = 0.2 * FM_TO_M                              # 2e-16 m
SIGMA_QCD_GEV2: Final[float] = 0.180                              # GeV^2 (lattice)
SIGMA_LATTICE_NATURAL: Final[float] = 0.51                        # dimensionless
SIGMA_QCD_SI: Final[float] = (
    SIGMA_QCD_GEV2 * GEV_TO_J ** 2 / (HBAR_SI * C_SI)
)
K_QCD_TARGET: Final[float] = SIGMA_QCD_SI / SIGMA_LATTICE_NATURAL  # ~2.87e5 Pa

# PDG 2024 heavy-quark current masses (MS-bar at 2 GeV for c,b; pole for t)
M_C_GEV: Final[float] = 1.27       # MS-bar
M_B_GEV: Final[float] = 4.18       # MS-bar
M_T_GEV: Final[float] = 173.0      # pole-ish


# ---------------------------------------------------------------------------
# Fitted running exponent (anchored at electron and QCD scales)
# ---------------------------------------------------------------------------

def fit_running_exponent_a() -> float:
    """Compute a in K(xi) = K_e * (xi_e / xi)^a from the two anchors.

    Algebra:
        K_QCD / K_e = (xi_e / xi_QCD)^a
        a = ln(K_QCD / K_e) / ln(xi_e / xi_QCD)

    Returns:
        Fitted a (negative when K decreases toward small xi).
    """
    return math.log(K_QCD_TARGET / K_ELECTRON) / math.log(XI_ELECTRON / XI_QCD)


def K_running(xi: float, a: float | None = None) -> float:
    """Power-law running stiffness K(xi) [Pa] anchored at xi_e.

    Args:
        xi: Length scale [m].
        a:  Optional override of the fitted exponent (default: anchored fit).

    Returns:
        K(xi) [Pa].
    """
    if a is None:
        a = fit_running_exponent_a()
    return K_ELECTRON * (XI_ELECTRON / xi) ** a


# ---------------------------------------------------------------------------
# Quark Compton scales
# ---------------------------------------------------------------------------

def xi_compton_GeV(m_GeV: float) -> float:
    """Compton wavelength xi = hbar c / (m c^2) for mass m [GeV].

    Args:
        m_GeV: Mass in GeV.

    Returns:
        xi [m].
    """
    return HBAR_SI * C_SI / (m_GeV * GEV_TO_J)


# ---------------------------------------------------------------------------
# Substrate predictions at heavy-quark Compton scales
# ---------------------------------------------------------------------------

@dataclass
class HeavyQuarkPrediction:
    """Substrate-derived quantities at one heavy-quark Compton scale.

    Attributes:
        name:                Quark name ("c", "b", "t").
        m_emp_GeV:           Empirical current-quark mass [GeV] (PDG).
        xi_q:                xi = hbar c / (m_q c^2) [m].
        K_at_xi_q:           K(xi_q) from the single-power running [Pa].
        sigma_at_xi_q_GeV2:  sigma(xi_q) = sigma_lat * K(xi_q) in GeV^2.
        sigma_ratio_to_QCD:  sigma(xi_q) / sigma_QCD.
        m_NG_GeV:            sqrt(sigma * hbar c) at xi_q, in GeV.
        m_kink_kin_GeV:      8 hbar/(c xi_q) = 8 m_q (kinematic, no K).
        m_NG_over_m_emp:     ratio m_NG / m_emp (1.0 = perfect prediction).
        K_emp_required_Pa:   K(xi_q) required to make m_NG = m_emp.
        a_implied_emp:       Implied exponent a' if K(xi_q) = K_emp_required.
    """

    name: str
    m_emp_GeV: float
    xi_q: float
    K_at_xi_q: float
    sigma_at_xi_q_GeV2: float
    sigma_ratio_to_QCD: float
    m_NG_GeV: float
    m_kink_kin_GeV: float
    m_NG_over_m_emp: float
    K_emp_required_Pa: float
    a_implied_emp: float


def predict_heavy_quark(name: str, m_emp_GeV: float, a: float | None = None
                        ) -> HeavyQuarkPrediction:
    """Compute substrate predictions for one heavy quark.

    Args:
        name:       Quark label.
        m_emp_GeV:  Empirical mass in GeV.
        a:          Optional override of running exponent.

    Returns:
        HeavyQuarkPrediction.
    """
    if a is None:
        a = fit_running_exponent_a()

    xi_q = xi_compton_GeV(m_emp_GeV)
    K_q = K_running(xi_q, a=a)

    sigma_SI = SIGMA_LATTICE_NATURAL * K_q
    sigma_GeV2 = sigma_SI * HBAR_SI * C_SI / GEV_TO_J ** 2
    sigma_ratio = sigma_GeV2 / SIGMA_QCD_GEV2

    # Nambu-Goto mass: m c^2 = sqrt(sigma_SI * hbar c)
    m_NG_J = math.sqrt(abs(sigma_SI) * HBAR_SI * C_SI)
    m_NG_GeV = m_NG_J / GEV_TO_J

    # Kinematic kink mass
    m_kink_kg = 8.0 * HBAR_SI / (C_SI * xi_q)
    m_kink_GeV = m_kink_kg * C_SI ** 2 / GEV_TO_J

    # K required at xi_q so that m_NG = m_emp:
    #   sigma_SI_req = (m_emp c^2)^2 / (hbar c)
    #   K_req = sigma_SI_req / sigma_lat
    sigma_req_SI = (m_emp_GeV * GEV_TO_J) ** 2 / (HBAR_SI * C_SI)
    K_req = sigma_req_SI / SIGMA_LATTICE_NATURAL

    # Implied exponent a' if K_e * (xi_e/xi_q)^{a'} = K_req
    if K_req > 0 and xi_q != XI_ELECTRON:
        a_imp = math.log(K_req / K_ELECTRON) / math.log(XI_ELECTRON / xi_q)
    else:
        a_imp = float("nan")

    return HeavyQuarkPrediction(
        name=name,
        m_emp_GeV=m_emp_GeV,
        xi_q=xi_q,
        K_at_xi_q=K_q,
        sigma_at_xi_q_GeV2=sigma_GeV2,
        sigma_ratio_to_QCD=sigma_ratio,
        m_NG_GeV=m_NG_GeV,
        m_kink_kin_GeV=m_kink_GeV,
        m_NG_over_m_emp=m_NG_GeV / m_emp_GeV,
        K_emp_required_Pa=K_req,
        a_implied_emp=a_imp,
    )


# ---------------------------------------------------------------------------
# Alternative regimes for K below xi_QCD
# ---------------------------------------------------------------------------

def K_saturated(xi: float, K_sat: float = K_QCD_TARGET,
                xi_sat: float = XI_QCD,
                a_above: float | None = None) -> float:
    """Saturated K: piecewise running.

    For xi >= xi_sat: K = K_running(xi, a_above)
    For xi <  xi_sat: K = K_sat (constant; no further softening below xi_QCD).

    Args:
        xi:        Length scale [m].
        K_sat:     Saturation value [Pa].
        xi_sat:    Saturation cutoff [m] (default xi_QCD).
        a_above:   Running exponent above the cutoff (default: anchored fit).

    Returns:
        K(xi) [Pa].
    """
    if xi >= xi_sat:
        return K_running(xi, a=a_above)
    return K_sat


def K_log_running(xi: float, K_QCD: float = K_QCD_TARGET,
                  xi_QCD_arg: float = XI_QCD,
                  b: float = 0.0) -> float:
    """Logarithmic running below xi_QCD: K(xi) = K_QCD * (1 + b * ln(xi_QCD/xi)).

    Mirrors a "QCD-like" beta function with slow logarithmic increase.
    For b=0 -> saturated.  For b>0 -> grows slowly toward small xi.

    Args:
        xi:          Length scale [m].
        K_QCD:       K at xi_QCD [Pa].
        xi_QCD_arg:  Reference xi [m].
        b:           Slope coefficient (dimensionless).

    Returns:
        K(xi) [Pa].
    """
    L = math.log(xi_QCD_arg / xi)  # > 0 if xi < xi_QCD
    return K_QCD * (1.0 + b * L)


def fit_b_log_to_anchor(xi_target: float, K_target: float,
                         K_QCD: float = K_QCD_TARGET,
                         xi_QCD_arg: float = XI_QCD) -> float:
    """Fit log-running coefficient b so K_log_running(xi_target) = K_target.

    Args:
        xi_target:    Desired xi [m].
        K_target:     Desired K(xi_target) [Pa].
        K_QCD:        Anchor K(xi_QCD) [Pa].
        xi_QCD_arg:   Anchor xi [m].

    Returns:
        b.
    """
    L = math.log(xi_QCD_arg / xi_target)
    if L == 0.0:
        return 0.0
    return ((K_target / K_QCD) - 1.0) / L


# ---------------------------------------------------------------------------
# Heavy-quark mass scaling: sqrt(sigma) and other dimensional combos
# ---------------------------------------------------------------------------

def m_from_sqrt_sigma(sigma_GeV2: float) -> float:
    """Mass scale from string tension: m c^2 = sqrt(sigma hbar c).

    Note: sigma is in GeV^2, so in natural units sqrt(sigma) is in GeV.

    Args:
        sigma_GeV2: String tension [GeV^2].

    Returns:
        m [GeV].
    """
    return math.sqrt(abs(sigma_GeV2))


def predict_m_q_via_running_only(m_emp_GeV: float,
                                 a: float | None = None
                                 ) -> tuple[float, float, float]:
    """Predict m_q from K(xi_q) only, using sqrt(sigma) scaling.

    Computes (m_NG, sigma_GeV2, K_at_xi_q).

    Args:
        m_emp_GeV: Empirical mass (used to fix the Compton scale).
        a:         Running exponent.

    Returns:
        (m_NG_GeV, sigma_GeV2, K_at_xi_q_Pa).
    """
    if a is None:
        a = fit_running_exponent_a()
    xi_q = xi_compton_GeV(m_emp_GeV)
    K_q = K_running(xi_q, a=a)
    sigma_SI = SIGMA_LATTICE_NATURAL * K_q
    sigma_GeV2 = sigma_SI * HBAR_SI * C_SI / GEV_TO_J ** 2
    m_NG_GeV = math.sqrt(abs(sigma_SI) * HBAR_SI * C_SI) / GEV_TO_J
    return m_NG_GeV, sigma_GeV2, K_q


# ---------------------------------------------------------------------------
# Summary container
# ---------------------------------------------------------------------------

@dataclass
class HeavyQuarkSummary:
    """Complete summary of substrate predictions for c, b, t.

    Attributes:
        a_fit:                       Fitted exponent a from electron-QCD anchors.
        predictions:                 List of HeavyQuarkPrediction (c, b, t).
        sigma_at_charm_to_QCD_ratio: sigma(xi_c) / sigma_QCD.
        regime_change_required:      True if sigma(xi_c) is more than 1000x off.
        verdict:                     Short verdict string.
        b_log_to_charm:              Log-running coefficient that would force charm.
        b_log_to_bottom:             Log-running coefficient that would force bottom.
        notes:                       Multi-line honest assessment text.
    """

    a_fit: float
    predictions: list[HeavyQuarkPrediction]
    sigma_at_charm_to_QCD_ratio: float
    regime_change_required: bool
    verdict: str
    b_log_to_charm: float
    b_log_to_bottom: float
    notes: list[str] = field(default_factory=list)


def run_full_analysis() -> HeavyQuarkSummary:
    """Run the full substrate heavy-quark mass analysis.

    Returns:
        HeavyQuarkSummary.
    """
    a = fit_running_exponent_a()

    preds: list[HeavyQuarkPrediction] = []
    for name, m_GeV in (("c", M_C_GEV), ("b", M_B_GEV), ("t", M_T_GEV)):
        preds.append(predict_heavy_quark(name, m_GeV, a=a))

    sigma_c_ratio = preds[0].sigma_ratio_to_QCD

    # The honest test is whether m_NG(xi_q) ~ m_emp, not just sigma at charm.
    # Heavy quarks fail badly even where sigma is "only" 4x off.
    worst_mass_ratio = max(abs(math.log10(p.m_NG_over_m_emp))
                           for p in preds if p.m_NG_over_m_emp > 0)
    regime_change = worst_mass_ratio > 0.5  # off by more than ~3x

    if not regime_change:
        verdict = (
            "K(xi) running predicts heavy-quark masses within a factor ~3: "
            f"worst |log10(m_NG/m_emp)| = {worst_mass_ratio:.2f}."
        )
    else:
        n_decades = worst_mass_ratio
        verdict = (
            "K(xi) single-power running does NOT predict heavy-quark masses: "
            f"m_NG/m_emp varies from {preds[0].m_NG_over_m_emp:.2e} (c) "
            f"to {preds[-1].m_NG_over_m_emp:.2e} (t); worst is off by "
            f"{n_decades:.1f} decades. The Nambu-Goto formula sqrt(sigma hbar c) "
            "fails for heavy quarks; their masses are a separate input."
        )

    # Log-running coefficient that would force m_NG(xi_c) = m_c (and m_b)
    K_req_c = preds[0].K_emp_required_Pa
    K_req_b = preds[1].K_emp_required_Pa
    b_c = fit_b_log_to_anchor(preds[0].xi_q, K_req_c)
    b_b = fit_b_log_to_anchor(preds[1].xi_q, K_req_b)

    notes = [
        "HONEST ASSESSMENT OF SUBSTRATE K(xi) AT HEAVY-QUARK SCALES",
        "=" * 60,
        "",
        f"Fitted running exponent a = {a:.4f} (NEGATIVE: K decreases toward small xi)",
        "Equivalently  K(xi) = K_e * (xi/xi_e)^{|a|}  with |a| ~ 5.69.",
        "",
        "PER-QUARK PREDICTIONS (single-power-law continuation):",
    ]
    for p in preds:
        notes += [
            f"  {p.name} (m_emp = {p.m_emp_GeV:.3f} GeV, xi_q = {p.xi_q:.3e} m):",
            f"    K(xi_q)        = {p.K_at_xi_q:.3e} Pa",
            f"    sigma(xi_q)    = {p.sigma_at_xi_q_GeV2:.3e} GeV^2  "
            f"(ratio to sigma_QCD = {p.sigma_ratio_to_QCD:.3e})",
            f"    m_NG(xi_q)     = {p.m_NG_GeV:.3e} GeV  "
            f"(observed {p.m_emp_GeV:.3f}, ratio m_NG/m_emp = {p.m_NG_over_m_emp:.3e})",
            f"    m_kink_kin     = {p.m_kink_kin_GeV:.3f} GeV  (= 8 m_emp; trivial)",
            f"    K required for m_NG = m_emp: {p.K_emp_required_Pa:.3e} Pa  "
            f"(implied a' = {p.a_implied_emp:.3f})",
        ]

    notes += [
        "",
        f"VERDICT: {verdict}",
        "",
        "ALTERNATIVE REGIME — LOG RUNNING BELOW xi_QCD:",
        "  The simplest fix is a regime change at xi_QCD:",
        "    K(xi) = K_QCD * (1 + b ln(xi_QCD/xi))   for xi < xi_QCD",
        f"  Forcing m_NG(xi_c) = m_c gives b = {b_c:.3e} (huge).",
        f"  Forcing m_NG(xi_b) = m_b gives b = {b_b:.3e}.",
        "  These b values are NOT consistent with each other (would not be a",
        "  single new beta-function but separate ad-hoc adjustments). So:",
        "  a single log-running cannot fix charm and bottom simultaneously",
        "  using sqrt(sigma hbar c) as the mass formula.",
        "",
        "WHAT WOULD BE NEEDED TO PREDICT HEAVY-QUARK MASSES IN THE SUBSTRATE",
        "  (a) A different mass formula above the QCD scale.  In real QCD,",
        "      heavy-quark masses are CURRENT (Lagrangian) masses set by",
        "      Yukawa couplings to the Higgs, not by the string tension.",
        "      The string tension only sets the BINDING energy of mesons,",
        "      not the heavy-quark mass itself.",
        "  (b) A piecewise K(xi) with different running on either side of xi_QCD.",
        "  (c) The Higgs-Yukawa sector, which is OUTSIDE the K-substrate model",
        "      as currently formulated.",
        "",
        "CONCLUSION (HONEST):",
        "  The single power-law K(xi) anchored at the electron Compton scale",
        "  and the QCD string tension does NOT predict m_c, m_b, m_t.",
        "  The substrate model as it stands captures the LOW-ENERGY hadronic",
        "  string tension and the lepton Compton hierarchy, but heavy-quark",
        "  current masses are an INDEPENDENT input that the substrate does",
        "  not produce from the K(xi) running alone.",
        "  This is a real boundary of the model, not a tuning failure.",
    ]

    return HeavyQuarkSummary(
        a_fit=a,
        predictions=preds,
        sigma_at_charm_to_QCD_ratio=sigma_c_ratio,
        regime_change_required=regime_change,
        verdict=verdict,
        b_log_to_charm=b_c,
        b_log_to_bottom=b_b,
        notes=notes,
    )
