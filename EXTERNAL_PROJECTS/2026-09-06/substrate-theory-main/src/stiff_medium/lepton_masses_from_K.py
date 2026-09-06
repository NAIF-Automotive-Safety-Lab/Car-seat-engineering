"""Lepton mass scale from the substrate K(xi) running.

Cross-scale consistency test for leptons (analogue of the meson Regge test).

Setup
-----
We use the SAME running substrate stiffness used for the meson trajectory:

    K(xi) = K_e * (xi_e / xi)^a,    a = 5.69

anchored at the electron Compton scale

    xi_e = hbar/(m_e c)         = 3.86e-13 m
    K_e  = hbar*c / xi_e^4      = 1.42e24 Pa

and at the QCD scale

    xi_QCD = 0.2 fm = 2e-16 m
    K_QCD  ~ 1.1e34 Pa  (so that sigma = sigma_lat * K_QCD = 0.180 GeV^2).

The Foot parametrization (Foot 1994) places ALL THREE charged leptons on
the Koide surface Q = 2/3 EXACTLY for any (M, delta) via

    m_n = M * (1 + sqrt(2) cos(2*pi*n/3 + delta))^2,   n = 0, 1, 2.

The empirical phase delta_emp ~ 1.404*pi sets the lepton mass HIERARCHY.
The empirical mass SCALE M_Foot is fixed by the geometric mean construction

    M_Foot = ((sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau)) / 3)^2.

QUESTION ANSWERED HERE
----------------------
Does the substrate K(xi) running produce M_Foot from first principles, or
is M_Foot an additional input not contained in the running stiffness alone?

We compute the following candidates and compare to M_Foot:

  (A) m_kink(xi_e)                   = 8 hbar/(c xi_e) = 8 m_e (~4.09 MeV).
  (B) sqrt(K_e * xi_e^3 * m_e c^2)  - dimensional combination at xi_e.
  (C) m_kink(xi*) such that m_kink(xi*) = M_Foot.  Solve for xi* and report
      whether xi* lies between xi_QCD (2e-16 m) and xi_e (3.86e-13 m).
      With the running K(xi) = K_e * (xi_e/xi)^a we have

          m_kink(xi) = 8 hbar/(c xi)         (kink-mass formula, model-independent of K)
                     = 8 K(xi) xi^3 / c^2    (after K xi^4 = hbar c is enforced).

      Note: enforcing the SUBSTRATE primitive hbar = K xi^4 / c at every
      scale would force K_e * xi_e^4 = K(xi) * xi^4, i.e. a UNIQUE running
      a = 4 (kinematic).  The phenomenological a = 5.69 is INCOMPATIBLE
      with hbar = K xi^4 / c at every scale; it is a *coarse-grained*
      stiffness, not the action quantum constraint.

      Therefore we have TWO inequivalent ways to extract m_kink(xi):

        (C1) Kinematic:  m_kink_kin(xi) = 8 hbar/(c xi)
                         (uses xi only; K plays no role).
        (C2) Phenomenological:  m_kink_phen(xi) = 8 K(xi) xi^3 / c^2
                         (uses the running K with a = 5.69; the result
                          differs from C1 except at the anchor xi_e).

      We report both, since the question is precisely whether the running
      K is "doing the work" or whether the kink mass is fixed by hbar alone.

  (D) Nambu-Goto-like at xi_tau = hbar/(m_tau c) = 0.111 fm:
            sigma(xi_tau) = sigma_lat * K(xi_tau)
            m_NG(xi_tau)  = sqrt(sigma(xi_tau) * hbar c)
      and compare to m_tau.

References
----------
- src/stiff_medium/substrate_primitives.py        (K_e, xi_e, hbar = K xi^4/c)
- src/stiff_medium/substrate_rg_running.py        (K(xi) running, a fit)
- src/stiff_medium/foot_phase_derivation.py       (Foot parametrization)
- src/stiff_medium/regge_spectrum.py              (sigma_QCD = 0.180 GeV^2)
- spec sections: 18.30 (stress-loaded vertex), 18.35 (m c^2 = hbar sqrt(kappa/I)).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final


# ---------------------------------------------------------------------------
# Physical constants (SI; PDG / CODATA)
# ---------------------------------------------------------------------------

C_SI: Final[float] = 2.99792458e8           # m/s
HBAR_SI: Final[float] = 1.054571817e-34     # J s
M_E_SI: Final[float] = 9.1093837015e-31     # kg
EV_TO_J: Final[float] = 1.602176634e-19     # J
MEV_TO_J: Final[float] = EV_TO_J * 1e6      # J
GEV_TO_J: Final[float] = EV_TO_J * 1e9      # J
FM_TO_M: Final[float] = 1e-15               # m

# PDG 2024 charged-lepton masses (MeV)
M_E_MEV: Final[float] = 0.51099895
M_MU_MEV: Final[float] = 105.6583755
M_TAU_MEV: Final[float] = 1776.86

# Substrate primitives at the electron Compton scale
XI_ELECTRON: Final[float] = HBAR_SI / (M_E_SI * C_SI)            # 3.8616e-13 m
K_ELECTRON: Final[float] = HBAR_SI * C_SI / XI_ELECTRON ** 4     # 1.42e24 Pa

# QCD anchor scale and stiffness target (so sigma_SI = sigma_lat * K_QCD = sigma_QCD)
XI_QCD: Final[float] = 0.2 * FM_TO_M                              # 2e-16 m
SIGMA_QCD_GEV2: Final[float] = 0.180                              # GeV^2 (lattice QCD)
SIGMA_LATTICE_NATURAL: Final[float] = 0.51                        # dimensionless, model
SIGMA_QCD_SI: Final[float] = (
    SIGMA_QCD_GEV2 * GEV_TO_J ** 2 / (HBAR_SI * C_SI)
)  # J/m
K_QCD_TARGET: Final[float] = SIGMA_QCD_SI / SIGMA_LATTICE_NATURAL  # ~1.1e34 Pa


# ---------------------------------------------------------------------------
# Running K(xi)
# ---------------------------------------------------------------------------

def fit_running_exponent() -> float:
    """Compute the power-law exponent a in K(xi) = K_e (xi_e / xi)^a.

    Fixed by anchoring K(xi_e) = K_e and K(xi_QCD) = K_QCD_TARGET.

    Returns:
        a (positive; K grows as xi shrinks).
    """
    # K_QCD = K_e * (xi_e / xi_QCD) ^ a
    return math.log(K_QCD_TARGET / K_ELECTRON) / math.log(XI_ELECTRON / XI_QCD)


def K_running(xi: float) -> float:
    """Return K(xi) [Pa] from the two-anchor power law.

    Args:
        xi: Length scale [m].

    Returns:
        K(xi) [Pa].
    """
    a = fit_running_exponent()
    return K_ELECTRON * (XI_ELECTRON / xi) ** a


# ---------------------------------------------------------------------------
# Foot parametrization scale M_Foot
# ---------------------------------------------------------------------------

def foot_scale_M(m_e_MeV: float = M_E_MEV,
                 m_mu_MeV: float = M_MU_MEV,
                 m_tau_MeV: float = M_TAU_MEV) -> float:
    """Empirical Foot scale M from the three lepton masses (MeV).

    M_Foot = ((sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau)) / 3)^2

    The geometric-mean construction is forced by the Foot parametrization:
    Sigma_n sqrt(m_n) = 3 sqrt(M_Foot), since the n-sum of (1 + sqrt(2) cos)
    equals 3.

    Args:
        m_e_MeV: m_e in MeV.
        m_mu_MeV: m_mu in MeV.
        m_tau_MeV: m_tau in MeV.

    Returns:
        M_Foot [MeV].
    """
    s = math.sqrt(m_e_MeV) + math.sqrt(m_mu_MeV) + math.sqrt(m_tau_MeV)
    return (s / 3.0) ** 2


# ---------------------------------------------------------------------------
# Mass candidates from substrate quantities
# ---------------------------------------------------------------------------

def m_kink_kinematic_MeV(xi: float) -> float:
    """Kinematic kink mass at scale xi:  m c^2 = 8 hbar c / xi.

    This formula uses NO information about K(xi); it follows from the
    sine-Gordon kink mass m_kink = 8 m_0 with m_0 c^2 = hbar c / xi.

    Args:
        xi: Length scale [m].

    Returns:
        m_kink [MeV].
    """
    m_kink_kg = 8.0 * HBAR_SI / (C_SI * xi)
    return m_kink_kg * C_SI ** 2 / MEV_TO_J


def m_kink_phenomenological_MeV(xi: float) -> float:
    """Phenomenological kink mass at scale xi using the running K:
       m c^2 = 8 K(xi) xi^3 / c^2 * c^2 = 8 K(xi) xi^3.

    Identity check:  if K xi^4 = hbar c (the substrate primitive holds),
    then 8 K xi^3 = 8 hbar c / xi -> identical to m_kink_kinematic.

    With the phenomenological running a = 5.69, K xi^4 is NOT constant
    away from the electron anchor, so the two formulas diverge.

    Args:
        xi: Length scale [m].

    Returns:
        m_kink [MeV].
    """
    K = K_running(xi)
    energy_J = 8.0 * K * xi ** 3
    return energy_J / MEV_TO_J


def M_dimensional_combination_MeV() -> float:
    """Candidate M = sqrt(K_e * xi_e^3 * m_e c^2)  -- dimensional combo at xi_e.

    Dimensional check:
        K_e: J/m^3
        xi_e^3: m^3
        m_e c^2: J
        product: J^2  -> sqrt: J  -> energy.

    Note that  K_e * xi_e^3 = (hbar c / xi_e^4) * xi_e^3 = hbar c / xi_e
    = m_e c^2,  so this candidate is identically  sqrt(m_e c^2 * m_e c^2)
    = m_e c^2 = 0.511 MeV.  We compute it to confirm the algebra.

    Returns:
        M [MeV].
    """
    energy_J_squared = K_ELECTRON * XI_ELECTRON ** 3 * (M_E_SI * C_SI ** 2)
    return math.sqrt(energy_J_squared) / MEV_TO_J


def find_xi_for_M_kinematic_MeV(M_target_MeV: float) -> float:
    """Solve  m_kink_kinematic(xi*) = M_target  -> xi*.

    Closed form:
        m c^2 = 8 hbar c / xi  =>  xi = 8 hbar c / (M_target * MeV_to_J).

    Args:
        M_target_MeV: Target mass [MeV].

    Returns:
        xi* [m].
    """
    M_target_J = M_target_MeV * MEV_TO_J
    return 8.0 * HBAR_SI * C_SI / M_target_J


def find_xi_for_M_phenomenological_MeV(M_target_MeV: float) -> float:
    """Solve  m_kink_phen(xi*) = M_target  -> xi*  using the running K.

    With K(xi) = K_e (xi_e / xi)^a:
        m c^2 [J] = 8 K(xi) xi^3
                  = 8 K_e xi_e^a xi^(3-a)
        =>  xi^(3-a) = (M_target_J) / (8 K_e xi_e^a)

    For a > 3 the exponent is negative, but xi can still be solved.

    Args:
        M_target_MeV: Target mass [MeV].

    Returns:
        xi* [m].
    """
    a = fit_running_exponent()
    M_target_J = M_target_MeV * MEV_TO_J
    rhs = M_target_J / (8.0 * K_ELECTRON * XI_ELECTRON ** a)
    return rhs ** (1.0 / (3.0 - a))


def nambu_goto_mass_at_xi_tau_MeV() -> tuple[float, float, float]:
    """Compute m_NG = sqrt(sigma(xi_tau) * hbar c) at xi_tau = hbar/(m_tau c).

    sigma(xi_tau) = sigma_lat * K(xi_tau)

    Returns:
        Tuple (xi_tau [m], sigma_GeV2, m_NG [MeV]).
    """
    m_tau_kg = M_TAU_MEV * MEV_TO_J / C_SI ** 2
    xi_tau = HBAR_SI / (m_tau_kg * C_SI)
    sigma_SI = SIGMA_LATTICE_NATURAL * K_running(xi_tau)
    sigma_GeV2 = sigma_SI * (HBAR_SI * C_SI) / GEV_TO_J ** 2
    m_NG_J = math.sqrt(abs(sigma_SI) * HBAR_SI * C_SI)
    m_NG_MeV = m_NG_J / MEV_TO_J
    return xi_tau, sigma_GeV2, m_NG_MeV


# ---------------------------------------------------------------------------
# Bundled report dataclass
# ---------------------------------------------------------------------------

@dataclass
class LeptonMassReport:
    """All numerical results from the cross-scale lepton-scale test.

    Attributes:
        a_running:                 Power-law exponent a in K(xi) = K_e (xi_e/xi)^a.
        M_foot_MeV:                Empirical Foot scale M [MeV].
        M_kink_xi_e_MeV:           Kinematic kink mass at xi_e [MeV] (= 8 m_e).
        M_dimensional_MeV:         sqrt(K_e xi_e^3 m_e c^2) [MeV] (= m_e identically).
        xi_star_kinematic_m:       Scale xi* where m_kink_kin = M_foot [m].
        xi_star_kinematic_fm:      xi* in fm.
        xi_star_kin_in_QCD_band:   bool: xi_QCD <= xi_star_kin <= xi_e?
        xi_star_phen_m:            Scale xi* with the running K [m].
        xi_star_phen_fm:           xi* (phenomenological) in fm.
        xi_star_phen_in_QCD_band:  bool for the phenomenological xi*.
        m_kink_phen_at_xi_QCD_MeV: m_kink_phen(xi_QCD) [MeV] (sanity check).
        m_kink_kin_at_xi_QCD_MeV:  m_kink_kin(xi_QCD) [MeV].
        xi_tau_m:                  hbar/(m_tau c) [m].
        sigma_at_xi_tau_GeV2:      sigma at xi_tau [GeV^2].
        m_NG_at_xi_tau_MeV:        Nambu-Goto mass sqrt(sigma * hbar c) at xi_tau [MeV].
        ratio_M_kink_to_M_foot:    M_kink(xi_e) / M_foot.
        ratio_m_NG_to_m_tau:       m_NG_at_xi_tau / m_tau.
    """

    a_running: float
    M_foot_MeV: float
    M_kink_xi_e_MeV: float
    M_dimensional_MeV: float
    xi_star_kinematic_m: float
    xi_star_kinematic_fm: float
    xi_star_kin_in_QCD_band: bool
    xi_star_phen_m: float
    xi_star_phen_fm: float
    xi_star_phen_in_QCD_band: bool
    m_kink_phen_at_xi_QCD_MeV: float
    m_kink_kin_at_xi_QCD_MeV: float
    xi_tau_m: float
    sigma_at_xi_tau_GeV2: float
    m_NG_at_xi_tau_MeV: float
    ratio_M_kink_to_M_foot: float
    ratio_m_NG_to_m_tau: float


def run_full_analysis() -> LeptonMassReport:
    """Run the complete K(xi) -> lepton-scale cross-test and return a report.

    Returns:
        LeptonMassReport with all numerical results.
    """
    a = fit_running_exponent()

    M_foot = foot_scale_M()                                        # MeV
    M_kink_e = m_kink_kinematic_MeV(XI_ELECTRON)                   # MeV
    M_dim = M_dimensional_combination_MeV()                        # MeV

    xi_star_kin = find_xi_for_M_kinematic_MeV(M_foot)              # m
    xi_star_phen = find_xi_for_M_phenomenological_MeV(M_foot)      # m

    in_band_kin = (XI_QCD <= xi_star_kin <= XI_ELECTRON)
    in_band_phen = (XI_QCD <= xi_star_phen <= XI_ELECTRON)

    # Sanity: kink mass at the QCD scale under both prescriptions
    m_kink_phen_qcd = m_kink_phenomenological_MeV(XI_QCD)
    m_kink_kin_qcd = m_kink_kinematic_MeV(XI_QCD)

    xi_tau, sigma_xi_tau_GeV2, m_NG_tau = nambu_goto_mass_at_xi_tau_MeV()

    return LeptonMassReport(
        a_running=a,
        M_foot_MeV=M_foot,
        M_kink_xi_e_MeV=M_kink_e,
        M_dimensional_MeV=M_dim,
        xi_star_kinematic_m=xi_star_kin,
        xi_star_kinematic_fm=xi_star_kin / FM_TO_M,
        xi_star_kin_in_QCD_band=in_band_kin,
        xi_star_phen_m=xi_star_phen,
        xi_star_phen_fm=xi_star_phen / FM_TO_M,
        xi_star_phen_in_QCD_band=in_band_phen,
        m_kink_phen_at_xi_QCD_MeV=m_kink_phen_qcd,
        m_kink_kin_at_xi_QCD_MeV=m_kink_kin_qcd,
        xi_tau_m=xi_tau,
        sigma_at_xi_tau_GeV2=sigma_xi_tau_GeV2,
        m_NG_at_xi_tau_MeV=m_NG_tau,
        ratio_M_kink_to_M_foot=M_kink_e / M_foot,
        ratio_m_NG_to_m_tau=m_NG_tau / M_TAU_MEV,
    )
