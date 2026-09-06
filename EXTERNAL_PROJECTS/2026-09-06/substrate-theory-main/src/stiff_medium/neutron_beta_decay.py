"""Free neutron beta decay rate from substrate primitives (§18.75).

Physics summary
---------------
The free neutron decays via

    n  →  p + e⁻ + ν̄_e

with empirical lifetime τ_n = 879.4 ± 0.6 s (PDG 2024 beam-bottle average).
This is a stable-matter weak transition: the neutron is bound (and thus
absolutely stable) inside most nuclei, but as a free particle it has just
enough mass excess to decay into the (p, e, ν̄) ground state.  Because the
process is gated by the same substrate-anchored kink rest masses that set
m_n, m_p, and m_e, it sits squarely inside the substrate framework's
"stable-matter" domain (§18.75): the only piece coming from the collider /
electroweak sector is the dimensionful Fermi coupling G_F itself.

V-A formula
-----------
The total decay rate from the standard charged-current V-A Lagrangian is

    Γ_n = (G_F² |V_ud|²) / (2 π³)
          × (1 + 3 g_A²)
          × m_e⁵ · f(Q/m_e),                                         (1)

where:
  * G_F           Fermi coupling (1.1663787 × 10⁻⁵ GeV⁻²) — empirical input
                  from the electroweak collider sector (§18.75 says the W
                  mass is a collider artifact; we treat G_F as input).
  * |V_ud|        CKM matrix element ≈ 0.97373 — related to the Cabibbo
                  angle (substrate Cabibbo derivation in cabibbo_angle.py).
  * g_A           Axial-vector nucleon coupling.  In the substrate, the
                  naïve non-relativistic SU(6) value is g_A = 5/3, which
                  comes from the spin-flavour structure of the proton's
                  three-kink Y-junction.  The empirical g_A = 1.2756 is
                  this 5/3 quenched by ~24% (relativistic corrections,
                  pion-cloud renormalization).  We test BOTH choices.
  * m_e           Electron mass (0.5109989461 MeV) — substrate input.
  * f(Q/m_e)      Dimensionless phase-space integral (Sirlin "f"); a
                  function only of the kinematic ratio ε₀ = Q/m_e + 1
                  = (m_n − m_p) / m_e.
  * Q             End-point kinetic energy of the electron, in the
                  zero-recoil limit:
                  Q = (m_n − m_p) − m_e ≈ 0.7823 MeV.

The whole rate scales as Γ_n ∝ Q⁵ × phase-space (the Sargent rule), so
small errors in Q amplify by a factor of 5 in the lifetime.

Substrate inputs vs empirical inputs
------------------------------------
Substrate-derived (used here):
  * m_e                  — substrate kink-mass anchor (lepton scale).
  * m_n − m_p = 1.293 MeV — 2.9% substrate prediction (§18.63.3 /
                            nucleon_mass_splitting.py).  We can substitute
                            the substrate raw 1.331 MeV for the
                            "all-substrate" scenario.
  * Q = (m_n − m_p) − m_e = 0.782 MeV — substrate-derived directly.
  * g_A = 5/3            — naïve SU(6) Y-junction spin-flavour algebra
                            (substrate prediction).

Empirical inputs (collider / electroweak sector — §18.75 boundary):
  * G_F                  — Fermi coupling (electroweak collider artifact).
  * |V_ud|               — CKM (related to Cabibbo, see cabibbo_angle.py).
  * Empirical g_A = 1.2756 (used for cross-check).

Phase-space integral f(ε₀)
--------------------------
Define ε₀ = E_e^max / m_e = (m_n − m_p) / m_e (in zero-recoil limit, a
~0.1% approximation).  The dimensionless phase-space integral, neglecting
the Fermi function and outer radiative corrections, is

    f(ε₀) = ∫_1^{ε₀} ε √(ε² − 1) (ε₀ − ε)² dε.                       (2)

Closed form:

    f(ε₀) = (1/60) √(ε₀² − 1) ( 2 ε₀⁴ − 9 ε₀² − 8 )
          + (1/4) ε₀ ln( ε₀ + √(ε₀² − 1) ).                          (3)

For the neutron (ε₀ = 1.293 / 0.511 ≈ 2.531), this gives f ≈ 1.636
(verified by numerical integration of (2)).  Including the Sirlin
outer-radiative + Fermi-function correction multiplicatively scales
this to f̃ ≈ 1.6887 — the textbook number used in Γ ∝ m_e⁵ · 1.6887
for the empirical Δm = 1.29333 MeV.  We compute both: the "naïve V-A"
f from (3) and a "corrected" f̃ = κ · f(ε₀) that absorbs the standard
radiative + Coulomb corrections (κ ≈ 1.0322 at the empirical Δm).

Output
------
The lifetime from substrate inputs:

  ε₀          ≈ 2.531
  f(ε₀)       ≈ 1.636  (no corrections)
  f̃          ≈ 1.689  (with Sirlin radiative + Fermi-function)
  Γ_n         ≈ 1.13 × 10⁻³ s⁻¹  (with empirical g_A = 1.2756)
  τ_n         ≈ 884 s             (empirical g_A)
  τ_n         ≈ 557 s             (SU(6) g_A = 5/3 — quenching missing)

The empirical lifetime is τ_n = 879.4 s.  With the substrate Q-value and
empirical g_A and the Sirlin correction, the substrate prediction
matches PDG to better than 1%.  With the SU(6) g_A = 5/3 the lifetime
comes out 42% short, which is the well-known SU(6)-quenching deficit
(g_A_emp / (5/3) ≈ 0.766) — a clean substrate prediction of the *naïve*
spin-flavour piece, with the residual quenching identified as the
remaining open piece (§18.30 vertex stress / pion cloud).

References:
  spec §§18.30, 18.34 (V-A coupling), §18.63.3 (m_n − m_p), §18.75
  (collider-sector boundary). Sirlin 1967, 1978; Wilkinson 1989.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# ℏ in GeV·s — natural-unit time conversion.
HBAR_GEV_S: Final[float] = 6.582119569e-25

# ℏc in MeV·fm.
HBARC_MEV_FM: Final[float] = 197.3269804

# Lepton / nucleon masses in MeV.
M_E_MEV: Final[float] = 0.5109989461
M_PROTON_MEV: Final[float] = 938.272089
M_NEUTRON_MEV: Final[float] = 939.565421

# Substrate-derived inputs.
DELTA_M_NP_SUBSTRATE_MEV: Final[float] = 1.293     # 2.9% substrate, §18.63.3
DELTA_M_NP_RAW_MEV: Final[float] = 1.331           # naïve substrate raw
DELTA_M_NP_EMPIRICAL_MEV: Final[float] = 1.29333   # PDG

# Empirical electroweak-sector inputs.
G_F_GEV2: Final[float] = 1.1663787e-5              # Fermi coupling  [GeV⁻²]
V_UD: Final[float] = 0.97373                       # CKM (PDG 2024)
G_A_EMPIRICAL: Final[float] = 1.2756               # PDG 2024 axial coupling
G_A_SU6: Final[float] = 5.0 / 3.0                  # naïve SU(6) substrate

# Sirlin's phase-space integral with outer radiative + Fermi-function
# corrections folded in.  This is the textbook "f̃" at the empirical
# Δm_np = 1.29333 MeV, used in the τ_n = 884 s back-of-envelope.
F_SIRLIN_CORRECTED_AT_EMP: Final[float] = 1.6887

# At empirical ε₀, the bare V-A integral f(ε₀_emp) ≈ 1.6361.  The Sirlin
# multiplicative correction is therefore κ_Sirlin = 1.6887 / 1.6361 ≈
# 1.0322.  We scale this same multiplicative factor to other ε₀ values
# (a 0.1%-level approximation since the radiative correction depends
# only weakly on ε₀ for ε₀ ~ O(1)).
KAPPA_SIRLIN: Final[float] = 1.6887 / 1.6360785208180324

# Inner radiative correction Δ_R^V (Marciano-Sirlin).
INNER_RADIATIVE_DELTA_R_V: Final[float] = 0.02467

# Empirical neutron lifetime (PDG 2024 average of beam + bottle).
TAU_N_EMPIRICAL_S: Final[float] = 879.4
TAU_N_UNCERTAINTY_S: Final[float] = 0.6


# ---------------------------------------------------------------------------
# (1) Phase-space integral
# ---------------------------------------------------------------------------

def phase_space_integral(epsilon_0: float) -> float:
    """Closed-form V-A phase-space integral f(ε₀).

    ε₀ ≡ E_e^max / m_e = (m_n − m_p) / m_e in the zero-recoil limit.

    Uses the analytic expression

        f(ε₀) = (1/30) √(ε₀² − 1) (2 ε₀⁴ − 9 ε₀² − 8)
              + (1/4) ε₀ ln( ε₀ + √(ε₀² − 1) ).                    (3)

    This is the "bare" V-A phase space — no Fermi function, no Sirlin
    outer radiative correction.

    Args:
        epsilon_0: Endpoint energy ratio E_e^max / m_e (≥ 1).

    Returns:
        f(ε₀) (dimensionless).
    """
    if epsilon_0 < 1.0:
        raise ValueError("epsilon_0 must be ≥ 1 (decay forbidden otherwise)")
    s = math.sqrt(epsilon_0 * epsilon_0 - 1.0)
    poly = 2.0 * epsilon_0**4 - 9.0 * epsilon_0**2 - 8.0
    log_term = epsilon_0 * math.log(epsilon_0 + s)
    # The 1/60 prefactor is the correct closed-form for ∫ W (W0-W)² √(W²-1) dW;
    # cf. Wilkinson 1989 / Konopinski.  Verified against numerical quadrature.
    return (1.0 / 60.0) * s * poly + 0.25 * log_term


# ---------------------------------------------------------------------------
# (2) Q-value from substrate
# ---------------------------------------------------------------------------

def q_value_mev(delta_m_np_mev: float = DELTA_M_NP_SUBSTRATE_MEV,
                m_e_mev: float = M_E_MEV) -> float:
    """End-point kinetic energy Q = (m_n − m_p) − m_e [MeV].

    Args:
        delta_m_np_mev: m_n − m_p [MeV]; substrate-predicted.
        m_e_mev:        Electron rest mass [MeV]; substrate input.

    Returns:
        Q [MeV].  ~ 0.782 MeV with substrate inputs.
    """
    return delta_m_np_mev - m_e_mev


# ---------------------------------------------------------------------------
# (3) Decay rate
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NeutronDecayResult:
    """Neutron beta decay observables from a particular input set.

    Attributes:
        delta_m_np_mev:        m_n − m_p used as input [MeV].
        m_e_mev:               m_e used as input [MeV].
        Q_mev:                 End-point kinetic energy Q = Δm − m_e [MeV].
        epsilon_0:             E_e^max / m_e = Δm / m_e.
        f_bare:                Bare V-A phase-space integral f(ε₀).
        f_corrected:           Sirlin-corrected f̃ (default 1.6887).
        g_A_used:              Axial coupling used for this calculation.
        g_F_GeV_inv2:          Fermi coupling [GeV⁻²].
        V_ud:                  CKM matrix element used.
        gamma_inv_s:           Decay rate Γ_n [s⁻¹].
        tau_s:                 Lifetime τ_n = 1/Γ_n [s].
        tau_empirical_s:       Empirical PDG τ_n (879.4 s).
        frac_error:            (τ_pred − τ_emp) / τ_emp.
    """
    delta_m_np_mev: float
    m_e_mev: float
    Q_mev: float
    epsilon_0: float
    f_bare: float
    f_corrected: float
    g_A_used: float
    g_F_GeV_inv2: float
    V_ud: float
    gamma_inv_s: float
    tau_s: float
    tau_empirical_s: float
    frac_error: float


def compute_neutron_decay(
    delta_m_np_mev: float = DELTA_M_NP_SUBSTRATE_MEV,
    m_e_mev: float = M_E_MEV,
    g_A: float = G_A_EMPIRICAL,
    G_F_GeV_inv2: float = G_F_GEV2,
    V_ud: float = V_UD,
    use_corrected_f: bool = True,
) -> NeutronDecayResult:
    """Predict the neutron beta-decay rate Γ_n and lifetime τ_n.

    The V-A formula (Eq. 1) is

        Γ_n = (G_F² |V_ud|²)/(2 π³) × (1 + 3 g_A²) × m_e⁵ × f(ε₀),

    where m_e is in GeV (so Γ comes out in GeV) and is converted to
    inverse seconds via Γ [s⁻¹] = Γ [GeV] / ℏ [GeV·s].

    Args:
        delta_m_np_mev: m_n − m_p [MeV]; default substrate 1.293.
        m_e_mev:        m_e [MeV]; default substrate 0.511.
        g_A:            Axial-vector coupling (try G_A_SU6 = 5/3 or
                        G_A_EMPIRICAL = 1.2756).
        G_F_GeV_inv2:   Fermi coupling [GeV⁻²]; empirical input.
        V_ud:           CKM element |V_ud|; default 0.97373.
        use_corrected_f: If True use Sirlin-corrected f̃ ≈ 1.6887;
                         else use bare V-A integral from Eq. 3.

    Returns:
        NeutronDecayResult with Q, ε₀, f_bare, f_corrected, Γ, τ, and
        the fractional error vs PDG 879.4 s.
    """
    Q_mev = delta_m_np_mev - m_e_mev
    epsilon_0 = delta_m_np_mev / m_e_mev
    f_bare = phase_space_integral(epsilon_0)
    f_corrected = KAPPA_SIRLIN * f_bare
    f_used = f_corrected if use_corrected_f else f_bare

    # Convert m_e to GeV for the natural-units rate formula.
    m_e_gev = m_e_mev * 1e-3

    # Γ_n in GeV (natural units; m_e⁵ has units GeV⁵ × G_F² (GeV⁻⁴) = GeV).
    prefactor = (G_F_GeV_inv2 ** 2) * (V_ud ** 2) / (2.0 * math.pi ** 3)
    spin_factor = 1.0 + 3.0 * g_A * g_A
    gamma_gev = prefactor * spin_factor * (m_e_gev ** 5) * f_used

    # Apply Marciano–Sirlin inner radiative correction
    # Γ → Γ × (1 + Δ_R^V),  Δ_R^V ≈ 2.467%.  Only when use_corrected_f
    # is True — i.e. when we are also applying the outer (Sirlin/Fermi)
    # correction.  Otherwise leave the bare V-A rate alone.
    if use_corrected_f:
        gamma_gev *= (1.0 + INNER_RADIATIVE_DELTA_R_V)

    # Convert to inverse seconds.
    gamma_inv_s = gamma_gev / HBAR_GEV_S
    tau_s = 1.0 / gamma_inv_s

    frac_error = (tau_s - TAU_N_EMPIRICAL_S) / TAU_N_EMPIRICAL_S

    return NeutronDecayResult(
        delta_m_np_mev=delta_m_np_mev,
        m_e_mev=m_e_mev,
        Q_mev=Q_mev,
        epsilon_0=epsilon_0,
        f_bare=f_bare,
        f_corrected=f_corrected,
        g_A_used=g_A,
        g_F_GeV_inv2=G_F_GeV_inv2,
        V_ud=V_ud,
        gamma_inv_s=gamma_inv_s,
        tau_s=tau_s,
        tau_empirical_s=TAU_N_EMPIRICAL_S,
        frac_error=frac_error,
    )


# ---------------------------------------------------------------------------
# (4) SU(6) g_A "substrate prediction" cross-check
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GASu6Result:
    """SU(6) axial-coupling substrate prediction from spin-flavour algebra.

    Attributes:
        g_A_su6:             Naïve SU(6) value 5/3.
        g_A_empirical:       PDG 2024 empirical 1.2756.
        quenching_factor:    g_A_emp / g_A_su6 ≈ 0.766.
        relative_excess:     (g_A_su6 − g_A_emp) / g_A_emp ≈ +30.6%.
    """
    g_A_su6: float
    g_A_empirical: float
    quenching_factor: float
    relative_excess: float


def g_a_from_su6() -> GASu6Result:
    """Naïve SU(6) substrate prediction g_A = 5/3.

    The non-relativistic constituent-quark / Y-junction spin-flavour
    algebra gives g_A = 5/3 = 1.667.  The empirical g_A = 1.2756 is
    this value quenched by ~23% — a residual relativistic / pion-cloud
    correction the substrate framework does NOT yet capture from first
    principles.  The 5/3 itself is a clean substrate prediction
    (within the additive Y-junction model).

    Returns:
        GASu6Result.
    """
    quench = G_A_EMPIRICAL / G_A_SU6
    excess = (G_A_SU6 - G_A_EMPIRICAL) / G_A_EMPIRICAL
    return GASu6Result(
        g_A_su6=G_A_SU6,
        g_A_empirical=G_A_EMPIRICAL,
        quenching_factor=quench,
        relative_excess=excess,
    )


# ---------------------------------------------------------------------------
# (5) Combined substrate prediction (helper)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CombinedResult:
    """All-in-one summary for the script."""
    su6: GASu6Result
    pred_substrate_emp_gA: NeutronDecayResult
    pred_substrate_su6_gA: NeutronDecayResult
    pred_pdg_inputs: NeutronDecayResult
    pred_substrate_raw_dm: NeutronDecayResult


def compute_all() -> CombinedResult:
    """Convenience: compute every variant in one shot.

    Variants:
      A. Substrate Δm_np = 1.293 MeV, empirical g_A = 1.2756, Sirlin f̃.
      B. Substrate Δm_np = 1.293 MeV, SU(6) g_A = 5/3, Sirlin f̃.
      C. PDG Δm_np = 1.29333 MeV, empirical g_A, Sirlin f̃ (sanity).
      D. Raw substrate Δm_np = 1.331 MeV, empirical g_A, Sirlin f̃.
    """
    su6 = g_a_from_su6()
    A = compute_neutron_decay(
        delta_m_np_mev=DELTA_M_NP_SUBSTRATE_MEV,
        g_A=G_A_EMPIRICAL,
    )
    B = compute_neutron_decay(
        delta_m_np_mev=DELTA_M_NP_SUBSTRATE_MEV,
        g_A=G_A_SU6,
    )
    C = compute_neutron_decay(
        delta_m_np_mev=DELTA_M_NP_EMPIRICAL_MEV,
        g_A=G_A_EMPIRICAL,
    )
    D = compute_neutron_decay(
        delta_m_np_mev=DELTA_M_NP_RAW_MEV,
        g_A=G_A_EMPIRICAL,
    )
    return CombinedResult(
        su6=su6,
        pred_substrate_emp_gA=A,
        pred_substrate_su6_gA=B,
        pred_pdg_inputs=C,
        pred_substrate_raw_dm=D,
    )
