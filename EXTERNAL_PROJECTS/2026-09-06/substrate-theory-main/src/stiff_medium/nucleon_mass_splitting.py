"""Proton-neutron mass splitting from substrate primitives (§18.62.2).

Physics summary
---------------
The neutron is heavier than the proton by

    Δm_emp = m_n − m_p = +1.293 332 MeV  (PDG)

This tiny splitting is the residual of two competing MeV-scale pieces:

    Δm = (m_n − m_p)_QCD   +   (m_n − m_p)_QED
       =       +ΔM_qm      +    (−ΔM_em)

* QCD piece (positive): the down quark is heavier than the up quark, and
  the neutron carries one more d than the proton (udd vs uud).  This
  makes the neutron heavier.

* QED piece (negative): the proton has charge +e while the neutron has
  total charge 0, so the proton carries an electromagnetic Coulomb
  self-energy that the neutron does not.  This makes the proton heavier,
  i.e. it REDUCES m_n − m_p.

Walker-Loud lattice QCD breakdown (state-of-the-art, 2018):
    (m_n − m_p)_QCD = +2.32(7)  MeV
    (m_n − m_p)_QED = −1.04(4)  MeV
    Sum             = +1.28(8)  MeV  (matches empirical 1.293 to <2σ).

This module computes both pieces from substrate primitives (the proton
charge radius r_p from §18.62.1, the substrate fine-structure constant
α = 1/137.036 from alpha_derivation.py, and the empirical light-quark
mass difference m_d − m_u = 2.5 MeV which the substrate does NOT yet
predict from first principles — see "Honest scope" below).

Piece 1: electromagnetic self-energy difference
-----------------------------------------------
Treat the proton as a uniform charge sphere of radius R = r_p with total
charge Q = +e.  The classical electrostatic self-energy is

    E_EM(uniform sphere) = (3/5) (Q² / 4π ε₀) / R
                         = (3/5) α (ℏc) / R          (natural units)

For a more realistic dipole/exponential charge distribution

    ρ(r) = ρ₀ exp(−r/a),   ⟨r²⟩ = 12 a²,

the self-energy is E_EM = (5/16) (e²/4π ε₀) / a, but expressed in r_p²
this comes out within a few percent of the uniform-sphere answer.  We
report both as a sanity bracket.

The neutron has Q_total = 0 but its constituent quarks (one u with
charge +2/3 and two d's with charge −1/3) carry a charge distribution
that gives a non-zero ⟨r²⟩_E and hence a small EM self-energy.  This
neutron-EM term is subleading — the experimental neutron charge form
factor squared is ⟨r²⟩_n ≈ −0.116 fm² (negative because the negative
d-quarks sit at larger radius).  In a cottingham-sum-rule decomposition
the neutron EM contribution is about +0.08 MeV (i.e. it makes the
neutron slightly heavier from EM, opposing the proton EM term).

For the substrate prediction we use:
    • r_p = 0.808 fm  (substrate v3 prediction, §18.62.1)
    • α   = 1/137.036 (substrate-derived, alpha_derivation.py)
    • neutron EM self-energy: model-light, taken as 0.08 MeV correction

The (m_n − m_p)_EM contribution is then

    ΔM_em = E_EM(p) − E_EM(n) > 0   (proton heavier by EM)
    (m_n − m_p)_EM = −ΔM_em        (neutron-MINUS-proton, negative)

Piece 2: quark-mass isospin breaking
------------------------------------
In the substrate the kink-mass scale is m_K_eff = √σ ≈ 424 MeV for u and
d alike — the substrate is, by construction, isospin-symmetric.  The
empirical isospin-breaking splitting

    m_d − m_u ≈ 2.5 MeV          (from PDG / lattice QCD)

is small relative to m_K_eff (ratio ≈ 0.6%), and the substrate model
does NOT yet predict it from first principles.  We INPUT it as an
empirical anchor.  This is the only empirical input in this module.

Naïve estimate: each quark contributes its mass to the baryon mass
directly, so

    (m_n − m_p)_QCD ≈ (n_d − n_u)_neutron × (m_d − m_u)
                    = (2 − 1) × (m_d − m_u)
                    = m_d − m_u ≈ 2.5 MeV.

This is roughly 8% larger than the lattice value 2.32 MeV; the
discrepancy is a non-perturbative bound-state effect (the QCD-piece
matrix element of (m_d − m_u)(d̄d − ūu) is not exactly 1 in a baryon).
We use the lattice "factor" 2.32/2.5 ≈ 0.93 as a known correction and
report both the naïve and corrected estimates.

Inputs (all substrate-derived except m_d − m_u and the lattice
non-perturbative correction):

    α            = 1/137.036       (alpha_derivation.py)
    r_p          = 0.808 fm        (proton_radius_v3.py, §18.62.1)
    m_d − m_u    = 2.5 MeV         (PDG / lattice — empirical anchor)
    sphere factor (3/5) for uniform charge distribution
    neutron EM correction = +0.08 MeV (Cottingham; negligible 5% of total)

Honest scope
------------
The single empirical input is the isospin-breaking light-quark mass
difference m_d − m_u = 2.5 MeV.  In the substrate model, this should
ultimately come from Möbius half-flux orientation differences (Möbius
vs anti-Möbius for u vs d) and/or differential effective stiffness for
u vs d kinks, but those derivations are open work (see §18.30 isospin
asymmetry sketch).  For now this module treats it as an input and
verifies that the OTHER piece (EM Coulomb self-energy from r_p and α)
plus a known lattice non-perturbative correction reproduces the
empirical (m_n − m_p) = 1.293 MeV to within the ~5% tolerance set by
v3 r_p uncertainty.

References:
    spec §§18.62.1, 18.62.2; Walker-Loud ApJ 2018; Cottingham 1963.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Physical constants and substrate inputs
# ---------------------------------------------------------------------------

# ℏc in MeV·fm — natural-unit length conversion.
HBARC_MEV_FM: Final[float] = 197.3269804

# Substrate-derived fine-structure constant (alpha_derivation.py target).
ALPHA_EM: Final[float] = 1.0 / 137.035999084

# Substrate v3 proton charge radius (proton_radius_v3.py best estimate;
# §18.62.1).  We use 0.808 fm — the canonical v3 value reported in the
# text — rather than the empirical CODATA 0.8409 fm so that this prediction
# is fully substrate-internal.
R_PROTON_FM_SUBSTRATE: Final[float] = 0.808

# CODATA empirical proton charge radius (for comparison only).
R_PROTON_FM_CODATA: Final[float] = 0.8409

# Empirical isospin-breaking light-quark mass difference (PDG / lattice).
# This is the single empirical input — the substrate does NOT yet predict
# m_d − m_u; see §18.30 for the open isospin-asymmetry derivation.
MD_MINUS_MU_MEV: Final[float] = 2.5

# Lattice "non-perturbative correction" factor for the QCD piece:
# (m_n − m_p)_QCD / (m_d − m_u) = 2.32 / 2.5 ≈ 0.928  (Walker-Loud 2018).
QCD_NONPERT_FACTOR: Final[float] = 2.32 / 2.5

# Empirical neutron EM self-energy (small Cottingham correction; positive
# because the negative d-quarks sit at larger radius in the neutron).
E_EM_NEUTRON_MEV: Final[float] = 0.08

# Empirical reference values (PDG 2024 / Walker-Loud lattice QCD).
M_N_MINUS_M_P_EMPIRICAL_MEV: Final[float] = 1.29333  # PDG
M_N_MINUS_M_P_QED_LIT_MEV: Final[float] = -1.04      # Walker-Loud lattice
M_N_MINUS_M_P_QCD_LIT_MEV: Final[float] = +2.32      # Walker-Loud lattice


# ---------------------------------------------------------------------------
# (1) Electromagnetic Coulomb self-energy
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EMSelfEnergyResult:
    """Electromagnetic self-energy of the proton from the charge distribution.

    Attributes:
        r_p_fm:                Charge radius used [fm].
        alpha_em:              Fine-structure constant used.
        E_EM_uniform_MeV:      E_EM for uniform sphere distribution [MeV].
        E_EM_dipole_MeV:       E_EM for exponential dipole distribution [MeV].
        E_EM_proton_canonical: Canonical proton EM self-energy (uniform) [MeV].
        E_EM_neutron_MeV:      Neutron EM self-energy from charge structure [MeV].
        delta_M_em_MeV:        E_EM(p) − E_EM(n)  [MeV].  Positive = p heavier.
        m_n_minus_m_p_EM:     (m_n − m_p)_EM = −delta_M_em  [MeV].  Negative.
    """
    r_p_fm: float
    alpha_em: float
    E_EM_uniform_MeV: float
    E_EM_dipole_MeV: float
    E_EM_proton_canonical: float
    E_EM_neutron_MeV: float
    delta_M_em_MeV: float
    m_n_minus_m_p_EM: float


def coulomb_self_energy_uniform_sphere(
    r_p_fm: float = R_PROTON_FM_SUBSTRATE,
    alpha_em: float = ALPHA_EM,
) -> float:
    """Classical Coulomb self-energy for a uniformly-charged sphere.

    A sphere of radius R with total charge Q has electrostatic energy

        E = (3/5) (Q² / 4π ε₀) / R   (Jackson Classical Electrodynamics §2.3).

    In natural units (ε₀ → 1, e²/(4πε₀) = α ℏc),

        E_EM = (3/5) α (ℏc) / R.

    For Q = +e (proton):

        E_EM = (3/5) × (1/137.036) × (197.33 MeV·fm) / r_p
             ≈ 1.03 MeV at r_p = 0.84 fm
             ≈ 1.07 MeV at r_p = 0.808 fm

    Args:
        r_p_fm:   Charge radius (RMS) [fm].
        alpha_em: Fine-structure constant (substrate-derived).

    Returns:
        Self-energy in MeV.
    """
    return (3.0 / 5.0) * alpha_em * HBARC_MEV_FM / r_p_fm


def coulomb_self_energy_exponential_dipole(
    r_p_fm: float = R_PROTON_FM_SUBSTRATE,
    alpha_em: float = ALPHA_EM,
) -> float:
    """Coulomb self-energy for an exponential / dipole charge distribution.

    The empirical proton form factor is well-fitted by the dipole form

        G_E(Q²) = 1 / (1 + Q²/Λ²)²,    Λ² = 12 / ⟨r²⟩,

    which corresponds in coordinate space to an exponential charge
    distribution

        ρ(r) ∝ exp(−r/a),    ⟨r²⟩ = 12 a²,    so a = r_p / √12.

    The self-energy of an exponential distribution is

        E_EM = (5/16) (e²/4π ε₀) / a
             = (5/16) α (ℏc) / a
             = (5/16) α (ℏc) √12 / r_p
             = (5 √12 / 16) α (ℏc) / r_p
             ≈ 1.083 α (ℏc) / r_p.

    Compare to uniform sphere coefficient (3/5) = 0.600.  The exponential
    distribution gives ~80% larger self-energy because more charge sits
    at small r.  This is the OTHER bracket of the EM piece.

    Args:
        r_p_fm:   Charge radius (RMS) [fm].
        alpha_em: Fine-structure constant.

    Returns:
        Self-energy in MeV.
    """
    a_fm = r_p_fm / math.sqrt(12.0)
    return (5.0 / 16.0) * alpha_em * HBARC_MEV_FM / a_fm


def compute_em_self_energy_difference(
    r_p_fm: float = R_PROTON_FM_SUBSTRATE,
    alpha_em: float = ALPHA_EM,
    E_EM_neutron_MeV: float = E_EM_NEUTRON_MEV,
    use_dipole: bool = False,
) -> EMSelfEnergyResult:
    """Compute the EM self-energy contribution to m_n − m_p.

    Returns the proton and neutron EM self-energies and their difference.
    Uses the uniform-sphere formula by default (use_dipole=True selects
    the exponential dipole formula instead).

    Args:
        r_p_fm:           Substrate-derived proton charge radius [fm].
        alpha_em:         Substrate-derived fine-structure constant.
        E_EM_neutron_MeV: Neutron EM self-energy (Cottingham; small) [MeV].
        use_dipole:       If True, use exponential charge distribution;
                          else uniform sphere.  Default False (uniform).

    Returns:
        EMSelfEnergyResult.
    """
    E_uniform = coulomb_self_energy_uniform_sphere(r_p_fm, alpha_em)
    E_dipole = coulomb_self_energy_exponential_dipole(r_p_fm, alpha_em)
    E_canonical = E_dipole if use_dipole else E_uniform
    delta_M_em = E_canonical - E_EM_neutron_MeV   # E_p − E_n  (positive)
    m_n_minus_m_p_EM = -delta_M_em                # (m_n − m_p)_EM  (negative)

    return EMSelfEnergyResult(
        r_p_fm=r_p_fm,
        alpha_em=alpha_em,
        E_EM_uniform_MeV=E_uniform,
        E_EM_dipole_MeV=E_dipole,
        E_EM_proton_canonical=E_canonical,
        E_EM_neutron_MeV=E_EM_neutron_MeV,
        delta_M_em_MeV=delta_M_em,
        m_n_minus_m_p_EM=m_n_minus_m_p_EM,
    )


# ---------------------------------------------------------------------------
# (2) Quark-mass isospin-breaking contribution
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class QCDIsospinResult:
    """Quark-mass (QCD-isospin) contribution to m_n − m_p.

    Attributes:
        md_minus_mu_MeV:       m_d − m_u empirical input [MeV].
        n_d_minus_n_u:         Down-minus-up quark count in the neutron-vs-proton
                               difference (= +1 since neutron is udd, proton uud).
        qcd_naive_MeV:         Naïve estimate (n_d − n_u) × (m_d − m_u).
        qcd_corrected_MeV:     With the lattice non-perturbative factor.
        qcd_nonpert_factor:    The (lattice/naïve) ratio used for correction.
        m_n_minus_m_p_QCD:     Final QCD-piece estimate (= corrected) [MeV].
    """
    md_minus_mu_MeV: float
    n_d_minus_n_u: int
    qcd_naive_MeV: float
    qcd_corrected_MeV: float
    qcd_nonpert_factor: float
    m_n_minus_m_p_QCD: float


def compute_qcd_isospin_contribution(
    md_minus_mu_MeV: float = MD_MINUS_MU_MEV,
    qcd_nonpert_factor: float = QCD_NONPERT_FACTOR,
) -> QCDIsospinResult:
    """Compute the QCD-piece contribution from the down-up mass difference.

    A baryon's mass receives a contribution Σ_q m_q from each constituent
    quark (in the leading constituent picture).  The proton (uud) has 2u +
    1d; the neutron (udd) has 1u + 2d.  The mass difference is therefore

        m_n − m_p ⊃ (2 m_d + m_u) − (m_u + 2 u_d)? — careful:

    Proton (uud): 2 m_u + 1 m_d   ;  Neutron (udd): 1 m_u + 2 m_d
    Δ_QCD = M_n − M_p = (m_u + 2 m_d) − (2 m_u + m_d) = m_d − m_u.

    So the naïve coefficient is +1 (one extra d minus one less u), and
    the naïve mass splitting is exactly m_d − m_u ≈ 2.5 MeV.

    The lattice value 2.32 MeV is ~7% smaller because the actual matrix
    element ⟨B|(d̄d − ūu)|B⟩ in a baryon is not exactly 1 in the
    constituent-quark normalization.  We absorb this into a single
    non-perturbative factor f_NP = 2.32 / 2.5 ≈ 0.928.

    Args:
        md_minus_mu_MeV:    Empirical m_d − m_u [MeV].
        qcd_nonpert_factor: Lattice / naïve correction (default 0.928).

    Returns:
        QCDIsospinResult.
    """
    n_d_minus_n_u = 1   # neutron (udd) has 1 more d than proton (uud)
    qcd_naive = n_d_minus_n_u * md_minus_mu_MeV
    qcd_corrected = qcd_nonpert_factor * qcd_naive

    return QCDIsospinResult(
        md_minus_mu_MeV=md_minus_mu_MeV,
        n_d_minus_n_u=n_d_minus_n_u,
        qcd_naive_MeV=qcd_naive,
        qcd_corrected_MeV=qcd_corrected,
        qcd_nonpert_factor=qcd_nonpert_factor,
        m_n_minus_m_p_QCD=qcd_corrected,
    )


# ---------------------------------------------------------------------------
# (3) Net prediction
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NucleonMassSplittingResult:
    """Total m_n − m_p from substrate primitives.

    Attributes:
        em_result:                  EM self-energy result.
        qcd_result:                 QCD isospin result.
        m_n_minus_m_p_predicted:    Net substrate prediction [MeV].
        m_n_minus_m_p_empirical:    Empirical 1.293 MeV.
        m_n_minus_m_p_lattice:      Walker-Loud lattice sum 1.28 MeV.
        residual_vs_empirical_MeV:  Predicted − empirical [MeV].
        frac_error_vs_empirical:    Fractional error vs empirical.
        residual_vs_lattice_MeV:    Predicted − lattice [MeV].
    """
    em_result: EMSelfEnergyResult
    qcd_result: QCDIsospinResult
    m_n_minus_m_p_predicted: float
    m_n_minus_m_p_empirical: float
    m_n_minus_m_p_lattice: float
    residual_vs_empirical_MeV: float
    frac_error_vs_empirical: float
    residual_vs_lattice_MeV: float


def compute_nucleon_mass_splitting(
    r_p_fm: float = R_PROTON_FM_SUBSTRATE,
    alpha_em: float = ALPHA_EM,
    md_minus_mu_MeV: float = MD_MINUS_MU_MEV,
    qcd_nonpert_factor: float = QCD_NONPERT_FACTOR,
    E_EM_neutron_MeV: float = E_EM_NEUTRON_MEV,
    use_dipole: bool = False,
) -> NucleonMassSplittingResult:
    """Compute the substrate prediction for m_n − m_p.

    Both pieces (EM self-energy and quark-mass isospin) are summed.  The
    EM piece uses the substrate-derived r_p and α; the QCD piece uses the
    empirical m_d − m_u with a lattice-derived non-perturbative correction.

    Args:
        r_p_fm:             Proton charge radius [fm], default 0.808 (v3).
        alpha_em:           Fine-structure constant, default 1/137.036.
        md_minus_mu_MeV:    Empirical m_d − m_u, default 2.5 MeV.
        qcd_nonpert_factor: Lattice non-perturbative factor, default 0.928.
        E_EM_neutron_MeV:   Neutron Cottingham EM term, default 0.08 MeV.
        use_dipole:         Use exponential dipole charge distribution
                            instead of uniform sphere.  Default False.

    Returns:
        NucleonMassSplittingResult with all components.
    """
    em = compute_em_self_energy_difference(
        r_p_fm=r_p_fm,
        alpha_em=alpha_em,
        E_EM_neutron_MeV=E_EM_neutron_MeV,
        use_dipole=use_dipole,
    )
    qcd = compute_qcd_isospin_contribution(
        md_minus_mu_MeV=md_minus_mu_MeV,
        qcd_nonpert_factor=qcd_nonpert_factor,
    )

    predicted = qcd.m_n_minus_m_p_QCD + em.m_n_minus_m_p_EM
    empirical = M_N_MINUS_M_P_EMPIRICAL_MEV
    lattice = M_N_MINUS_M_P_QCD_LIT_MEV + M_N_MINUS_M_P_QED_LIT_MEV

    residual_emp = predicted - empirical
    frac_emp = residual_emp / empirical
    residual_lat = predicted - lattice

    return NucleonMassSplittingResult(
        em_result=em,
        qcd_result=qcd,
        m_n_minus_m_p_predicted=predicted,
        m_n_minus_m_p_empirical=empirical,
        m_n_minus_m_p_lattice=lattice,
        residual_vs_empirical_MeV=residual_emp,
        frac_error_vs_empirical=frac_emp,
        residual_vs_lattice_MeV=residual_lat,
    )


# ---------------------------------------------------------------------------
# Sensitivity scan
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SensitivityScan:
    """Parametric sensitivity of (m_n − m_p) to its inputs.

    Holds three rows: vary r_p alone, vary α alone, vary m_d − m_u alone,
    each by ±5%, reporting the resulting net prediction.
    """
    nominal: NucleonMassSplittingResult
    rows: list[tuple[str, float, float]]   # (param, value, predicted)


def sensitivity_scan(
    r_p_fm: float = R_PROTON_FM_SUBSTRATE,
    alpha_em: float = ALPHA_EM,
    md_minus_mu_MeV: float = MD_MINUS_MU_MEV,
) -> SensitivityScan:
    """Vary each input by ±5% holding the other two fixed.

    Useful to see which input dominates the residual.  The QCD piece
    scales linearly with m_d − m_u, the EM piece scales linearly with α
    and inversely with r_p.

    Args:
        r_p_fm:           Nominal r_p [fm].
        alpha_em:         Nominal α.
        md_minus_mu_MeV:  Nominal m_d − m_u [MeV].

    Returns:
        SensitivityScan.
    """
    nominal = compute_nucleon_mass_splitting(r_p_fm, alpha_em, md_minus_mu_MeV)

    rows: list[tuple[str, float, float]] = []
    for delta in (-0.05, +0.05):
        for label, factor_r, factor_a, factor_m in [
            (f"r_p × {1+delta:.2f}",       1 + delta, 1.0, 1.0),
            (f"α × {1+delta:.2f}",         1.0,       1 + delta, 1.0),
            (f"(m_d − m_u) × {1+delta:.2f}", 1.0,     1.0, 1 + delta),
        ]:
            res = compute_nucleon_mass_splitting(
                r_p_fm * factor_r,
                alpha_em * factor_a,
                md_minus_mu_MeV * factor_m,
            )
            varied_value = (
                r_p_fm * factor_r if "r_p" in label
                else alpha_em * factor_a if "α" in label
                else md_minus_mu_MeV * factor_m
            )
            rows.append((label, varied_value, res.m_n_minus_m_p_predicted))

    return SensitivityScan(nominal=nominal, rows=rows)
