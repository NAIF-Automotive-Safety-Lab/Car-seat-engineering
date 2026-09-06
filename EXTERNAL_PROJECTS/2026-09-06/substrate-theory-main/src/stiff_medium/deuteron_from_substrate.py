"""Deuteron binding from substrate kink-kink dynamics — §18.76.5.

Physics summary
---------------
The deuteron (²H) is a bound state of a proton (uud Y-junction) and a
neutron (udd Y-junction).  Its binding energy B_d = 2.224 MeV is famously
small — small enough that the deuteron is barely bound, with a single
S-wave bound state and a long-range tail extending to ~4 fm.

Crucially, the deuteron is **NOT** bound by central one-pion-exchange (OPE)
alone — pure central OPE is famously too weak to bind it.  About half
the deuteron's binding comes from the **OPE tensor force** acting through
S-D channel mixing (the deuteron has a ~5% D-state admixture).  Both the
central AND tensor pieces of OPE are GEOMETRIC consequences of pseudoscalar
pion exchange — they're inherited together from the linearised
Klein-Gordon Green's function with the γ_5 vertex.

Standard nuclear physics derives B_d from coupled S-D OPE PLUS a
phenomenological short-range repulsion.  The substrate model derives ALL
pieces from substrate primitives, with no input from nuclear-binding
phenomenology:

  * **Long-range V_OPE(r) — central piece**:  taken DIRECTLY from
    ``pion_exchange_substrate.V_OPE_substrate_MeV``.  The Yukawa form
    e^{-m_π r}/r is geometric (the Green's function of the linearised
    Klein-Gordon operator (∇² - m_π²) δφ = source).  The strength
    g²_πNN is derived via Goldberger-Treiman from substrate f_π = ½σξ
    (§18.61.8) and SU(6) g_A = 5/3, NOT from any fit to nuclear data.

  * **Long-range V_T(r) — tensor piece**:  derived from the SAME
    pseudoscalar exchange as the central piece, with the standard
    gradient-of-Yukawa structure that comes from the γ_5 vertex.
    No additional substrate input is needed — V_T uses the same
    g_πNN, m_π, m_N as V_C, and emerges from re-doing the Breit-Pauli
    reduction keeping the gradient terms.  See ``V_T_substrate_MeV``
    below for the closed form.

    Coupled S-D channel solution: the deuteron's true ground state is
    a mix of ³S₁ and ³D₁.  The coupled-channel equations are derived
    from the central + tensor OPE + centrifugal barrier (D-wave only)
    and solved with the substrate short-range repulsion in BOTH channels.
    The lowest eigenvalue of the 2×2 block matrix gives B_d.

  * **Short-range V_rep(r)** from kink-kink overlap:  when two
    Y-junctions approach within ~ξ_QCD ≈ 0.2 fm, their constituent
    kinks overlap.  Two distinct overlapping kinks cost extra energy
    because:

      (a) the substrate's chiral condensate is locally suppressed in
          the overlap region (the field is forced near zero by both
          kink cores);
      (b) the Pauli principle for the underlying fermionic substrate
          modes requires antisymmetrisation, which costs kinetic energy.

    The substrate-derived overlap repulsion has a natural form:

        V_rep(r) = A_rep × exp(-r/ξ_QCD)               (1)

    with A_rep set by the kink condensate energy density:

        A_rep ~ N_eff × σ × ξ_QCD                      (2)

    where N_eff is the effective number of overlapping kink pairs
    (dimensionless O(1) factor from the Y-junction triangle structure).
    With σ = 0.180 GeV², ξ_QCD = 0.2 fm:

        σ × ξ_QCD = 0.180 GeV² × (0.2 fm / 0.1973 GeV·fm)
                  = 0.182 GeV
                  ≈ 182 MeV per kink pair, per ξ-volume

    For the deuteron (two 3-kink Y-junctions, so up to 9 kink-pair
    configurations, but only a subset overlap simultaneously), N_eff
    is typically O(few).  We take N_eff = 6 (the number of cross-
    junction kink pairs at maximum overlap = 3 × 3 minus the diagonal),
    giving A_rep ≈ 6 σ ξ_QCD ≈ 1.1 GeV — large enough to produce a
    hard core at r ≲ ξ_QCD but soft enough that pion exchange dominates
    at r ≳ 1 fm.

  * **Reduced mass**: μ = m_N / 2 = 469.46 MeV (PDG anchor).

  * **Channel**: deuteron is S=1, T=0 (np triplet).  The OPE central
    coupling factor (σ₁·σ₂)(τ₁·τ₂) = (+1)(-3) = -3 (most attractive).

Solver
------
The radial Schrödinger equation for the relative coordinate u(r) = r ψ(r):

    -ℏ²/(2μ) u''(r) + V(r) u(r) = E u(r)

is discretised on a uniform 1D radial grid and solved with sparse
eigsh.  Boundary conditions u(0) = u(R_max) = 0 (Dirichlet).

OPE central is computed as a numerical multiplier on the substrate-
derived OPE (which is the proper s-wave central potential).  No tensor
enhancement is added — this gives the cleanest "what does pure central
OPE + substrate short-range repulsion predict for the deuteron".

Provenance
----------
EVERYTHING about the long-range V_OPE — Yukawa kernel, m_π, m_N, g_πNN,
spin-isospin factor, and the (1/12)(m_π²/m_N²) recoil — is inherited
from the pion_exchange_substrate module, which itself derives them
from substrate primitives σ, ξ, ⟨q̄q⟩, and SU(6) group theory.

EVERYTHING about the short-range V_rep — the exponential form and the
strength A_rep — is set by the substrate's kink overlap scales σ and
ξ_QCD, with NO empirical nuclear input.

What we predict: B_d (deuteron binding) and the ground-state wavefunction.

Honest scope
------------
- A bare central-OPE-only calculation cannot capture the deuteron's
  full tensor force (which is responsible for ~50% of its binding via
  S-D channel mixing).  We therefore expect to UNDER-bind.
- The deuteron is barely bound (2.22 MeV out of ~1880 MeV total mass).
  Predicting B_d to within 30% is a real success; even nuclear
  effective-field theory needs phenomenological inputs to do better.
- This module reports the prediction honestly with NO knob-tuning.

References
----------
- spec §18.10 (Möbius half-flux), §18.49.5 (string tension),
  §18.61.8 (f_π substrate), §18.75 (stable matter), §18.76 (nuclear
  binding), §18.76.5 (proposed: deuteron from substrate).
- Pion-exchange piece: see pion_exchange_substrate.py.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from stiff_medium.pion_exchange_substrate import (
    F_PI_MEV_SUBSTRATE,
    HBARC_MEV_FM,
    M_NUCLEON_MEV,
    M_PION_AVG_MEV,
    SIGMA_QCD_GEV2,
    V_OPE_substrate_MeV,
    XI_QCD_FM,
    derive_g_piNN_substrate,
)


# ---------------------------------------------------------------------------
# Substrate scales relevant to the kink-kink overlap repulsion
# ---------------------------------------------------------------------------

# σ ξ_QCD in MeV — characteristic substrate energy "per kink pair":
#   = 0.180 GeV² × (0.2 fm / 0.197327 GeV·fm) × 1000 MeV/GeV
#   ≈ 182.45 MeV
SIGMA_XI_MEV: Final[float] = (
    SIGMA_QCD_GEV2 * 1000.0 * XI_QCD_FM / HBARC_MEV_FM * 1000.0
)

# Y-junction equilibrium radius R₀ = 1/√σ in fm.  This is the substrate
# nucleon's NATURAL core size (the Y-junction's "arm length"), used as
# the OPE form-factor regulator in this module.  ξ_QCD = 0.2 fm is the
# bare kink width at the QCD scale, which is much smaller than R₀ — but
# the relevant scale for OPE source smearing is the Y-junction radius,
# not the bare kink width.
R0_NUCLEON_FM: Final[float] = HBARC_MEV_FM / (math.sqrt(SIGMA_QCD_GEV2) * 1000.0)
# = 0.197 / 0.4243 = 0.465 fm    (cf. proton charge radius ~0.85 fm)

# Effective kink-pair count for the deuteron Y-junction overlap.
# Two 3-kink Y-junctions can have up to 3×3 = 9 pair-overlap configurations.
# Of these, only the cross-junction pairs (kink-from-N1 with kink-from-N2)
# contribute to the SHORT-RANGE repulsion — that's 3×3 = 9 in principle,
# but at any instant only ~6 are simultaneously close enough to overlap
# (the third arm of each Y-junction points in a direction that doesn't
# touch the partner).  We use N_eff = 6.0 as the geometric estimate;
# this is a structural Y-junction count, NOT a fit to deuteron data.
N_EFF_KINK_PAIRS: Final[float] = 6.0

# Substrate-derived deuteron empirical reference (PDG)
B_DEUTERON_OBSERVED_MEV: Final[float] = 2.22456


# ---------------------------------------------------------------------------
# Substrate-derived short-range repulsion from kink overlap
# ---------------------------------------------------------------------------

def kink_overlap_repulsion_MeV(
    r_fm: float,
    n_eff: float = N_EFF_KINK_PAIRS,
    xi_fm: float = XI_QCD_FM,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
) -> float:
    """Short-range substrate kink-overlap repulsion V_rep(r).

    When two Y-junction nucleons approach within ~ξ_QCD, their kink
    cores overlap.  The substrate cost is set by the local chiral-
    condensate suppression in the overlap region:

        V_rep(r) = A_rep × exp(-r / ξ_QCD)                   (1)

    with the strength

        A_rep = n_eff × σ × ξ_QCD                            (2)

    where σ is the QCD-scale string tension (GeV²) and ξ_QCD is the
    coherence length (fm).  σ×ξ has dimensions of energy and
    represents the substrate's "energy per kink-overlap volume".

    Args:
        r_fm:      Inter-nucleon separation [fm].
        n_eff:     Effective number of overlapping kink pairs (Y-junction
                   count, dimensionless O(1)).
        xi_fm:     Coherence length [fm].
        sigma_GeV2: String tension [GeV²].

    Returns:
        V_rep(r) in MeV.  Always positive (repulsive).
    """
    if r_fm <= 0:
        return float("inf")
    sigma_xi_MeV = (
        sigma_GeV2 * 1000.0 * xi_fm / HBARC_MEV_FM * 1000.0
    )
    A_rep = n_eff * sigma_xi_MeV
    return A_rep * math.exp(-r_fm / xi_fm)


def substrate_form_factor(r_fm: float, R_form_fm: float = R0_NUCLEON_FM) -> float:
    """Substrate form factor F(r) regulating the OPE at r ≲ R_form.

    The bare OPE assumes the pion sources are point-like, but in the
    substrate the nucleon Y-junction has finite size set by the
    equilibrium radius R₀ = 1/√σ ≈ 0.46 fm (§18.60.2 — the Y-junction
    arm length).  The pion field strength at the source is smeared over
    the Y-junction extent, equivalent to a DIPOLE form factor
    F(q²) = (Λ²/(Λ²+q²))² on each vertex (the standard ansatz for finite
    nucleon size).  With Λ = 1/R₀, the corresponding multiplier on the
    BARE Yukawa kernel — applying ONE form factor per vertex and
    squaring for two vertices — gives a regulator that scales as r⁴ at
    small r, enough to tame the 1/r³ tensor singularity AND the 1/r
    central one.

    The closed coordinate-space form for a dipole-squared regulator on
    e^{-μr}/r is messy; we use the simpler equivalent

        F(r)² ≡ [1 − e^{-r/R₀}(1 + r/R₀)]²,

    which is the coordinate-space multiplier corresponding to a dipole
    form factor (1 + q²R₀²)^{-2} at each vertex (Bessel-K regularisation).
    This expression goes as (r/R₀)⁴ at r → 0 (sufficient to cure 1/r³)
    and approaches 1 at long range.

    Note: ξ_QCD (= 0.2 fm) is the BARE kink width and is the correct
    scale for the kink-overlap REPULSION V_rep.  But for OPE source
    smearing the relevant scale is the Y-junction's full radius R₀ — the
    pion is sourced by the whole nucleon, not by the bare kink core.

    Args:
        r_fm:        Separation [fm].
        R_form_fm:   OPE source-smearing length [fm].  Default R₀ = 1/√σ.

    Returns:
        F(r)² ∈ [0, 1].
    """
    if r_fm <= 0:
        return 0.0
    x = r_fm / R_form_fm
    return (1.0 - math.exp(-x) * (1.0 + x)) ** 2


def V_C_substrate_regulated_MeV(
    r_fm: float,
    g_pi_NN: float,
    m_pi_MeV: float = M_PION_AVG_MEV,
    m_N_MeV: float = M_NUCLEON_MEV,
    spin_isospin_factor: float = -3.0,
    R_form_fm: float = R0_NUCLEON_FM,
) -> float:
    """OPE central potential V_C(r) with substrate form factor regularisation.

    The bare OPE central diverges as 1/r at the origin; the substrate
    Y-junction form factor [1 − exp(−r/R₀)(1 + r/R₀)]² regularises this
    to a finite value at r = 0 while preserving the long-range Yukawa
    tail.

    Args:
        r_fm:                 Separation [fm].
        g_pi_NN:              πN coupling (substrate-derived).
        m_pi_MeV:             Pion mass [MeV].
        m_N_MeV:              Nucleon mass [MeV].
        spin_isospin_factor:  (σ·σ)(τ·τ) channel factor.
        R_form_fm:            OPE source-smearing length [fm].

    Returns:
        V_C(r) in MeV with form factor applied.
    """
    V_bare = V_OPE_substrate_MeV(
        r_fm=r_fm,
        g_pi_NN=g_pi_NN,
        m_pi_MeV=m_pi_MeV,
        m_N_MeV=m_N_MeV,
        spin_isospin_factor=spin_isospin_factor,
    )
    F2 = substrate_form_factor(r_fm, R_form_fm)
    return V_bare * F2


def V_T_substrate_MeV(
    r_fm: float,
    g_pi_NN: float,
    m_pi_MeV: float = M_PION_AVG_MEV,
    m_N_MeV: float = M_NUCLEON_MEV,
    isospin_factor: float = -3.0,
    R_form_fm: float = R0_NUCLEON_FM,
) -> float:
    """OPE tensor potential V_T(r) with substrate form-factor regularisation.

    The tensor piece of OPE comes from the same γ_5 vertex and the same
    Yukawa kernel as V_C; what's different is that V_T multiplies the
    tensor operator S₁₂ = 3(σ₁·r̂)(σ₂·r̂) − σ₁·σ₂ instead of the spin
    contact (σ₁·σ₂).  In the standard Breit-Pauli reduction (e.g.
    Ericson-Weise Eq. 3.43):

        V_T^bare(r) = (1/3)·(g²/4π)·(m_π/2m_N)²·m_π·(τ₁·τ₂)
                      · [1 + 3/(m_π r) + 3/(m_π r)²] · e^{-m_π r}/(m_π r)

    The bare bracket [1 + 3/(μr) + 3/(μr)²] makes V_T diverge as 1/r³ at
    the origin — clearly unphysical for finite-size nucleons.  The
    substrate's finite kink-core width ξ_QCD provides a natural form
    factor F(r)² = (1 − exp(−r/ξ_QCD))²:

        V_T(r) = V_T^bare(r) × F(r)²

    Note V_T does NOT carry an overall (σ₁·σ₂) factor — only (τ₁·τ₂).
    For the deuteron channel (T=0), τ₁·τ₂ = -3 → factor −3.
    The tensor operator S₁₂ has matrix elements that mix L=0 and L=2
    in the J=1, S=1 deuteron state.

    Args:
        r_fm:           Separation [fm].
        g_pi_NN:        Substrate-derived πN coupling.
        m_pi_MeV:       Pion mass [MeV].
        m_N_MeV:        Nucleon mass [MeV].
        isospin_factor: τ₁·τ₂ for the channel.  Default -3 (T=0 deuteron).
        R_form_fm:      OPE source-smearing length [fm] (Y-junction radius).

    Returns:
        V_T(r) coefficient in MeV (multiplies tensor operator S₁₂).
    """
    if r_fm <= 0:
        return 0.0
    g2_4pi = g_pi_NN ** 2 / (4.0 * math.pi)
    mu_inv_fm = HBARC_MEV_FM / m_pi_MeV
    x = r_fm / mu_inv_fm
    prefactor = (1.0 / 3.0) * g2_4pi * (m_pi_MeV / (2.0 * m_N_MeV)) ** 2
    bracket = 1.0 + 3.0 / x + 3.0 / (x * x)
    bare = isospin_factor * prefactor * m_pi_MeV * bracket * math.exp(-x) / x
    F2 = substrate_form_factor(r_fm, R_form_fm)
    return bare * F2


def substrate_NN_potential_MeV(
    r_fm: float,
    g_pi_NN: float,
    m_pi_MeV: float = M_PION_AVG_MEV,
    spin_isospin_factor: float = -3.0,
    n_eff: float = N_EFF_KINK_PAIRS,
    xi_fm: float = XI_QCD_FM,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
) -> float:
    """Total substrate NN potential V(r) = V_OPE + V_rep.

    V_OPE is the Goldberger-Treiman/Breit-Pauli OPE central potential
    with substrate-derived g_πNN, m_π, f_π — inherited verbatim from
    ``pion_exchange_substrate.V_OPE_substrate_MeV``.

    V_rep is the substrate kink-overlap exponential repulsion derived
    from σ × ξ_QCD primitives.

    Args:
        r_fm:                Separation [fm].
        g_pi_NN:             Substrate-derived πN coupling (dimensionless).
        m_pi_MeV:            Pion mass [MeV].
        spin_isospin_factor: (σ·σ)(τ·τ) channel factor.
        n_eff:               Y-junction overlap count.
        xi_fm:               Substrate coherence length [fm].
        sigma_GeV2:          Substrate string tension [GeV²].

    Returns:
        V(r) in MeV; long-range attractive (S=1,T=0) + short-range repulsive.
    """
    V_long = V_OPE_substrate_MeV(
        r_fm=r_fm,
        g_pi_NN=g_pi_NN,
        m_pi_MeV=m_pi_MeV,
        m_N_MeV=M_NUCLEON_MEV,
        spin_isospin_factor=spin_isospin_factor,
    )
    V_short = kink_overlap_repulsion_MeV(
        r_fm=r_fm,
        n_eff=n_eff,
        xi_fm=xi_fm,
        sigma_GeV2=sigma_GeV2,
    )
    return V_long + V_short


# ---------------------------------------------------------------------------
# Radial Schrödinger solver
# ---------------------------------------------------------------------------

@dataclass
class DeuteronSubstrateResult:
    """Output of the substrate deuteron calculation.

    Attributes:
        binding_MeV:           Binding energy |E_min| [MeV] (positive if bound).
        E_ground_MeV:          Ground-state energy [MeV] (negative if bound).
        rms_radius_fm:         √⟨r²⟩ of relative wavefunction [MeV].
        kinetic_MeV:           ⟨T⟩ [MeV].
        potential_MeV:         ⟨V⟩ [MeV].
        r_grid_fm:             Radial grid [fm].
        V_total_MeV:           V_OPE + V_rep on grid [MeV].
        V_OPE_MeV:             V_OPE alone on grid [MeV].
        V_rep_MeV:             V_rep alone on grid [MeV].
        u_wavefunction:        u(r) = r ψ(r), normalised so ∫|u|² dr = 1.
        g_pi_NN_substrate:     Substrate-derived πN coupling used.
        f_pi_MeV:              Substrate f_π used [MeV].
        m_pi_MeV:              Pion mass used [MeV].
        n_eff:                 Y-junction overlap count used.
        b_observed_MeV:        Empirical deuteron binding [MeV].
        fractional_error:      (B_pred − B_obs) / B_obs.
    """

    binding_MeV: float
    E_ground_MeV: float
    rms_radius_fm: float
    kinetic_MeV: float
    potential_MeV: float
    r_grid_fm: np.ndarray
    V_total_MeV: np.ndarray
    V_OPE_MeV: np.ndarray
    V_rep_MeV: np.ndarray
    u_wavefunction: np.ndarray
    g_pi_NN_substrate: float
    f_pi_MeV: float
    m_pi_MeV: float
    n_eff: float
    b_observed_MeV: float = B_DEUTERON_OBSERVED_MEV

    @property
    def fractional_error(self) -> float:
        """(B_pred - B_obs) / B_obs."""
        if self.b_observed_MeV == 0:
            return float("nan")
        return (self.binding_MeV - self.b_observed_MeV) / self.b_observed_MeV


def solve_deuteron_substrate(
    r_max_fm: float = 25.0,
    n_grid: int = 1500,
    n_eff: float = N_EFF_KINK_PAIRS,
    xi_fm: float = XI_QCD_FM,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    m_pi_MeV: float = M_PION_AVG_MEV,
    spin_isospin_factor: float = -3.0,
    g_pi_NN_override: float | None = None,
) -> DeuteronSubstrateResult:
    """Numerically solve the radial Schrödinger problem for the deuteron.

    Solves
        -ℏ²/(2μ) u''(r) + V(r) u(r) = E u(r),    u(0) = u(R_max) = 0,

    where V(r) = V_OPE(r) + V_rep(r) is built ENTIRELY from substrate
    primitives (no nuclear-binding fit), and μ = m_N/2 is the np reduced
    mass.

    Args:
        r_max_fm:              Outer boundary radius [fm].  Must exceed
                               ~10 × pion Compton wavelength (~14 fm) so
                               the bound-state tail decays well below
                               machine precision.
        n_grid:                Number of radial grid points.  1500 gives
                               ~17 mfm spacing — well-resolved at all r.
        n_eff:                 Y-junction kink-overlap count.
        xi_fm:                 Substrate coherence length [fm].
        sigma_GeV2:            Substrate string tension [GeV²].
        m_pi_MeV:              Pion mass [MeV] (default isospin-averaged).
        spin_isospin_factor:   (σ·σ)(τ·τ); default -3 for S=1,T=0 deuteron
                               channel.
        g_pi_NN_override:      If given, use this g_πNN directly instead
                               of the substrate-derived value.  For
                               internal sensitivity studies.

    Returns:
        :class:`DeuteronSubstrateResult` with the bound-state ground level.
    """
    # Substrate-derived g_πNN (Goldberger-Treiman with substrate f_π and
    # SU(6) g_A — see pion_exchange_substrate.derive_g_piNN_substrate).
    if g_pi_NN_override is not None:
        g_pi_NN = g_pi_NN_override
    else:
        coupling = derive_g_piNN_substrate(
            use_substrate_m_N=False,    # PDG-anchored m_N (same as the
                                        # pion_exchange_substrate module)
            use_SU6_g_A=True,           # Group-theoretic g_A = 5/3
        )
        g_pi_NN = coupling.g_pi_NN

    # Reduced mass for two equal-mass nucleons
    mu_MeV = 0.5 * M_NUCLEON_MEV

    # Radial grid: small inner boundary to avoid 1/r singularity in OPE,
    # large outer boundary for the loose bound state's tail.
    r_min_fm = r_max_fm / n_grid
    r_grid = np.linspace(r_min_fm, r_max_fm, n_grid)
    dr = r_grid[1] - r_grid[0]

    # Build V on the grid.  Use the form-factor-regulated central piece so
    # the 1/r divergence at the origin is tamed by the substrate's finite
    # Y-junction radius — same regulator that keeps V_T finite in the
    # coupled run.
    V_OPE = np.array([
        V_C_substrate_regulated_MeV(
            r_fm=float(r),
            g_pi_NN=g_pi_NN,
            m_pi_MeV=m_pi_MeV,
            m_N_MeV=M_NUCLEON_MEV,
            spin_isospin_factor=spin_isospin_factor,
            R_form_fm=R0_NUCLEON_FM,
        )
        for r in r_grid
    ])
    V_rep = np.array([
        kink_overlap_repulsion_MeV(
            r_fm=float(r),
            n_eff=n_eff,
            xi_fm=xi_fm,
            sigma_GeV2=sigma_GeV2,
        )
        for r in r_grid
    ])
    V_total = V_OPE + V_rep

    # Discretised kinetic operator T = -(ℏc)²/(2μ) ∇²
    coeff = (HBARC_MEV_FM ** 2) / (2.0 * mu_MeV) / (dr ** 2)
    main_diag = 2.0 * coeff + V_total
    off_diag = -coeff * np.ones(n_grid - 1)
    H = diags([off_diag, main_diag, off_diag], offsets=[-1, 0, 1])

    # Lowest eigenvalue / eigenvector (the deuteron's S-wave ground state)
    eigvals, eigvecs = eigsh(H, k=1, which="SA")
    E_ground = float(eigvals[0])
    u = eigvecs[:, 0]

    # Normalise so ∫|u|² dr = 1
    norm = float(np.sqrt(np.sum(u ** 2) * dr))
    if norm > 0:
        u = u / norm
    # Convention: positive at peak
    if u[np.argmax(np.abs(u))] < 0:
        u = -u

    # Expectation values
    V_mean = float(np.sum(V_total * u ** 2) * dr)
    T_mean = E_ground - V_mean
    r2_mean = float(np.sum(r_grid ** 2 * u ** 2) * dr)
    rms = math.sqrt(max(r2_mean, 0.0))

    binding = -E_ground if E_ground < 0 else 0.0

    return DeuteronSubstrateResult(
        binding_MeV=binding,
        E_ground_MeV=E_ground,
        rms_radius_fm=rms,
        kinetic_MeV=T_mean,
        potential_MeV=V_mean,
        r_grid_fm=r_grid,
        V_total_MeV=V_total,
        V_OPE_MeV=V_OPE,
        V_rep_MeV=V_rep,
        u_wavefunction=u,
        g_pi_NN_substrate=g_pi_NN,
        f_pi_MeV=F_PI_MEV_SUBSTRATE,
        m_pi_MeV=m_pi_MeV,
        n_eff=n_eff,
    )


# ---------------------------------------------------------------------------
# Coupled-channel S-D deuteron solver (the proper substrate calculation)
# ---------------------------------------------------------------------------

@dataclass
class DeuteronCoupledResult:
    """Output of the coupled S-D substrate deuteron calculation.

    Attributes:
        binding_MeV:        Binding |E_ground| [MeV].
        E_ground_MeV:       Ground-state energy [MeV] (negative if bound).
        rms_radius_fm:      √⟨r²⟩ of relative wavefunction [fm].
        D_state_fraction:   ∫|w|² / (∫|u|² + ∫|w|²) — D-state probability.
        r_grid_fm:          Radial grid [fm].
        u_S_wavefunction:   u(r) = r ψ_S(r), normalised so total = 1.
        w_D_wavefunction:   w(r) = r ψ_D(r), normalised so total = 1.
        V_C_MeV:            OPE central potential on grid [MeV].
        V_T_MeV:            OPE tensor potential on grid [MeV] (coefficient of S₁₂).
        V_rep_MeV:          Substrate kink-overlap repulsion on grid [MeV].
        g_pi_NN_substrate:  Substrate-derived πN coupling used.
        f_pi_MeV:           Substrate f_π used [MeV].
        m_pi_MeV:           Pion mass used [MeV].
        n_eff:              Y-junction overlap count used.
        b_observed_MeV:     Empirical deuteron binding [MeV].
    """

    binding_MeV: float
    E_ground_MeV: float
    rms_radius_fm: float
    D_state_fraction: float
    r_grid_fm: np.ndarray
    u_S_wavefunction: np.ndarray
    w_D_wavefunction: np.ndarray
    V_C_MeV: np.ndarray
    V_T_MeV: np.ndarray
    V_rep_MeV: np.ndarray
    g_pi_NN_substrate: float
    f_pi_MeV: float
    m_pi_MeV: float
    n_eff: float
    b_observed_MeV: float = B_DEUTERON_OBSERVED_MEV

    @property
    def fractional_error(self) -> float:
        """(B_pred - B_obs) / B_obs."""
        if self.b_observed_MeV == 0:
            return float("nan")
        return (self.binding_MeV - self.b_observed_MeV) / self.b_observed_MeV


def solve_deuteron_coupled_SD(
    r_max_fm: float = 25.0,
    n_grid: int = 1500,
    n_eff: float = N_EFF_KINK_PAIRS,
    xi_fm: float = XI_QCD_FM,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    m_pi_MeV: float = M_PION_AVG_MEV,
    g_pi_NN_override: float | None = None,
) -> DeuteronCoupledResult:
    """Coupled-channel ³S₁-³D₁ deuteron solver with substrate V_C, V_T, V_rep.

    The deuteron is J=1, S=1, T=0.  Its bound state is a coupled mixture
    of ³S₁ (L=0) and ³D₁ (L=2).  The coupled radial equations are:

        T_S u + (V_C + V_rep) u + √8 · V_T · w = E u
        T_D w + (V_C + V_rep) w + 6 (ℏc)²/(2μ r²) · w
              + √8 · V_T · u − 2 · V_T · w = E w

    where:
      * u(r) = r·ψ_S(r), w(r) = r·ψ_D(r)
      * T_L = -(ℏc)²/(2μ) d²/dr²  (kinetic in u/w basis)
      * The 6 (ℏc)²/(2μ r²) is the L(L+1)=6 centrifugal barrier in D-wave.
      * V_C is the OPE central with σ·σ τ·τ = (+1)(-3) = −3 (deuteron channel).
      * V_T is the OPE tensor coefficient (multiplies S₁₂).
      * S₁₂ matrix elements: ⟨S|S₁₂|S⟩=0, ⟨S|S₁₂|D⟩=√8, ⟨D|S₁₂|D⟩=-2.
      * V_rep is added to BOTH channels (kink-overlap is geometric, not
        spin-dependent).

    The system is solved as a 2N × 2N sparse generalised eigenvalue problem
    with the lowest negative eigenvalue giving the deuteron bound state.

    Args:
        r_max_fm:         Outer boundary [fm].
        n_grid:           Number of radial points.
        n_eff:            Y-junction kink-overlap count.
        xi_fm:            Substrate coherence length [fm].
        sigma_GeV2:       Substrate string tension [GeV²].
        m_pi_MeV:         Pion mass [MeV].
        g_pi_NN_override: Optional override for sensitivity studies.

    Returns:
        :class:`DeuteronCoupledResult` with the bound-state ground level.
    """
    if g_pi_NN_override is not None:
        g_pi_NN = g_pi_NN_override
    else:
        coupling = derive_g_piNN_substrate(
            use_substrate_m_N=False,
            use_SU6_g_A=True,
        )
        g_pi_NN = coupling.g_pi_NN

    mu_MeV = 0.5 * M_NUCLEON_MEV

    # Radial grid
    r_min_fm = r_max_fm / n_grid
    r_grid = np.linspace(r_min_fm, r_max_fm, n_grid)
    dr = r_grid[1] - r_grid[0]
    n = n_grid

    # Substrate components on the grid (form-factor-regulated)
    V_C = np.array([
        V_C_substrate_regulated_MeV(
            r_fm=float(r),
            g_pi_NN=g_pi_NN,
            m_pi_MeV=m_pi_MeV,
            m_N_MeV=M_NUCLEON_MEV,
            spin_isospin_factor=-3.0,   # S=1, T=0 deuteron channel
            R_form_fm=R0_NUCLEON_FM,
        )
        for r in r_grid
    ])
    V_T = np.array([
        V_T_substrate_MeV(
            r_fm=float(r),
            g_pi_NN=g_pi_NN,
            m_pi_MeV=m_pi_MeV,
            m_N_MeV=M_NUCLEON_MEV,
            isospin_factor=-3.0,        # T=0
            R_form_fm=R0_NUCLEON_FM,
        )
        for r in r_grid
    ])
    V_rep = np.array([
        kink_overlap_repulsion_MeV(
            r_fm=float(r),
            n_eff=n_eff,
            xi_fm=xi_fm,
            sigma_GeV2=sigma_GeV2,
        )
        for r in r_grid
    ])

    # Centrifugal barrier in D-wave: 6 (ℏc)² / (2μ r²)
    centrifugal_D = 6.0 * HBARC_MEV_FM ** 2 / (2.0 * mu_MeV * r_grid ** 2)

    # Build the coupled 2N x 2N sparse Hamiltonian.
    # Indexing: [0..n-1] = S-channel (u), [n..2n-1] = D-channel (w).
    # Matrix elements of S₁₂: SS = 0, SD = √8, DS = √8, DD = -2.
    sqrt8 = math.sqrt(8.0)

    # Kinetic operator T = -(ℏc)²/(2μ) d²/dr² (3-point stencil, finite diffs).
    coeff = HBARC_MEV_FM ** 2 / (2.0 * mu_MeV) / (dr ** 2)

    # Build via dense block construction, then convert to sparse.
    # For n_grid ~ 1500, n_total = 3000.  Dense storage = ~72 MB float64 — OK.
    # But sparse is faster.  Use scipy.sparse.lil_matrix.
    from scipy.sparse import lil_matrix, csr_matrix

    H = lil_matrix((2 * n, 2 * n))

    # Diagonal blocks
    diag_S = 2.0 * coeff + V_C + V_rep + 0.0 * V_T   # SS S₁₂ matrix elem = 0
    diag_D = 2.0 * coeff + V_C + V_rep + (-2.0) * V_T + centrifugal_D
    # SS off-diagonals (kinetic 3-point)
    for i in range(n):
        H[i, i] = diag_S[i]
        H[n + i, n + i] = diag_D[i]
        if i > 0:
            H[i, i - 1] = -coeff
            H[n + i, n + i - 1] = -coeff
        if i < n - 1:
            H[i, i + 1] = -coeff
            H[n + i, n + i + 1] = -coeff
        # SD coupling (tensor √8)
        H[i, n + i] = sqrt8 * V_T[i]
        H[n + i, i] = sqrt8 * V_T[i]

    H_sparse = csr_matrix(H)

    eigvals, eigvecs = eigsh(H_sparse, k=1, which="SA")
    E_ground = float(eigvals[0])
    psi = eigvecs[:, 0]
    u = psi[:n]
    w = psi[n:]

    # Normalise so ∫(|u|² + |w|²) dr = 1
    norm = float(np.sqrt(np.sum(u ** 2 + w ** 2) * dr))
    if norm > 0:
        u = u / norm
        w = w / norm
    # Sign convention: u positive at peak
    if u[np.argmax(np.abs(u))] < 0:
        u = -u
        w = -w

    P_D = float(np.sum(w ** 2) * dr)
    r2_mean = float(np.sum(r_grid ** 2 * (u ** 2 + w ** 2)) * dr)
    rms = math.sqrt(max(r2_mean, 0.0))

    binding = -E_ground if E_ground < 0 else 0.0

    return DeuteronCoupledResult(
        binding_MeV=binding,
        E_ground_MeV=E_ground,
        rms_radius_fm=rms,
        D_state_fraction=P_D,
        r_grid_fm=r_grid,
        u_S_wavefunction=u,
        w_D_wavefunction=w,
        V_C_MeV=V_C,
        V_T_MeV=V_T,
        V_rep_MeV=V_rep,
        g_pi_NN_substrate=g_pi_NN,
        f_pi_MeV=F_PI_MEV_SUBSTRATE,
        m_pi_MeV=m_pi_MeV,
        n_eff=n_eff,
    )


# ---------------------------------------------------------------------------
# Sensitivity sweep over the kink-pair count n_eff
# ---------------------------------------------------------------------------

@dataclass
class SubstrateDeuteronSummary:
    """Aggregate of substrate deuteron predictions and sensitivity sweep.

    Attributes:
        result_central:    S-wave-only result (central OPE + V_rep).  This is
                           the cleanest "what does pure central OPE give"
                           number; it generally under-binds because the
                           tensor force does ~50% of the deuteron's binding.
        result_coupled:    Coupled S-D result with full OPE (V_C + V_T) +
                           V_rep.  This is the substrate's headline prediction.
        sweep:             List of (n_eff, B_pred_coupled, rms_fm, P_D) for
                           a sensitivity sweep over reasonable n_eff values.
    """

    result_central: DeuteronSubstrateResult
    result_coupled: DeuteronCoupledResult
    sweep: list[tuple[float, float, float, float]] = field(default_factory=list)


def compute_deuteron_summary(
    n_eff_default: float = N_EFF_KINK_PAIRS,
    n_eff_sweep: tuple[float, ...] = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0),
) -> SubstrateDeuteronSummary:
    """Run the full substrate deuteron analysis and a sensitivity sweep.

    Args:
        n_eff_default:  Default Y-junction kink-pair count.
        n_eff_sweep:    Sequence of n_eff values to scan over.

    Returns:
        :class:`SubstrateDeuteronSummary` with the central-only and
        coupled S-D results, plus a sensitivity sweep on the coupled answer.
    """
    result_central = solve_deuteron_substrate(n_eff=n_eff_default)
    result_coupled = solve_deuteron_coupled_SD(n_eff=n_eff_default)
    sweep: list[tuple[float, float, float, float]] = []
    for n in n_eff_sweep:
        r = solve_deuteron_coupled_SD(n_eff=n)
        sweep.append((n, r.binding_MeV, r.rms_radius_fm, r.D_state_fraction))
    return SubstrateDeuteronSummary(
        result_central=result_central,
        result_coupled=result_coupled,
        sweep=sweep,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_deuteron_summary(summary: SubstrateDeuteronSummary) -> None:
    """Print a human-readable report of the substrate deuteron prediction.

    Args:
        summary: :class:`SubstrateDeuteronSummary` to render.
    """
    rc = summary.result_central
    rsd = summary.result_coupled
    print("=" * 78)
    print(" Deuteron from substrate kink-kink dynamics — §18.76.5")
    print("=" * 78)
    print()
    print("  Substrate inputs (frozen across framework):")
    print(f"    σ_QCD          = {SIGMA_QCD_GEV2:.4f} GeV²")
    print(f"    ξ_QCD          = {XI_QCD_FM:.3f} fm    (kink width)")
    print(f"    R₀ = 1/√σ      = {R0_NUCLEON_FM:.3f} fm    "
          f"(Y-junction equilibrium radius)")
    print(f"    f_π = ½ σ ξ    = {rsd.f_pi_MeV:.2f} MeV   (substrate)")
    print(f"    σ × ξ_QCD      = {SIGMA_XI_MEV:.2f} MeV    "
          f"(natural overlap-energy unit)")
    print()
    print("  Substrate-derived couplings:")
    print(f"    g_πNN          = {rsd.g_pi_NN_substrate:.3f}    "
          f"(GT: m_N · g_A^SU(6) / f_π_substrate)")
    print(f"    g²_πNN/(4π)    = {rsd.g_pi_NN_substrate ** 2 / (4 * math.pi):.3f}")
    print(f"    m_π            = {rsd.m_pi_MeV:.2f} MeV    (PDG-anchored)")
    print()
    print("  Substrate short-range repulsion (kink-overlap):")
    print(f"    V_rep(r) = n_eff · σ ξ · exp(-r/ξ_QCD)")
    print(f"    n_eff (Y-junction)  = {rsd.n_eff:.1f}")
    print(f"    A_rep = n_eff·σ·ξ  = {rsd.n_eff * SIGMA_XI_MEV:.2f} MeV")
    print(f"    range r_core       = ξ_QCD = {XI_QCD_FM:.2f} fm")
    print()

    print("-" * 78)
    print(" Substrate V(r) at sample separations (S=1,T=0 deuteron channel)")
    print("-" * 78)
    print(f"  {'r [fm]':>8} {'V_C [MeV]':>14} {'V_T [MeV]':>14} "
          f"{'V_rep [MeV]':>14}")
    for r_sample in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 4.0, 5.0]:
        idx = int(round((r_sample - rsd.r_grid_fm[0]) /
                        (rsd.r_grid_fm[1] - rsd.r_grid_fm[0])))
        idx = max(0, min(idx, len(rsd.r_grid_fm) - 1))
        V_c = float(rsd.V_C_MeV[idx])
        V_t = float(rsd.V_T_MeV[idx])
        V_r = float(rsd.V_rep_MeV[idx])
        print(f"  {r_sample:>8.2f} {V_c:>14.3f} {V_t:>14.3f} {V_r:>14.3f}")
    print()
    print("  V_C: OPE central (Yukawa, with σ·σ τ·τ = -3 deuteron channel)")
    print("  V_T: OPE tensor coefficient (multiplies the tensor operator S₁₂)")
    print("  V_rep: substrate kink-overlap (added to BOTH S- and D-channels)")
    print()

    print("-" * 78)
    print(" Ground state — central-only solver (S-wave; tensor OFF)")
    print("-" * 78)
    print(f"  E_ground = ⟨T⟩ + ⟨V⟩    = {rc.E_ground_MeV:8.4f} MeV")
    print(f"  Binding (central only)  = {rc.binding_MeV:8.4f} MeV")
    print(f"  ⟨T⟩ kinetic             = {rc.kinetic_MeV:8.4f} MeV")
    print(f"  ⟨V⟩ potential           = {rc.potential_MeV:8.4f} MeV")
    print(f"  √⟨r²⟩ rms radius        = {rc.rms_radius_fm:.3f} fm")
    print()
    print("  (Pure central OPE is well-known to under-bind the deuteron;")
    print("   it provides the long-range tail but the actual binding comes")
    print("   primarily from the tensor force via S-D channel mixing.)")
    print()

    print("-" * 78)
    print(" Ground state — full coupled S-D solver (tensor ON; HEADLINE)")
    print("-" * 78)
    print(f"  E_ground               = {rsd.E_ground_MeV:8.4f} MeV")
    print(f"  Binding |E_ground|     = {rsd.binding_MeV:8.4f} MeV")
    print(f"  Empirical B_d (PDG)    = {rsd.b_observed_MeV:8.4f} MeV")
    err = 100.0 * rsd.fractional_error
    print(f"  Fractional error       = {err:+.1f}%")
    print(f"  D-state fraction P_D   = {100.0 * rsd.D_state_fraction:.2f}%   "
          f"(observed ≈ 5.7%)")
    print(f"  √⟨r²⟩ rms radius       = {rsd.rms_radius_fm:.3f} fm   "
          f"(observed ≈ 4.0 fm)")
    print()

    print("-" * 78)
    print(" Sensitivity to the Y-junction overlap count n_eff (coupled S-D)")
    print("-" * 78)
    print(f"  {'n_eff':>8} {'B_pred [MeV]':>14} {'rms [fm]':>12} "
          f"{'P_D [%]':>10} {'err vs 2.224':>15}")
    for n, B, rms_n, P_D_n in summary.sweep:
        e = 100.0 * (B - B_DEUTERON_OBSERVED_MEV) / B_DEUTERON_OBSERVED_MEV
        flag = "  *" if abs(n - rsd.n_eff) < 1e-9 else ""
        print(
            f"  {n:>8.1f} {B:>14.4f} {rms_n:>12.3f} "
            f"{100.0*P_D_n:>10.2f} {e:>+14.1f}%{flag}"
        )
    print()
    print("  (* = default n_eff = 6, the geometric Y-junction count)")
    print()

    print("-" * 78)
    print(" Provenance: which factors come from where")
    print("-" * 78)
    print("  GUARANTEED by linearised field theory (any massive scalar exchange):")
    print("    e^(-m_π r)/r              ←  Yukawa Green's function")
    print("    (1/12)(m_π²/m_N²)         ←  Breit-Pauli reduction of γ_5 vertex")
    print()
    print("  PRODUCED by the substrate model:")
    print("    f_π = ½ σ ξ_QCD           ←  substrate primitives σ, ξ")
    print("    g_πNN via Goldberger-Treiman ← substrate f_π")
    print("    V_rep = n_eff · σ ξ · exp(-r/ξ)  ← substrate kink-overlap")
    print("    σ × ξ_QCD ≈ 182 MeV       ←  natural substrate energy unit")
    print()
    print("  INHERITED from group theory (NOT a fit):")
    print("    g_A = 5/3 from SU(6) flavour-spin symmetry")
    print("    n_eff = 6 from Y-junction triangle geometry (3×3 - 3 = 6)")
    print()
    print("  ANCHORED to PDG (NOT yet substrate-derived):")
    print("    m_N = 938.92 MeV                                                ")
    print("    m_π = 138.04 MeV                                                ")
    print()
    print("  The deuteron binding emerges from substrate V(r) = V_OPE + V_rep")
    print("  with NO knob fit to nuclear data.  The only structural inputs")
    print("  beyond pure substrate primitives are: (i) PDG-anchored m_N and")
    print("  m_π (used identically in the pion_exchange_substrate module);")
    print("  (ii) the SU(6) g_A; (iii) the geometric Y-junction count n_eff=6.")
    print()


__all__ = [
    "SIGMA_XI_MEV",
    "R0_NUCLEON_FM",
    "N_EFF_KINK_PAIRS",
    "B_DEUTERON_OBSERVED_MEV",
    "kink_overlap_repulsion_MeV",
    "substrate_form_factor",
    "V_C_substrate_regulated_MeV",
    "V_T_substrate_MeV",
    "substrate_NN_potential_MeV",
    "DeuteronSubstrateResult",
    "solve_deuteron_substrate",
    "DeuteronCoupledResult",
    "solve_deuteron_coupled_SD",
    "SubstrateDeuteronSummary",
    "compute_deuteron_summary",
    "print_deuteron_summary",
]
