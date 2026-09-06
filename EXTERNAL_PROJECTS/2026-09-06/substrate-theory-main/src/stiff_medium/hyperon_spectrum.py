"""Hyperon mass spectrum from substrate primitives — §18.64.

Physics summary
---------------
Hyperons are baryons that contain at least one strange quark.  In the
substrate model, each "quark" is a Möbius half-flux kink bound in a Y-junction
geometry of length R₀ ≈ 1/√σ (§18.60.2).  The light (u, d) constituent kink
mass is set by the geometric scale m_K_eff = √σ ≈ 424 MeV (§18.61.1) and
absorbs all the structural mass into a single number.  The strange kink
differs from the light kink in two SUBSTRATE-NATURAL ways:

  (1) it lives on a slightly larger pocket of the running stiffness K(ξ)
      (i.e. its effective coherence length ξ_s > ξ_QCD), and

  (2) its constituent inertia is correspondingly larger,
          m_s_const = m_K_eff + Δ_s,
      with Δ_s the only flavour input and SET BY ONE HYPERON ANCHOR.

Below we adopt the simplest substrate-natural decomposition:

    m_baryon = Σᵢ m_qᵢ_const + ΔE_chromomagnetic(spin, flavour)              (1)

with the chromomagnetic shift coming from the Möbius spin-spin contact term
(§18.61.1):

    ΔE_chromo = K_substrate × Σ_{i<j} (S_i · S_j) / (m_qᵢ m_qⱼ)               (2)

    K_substrate = (8π/3) α_M C_F |ψ(0)|²
                = (8/3) σ ξ_QCD² σ^{3/2}        (collapsed form)             (3)

For a baryon with three identical light quarks (uud / udd) eq. (2) reduces to

    ΔE_chromo(N) = K_substrate × (S² - 9/4)/2 / m_q²                        (4)
                 = K_substrate × (-3/4) / m_q²    (J = 1/2)

and the N-Δ splitting (§18.61.1) is

    Δm_NΔ = K_substrate × (3/2) / m_q² = 4 ξ² σ^{3/2}.                      (5)

Inverting (5) gives the chromomag pair coupling:

    (chromo per pair) = (2/3) × Δm_NΔ = (8/3) ξ² σ^{3/2}.                   (6)

Anchoring the absolute baryon-mass scale to the proton (m_p = 939 MeV):

    m_p = 3 m_q_struct - (3/4) × (chromo per pair) / m_q_struct²            (7)

This implicit equation determines the "structure mass" m_q_struct ≈ 363 MeV
that absorbs the proton's confinement binding.  It differs slightly from the
chromomag mass m_K_eff = √σ used INSIDE the spin-spin contact formula —
the structure mass is a hadronisation parameter (kinetic plus residual
binding) while the chromomag mass is the geometric ℏc/R₀.  Both are
substrate-derived from σ alone.

Hyperons in the same scheme
---------------------------
For a baryon with strange content, eq. (1) generalises to

    m_baryon = N_q_light × m_q_struct + N_s × m_s_struct + ΔE_chromo,       (8)

where m_s_struct = m_q_struct + Δ_s sets the strange-kink constituent mass.
The chromomag shift (2) becomes flavour-dependent:

    ΔE_chromo = K_subst × [
        Σ_{i<j ∈ light}    (S_i·S_j) / m_q²
      + Σ_{i ∈ light, j=s} (S_i·S_j) / (m_q m_s)
      + Σ_{i<j ∈ s}        (S_i·S_j) / m_s²
    ].                                                                       (9)

The three flavour-dependent pair-coupling buckets are computed from m_q,
m_s, and the ground-state spin coupling appropriate to each baryon
(SU(3) Clebsch-Gordan from §18.49 inheritance).

Spin-flavour coupling matrix
----------------------------
For the J=1/2 octet, the standard SU(6) result for ⟨Σ S_i·S_j /(m_i m_j)⟩
in a baryon B is the Sakharov-like sum:

  N (uud), J=1/2:    -3/4 × 1/m_q²
  Σ⁺ (uus), J=1/2:  +1/4 × 1/m_q²  -  1   × 1/(m_q m_s)
  Σ⁰/Σ⁻  same form (uds with mass-symmetric average)
  Λ (uds), J=1/2:   -3/4 × 1/m_q²    (pure light pair, s spectator)
  Ξ⁰ (uss), J=1/2:  +1/4 × 1/m_s²  -  1   × 1/(m_q m_s)
  Ω⁻ (sss), J=3/2:  +3/4 × 1/m_s²

These coefficients come from the standard quark model (e.g. De Rújula,
Georgi, Glashow 1975) and are inherited identically by the substrate via
§18.49 SU(3) gauge inheritance plus spin-½ kink statistics (§18.10).

Anchor and predictions
----------------------
We FIX m_s_struct from the Λ hyperon (cleanest because Λ has a pure light
spin singlet — only the s spectator differs from N's spin structure), then
PREDICT all other hyperons.  In an SU(3)-symmetric model with one extra
input (Δ_s), the resulting predictions are zero-fit checks of the
substrate's flavour algebra.

Honest scope
------------
* The substrate does NOT yet predict m_s_struct from first principles.
  The strange-quark-running of K(ξ) at ξ_s > ξ_QCD would in principle
  pin Δ_s but we have not derived that running independently.  We thus
  treat m_s_struct as ONE flavour anchor.

* SU(3) symmetry breaking beyond the Δ_s shift (e.g. Σ-Λ mass difference,
  Σ⁺-Σ⁻ EM splitting) is not modelled here; expect 1-2% residual errors.

* The Σ-Λ splitting is a non-trivial test: both are uds with J=1/2, and
  the only difference is the spin coupling of the light ud pair (singlet
  for Λ, triplet for Σ).  Eq. (9) captures this through the spin
  coefficients without any extra input.

References
----------
spec §18.49 (SU(3) inheritance), §18.49.5 (string tension), §18.60.2
(R₀ = 1/√σ), §18.61.1 (chromomag N-Δ splitting), §18.62.2 (n-p splitting),
§18.64 (hyperon spectrum).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Substrate primitives at the QCD scale (inherited from chromomagnetic_substrate)
# ---------------------------------------------------------------------------

SIGMA_QCD_GEV2: Final[float] = 0.180
XI_QCD_FM: Final[float] = 0.2
HBARC_GEVFM: Final[float] = 0.197326980
XI_QCD_INV_GEV: Final[float] = XI_QCD_FM / HBARC_GEVFM     # ≈ 1.0135 GeV⁻¹

N_COLOURS: Final[int] = 3
C_F: Final[float] = (N_COLOURS ** 2 - 1.0) / (2.0 * N_COLOURS)   # 4/3

# Substrate constituent kink mass for the chromomagnetic propagator
M_K_CHROMO_GEV: Final[float] = math.sqrt(SIGMA_QCD_GEV2)   # √σ ≈ 0.4243 GeV

# Empirical baryon masses (PDG 2024) — kept for comparison only
M_PROTON_MEV: Final[float] = 938.272
M_NEUTRON_MEV: Final[float] = 939.565
M_NUCLEON_MEAN_MEV: Final[float] = 0.5 * (M_PROTON_MEV + M_NEUTRON_MEV)   # 938.92
M_DELTA_MEV: Final[float] = 1232.0
DELTA_M_NDELTA_MEV: Final[float] = M_DELTA_MEV - M_NUCLEON_MEAN_MEV       # ≈ 293.1 MeV

M_LAMBDA_MEV: Final[float] = 1115.683
M_SIGMA_PLUS_MEV: Final[float] = 1189.37
M_SIGMA_ZERO_MEV: Final[float] = 1192.642
M_SIGMA_MINUS_MEV: Final[float] = 1197.449
M_SIGMA_MEAN_MEV: Final[float] = (
    M_SIGMA_PLUS_MEV + M_SIGMA_ZERO_MEV + M_SIGMA_MINUS_MEV
) / 3.0     # 1193.15
M_XI_ZERO_MEV: Final[float] = 1314.86
M_XI_MINUS_MEV: Final[float] = 1321.71
M_XI_MEAN_MEV: Final[float] = 0.5 * (M_XI_ZERO_MEV + M_XI_MINUS_MEV)      # 1318.29
M_OMEGA_MEV: Final[float] = 1672.45


# ---------------------------------------------------------------------------
# K_substrate from substrate primitives — §18.61.1 contact coefficient
# ---------------------------------------------------------------------------

def K_substrate_GeV3(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
) -> float:
    """K_substrate = (8π/3) α_M C_F |ψ(0)|² = (8/3) σ ξ² σ^{3/2}.

    This is the chromomagnetic contact-coefficient that controls the
    spin-spin splittings.  Only σ and ξ enter; numerical prefactor (8/3)
    arises from collapsing (8π/3) × (4/3) × (3/4π) → 8/3.

    Args:
        sigma_GeV2: String tension [GeV²].
        xi_inv_GeV: Coherence length expressed in GeV⁻¹.

    Returns:
        K_substrate in GeV³.
    """
    return (8.0 / 3.0) * sigma_GeV2 * xi_inv_GeV ** 2 * sigma_GeV2 ** 1.5


def chromomag_pair_coupling_MeV(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
) -> float:
    """Per-pair chromomag coupling, evaluated at light-quark mass m_q = √σ.

    K_substrate / m_q² is the "natural" chromomag pair coupling for the
    nucleon-Δ system (see eq. 5: Δm_NΔ = (3/2) × K/m_q²).  Equivalently:

        K / m_q² = (8/3) σ^{3/2} ξ² × (1/σ) = (8/3) ξ² σ^{1/2}
                 = (2/3) × Δm_NΔ

    Returns:
        Per-pair chromomag coupling in MeV.
    """
    K = K_substrate_GeV3(sigma_GeV2, xi_inv_GeV)
    m_q2 = M_K_CHROMO_GEV ** 2
    return 1000.0 * K / m_q2


# ---------------------------------------------------------------------------
# Spin-flavour coefficient table (J = 1/2 octet + Ω⁻ J = 3/2 decuplet)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HyperonSpinFlavour:
    """Pair-coupling decomposition for one baryon.

    Attributes:
        name:           PDG-style name.
        n_light:        Number of light (u, d) constituent kinks.
        n_strange:      Number of strange constituent kinks.
        c_qq:           Coefficient of 1/m_q² (light-light pair sum, ⟨S·S⟩).
        c_qs:           Coefficient of 1/(m_q m_s) (light-strange pair sum).
        c_ss:           Coefficient of 1/m_s² (strange-strange pair sum).
        J:              Total spin (½ for octet, 3/2 for Ω⁻).
        mass_pdg_MeV:   Empirical mass [MeV].
    """
    name: str
    n_light: int
    n_strange: int
    c_qq: float
    c_qs: float
    c_ss: float
    J: float
    mass_pdg_MeV: float


# Spin coefficients from De Rújula-Georgi-Glashow analysis of the J=1/2
# octet and the Δ/Ω decuplet.  Same Clebsch-Gordan algebra as standard
# constituent quark model — substrate inherits via §18.49 SU(3).
#
# Conventions:
#   * Each coefficient is the matrix element of the spin pair sum
#     Σ_{i<j} (S_i · S_j) restricted to that pair-flavour bucket.
#   * For mixed quark masses these are evaluated using the Sakharov
#     spin-symmetry expansion (J(J+1) - 3/4 N decomposition projected
#     onto SU(3) flavour subgroups).
#
# Reference values (textbook, e.g. Halzen-Martin §6.7):
#   N (uud)   J=1/2:  c_qq = -3/4
#   Σ (uus / dds): c_qq = +1/4 (light pair triplet),  c_qs = -1
#   Λ (uds)  J=1/2:  c_qq = -3/4 (light pair singlet), c_qs = 0 (s spectator)
#   Σ⁰ (uds): same coefficients as Σ⁺ since spin structure is Σ-like
#   Ξ (uss / dss): c_qs = -1, c_ss = +1/4
#   Ω⁻ (sss) J=3/2:  c_ss = +3/4

HYPERON_TABLE: list[HyperonSpinFlavour] = [
    HyperonSpinFlavour(
        name="N(939)",
        n_light=3, n_strange=0,
        c_qq=-0.75, c_qs=0.0, c_ss=0.0,
        J=0.5, mass_pdg_MeV=M_NUCLEON_MEAN_MEV,
    ),
    HyperonSpinFlavour(
        name="Λ",
        n_light=2, n_strange=1,
        c_qq=-0.75, c_qs=0.0, c_ss=0.0,
        J=0.5, mass_pdg_MeV=M_LAMBDA_MEV,
    ),
    HyperonSpinFlavour(
        name="Σ",
        n_light=2, n_strange=1,
        c_qq=+0.25, c_qs=-1.0, c_ss=0.0,
        J=0.5, mass_pdg_MeV=M_SIGMA_MEAN_MEV,
    ),
    HyperonSpinFlavour(
        name="Ξ",
        n_light=1, n_strange=2,
        c_qq=0.0, c_qs=-1.0, c_ss=+0.25,
        J=0.5, mass_pdg_MeV=M_XI_MEAN_MEV,
    ),
    HyperonSpinFlavour(
        name="Ω⁻",
        n_light=0, n_strange=3,
        c_qq=0.0, c_qs=0.0, c_ss=+0.75,
        J=1.5, mass_pdg_MeV=M_OMEGA_MEV,
    ),
    HyperonSpinFlavour(
        name="Δ(1232)",
        n_light=3, n_strange=0,
        c_qq=+0.75, c_qs=0.0, c_ss=0.0,
        J=1.5, mass_pdg_MeV=M_DELTA_MEV,
    ),
]


# ---------------------------------------------------------------------------
# Structure mass m_q_struct from the proton anchor
# ---------------------------------------------------------------------------

def solve_m_q_struct_MeV(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
    m_target_MeV: float = M_NUCLEON_MEAN_MEV,
) -> float:
    """Solve m_p = 3 m_q_struct + ΔE_chromo(N) for the structure mass.

    The chromomag pair coupling is evaluated at the GEOMETRIC mass m_K_chromo = √σ
    (this is the natural propagator scale — the substrate's R₀-derived inertia).
    The STRUCTURE mass m_q_struct is what appears as the LINEAR Σ m_q in the
    baryon mass formula and absorbs the residual confinement binding.

    Equation:
        m_target = 3 m_q_struct + K_subst × (-3/4) / m_K_chromo²
        ⇒ m_q_struct = (m_target + (3/4) K_subst / m_K_chromo²) / 3
                     = (m_target + (3/4) × Δm_pair) / 3

    With m_target = 938.92 MeV, Δm_pair = chromomag pair coupling ≈ 209 MeV,
    we find m_q_struct ≈ 365 MeV.

    Args:
        sigma_GeV2: String tension [GeV²].
        xi_inv_GeV: Coherence length [GeV⁻¹].
        m_target_MeV: Anchor mass [MeV] (default: nucleon mean).

    Returns:
        m_q_struct in MeV.
    """
    chromo_pair_MeV = chromomag_pair_coupling_MeV(sigma_GeV2, xi_inv_GeV)
    # Nucleon: ΔE = -3/4 × chromo_pair
    delta_E_N_MeV = -0.75 * chromo_pair_MeV
    return (m_target_MeV - delta_E_N_MeV) / 3.0


# ---------------------------------------------------------------------------
# Hyperon mass formula
# ---------------------------------------------------------------------------

def baryon_mass_MeV(
    spec: HyperonSpinFlavour,
    m_q_struct_MeV: float,
    m_s_struct_MeV: float,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
    use_struct_mass_in_chromo: bool = False,
) -> float:
    """Baryon mass from substrate structure + chromomag spin-spin.

    Mass formula (§18.64):

        m_B = N_q × m_q_struct + N_s × m_s_struct + K × (
                 c_qq / m_q² + c_qs / (m_q m_s) + c_ss / m_s² )

    The chromomag mass that enters the propagator (1/m_q m_s) defaults to
    the GEOMETRIC scale (m_K_chromo = √σ for light, m_s_chromo = m_K_chromo
    × m_s_struct/m_q_struct for strange — same ratio as the structural
    masses).  Setting use_struct_mass_in_chromo=True uses m_q_struct
    directly in the chromomag denominator (equivalent to the QCD constituent
    quark identification).

    Args:
        spec: HyperonSpinFlavour entry for this baryon.
        m_q_struct_MeV: Light-quark structure mass [MeV].
        m_s_struct_MeV: Strange-quark structure mass [MeV].
        sigma_GeV2: String tension [GeV²].
        xi_inv_GeV: Coherence length [GeV⁻¹].
        use_struct_mass_in_chromo: If True, use the structure masses
            directly in the chromomag denominator (QCD-style).

    Returns:
        Baryon mass in MeV.
    """
    K_GeV3 = K_substrate_GeV3(sigma_GeV2, xi_inv_GeV)
    K_MeV3 = K_GeV3 * (1000.0 ** 3)

    if use_struct_mass_in_chromo:
        m_q_chromo = m_q_struct_MeV
        m_s_chromo = m_s_struct_MeV
    else:
        m_q_chromo = 1000.0 * M_K_CHROMO_GEV  # = 1000 √σ ≈ 424 MeV
        # Same ratio s/q as the structure masses, anchored at √σ for light
        m_s_chromo = m_q_chromo * (m_s_struct_MeV / m_q_struct_MeV)

    chromo_qq = spec.c_qq / (m_q_chromo ** 2)
    chromo_qs = spec.c_qs / (m_q_chromo * m_s_chromo)
    chromo_ss = spec.c_ss / (m_s_chromo ** 2)
    delta_E_MeV = K_MeV3 * (chromo_qq + chromo_qs + chromo_ss)

    return (
        spec.n_light * m_q_struct_MeV
        + spec.n_strange * m_s_struct_MeV
        + delta_E_MeV
    )


# ---------------------------------------------------------------------------
# Anchor m_s from one hyperon, predict the rest
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HyperonPrediction:
    """One predicted/observed pair for the hyperon spectrum.

    Attributes:
        name:           Baryon name.
        mass_pred_MeV:  Substrate prediction.
        mass_pdg_MeV:   PDG empirical mass.
        frac_err:       (predicted − observed) / observed.
        anchor:         True if this baryon was used to anchor m_s.
    """
    name: str
    mass_pred_MeV: float
    mass_pdg_MeV: float
    frac_err: float
    anchor: bool = False


def solve_m_s_struct_from_anchor_MeV(
    anchor_name: str,
    m_q_struct_MeV: float,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
    use_struct_mass_in_chromo: bool = False,
) -> float:
    """Solve m_s_struct so that the named hyperon hits its PDG mass.

    This is a 1-D root find (the chromomag part depends on m_s through
    1/(m_q m_s) and 1/m_s², which makes the equation nonlinear).  We use
    a fast Newton iteration — typically 3-5 iterations are enough.

    Args:
        anchor_name: Name from HYPERON_TABLE used as the m_s anchor.
        m_q_struct_MeV: Light-quark structure mass [MeV].
        sigma_GeV2: String tension [GeV²].
        xi_inv_GeV: Coherence length [GeV⁻¹].
        use_struct_mass_in_chromo: Pass-through to baryon_mass_MeV.

    Returns:
        m_s_struct in MeV.
    """
    spec = next(h for h in HYPERON_TABLE if h.name == anchor_name)
    target = spec.mass_pdg_MeV

    # Bracket: m_s should be 50-200 MeV heavier than m_q
    m_lo = m_q_struct_MeV
    m_hi = m_q_struct_MeV + 800.0
    for _ in range(80):
        m_mid = 0.5 * (m_lo + m_hi)
        m_pred = baryon_mass_MeV(
            spec, m_q_struct_MeV, m_mid,
            sigma_GeV2=sigma_GeV2, xi_inv_GeV=xi_inv_GeV,
            use_struct_mass_in_chromo=use_struct_mass_in_chromo,
        )
        if m_pred > target:
            m_hi = m_mid
        else:
            m_lo = m_mid
        if abs(m_hi - m_lo) < 1e-6:
            break
    return 0.5 * (m_lo + m_hi)


@dataclass
class HyperonSpectrumResult:
    """Aggregate substrate prediction of the J=1/2 octet + Ω⁻.

    Attributes:
        sigma_GeV2:           String tension used.
        xi_inv_GeV:           Coherence length used.
        K_substrate_GeV3:     Chromomag contact coefficient.
        chromo_pair_MeV:      K / m_q² (per-pair light coupling).
        m_q_struct_MeV:       Light-quark structure mass.
        m_s_struct_MeV:       Strange-quark structure mass.
        delta_s_MeV:          m_s − m_q (the only flavour input).
        anchor_name:          Hyperon used to fix m_s.
        predictions:          List of HyperonPrediction (one per baryon).
        max_abs_err:          Max |frac_err| over non-anchor baryons.
    """
    sigma_GeV2: float
    xi_inv_GeV: float
    K_substrate_GeV3: float
    chromo_pair_MeV: float
    m_q_struct_MeV: float
    m_s_struct_MeV: float
    delta_s_MeV: float
    anchor_name: str
    predictions: list[HyperonPrediction]
    max_abs_err: float


def compute_hyperon_spectrum(
    anchor_name: str = "Λ",
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
    use_struct_mass_in_chromo: bool = False,
) -> HyperonSpectrumResult:
    """Full substrate prediction of N, Λ, Σ, Ξ, Ω⁻, Δ.

    Step 1: solve m_q_struct from the proton anchor (this is fixed by the
            absolute baryon-mass scale; the substrate gives ≈ 365 MeV).
    Step 2: solve m_s_struct so the chosen hyperon (default Λ) hits PDG.
    Step 3: predict the others from the same m_s_struct.

    Args:
        anchor_name: Hyperon used to fix m_s.  Must be in HYPERON_TABLE.
        sigma_GeV2: String tension [GeV²].
        xi_inv_GeV: Coherence length [GeV⁻¹].
        use_struct_mass_in_chromo: If True, evaluate chromomag at the
            structure masses (QCD-style).  Default False uses the
            geometric chromomag mass m_K_chromo = √σ.

    Returns:
        HyperonSpectrumResult with the full per-baryon table.
    """
    K_GeV3 = K_substrate_GeV3(sigma_GeV2, xi_inv_GeV)
    chromo_pair_MeV = chromomag_pair_coupling_MeV(sigma_GeV2, xi_inv_GeV)
    m_q_struct = solve_m_q_struct_MeV(sigma_GeV2, xi_inv_GeV)
    m_s_struct = solve_m_s_struct_from_anchor_MeV(
        anchor_name, m_q_struct,
        sigma_GeV2=sigma_GeV2, xi_inv_GeV=xi_inv_GeV,
        use_struct_mass_in_chromo=use_struct_mass_in_chromo,
    )

    predictions: list[HyperonPrediction] = []
    for h in HYPERON_TABLE:
        m_pred = baryon_mass_MeV(
            h, m_q_struct, m_s_struct,
            sigma_GeV2=sigma_GeV2, xi_inv_GeV=xi_inv_GeV,
            use_struct_mass_in_chromo=use_struct_mass_in_chromo,
        )
        frac_err = (m_pred - h.mass_pdg_MeV) / h.mass_pdg_MeV
        predictions.append(
            HyperonPrediction(
                name=h.name,
                mass_pred_MeV=m_pred,
                mass_pdg_MeV=h.mass_pdg_MeV,
                frac_err=frac_err,
                anchor=(h.name == anchor_name),
            )
        )

    non_anchor_errs = [
        abs(p.frac_err) for p in predictions
        if not p.anchor and p.name != "N(939)"
    ]
    max_abs_err = max(non_anchor_errs) if non_anchor_errs else 0.0

    return HyperonSpectrumResult(
        sigma_GeV2=sigma_GeV2,
        xi_inv_GeV=xi_inv_GeV,
        K_substrate_GeV3=K_GeV3,
        chromo_pair_MeV=chromo_pair_MeV,
        m_q_struct_MeV=m_q_struct,
        m_s_struct_MeV=m_s_struct,
        delta_s_MeV=m_s_struct - m_q_struct,
        anchor_name=anchor_name,
        predictions=predictions,
        max_abs_err=max_abs_err,
    )


# ---------------------------------------------------------------------------
# Sensitivity sweep
# ---------------------------------------------------------------------------

def sensitivity_sweep() -> list[tuple[str, HyperonSpectrumResult]]:
    """Sweep σ, ξ, anchor choice, and chromomag-mass convention.

    Returns:
        List of (label, HyperonSpectrumResult) tuples.
    """
    runs: list[tuple[str, HyperonSpectrumResult]] = []
    runs.append((
        "default (σ=0.180 GeV², ξ=0.2 fm, anchor=Λ, chromo m=√σ)",
        compute_hyperon_spectrum(),
    ))
    runs.append((
        "anchor = Σ instead of Λ",
        compute_hyperon_spectrum(anchor_name="Σ"),
    ))
    runs.append((
        "anchor = Ξ instead of Λ",
        compute_hyperon_spectrum(anchor_name="Ξ"),
    ))
    runs.append((
        "σ × 1.05 (lattice +5%)",
        compute_hyperon_spectrum(sigma_GeV2=0.180 * 1.05),
    ))
    runs.append((
        "σ × 0.95 (lattice −5%)",
        compute_hyperon_spectrum(sigma_GeV2=0.180 * 0.95),
    ))
    for xi_fm in (0.18, 0.22):
        runs.append((
            f"ξ_QCD = {xi_fm} fm",
            compute_hyperon_spectrum(xi_inv_GeV=xi_fm / HBARC_GEVFM),
        ))
    runs.append((
        "chromo m = m_struct (QCD-style propagator)",
        compute_hyperon_spectrum(use_struct_mass_in_chromo=True),
    ))
    return runs


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

def format_result(result: HyperonSpectrumResult) -> str:
    """Multi-line report for one HyperonSpectrumResult."""
    lines = [
        f"  σ                    = {result.sigma_GeV2:.4f} GeV²",
        f"  ξ⁻¹                  = {result.xi_inv_GeV:.4f} GeV⁻¹  "
        f"(σ ξ² = {result.sigma_GeV2 * result.xi_inv_GeV ** 2:.4f})",
        f"  K_substrate          = {result.K_substrate_GeV3:.4e} GeV³",
        f"  chromo pair coupling = {result.chromo_pair_MeV:.2f} MeV  "
        f"(= (8/3) ξ² σ^{{1/2}})",
        f"  m_q_struct           = {result.m_q_struct_MeV:.2f} MeV "
        f"(from p anchor)",
        f"  m_s_struct           = {result.m_s_struct_MeV:.2f} MeV "
        f"(from {result.anchor_name} anchor)",
        f"  Δ_s = m_s − m_q      = {result.delta_s_MeV:+.2f} MeV  "
        f"(only flavour input)",
        f"  max |err| (non-anchor) = {100.0 * result.max_abs_err:+.2f}%",
        "",
        "  Per-baryon predictions:",
        f"    {'baryon':<8} {'pred [MeV]':>11} {'PDG [MeV]':>11}"
        f" {'err':>9}",
    ]
    for p in result.predictions:
        tag = "  ← anchor" if p.anchor else ""
        if p.name == "N(939)":
            tag = "  (proton-anchor)"
        lines.append(
            f"    {p.name:<8} {p.mass_pred_MeV:>11.2f} {p.mass_pdg_MeV:>11.2f}"
            f" {100.0 * p.frac_err:>+8.2f}%{tag}"
        )
    return "\n".join(lines)


def print_full_report() -> None:
    """Print the default prediction plus the sensitivity sweep."""
    print("=" * 80)
    print(" Substrate prediction of the hyperon mass spectrum (§18.64)")
    print(" (J=1/2 octet + J=3/2 Ω⁻; flavour-dependent chromomagnetic spin-spin)")
    print("=" * 80)
    runs = sensitivity_sweep()
    default_label, default_result = runs[0]

    print()
    print(f"DEFAULT — {default_label}")
    print(format_result(default_result))
    print()

    print("-" * 80)
    print(" Sensitivity sweep — max |err| over predicted (non-anchor) baryons")
    print("-" * 80)
    print(f"  {'variant':<48}{'m_s [MeV]':>13}{'max |err|':>13}")
    for label, res in runs:
        print(
            f"  {label:<48}{res.m_s_struct_MeV:>13.2f}"
            f"{100.0 * res.max_abs_err:>+12.2f}%"
        )
    print()

    print("-" * 80)
    print(" Interpretation")
    print("-" * 80)
    print("  * The ONLY flavour input is the strange-quark structure mass m_s,")
    print("    fixed by ONE hyperon anchor (default Λ at 1115.68 MeV).  All other")
    print("    hyperon masses are PREDICTED from the same substrate K_subst, σ, ξ.")
    print("  * The Σ-Λ splitting (≈75 MeV empirical, both uds with J=1/2) is a")
    print("    pure prediction — it comes from the spin-coefficient difference")
    print("    c_qq=-3/4 (Λ, light singlet) vs c_qq=+1/4 with c_qs=-1 (Σ, light")
    print("    triplet).  No extra parameter.")
    print("  * The Ω⁻ mass (sss, J=3/2) is the cleanest prediction: pure strange")
    print("    triplet with no flavour-mixed pairs and no light contribution.")


if __name__ == "__main__":
    print_full_report()
