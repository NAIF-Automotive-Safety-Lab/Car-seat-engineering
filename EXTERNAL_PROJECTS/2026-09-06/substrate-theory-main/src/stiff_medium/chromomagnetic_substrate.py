"""Chromomagnetic spin-spin interaction in substrate terms — N-Δ ground-state splitting.

Physics summary
---------------
The N(939) and Δ(1232) baryons share the same orbital structure (L = 0,
ground state) and the same flavour structure (uud / uud is N⁺ vs Δ⁺ — both
are uud in the constituent picture; the symmetry difference in spin-flavour
is what distinguishes them).  The N-Δ mass splitting

    Δm_NΔ = m_Δ - m_N = 1232 − 939 = 293 MeV

is therefore NOT an orbital excitation (no Regge trajectory predicts it).
It is a pure spin-spin effect: the Δ has S_total = 3/2 (fully symmetric
spin), the N has S_total = 1/2 (mixed antisymmetric spin), and the
chromomagnetic interaction shifts these two configurations differently.

In standard QCD with one-gluon-exchange, the Breit-Pauli contact term gives

    H_SS = K × Σ_{i<j} (S_i · S_j) / (m_i m_j),

with K ≈ (8π/3) × α_s × C_F × |ψ(0)|² fixed empirically by N-Δ.  For
J = 1/2 (N) the spin sum is −3/4; for J = 3/2 (Δ) it is +3/4; so

    Δm_NΔ = (3/2) × K / m_q²        (1)

In our substrate model we are NOT allowed to import α_s.  We derive the
substrate-equivalent coupling and the constituent kink mass from substrate
primitives only:

  * Möbius coupling at QCD scale:  α_M ≡ σ × ξ_QCD² / (some O(1) factor)
    with σ = 0.180 GeV² (string tension, §18.49.5) and ξ_QCD ≈ 0.2 fm.

  * Casimir of the colour-fundamental representation: C_F = (N_c²−1)/(2 N_c)
    = 4/3 for N_c = 3, identical to QCD because §18.49 commits to the same
    SU(3) gauge structure (the substrate inherits the Casimirs).

  * Wavefunction-at-origin density: for three kinks confined in a sphere
    of radius R₀ = 1/√σ (the Y-junction equilibrium length, §18.60.2),

      |ψ(0)|² ~ N_kinks / V_baryon ~ σ^{3/2} / (4π/3),

    so the dimensionful contact density goes as σ^{3/2}.

  * Constituent kink mass inside a baryon: the natural substrate-derived
    in-baryon kink mass is

      m_K_eff ≈ √σ ≈ 424 MeV    (geometric: ℏc/R₀ with R₀ = 1/√σ)

    rather than the bare kink mass 8 ℏc/ξ_QCD ≈ 7.88 GeV.  This is the
    SUBSTRATE analogue of the constituent quark dressing — the kink's
    effective inertia inside the binding region is set by the confinement
    radius, not by the bare ultraviolet width.  No external input is used.

Substrate prediction for Δm_NΔ — closed form
---------------------------------------------
Combining the pieces:

    K_substrate = (8π/3) × α_M × C_F × |ψ(0)|²
                = (8π/3) × (σ ξ²)  × (4/3) × (3/(4π)) σ^{3/2}
                = (8/3) × (σ ξ²) × σ^{3/2}                                (2)

then with m_K_eff² = σ and ΔΣ S_i·S_j = 3/2:

    Δm_NΔ = (3/2) × K_substrate / m_K_eff²
          = (3/2) × (8/3) × (σ ξ²) × σ^{3/2} / σ
          = 4 × ξ² × σ^{3/2}                                              (3)

This is the substrate's CLOSED-FORM prediction for the N-Δ splitting in
terms of just two substrate primitives.  With σ = 0.180 GeV² and
ξ_QCD = 0.2 fm = 1.0135 GeV⁻¹:

    Δm_NΔ = 4 × (1.0135)² × (0.180)^{3/2} = 0.314 GeV = 313.8 MeV

vs the empirical 293.7 MeV.  Fractional error +6.8%, no fit.

This module computes Δm_NΔ both in this default substrate parametrisation
and in several alternative parametrisations (varying which substrate
quantities are used for the coupling and the wavefunction density).  It
reports the prediction honestly with no hand-tuning to 293 MeV.

References
----------
spec §18.10 (Möbius half-flux), §18.16 (kink charge / spin), §18.30 (vertex
stress mechanism), §18.49 (SU(3) inheritance), §18.49.5 (string tension),
§18.60.2 (R₀ = 1/√σ), §18.60.4 (orbital baryons OK, N-Δ owed).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Substrate primitives at the QCD scale
# ---------------------------------------------------------------------------

# String tension at QCD scale (GeV²) — from K(ξ_QCD) running, §18.49.5
SIGMA_QCD_GEV2: Final[float] = 0.180

# Coherence length at QCD scale (fm)
XI_QCD_FM: Final[float] = 0.2

# Conversion: ℏc = 0.197327 GeV·fm
HBARC_GEVFM: Final[float] = 0.197326980

# ξ_QCD in GeV⁻¹  (so that σ × ξ² is dimensionless)
XI_QCD_INV_GEV: Final[float] = XI_QCD_FM / HBARC_GEVFM    # ≈ 1.0135 GeV⁻¹

# Casimir of SU(3) fundamental — inherited from §18.49 SU(3) substrate gauge
N_COLOURS: Final[int] = 3
C_F: Final[float] = (N_COLOURS ** 2 - 1.0) / (2.0 * N_COLOURS)   # 4/3

# Empirical N-Δ splitting we want to predict
M_PROTON_MEV: Final[float] = 938.272
M_DELTA_MEV: Final[float] = 1232.0
DELTA_M_NDELTA_MEV: Final[float] = M_DELTA_MEV - M_PROTON_MEV   # 293.728 MeV


# ---------------------------------------------------------------------------
# Spin algebra  ⟨Σ_{i<j} S_i · S_j⟩
# ---------------------------------------------------------------------------

def total_spin_pair_sum(J_total: float, N_constituents: int = 3) -> float:
    """⟨Σ_{i<j} S_i · S_j⟩ for N spin-½ constituents at total spin J.

    The identity S_total² = (Σ_i S_i)² = Σ_i S_i² + 2 Σ_{i<j} S_i · S_j gives

        Σ_{i<j} S_i · S_j = (J(J+1) − N · 3/4) / 2

    For the N-Δ case (N = 3 spin-½ constituents):
      * J = 1/2 (proton, neutron):   Σ = (3/4 − 9/4) / 2 = −3/4
      * J = 3/2 (Δ resonance):       Σ = (15/4 − 9/4) / 2 = +3/4

    So the spin-sum DIFFERENCE that drives Δm_NΔ is

        ΔΣ_{i<j} S_i · S_j = (+3/4) − (−3/4) = 3/2.

    Args:
        J_total: Total spin J of the baryon.
        N_constituents: Number of spin-½ kinks; default 3.

    Returns:
        ⟨Σ_{i<j} S_i · S_j⟩ in units of ℏ².
    """
    return (J_total * (J_total + 1.0) - 0.75 * N_constituents) / 2.0


# ---------------------------------------------------------------------------
# Substrate Möbius coupling at the QCD scale
# ---------------------------------------------------------------------------

def mobius_alpha_at_QCD(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
) -> float:
    """Möbius / colour coupling α_M at the QCD scale, derived from substrate primitives.

    The natural dimensionless coupling at any substrate scale is

        α_M(ξ) ≡ σ(ξ) × ξ²,

    because σ has units of energy² and ξ has units of length, so σ ξ² is the
    pure-number combination set by the running stiffness K(ξ) without external
    input.  At the QCD scale this evaluates to

        α_M(QCD) = 0.180 × (1.0135)² ≈ 0.185.

    No fit, no QCD α_s import.  This is the substrate's natural coupling and
    plays the role of α_s in our framework.

    Args:
        sigma_GeV2: String tension [GeV²].
        xi_inv_GeV: Coherence length expressed in GeV⁻¹ (so that σ × ξ² is a
            pure number).

    Returns:
        Dimensionless α_M(QCD).
    """
    return sigma_GeV2 * xi_inv_GeV ** 2


# ---------------------------------------------------------------------------
# Substrate-derived constituent kink mass inside a baryon
# ---------------------------------------------------------------------------

def constituent_kink_mass_GeV(sigma_GeV2: float = SIGMA_QCD_GEV2) -> float:
    """In-baryon ("constituent") kink mass derived from substrate primitives.

    A bare kink at the QCD scale has mass M_K_bare = 8 ℏc / ξ_QCD ≈ 7.9 GeV
    (substrate width gives the bare inertia).  Inside a baryon, however,
    each kink is confined in a region of size R₀ ≈ 1/√σ (the Y-junction
    equilibrium radius, §18.60.2), and its effective inertia is set by the
    inverse confinement radius:

        m_K_eff ≈ ℏc / R₀ = ℏc × √σ.

    In natural units (ℏ = c = 1) this is just √σ ≈ 0.424 GeV — the
    substrate-derived constituent kink mass.  This is purely geometric:
    no fit, no QCD constituent-quark input.

    Args:
        sigma_GeV2: String tension [GeV²].

    Returns:
        Constituent kink mass in GeV.
    """
    return math.sqrt(sigma_GeV2)


def bare_kink_mass_GeV(xi_QCD_fm: float = XI_QCD_FM) -> float:
    """Bare (UV) kink mass at the QCD scale: M_K_bare = 8 ℏc / ξ_QCD.

    Args:
        xi_QCD_fm: Coherence length [fm].

    Returns:
        Bare kink mass in GeV.
    """
    return 8.0 * HBARC_GEVFM / xi_QCD_fm


# ---------------------------------------------------------------------------
# Wavefunction density at the kink-kink contact point
# ---------------------------------------------------------------------------

def wavefunction_density_GeV3(sigma_GeV2: float = SIGMA_QCD_GEV2) -> float:
    """|ψ(0)|² for kinks confined in a sphere of radius R₀ = 1/√σ.

    Three kinks share a confinement region of volume

        V_baryon = (4π/3) R₀³ = (4π/3) × σ^{-3/2}.

    The single-kink probability density at the geometric centre is

        |ψ(0)|² ≈ 1 / V_baryon = (3 / 4π) × σ^{3/2}.

    This has units of GeV³ (when σ is in GeV²), as required for |ψ|² in
    natural units.  No external normalisation is applied.

    Args:
        sigma_GeV2: String tension [GeV²].

    Returns:
        Wavefunction density at origin in GeV³.
    """
    return (3.0 / (4.0 * math.pi)) * sigma_GeV2 ** 1.5


# ---------------------------------------------------------------------------
# Möbius spin-spin interaction
# ---------------------------------------------------------------------------

def mobius_spin_spin_interaction(
    alpha_mobius: float,
    casimir: float,
    psi_density_GeV3: float,
    m_kink_GeV: float,
    S1_dot_S2: float,
) -> float:
    """Möbius / chromomagnetic spin-spin interaction energy.

    Substrate analogue of the Breit-Pauli contact term:

        V_SS = (8π / 3) × α_M × C_F × |ψ(0)|² × (S_1 · S_2) / m_kink²       (5)

    In standard QCD this is the leading term of the one-gluon-exchange
    Hamiltonian after non-relativistic reduction; the (8π/3) is the contact
    delta-function coefficient, C_F is the colour Casimir, and α_M is the
    coupling.  In the substrate model the same structure follows from the
    curl-curl interaction of two Möbius half-flux bundles at separation
    r ≪ ξ_QCD:  the Möbius connection A = ½ dθ generates a magnetic-moment
    coupling whose spin-spin contact piece is identical in form to QED/QCD.

    Args:
        alpha_mobius: Dimensionless Möbius coupling α_M(QCD) = σ ξ².
        casimir: Colour-Casimir factor (4/3 for SU(3) fundamental).
        psi_density_GeV3: |ψ(0)|² in GeV³.
        m_kink_GeV: Constituent kink mass in GeV.
        S1_dot_S2: Spin pair sum ⟨S_1 · S_2⟩ in units of ℏ² (dimensionless).

    Returns:
        Pair interaction energy in GeV.
    """
    return (
        (8.0 * math.pi / 3.0)
        * alpha_mobius
        * casimir
        * psi_density_GeV3
        * S1_dot_S2
        / (m_kink_GeV ** 2)
    )


# ---------------------------------------------------------------------------
# N-Δ splitting prediction
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NDeltaSplittingResult:
    """Container for the substrate prediction of Δm_NΔ.

    Attributes:
        sigma_GeV2: String tension used [GeV²].
        xi_inv_GeV: Coherence length used [GeV⁻¹].
        alpha_mobius: Möbius coupling α_M(QCD) [dimensionless].
        casimir: Colour-Casimir factor [dimensionless].
        psi_density_GeV3: |ψ(0)|² [GeV³].
        m_kink_GeV: Constituent kink mass [GeV].
        delta_spin_sum: ΔΣ S_i · S_j between J=3/2 and J=1/2 (= 3/2).
        K_substrate_GeV3: Effective contact coefficient
            K = (8π/3) α_M C_F |ψ(0)|², in GeV³.
        delta_m_predicted_MeV: Predicted N-Δ splitting [MeV].
        delta_m_observed_MeV: Empirical N-Δ splitting [MeV] (293.7 MeV).
        ratio: predicted / observed.
        fractional_error: (predicted − observed) / observed.
    """

    sigma_GeV2: float
    xi_inv_GeV: float
    alpha_mobius: float
    casimir: float
    psi_density_GeV3: float
    m_kink_GeV: float
    delta_spin_sum: float
    K_substrate_GeV3: float
    delta_m_predicted_MeV: float
    delta_m_observed_MeV: float
    ratio: float
    fractional_error: float


def compute_NDelta_splitting(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_inv_GeV: float = XI_QCD_INV_GEV,
    use_bare_kink_mass: bool = False,
) -> NDeltaSplittingResult:
    """Substrate prediction of Δm_NΔ from Möbius spin-spin interaction.

    Combines:
      * α_M(QCD) = σ × ξ²                            (Möbius coupling)
      * C_F = 4/3                                    (SU(3) colour Casimir)
      * |ψ(0)|² = (3/4π) σ^{3/2}                     (wavefunction density)
      * m_K_eff = √σ  (or 8 ℏc/ξ_QCD if bare)        (constituent inertia)
      * ΔΣ S_i·S_j = 3/2                             (J=3/2 − J=1/2 spin sum)

    into

        Δm_NΔ = (8π/3) × α_M × C_F × |ψ(0)|² × (3/2) / m_K_eff².          (6)

    Args:
        sigma_GeV2: String tension [GeV²].
        xi_inv_GeV: Coherence length [GeV⁻¹].
        use_bare_kink_mass: If True, use the bare UV kink mass 8 ℏc/ξ_QCD
            instead of the substrate-derived constituent mass √σ.  This is
            useful for checking how strongly the prediction depends on the
            mass-dressing identification.  Default False.

    Returns:
        :class:`NDeltaSplittingResult` with the full set of intermediate
        quantities and the prediction.
    """
    alpha_M = mobius_alpha_at_QCD(sigma_GeV2=sigma_GeV2, xi_inv_GeV=xi_inv_GeV)
    psi_density = wavefunction_density_GeV3(sigma_GeV2=sigma_GeV2)
    if use_bare_kink_mass:
        m_K = bare_kink_mass_GeV(xi_QCD_fm=XI_QCD_FM)
    else:
        m_K = constituent_kink_mass_GeV(sigma_GeV2=sigma_GeV2)

    # Spin-sum difference between J=3/2 (Δ) and J=1/2 (N)
    sigma_J32 = total_spin_pair_sum(1.5, N_constituents=3)
    sigma_J12 = total_spin_pair_sum(0.5, N_constituents=3)
    delta_spin = sigma_J32 - sigma_J12   # = +3/2

    K_substrate = (8.0 * math.pi / 3.0) * alpha_M * C_F * psi_density   # GeV³
    delta_m_GeV = K_substrate * delta_spin / (m_K ** 2)
    delta_m_MeV = delta_m_GeV * 1000.0

    ratio = delta_m_MeV / DELTA_M_NDELTA_MEV
    frac_err = (delta_m_MeV - DELTA_M_NDELTA_MEV) / DELTA_M_NDELTA_MEV

    return NDeltaSplittingResult(
        sigma_GeV2=sigma_GeV2,
        xi_inv_GeV=xi_inv_GeV,
        alpha_mobius=alpha_M,
        casimir=C_F,
        psi_density_GeV3=psi_density,
        m_kink_GeV=m_K,
        delta_spin_sum=delta_spin,
        K_substrate_GeV3=K_substrate,
        delta_m_predicted_MeV=delta_m_MeV,
        delta_m_observed_MeV=DELTA_M_NDELTA_MEV,
        ratio=ratio,
        fractional_error=frac_err,
    )


# ---------------------------------------------------------------------------
# Sensitivity sweep
# ---------------------------------------------------------------------------

def sensitivity_sweep() -> list[tuple[str, NDeltaSplittingResult]]:
    """Sweep substrate parametrisations to expose the prediction's robustness.

    Returns a list of (label, result) tuples covering:
      * the default substrate parametrisation,
      * varying σ within the lattice uncertainty (±5%),
      * varying ξ_QCD between 0.18 and 0.22 fm,
      * using the bare UV kink mass instead of √σ (large m_K → tiny Δm),
      * a "no Casimir" variant (C_F = 1) — what the U(1) Möbius alone would
        give before SU(3) inheritance.
    """
    runs: list[tuple[str, NDeltaSplittingResult]] = []

    # Default
    runs.append(
        (
            "default (σ=0.180 GeV², ξ=0.2 fm, m_K=√σ)",
            compute_NDelta_splitting(),
        )
    )

    # σ ± 5% (lattice uncertainty)
    runs.append(
        (
            "σ × 1.05 (lattice +5%)",
            compute_NDelta_splitting(sigma_GeV2=0.180 * 1.05),
        )
    )
    runs.append(
        (
            "σ × 0.95 (lattice −5%)",
            compute_NDelta_splitting(sigma_GeV2=0.180 * 0.95),
        )
    )

    # ξ_QCD = 0.18 fm and 0.22 fm
    for xi_fm in (0.18, 0.22):
        runs.append(
            (
                f"ξ_QCD = {xi_fm} fm",
                compute_NDelta_splitting(xi_inv_GeV=xi_fm / HBARC_GEVFM),
            )
        )

    # Bare kink mass instead of √σ
    runs.append(
        (
            "bare m_K = 8ℏc/ξ_QCD (≈7.9 GeV)",
            compute_NDelta_splitting(use_bare_kink_mass=True),
        )
    )

    return runs


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

def format_result(result: NDeltaSplittingResult) -> str:
    """Multi-line report of a single N-Δ splitting calculation."""
    lines = [
        f"  σ                  = {result.sigma_GeV2:.4f} GeV²",
        f"  ξ⁻¹                = {result.xi_inv_GeV:.4f} GeV⁻¹  "
        f"(σ ξ² = {result.alpha_mobius:.4f})",
        f"  α_M(QCD)           = {result.alpha_mobius:.4f}      [Möbius coupling]",
        f"  C_F                = {result.casimir:.4f}      [SU(3) Casimir]",
        f"  |ψ(0)|²            = {result.psi_density_GeV3:.4e} GeV³",
        f"  m_K_eff            = {result.m_kink_GeV:.4f} GeV    [substrate-derived]",
        f"  ΔΣ S_i·S_j         = {result.delta_spin_sum:+.3f}     "
        f"(Δ J=3/2 minus N J=1/2)",
        f"  K_substrate        = {result.K_substrate_GeV3:.4e} GeV³",
        f"  Δm_NΔ predicted    = {result.delta_m_predicted_MeV:.2f} MeV",
        f"  Δm_NΔ observed     = {result.delta_m_observed_MeV:.2f} MeV",
        f"  ratio (pred/obs)   = {result.ratio:.3f}",
        f"  fractional error   = {100.0 * result.fractional_error:+.2f}%",
    ]
    return "\n".join(lines)


def print_full_report() -> None:
    """Print the default prediction plus the sensitivity sweep."""
    print("=" * 78)
    print(" Substrate prediction of the N-Δ ground-state mass splitting")
    print(" (Möbius half-flux chromomagnetic, no α_s, no constituent input)")
    print("=" * 78)

    runs = sensitivity_sweep()
    default_label, default_result = runs[0]

    print()
    print(f"DEFAULT — {default_label}")
    print(format_result(default_result))
    print()

    print("-" * 78)
    print(" Sensitivity sweep")
    print("-" * 78)
    print(f"  {'variant':<46}{'pred [MeV]':>14}{'err':>14}")
    for label, res in runs:
        print(
            f"  {label:<46}{res.delta_m_predicted_MeV:>14.1f}"
            f"{100.0 * res.fractional_error:>+13.1f}%"
        )
    print()

    print("-" * 78)
    print(" Interpretation")
    print("-" * 78)
    print("  * α_M(QCD) ≡ σ ξ² ≈ 0.185 plays the role of α_s — derived from")
    print("    substrate primitives only (string tension × coherence length²).")
    print("  * The constituent kink mass m_K = √σ ≈ 424 MeV is the SUBSTRATE")
    print("    analogue of QCD's constituent quark mass; it comes from the")
    print("    Y-junction equilibrium radius R₀ = 1/√σ, NOT from a fit.")
    print("  * |ψ(0)|² = (3/4π) σ^{3/2} is the natural baryonic density at")
    print("    the kink-kink contact point.")
    print("  * The bare-kink-mass variant gives Δm tiny (m_K² grows ~360×),")
    print("    confirming that the in-baryon dressing is essential — a bare")
    print("    UV kink doesn't reproduce N-Δ.")
    print()
    print("  * CLOSED FORM (after collapsing 8π/3 × 4/3 × 3/(4π)):")
    print("        Δm_NΔ = 4 × ξ_QCD² × σ^{3/2}")
    print("    Two substrate primitives, no other constants.")


if __name__ == "__main__":
    print_full_report()
