"""Proton charge radius from the substrate Y-junction picture.

Physics summary
---------------
In the substrate model the proton is a Y-junction of three quark kinks bound
at a vertex by three strings of tension

    σ = σ_QCD = 0.180 GeV²        (from the K(ξ) running, see regge_spectrum.py)

Each kink has an intrinsic transverse width set by the QCD coherence scale

    ξ_QCD ≈ 0.2 fm.

The proton's mean-square charge radius receives two physically distinct
contributions:

    ⟨r²⟩_p = ⟨r²⟩_position + ⟨r²⟩_intrinsic

where ⟨r²⟩_position is the spread of the three quark CENTERS about the
Y-junction, and ⟨r²⟩_intrinsic is the rms charge spread WITHIN each kink.
We evaluate both from first principles (no hand-tuning).

The position contribution is set by the Y-junction equilibrium length R₀.
Two physically motivated estimates of R₀ exist; we report both:

    (A) Dimensional / variational: a relativistic kink in a linear
        confining potential V(r) = σ r has ground-state ⟨r⟩ ≈ 1/√σ.
        This is the natural length scale for a single string of tension σ.

    (B) Rotating Y-junction: for the proton mass M_p anchored at 0.938 GeV
        and three Nambu-Goto strings rotating at the speed of light,
            M_p = (3π/2) σ L  ⇒  L = 2 M_p / (3π σ)
        L is the length from vertex (center) to each quark.

Estimate (A) is parameter-free in our model and is the cleanest "prediction"
from σ alone; estimate (B) uses the empirical proton mass as an extra anchor.

For the intrinsic kink width contribution, a 1D sine-Gordon kink φ(x) =
4 arctan(exp(x/ξ)) has charge density d_x φ ∝ sech²(x/ξ).  The mean-square
spread of this profile is ⟨x²⟩ = (π² / 12) ξ² ≈ 0.822 ξ².  Taking the
3D smearing as three independent transverse copies, ⟨r²⟩_intrinsic =
3 × ⟨x²⟩ for a fully isotropic kink, but the proton kinks lie along the
(approximately planar) Y-arms, so we use the 1D value as a conservative
lower bound and the 3D value as an upper bound.

The total prediction is then:

    r_p² = ⟨r²⟩_position + ⟨r²⟩_intrinsic
         = c_geom × R₀² + c_kink × ξ_QCD²

The geometric factor c_geom is exactly the charge-weighted average of the
quark positions for a uud proton with q_u = +2/3, q_d = -1/3:

    c_geom = (2 q_u + q_d) / Q_p = (4/3 - 1/3) / 1 = 1

(every quark is at the same distance R₀ from the center, so the weighted
mean-square distance is just R₀² for any consistent set of charges).

Honest assessment
-----------------
- σ_QCD = 0.180 GeV² and ξ_QCD = 0.2 fm are taken from the K(ξ) substrate
  running (regge_spectrum.py / substrate_rg_running.py).  They are NOT
  fitted to the charge radius.

- The Y-junction geometry (3 strings meeting at a vertex, 3 quarks at the
  far ends) is the simplest baryon configuration and is what gives the
  correct Regge slope α' = 1/(2πσ) for nucleons.  No charge-radius input
  is used to set this up.

- The geometric factor c_geom = 1 follows from the equilateral triangle
  charge distribution; this is not adjustable.

- The kink-width factor c_kink lies in [π²/12, π²/4] depending on whether
  the kink width adds in 1D, 2D, or 3D.  We report all three.

- We do NOT hand-tune to match 0.84 fm.  The reported numbers come from
  the formulae above with σ_QCD and ξ_QCD inserted; the 5% sensitivity
  scans probe the parametric dependence.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Physical constants and substrate inputs
# ---------------------------------------------------------------------------

# ℏc in MeV·fm — converts inverse-energy lengths to fm.
HBARC_MEV_FM: float = 197.3269804

# ℏc in GeV·fm
HBARC_GEV_FM: float = HBARC_MEV_FM * 1e-3  # = 0.1973269804

# String tension at QCD scale (from K(ξ) running).
SIGMA_QCD_GEV2: float = 0.180  # GeV²

# Coherence length at QCD scale.
XI_QCD_FM: float = 0.2  # fm

# Empirical proton charge radius (CODATA 2018 / muonic hydrogen).
R_PROTON_CHARGE_EXP_FM: float = 0.8409  # ± 0.0004 fm

# Empirical proton mass.
M_PROTON_GEV: float = 0.938272  # GeV


# ---------------------------------------------------------------------------
# Y-junction equilibrium radius
# ---------------------------------------------------------------------------

def equilibrium_radius_dimensional(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
) -> float:
    """Y-junction equilibrium length from dimensional / variational scaling.

    A relativistic massless quark in the linear confining potential
    V(r) = σ r has the natural length scale R₀ = 1/√σ (in natural units).
    This is the only length that can be built out of σ alone, so it sets
    the size of the bound state by dimensional analysis.

    Variationally, a Gaussian trial wave function ψ(r) ∝ exp(-r²/2 a²) in
    the potential V = σ r gives ⟨H⟩ minimized at a ≈ 0.94 / √σ, and
    ⟨r⟩ ≈ √(8/π) a ≈ 1.50 / √σ.  We use the simplest dimensional value
    R₀ = 1/√σ as the cleanest first-principles estimate.

    Args:
        sigma_GeV2: String tension [GeV²].

    Returns:
        R₀ in fm.
    """
    R0_GeV_inv = 1.0 / math.sqrt(sigma_GeV2)  # GeV⁻¹
    return R0_GeV_inv * HBARC_GEV_FM           # fm


def equilibrium_radius_rotating_Y(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    M_proton_GeV: float = M_PROTON_GEV,
) -> float:
    """Y-junction equilibrium length from rotating-Y mass relation.

    For three relativistic Nambu-Goto strings rotating about a common
    vertex with the quark ends at the speed of light, the classical mass
    of each string is M_s = (π/2) σ L, where L is the vertex-to-quark
    length.  Three strings give M_p = (3π/2) σ L, hence

        L = 2 M_p / (3 π σ).

    This is the standard rotating-Y derivation; it uses the empirical
    proton mass as an additional anchor (so this estimate is NOT
    parameter-free in the way the dimensional one is).

    Args:
        sigma_GeV2:    String tension [GeV²].
        M_proton_GeV:  Proton mass [GeV].

    Returns:
        L in fm.
    """
    L_GeV_inv = 2.0 * M_proton_GeV / (3.0 * math.pi * sigma_GeV2)
    return L_GeV_inv * HBARC_GEV_FM


# ---------------------------------------------------------------------------
# Intrinsic kink width contribution to ⟨r²⟩
# ---------------------------------------------------------------------------

def intrinsic_r2_kink_1D(xi_fm: float = XI_QCD_FM) -> float:
    """⟨x²⟩ of the sine-Gordon kink charge profile (1D).

    The kink profile φ(x) = 4 arctan(exp(x/ξ)) has charge density
    ρ(x) = (1/π) sech(x/ξ) (the topological charge density), normalised
    so ∫ ρ dx = 1.  Its mean-square spread is

        ⟨x²⟩ = (π²/4) ξ²  ≈ 2.467 ξ²    (sech profile)

    For the canonical sine-Gordon kink the charge density along the
    string axis is dφ/dx ∝ sech²(x/ξ), normalised as the energy density;
    the corresponding ⟨x²⟩ is

        ⟨x²⟩ = (π²/12) ξ²  ≈ 0.822 ξ²   (sech² profile)

    We use the sech² (energy-density) value as the standard kink width.

    Args:
        xi_fm: Kink width [fm].

    Returns:
        ⟨x²⟩ in fm².
    """
    return (math.pi ** 2 / 12.0) * xi_fm ** 2


def intrinsic_r2_kink_3D(xi_fm: float = XI_QCD_FM) -> float:
    """⟨r²⟩ for a fully isotropic 3D kink charge cloud.

    Treating the kink width as transverse to the string in all three
    directions gives ⟨r²⟩ = ⟨x²⟩ + ⟨y²⟩ + ⟨z²⟩ = 3 ⟨x²⟩.

    Args:
        xi_fm: Kink width [fm].

    Returns:
        ⟨r²⟩ in fm².
    """
    return 3.0 * intrinsic_r2_kink_1D(xi_fm)


# ---------------------------------------------------------------------------
# Charge radius from Y-junction configuration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ProtonRadiusEstimate:
    """One charge-radius estimate with its physical inputs.

    Attributes:
        label:           Human-readable name of the estimate.
        R0_fm:           Y-junction equilibrium radius [fm].
        r2_position:     ⟨r²⟩_position contribution [fm²].
        r2_intrinsic:    ⟨r²⟩_intrinsic contribution [fm²].
        r_p_fm:          Total predicted r_p [fm].
        frac_error:      (r_p_pred - r_p_exp) / r_p_exp.
    """
    label: str
    R0_fm: float
    r2_position: float
    r2_intrinsic: float
    r_p_fm: float
    frac_error: float


def proton_charge_radius(
    R0_fm: float,
    xi_fm: float = XI_QCD_FM,
    c_geom: float = 1.0,
    kink_dim: int = 3,
    label: str = "estimate",
) -> ProtonRadiusEstimate:
    """Compute the proton charge radius from a given Y-junction radius and kink width.

    Formula:

        r_p² = c_geom × R₀² + ⟨r²⟩_intrinsic(ξ, kink_dim)

    where c_geom = 1 for an equilateral uud Y-junction (charge-weighted
    quark positions all sit at the same distance R₀ from the center).

    Args:
        R0_fm:    Equilibrium quark-to-vertex distance [fm].
        xi_fm:    Kink coherence width [fm].
        c_geom:   Geometric coefficient on R₀² (default 1 for equilateral).
        kink_dim: Dimensionality of the kink width contribution (1 or 3).
        label:    Tag for the estimate.

    Returns:
        ProtonRadiusEstimate with all components and the fractional error
        relative to the empirical r_p^charge = 0.8409 fm.
    """
    if kink_dim == 1:
        r2_intrinsic = intrinsic_r2_kink_1D(xi_fm)
    elif kink_dim == 3:
        r2_intrinsic = intrinsic_r2_kink_3D(xi_fm)
    else:
        raise ValueError(f"kink_dim must be 1 or 3, got {kink_dim}")

    r2_position = c_geom * R0_fm ** 2
    r2_total = r2_position + r2_intrinsic
    r_p_fm = math.sqrt(r2_total)
    frac_error = (r_p_fm - R_PROTON_CHARGE_EXP_FM) / R_PROTON_CHARGE_EXP_FM

    return ProtonRadiusEstimate(
        label=label,
        R0_fm=R0_fm,
        r2_position=r2_position,
        r2_intrinsic=r2_intrinsic,
        r_p_fm=r_p_fm,
        frac_error=frac_error,
    )


# ---------------------------------------------------------------------------
# Sensitivity analysis
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SensitivityScan:
    """Sensitivity of r_p to ±5% variation in σ and ξ.

    Attributes:
        label:           Estimate this scan was run on.
        nominal_r_p:     r_p at central σ, ξ [fm].
        delta_sigma_pct: dr_p/r_p (%) for +5% σ variation.
        delta_xi_pct:    dr_p/r_p (%) for +5% ξ variation.
    """
    label: str
    nominal_r_p: float
    delta_sigma_plus_pct: float
    delta_sigma_minus_pct: float
    delta_xi_plus_pct: float
    delta_xi_minus_pct: float


def sensitivity_scan(
    estimate_fn,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_fm: float = XI_QCD_FM,
    label: str = "scan",
) -> SensitivityScan:
    """Compute fractional shift in r_p for ±5% variations in σ and ξ.

    Args:
        estimate_fn: Callable (sigma, xi) -> ProtonRadiusEstimate.
        sigma_GeV2:  Central value of σ [GeV²].
        xi_fm:       Central value of ξ [fm].
        label:       Tag.

    Returns:
        SensitivityScan with all four perturbation results.
    """
    nominal = estimate_fn(sigma_GeV2, xi_fm).r_p_fm

    def pct(r_p_perturbed: float) -> float:
        return 100.0 * (r_p_perturbed - nominal) / nominal

    return SensitivityScan(
        label=label,
        nominal_r_p=nominal,
        delta_sigma_plus_pct=pct(estimate_fn(sigma_GeV2 * 1.05, xi_fm).r_p_fm),
        delta_sigma_minus_pct=pct(estimate_fn(sigma_GeV2 * 0.95, xi_fm).r_p_fm),
        delta_xi_plus_pct=pct(estimate_fn(sigma_GeV2, xi_fm * 1.05).r_p_fm),
        delta_xi_minus_pct=pct(estimate_fn(sigma_GeV2, xi_fm * 0.95).r_p_fm),
    )


# ---------------------------------------------------------------------------
# Top-level summary
# ---------------------------------------------------------------------------

@dataclass
class ProtonRadiusSummary:
    """All proton-radius estimates and sensitivity scans."""
    sigma_GeV2: float
    xi_fm: float
    estimates: list[ProtonRadiusEstimate]
    sensitivities: list[SensitivityScan]


def compute_proton_radius_summary(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_fm: float = XI_QCD_FM,
) -> ProtonRadiusSummary:
    """Compute all proton-radius estimates from substrate inputs.

    We report four estimates that combine the two R₀ choices (dimensional
    vs rotating-Y) with the two kink-width regimes (1D vs 3D).

    Args:
        sigma_GeV2: String tension [GeV²].
        xi_fm:      Kink coherence width [fm].

    Returns:
        ProtonRadiusSummary with all four estimates and four sensitivity scans.
    """
    R0_dim = equilibrium_radius_dimensional(sigma_GeV2)
    R0_rot = equilibrium_radius_rotating_Y(sigma_GeV2, M_PROTON_GEV)

    estimates = [
        proton_charge_radius(
            R0_dim, xi_fm, kink_dim=3,
            label=f"R0 = 1/√σ (dimensional), 3D kink",
        ),
        proton_charge_radius(
            R0_dim, xi_fm, kink_dim=1,
            label=f"R0 = 1/√σ (dimensional), 1D kink",
        ),
        proton_charge_radius(
            R0_rot, xi_fm, kink_dim=3,
            label=f"R0 = 2 M_p/(3πσ) (rotating Y), 3D kink",
        ),
        proton_charge_radius(
            R0_rot, xi_fm, kink_dim=1,
            label=f"R0 = 2 M_p/(3πσ) (rotating Y), 1D kink",
        ),
    ]

    def make_dim3(s, x):
        return proton_charge_radius(
            equilibrium_radius_dimensional(s), x, kink_dim=3,
            label="dim+3D",
        )

    def make_dim1(s, x):
        return proton_charge_radius(
            equilibrium_radius_dimensional(s), x, kink_dim=1,
            label="dim+1D",
        )

    def make_rot3(s, x):
        return proton_charge_radius(
            equilibrium_radius_rotating_Y(s, M_PROTON_GEV), x, kink_dim=3,
            label="rot+3D",
        )

    def make_rot1(s, x):
        return proton_charge_radius(
            equilibrium_radius_rotating_Y(s, M_PROTON_GEV), x, kink_dim=1,
            label="rot+1D",
        )

    sensitivities = [
        sensitivity_scan(make_dim3, sigma_GeV2, xi_fm, "dimensional + 3D kink"),
        sensitivity_scan(make_dim1, sigma_GeV2, xi_fm, "dimensional + 1D kink"),
        sensitivity_scan(make_rot3, sigma_GeV2, xi_fm, "rotating Y + 3D kink"),
        sensitivity_scan(make_rot1, sigma_GeV2, xi_fm, "rotating Y + 1D kink"),
    ]

    return ProtonRadiusSummary(
        sigma_GeV2=sigma_GeV2,
        xi_fm=xi_fm,
        estimates=estimates,
        sensitivities=sensitivities,
    )
