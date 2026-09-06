"""Proton charge radius v2: add string-fluctuation broadening (§18.60.2 follow-up).

Physics summary
---------------
The bare Y-junction + sine-Gordon kink picture (proton_radius.py, §18.60.2)
predicts r_p ≈ 0.47-0.56 fm, missing the empirical 0.8409 fm by ~33%. The
spec already flagged the missing physics: **transverse string fluctuations
of the Y-arms broaden the effective charge distribution**. This module
computes that contribution from first principles.

For a relativistic Nambu-Goto string of length L with one end fixed at the
Y-junction vertex and the other at a quark, transverse fluctuations of the
string have a logarithmic UV divergence regularised by the kink width
ξ = ξ_QCD. The standard string-theory result for the midpoint of an
open string with both endpoints fixed (the Lüscher formula) is

    ⟨X⊥²(L/2)⟩_per_dim = (1/(2π σ)) × ln(L / ξ_UV)         [midpoint, fixed-fixed]

For a free endpoint (massless quark), the rms displacement is larger by
a factor of 2 from the standard mode-sum (the free-endpoint Dirichlet
boundary doubles the lowest-mode amplitude):

    ⟨X⊥²(L)⟩_per_dim = (1/(π σ)) × ln(L / ξ_UV)            [endpoint, free]

We report BOTH so the reader can see the size of this modelling choice.
The endpoint formula is the physically correct one for the Y-junction with
massless / very-light quarks, but the bare midpoint formula is the
conservative lower bound and matches the prompt's wording.

For three quarks at the ends of three Y-arms, the total contribution to
the charge-radius² is the average of the per-quark fluctuations
(charge-weighted, but since the equilateral uud has c_geom = 1, the
average is just the per-quark value). Two transverse dimensions per
string contribute, so

    ⟨δr²⟩_string_fluct = 2 × ⟨X⊥²(L)⟩_per_dim
                       = (1/(π σ)) × ln(R₀ / ξ_QCD)        [midpoint formula]
                       = (2/(π σ)) × ln(R₀ / ξ_QCD)        [endpoint formula]

Total decomposition
-------------------
The full proton charge-radius² is

    r_p² = R₀²                      [Y-junction equilibrium positions]
         + ⟨δr²⟩_string_fluct       [transverse string fluctuations]
         + ⟨r²⟩_intrinsic_kink      [sine-Gordon sech² profile]
         + ⟨δr²⟩_pion_cloud         [optional: meson cloud, not in bare model]

with σ_QCD = 0.180 GeV², ξ_QCD = 0.2 fm, R₀ = 1/√σ ≈ 0.465 fm.

Pion-cloud caveat
-----------------
In conventional QCD the pion cloud contributes ~+0.05-0.10 fm² to ⟨r_p²⟩
(see e.g. Hammer & Meißner reviews). In the substrate model pions are
kinks of the chiral-condensate field; their CONTRIBUTION to ⟨r_p²⟩ is
dominated by transverse string fluctuations of the Y-arms (which is
already what we are calculating). To avoid double-counting we report
the pion-cloud literature value SEPARATELY, as a cross-check, and not
as an additive term in the substrate prediction.

Honest accounting
-----------------
We report a 4-component decomposition for each estimate:
  • r_p²_position    = R₀²
  • r_p²_string      = ⟨δr²⟩_string_fluct
  • r_p²_kink        = ⟨r²⟩_intrinsic_kink
  • r_p²_pion(lit)   = literature pion-cloud value (cross-check only)

The substrate prediction is r_p² = position + string + kink. The pion-cloud
column lets us see whether adding the literature value would over- or
under-shoot empirical 0.7071 fm² = (0.8409 fm)².

NO inputs are fitted to r_p. σ comes from K(ξ) running at ξ_QCD;
ξ_QCD = 0.2 fm is the substrate coherence length.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from stiff_medium.proton_radius import (
    HBARC_GEV_FM,
    M_PROTON_GEV,
    R_PROTON_CHARGE_EXP_FM,
    SIGMA_QCD_GEV2,
    XI_QCD_FM,
    equilibrium_radius_dimensional,
    equilibrium_radius_rotating_Y,
    intrinsic_r2_kink_1D,
    intrinsic_r2_kink_3D,
)

# ---------------------------------------------------------------------------
# Pion-cloud literature value (cross-check only, not additive in our model)
# ---------------------------------------------------------------------------

# Conventional QCD: pion-cloud contributes Δ⟨r²⟩ ≈ 0.05-0.10 fm² (Hammer &
# Meißner, "Nucleon Form Factors", PPNP 2020). We use the central value.
R2_PION_CLOUD_LIT_FM2: float = 0.075  # fm²


# ---------------------------------------------------------------------------
# String-fluctuation contributions
# ---------------------------------------------------------------------------

def string_fluctuation_per_dim_midpoint(
    sigma_GeV2: float,
    R0_fm: float,
    xi_fm: float,
) -> float:
    """Transverse string fluctuation per dimension at the midpoint (Lüscher).

    Formula (Lüscher / Polyakov, both endpoints fixed):

        ⟨X⊥²(L/2)⟩_per_dim = (1/(2π σ)) × ln(L / ξ_UV)

    With σ in GeV² and lengths in fm we convert via (ℏc)² = 0.0389 GeV²·fm².

    Args:
        sigma_GeV2: String tension [GeV²].
        R0_fm:      String length L = R₀ [fm].
        xi_fm:      UV cutoff ξ = ξ_QCD [fm].

    Returns:
        ⟨X⊥²⟩ in fm² per transverse dimension.
    """
    if R0_fm <= xi_fm:
        return 0.0  # log argument < 1 — no fluctuation enhancement
    log_arg = math.log(R0_fm / xi_fm)
    rms_GeV_inv2 = log_arg / (2.0 * math.pi * sigma_GeV2)  # GeV⁻²
    return rms_GeV_inv2 * HBARC_GEV_FM ** 2                 # fm²


def string_fluctuation_per_dim_endpoint(
    sigma_GeV2: float,
    R0_fm: float,
    xi_fm: float,
) -> float:
    """Transverse string fluctuation per dimension at a free endpoint.

    For a string with one end fixed (Y-junction vertex) and one free end
    (massless quark), the standard mode-sum gives a factor of 2 enhancement
    relative to the fixed-fixed midpoint:

        ⟨X⊥²(L)⟩_per_dim = (1/(π σ)) × ln(L / ξ_UV)

    This is the physically appropriate formula for the Y-junction's quark
    end. We report it alongside the midpoint formula to bracket the
    modelling choice.

    Args:
        sigma_GeV2: String tension [GeV²].
        R0_fm:      String length L = R₀ [fm].
        xi_fm:      UV cutoff ξ = ξ_QCD [fm].

    Returns:
        ⟨X⊥²⟩ in fm² per transverse dimension.
    """
    return 2.0 * string_fluctuation_per_dim_midpoint(sigma_GeV2, R0_fm, xi_fm)


def string_fluctuation_contribution(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    R0_fm: float | None = None,
    xi_fm: float = XI_QCD_FM,
    endpoint: bool = False,
) -> float:
    """Total ⟨δr²⟩_string_fluct contribution to r_p² (3D, both transverse dims).

    Each Y-arm carries a quark whose position fluctuates transversely.
    There are 2 transverse dimensions per string, so:

        ⟨δr²⟩_per_quark = 2 × ⟨X⊥²⟩_per_dim

    In an equilateral uud Y-junction the 3 quarks all have the same
    fluctuation amplitude (geometry symmetry); the charge-weighted
    mean-square deviation is just the per-quark value, since each quark
    has charge magnitude similar enough that c_geom = 1 for the position
    sum and c_geom = 1 for the fluctuation sum (independent fluctuations
    average independently).

    Args:
        sigma_GeV2: String tension [GeV²]. Default σ_QCD.
        R0_fm:      Y-arm length [fm]. Default 1/√σ.
        xi_fm:      UV cutoff [fm]. Default ξ_QCD.
        endpoint:   If True, use free-endpoint formula (factor 2 larger).

    Returns:
        Contribution to r_p² in fm².
    """
    if R0_fm is None:
        R0_fm = equilibrium_radius_dimensional(sigma_GeV2)
    fn = (
        string_fluctuation_per_dim_endpoint
        if endpoint
        else string_fluctuation_per_dim_midpoint
    )
    return 2.0 * fn(sigma_GeV2, R0_fm, xi_fm)


# ---------------------------------------------------------------------------
# Full r_p² with decomposition
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ProtonRadiusV2:
    """Full charge-radius prediction with all components.

    Attributes:
        label:           Configuration label.
        R0_fm:           Y-junction equilibrium radius [fm].
        r2_position:     Position contribution R₀² [fm²].
        r2_string:       String-fluctuation contribution [fm²].
        r2_kink:         Intrinsic kink-width contribution [fm²].
        r2_total:        Substrate prediction r_p² [fm²].
        r_p_fm:          Substrate prediction r_p [fm].
        r2_pion_lit:     Literature pion-cloud value (cross-check) [fm²].
        r_p_with_pion:   r_p including pion-cloud literature value [fm].
        frac_error:      (r_p_pred - r_p_exp) / r_p_exp.
    """
    label: str
    R0_fm: float
    r2_position: float
    r2_string: float
    r2_kink: float
    r2_total: float
    r_p_fm: float
    r2_pion_lit: float
    r_p_with_pion: float
    frac_error: float


def proton_radius_full(
    R0_fm: float,
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_fm: float = XI_QCD_FM,
    kink_dim: int = 3,
    endpoint: bool = False,
    label: str = "estimate",
) -> ProtonRadiusV2:
    """Full r_p² prediction with explicit 4-component decomposition.

        r_p² = R₀² + ⟨δr²⟩_string_fluct + ⟨r²⟩_intrinsic_kink

    The pion-cloud literature value is reported as a separate column
    (cross-check only — not added to the substrate prediction, since the
    substrate's pion = chiral kink IS the string-fluctuation contribution).

    Args:
        R0_fm:      Y-junction equilibrium radius [fm].
        sigma_GeV2: String tension [GeV²].
        xi_fm:      Kink coherence width [fm].
        kink_dim:   1 or 3 — dimensionality of the kink-width contribution.
        endpoint:   If True, use free-endpoint string-fluctuation formula.
        label:      Tag for the estimate.

    Returns:
        ProtonRadiusV2 with all components and fractional error.
    """
    if kink_dim == 1:
        r2_kink = intrinsic_r2_kink_1D(xi_fm)
    elif kink_dim == 3:
        r2_kink = intrinsic_r2_kink_3D(xi_fm)
    else:
        raise ValueError(f"kink_dim must be 1 or 3, got {kink_dim}")

    r2_position = R0_fm ** 2
    r2_string = string_fluctuation_contribution(
        sigma_GeV2=sigma_GeV2, R0_fm=R0_fm, xi_fm=xi_fm, endpoint=endpoint
    )

    r2_total = r2_position + r2_string + r2_kink
    r_p_fm = math.sqrt(r2_total)

    r2_pion_lit = R2_PION_CLOUD_LIT_FM2
    r_p_with_pion = math.sqrt(r2_total + r2_pion_lit)

    frac_error = (r_p_fm - R_PROTON_CHARGE_EXP_FM) / R_PROTON_CHARGE_EXP_FM

    return ProtonRadiusV2(
        label=label,
        R0_fm=R0_fm,
        r2_position=r2_position,
        r2_string=r2_string,
        r2_kink=r2_kink,
        r2_total=r2_total,
        r_p_fm=r_p_fm,
        r2_pion_lit=r2_pion_lit,
        r_p_with_pion=r_p_with_pion,
        frac_error=frac_error,
    )


# ---------------------------------------------------------------------------
# Aggregate summary
# ---------------------------------------------------------------------------

@dataclass
class ProtonRadiusV2Summary:
    """Aggregate of all v2 r_p estimates."""
    sigma_GeV2: float
    xi_fm: float
    estimates: list[ProtonRadiusV2]


def compute_proton_radius_v2_summary(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_fm: float = XI_QCD_FM,
) -> ProtonRadiusV2Summary:
    """Compute v2 estimates with string-fluctuation contribution.

    We report 4 estimates that combine:
      • R₀ choice: dimensional (1/√σ) vs rotating-Y (2 M_p / 3πσ).
      • String-fluctuation formula: midpoint (fixed-fixed) vs endpoint (free).

    All use the 3D kink-width contribution since the bare 3D estimate is
    the most physically motivated (the kink is a localised soliton, smeared
    isotropically by the substrate primitives).

    Args:
        sigma_GeV2: String tension [GeV²].
        xi_fm:      Kink coherence width [fm].

    Returns:
        ProtonRadiusV2Summary with all four estimates.
    """
    R0_dim = equilibrium_radius_dimensional(sigma_GeV2)
    R0_rot = equilibrium_radius_rotating_Y(sigma_GeV2, M_PROTON_GEV)

    estimates = [
        proton_radius_full(
            R0_dim, sigma_GeV2, xi_fm, kink_dim=3, endpoint=False,
            label="dimensional R₀ + midpoint fluct",
        ),
        proton_radius_full(
            R0_dim, sigma_GeV2, xi_fm, kink_dim=3, endpoint=True,
            label="dimensional R₀ + endpoint fluct",
        ),
        proton_radius_full(
            R0_rot, sigma_GeV2, xi_fm, kink_dim=3, endpoint=False,
            label="rotating-Y R₀ + midpoint fluct",
        ),
        proton_radius_full(
            R0_rot, sigma_GeV2, xi_fm, kink_dim=3, endpoint=True,
            label="rotating-Y R₀ + endpoint fluct",
        ),
    ]

    return ProtonRadiusV2Summary(
        sigma_GeV2=sigma_GeV2,
        xi_fm=xi_fm,
        estimates=estimates,
    )
