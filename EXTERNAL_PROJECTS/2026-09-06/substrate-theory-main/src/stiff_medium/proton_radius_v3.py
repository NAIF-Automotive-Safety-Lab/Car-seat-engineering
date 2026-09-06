"""Proton charge radius v3: full 3-body wavefunction in the Y-junction (§18.61.4).

Physics summary
---------------
v2 (proton_radius_v2.py) used R₀ = 1/√σ (single-string equilibrium) for the
position-spread term and added Lüscher transverse string fluctuations on top.
That closed the gap from −33% (bare) to about −22% (canonical v2 with
endpoint formula). The remaining gap is the difference between the
single-string equilibrium length and the actual 3-quark wave-function
SPREAD inside the Y-junction.

This module computes that spread from first principles by solving the
Schrödinger problem for one effective constituent in the linear potential
V(r) = σ r against the Y-vertex (the COM frame of the 3 quarks is the
vertex itself, so the relative coordinate from each quark to the vertex
is the natural radial variable). The constituent inertia is the substrate-
derived m_eff = √σ ≈ 424 MeV (chromomagnetic_substrate.py).

Two independent calculations are performed:

    (1) Variational Gaussian
            ψ(r) = (α/π)^(3/4) exp(-α r²/2)
            ⟨T⟩ = (3/4) α ℏ²/μ
            ⟨V⟩ = σ × 2/√(πα)
            E(α) minimised in α  ⇒  α_min = [(4 μ σ)/(3 √π ℏ²)]^(2/3)
            ⟨r²⟩ = 3/(2 α_min)

    (2) Numerical sparse-matrix radial Schrödinger
            -ℏ²/(2μ) [u''(r)]  +  σ r u(r)  =  E u(r),  u = r ψ
            Discretised on a finite grid r ∈ (0, R_max], 200 points,
            Dirichlet boundary u(0)=u(R_max)=0, ground-state extracted
            via scipy.sparse.linalg.eigsh.

The variational answer is the analytic anchor; the numerical answer is
the cross-check (the true ground state is NOT an exact Gaussian — it is
proportional to Ai((r − E/σ)·(2μσ/ℏ²)^(1/3)) for s-wave, an Airy function).

Total decomposition (v3)
------------------------
The full proton charge-radius² is

    r_p² = ⟨r²⟩_3body + ⟨X⊥²⟩_string_fluct + r²_kink_width

with each term physically distinct:
    • ⟨r²⟩_3body         — 3-quark wavefunction spread (this module)
    • ⟨X⊥²⟩_string_fluct — Lüscher transverse string fluctuations
                            (proton_radius_v2.string_fluctuation_contribution)
    • r²_kink_width      — intrinsic sech² kink width
                            (proton_radius.intrinsic_r2_kink_3D)

For the string-fluctuation contribution we use R₀_eff = √⟨r²⟩_3body as the
characteristic Y-arm length (this is the correct length scale for the
log enhancement, since the actual quark sits at rms distance √⟨r²⟩ from
the vertex, NOT at 1/√σ).

Inputs (NONE fitted to r_p):
    σ_QCD   = 0.180 GeV²  (K(ξ) running, regge_spectrum.py)
    ξ_QCD   = 0.2 fm      (substrate coherence length)
    m_eff   = √σ ≈ 0.424 GeV  (chromomagnetic_substrate.py)

Honest accounting
-----------------
The variational Gaussian is provably ABOVE the true ground-state energy.
For ⟨r²⟩ specifically, the variational Gaussian under-estimates the spread
slightly because the true Airy-function ground state has a longer outer
tail (V grows linearly forever — there is no exponential cutoff). The
numerical eigensolver corrects this and gives the unbiased answer.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from stiff_medium.proton_radius import (
    HBARC_GEV_FM,
    R_PROTON_CHARGE_EXP_FM,
    SIGMA_QCD_GEV2,
    XI_QCD_FM,
    intrinsic_r2_kink_3D,
)
from stiff_medium.proton_radius_v2 import (
    R2_PION_CLOUD_LIT_FM2,
    string_fluctuation_per_dim_endpoint,
    string_fluctuation_per_dim_midpoint,
)

# ---------------------------------------------------------------------------
# Substrate-derived constituent inertia
# ---------------------------------------------------------------------------

def m_eff_from_sigma(sigma_GeV2: float = SIGMA_QCD_GEV2) -> float:
    """Effective constituent mass m_eff = √σ.

    Chromomagnetic substrate result (chromomagnetic_substrate.py §18.62):
    the kink localisation scale R₀ = 1/√σ implies a kinetic inertia
    m_eff = ℏc / (R₀ × c²) → in natural units, m_eff = √σ.

    Args:
        sigma_GeV2: String tension [GeV²].

    Returns:
        m_eff in GeV.
    """
    return math.sqrt(sigma_GeV2)


# ---------------------------------------------------------------------------
# (1) Variational Gaussian
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class VariationalResult:
    """Variational ground state of -ℏ²/(2μ) ∇² + σ r in 3D.

    Attributes:
        alpha_min_GeV2:  Optimal Gaussian width parameter [GeV²].
        E_min_GeV:       Variational ground-state energy [GeV].
        r2_GeV_inv2:     ⟨r²⟩ in GeV⁻².
        r2_fm2:          ⟨r²⟩ in fm² (after × (ℏc)²).
        r_rms_fm:        √⟨r²⟩ in fm.
    """
    alpha_min_GeV2: float
    E_min_GeV: float
    r2_GeV_inv2: float
    r2_fm2: float
    r_rms_fm: float


def solve_variational_quark_wavefunction(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    m_eff_GeV: float | None = None,
) -> VariationalResult:
    """Minimise ⟨H⟩ for ψ(r) = (α/π)^(3/4) exp(-α r²/2) in V = σ r.

    Energy expectation (natural units ℏ = 1):

        E(α) = (3 α) / (4 μ) + 2 σ / √(π α)

    Stationarity dE/dα = 0:

        3 / (4 μ)  =  σ / √π · α^(-3/2)
        ⇒ α_min^(3/2) = (4 μ σ) / (3 √π)
        ⇒ α_min = [(4 μ σ) / (3 √π)]^(2/3)

    Then

        E_min   = (3 α_min) / (4 μ) + 2 σ / √(π α_min)
        ⟨r²⟩   = 3 / (2 α_min)

    α has units GeV² (so r² ~ 1/α has units GeV⁻²); convert via
    (ℏc)² = 0.0389 GeV²·fm².

    Args:
        sigma_GeV2: String tension σ [GeV²].
        m_eff_GeV:  Constituent inertia μ [GeV]; default √σ.

    Returns:
        VariationalResult.
    """
    if m_eff_GeV is None:
        m_eff_GeV = m_eff_from_sigma(sigma_GeV2)

    mu = m_eff_GeV
    sig = sigma_GeV2

    alpha_pow_3_2 = (4.0 * mu * sig) / (3.0 * math.sqrt(math.pi))   # GeV³
    alpha_min = alpha_pow_3_2 ** (2.0 / 3.0)                         # GeV²

    E_min = (3.0 * alpha_min) / (4.0 * mu) \
        + 2.0 * sig / math.sqrt(math.pi * alpha_min)                 # GeV

    r2_GeV_inv2 = 3.0 / (2.0 * alpha_min)                            # GeV⁻²
    r2_fm2 = r2_GeV_inv2 * HBARC_GEV_FM ** 2                         # fm²
    r_rms_fm = math.sqrt(r2_fm2)

    return VariationalResult(
        alpha_min_GeV2=alpha_min,
        E_min_GeV=E_min,
        r2_GeV_inv2=r2_GeV_inv2,
        r2_fm2=r2_fm2,
        r_rms_fm=r_rms_fm,
    )


# ---------------------------------------------------------------------------
# (2) Numerical sparse-matrix radial Schrödinger
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NumericalResult:
    """Numerical ground state of the s-wave radial equation.

    The reduced radial wavefunction u(r) = r ψ(r) satisfies

        -ℏ²/(2μ) u''(r) + σ r u(r) = E u(r),  u(0) = u(R_max) = 0.

    Attributes:
        E_GeV:         Ground-state energy [GeV].
        r2_GeV_inv2:   ⟨r²⟩ from numerical eigenvector [GeV⁻²].
        r2_fm2:        ⟨r²⟩ in fm² (after × (ℏc)²).
        r_rms_fm:      √⟨r²⟩ in fm.
        N_r:           Number of grid points used.
        R_max_GeV_inv: Radial cutoff used [GeV⁻¹].
    """
    E_GeV: float
    r2_GeV_inv2: float
    r2_fm2: float
    r_rms_fm: float
    N_r: int
    R_max_GeV_inv: float


def solve_radial_schrodinger_numerical(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    m_eff_GeV: float | None = None,
    N_r: int = 200,
    R_max_GeV_inv: float = 30.0,
) -> NumericalResult:
    """Solve the s-wave radial equation by sparse matrix eigendecomposition.

    Discretisation:
        Uniform grid r_i = i·dr, i = 1, …, N_r,  dr = R_max / (N_r + 1)
        (i = 0 and i = N_r+1 are boundary points where u = 0.)

        Kinetic operator (3-point stencil, second derivative):
            u''(r_i) ≈ [u_{i+1} − 2 u_i + u_{i−1}] / dr²

        Hamiltonian:
            H_ii    = +ℏ²/(μ dr²) + σ r_i
            H_i,i±1 = −ℏ²/(2μ dr²)

        Lowest eigenvalue → E; eigenvector u → ψ(r) = u(r)/r normalised
        so ∫|ψ|² 4π r² dr = 1, i.e. ∫|u|² 4π dr = 1 / r factors? Actually
        for u = r ψ in 3D, the normalisation of ψ(r) reads

            ∫|ψ|² d³r = 4π ∫|u|² dr  (since |u|² = r²|ψ|²)

        So ∫|u|² dr = 1/(4π).  Then

            ⟨r²⟩ = ∫|ψ|² r² 4π r² dr = 4π ∫|u|² r² dr.

    Memory budget: N_r × N_r tridiagonal matrix in CSR; ~3 N_r doubles ≈
    5 KB for N_r = 200. Trivial.

    Args:
        sigma_GeV2:      String tension [GeV²].
        m_eff_GeV:       Constituent mass [GeV]; default √σ.
        N_r:             Number of interior grid points (200 is plenty).
        R_max_GeV_inv:   Radial cutoff in GeV⁻¹ (≈ 5.9 fm at 30/GeV).

    Returns:
        NumericalResult.
    """
    if m_eff_GeV is None:
        m_eff_GeV = m_eff_from_sigma(sigma_GeV2)

    mu = m_eff_GeV
    sig = sigma_GeV2

    # Grid: r_i = i·dr for i = 1..N_r (i = 0 and i = N_r+1 are boundary).
    dr = R_max_GeV_inv / (N_r + 1)             # GeV⁻¹
    r = np.arange(1, N_r + 1) * dr             # GeV⁻¹

    # Kinetic prefactor 1/(2μ dr²) [GeV]; ℏ = 1.
    kin = 1.0 / (2.0 * mu * dr ** 2)           # GeV

    # Diagonal: +2 × kin + σ r
    main_diag = 2.0 * kin + sig * r            # GeV
    # Off-diagonals: −kin
    off_diag = -kin * np.ones(N_r - 1)         # GeV

    H = diags(
        diagonals=[off_diag, main_diag, off_diag],
        offsets=[-1, 0, 1],
        format="csr",
    )

    # Lowest eigenvalue / -vector (sigma=0 shift-invert is faster, but
    # for N_r = 200 the basic eigsh on small matrix is fine; we use
    # which='SA' (smallest algebraic) for robustness.)
    E_arr, U_arr = eigsh(H, k=1, which="SA")
    E = float(E_arr[0])                         # GeV
    u = U_arr[:, 0]                             # un-normalised

    # Normalise ψ in 3D: ∫|ψ|² d³r = 4π ∫|u|² dr = 1.
    norm_u2 = np.sum(u ** 2) * dr               # GeV⁻¹  (∫|u|² dr)
    # Want 4π × norm_u2 = 1  ⇒  rescale u by 1/√(4π × norm_u2).
    scale = 1.0 / math.sqrt(4.0 * math.pi * norm_u2)
    u_norm = u * scale

    # ⟨r²⟩ = 4π ∫|u|² r² dr.
    r2 = 4.0 * math.pi * np.sum(u_norm ** 2 * r ** 2) * dr   # GeV⁻²
    r2_fm2 = r2 * HBARC_GEV_FM ** 2

    return NumericalResult(
        E_GeV=E,
        r2_GeV_inv2=r2,
        r2_fm2=r2_fm2,
        r_rms_fm=math.sqrt(r2_fm2),
        N_r=N_r,
        R_max_GeV_inv=R_max_GeV_inv,
    )


# ---------------------------------------------------------------------------
# Full v3 prediction
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ProtonRadiusV3:
    """Full v3 charge-radius prediction.

    Attributes:
        method:          'variational' or 'numerical'.
        sigma_GeV2:      String tension used.
        xi_fm:           Kink coherence width used.
        m_eff_GeV:       Constituent mass used.
        r2_3body:        3-body wavefunction spread [fm²].
        r_rms_3body:     √⟨r²⟩_3body [fm] (= R₀_eff for string log).
        r2_string:       Lüscher string-fluct contribution [fm²].
        r2_kink:         Intrinsic kink width [fm²].
        r2_total:        Total r_p² [fm²].
        r_p_fm:          Total r_p [fm].
        frac_error:      (r_p − r_exp) / r_exp.
        E_GeV:           Ground-state energy [GeV] (cross-check).
        r2_pion_lit:     Literature pion-cloud value (cross-check) [fm²].
        r_p_with_pion:   r_p including pion-cloud literature value [fm].
    """
    method: str
    sigma_GeV2: float
    xi_fm: float
    m_eff_GeV: float
    r2_3body: float
    r_rms_3body: float
    r2_string: float
    r2_kink: float
    r2_total: float
    r_p_fm: float
    frac_error: float
    E_GeV: float
    r2_pion_lit: float
    r_p_with_pion: float


def proton_radius_full_v3(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_fm: float = XI_QCD_FM,
    m_eff_GeV: float | None = None,
    method: str = "variational",
    endpoint: bool = True,
    N_r: int = 200,
    R_max_GeV_inv: float = 30.0,
) -> ProtonRadiusV3:
    """Compute r_p² = ⟨r²⟩_3body + ⟨δr²⟩_string + r²_kink (v3).

    The 3-body wavefunction term replaces the v2 R₀² term. The
    string-fluctuation term uses R₀_eff = √⟨r²⟩_3body as the Y-arm length
    in the Lüscher log. The kink-width term is unchanged from v2.

    Args:
        sigma_GeV2:    String tension [GeV²].
        xi_fm:         Kink coherence width [fm].
        m_eff_GeV:     Constituent mass [GeV]; default √σ.
        method:        'variational' (default) or 'numerical'.
        endpoint:      Use free-endpoint Lüscher formula (default True
                       — physically correct for free quark ends).
        N_r:           Numerical grid size (numerical method only).
        R_max_GeV_inv: Numerical radial cutoff (numerical method only).

    Returns:
        ProtonRadiusV3.
    """
    if m_eff_GeV is None:
        m_eff_GeV = m_eff_from_sigma(sigma_GeV2)

    if method == "variational":
        v = solve_variational_quark_wavefunction(sigma_GeV2, m_eff_GeV)
        r2_3body = v.r2_fm2
        r_rms_3body = v.r_rms_fm
        E_GeV = v.E_min_GeV
    elif method == "numerical":
        n = solve_radial_schrodinger_numerical(
            sigma_GeV2, m_eff_GeV, N_r=N_r, R_max_GeV_inv=R_max_GeV_inv
        )
        r2_3body = n.r2_fm2
        r_rms_3body = n.r_rms_fm
        E_GeV = n.E_GeV
    else:
        raise ValueError(f"method must be 'variational' or 'numerical', got {method}")

    # String-fluctuation: use R₀_eff = √⟨r²⟩_3body in the log.
    if endpoint:
        rms_per_dim_fm2 = string_fluctuation_per_dim_endpoint(
            sigma_GeV2, r_rms_3body, xi_fm
        )
    else:
        rms_per_dim_fm2 = string_fluctuation_per_dim_midpoint(
            sigma_GeV2, r_rms_3body, xi_fm
        )
    r2_string = 2.0 * rms_per_dim_fm2     # 2 transverse dims

    r2_kink = intrinsic_r2_kink_3D(xi_fm)

    r2_total = r2_3body + r2_string + r2_kink
    r_p_fm = math.sqrt(r2_total)
    frac_error = (r_p_fm - R_PROTON_CHARGE_EXP_FM) / R_PROTON_CHARGE_EXP_FM

    r2_pion_lit = R2_PION_CLOUD_LIT_FM2
    r_p_with_pion = math.sqrt(r2_total + r2_pion_lit)

    return ProtonRadiusV3(
        method=method,
        sigma_GeV2=sigma_GeV2,
        xi_fm=xi_fm,
        m_eff_GeV=m_eff_GeV,
        r2_3body=r2_3body,
        r_rms_3body=r_rms_3body,
        r2_string=r2_string,
        r2_kink=r2_kink,
        r2_total=r2_total,
        r_p_fm=r_p_fm,
        frac_error=frac_error,
        E_GeV=E_GeV,
        r2_pion_lit=r2_pion_lit,
        r_p_with_pion=r_p_with_pion,
    )


# ---------------------------------------------------------------------------
# Aggregate summary
# ---------------------------------------------------------------------------

@dataclass
class ProtonRadiusV3Summary:
    """Aggregate of v3 estimates: variational vs numerical, midpoint vs endpoint."""
    sigma_GeV2: float
    xi_fm: float
    m_eff_GeV: float
    estimates: list[ProtonRadiusV3]


def compute_proton_radius_v3_summary(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    xi_fm: float = XI_QCD_FM,
    m_eff_GeV: float | None = None,
) -> ProtonRadiusV3Summary:
    """Compute v3 summary: variational and numerical, both Lüscher choices.

    Args:
        sigma_GeV2: String tension [GeV²].
        xi_fm:      Kink coherence width [fm].
        m_eff_GeV:  Constituent mass [GeV]; default √σ.

    Returns:
        ProtonRadiusV3Summary with 4 estimates.
    """
    if m_eff_GeV is None:
        m_eff_GeV = m_eff_from_sigma(sigma_GeV2)

    estimates = [
        proton_radius_full_v3(
            sigma_GeV2, xi_fm, m_eff_GeV,
            method="variational", endpoint=False,
        ),
        proton_radius_full_v3(
            sigma_GeV2, xi_fm, m_eff_GeV,
            method="variational", endpoint=True,
        ),
        proton_radius_full_v3(
            sigma_GeV2, xi_fm, m_eff_GeV,
            method="numerical", endpoint=False,
        ),
        proton_radius_full_v3(
            sigma_GeV2, xi_fm, m_eff_GeV,
            method="numerical", endpoint=True,
        ),
    ]

    return ProtonRadiusV3Summary(
        sigma_GeV2=sigma_GeV2,
        xi_fm=xi_fm,
        m_eff_GeV=m_eff_GeV,
        estimates=estimates,
    )
