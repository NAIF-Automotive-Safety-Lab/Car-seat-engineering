"""Extended Möbius-Dirac vertex eigenproblems — beyond the §18.61.4 minimal model.

Spec context (§§18.10, 18.30, 18.35, 18.48):
- §18.61.4: minimal Möbius half-flux + Z₃ vertex topology gave κ_2/κ_0 ≈ 1.084
  (8% spread) → mass ratios m_2/m_0 ≈ 1.04, far below the empirical (×207, ×3477).
- §18.48.2: empirical fit κ_n ∝ (n+1)^k with k ≈ 15 reproduces the ratios
  individually but breaks Koide Q = 2/3.
- The "stress-loading" picture (§18.30) suggests profile-specific stiffness:
  loaded vertex configurations have higher local stress → larger κ.

GOAL OF THIS MODULE
===================
Test three physically motivated extensions to the minimal §18.30 topology
that might produce the exponential lepton mass hierarchy and δ ≈ 1.404π:

  Variant A (PROFILE-SPECIFIC STIFFNESS): replace constant kink mass M
            by m(r) = M · cosh(r/ξ)^p  for stress-loading exponents p ∈ {1, 2, 3}.
            This is the Pöschl-Teller-like profile with explicit stress
            growth at large r — the stress-loaded vertex sees an effective
            mass that grows with the cosh-power.

  Variant B (CONE CURVATURE): include the actual 45° cone metric on the
            angular sector. The 2D Dirac on a cone has an explicit
            curvature term. With Möbius half-flux ℓ_eff = m_z + 1/2,
            the eigenvalue problem on the cone is
                [-∂_r² - (1/r)∂_r + ℓ_eff²/r² + m(r)² + m'(r)] ψ = E² ψ
            (the 1/r kinetic term comes from the cone-Laplacian, NOT
            from the centrifugal barrier alone).

  Variant C (NESTED Z₃ TREE): two-level Z₃ topology. Outer triplet of
            stiffnesses κ^outer_a (a=0,1,2), each with inner triplet
            κ^inner_b (b=0,1,2). The 9 vertex eigenvalues come in three
            "tiers" with hierarchical scale set by the inner/outer ratio.
            Take the lowest of each tier as the (e, μ, τ) candidate.

For each variant we solve for κ_n and report:
  (i)   eigenvalue ratio match: κ_2/κ_0 vs target ≈ 1.21e7 (i.e. m_τ²/m_e²)
        and √(κ ratio) vs target (×207, ×3477)
  (ii)  Koide Q (target 2/3)
  (iii) implied δ via Foot inversion
  (iv)  whether the §18.48.2 power-law κ_n ∝ (n+1)^k is recovered with
        k ≈ 15

CRITICAL: We do NOT hand-tune the exponent p, cone angle, or topology to
fit the ratios. We try p ∈ {1, 2, 3, 5, 10, 15} as physically motivated
integer / half-integer values, and report them all honestly.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize_scalar

from stiff_medium.mobius_dirac_vertex import (
    DELTA_EMP_PI,
    EMPIRICAL_RATIOS,
    KappaEigenvalues,
    M_E_MEV,
    M_FOOT_MEV,
    M_MU_MEV,
    M_TAU_MEV,
    _foot_masses_array,
    _koide_q,
    compare_to_foot,
)


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------


@dataclass
class VariantResult:
    """Result of a single extended-Möbius variant.

    Attributes:
        name: Variant name (A, B, or C plus parameter label).
        kappa: Three sorted κ eigenvalues (or NaN if failed).
        E_n: √κ_n (Dirac energies).
        ratio_mu_e: κ_1/κ_0 (or for the mass scale, √(κ_1/κ_0)).
        ratio_tau_e: κ_2/κ_0 (or √(κ_2/κ_0)).
        ratio_mu_e_target: empirical √(κ_μ/κ_e) ≈ 207.
        ratio_tau_e_target: empirical √(κ_τ/κ_e) ≈ 3477.
        koide_q: Koide invariant for the predicted three masses.
        delta_pi: best-fit Foot phase δ/π.
        residual_delta_pi: |δ_fit − 1.404|π (modulo 2).
        power_law_k: best-fit k for κ_n ∝ (n+1)^k.
        successful: Whether eigensolver converged.
        notes: Honest description of caveats.
        params: Dict of parameters used (p, cone angle, ratio, etc.).
    """

    name: str
    kappa: np.ndarray
    E_n: np.ndarray
    ratio_mu_e: float
    ratio_tau_e: float
    ratio_mu_e_target: float
    ratio_tau_e_target: float
    koide_q: float
    delta_pi: float
    residual_delta_pi: float
    power_law_k: float
    successful: bool
    notes: str
    params: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Common building blocks
# ---------------------------------------------------------------------------


def _radial_grid(N: int, R_max: float) -> tuple[np.ndarray, float]:
    dr = R_max / float(N)
    return np.linspace(dr, R_max, N), dr


def _kinetic_1d(N: int, dr: float) -> sp.csr_matrix:
    """Standard -∂_r² with Dirichlet BC (3-pt central difference)."""
    diag = np.full(N, 2.0 / dr**2)
    off = np.full(N - 1, -1.0 / dr**2)
    return sp.diags([off, diag, off], offsets=[-1, 0, 1], format="csr")


def _kinetic_cone_1d(r: np.ndarray, dr: float) -> sp.csr_matrix:
    """Cone-Laplacian radial part: -∂_r² - (1/r)∂_r.

    On a 2D cone with metric ds² = dr² + r² dθ² (i.e., flat polar metric),
    the radial Laplacian is -∂_r² - (1/r)∂_r. The 45°-cone deficit angle
    is absorbed into the angular eigenvalue (modifying m_z range), but
    the radial structure is the same as flat 2D polar — that's the key
    observation: a CONE is locally flat AWAY from the apex, so it
    differs from a plane only in the global angular periodicity.
    """
    N = r.size
    # Second derivative
    d2_diag = np.full(N, 2.0 / dr**2)
    d2_off = np.full(N - 1, -1.0 / dr**2)
    T = sp.diags([d2_off, d2_diag, d2_off], offsets=[-1, 0, 1], format="csr")
    # First derivative: (1/r) × [(ψ_{i+1} - ψ_{i-1})/(2 dr)]
    inv_r = 1.0 / r
    upper = -inv_r[:-1] / (2.0 * dr)
    lower = inv_r[1:] / (2.0 * dr)
    D1 = sp.diags([lower, np.zeros(N), upper], offsets=[-1, 0, 1], format="csr")
    return T - D1


def _solve_lowest(H: sp.csr_matrix, N: int, n_eigs: int = 6) -> np.ndarray:
    """Find the lowest n_eigs eigenvalues of a symmetric sparse H."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ev = spla.eigsh(
                H,
                k=min(n_eigs, N - 2),
                sigma=-1e-3,
                which="LM",
                return_eigenvectors=False,
                tol=1e-9,
            )
        return np.sort(ev)
    except Exception:
        return np.array([np.nan] * n_eigs)


# ---------------------------------------------------------------------------
# Variant A: profile-specific stiffness (cosh^p mass profile)
# ---------------------------------------------------------------------------


def variant_a_profile_stiffness(
    p: float,
    N_radial: int = 300,
    M_kink: float = 1.0,
    xi: float = 1.0,
    R_max: float = 8.0,
) -> VariantResult:
    """Variant A: m(r) = M · tanh(r/ξ) · cosh(r/ξ)^p — stress-loaded profile.

    Physical picture: the stress-loaded vertex experiences a mass
    background that grows as the kink "breathes" out at large r. The
    cosh^p factor models substrate stiffness rising with displacement
    (cf. §18.48 stress-loading).

    For p=0 we recover the minimal §18.61.4 case.
    For p>0 the effective mass squared term in the squared-Dirac
    operator grows as cosh^(2p), pushing eigenvalues up exponentially.

    NOTE: cosh grows exponentially, so we must use a SMALLER R_max
    when p is large to avoid the asymptotic mass term diverging.
    We use R_max ≈ 5 ξ for p ≤ 5 and R_max ≈ 3 ξ for p > 5.

    Args:
        p: Stress-loading exponent. Try integer values 1, 2, 3, 5, 10, 15.
        N_radial: Grid size.
        M_kink: Mass scale.
        xi: Kink width.
        R_max: Outer cutoff (auto-scaled for large p).

    Returns:
        VariantResult.
    """
    # Auto-scale R_max to avoid cosh^p exploding past machine precision
    if p >= 10:
        R_max = min(R_max, 2.5 * xi)
    elif p >= 5:
        R_max = min(R_max, 4.0 * xi)
    elif p >= 2:
        R_max = min(R_max, 5.0 * xi)

    r, dr = _radial_grid(N_radial, R_max)
    # Stress-loaded mass profile: m(r) = M · tanh(r/ξ) · cosh(r/ξ)^p
    # tanh keeps the kink structure; cosh^p adds the stress-loading growth.
    u = r / xi
    cosh_u = np.cosh(u)
    sinh_u = np.sinh(u)
    tanh_u = np.tanh(u)
    sech2 = 1.0 / cosh_u**2

    m = M_kink * tanh_u * cosh_u**p
    # m'(r) = M·[ sech² · cosh^p + tanh · p · cosh^(p-1) · sinh ] / ξ
    #       = M/ξ · [ cosh^(p-2) + p · sinh² · cosh^(p-2) ]
    #       = M/ξ · cosh^(p-2) · (1 + p sinh²)
    m_prime = (M_kink / xi) * cosh_u**(p - 2) * (1.0 + p * sinh_u**2)

    energies = np.zeros(3)
    notes_parts = []

    # Same Z₃ × Möbius half-flux sectors as the minimal version
    ell_eff_vals = np.array([-0.5, 0.5, 1.5])

    T = _kinetic_1d(N_radial, dr)

    for i, ell_eff in enumerate(ell_eff_vals):
        V_cent = ell_eff * (ell_eff + 1.0) / r**2
        V_total = V_cent + m**2 + m_prime
        H = T + sp.diags(V_total, format="csr")
        ev = _solve_lowest(H, N_radial)
        # Filter: positive, real, finite
        physical = ev[(ev > 1e-12) & np.isfinite(ev)]
        if len(physical) == 0:
            energies[i] = np.nan
            notes_parts.append(f"sector ℓ={ell_eff}: no positive eigenvalue")
        else:
            energies[i] = float(np.sqrt(physical[0]))

    valid = ~np.isnan(energies)
    successful = bool(valid.all() and np.all(energies > 0))

    if successful:
        kappa = np.sort(energies**2)
        E_n = np.sqrt(kappa)
    else:
        kappa = energies**2
        E_n = energies

    return _build_variant_result(
        name=f"A_p={p}",
        kappa=kappa,
        E_n=E_n,
        successful=successful,
        notes="; ".join(notes_parts) if notes_parts else "all sectors converged",
        params={"p": p, "N_radial": N_radial, "M_kink": M_kink, "xi": xi, "R_max": R_max},
    )


# ---------------------------------------------------------------------------
# Variant B: cone curvature corrections
# ---------------------------------------------------------------------------


def variant_b_cone_curvature(
    cone_angle_deg: float = 45.0,
    N_radial: int = 300,
    M_kink: float = 1.0,
    xi: float = 1.0,
    R_max: float = 12.0,
) -> VariantResult:
    """Variant B: Dirac on a 45° cone with the actual cone-Laplacian.

    On a 2D cone with metric ds² = dr² + α²r² dθ² (deficit angle 2π(1-α)),
    the radial Laplacian for an angular mode m_z (Möbius shifted by 1/2)
    is

        Δ_r = ∂_r² + (1/r)∂_r - ℓ_eff²/(α²r²)

    where ℓ_eff = m_z + 1/2. Note that the EFFECTIVE angular momentum
    is ℓ_eff/α — the deficit MULTIPLIES the centrifugal term.

    For a 45° cone, α = 45°/360° × 4 = 1/2 (the cone wraps half the
    angular range of a flat plane). So ℓ_eff/α = 2·ℓ_eff doubles the
    centrifugal scale.

    Wait — let me be more careful. The "45° cone" of §18.10 refers to
    the OPENING angle of the cone embedded in 3D. For a cone with
    half-opening angle θ_0 (measured from the cone axis), the induced
    metric on the cone (in cone-intrinsic radial coordinate r and
    azimuth φ ∈ [0, 2π)) is

        ds² = dr² + (sin θ_0)² · r² dφ²

    so α = sin(θ_0). For θ_0 = 45°, α = sin(45°) = 1/√2, and
    ℓ_eff/α = √2 · ℓ_eff. The centrifugal term scales as 2 ℓ_eff².

    The squared-Dirac operator on the cone reads
        [-∂_r² - (1/r)∂_r + ℓ_eff²/(α²r²) + m(r)² + m'(r)] ψ = E² ψ

    Args:
        cone_angle_deg: Half-opening angle in degrees (default 45°).
        N_radial: Grid size.
        M_kink: Mass scale.
        xi: Kink width.
        R_max: Outer cutoff.

    Returns:
        VariantResult.
    """
    alpha = np.sin(np.deg2rad(cone_angle_deg))
    r, dr = _radial_grid(N_radial, R_max)

    # Cone-Laplacian radial part
    T_cone = _kinetic_cone_1d(r, dr)

    m = M_kink * np.tanh(r / xi)
    m_prime = (M_kink / xi) / np.cosh(r / xi)**2

    energies = np.zeros(3)
    notes_parts = []

    # Möbius half-flux: ell_eff = m_z + 1/2 for m_z ∈ {-1, 0, +1}
    ell_eff_vals = np.array([-0.5, 0.5, 1.5])

    for i, ell_eff in enumerate(ell_eff_vals):
        # Centrifugal on cone: ℓ_eff² / (α² r²)
        V_cent = ell_eff**2 / (alpha**2 * r**2)
        V_total = V_cent + m**2 + m_prime
        H = T_cone + sp.diags(V_total, format="csr")
        # Make symmetric by averaging with transpose (cone-Laplacian's
        # 1/r drift is not exactly symmetric in the discretization).
        H_sym = 0.5 * (H + H.T)
        ev = _solve_lowest(H_sym, N_radial)
        physical = ev[(ev > 1e-12) & np.isfinite(ev)]
        if len(physical) == 0:
            energies[i] = np.nan
            notes_parts.append(f"sector ℓ={ell_eff}: no positive eigenvalue")
        else:
            energies[i] = float(np.sqrt(physical[0]))

    valid = ~np.isnan(energies)
    successful = bool(valid.all() and np.all(energies > 0))

    if successful:
        kappa = np.sort(energies**2)
        E_n = np.sqrt(kappa)
    else:
        kappa = energies**2
        E_n = energies

    return _build_variant_result(
        name=f"B_cone={cone_angle_deg}deg",
        kappa=kappa,
        E_n=E_n,
        successful=successful,
        notes="; ".join(notes_parts) if notes_parts else "all sectors converged",
        params={
            "cone_angle_deg": cone_angle_deg,
            "alpha_sin": alpha,
            "N_radial": N_radial,
            "M_kink": M_kink,
            "xi": xi,
            "R_max": R_max,
        },
    )


# ---------------------------------------------------------------------------
# Variant C: nested Z₃ × Z₃ vertex tree
# ---------------------------------------------------------------------------


def variant_c_nested_z3(
    inner_outer_ratio: float = 3.0,
    N_radial: int = 250,
    M_kink: float = 1.0,
    xi: float = 1.0,
    R_max: float = 12.0,
) -> VariantResult:
    """Variant C: nested Z₃ × Z₃ vertex topology.

    Picture: outer triplet of vertex states (m_z^outer ∈ {-1, 0, +1}),
    each carrying an inner triplet (m_z^inner ∈ {-1, 0, +1}) with its
    own Möbius half-flux. The 9 effective angular momenta are
        ℓ_eff(a, b) = (m_z^outer + 1/2) + ratio · (m_z^inner + 1/2)
    where `ratio` is the inner/outer scale ratio (set by the geometry
    of nested Möbius bundles — physically the inner bundle has higher
    curvature so its half-flux contribution is amplified).

    For each (a, b) sector we solve the radial Dirac and harvest the
    lowest bound-state energy. The 9 energies cluster into three "tiers"
    by outer label `a`. Take the LOWEST of each tier as the (e, μ, τ)
    candidate.

    Args:
        inner_outer_ratio: Scale factor for the inner bundle's
            angular contribution. Try {1, 2, 3, 5, 10}.
        N_radial: Grid size.
        M_kink: Mass scale.
        xi: Kink width.
        R_max: Outer cutoff.

    Returns:
        VariantResult.
    """
    r, dr = _radial_grid(N_radial, R_max)
    m = M_kink * np.tanh(r / xi)
    m_prime = (M_kink / xi) / np.cosh(r / xi)**2

    T = _kinetic_1d(N_radial, dr)

    # 9 sectors: outer m_z^o ∈ {-1,0,1}, inner m_z^i ∈ {-1,0,1}
    energies_grid = np.zeros((3, 3))

    for a, mz_o in enumerate([-1, 0, 1]):
        for b, mz_i in enumerate([-1, 0, 1]):
            ell_outer = mz_o + 0.5
            ell_inner = mz_i + 0.5
            ell_eff_total = ell_outer + inner_outer_ratio * ell_inner

            V_cent = ell_eff_total * (ell_eff_total + 1.0) / r**2
            # For very negative ell_eff the centrifugal can be wildly
            # attractive; cap it sensibly to avoid spurious -∞ states.
            V_total = V_cent + m**2 + m_prime
            H = T + sp.diags(V_total, format="csr")
            ev = _solve_lowest(H, N_radial)
            physical = ev[(ev > 1e-12) & np.isfinite(ev)]
            if len(physical) == 0:
                energies_grid[a, b] = np.nan
            else:
                energies_grid[a, b] = float(np.sqrt(physical[0]))

    notes_parts = []
    # For each outer tier, take the lowest energy across inner triplet
    tier_energies = np.zeros(3)
    for a in range(3):
        row = energies_grid[a, :]
        valid_row = row[~np.isnan(row) & (row > 0)]
        if len(valid_row) == 0:
            tier_energies[a] = np.nan
            notes_parts.append(f"outer tier a={a}: all sectors failed")
        else:
            tier_energies[a] = float(np.min(valid_row))

    valid = ~np.isnan(tier_energies)
    successful = bool(valid.all() and np.all(tier_energies > 0))

    if successful:
        kappa = np.sort(tier_energies**2)
        E_n = np.sqrt(kappa)
    else:
        kappa = tier_energies**2
        E_n = tier_energies

    return _build_variant_result(
        name=f"C_ratio={inner_outer_ratio}",
        kappa=kappa,
        E_n=E_n,
        successful=successful,
        notes="; ".join(notes_parts) if notes_parts else "all 9 sectors converged",
        params={
            "inner_outer_ratio": inner_outer_ratio,
            "N_radial": N_radial,
            "M_kink": M_kink,
            "xi": xi,
            "R_max": R_max,
            "energies_grid": energies_grid.tolist(),
        },
    )


# ---------------------------------------------------------------------------
# Hierarchical κ ansatz (analytic check on §18.48.2)
# ---------------------------------------------------------------------------


def variant_c_hierarchical_powerlaw(
    k_exponent: float,
    M_kink: float = 1.0,
) -> VariantResult:
    """Variant C analytic form: κ_n = κ_0 · (n+1)^k.

    This is the §18.48.2 empirical fit form, plugged in directly. We check
    whether ANY topological mechanism (e.g. nested Möbius) could produce
    this k. The answer comes from numerics; this is a sanity reference.

    Args:
        k_exponent: Power-law exponent (try k ∈ {2, 5, 10, 14.84, 15, 15.38}).
        M_kink: Overall mass scale.

    Returns:
        VariantResult.
    """
    n_arr = np.arange(3)  # 0, 1, 2
    kappa = M_kink**2 * (n_arr + 1.0) ** k_exponent
    E_n = np.sqrt(kappa)

    return _build_variant_result(
        name=f"C_powerlaw_k={k_exponent}",
        kappa=kappa,
        E_n=E_n,
        successful=True,
        notes=f"analytic ansatz κ_n ∝ (n+1)^{k_exponent}",
        params={"k_exponent": k_exponent, "M_kink": M_kink},
    )


# ---------------------------------------------------------------------------
# Common scoring / Foot-fit logic
# ---------------------------------------------------------------------------


def _power_law_fit(kappa: np.ndarray) -> float:
    """Fit κ_n = κ_0 · (n+1)^k by least-squares on log(κ_n) vs log(n+1).

    Args:
        kappa: Three sorted κ values.

    Returns:
        Best-fit k.
    """
    if np.any(kappa <= 0) or np.any(~np.isfinite(kappa)):
        return float("nan")
    n_arr = np.arange(3)
    x = np.log(n_arr + 1.0)
    y = np.log(kappa)
    # Linear fit y = log(κ_0) + k · x
    # (subtract first value to remove κ_0)
    if np.all(x == 0):
        return float("nan")
    # Use slope from indices 1, 2 only (x_0 = 0 makes least-squares trivial)
    A = np.vstack([x, np.ones_like(x)]).T
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(sol[0])


def _build_variant_result(
    name: str,
    kappa: np.ndarray,
    E_n: np.ndarray,
    successful: bool,
    notes: str,
    params: dict,
) -> VariantResult:
    """Build a VariantResult with all derived quantities."""
    if not successful or np.any(~np.isfinite(kappa)) or np.any(kappa <= 0):
        return VariantResult(
            name=name,
            kappa=kappa,
            E_n=E_n,
            ratio_mu_e=float("nan"),
            ratio_tau_e=float("nan"),
            ratio_mu_e_target=M_MU_MEV / M_E_MEV,
            ratio_tau_e_target=M_TAU_MEV / M_E_MEV,
            koide_q=float("nan"),
            delta_pi=float("nan"),
            residual_delta_pi=float("nan"),
            power_law_k=float("nan"),
            successful=False,
            notes=notes,
            params=params,
        )

    # Mass = √κ, normalized so smallest = 1
    masses = np.sqrt(kappa)
    masses_sorted = np.sort(masses)

    ratio_mu_e = float(masses_sorted[1] / masses_sorted[0])
    ratio_tau_e = float(masses_sorted[2] / masses_sorted[0])
    koide_q = _koide_q(masses_sorted)

    # Foot phase fit
    try:
        comparison = compare_to_foot(kappa)
        delta_pi = comparison.delta_fit_pi
        residual_delta_pi = comparison.residual_delta_pi
    except Exception:
        delta_pi = float("nan")
        residual_delta_pi = float("nan")

    # Power-law k fit on κ_n
    k_fit = _power_law_fit(np.sort(kappa))

    return VariantResult(
        name=name,
        kappa=np.sort(kappa),
        E_n=np.sort(E_n),
        ratio_mu_e=ratio_mu_e,
        ratio_tau_e=ratio_tau_e,
        ratio_mu_e_target=M_MU_MEV / M_E_MEV,
        ratio_tau_e_target=M_TAU_MEV / M_E_MEV,
        koide_q=koide_q,
        delta_pi=delta_pi,
        residual_delta_pi=residual_delta_pi,
        power_law_k=k_fit,
        successful=True,
        notes=notes,
        params=params,
    )


# ---------------------------------------------------------------------------
# Combined scoring
# ---------------------------------------------------------------------------


def score_variant(result: VariantResult) -> dict:
    """Score a variant on three axes.

    Score components (lower is better):
      ratio_score: log10-distance from empirical (m_μ/m_e, m_τ/m_e)
      koide_score: |Q − 2/3|
      delta_score: |δ_fit − 1.404π| / π

    Args:
        result: VariantResult.

    Returns:
        Dict with score breakdown.
    """
    if not result.successful:
        return {
            "ratio_score": float("inf"),
            "koide_score": float("inf"),
            "delta_score": float("inf"),
            "total_score": float("inf"),
        }

    target_mu_e = M_MU_MEV / M_E_MEV
    target_tau_e = M_TAU_MEV / M_E_MEV

    if result.ratio_mu_e <= 0 or result.ratio_tau_e <= 0:
        ratio_score = float("inf")
    else:
        ratio_score = float(
            abs(np.log10(result.ratio_mu_e / target_mu_e))
            + abs(np.log10(result.ratio_tau_e / target_tau_e))
        )

    koide_score = float(abs(result.koide_q - 2.0 / 3.0))
    delta_score = float(result.residual_delta_pi)

    total = ratio_score + koide_score + delta_score
    return {
        "ratio_score": ratio_score,
        "koide_score": koide_score,
        "delta_score": delta_score,
        "total_score": total,
    }


# ---------------------------------------------------------------------------
# High-level driver: scan all variants
# ---------------------------------------------------------------------------


def run_all_extended_variants(
    p_values: Sequence[float] = (1, 2, 3, 5, 10, 15),
    cone_angles_deg: Sequence[float] = (45.0, 30.0, 60.0, 90.0),
    inner_outer_ratios: Sequence[float] = (1.0, 2.0, 3.0, 5.0, 10.0),
    k_powerlaw_values: Sequence[float] = (2.0, 5.0, 10.0, 14.84, 15.0, 15.38),
    N_radial: int = 250,
    M_kink: float = 1.0,
    xi: float = 1.0,
) -> dict:
    """Run all three variants with a sweep of parameters.

    Returns:
        Dict with keys 'variant_a', 'variant_b', 'variant_c',
        'powerlaw_reference', each containing a list of VariantResult.
    """
    results = {
        "variant_a": [],
        "variant_b": [],
        "variant_c_numerical": [],
        "powerlaw_reference": [],
    }

    for p in p_values:
        try:
            r = variant_a_profile_stiffness(p=float(p), N_radial=N_radial,
                                            M_kink=M_kink, xi=xi)
            results["variant_a"].append(r)
        except Exception as exc:
            results["variant_a"].append(
                _build_variant_result(
                    name=f"A_p={p}",
                    kappa=np.array([np.nan, np.nan, np.nan]),
                    E_n=np.array([np.nan, np.nan, np.nan]),
                    successful=False,
                    notes=f"exception: {exc}",
                    params={"p": p},
                )
            )

    for angle in cone_angles_deg:
        try:
            r = variant_b_cone_curvature(cone_angle_deg=float(angle),
                                         N_radial=N_radial, M_kink=M_kink,
                                         xi=xi)
            results["variant_b"].append(r)
        except Exception as exc:
            results["variant_b"].append(
                _build_variant_result(
                    name=f"B_cone={angle}deg",
                    kappa=np.array([np.nan, np.nan, np.nan]),
                    E_n=np.array([np.nan, np.nan, np.nan]),
                    successful=False,
                    notes=f"exception: {exc}",
                    params={"cone_angle_deg": angle},
                )
            )

    for ratio in inner_outer_ratios:
        try:
            r = variant_c_nested_z3(inner_outer_ratio=float(ratio),
                                    N_radial=min(N_radial, 200),
                                    M_kink=M_kink, xi=xi)
            results["variant_c_numerical"].append(r)
        except Exception as exc:
            results["variant_c_numerical"].append(
                _build_variant_result(
                    name=f"C_ratio={ratio}",
                    kappa=np.array([np.nan, np.nan, np.nan]),
                    E_n=np.array([np.nan, np.nan, np.nan]),
                    successful=False,
                    notes=f"exception: {exc}",
                    params={"inner_outer_ratio": ratio},
                )
            )

    for k in k_powerlaw_values:
        try:
            r = variant_c_hierarchical_powerlaw(k_exponent=float(k),
                                                M_kink=M_kink)
            results["powerlaw_reference"].append(r)
        except Exception as exc:
            results["powerlaw_reference"].append(
                _build_variant_result(
                    name=f"C_powerlaw_k={k}",
                    kappa=np.array([np.nan, np.nan, np.nan]),
                    E_n=np.array([np.nan, np.nan, np.nan]),
                    successful=False,
                    notes=f"exception: {exc}",
                    params={"k_exponent": k},
                )
            )

    return results


def find_best_variant(all_results: dict) -> tuple[str, VariantResult, dict]:
    """Find the variant with the lowest combined score.

    Args:
        all_results: Output of run_all_extended_variants.

    Returns:
        (variant_name, best_result, score_dict).
    """
    best_score = float("inf")
    best_name = None
    best_result = None
    best_score_dict = None

    for category, lst in all_results.items():
        for r in lst:
            sc = score_variant(r)
            if sc["total_score"] < best_score:
                best_score = sc["total_score"]
                best_name = f"{category}/{r.name}"
                best_result = r
                best_score_dict = sc

    return best_name, best_result, best_score_dict
