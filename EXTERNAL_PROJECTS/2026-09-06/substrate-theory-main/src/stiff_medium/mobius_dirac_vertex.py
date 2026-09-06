"""Möbius-Dirac vertex eigenproblem for lepton stiffness eigenvalues κ_n.

Spec context (§§18.4, 18.10, 18.30, 18.35):
- §18.30: Three lepton generations from vertex closure (Z₃ symmetry).
- §18.35: m c² = ℏ √(κ/I) — masses from stiffness eigenvalues.
- §18.10: Möbius half-flux U(1) connection A = (1/2) dθ.
- §18.4:  Half-flux holonomy ψ(θ + 2π) = -ψ(θ) → fermionic statistics.

GOAL OF THIS MODULE
===================
Solve the 1D radial Dirac equation on a discrete approximation of the
45° cone with a Möbius half-flux connection, sectored by Z₃ angular
symmetry.  Extract the three lowest stiffness eigenvalues κ_n and ask:

    Do the resulting m_n / M ratios match the Foot parametrization
        m_n = M (1 + √2 cos(2πn/3 + δ))²
    with δ close to the empirical 1.404π?

If yes, §18.30 is structurally complete.  If no, the topology
specified there is incomplete and additional ingredients (e.g. a
specific kink potential, a k-running stiffness model) are required.

GEOMETRY
========
The 45° cone has angular structure carried by a U(1) bundle with
connection A = (1/2) dθ.  In the ANGULAR sector, the Möbius half-flux
phase combines with an integer m_z angular quantum number to give
effective angular momentum ℓ_eff = m_z + 1/2 (half-integer).

We project onto the THREE lowest angular sectors of the Z₃-vertex
closure, namely m_z ∈ {−1, 0, +1} (or equivalently the Z₃-eigenstates
of phase 0, 2π/3, 4π/3).  For each sector, we solve the radial Dirac
equation with a kink mass profile m(r) = M tanh(r/ξ) and harvest the
lowest bound-state energy.  These three energies are the eigenvalues
of the vertex stiffness operator.

DISCRETIZATION
==============
Radial coordinate r ∈ (0, R_max] sampled at N_radial points (≤500 to
keep the matrix small).  We use the Wilson-fermion squared-Dirac form
for the upper component:

    H_+ = -∂_r² + ℓ_eff(ℓ_eff + 1)/r² + m(r)² + m'(r)

The κ eigenvalue is identified with E² (the squared Dirac mass), so
that m_n c² = ℏ √(κ_n / I) reads κ_n ≡ E_n² in natural units (ℏ = c = I = 1).

For the Z₃-symmetry projection, we use ℓ_eff = j + 1/2 for j = −1, 0, +1
giving ℓ_eff ∈ {−1/2, +1/2, +3/2}.  The half-integer values come from
the Möbius half-flux holonomy.  The squared centrifugal term ℓ(ℓ+1)
gives {−1/4, +3/4, +15/4} — an irregular pattern, exactly the kind
of three-tier splitting we want.

Honest expectation: this minimal model will NOT give δ ≈ 1.404π out of
the box.  It will give *some* δ, and we will report it honestly.  The
Foot Q = 2/3 is structural to the parametrization, so the three masses
will line up on the Koide surface AT the implied δ.  But matching the
PHASE δ to 1.404π requires more than this minimal topology.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize_scalar

# ---------------------------------------------------------------------------
# Foot calibration targets (PDG 2024 lepton masses)
# ---------------------------------------------------------------------------

M_E_MEV: float = 0.51099895
M_MU_MEV: float = 105.6583755
M_TAU_MEV: float = 1776.86

# Foot calibration: M_FOOT and δ_emp such that
#   m_n = M_FOOT (1 + √2 cos(2πn/3 + δ_emp))²
M_FOOT_MEV: float = 313.86
DELTA_EMP_PI: float = 1.404  # δ_emp / π

EMPIRICAL_RATIOS = (
    M_E_MEV / M_FOOT_MEV,
    M_MU_MEV / M_FOOT_MEV,
    M_TAU_MEV / M_FOOT_MEV,
)  # ≈ (0.001628, 0.3367, 5.661)


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------


@dataclass
class KappaEigenvalues:
    """The three stiffness eigenvalues from the Möbius-Dirac vertex.

    Attributes:
        kappa: Array of three κ values (sorted ascending).
        E: Array of three √κ values (Dirac energies).
        ell_eff: Effective half-integer ℓ values used per sector.
        N_radial: Number of radial grid points used.
        R_max: Outer radial cutoff.
        M_kink: Asymptotic kink mass parameter.
        xi: Kink width.
        successful: Whether eigensolver converged sensibly.
        notes: Honest description of any caveats.
    """

    kappa: np.ndarray
    E: np.ndarray
    ell_eff: np.ndarray
    N_radial: int
    R_max: float
    M_kink: float
    xi: float
    successful: bool
    notes: str


@dataclass
class FootComparison:
    """Comparison of computed κ_n to Foot parametrization.

    Attributes:
        kappa_normalized: κ_n / κ_min (so smallest = 1).
        m_predicted_over_M: m_n / M_predicted, with M_predicted from fit.
        m_empirical_over_M: (m_e, m_μ, m_τ) / M_FOOT_MEV.
        delta_fit_pi: best-fit δ in units of π.
        M_fit: best-fit overall scale.
        delta_emp_pi: empirical δ_emp/π = 1.404 for reference.
        koide_q: Koide invariant of the predicted three masses.
        koide_target: 2/3.
        ratio_mu_e_predicted: m_2/m_1 from κ_n.
        ratio_tau_e_predicted: m_3/m_1 from κ_n.
        ratio_mu_e_target: empirical m_μ/m_e ≈ 207.
        ratio_tau_e_target: empirical m_τ/m_e ≈ 3477.
        residual_delta_pi: |δ_fit − 1.404|π.
        verdict: honest summary text.
    """

    kappa_normalized: np.ndarray
    m_predicted_over_M: np.ndarray
    m_empirical_over_M: np.ndarray
    delta_fit_pi: float
    M_fit: float
    delta_emp_pi: float
    koide_q: float
    koide_target: float
    ratio_mu_e_predicted: float
    ratio_tau_e_predicted: float
    ratio_mu_e_target: float
    ratio_tau_e_target: float
    residual_delta_pi: float
    verdict: str


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------


def _radial_grid(N_radial: int, R_max: float) -> tuple[np.ndarray, float]:
    """Build a uniform radial grid on (0, R_max].

    The grid avoids r = 0 (where the centrifugal term diverges) by
    starting at r_1 = dr.

    Args:
        N_radial: Number of radial grid points.
        R_max: Outer cutoff radius.

    Returns:
        (r_grid, dr) — the array of radii and the grid spacing.
    """
    dr = R_max / float(N_radial)
    r = np.linspace(dr, R_max, N_radial)
    return r, dr


def _kink_mass_profile(
    r: np.ndarray,
    M_kink: float,
    xi: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Pöschl-Teller kink mass m(r) = M tanh(r/ξ) and its derivative.

    The kink is centered at r = 0 (the cone tip).  m(r) → M as r → ∞.

    Args:
        r: Radial grid (positive).
        M_kink: Asymptotic mass.
        xi: Kink width.

    Returns:
        (m, m_prime) — mass and its derivative on the grid.
    """
    m = M_kink * np.tanh(r / xi)
    m_prime = (M_kink / xi) / np.cosh(r / xi) ** 2
    return m, m_prime


def _build_radial_kinetic(N_radial: int, dr: float) -> sp.csr_matrix:
    """Build the 1D second-derivative operator -∂_r² as sparse matrix.

    Three-point central difference with Dirichlet BC at r = 0 and r = R_max.

    Args:
        N_radial: Grid size.
        dr: Grid spacing.

    Returns:
        Sparse N×N matrix for -∂_r².
    """
    diag = np.full(N_radial, 2.0 / dr**2)
    off = np.full(N_radial - 1, -1.0 / dr**2)
    return sp.diags([off, diag, off], offsets=[-1, 0, 1], format="csr")


def _centrifugal_potential(r: np.ndarray, ell_eff: float) -> np.ndarray:
    """Centrifugal barrier ℓ_eff (ℓ_eff + 1) / r² on the grid.

    For half-integer ℓ_eff this can be NEGATIVE (e.g. ℓ_eff = −1/2 gives
    −1/4 / r² — a weak attractive singularity).  We regularize at r → 0
    by the grid offset (smallest r is dr, not 0).

    Args:
        r: Radial grid.
        ell_eff: Effective angular quantum number (half-integer).

    Returns:
        Diagonal potential array.
    """
    return ell_eff * (ell_eff + 1.0) / (r**2)


# ---------------------------------------------------------------------------
# Main API: build operator and solve
# ---------------------------------------------------------------------------


def build_mobius_dirac_operator(
    N_radial: int,
    M_kink: float,
    xi: float,
    R_max: float,
    ell_eff: float,
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Build the squared-Dirac operator for a single Möbius angular sector.

    For the squared upper-component Dirac operator on the cone with a
    radial kink mass m(r) and effective half-integer angular momentum
    ℓ_eff, the eigenvalue equation is

        [ -∂_r² + ℓ_eff(ℓ_eff+1)/r² + m(r)² + m'(r) ] ψ = E² ψ

    Args:
        N_radial: Number of radial grid points (≤ 500).
        M_kink: Asymptotic mass (sets the energy scale).
        xi: Kink width.
        R_max: Outer cutoff.
        ell_eff: Effective ℓ for this sector (half-integer).

    Returns:
        (H, r) — sparse Hamiltonian and the radial grid.
    """
    r, dr = _radial_grid(N_radial, R_max)
    m, m_prime = _kink_mass_profile(r, M_kink, xi)
    V_cent = _centrifugal_potential(r, ell_eff)
    V_total = V_cent + m**2 + m_prime
    T = _build_radial_kinetic(N_radial, dr)
    H = T + sp.diags(V_total, format="csr")
    return H, r


def solve_kappa_eigenvalues(
    N_radial: int = 400,
    M_kink: float = 1.0,
    xi: float = 1.0,
    R_max: Optional[float] = None,
    z3_sector_offset: int = 0,
) -> KappaEigenvalues:
    """Solve for the three κ eigenvalues from Möbius-Dirac vertex.

    Strategy: solve three separate radial eigenvalue problems, one per
    Z₃-vertex sector, with effective half-integer angular momentum:

        m_z ∈ { -1,  0, +1 }  (Z₃ sectors)
        ℓ_eff = m_z + 1/2     (Möbius half-flux shift)
              ∈ {-1/2, +1/2, +3/2}

    Take the LOWEST bound-state energy from each sector.  The squared
    energies are the κ eigenvalues.

    Sector-offset note: pure Z₃ symmetry around the Möbius bundle gives
    m_z values modulo 3.  The "natural" three sectors are {−1, 0, +1}
    centered on m_z = 0 (corresponding to z3_sector_offset = 0).

    Args:
        N_radial: Radial grid size (≤ 500).
        M_kink: Mass scale.
        xi: Kink width.
        R_max: Outer cutoff.  Defaults to 12 ξ.
        z3_sector_offset: Integer shift of the m_z labelling.

    Returns:
        KappaEigenvalues with three κ values and metadata.
    """
    if R_max is None:
        R_max = 12.0 * xi

    if N_radial > 500:
        raise ValueError(
            f"N_radial = {N_radial} exceeds memory budget; use ≤ 500."
        )

    # The three Z₃ sectors with Möbius half-flux shift
    m_z_vals = np.array([-1, 0, +1]) + z3_sector_offset
    ell_eff_vals = m_z_vals + 0.5

    energies = np.zeros(3)
    notes_parts = []

    for i, ell_eff in enumerate(ell_eff_vals):
        H, r = build_mobius_dirac_operator(
            N_radial=N_radial,
            M_kink=M_kink,
            xi=xi,
            R_max=R_max,
            ell_eff=ell_eff,
        )
        # Sparse eigensolver: smallest 4 eigenvalues (use shift-invert)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                eigvals = spla.eigsh(
                    H,
                    k=min(6, N_radial - 2),
                    sigma=-1e-3,  # shift below ground-state to find lowest
                    which="LM",
                    return_eigenvectors=False,
                    tol=1e-9,
                )
            eigvals = np.sort(eigvals)
        except Exception as exc:  # pragma: no cover
            notes_parts.append(f"sector ℓ={ell_eff}: solver failed ({exc})")
            energies[i] = np.nan
            continue

        # Filter physical (positive, real, below threshold M²)
        threshold = M_kink**2 * 0.999  # asymptote of m(r)²
        physical = eigvals[(eigvals > 1e-12) & (eigvals < threshold)]
        if len(physical) == 0:
            # Some sectors might have NO bound state below threshold
            # (e.g. ℓ_eff = -1/2 with weak attractive centrifugal can have
            # only continuum); fall back to the smallest positive eigenvalue
            physical = eigvals[eigvals > 0]
            notes_parts.append(
                f"sector ℓ_eff={ell_eff}: no truly bound state, "
                "using lowest positive eigenvalue"
            )

        if len(physical) == 0:
            energies[i] = np.nan
            notes_parts.append(f"sector ℓ_eff={ell_eff}: no positive eigenvalue")
        else:
            E_lowest = float(np.sqrt(physical[0]))
            energies[i] = E_lowest

    # κ_n = E_n²  (in natural units where ℏ = I = 1)
    valid = ~np.isnan(energies)
    successful = bool(valid.all() and np.all(energies > 0))

    if successful:
        kappa = energies**2
        # Sort ascending (so κ_0 is smallest)
        order = np.argsort(kappa)
        kappa = kappa[order]
        energies = energies[order]
        ell_eff_sorted = ell_eff_vals[order]
    else:
        kappa = energies**2
        ell_eff_sorted = ell_eff_vals

    notes = "; ".join(notes_parts) if notes_parts else "all sectors converged"

    return KappaEigenvalues(
        kappa=kappa,
        E=energies,
        ell_eff=np.asarray(ell_eff_sorted, dtype=float),
        N_radial=N_radial,
        R_max=R_max,
        M_kink=M_kink,
        xi=xi,
        successful=successful,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Foot comparison
# ---------------------------------------------------------------------------


def _foot_masses_array(delta: float, M: float = 1.0) -> np.ndarray:
    """Foot parametrization: m_n = M (1 + √2 cos(2πn/3 + δ))², sorted."""
    ns = np.array([0.0, 1.0, 2.0])
    x = 1.0 + np.sqrt(2.0) * np.cos(2.0 * np.pi * ns / 3.0 + delta)
    m = M * x**2
    return np.sort(m)


def _koide_q(masses: np.ndarray) -> float:
    """Koide invariant for three masses."""
    if np.any(masses <= 0):
        return float("nan")
    return float(np.sum(masses) / np.sum(np.sqrt(masses)) ** 2)


def _fit_foot_to_kappa(
    kappa: np.ndarray,
) -> tuple[float, float, float]:
    """Find (δ, M) such that Foot masses best match m_n = √(κ_n/I).

    We treat the three κ_n as proportional to the three m_n² via
    m_n c² = ℏ √(κ_n / I).  In natural units this is simply m_n = √κ_n.
    Then we fit the Foot family m_n = M (1 + √2 cos(2πn/3 + δ))² to the
    sorted √κ values.

    Args:
        kappa: Sorted (ascending) array of three κ eigenvalues.

    Returns:
        (delta_pi, M_fit, residual) where delta_pi = δ/π.
    """
    # Convert κ to "m" via m = √κ
    m_target = np.sqrt(np.maximum(kappa, 0.0))
    m_target_sorted = np.sort(m_target)

    def loss(delta: float) -> float:
        m_foot = _foot_masses_array(delta, M=1.0)
        # The three m_foot are determined up to overall scale; fit M
        # that minimizes log-ratio² error.
        if np.any(m_foot <= 1e-30):
            return 1e12
        # Only ratios matter.
        # Use least-squares of log(m_target_n / m_foot_n) → constant
        log_diff = np.log(m_target_sorted) - np.log(m_foot)
        # Best M is exp(mean(log_diff)); residual is variance of log_diff
        return float(np.var(log_diff))

    # Global scan + refine
    deltas = np.linspace(0.0, 2.0 * np.pi, 4000)
    losses = np.array([loss(d) for d in deltas])
    best_idx = int(np.argmin(losses))
    d_coarse = float(deltas[best_idx])

    res = minimize_scalar(
        loss,
        bounds=(d_coarse - 0.05, d_coarse + 0.05),
        method="bounded",
        options={"xatol": 1e-10},
    )
    delta_best = float(res.x % (2.0 * np.pi))

    # Best M
    m_foot = _foot_masses_array(delta_best, M=1.0)
    log_diff = np.log(m_target_sorted) - np.log(m_foot)
    log_M = float(np.mean(log_diff))
    M_best = float(np.exp(log_M))

    residual = float(res.fun) ** 0.5
    return delta_best / np.pi, M_best, residual


def compare_to_foot(kappa_n: np.ndarray) -> FootComparison:
    """Extract δ, M from κ_n eigenvalues; compare to empirical Foot.

    Treats m_n = √κ_n (in natural units) and fits the Foot family.

    Args:
        kappa_n: Array of three κ eigenvalues (will be sorted internally).

    Returns:
        FootComparison with fit parameters and verdict.
    """
    kappa_sorted = np.sort(np.asarray(kappa_n, dtype=float))
    kappa_min = kappa_sorted[0] if kappa_sorted[0] > 0 else 1e-30
    kappa_norm = kappa_sorted / kappa_min

    delta_fit_pi, M_fit, residual_log = _fit_foot_to_kappa(kappa_sorted)

    # Predicted masses (in natural units of √κ)
    m_predicted = _foot_masses_array(delta_fit_pi * np.pi, M=M_fit)

    # Normalize to fit-M for "m_n / M_predicted"
    m_pred_over_M = m_predicted / M_fit

    # Empirical reference
    m_emp_over_M = np.array(EMPIRICAL_RATIOS)

    q_pred = _koide_q(m_predicted)

    # Mass ratios (m_2/m_1 ≈ 207, m_3/m_1 ≈ 3477 if matching empirical)
    if m_predicted[0] > 0:
        r10 = float(m_predicted[1] / m_predicted[0])
        r20 = float(m_predicted[2] / m_predicted[0])
    else:
        r10 = float("nan")
        r20 = float("nan")
    r10_target = M_MU_MEV / M_E_MEV
    r20_target = M_TAU_MEV / M_E_MEV

    # Distance from empirical δ
    residual_delta_pi = float(min(
        abs(delta_fit_pi - DELTA_EMP_PI),
        abs(delta_fit_pi - DELTA_EMP_PI - 2.0),
        abs(delta_fit_pi - DELTA_EMP_PI + 2.0),
    ))

    # Verdict
    if residual_delta_pi < 0.02:
        verdict = (
            f"PASS: δ = {delta_fit_pi:.4f}π is within 0.02π of "
            f"empirical {DELTA_EMP_PI}π."
        )
    elif residual_delta_pi < 0.1:
        verdict = (
            f"NEAR: δ = {delta_fit_pi:.4f}π differs from empirical "
            f"{DELTA_EMP_PI}π by {residual_delta_pi:.3f}π — close but "
            "not exact."
        )
    else:
        verdict = (
            f"FAIL: δ = {delta_fit_pi:.4f}π differs from empirical "
            f"{DELTA_EMP_PI}π by {residual_delta_pi:.3f}π. The minimal "
            "Möbius half-flux + Z₃ topology does NOT predict the empirical "
            "Foot phase. The §18.30 specification is incomplete: additional "
            "ingredients (kink-profile-specific stiffness, k-running, or "
            "non-trivial cone curvature) are required."
        )

    return FootComparison(
        kappa_normalized=kappa_norm,
        m_predicted_over_M=m_pred_over_M,
        m_empirical_over_M=m_emp_over_M,
        delta_fit_pi=delta_fit_pi,
        M_fit=M_fit,
        delta_emp_pi=DELTA_EMP_PI,
        koide_q=q_pred,
        koide_target=2.0 / 3.0,
        ratio_mu_e_predicted=r10,
        ratio_tau_e_predicted=r20,
        ratio_mu_e_target=r10_target,
        ratio_tau_e_target=r20_target,
        residual_delta_pi=residual_delta_pi,
        verdict=verdict,
    )


# ---------------------------------------------------------------------------
# Convenience high-level driver
# ---------------------------------------------------------------------------


def run_mobius_dirac_vertex(
    N_radial: int = 400,
    M_kink: float = 1.0,
    xi: float = 1.0,
    R_max: Optional[float] = None,
    z3_sector_offset: int = 0,
) -> dict:
    """Run the full Möbius-Dirac vertex calculation.

    Args:
        N_radial: Radial grid size.
        M_kink: Mass scale.
        xi: Kink width.
        R_max: Outer cutoff (default 12ξ).
        z3_sector_offset: Z₃ sector shift (default 0 → m_z ∈ {-1,0,+1}).

    Returns:
        Dict with 'kappa', 'comparison', 'parameters'.
    """
    kappa_result = solve_kappa_eigenvalues(
        N_radial=N_radial,
        M_kink=M_kink,
        xi=xi,
        R_max=R_max,
        z3_sector_offset=z3_sector_offset,
    )
    if kappa_result.successful:
        comparison = compare_to_foot(kappa_result.kappa)
    else:
        comparison = None

    return {
        "kappa": kappa_result,
        "comparison": comparison,
        "parameters": {
            "N_radial": N_radial,
            "M_kink": M_kink,
            "xi": xi,
            "R_max": kappa_result.R_max,
            "z3_sector_offset": z3_sector_offset,
        },
    }
