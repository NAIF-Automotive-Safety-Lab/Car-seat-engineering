"""Confinement potential derived from substrate field dynamics (§18.49).

KEY PHYSICS: We derive — not impose — the inter-kink potential by solving for
the minimum-energy static field configuration between two kinks placed at
separation R in the 3D substrate.  If the medium has the right nonlinear
elastic response (sine-Gordon potential), the flux tube that bridges the two
kinks carries a constant energy density per unit length → linear confining
potential V(R) = σR.

DERIVATION
----------
The substrate field φ(x) satisfies the static (Euclidean) EOM:

    ∇²φ = (1/ξ²) sin(φ)        [static sine-Gordon]

A "kink-antikink" pair separated by R along the z-axis corresponds to
the field configuration that winds by +2π at z₁ and by -2π at z₂, with
vacuum φ=0 at both z → ±∞.  In 1D the exact solution is:

    φ_1D(z) = 4 arctan(exp((z−z₁)/ξ)) − 4 arctan(exp((z−z₂)/ξ))

This configuration has φ ≈ 2π in the "tube" z₁ < z < z₂ and φ ≈ 0
outside.  The tube energy per unit length is σ = K/ξ (exact in 1D).

In 3D the field disperses transversely.  The minimum-energy configuration
is found by relaxation (gradient descent on the energy functional).  The
central z-column is pinned to the 1D profile to enforce the topological
winding number; without this, the 3D field collapses to φ=0 everywhere
because the Dirichlet vacuum BCs do not protect the winding in 3D.

After relaxation, the total energy is integrated:

    E(R) = ∫ [ ½|∇φ|² + (K/ξ²)(1 − cos φ) ] d³x

and fitted to E(R) = σR + α_C/R + E₀ to extract σ.

TOPOLOGICAL NOTE
----------------
In genuine 1+1D QCD, kink topological charge is conserved (π₁(S¹) = Z).
In a 3+1D substrate the global topology is trivial and the tube CAN shrink
to zero unless kink positions are fixed.  We pin the field on a small
cylindrical core (radius ∼ pin_rad ξ) to the exact 1D profile.  This is
the 3D analogue of fixing quark positions in a lattice QCD Wilson loop
calculation.  The energy extracted is then E(R) − E(0), measuring the
"string" contribution from the field between the kinks.

PHYSICAL CONSTANTS (SI)
-----------------------
All SI conversion uses:
    K_SI = ℏc/ξ⁴   [from SubstratePrimitives, Solution A: ξ = λ_C(electron)]
    σ_SI = σ_lattice × K_SI × ξ_SI
    σ_GeV² = σ_SI × ℏc / (GeV²)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np
from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

C_SI: Final[float] = 2.99792458e8
HBAR_SI: Final[float] = 1.054571817e-34
M_E_SI: Final[float] = 9.1093837015e-31
M_P_SI: Final[float] = 1.67262192369e-27
M_P_MEV: Final[float] = 938.272
M_E_MEV: Final[float] = 0.51099895
PROTON_TO_ELECTRON: Final[float] = M_P_SI / M_E_SI
EV_TO_J: Final[float] = 1.602176634e-19
GEV_TO_J: Final[float] = EV_TO_J * 1e9
FM_TO_M: Final[float] = 1e-15
HBAR_C_GEV_FM: Final[float] = 0.19732698  # GeV·fm
SIGMA_QCD_GEV2: Final[float] = 0.18        # GeV²  (lattice QCD)


# ---------------------------------------------------------------------------
# Lattice parameters
# ---------------------------------------------------------------------------

@dataclass
class LatticeParams:
    """Parameters for the 3D relaxation lattice.

    Attributes:
        N:        Grid size (N × N × N).
        dx:       Lattice spacing in natural units (ξ = 1).
        eta:      Gradient descent step size.
        max_iter: Maximum relaxation iterations.
        tol:      Convergence tolerance on max |residual|.
        xi:       Coherence length (= 1.0 in natural units).
        pin_rad:  Radius (lattice sites) of the topological pinning core.
    """

    N: int = 40
    dx: float = 1.0
    eta: float = 0.05
    max_iter: int = 3000
    tol: float = 1e-5
    xi: float = 1.0
    pin_rad: int = 2


# ---------------------------------------------------------------------------
# 1D sine-Gordon kink profile
# ---------------------------------------------------------------------------

def kink_profile_1d(
    z_arr: NDArray[np.float64],
    z0: float,
    xi: float,
) -> NDArray[np.float64]:
    """Exact 1D sine-Gordon kink profile centred at z₀.

    φ_k(z) = 4 arctan(exp((z − z₀)/ξ))

    This goes from 0 at z → −∞ to 2π at z → +∞.

    Args:
        z_arr: Array of z coordinates.
        z0:    Kink centre position.
        xi:    Coherence length.

    Returns:
        Profile values at each z.
    """
    return 4.0 * np.arctan(np.exp((z_arr - z0) / xi))


def kink_antikink_1d(
    z_arr: NDArray[np.float64],
    z1: float,
    z2: float,
    xi: float,
) -> NDArray[np.float64]:
    """1D kink–antikink profile: φ = 0 outside, φ ≈ 2π in the tube.

    φ(z) = 4 arctan(exp((z − z₁)/ξ)) − 4 arctan(exp((z − z₂)/ξ))

    At z₁: +kink (φ rises through π).
    At z₂: −antikink (φ falls through π, returning to 0).
    In the tube z₁ < z < z₂: φ ≈ 2π.

    Args:
        z_arr: Array of z coordinates.
        z1:    First kink position (φ rises from 0 toward 2π).
        z2:    Second kink position (φ falls from 2π back to 0).
        xi:    Coherence length.

    Returns:
        1D field profile.
    """
    return (
        4.0 * np.arctan(np.exp((z_arr - z1) / xi))
        - 4.0 * np.arctan(np.exp((z_arr - z2) / xi))
    )


# ---------------------------------------------------------------------------
# 3D field construction
# ---------------------------------------------------------------------------

def build_initial_field(
    N: int,
    xi: float,
    z1: int,
    z2: int,
    center: int,
    tube_radius: float = 2.5,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Build the 3D initial field and the 1D pinning profile.

    Constructs φ(x,y,z) = φ_1D(z) × exp(−ρ²/(2 r_tube²)) where φ_1D is the
    exact 1D kink–antikink profile and ρ = sqrt((x−x_c)² + (y−y_c)²) is the
    transverse distance from the tube axis.

    Args:
        N:           Grid size.
        xi:          Coherence length (natural units).
        z1:          z-index of the first kink.
        z2:          z-index of the second kink.
        center:      x, y centre index.
        tube_radius: Initial Gaussian envelope radius (in units of ξ).

    Returns:
        Tuple of (phi_3d, phi_1d_profile).
        phi_1d_profile is used later for pinning.
    """
    z_idx = np.arange(N, dtype=float)
    phi_1d = kink_antikink_1d(z_idx, float(z1), float(z2), xi)

    phi = np.zeros((N, N, N))
    for ix in range(N):
        for iy in range(N):
            rho2 = (ix - center) ** 2 + (iy - center) ** 2
            envelope = math.exp(-rho2 / (2.0 * tube_radius ** 2))
            phi[ix, iy, :] = phi_1d * envelope

    # Vacuum BCs
    phi[0, :, :] = 0.0
    phi[-1, :, :] = 0.0
    phi[:, 0, :] = 0.0
    phi[:, -1, :] = 0.0
    phi[:, :, 0] = 0.0
    phi[:, :, -1] = 0.0

    return phi, phi_1d


# ---------------------------------------------------------------------------
# PDE solver
# ---------------------------------------------------------------------------

def _laplacian_7pt(
    phi: NDArray[np.float64],
    dx: float,
) -> NDArray[np.float64]:
    """7-point Laplacian stencil on interior points."""
    lap = np.zeros_like(phi)
    inv_dx2 = 1.0 / dx ** 2
    lap[1:-1, 1:-1, 1:-1] = inv_dx2 * (
        phi[2:, 1:-1, 1:-1]
        + phi[:-2, 1:-1, 1:-1]
        + phi[1:-1, 2:, 1:-1]
        + phi[1:-1, :-2, 1:-1]
        + phi[1:-1, 1:-1, 2:]
        + phi[1:-1, 1:-1, :-2]
        - 6.0 * phi[1:-1, 1:-1, 1:-1]
    )
    return lap


def _residual_sine_gordon(
    phi: NDArray[np.float64],
    dx: float,
    xi: float,
) -> NDArray[np.float64]:
    """Sine-Gordon residual R = −∇²φ + sin(φ)/ξ²."""
    return -_laplacian_7pt(phi, dx) + np.sin(phi) / xi ** 2


def _apply_vacuum_bc(phi: NDArray[np.float64]) -> None:
    """Set outer-boundary faces to zero (vacuum) in-place."""
    phi[0, :, :] = 0.0
    phi[-1, :, :] = 0.0
    phi[:, 0, :] = 0.0
    phi[:, -1, :] = 0.0
    phi[:, :, 0] = 0.0
    phi[:, :, -1] = 0.0


def _apply_pin(
    phi: NDArray[np.float64],
    phi_1d: NDArray[np.float64],
    center: int,
    pin_rad: int,
) -> None:
    """Pin the cylindrical core around the z-axis to the 1D profile.

    This enforces topological protection: the string cannot collapse to zero
    because the field on the axis is fixed to the exact kink–antikink profile.

    Args:
        phi:     3D field array (modified in place).
        phi_1d:  1D kink–antikink profile array of length N.
        center:  x, y index of the tube axis.
        pin_rad: Radius of the pinned cylinder.
    """
    N = phi.shape[0]
    for di in range(-pin_rad, pin_rad + 1):
        for dj in range(-pin_rad, pin_rad + 1):
            if di ** 2 + dj ** 2 <= pin_rad ** 2:
                ix, iy = center + di, center + dj
                if 0 <= ix < N and 0 <= iy < N:
                    phi[ix, iy, :] = phi_1d


def relax_field(
    phi: NDArray[np.float64],
    phi_1d: NDArray[np.float64],
    params: LatticeParams,
    center: int,
    verbose: bool = False,
) -> tuple[NDArray[np.float64], list[float], int]:
    """Relax the field to the sine-Gordon minimum via gradient descent.

    The topological pinning core is enforced after every update step to
    prevent the string from collapsing.

    Args:
        phi:      Initial field (not modified; a copy is made).
        phi_1d:   1D profile for pinning (shape N,).
        params:   Lattice and solver parameters.
        center:   x, y axis index.
        verbose:  Print convergence progress every 500 iterations.

    Returns:
        (relaxed_phi, residual_history, iterations_used)
    """
    phi = phi.copy()
    residuals: list[float] = []

    for iteration in range(params.max_iter):
        res = _residual_sine_gordon(phi, params.dx, params.xi)
        max_res = float(np.max(np.abs(res[1:-1, 1:-1, 1:-1])))
        residuals.append(max_res)

        # Gradient descent on interior
        phi[1:-1, 1:-1, 1:-1] -= params.eta * res[1:-1, 1:-1, 1:-1]

        # Enforce BCs and topological pin
        _apply_vacuum_bc(phi)
        _apply_pin(phi, phi_1d, center, params.pin_rad)

        if verbose and iteration % 500 == 0:
            print(f"  iter {iteration:5d}: max_res = {max_res:.3e}")

        if max_res < params.tol:
            if verbose:
                print(f"  Converged at iter {iteration} (max_res={max_res:.2e})")
            return phi, residuals, iteration

    if verbose:
        print(f"  Max iters reached; final max_res = {max_res:.2e}")
    return phi, residuals, params.max_iter


# ---------------------------------------------------------------------------
# Energy functional
# ---------------------------------------------------------------------------

def compute_field_energy(
    phi: NDArray[np.float64],
    dx: float,
    xi: float,
) -> float:
    """Integrate the sine-Gordon energy density over the lattice.

    ε(x) = ½ |∇φ|² + (1/ξ²)(1 − cos φ)

    The gradient uses central differences on interior points.

    Args:
        phi: Field configuration.
        dx:  Lattice spacing.
        xi:  Coherence length.

    Returns:
        Total dimensionless energy.
    """
    inv_dx = 0.5 / dx
    dphi_dx = (phi[2:, 1:-1, 1:-1] - phi[:-2, 1:-1, 1:-1]) * inv_dx
    dphi_dy = (phi[1:-1, 2:, 1:-1] - phi[1:-1, :-2, 1:-1]) * inv_dx
    dphi_dz = (phi[1:-1, 1:-1, 2:] - phi[1:-1, 1:-1, :-2]) * inv_dx
    grad2 = dphi_dx ** 2 + dphi_dy ** 2 + dphi_dz ** 2
    phi_int = phi[1:-1, 1:-1, 1:-1]
    potential = (1.0 / xi ** 2) * (1.0 - np.cos(phi_int))
    return float(np.sum(0.5 * grad2 + potential)) * dx ** 3


# ---------------------------------------------------------------------------
# String tension extraction
# ---------------------------------------------------------------------------

@dataclass
class ConfinementResult:
    """Results of the confinement potential calculation.

    Attributes:
        R_values:               Kink separations (lattice units).
        E_values:               Total energies at each R.
        sigma_lattice:          String tension (lattice units, dE/dR slope).
        alpha_C:                Coulomb-like 1/R coefficient.
        E_offset:               Constant energy offset.
        fit_quality:            R² of the σR + α_C/R + E₀ fit.
        sigma_SI:               String tension in J/m.
        sigma_GeV2:             String tension in GeV².
        proton_mass_MeV:        Predicted proton mass from Nambu-Goto formula.
        proton_to_electron_ratio: mp/me from the model.
        notes:                  Derivation notes.
    """

    R_values: list[float]
    E_values: list[float]
    sigma_lattice: float
    alpha_C: float
    E_offset: float
    fit_quality: float
    sigma_SI: float = 0.0
    sigma_GeV2: float = 0.0
    proton_mass_MeV: float = 0.0
    proton_to_electron_ratio: float = 0.0
    notes: list[str] = field(default_factory=list)


def fit_confinement_curve(
    R_arr: NDArray[np.float64],
    E_arr: NDArray[np.float64],
) -> tuple[float, float, float, float]:
    """Fit E(R) = σR + α_C/R + E₀ via linear least squares.

    Args:
        R_arr: Separations.
        E_arr: Energies.

    Returns:
        (sigma, alpha_C, E_offset, R_squared)
    """
    A = np.column_stack([R_arr, 1.0 / R_arr, np.ones_like(R_arr)])
    coeffs, _, _, _ = np.linalg.lstsq(A, E_arr, rcond=None)
    sigma, alpha_C, E0 = float(coeffs[0]), float(coeffs[1]), float(coeffs[2])
    E_fit = sigma * R_arr + alpha_C / R_arr + E0
    ss_res = float(np.sum((E_arr - E_fit) ** 2))
    ss_tot = float(np.sum((E_arr - float(np.mean(E_arr))) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return sigma, alpha_C, E0, r2


def scan_separations(
    R_list: list[int],
    params: LatticeParams,
    verbose: bool = True,
) -> ConfinementResult:
    """Compute E(R) for a range of kink separations and extract σ.

    For each R:
        1. Build the 3D field with the exact 1D kink–antikink profile.
        2. Relax to the static sine-Gordon minimum (with topological pinning).
        3. Integrate the energy functional.
        4. Fit E(R) = σR + α_C/R + E₀.

    Args:
        R_list:  Integer separations in lattice units.
        params:  Lattice and solver parameters.
        verbose: Progress printing.

    Returns:
        ConfinementResult with extracted σ.
    """
    N = params.N
    center = N // 2
    dx = params.dx
    xi = params.xi

    R_values: list[float] = []
    E_values: list[float] = []

    for R in R_list:
        if R < 2:
            continue
        z1 = center - R // 2
        z2 = center + R // 2
        if z1 < params.pin_rad + 2 or z2 > N - params.pin_rad - 3:
            if verbose:
                print(f"  R={R}: kinks too close to boundary — skipping.")
            continue

        if verbose:
            print(f"\n[R = {R:3d}]  kinks at z={z1}, z={z2}")

        phi0, phi_1d = build_initial_field(N, xi, z1, z2, center)
        phi_r, residuals, n_iter = relax_field(
            phi0, phi_1d, params, center, verbose=verbose
        )
        E = compute_field_energy(phi_r, dx, xi)

        if verbose:
            final_res = residuals[-1]
            print(f"  n_iter={n_iter}, final_res={final_res:.2e}, E={E:.5f}")

        R_values.append(float(R) * dx)
        E_values.append(E)

    R_arr = np.array(R_values)
    E_arr = np.array(E_values)
    sigma, alpha_C, E0, r2 = fit_confinement_curve(R_arr, E_arr)

    return ConfinementResult(
        R_values=list(R_arr),
        E_values=list(E_arr),
        sigma_lattice=sigma,
        alpha_C=alpha_C,
        E_offset=E0,
        fit_quality=r2,
        notes=[
            f"E(R) = {sigma:.4f}·R + {alpha_C:.4f}/R + {E0:.4f}",
            f"R² = {r2:.5f}",
        ],
    )


# ---------------------------------------------------------------------------
# SI conversion and physical predictions
# ---------------------------------------------------------------------------

def convert_to_si(
    result: ConfinementResult,
    K_SI: float,
    xi_SI: float,
) -> ConfinementResult:
    """Convert lattice σ to SI and predict physical observables.

    Unit derivation (careful):
        In natural units dx = ξ = 1, the energy integrand is
            ε_nat = ½|∇φ|² + (1/ξ²)(1 − cos φ)
        where gradients ∂φ/∂x are taken in lattice units, and ξ = 1.
        In physical units each term has a factor K_SI (energy density scale):
            ε_phys [J/m³] = K_SI × ε_nat
        The volume element is (ξ_SI)³, so:
            E_phys [J] = K_SI × ξ_SI³ × sum(ε_nat) = K_SI × ξ_SI × E_nat
        where E_nat = sum(ε_nat) × dx³ and dx = 1 (one lattice unit = ξ_SI).
        Then:
            σ_SI [J/m] = dE_phys/dR_phys
                       = (dE_nat/dR_nat) × (K_SI × ξ_SI) / ξ_SI
                       = σ_nat × K_SI
        The ξ_SI factors cancel exactly.  Hence:
            σ_SI = σ_lattice × K_SI

    Proton mass predictions:
        Method 2 (Nambu-Goto):  m_p c² = √(σ_SI × ℏc)
        Method 3 (3 kinks):     m_p c² = 3 × m_kink c²
                                m_kink [kg] = 8ℏ/(c·ξ),  m_kink c² [MeV] = 8m_e c²

    Args:
        result: ConfinementResult with lattice σ.
        K_SI:   Substrate stiffness [J/m³].
        xi_SI:  Coherence length [m].

    Returns:
        Updated ConfinementResult with SI fields and proton predictions.
    """
    # Correct unit conversion: σ_SI = σ_nat × K_SI  (ξ cancels in dE/dR)
    sigma_SI = result.sigma_lattice * K_SI
    hbar_c_SI = HBAR_SI * C_SI
    sigma_GeV2 = sigma_SI * hbar_c_SI / GEV_TO_J ** 2

    # Nambu-Goto / Regge: m_p c² = √(σ_SI × ℏc)  [J → MeV]
    mp_m2_J = math.sqrt(abs(sigma_SI) * hbar_c_SI)
    mp_m2_MeV = mp_m2_J / (EV_TO_J * 1e6)
    mp_m2_ratio = mp_m2_MeV / M_E_MEV

    # 3 × m_kink: m_kink [kg] = 8ℏ/(c·ξ), then multiply by c² to get energy
    m_kink_kg = 8.0 * HBAR_SI / (C_SI * xi_SI)         # kg
    m_kink_J_energy = m_kink_kg * C_SI ** 2             # J  (rest energy)
    m_kink_MeV = m_kink_J_energy / (EV_TO_J * 1e6)      # MeV
    mp_m3_MeV = 3.0 * m_kink_MeV

    # Best proton mass estimate: Nambu-Goto (physically motivated by string model)
    proton_mass_MeV = mp_m2_MeV
    proton_to_electron = proton_mass_MeV / M_E_MEV

    notes = list(result.notes) + [
        f"K_SI = {K_SI:.4e} J/m³,  xi_SI = {xi_SI:.4e} m",
        f"σ_SI = σ_nat × K_SI = {sigma_SI:.4e} J/m  [ξ cancels in dE/dR]",
        f"σ_GeV² = {sigma_GeV2:.4e}  (QCD lattice: {SIGMA_QCD_GEV2:.4f} GeV²)",
        f"σ ratio model/QCD = {sigma_GeV2 / SIGMA_QCD_GEV2:.3e}",
        f"m_p (Nambu-Goto √(σℏc)) = {mp_m2_MeV:.2f} MeV, mp/me = {mp_m2_ratio:.1f}",
        f"m_kink = 8ℏc/(c²ξ) = 8m_e = {m_kink_MeV:.4f} MeV  (ξ=λ_C(e))",
        f"m_p (3×m_kink) = {mp_m3_MeV:.4f} MeV",
    ]

    return ConfinementResult(
        R_values=result.R_values,
        E_values=result.E_values,
        sigma_lattice=result.sigma_lattice,
        alpha_C=result.alpha_C,
        E_offset=result.E_offset,
        fit_quality=result.fit_quality,
        sigma_SI=sigma_SI,
        sigma_GeV2=sigma_GeV2,
        proton_mass_MeV=proton_mass_MeV,
        proton_to_electron_ratio=proton_to_electron,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Flux tube geometry analysis
# ---------------------------------------------------------------------------

def analyse_flux_tube(
    phi: NDArray[np.float64],
    z1: int,
    z2: int,
    dx: float,
    xi: float,
) -> dict[str, object]:
    """Characterise the emergent flux tube geometry.

    Computes:
    - Energy density on the tube axis (x = y = center, varying z)
    - Transverse RMS width at the tube mid-point
    - Mid-axis uniformity (1 = perfectly constant energy density → string-like)

    Args:
        phi:  Relaxed field.
        z1:   First kink z-index.
        z2:   Second kink z-index.
        dx:   Lattice spacing.
        xi:   Coherence length.

    Returns:
        Dict with geometry metrics.
    """
    N = phi.shape[0]
    center = N // 2
    inv_dx = 1.0 / dx

    # Axis energy density profile
    axis_energy: list[float] = []
    for iz in range(N):
        ix = iy = center
        if any(v == 0 or v == N - 1 for v in [ix, iy, iz]):
            axis_energy.append(0.0)
            continue
        grad2 = (
            ((phi[ix + 1, iy, iz] - phi[ix - 1, iy, iz]) * 0.5 * inv_dx) ** 2
            + ((phi[ix, iy + 1, iz] - phi[ix, iy - 1, iz]) * 0.5 * inv_dx) ** 2
            + ((phi[ix, iy, iz + 1] - phi[ix, iy, iz - 1]) * 0.5 * inv_dx) ** 2
        )
        pot = (1.0 / xi ** 2) * (1.0 - math.cos(phi[ix, iy, iz]))
        axis_energy.append(0.5 * grad2 + pot)

    # Transverse profile at tube mid-point
    z_mid = (z1 + z2) // 2
    transverse_profile: list[tuple[float, float]] = []
    for ix in range(N):
        rho = abs(ix - center) * dx
        ix_m = max(ix - 1, 0)
        ix_p = min(ix + 1, N - 1)
        grad2_approx = (
            ((phi[ix_p, center, z_mid] - phi[ix_m, center, z_mid]) * 0.5 * inv_dx) ** 2
        )
        pot = (1.0 / xi ** 2) * (1.0 - math.cos(phi[ix, center, z_mid]))
        transverse_profile.append((rho, 0.5 * grad2_approx + pot))

    rho_arr = np.array([p[0] for p in transverse_profile])
    eps_arr = np.array([p[1] for p in transverse_profile])
    total_eps = float(np.sum(eps_arr))
    rms_width = (
        float(np.sqrt(np.sum(rho_arr ** 2 * eps_arr) / total_eps))
        if total_eps > 0
        else 0.0
    )

    # Mid-axis uniformity
    mid_start = z1 + 2
    mid_end = z2 - 2
    if mid_end > mid_start:
        mid_e = [axis_energy[z] for z in range(mid_start, mid_end)]
        mean_e = float(np.mean(mid_e))
        std_e = float(np.std(mid_e))
        uniformity = 1.0 - std_e / (mean_e + 1e-12)
    else:
        uniformity = 0.0

    return {
        "axis_energy_profile": axis_energy,
        "transverse_profile": transverse_profile,
        "tube_rms_width_lattice": rms_width / dx,
        "tube_rms_width_xi": rms_width / xi,
        "mid_axis_uniformity": uniformity,
        "is_string_like": uniformity > 0.7,
    }
