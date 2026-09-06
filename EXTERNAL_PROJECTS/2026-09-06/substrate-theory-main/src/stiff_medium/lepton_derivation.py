"""Lepton mass spectrum derivation from 3D Dirac equation in kink backgrounds.

Spec context (§18.30, §18.35, §18.48):
- §18.30: muon and tau are the same electron in stress-loaded excited states.
  2 free excitation energies (Δ_1, Δ_2), fewer than SM's 3 Yukawa couplings.
- §18.35: m c² = ℏ √(κ/I), where κ is the directional stiffness of the
  medium. Stress-loading redistributes κ across vertex quantum numbers.
- §18.48.2: empirical power-law κ_n ∝ (n+1)^k with k ≈ 15 fits both
  ratios within 4% individually, but gives Koide Q ≈ 0.91, not 2/3.

This module ATTEMPTS to derive lepton mass ratios from first principles
using three ansatzes:
  (a) Poschl-Teller-like m(r) = M tanh(r/xi) — radial 3D Dirac
  (b) Multi-kink superposition profiles (kink-antikink-kink chain)
  (c) Vertex-loaded profiles with localized stress bumps

Honest status (inherited from lepton_dirac_solver.py and multi_kink_dirac.py):
  - Single-kink Poschl-Teller: ratios bounded by √2 ≈ 1.41, far below 207.
  - Multi-kink chains: still O(few) — not O(200).
  - Stress-loaded bumps: bumps produce new localized states near the bump,
    not excited states of the same confinement well; ratios remain O(few).
  - Power-law κ_n ∝ (n+1)^k: reproduces EACH ratio at k ≈ 15 individually
    but violates Koide Q = 2/3 (gives Q ≈ 0.91).
  - No single profile found that simultaneously gives m_μ/m_e ≈ 207,
    m_τ/m_e ≈ 3477, and Koide Q = 2/3.

This module provides rigorous numerics for all three ansatzes plus the
§18.48.2 power-law search, so the caller can assess what the model can
and cannot derive.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional
import numpy as np
from scipy.linalg import eigh

# ---------------------------------------------------------------------------
# Physical targets (PDG 2024)
# ---------------------------------------------------------------------------

M_E_GEV: float = 0.000511  # electron mass GeV
M_MU_GEV: float = 0.105658  # muon mass GeV
M_TAU_GEV: float = 1.776860  # tau mass GeV

RATIO_MU_E: float = M_MU_GEV / M_E_GEV    # ≈ 206.77
RATIO_TAU_E: float = M_TAU_GEV / M_E_GEV  # ≈ 3476.82

# Koide constant Q = (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)²
# Measured to be 2/3 within 1e-5
KOIDE_TARGET: float = 2.0 / 3.0


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class BoundSpectrum:
    """Eigenvalues from a radial Dirac / Schrodinger-like solve.

    Attributes:
        energies: Sorted array of bound-state energies E_n, n=0,1,2,...
        ansatz: Human-readable label for the potential used.
        m_inf: Asymptotic mass parameter (in natural units).
        n_bound: Total number of bound states found.
    """

    energies: np.ndarray
    ansatz: str
    m_inf: float
    n_bound: int

    @property
    def ratio_10(self) -> Optional[float]:
        """E_1 / E_0 (analogous to m_μ / m_e if n=0 → electron)."""
        if self.n_bound >= 2:
            return float(self.energies[1] / self.energies[0])
        return None

    @property
    def ratio_20(self) -> Optional[float]:
        """E_2 / E_0 (analogous to m_τ / m_e if n=0 → electron)."""
        if self.n_bound >= 3:
            return float(self.energies[2] / self.energies[0])
        return None

    def koide_q(self) -> Optional[float]:
        """Koide invariant Q = (Σ E_n) / (Σ √E_n)² for the three lowest states.

        Returns None if fewer than 3 bound states exist.
        """
        if self.n_bound < 3:
            return None
        e = self.energies[:3]
        return float(np.sum(e) / np.sum(np.sqrt(e)) ** 2)

    def summary(self) -> str:
        """One-line human-readable summary."""
        parts = [f"ansatz={self.ansatz!r}", f"m_inf={self.m_inf:.3g}",
                 f"n_bound={self.n_bound}"]
        r10 = self.ratio_10
        r20 = self.ratio_20
        if r10 is not None:
            parts.append(f"r10={r10:.4g}")
        if r20 is not None:
            parts.append(f"r20={r20:.4g}")
        q = self.koide_q()
        if q is not None:
            parts.append(f"Q={q:.4f}")
        return "  ".join(parts)


# ---------------------------------------------------------------------------
# Core solver: squared Dirac Hamiltonian on 1D/radial grid
# ---------------------------------------------------------------------------

def _build_laplacian(N: int, dx: float) -> np.ndarray:
    """Build the 1D finite-difference second-derivative matrix (N×N).

    Uses three-point stencil with Dirichlet boundary conditions.

    Args:
        N: Grid points.
        dx: Grid spacing.

    Returns:
        Dense N×N matrix for -d²/dx².
    """
    diag = np.full(N, 2.0 / dx**2)
    off = np.full(N - 1, -1.0 / dx**2)
    mat = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    return mat


def solve_squared_dirac(
    mass_profile: np.ndarray,
    mass_prime: np.ndarray,
    m_inf: float,
    dx: float,
    n_states: int = 8,
) -> np.ndarray:
    """Solve the squared Dirac equation H ψ = E² ψ.

    The Dirac equation in background m(x) for the upper component ψ_+
    reduces (via squaring) to:

        H_+ = -∂_x² + m(x)² + m'(x)

    Bound states satisfy 0 < E² < m_inf².

    Args:
        mass_profile: Array of m(x) values on the grid.
        mass_prime: Array of m'(x) values on the grid.
        m_inf: Asymptotic |m(x → ±∞)|, used as the continuum threshold.
        dx: Grid spacing.
        n_states: Maximum bound states to return.

    Returns:
        Sorted array of bound-state energies E_n (not E²).
    """
    N = len(mass_profile)
    lap = _build_laplacian(N, dx)
    V = mass_profile**2 + mass_prime
    H = lap + np.diag(V)

    # Only request the lowest eigenvalues to keep cost manageable
    n_req = min(n_states + 20, N - 2)
    eigvals = eigh(H, subset_by_index=[0, n_req - 1], eigvals_only=True)
    eigvals = np.sort(eigvals)

    bound_mask = (eigvals > 1e-12) & (eigvals < m_inf**2)
    e_sq_bound = eigvals[bound_mask]
    return np.sqrt(e_sq_bound)[:n_states]


# ---------------------------------------------------------------------------
# Ansatz (a): Poschl-Teller-like tanh profile
# ---------------------------------------------------------------------------

def poschl_teller_profile(
    x: np.ndarray,
    m_inf: float,
    xi: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Poschl-Teller kink mass profile m(x) = m_inf * tanh(x/xi).

    Args:
        x: Grid coordinates.
        m_inf: Asymptotic mass.
        xi: Kink width (coherence length).

    Returns:
        Tuple (m, m_prime) of mass profile and its derivative.
    """
    m = m_inf * np.tanh(x / xi)
    m_prime = (m_inf / xi) / np.cosh(x / xi) ** 2
    return m, m_prime


def solve_poschl_teller(
    m_inf: float,
    xi: float = 1.0,
    L_box: float = 30.0,
    N: int = 801,
    n_states: int = 6,
) -> BoundSpectrum:
    """Solve the Dirac equation in a single Poschl-Teller kink potential.

    Per §18.11/§18.12: the electron corresponds to the lowest non-zero
    bound state (n=0 is the Jackiw-Rebbi zero-mode = neutrino). The next
    two excited states would be muon and tau if ratios worked out.

    Result (as already established in lepton_dirac_solver.py): the ratios
    are bounded by ~√2 for large m_inf*xi, far below the observed 207.

    Args:
        m_inf: Asymptotic mass parameter in natural units.
        xi: Kink width.
        L_box: Half-box size.
        N: Grid resolution.
        n_states: States to return.

    Returns:
        BoundSpectrum for this ansatz.
    """
    x = np.linspace(-L_box, L_box, N)
    dx = x[1] - x[0]
    m, mp = poschl_teller_profile(x, m_inf, xi)
    energies = solve_squared_dirac(m, mp, m_inf, dx, n_states)
    return BoundSpectrum(
        energies=energies,
        ansatz=f"poschl_teller(m_inf={m_inf:.3g},xi={xi:.3g})",
        m_inf=m_inf,
        n_bound=len(energies),
    )


def scan_poschl_teller(
    m_inf_values: Optional[list[float]] = None,
    xi: float = 1.0,
    L_box: float = 40.0,
    N: int = 801,
) -> list[BoundSpectrum]:
    """Scan over m_inf values for the Poschl-Teller ansatz.

    The parameter k = m_inf * xi controls the number of bound states and
    the ratio saturation. For large k the ratios saturate near √(2n+1).

    Args:
        m_inf_values: List of m_inf values. Defaults to broad sweep.
        xi: Kink width (fixed).
        L_box: Half-box size.
        N: Grid resolution.

    Returns:
        List of BoundSpectrum results.
    """
    if m_inf_values is None:
        m_inf_values = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
    return [
        solve_poschl_teller(m_inf=m, xi=xi, L_box=L_box, N=N)
        for m in m_inf_values
    ]


# ---------------------------------------------------------------------------
# Ansatz (b): Multi-kink superposition
# ---------------------------------------------------------------------------

def multi_kink_profile(
    x: np.ndarray,
    kink_positions: list[float],
    kink_signs: list[float],
    m_inf: float,
    xi: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Mass profile as superposition of signed tanh kinks.

    Each kink at x_k with sign s_k:
        m(x) = Σ_k  s_k * m_inf * tanh((x - x_k)/xi)

    m_inf is the amplitude per kink; the effective asymptote depends on
    the sign arrangement.

    Args:
        x: Grid coordinates.
        kink_positions: List of kink center positions.
        kink_signs: +1 (kink) or -1 (antikink) per element.
        m_inf: Mass amplitude per kink.
        xi: Kink width.

    Returns:
        Tuple (m, m_prime).
    """
    m = np.zeros_like(x)
    m_prime = np.zeros_like(x)
    for x_k, s_k in zip(kink_positions, kink_signs):
        u = (x - x_k) / xi
        m += s_k * m_inf * np.tanh(u)
        m_prime += s_k * (m_inf / xi) / np.cosh(u) ** 2
    return m, m_prime


@dataclass
class MultiKinkConfig:
    """Configuration for a multi-kink Dirac solve.

    Attributes:
        positions: Kink center positions.
        signs: +1 (kink) or -1 (antikink) for each center.
        m_inf: Mass amplitude per kink.
        xi: Kink width.
        label: Human-readable name.
    """

    positions: list[float]
    signs: list[float]
    m_inf: float
    xi: float
    label: str = "multi_kink"


def solve_multi_kink(
    config: MultiKinkConfig,
    L_box: float = 40.0,
    N: int = 1001,
    n_states: int = 8,
) -> BoundSpectrum:
    """Solve the Dirac equation in a multi-kink superposition potential.

    The effective asymptotic mass is estimated from the boundary values of
    |m(x)| at the grid edges.

    Args:
        config: MultiKinkConfig describing the kink arrangement.
        L_box: Half-box size.
        N: Grid resolution.
        n_states: States to return.

    Returns:
        BoundSpectrum for this configuration.
    """
    x = np.linspace(-L_box, L_box, N)
    dx = x[1] - x[0]
    m, mp = multi_kink_profile(x, config.positions, config.signs,
                                config.m_inf, config.xi)
    # Effective threshold: |m(boundary)| to catch both symmetric and
    # asymmetric configurations
    m_eff_inf = max(abs(float(m[0])), abs(float(m[-1])), 1e-12)
    energies = solve_squared_dirac(m, mp, m_eff_inf, dx, n_states)
    return BoundSpectrum(
        energies=energies,
        ansatz=config.label,
        m_inf=m_eff_inf,
        n_bound=len(energies),
    )


def standard_multi_kink_configs(m_inf: float = 6.0, xi: float = 1.0) -> list[MultiKinkConfig]:
    """Return a canonical set of multi-kink configurations to test.

    Includes kink-antikink (meson-like), kink-antikink-kink (baryon-like),
    and same-sign double kink, mirroring §18.30 stress-loading topologies.

    Args:
        m_inf: Mass amplitude per kink.
        xi: Kink width.

    Returns:
        List of MultiKinkConfig objects.
    """
    sep = 3.0 * xi
    return [
        MultiKinkConfig([-sep, sep], [+1, -1], m_inf, xi,
                        f"kink_antikink(sep={2*sep:.1f})"),
        MultiKinkConfig([-sep, 0.0, sep], [+1, -1, +1], m_inf, xi,
                        f"kink_antikink_kink(sep={sep:.1f})"),
        MultiKinkConfig([-sep, sep], [+1, +1], m_inf, xi,
                        f"double_kink(sep={2*sep:.1f})"),
        MultiKinkConfig([-2*sep, 0.0, 2*sep], [-1, +1, -1], m_inf, xi,
                        f"antikink_kink_antikink(sep={2*sep:.1f})"),
    ]


# ---------------------------------------------------------------------------
# Ansatz (c): Vertex-loaded profiles with localized stress
# ---------------------------------------------------------------------------

def vertex_loaded_profile(
    x: np.ndarray,
    m_inf: float,
    xi: float,
    stress_amplitudes: list[float],
    stress_positions: list[float],
    stress_sigma: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Single Poschl-Teller kink with Gaussian stress bumps (vertex quanta).

    Per §18.30: muon and tau arise from loading discrete stress quanta onto
    the vertex. This is modeled as additional Gaussian mass bumps on top of
    the base kink profile.

    m(x) = m_inf * tanh(x/xi) + Σ_i A_i * Gauss(x; x_i, sigma)

    Args:
        x: Grid coordinates.
        m_inf: Base kink asymptotic mass.
        xi: Base kink width.
        stress_amplitudes: Amplitude of each stress quantum bump.
        stress_positions: Position of each stress quantum.
        stress_sigma: Width of each Gaussian bump.

    Returns:
        Tuple (m, m_prime).
    """
    m = m_inf * np.tanh(x / xi)
    m_prime = (m_inf / xi) / np.cosh(x / xi) ** 2
    norm = stress_sigma * np.sqrt(2 * np.pi)
    for A, x_i in zip(stress_amplitudes, stress_positions):
        bump = A * np.exp(-0.5 * ((x - x_i) / stress_sigma) ** 2) / norm
        dbump = bump * (-(x - x_i) / stress_sigma ** 2)
        m += bump
        m_prime += dbump
    return m, m_prime


def solve_vertex_loaded(
    m_inf: float,
    xi: float,
    n_quanta: int,
    stress_amplitude: float,
    stress_sigma: float = 0.5,
    L_box: float = 40.0,
    N: int = 1001,
    n_states: int = 6,
) -> BoundSpectrum:
    """Solve Dirac with n_quanta stress bumps placed symmetrically near x=0.

    Per §18.30: n_quanta = 0 → electron, 1 → muon, 2 → tau.

    Stress bumps are placed at x = 0 for n=1, at x = ±xi for n=2,
    modeling the vertex angular-momentum redistribution.

    Args:
        m_inf: Asymptotic mass of the bare kink.
        xi: Kink width.
        n_quanta: Number of stress quanta (0, 1, or 2).
        stress_amplitude: Amplitude A per quantum bump.
        stress_sigma: Gaussian half-width for each bump.
        L_box: Half-box size.
        N: Grid resolution.
        n_states: Maximum bound states to return.

    Returns:
        BoundSpectrum for this configuration.
    """
    if n_quanta == 0:
        positions: list[float] = []
        amplitudes: list[float] = []
    elif n_quanta == 1:
        positions = [0.0]
        amplitudes = [stress_amplitude]
    else:  # n_quanta == 2
        positions = [-xi, xi]
        amplitudes = [stress_amplitude, stress_amplitude]

    x = np.linspace(-L_box, L_box, N)
    dx = x[1] - x[0]
    m, mp = vertex_loaded_profile(x, m_inf, xi, amplitudes, positions, stress_sigma)
    m_eff_inf = max(abs(float(m[0])), abs(float(m[-1])), 1e-12)
    energies = solve_squared_dirac(m, mp, m_eff_inf, dx, n_states)
    label = (f"vertex_loaded(m_inf={m_inf:.3g},A={stress_amplitude:.3g},"
             f"n={n_quanta})")
    return BoundSpectrum(
        energies=energies,
        ansatz=label,
        m_inf=m_eff_inf,
        n_bound=len(energies),
    )


def scan_vertex_loaded(
    m_inf: float = 4.0,
    xi: float = 1.0,
    stress_amplitudes: Optional[list[float]] = None,
    L_box: float = 40.0,
    N: int = 1001,
) -> list[dict]:
    """Scan stress amplitudes and extract three-state ratios.

    For each stress amplitude A, we solve n_quanta = 0, 1, 2 separately
    and record the lowest bound-state energy for each as a proxy for
    m_e, m_μ, m_τ excitation energies.

    Args:
        m_inf: Base kink mass.
        xi: Kink width.
        stress_amplitudes: Amplitudes to scan.
        L_box: Half-box size.
        N: Grid resolution.

    Returns:
        List of dicts with keys: stress_A, E0, E1, E2, ratio_10, ratio_20, Q.
    """
    if stress_amplitudes is None:
        stress_amplitudes = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    results = []
    for A in stress_amplitudes:
        spectra = [
            solve_vertex_loaded(m_inf, xi, nq, A, L_box=L_box, N=N, n_states=4)
            for nq in (0, 1, 2)
        ]
        E0 = spectra[0].energies[0] if spectra[0].n_bound >= 1 else None
        E1 = spectra[1].energies[0] if spectra[1].n_bound >= 1 else None
        E2 = spectra[2].energies[0] if spectra[2].n_bound >= 1 else None
        row: dict = {"stress_A": A, "E0": E0, "E1": E1, "E2": E2,
                     "ratio_10": None, "ratio_20": None, "Q": None}
        if E0 is not None and E1 is not None and E0 > 1e-15:
            row["ratio_10"] = float(E1 / E0)
        if E0 is not None and E2 is not None and E0 > 1e-15:
            row["ratio_20"] = float(E2 / E0)
        if E0 is not None and E1 is not None and E2 is not None:
            es = np.array([E0, E1, E2])
            q_denom = np.sum(np.sqrt(es)) ** 2
            row["Q"] = float(np.sum(es) / q_denom) if q_denom > 0 else None
        results.append(row)
    return results


# ---------------------------------------------------------------------------
# §18.48.2: Power-law κ model (cone-bouncing formulation)
# ---------------------------------------------------------------------------

def power_law_kappa_spectrum(
    k_exp: float,
    kappa_0: float = 1.0,
    I: float = 1.0,
    hbar: float = 1.0,
) -> tuple[np.ndarray, Optional[float]]:
    """Compute lepton masses from the §18.48.2 power-law stiffness model.

    Per §18.35: m_n c² = hbar * sqrt(kappa_n / I).
    Per §18.48.2 empirical fit: kappa_n = kappa_0 * (n+1)^k_exp.

    So m_n = hbar * sqrt(kappa_0/I) * (n+1)^(k_exp/2).

    This gives ratio m_1/m_0 = 2^(k_exp/2) and m_2/m_0 = 3^(k_exp/2).
    To match m_mu/m_e = 207: k_exp = 2*log(207)/log(2) ≈ 15.38
    To match m_tau/m_e = 3477: k_exp = 2*log(3477)/log(3) ≈ 14.84

    Args:
        k_exp: Power-law exponent k.
        kappa_0: Base stiffness for n=0 state.
        I: Effective inertia.
        hbar: Reduced Planck constant (natural units).

    Returns:
        Tuple of (masses_array[3], koide_Q). Q is None if masses are degenerate.
    """
    masses = np.array([
        hbar * np.sqrt(kappa_0 / I) * (n + 1) ** (k_exp / 2)
        for n in range(3)
    ])
    denom = np.sum(np.sqrt(masses)) ** 2
    if denom < 1e-30:
        return masses, None
    q = float(np.sum(masses) / denom)
    return masses, q


def scan_power_law_kappa(
    k_values: Optional[list[float]] = None,
) -> list[dict]:
    """Scan the power-law exponent k and report ratios + Koide Q.

    Args:
        k_values: List of k_exp values to scan.

    Returns:
        List of dicts: k_exp, ratio_10, ratio_20, Q, Q_deviation_from_2over3.
    """
    if k_values is None:
        k_values = [10.0, 12.0, 13.0, 14.0, 14.84, 15.0, 15.38, 16.0, 18.0]
    results = []
    for k in k_values:
        masses, q = power_law_kappa_spectrum(k_exp=k)
        r10 = float(masses[1] / masses[0]) if masses[0] > 0 else None
        r20 = float(masses[2] / masses[0]) if masses[0] > 0 else None
        q_dev = abs(q - KOIDE_TARGET) / KOIDE_TARGET if q is not None else None
        results.append({
            "k_exp": k,
            "ratio_10": r10,
            "ratio_20": r20,
            "Q": q,
            "Q_deviation_pct": 100.0 * q_dev if q_dev is not None else None,
        })
    return results


# ---------------------------------------------------------------------------
# Koide constraint analysis
# ---------------------------------------------------------------------------

def koide_from_ratios(r1: float, r2: float) -> float:
    """Compute Koide Q from mass ratios r1 = m1/m0, r2 = m2/m0.

    Q = (1 + r1 + r2) / (1 + sqrt(r1) + sqrt(r2))^2

    Args:
        r1: m_1 / m_0 ratio.
        r2: m_2 / m_0 ratio.

    Returns:
        Koide invariant Q.
    """
    numerator = 1.0 + r1 + r2
    denominator = (1.0 + np.sqrt(r1) + np.sqrt(r2)) ** 2
    if denominator < 1e-30:
        return float("nan")
    return float(numerator / denominator)


def koide_from_masses(m0: float, m1: float, m2: float) -> float:
    """Compute Koide Q from three absolute masses.

    Args:
        m0: Ground state mass.
        m1: First excitation mass.
        m2: Second excitation mass.

    Returns:
        Koide invariant Q.
    """
    total = m0 + m1 + m2
    sqrt_sum = np.sqrt(m0) + np.sqrt(m1) + np.sqrt(m2)
    if sqrt_sum**2 < 1e-30:
        return float("nan")
    return float(total / sqrt_sum**2)


def measured_koide() -> float:
    """Return Koide Q from PDG lepton masses (sanity check = 2/3)."""
    return koide_from_masses(M_E_GEV, M_MU_GEV, M_TAU_GEV)


# ---------------------------------------------------------------------------
# Best-match finder: what (r1, r2) pair minimizes distance to (207, 3477)
# while staying on the Koide surface Q = 2/3?
# ---------------------------------------------------------------------------

def koide_surface_ratios(r1: float) -> tuple[float, float]:
    """Given r1 = m1/m0, find r2 = m2/m0 such that Koide Q = 2/3.

    Derive from Q = 2/3, letting s1 = sqrt(r1), s2 = sqrt(r2):

        (1 + s1^2 + s2^2) / (1 + s1 + s2)^2 = 2/3

    Cross-multiplying and rearranging yields a quadratic in s2:

        s2^2 - (4*s1 + 4)*s2 + (s1^2 - 4*s1 + 1) = 0

    This has two positive roots (for physical r1 > 0), corresponding to
    the two branches of the Koide surface.  The LARGE root gives the
    physical lepton tau (r2 >> r1 >> 1); the small root is a low-mass
    satellite solution.

    Args:
        r1: m_1 / m_0 (must be > 0).

    Returns:
        Tuple (r2_large, r2_small) where r2_large > r1 > r2_small.
        Returns (nan, nan) if no real positive solutions exist.
    """
    if r1 <= 0:
        return float("nan"), float("nan")
    s1 = np.sqrt(r1)
    a_coeff = 1.0
    b_coeff = -(4.0 * s1 + 4.0)
    c_coeff = s1**2 - 4.0 * s1 + 1.0
    discriminant = b_coeff**2 - 4.0 * a_coeff * c_coeff
    if discriminant < 0:
        return float("nan"), float("nan")
    sqrt_disc = np.sqrt(discriminant)
    s2_large = (-b_coeff + sqrt_disc) / (2.0 * a_coeff)
    s2_small = (-b_coeff - sqrt_disc) / (2.0 * a_coeff)
    r2_large = float(s2_large**2) if s2_large >= 0 else float("nan")
    r2_small = float(s2_small**2) if s2_small >= 0 else float("nan")
    return r2_large, r2_small


def koide_surface_prediction() -> dict:
    """Given m_mu/m_e = 207, what does Koide demand for m_tau/m_e?

    Uses the corrected quadratic derivation in koide_surface_ratios().
    There are TWO solutions for r2 given r1; the large root is the
    physical one (tau >> muon >> electron).

    Args: none (uses PDG targets).

    Returns:
        Dict with keys:
          r1_input, r2_koide_large, r2_koide_small, r2_target,
          r2_error_pct, koide_internally_consistent,
          verification_Q_at_solution.
    """
    r2_large, r2_small = koide_surface_ratios(RATIO_MU_E)

    # Verify Q at the large root (should be exactly 2/3)
    q_check = koide_from_ratios(RATIO_MU_E, r2_large) if np.isfinite(r2_large) else float("nan")
    r2_error_pct = abs(r2_large - RATIO_TAU_E) / RATIO_TAU_E * 100.0 if np.isfinite(r2_large) else float("nan")

    # Reverse: use r2=3477 to predict r1 (by symmetry of the Koide relation,
    # the same quadratic applies with r1 <-> r2 role swap).
    r1_large, r1_small = koide_surface_ratios(RATIO_TAU_E)
    # The SMALL root gives the physical muon/electron ratio (r1 << r2)
    r1_pred = r1_small if np.isfinite(r1_small) else r1_large
    r1_error_pct = abs(r1_pred - RATIO_MU_E) / RATIO_MU_E * 100.0

    consistent = bool(r2_error_pct < 0.1 and r1_error_pct < 0.1)
    return {
        "given_r1_207_koide_r2_large": r2_large,
        "given_r1_207_koide_r2_small": r2_small,
        "target_r2_3477": RATIO_TAU_E,
        "r2_error_pct": r2_error_pct,
        "given_r2_3477_koide_r1_pred": r1_pred,
        "target_r1_207": RATIO_MU_E,
        "r1_error_pct": r1_error_pct,
        "verification_Q": q_check,
        "koide_internally_consistent": consistent,
    }


# ---------------------------------------------------------------------------
# High-level driver: run all ansatzes and collect results
# ---------------------------------------------------------------------------

def run_all_ansatzes(
    n_grid: int = 801,
    L_box: float = 40.0,
    verbose: bool = False,
) -> dict:
    """Run all three ansatzes and the power-law κ scan.

    Returns a dict with four keys, each containing a list of results:
      - 'poschl_teller': list of BoundSpectrum
      - 'multi_kink': list of BoundSpectrum
      - 'vertex_loaded': list of dicts (per stress-amplitude scan)
      - 'power_law_kappa': list of dicts (per k_exp value)
      - 'koide_surface': dict (Koide constraint cross-check)
      - 'best_ratios': dict summarising peak ratios achieved

    Args:
        n_grid: Grid resolution for eigenvalue solves.
        L_box: Half-box size.
        verbose: If True, print progress lines.

    Returns:
        Nested dict of all numerical results.
    """
    if verbose:
        print("[ansatz a] Poschl-Teller scan ...")
    pt_results = scan_poschl_teller(
        m_inf_values=[2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0],
        xi=1.0, L_box=L_box, N=n_grid,
    )

    if verbose:
        print("[ansatz b] Multi-kink scan ...")
    mk_configs = []
    for m_inf in [4.0, 8.0, 16.0]:
        mk_configs.extend(standard_multi_kink_configs(m_inf=m_inf, xi=1.0))
    mk_results = [
        solve_multi_kink(cfg, L_box=L_box, N=n_grid)
        for cfg in mk_configs
    ]

    if verbose:
        print("[ansatz c] Vertex-loaded scan ...")
    vl_results = scan_vertex_loaded(
        m_inf=4.0, xi=1.0,
        stress_amplitudes=[0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0],
        L_box=L_box, N=n_grid,
    )

    if verbose:
        print("[§18.48.2] Power-law κ scan ...")
    kl_results = scan_power_law_kappa()

    if verbose:
        print("[Koide surface] Cross-check ...")
    koide_surface = koide_surface_prediction()

    # Collect best ratios achieved
    all_r10: list[float] = []
    all_r20: list[float] = []

    for sp in pt_results + mk_results:
        if sp.ratio_10 is not None:
            all_r10.append(sp.ratio_10)
        if sp.ratio_20 is not None:
            all_r20.append(sp.ratio_20)

    for row in vl_results:
        if row["ratio_10"] is not None:
            all_r10.append(row["ratio_10"])
        if row["ratio_20"] is not None:
            all_r20.append(row["ratio_20"])

    best_ratio_10 = max(all_r10) if all_r10 else None
    best_ratio_20 = max(all_r20) if all_r20 else None

    return {
        "poschl_teller": pt_results,
        "multi_kink": mk_results,
        "vertex_loaded": vl_results,
        "power_law_kappa": kl_results,
        "koide_surface": koide_surface,
        "best_ratios": {
            "best_r10_dirac": best_ratio_10,
            "best_r20_dirac": best_ratio_20,
            "target_r10": RATIO_MU_E,
            "target_r20": RATIO_TAU_E,
        },
    }
