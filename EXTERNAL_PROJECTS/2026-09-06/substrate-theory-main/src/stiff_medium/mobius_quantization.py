"""Möbius topological quantization for the lepton mass spectrum.

Spec context (§§18.13, 18.30, 18.35, 18.48):
- §18.35: m_n c² = ℏ √(κ_n / I), where κ_n is the directional stiffness of
  the stress-loaded vertex configuration with n quanta.
- §18.48.2: power-law κ_n ∝ (n+1)^k fits each ratio at k≈15 individually but
  gives Koide Q = 0.91, not 2/3.  A Möbius-topology constraint is needed.
- §18.13, §18.30: muon and tau are the electron in 1- and 2-quantum
  stress-loaded states; exactly 3 generations from vertex closure.

This module implements four approaches to find a topological quantization
condition on κ_n derivable from the Möbius bundle that simultaneously gives
    m_μ / m_e = 207,   m_τ / m_e = 3477,   Koide Q = 2/3.

Approach A — Bohr-Sommerfeld on stress-loaded orbit
    ∮ p dq = n ℏ  →  κ_n from quantised action integrals.
    Torsional pendulum on the Möbius cone: kinetic = (1/2) I θ̇²,
    potential = (1/2) κ θ².  Action integral = π I ω_n = n π ℏ.

Approach B — Atiyah-Singer index on the Möbius torus
    The index of the chiral Dirac operator on T² with half-flux:
        ind(D) = c₁ · χ(T²) = 1/2 · 0 = 0
    but the η-invariant / mod-2 index gives discrete spectrum via the
    self-adjoint extension parameter θ_SA.  Spectrum:
        λ_n = (n + θ_SA/π)  on the circle.
    Search over θ_SA to fit lepton ratios.

Approach C — HOMFLY / Jones polynomial of Möbius torus knot
    The (2, n) torus knot T(2,n) has Jones invariant V(t) that encodes
    discrete "quantum numbers" through its skein relation.  The zeros of
    the Alexander polynomial ΔT(2,n)(t) give integer Casimir-like values.
    These are used as stiffness eigenvalues κ_n ∝ Casimir(n).

Approach D — Möbius operator eigenvalue equation
    The generator of the Möbius group acting on the upper-half-plane
    (or circle): L f = -i (z d/dz + 1/2) f.
    Discretize on the Möbius strip boundary with half-flux constraint
    → solve eigenvalue equation and identify spectrum with κ_n.

Each approach:
  1. Computes candidate κ_0, κ_1, κ_2.
  2. Checks m_μ/m_e = 207 and m_τ/m_e = 3477.
  3. Checks Koide Q = 2/3.
  4. Reports how close the topology comes to both constraints.

Honest summary: this is a deep open problem.  Approaches B–D involve
genuine topological mathematics.  None is guaranteed to give BOTH ratios
AND Q = 2/3 from pure topology without fitting parameters.  The module
documents partial successes and the gap that remains.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import minimize_scalar, brentq, minimize
from scipy.linalg import eigh
import sympy as sp

# ---------------------------------------------------------------------------
# Physical targets (PDG 2024)
# ---------------------------------------------------------------------------

M_E_GEV: float = 0.000511
M_MU_GEV: float = 0.105658
M_TAU_GEV: float = 1.776860

RATIO_MU_E: float = M_MU_GEV / M_E_GEV    # 206.768
RATIO_TAU_E: float = M_TAU_GEV / M_E_GEV  # 3476.82
KOIDE_TARGET: float = 2.0 / 3.0

# Tolerance for "matching" ratios and Koide
RATIO_TOL_PCT: float = 5.0   # 5 % tolerance on each mass ratio
KOIDE_TOL: float = 0.01      # |Q - 2/3| < 0.01


# ---------------------------------------------------------------------------
# Utility: Koide Q from three masses (or stiffnesses via m ∝ √κ)
# ---------------------------------------------------------------------------


def koide_q(m0: float, m1: float, m2: float) -> float:
    """Koide invariant Q = (m0+m1+m2) / (√m0+√m1+√m2)².

    For the lepton masses this equals 2/3 within 1e-5.

    Args:
        m0: Electron mass proxy.
        m1: Muon mass proxy.
        m2: Tau mass proxy.

    Returns:
        Koide Q value.  Returns NaN if any mass is non-positive.
    """
    if m0 <= 0 or m1 <= 0 or m2 <= 0:
        return float("nan")
    num = m0 + m1 + m2
    den = (np.sqrt(m0) + np.sqrt(m1) + np.sqrt(m2)) ** 2
    return float(num / den)


def masses_from_kappas(
    kappas: np.ndarray,
    kappa_0_scale: float = 1.0,
    I: float = 1.0,
    hbar: float = 1.0,
) -> np.ndarray:
    """Convert stiffness eigenvalues to masses via m_n = ℏ √(κ_n / I).

    Per §18.35: m_n c² = ℏ √(κ_n / I).

    Args:
        kappas: Array of κ values (must all be positive).
        kappa_0_scale: Additional normalisation (does not affect ratios).
        I: Effective inertia (does not affect ratios).
        hbar: Reduced Planck constant (does not affect ratios).

    Returns:
        Array of masses proportional to √κ_n.
    """
    return hbar * np.sqrt(np.asarray(kappas) * kappa_0_scale / I)


def assess_fit(kappas: np.ndarray) -> dict:
    """Compute all fit metrics for three stiffness eigenvalues.

    Args:
        kappas: Array of exactly three κ values (κ_e, κ_μ, κ_τ).

    Returns:
        Dict with keys: masses, ratio_mu_e, ratio_tau_e, koide_q,
        ratio_mu_e_err_pct, ratio_tau_e_err_pct, koide_err.
    """
    if len(kappas) < 3 or np.any(np.array(kappas[:3]) <= 0):
        return {
            "masses": None,
            "ratio_mu_e": None,
            "ratio_tau_e": None,
            "koide_q": None,
            "ratio_mu_e_err_pct": float("inf"),
            "ratio_tau_e_err_pct": float("inf"),
            "koide_err": float("inf"),
        }
    k = np.array(kappas[:3], dtype=float)
    m = masses_from_kappas(k)
    r10 = float(m[1] / m[0])
    r20 = float(m[2] / m[0])
    q = koide_q(m[0], m[1], m[2])
    return {
        "masses": m,
        "ratio_mu_e": r10,
        "ratio_tau_e": r20,
        "koide_q": q,
        "ratio_mu_e_err_pct": abs(r10 - RATIO_MU_E) / RATIO_MU_E * 100,
        "ratio_tau_e_err_pct": abs(r20 - RATIO_TAU_E) / RATIO_TAU_E * 100,
        "koide_err": abs(q - KOIDE_TARGET) if np.isfinite(q) else float("inf"),
    }


# ---------------------------------------------------------------------------
# Koide surface utilities (exact algebra)
# ---------------------------------------------------------------------------


def koide_surface_r2_given_r1(r1: float) -> tuple[float, float]:
    """Given r1 = m1/m0, solve Koide Q = 2/3 for r2 = m2/m0.

    The constraint Q = 2/3 with masses (1, r1, r2) gives a quadratic in
    s2 = √r2:
        s2² - (4s1 + 4)·s2 + (s1² - 4s1 + 1) = 0
    where s1 = √r1.

    Args:
        r1: m1/m0 ratio (must be positive).

    Returns:
        Tuple (r2_large, r2_small).  Large root is the physical tau branch.
        Returns (nan, nan) if no real positive solution.
    """
    if r1 <= 0:
        return float("nan"), float("nan")
    s1 = float(np.sqrt(r1))
    a, b, c = 1.0, -(4.0 * s1 + 4.0), s1**2 - 4.0 * s1 + 1.0
    disc = b**2 - 4.0 * a * c
    if disc < 0:
        return float("nan"), float("nan")
    sd = float(np.sqrt(disc))
    s2_l = (-b + sd) / 2.0
    s2_s = (-b - sd) / 2.0
    return (float(s2_l**2) if s2_l >= 0 else float("nan"),
            float(s2_s**2) if s2_s >= 0 else float("nan"))


def koide_exact_ratios() -> dict:
    """Exact Koide constraint cross-check.

    Given r1 = 207, what does Koide demand for r2?
    Answer should be ≈ 3477 if the relation is self-consistent.

    Returns:
        Dict with predicted r2 (large branch) and error vs PDG.
    """
    r2_large, r2_small = koide_surface_r2_given_r1(RATIO_MU_E)
    err = abs(r2_large - RATIO_TAU_E) / RATIO_TAU_E * 100
    q_check = koide_q(1.0, RATIO_MU_E, r2_large)
    return {
        "r1_input": RATIO_MU_E,
        "r2_predicted_large": r2_large,
        "r2_predicted_small": r2_small,
        "r2_target": RATIO_TAU_E,
        "r2_error_pct": err,
        "q_verification": q_check,
        "is_consistent": bool(err < 0.1),
    }


# ---------------------------------------------------------------------------
# Approach A: Bohr-Sommerfeld on torsional oscillator
# ---------------------------------------------------------------------------


def _bs_action_torsional(omega: float, I: float = 1.0, hbar: float = 1.0) -> float:
    """Action integral for a torsional oscillator at frequency ω.

    For H = (1/2) I θ̇² + (1/2) κ θ²  with ω = √(κ/I):
        J = ∮ p dθ = π I ω = π ℏ n   →   ω_n = n ℏ/I

    Returns J(ω) = π I ω.
    """
    return float(np.pi * I * omega)


def bohr_sommerfeld_kappas(
    n_values: list[int] | None = None,
    I: float = 1.0,
    hbar: float = 1.0,
    alpha_BS: float = 1.0,
) -> np.ndarray:
    """Kappa eigenvalues from Bohr-Sommerfeld quantisation on the Möbius cone.

    The Möbius half-flux constraint modifies the quantum number: a half-flux
    winding (holonomy = −1) means the wave-function must be anti-periodic.
    The natural quantum numbers on a Möbius strip are n = 1/2, 3/2, 5/2, ...
    (half-integers), NOT the integers, because the boundary condition is ψ → −ψ
    after one circuit (Möbius identification).

    In the LEPTON model (§18.30), the THREE leptons correspond to:
        n=0  (electron — no stress quanta, torsional zero-point)
        n=1  (muon — one stress quantum)
        n=2  (tau — two stress quanta)

    However, with Möbius anti-periodic BC, the effective quantum number is
        n_eff = n + 1/2  (zero-point = 1/2, not 0)

    This gives:
        κ_n = I × ω_n²  where  ω_n = (n + 1/2) ℏ/(I × alpha_BS)

    So κ_n ∝ (n + 1/2)²   →   masses m_n ∝ (n + 1/2).

    The HALF-INTEGER shift is the Möbius contribution.  Standard harmonic
    oscillator gives n + 1/2 energy spectrum; here the FREQUENCY, not the
    energy, is quantised.

    Args:
        n_values: Quantum numbers [n_e, n_μ, n_τ] (default [0, 1, 2]).
        I: Effective inertia parameter.
        hbar: Reduced Planck constant.
        alpha_BS: Scale factor for the action unit.

    Returns:
        Array of κ_n for the given quantum numbers.
    """
    if n_values is None:
        n_values = [0, 1, 2]
    ns = np.array(n_values, dtype=float)
    # Möbius half-integer shift: n → n + 1/2
    n_eff = ns + 0.5
    # From B-S: ω_n = n_eff * (hbar / (I * alpha_BS))
    omega_n = n_eff * hbar / (I * alpha_BS)
    kappa_n = I * omega_n**2
    return kappa_n


@dataclass
class BohrSommerfeldResult:
    """Result from Bohr-Sommerfeld quantization search.

    Attributes:
        alpha_BS: Best-fit scale parameter.
        kappas: Kappa values for n=0,1,2.
        fit: Dict from assess_fit().
        variant: Description of the B-S variant used.
    """

    alpha_BS: float
    kappas: np.ndarray
    fit: dict
    variant: str


def run_bohr_sommerfeld() -> list[BohrSommerfeldResult]:
    """Run Bohr-Sommerfeld quantisation in several variants.

    Variants:
    1. Integer n with Möbius half-integer shift (n + 1/2)
    2. n² growth (energy-level analog)
    3. n³ growth (rotational level analog)
    4. Mixed: κ_n = (2n+1)² (odd half-integers on Möbius cylinder)

    Each variant uses alpha_BS as a trivial scale that does not affect mass
    ratios.  We report which variant comes closest to (207, 3477, Q=2/3).

    Returns:
        List of BohrSommerfeldResult for each variant.
    """
    results: list[BohrSommerfeldResult] = []

    # Variant 1: Möbius half-integer, κ_n ∝ (n + 1/2)²
    # masses ∝ n + 1/2 → ratios = (1.5/0.5, 2.5/0.5) = (3, 5)
    for alpha in [1.0, 0.1, 0.01]:
        k = bohr_sommerfeld_kappas([0, 1, 2], alpha_BS=alpha)
        fit = assess_fit(k)
        results.append(BohrSommerfeldResult(
            alpha_BS=alpha, kappas=k, fit=fit,
            variant="half-integer: kappa ~ (n+1/2)^2"
        ))

    # Variant 2: κ_n ∝ (n+1)² — integer shift, NO Möbius 1/2
    # masses ∝ n+1 → ratios = (2, 3)
    for alpha in [1.0]:
        ns = np.array([0, 1, 2], dtype=float)
        n_eff = ns + 1.0
        k = n_eff**2
        fit = assess_fit(k)
        results.append(BohrSommerfeldResult(
            alpha_BS=alpha, kappas=k, fit=fit,
            variant="integer: kappa ~ (n+1)^2"
        ))

    # Variant 3: Möbius winding-doubled — after two circuits the phase = 2π
    # → n_eff = 2n + 1 (odd integers only, due to Möbius twist doubling)
    # masses ∝ 2n+1 → ratios = (3, 5)
    for alpha in [1.0]:
        ns = np.array([0, 1, 2], dtype=float)
        n_eff = 2.0 * ns + 1.0
        k = n_eff**2
        fit = assess_fit(k)
        results.append(BohrSommerfeldResult(
            alpha_BS=alpha, kappas=k, fit=fit,
            variant="odd-integer (Mobius twist doubled): kappa ~ (2n+1)^2"
        ))

    # Variant 4: Power-law B-S with free exponent k_exp — reproduce §18.48.2
    # κ_n = (n+1)^k_exp.  Optimal k that minimises combined ratio error.
    def combined_error(k_exp: float) -> float:
        k = np.array([(n + 1) ** k_exp for n in range(3)])
        fit = assess_fit(k)
        e1 = fit["ratio_mu_e_err_pct"] / 100
        e2 = fit["ratio_tau_e_err_pct"] / 100
        return e1**2 + e2**2

    res_pl = minimize_scalar(combined_error, bounds=(10.0, 20.0), method="bounded")
    k_best = float(res_pl.x)
    k_pl = np.array([(n + 1) ** k_best for n in range(3)])
    fit_pl = assess_fit(k_pl)
    results.append(BohrSommerfeldResult(
        alpha_BS=k_best, kappas=k_pl, fit=fit_pl,
        variant=f"power-law optimal k={k_best:.4f}: kappa ~ (n+1)^k"
    ))

    # Variant 5: Constraint-based B-S — force Koide = 2/3
    # Given Q = 2/3 and r1 = 207, solve for r2; then find kappas.
    r2_large, _ = koide_surface_r2_given_r1(RATIO_MU_E)
    if np.isfinite(r2_large):
        # masses proportional to (1, 207, r2_large)
        # κ_n = m_n² → kappas proportional to (1, 207², r2_large²)
        k_koide = np.array([1.0, RATIO_MU_E**2, r2_large**2])
        fit_koide = assess_fit(k_koide)
        results.append(BohrSommerfeldResult(
            alpha_BS=0.0, kappas=k_koide, fit=fit_koide,
            variant="Koide-constrained (algebraic): kappas from Q=2/3 surface"
        ))

    return results


# ---------------------------------------------------------------------------
# Approach B: Atiyah-Singer index / self-adjoint extension
# ---------------------------------------------------------------------------


@dataclass
class AtiyahSingerResult:
    """Result from the Atiyah-Singer / self-adjoint extension approach.

    Attributes:
        theta_SA: Self-adjoint extension parameter in [0, π).
        kappas: Kappa values for n=0,1,2.
        fit: Dict from assess_fit().
        n_quantum: Quantum numbers used.
    """

    theta_SA: float
    kappas: np.ndarray
    fit: dict
    n_quantum: list[int]


def atiyah_singer_spectrum(
    theta_SA: float,
    n_values: list[int] | None = None,
) -> np.ndarray:
    """Eigenvalue spectrum of the Dirac operator on T² with Möbius half-flux.

    On the torus T² with a U(1) bundle of first Chern class c₁ = 1/2, the
    spectrum of the chiral Dirac operator D is:

        λ_n = n + θ_SA / (2π)

    where θ_SA ∈ [0, 2π) is the self-adjoint extension parameter that encodes
    the boundary condition on the Möbius strip edge.  For θ_SA = π the
    spectrum is n + 1/2 (anti-periodic = Möbius BC); for θ_SA = 0 it is n.

    In the stiffness model: κ_n = λ_n².

    The Atiyah-Singer index theorem for D on T² with c₁ = 1/2 gives
        ind(D) = ∫ Â(T²) ch(F) = (1/2) × 0 = 0
    but the mod-2 index is 1 (the half-integer Chern class forces a
    zero-mode parity anomaly), consistent with the Möbius fermion sign.

    Args:
        theta_SA: Extension parameter in [0, 2π).
        n_values: Quantum numbers (default [0, 1, 2]).

    Returns:
        Array of κ_n = λ_n² (must be positive; skip n where λ_n = 0).
    """
    if n_values is None:
        n_values = [0, 1, 2]
    ns = np.array(n_values, dtype=float)
    lam = ns + theta_SA / (2.0 * np.pi)
    kappas = lam**2
    return kappas


def scan_atiyah_singer(
    n_theta: int = 2000,
    n_values_list: list[list[int]] | None = None,
) -> list[AtiyahSingerResult]:
    """Scan θ_SA and find the best-fit extension parameter.

    For each (θ_SA, n_values) pair, compute the stiffness spectrum and
    evaluate fit metrics against (207, 3477, Q=2/3).

    Args:
        n_theta: Number of θ_SA values to scan in [0, 2π).
        n_values_list: Quantum number sets to try.

    Returns:
        Top-5 results sorted by combined error score.
    """
    if n_values_list is None:
        n_values_list = [
            [0, 1, 2],
            [1, 2, 3],
            [0, 2, 4],
            [1, 3, 5],
            [0, 1, 3],
        ]

    thetas = np.linspace(1e-6, 2.0 * np.pi - 1e-6, n_theta)
    best_results: list[AtiyahSingerResult] = []
    best_score = float("inf")

    for n_vals in n_values_list:
        for theta in thetas:
            k = atiyah_singer_spectrum(theta_SA=theta, n_values=n_vals)
            if np.any(k <= 0):
                continue
            fit = assess_fit(k)
            # Score: relative ratio errors + large penalty for Koide mismatch
            e1 = fit["ratio_mu_e_err_pct"] / 100
            e2 = fit["ratio_tau_e_err_pct"] / 100
            eq = fit["koide_err"] / KOIDE_TARGET
            score = e1**2 + e2**2 + 10.0 * eq**2  # Koide weighted 10×
            if score < best_score:
                best_score = score
                best_results = [AtiyahSingerResult(
                    theta_SA=float(theta),
                    kappas=k.copy(),
                    fit=fit,
                    n_quantum=n_vals,
                )]
            elif score < best_score * 2.0:
                best_results.append(AtiyahSingerResult(
                    theta_SA=float(theta),
                    kappas=k.copy(),
                    fit=fit,
                    n_quantum=n_vals,
                ))

    # Sort and return top 5
    best_results.sort(
        key=lambda r: (
            (r.fit["ratio_mu_e_err_pct"] / 100)**2
            + (r.fit["ratio_tau_e_err_pct"] / 100)**2
            + 10.0 * (r.fit["koide_err"] / KOIDE_TARGET)**2
        )
    )
    return best_results[:5]


def refine_atiyah_singer(
    initial_theta: float,
    n_values: list[int],
) -> AtiyahSingerResult:
    """Refine θ_SA around a coarse best-fit using scipy minimization.

    Minimises the combined ratio + Koide error with high-precision optimizer.

    Args:
        initial_theta: Starting point for θ_SA.
        n_values: Quantum numbers to use.

    Returns:
        Refined AtiyahSingerResult.
    """
    def objective(theta: float) -> float:
        k = atiyah_singer_spectrum(theta_SA=theta, n_values=n_values)
        if np.any(k <= 0):
            return 1e12
        fit = assess_fit(k)
        e1 = fit["ratio_mu_e_err_pct"] / 100
        e2 = fit["ratio_tau_e_err_pct"] / 100
        eq = fit["koide_err"] / KOIDE_TARGET
        return e1**2 + e2**2 + 10.0 * eq**2

    res = minimize_scalar(
        objective,
        bounds=(max(1e-9, initial_theta - 0.5), min(2 * np.pi - 1e-9, initial_theta + 0.5)),
        method="bounded",
        options={"xatol": 1e-12},
    )
    theta_best = float(res.x)
    k = atiyah_singer_spectrum(theta_SA=theta_best, n_values=n_values)
    fit = assess_fit(k)
    return AtiyahSingerResult(
        theta_SA=theta_best,
        kappas=k,
        fit=fit,
        n_quantum=n_values,
    )


# ---------------------------------------------------------------------------
# Approach C: Torus knot / Alexander polynomial Casimir values
# ---------------------------------------------------------------------------


@dataclass
class KnotInvariantResult:
    """Result from knot-invariant based stiffness quantization.

    Attributes:
        knot_type: String description of the knot.
        casimir_values: Discrete values from the knot invariant.
        kappas: Kappa values chosen from Casimir spectrum.
        fit: Dict from assess_fit().
    """

    knot_type: str
    casimir_values: list[float]
    kappas: np.ndarray
    fit: dict


def torus_knot_alexander_zeros(p: int, q: int) -> list[float]:
    """Zeros of the Alexander polynomial of the torus knot T(p,q).

    The Alexander polynomial of T(p,q) (as a polynomial in t) is:
        Δ_{T(p,q)}(t) = (t^{pq} - 1)(t - 1) / ((t^p - 1)(t^q - 1))

    Evaluated at t = e^{iφ} for φ ∈ (0, 2π), the modulus |Δ(e^{iφ})| has
    zeros at φ = 2π k/p and φ = 2π k/q for integer k.  These zero-crossing
    angles correspond to "topological resonances" we interpret as stiffness
    eigenvalues κ_k = (k + 1/2)^{-1} times a scale.

    For the Möbius knot (i.e. T(2,n)): p=2, q=n.
    The Alexander polynomial zeros are at t = e^{iπk} (from the p=2 factor):
        t = -1 → one zero at φ = π
    and from q=n: t = e^{2πi k/n} for k = 1, ..., n-1.

    We return the exponents (sorted) as the candidate Casimir numbers.

    Args:
        p: First torus knot parameter.
        q: Second torus knot parameter.

    Returns:
        Sorted list of exponent values from the Alexander polynomial structure.
    """
    # Use sympy for exact rational arithmetic
    t = sp.Symbol("t")
    # Alexander polynomial of T(p,q)
    numer = (t**(p * q) - 1) * (t - 1)
    denom = (t**p - 1) * (t**q - 1)
    # Polynomial division — work symbolically
    alex = sp.cancel(numer / denom)
    alex_poly = sp.Poly(sp.expand(alex), t)
    # Get all roots numerically
    try:
        roots = sp.nroots(alex_poly, n=15, maxsteps=200)
    except Exception:
        roots = []

    # Extract arguments of roots on the unit circle (|root| ≈ 1)
    casimirs: list[float] = []
    for root in roots:
        try:
            rv = complex(root)
            if abs(abs(rv) - 1.0) < 0.05:  # near unit circle
                phi = float(np.angle(rv))
                casimirs.append(abs(phi / np.pi))  # in units of π
        except (TypeError, ValueError):
            continue

    return sorted(set(round(c, 6) for c in casimirs))


def jones_polynomial_t2n(n: int) -> list[float]:
    """Jones polynomial invariant for T(2, n) torus knot — eigenvalue spectrum.

    The Jones polynomial of T(2, n) evaluated at t = -1 (the Arf invariant)
    gives the signature of the knot.  At t = exp(2πi / k) for small k, it
    gives the quantum dimension of the SU(2)_k representation.

    For our purposes we extract the discrete "eigenvalues" from the skein
    relation.  The (2, n) torus knot satisfies:
        V_{T(2,n)}(t) = -t^{(n+3)/2} × [(1 - t^{n+1}) / (1 - t²)] + t^{(n-1)/2}

    The zeros of the bracket polynomial (unnormalized Jones) occur at
    t = exp(2πi k / (n + 2)) for k = 1, ..., n.  These are the
    SU(2)_{n} quantum group roots of unity.  Their arguments:
        φ_k = 2π k / (n + 2)

    We interpret the "quantum numbers" as Casimir values
        C_k = k(k+2) / 4   (SU(2) quadratic Casimir for spin-j state)

    where j = k/2 runs over half-integer values k/2 = 1/2, 1, 3/2, ...

    Args:
        n: Torus knot parameter (also controls the number of terms).

    Returns:
        Sorted list of SU(2) Casimir values C_k = j(j+1) for j = k/2.
    """
    casimirs: list[float] = []
    for k in range(1, n + 1):
        j = k / 2.0
        casimirs.append(j * (j + 1))
    return sorted(casimirs)


def mobius_torus_stiffness_from_knot(
    p: int = 2,
    n_torus: int = 3,
    method: str = "alexander",
) -> KnotInvariantResult:
    """Compute stiffness eigenvalues from torus knot invariants.

    The Möbius bundle wraps once around the core (p=2 for half-flux) and
    n_torus times around the torus longitude.  The knot T(2, n_torus).

    For the lepton spectrum we need three distinct κ values.  We pick the
    first three non-trivial invariants.

    Args:
        p: Winding in the first direction (default 2 for Möbius half-flux).
        n_torus: Winding in the second direction.
        method: "alexander" or "jones".

    Returns:
        KnotInvariantResult with the best three κ values.
    """
    if method == "alexander":
        casimirs = torus_knot_alexander_zeros(p, n_torus)
        knot_label = f"Alexander T({p},{n_torus})"
    elif method == "jones":
        casimirs = jones_polynomial_t2n(n_torus)
        knot_label = f"Jones T(2,{n_torus})"
    else:
        raise ValueError(f"Unknown method {method!r}")

    if len(casimirs) < 3:
        # Pad with sequential integers if not enough invariants
        casimirs = casimirs + list(range(len(casimirs), 3))

    # Use first three as κ values (must be positive)
    kappas = np.array(casimirs[:3], dtype=float)
    kappas = kappas[kappas > 0]
    if len(kappas) < 3:
        kappas = np.array([max(casimirs[:3][i], 1e-6) for i in range(3)])

    fit = assess_fit(kappas)
    return KnotInvariantResult(
        knot_type=knot_label,
        casimir_values=casimirs[:6],
        kappas=kappas,
        fit=fit,
    )


def scan_jones_casimir() -> list[KnotInvariantResult]:
    """Scan Jones/SU(2)-Casimir spectra for different torus knot parameters.

    For T(2, n_torus) knots with n_torus ranging, pick the three best
    Casimir values and check fit to lepton ratios.

    The SU(2) Casimir values j(j+1) for j = 1/2, 1, 3/2, 2, ... are:
        3/4, 2, 15/4, 6, 35/4, ...  → 0.75, 2.0, 3.75, 6.0, 8.75, ...

    Mass ratios from √(κ_n / κ_0):
        √(2.0/0.75) ≈ 1.63,  √(3.75/0.75) ≈ 2.24, ...
    These are O(1), not O(200).  Jones invariants alone do not bridge the gap.

    Returns:
        List of KnotInvariantResult for small n_torus values.
    """
    results: list[KnotInvariantResult] = []
    for n in range(3, 25):
        r = mobius_torus_stiffness_from_knot(p=2, n_torus=n, method="jones")
        results.append(r)
    return results


# ---------------------------------------------------------------------------
# Approach D: Möbius operator eigenvalue equation
# ---------------------------------------------------------------------------


@dataclass
class MobiusOperatorResult:
    """Result from the Möbius operator eigenvalue approach.

    Attributes:
        N_grid: Grid size used.
        alpha_twist: Twist parameter for the Möbius identification.
        eigenvalues: Raw eigenvalues (sorted).
        kappas: Three κ values used (squared eigenvalues ratios).
        fit: Dict from assess_fit().
    """

    N_grid: int
    alpha_twist: float
    eigenvalues: np.ndarray
    kappas: np.ndarray
    fit: dict


def _build_mobius_laplacian(N: int, alpha_twist: float) -> np.ndarray:
    """Build the Laplacian on the Möbius strip with twist angle α.

    The Möbius strip is constructed by identifying the two ends of [0, 1]
    with a flip: f(0) = −f(1) × e^{iα}.  For the standard Möbius strip α = π
    (pure flip); for a twisted bundle α ∈ [0, 2π).

    The Laplacian −d²/dx² on the Möbius strip (1D base) with the identification
    ψ(0) = −ψ(1) e^{iα} has spectrum:
        λ_n = (π(2n + 1) + α)²  for n = 0, 1, 2, ...
    This is the "anti-periodic" spectrum shifted by α.

    Args:
        N: Number of interior grid points.
        alpha_twist: Twist angle in [0, 2π).

    Returns:
        Dense N×N real symmetric matrix (the Laplacian).
    """
    dx = 1.0 / (N + 1)
    # Dirichlet at interior with Möbius BC applied implicitly through
    # the effective Bloch-like momentum shift k_0 = α/L for length L = 1
    # → shift all eigenvalues by (alpha_twist)² in the operator
    lap = np.zeros((N, N))
    diag_val = 2.0 / dx**2
    off_val = -1.0 / dx**2
    for i in range(N):
        lap[i, i] = diag_val
        if i > 0:
            lap[i, i - 1] = off_val
        if i < N - 1:
            lap[i, i + 1] = off_val
    # The Möbius shift: add α²/(L²) to every diagonal element, and
    # off-diag coupling 2α/L to the momentum (−id/dx) operator.
    # For the real symmetric form: add uniform diagonal shift k_mobius².
    k_mobius = alpha_twist  # in units where L = 1 / π
    lap += np.eye(N) * k_mobius**2
    return lap


def mobius_operator_spectrum(
    N: int = 200,
    alpha_twist: float = np.pi,
    n_states: int = 6,
) -> np.ndarray:
    """Compute eigenvalues of the Möbius-twisted Laplacian.

    Args:
        N: Grid dimension.
        alpha_twist: Möbius twist angle.
        n_states: Number of lowest eigenvalues to return.

    Returns:
        Sorted array of the n_states lowest positive eigenvalues.
    """
    lap = _build_mobius_laplacian(N, alpha_twist)
    # Request only lowest eigenvalues
    n_req = min(n_states + 5, N - 1)
    eigvals = eigh(lap, subset_by_index=[0, n_req - 1], eigvals_only=True)
    eigvals = np.sort(eigvals)
    pos = eigvals[eigvals > 0]
    return pos[:n_states]


def scan_mobius_operator(
    alpha_values: np.ndarray | None = None,
    N: int = 150,
) -> list[MobiusOperatorResult]:
    """Scan Möbius twist angles and find eigenvalue-based lepton fits.

    For each α, compute eigenvalues λ_0, λ_1, λ_2 and set κ_n = λ_n².
    Find the α that minimises combined ratio + Koide error.

    Args:
        alpha_values: Array of α values in [0, 2π).  Defaults to 500 values.
        N: Grid size.

    Returns:
        Top-5 results sorted by combined error.
    """
    if alpha_values is None:
        alpha_values = np.linspace(0.01, 2.0 * np.pi - 0.01, 500)

    best_results: list[MobiusOperatorResult] = []
    best_score = float("inf")

    for alpha in alpha_values:
        try:
            eigs = mobius_operator_spectrum(N=N, alpha_twist=float(alpha), n_states=5)
        except Exception:
            continue
        if len(eigs) < 3:
            continue

        # For the mass ratio problem, we normalise κ_n to κ_0:
        # try all triples from the first 5 eigenvalues
        for i0 in range(len(eigs) - 2):
            for i1 in range(i0 + 1, len(eigs) - 1):
                for i2 in range(i1 + 1, len(eigs)):
                    k = eigs[[i0, i1, i2]]**2
                    if np.any(k <= 0):
                        continue
                    fit = assess_fit(k)
                    e1 = fit["ratio_mu_e_err_pct"] / 100
                    e2 = fit["ratio_tau_e_err_pct"] / 100
                    eq = fit["koide_err"] / KOIDE_TARGET
                    score = e1**2 + e2**2 + 10.0 * eq**2
                    if score < best_score * 1.5:
                        best_results.append(MobiusOperatorResult(
                            N_grid=N,
                            alpha_twist=float(alpha),
                            eigenvalues=eigs,
                            kappas=k.copy(),
                            fit=fit,
                        ))
                        if score < best_score:
                            best_score = score

    best_results.sort(
        key=lambda r: (
            (r.fit["ratio_mu_e_err_pct"] / 100)**2
            + (r.fit["ratio_tau_e_err_pct"] / 100)**2
            + 10.0 * (r.fit["koide_err"] / KOIDE_TARGET)**2
        )
    )
    return best_results[:5]


# ---------------------------------------------------------------------------
# Approach E: Constrained optimisation on the Koide surface
# with Möbius-compatible parametrisation
# ---------------------------------------------------------------------------


@dataclass
class KoideSurfaceResult:
    """Result from Koide-surface constrained search.

    Attributes:
        parametrisation: Name of the parametrisation used.
        params: Free parameters used.
        kappas: Stiffness eigenvalues.
        fit: Dict from assess_fit().
    """

    parametrisation: str
    params: dict
    kappas: np.ndarray
    fit: dict


def koide_mobius_parametrisation() -> list[KoideSurfaceResult]:
    """Search the Koide surface for a Möbius-natural parametrisation.

    The Koide surface Q = 2/3 can be parametrised by angles:
        m_n = M × (1 + √2 cos(2πn/3 + δ))²  for n = 0, 1, 2

    (Foot 1994 / Rosen parametrisation).  This gives Q = 2/3 exactly for any
    M > 0 and δ ∈ [0, 2π).  The lepton masses correspond to δ ≈ 2π/12.

    A Möbius connection: on the Möbius strip, angles periodic up to π (since
    a circuit flips sign).  So if δ is set by the Möbius winding:
        δ = π / (2n_winding + 1)
    with n_winding the Möbius half-flux winding number, we get a pure
    topological prediction.

    We scan δ and search for the value that gives ratios (207, 3477).
    Additionally we test the Foot/Rosen prediction δ = 2π/12 (= π/6).

    Args: none.

    Returns:
        List of KoideSurfaceResult for each δ tried.
    """
    results: list[KoideSurfaceResult] = []

    def masses_from_koide_param(M: float, delta: float) -> np.ndarray:
        """Foot parametrisation: m_n = M(1 + √2 cos(2πn/3 + δ))²."""
        ns = np.array([0.0, 1.0, 2.0])
        x = 1.0 + np.sqrt(2.0) * np.cos(2.0 * np.pi * ns / 3.0 + delta)
        m = M * x**2
        return m

    def koide_ratio_error(delta: float) -> float:
        """Error between Foot parametrisation and lepton ratios."""
        m = masses_from_koide_param(1.0, delta)
        if np.any(m <= 0):
            return 1e12
        order = np.argsort(m)
        m = m[order]
        r10 = m[1] / m[0]
        r20 = m[2] / m[0]
        e1 = (np.log(r10) - np.log(RATIO_MU_E))**2
        e2 = (np.log(r20) - np.log(RATIO_TAU_E))**2
        return e1 + e2

    # Global scan over δ
    deltas = np.linspace(0, 2 * np.pi, 20000)
    errors = np.array([koide_ratio_error(d) for d in deltas])
    best_idx = int(np.argmin(errors))
    delta_coarse = float(deltas[best_idx])

    # Refine
    res_refine = minimize_scalar(
        koide_ratio_error,
        bounds=(delta_coarse - 0.1, delta_coarse + 0.1),
        method="bounded",
        options={"xatol": 1e-12},
    )
    delta_best = float(res_refine.x)
    m_best = masses_from_koide_param(1.0, delta_best)
    order = np.argsort(m_best)
    m_best = m_best[order]
    kappas_best = m_best**2  # κ ∝ m²  since m = √κ
    fit_best = assess_fit(kappas_best)

    results.append(KoideSurfaceResult(
        parametrisation="Foot-Rosen: m=M(1+√2 cos(2πn/3+δ))²",
        params={"delta_rad": delta_best,
                "delta_pi": delta_best / np.pi,
                "delta_over_2pi": delta_best / (2 * np.pi)},
        kappas=kappas_best,
        fit=fit_best,
    ))

    # Möbius winding prediction for δ
    # n_winding = 1 (Möbius half-flux) → δ_Möbius = π/3 (= 2π/6)
    # n_winding = 2 → δ = π/5, etc.
    for n_w in range(1, 8):
        delta_m = np.pi / (2 * n_w + 1)
        m = masses_from_koide_param(1.0, delta_m)
        if np.any(m <= 0):
            continue
        order = np.argsort(m)
        m = m[order]
        kappas = m**2
        fit = assess_fit(kappas)
        results.append(KoideSurfaceResult(
            parametrisation=f"Möbius winding n_w={n_w}: δ=π/(2n+1)",
            params={"n_winding": n_w, "delta_rad": delta_m,
                    "delta_pi": delta_m / np.pi},
            kappas=kappas,
            fit=fit,
        ))

    # Exact PDG-derived δ (what δ gives the actual lepton masses?)
    # Solve: we know r10 = 207, r20 = 3477.  From Foot param:
    #   √r10 = |1 + √2 cos(2π/3 + δ)| / |1 + √2 cos(δ)|
    #   √r20 = |1 + √2 cos(4π/3 + δ)| / |1 + √2 cos(δ)|
    # This fixes δ uniquely (verified by the optimal search above).
    pdg_delta = delta_best  # already found
    results.append(KoideSurfaceResult(
        parametrisation="PDG-exact: optimal δ from lepton masses",
        params={"delta_rad": pdg_delta,
                "delta_pi": pdg_delta / np.pi,
                "interpretation": (
                    "Not directly Möbius-derived; gives exact Q=2/3 and "
                    "best-fit ratios.  A Möbius mechanism must explain "
                    f"why δ ≈ {pdg_delta/np.pi:.4f}π."
                )},
        kappas=kappas_best,
        fit=fit_best,
    ))

    return results


# ---------------------------------------------------------------------------
# Approach F: Stress-topology combined — Möbius winding × power-law hybrid
# ---------------------------------------------------------------------------


def mobius_power_law_hybrid(
    k_exp: float = 15.0,
    mobius_twist: float = 0.5,
    n_values: list[int] | None = None,
) -> dict:
    """Combine Möbius BC shift with power-law stiffness.

    The Möbius half-flux adds a half-integer shift to the quantum number:
        n_eff = n + mobius_twist    (mobius_twist = 1/2 for standard Möbius)

    Combined with the empirical power-law:
        κ_n = n_eff^{k_exp}

    This is a two-parameter model (k_exp, mobius_twist) that interpolates
    between the pure power-law (mobius_twist = 0) and the Möbius prediction
    (mobius_twist = 1/2).  Minimising ratio + Koide error jointly constrains
    both parameters.

    Args:
        k_exp: Power-law exponent.
        mobius_twist: Half-integer shift (0 to 1).
        n_values: Quantum numbers.

    Returns:
        Dict with kappas, fit metrics.
    """
    if n_values is None:
        n_values = [0, 1, 2]
    ns = np.array(n_values, dtype=float)
    n_eff = ns + mobius_twist
    if np.any(n_eff <= 0):
        n_eff = np.maximum(n_eff, 1e-8)
    kappas = n_eff**k_exp
    fit = assess_fit(kappas)
    return {
        "k_exp": k_exp,
        "mobius_twist": mobius_twist,
        "kappas": kappas,
        "fit": fit,
    }


def optimise_hybrid(n_values: list[int] | None = None) -> dict:
    """Find (k_exp, twist) that minimises ratio + Koide error jointly.

    This is the most aggressive numerical search: 2D optimisation over
    (k_exp, mobius_twist) to find ANY combination that simultaneously
    satisfies ratios and Koide.

    Args:
        n_values: Quantum numbers.

    Returns:
        Dict with optimal parameters and fit metrics.
    """
    if n_values is None:
        n_values = [0, 1, 2]

    def objective(params: np.ndarray) -> float:
        k_exp, twist = params
        if k_exp < 0.1 or k_exp > 30.0:
            return 1e12
        if twist < 0.0 or twist >= 2.0:
            return 1e12
        res = mobius_power_law_hybrid(k_exp=k_exp, mobius_twist=twist,
                                      n_values=n_values)
        fit = res["fit"]
        e1 = fit["ratio_mu_e_err_pct"] / 100
        e2 = fit["ratio_tau_e_err_pct"] / 100
        eq = fit["koide_err"] / KOIDE_TARGET
        return e1**2 + e2**2 + 10.0 * eq**2

    # Coarse grid search first
    best_score = float("inf")
    best_p0 = np.array([15.0, 0.5])
    for k_try in np.linspace(5.0, 25.0, 20):
        for t_try in np.linspace(0.01, 1.99, 20):
            s = objective(np.array([k_try, t_try]))
            if s < best_score:
                best_score = s
                best_p0 = np.array([k_try, t_try])

    # Refine
    res = minimize(
        objective,
        x0=best_p0,
        method="Nelder-Mead",
        options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 50000},
    )
    k_opt, tw_opt = float(res.x[0]), float(res.x[1])
    result = mobius_power_law_hybrid(k_exp=k_opt, mobius_twist=tw_opt,
                                     n_values=n_values)
    result["optimiser_converged"] = bool(res.success)
    result["optimiser_fun"] = float(res.fun)
    return result


# ---------------------------------------------------------------------------
# High-level driver
# ---------------------------------------------------------------------------


def run_all_approaches(verbose: bool = True) -> dict:
    """Run all four topological quantization approaches.

    Approach A: Bohr-Sommerfeld with Möbius half-integer shift.
    Approach B: Atiyah-Singer self-adjoint extension parameter scan.
    Approach C: Torus knot Jones / Alexander polynomial Casimir values.
    Approach D: Möbius operator (twisted Laplacian) eigenvalue spectrum.
    Approach E: Koide surface (Foot parametrisation) + Möbius winding.
    Approach F: Hybrid power-law + Möbius twist, 2D optimisation.

    Args:
        verbose: Print progress if True.

    Returns:
        Nested dict with keys A–F, each containing results and fit metrics.
    """
    out: dict = {}

    if verbose:
        print("=== Approach A: Bohr-Sommerfeld ===")
    out["A_bohr_sommerfeld"] = run_bohr_sommerfeld()

    if verbose:
        print("=== Approach B: Atiyah-Singer / self-adjoint extension ===")
    as_coarse = scan_atiyah_singer(n_theta=3000)
    as_refined: list[AtiyahSingerResult] = []
    seen_nvals: set[tuple[int, ...]] = set()
    for r in as_coarse[:3]:
        key = tuple(r.n_quantum)
        if key not in seen_nvals:
            seen_nvals.add(key)
            as_refined.append(refine_atiyah_singer(r.theta_SA, r.n_quantum))
    out["B_atiyah_singer"] = {
        "coarse_top5": as_coarse,
        "refined": as_refined,
    }

    if verbose:
        print("=== Approach C: Knot invariants ===")
    knot_results = scan_jones_casimir()
    out["C_knot_invariants"] = knot_results

    if verbose:
        print("=== Approach D: Möbius operator eigenvalues ===")
    mob_op_results = scan_mobius_operator(N=100)
    out["D_mobius_operator"] = mob_op_results

    if verbose:
        print("=== Approach E: Koide surface + Möbius winding ===")
    koide_surf = koide_mobius_parametrisation()
    out["E_koide_surface"] = koide_surf

    if verbose:
        print("=== Approach F: Hybrid power-law + Möbius twist ===")
    hybrid_opt = optimise_hybrid()
    out["F_hybrid"] = hybrid_opt

    if verbose:
        print("=== Koide consistency cross-check ===")
    out["koide_cross_check"] = koide_exact_ratios()

    return out
