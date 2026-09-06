"""Hypothesis tests for the current missing pieces in the substrate model.

This module is deliberately conservative.  It does not treat numerical
near-matches as derivations.  It asks a narrower question:

    If the open sectors need one extra mechanism, what numerical target must
    that mechanism hit, and do simple substrate-shaped candidates survive a
    first sanity check?

Covered gaps:
    - UV/Planck suppression chi_UV = l_P / xi_e.
    - Lepton Foot/Koide phase selection and branch conventions.
    - CKM/Cabibbo rational-angle selector degeneracy.
    - Pre-CMB matter/radiation transfer-window hiding.
    - Dark-matter dimer geometric self-interaction scale.
    - Matter/antimatter orientation bias needed without one-shot baryogenesis.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

from .substrate_planck_scale import (
    L_PLANCK_OBS,
    SIGMA_LAT,
    SIGMA_MAX,
    XI_E,
)


# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

M_E_MEV = 0.51099895
M_MU_MEV = 105.6583755
M_TAU_MEV = 1776.86
RATIO_MU_E = M_MU_MEV / M_E_MEV
RATIO_TAU_E = M_TAU_MEV / M_E_MEV
KOIDE_TARGET = 2.0 / 3.0

THETA_C_DEG = 13.04
THETA_C_RAD = math.radians(THETA_C_DEG)

GEV_TO_G = 1.78266192e-24
FM_TO_CM = 1.0e-13


def _pct_error(value: float, target: float) -> float:
    return (value - target) / target * 100.0


def _abs_pct_error(value: float, target: float) -> float:
    return abs(_pct_error(value, target))


# ---------------------------------------------------------------------------
# 1. UV / Planck suppression
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class UVSuppressionCandidate:
    """A dimensionless candidate for chi_UV = l_P / xi_e."""

    formula: str
    n_inst: int
    p_sigma_lat: int
    q_sigma_max: int
    r_twopi: int
    value: float
    target: float
    relative_error: float
    log_error: float
    complexity: int
    verdict: str


def uv_target_ratio() -> float:
    """Return the required dimensionless Planck/electron ratio."""

    return L_PLANCK_OBS / XI_E


def evaluate_uv_candidate(
    n_inst: int,
    p_sigma_lat: int = 0,
    q_sigma_max: int = 0,
    r_twopi: int = 0,
) -> UVSuppressionCandidate:
    """Evaluate exp(-n*pi) sigma_lat^p sigma_max^q (2*pi)^r.

    The exp(-n*pi) form is included because UV instanton/barrier actions are
    one of the few natural ways to generate exponentially small numbers.  A
    close match is still only numerology until the model derives the action
    integer and the prefactors from a substrate saddle.
    """

    target = uv_target_ratio()
    value = (
        math.exp(-n_inst * math.pi)
        * (SIGMA_LAT ** p_sigma_lat)
        * (SIGMA_MAX ** q_sigma_max)
        * ((2.0 * math.pi) ** r_twopi)
    )
    rel = value / target - 1.0
    log_error = math.log(value / target)
    pieces = [f"exp(-{n_inst}*pi)"]
    if p_sigma_lat:
        pieces.append(f"sigma_lat^{p_sigma_lat}")
    if q_sigma_max:
        pieces.append(f"sigma_max^{q_sigma_max}")
    if r_twopi:
        pieces.append(f"(2*pi)^{r_twopi}")
    complexity = (
        abs(n_inst)
        + abs(p_sigma_lat)
        + abs(q_sigma_max)
        + abs(r_twopi)
        + sum(1 for x in (p_sigma_lat, q_sigma_max, r_twopi) if x)
    )
    if abs(rel) < 0.1:
        verdict = "near numeric hit; needs a real UV action mechanism"
    elif abs(rel) < 1.0:
        verdict = "same order; still not a derivation"
    else:
        verdict = "miss"
    return UVSuppressionCandidate(
        formula=" * ".join(pieces),
        n_inst=n_inst,
        p_sigma_lat=p_sigma_lat,
        q_sigma_max=q_sigma_max,
        r_twopi=r_twopi,
        value=value,
        target=target,
        relative_error=rel,
        log_error=log_error,
        complexity=complexity,
        verdict=verdict,
    )


def search_uv_suppressions(
    max_inst: int = 24,
    max_sigma_power: int = 8,
    twopi_power_range: range = range(-2, 3),
    top_n: int = 12,
) -> list[UVSuppressionCandidate]:
    """Search a small family of exponential UV-suppression candidates."""

    candidates: list[UVSuppressionCandidate] = []
    for n_inst in range(1, max_inst + 1):
        for p_lat in range(0, max_sigma_power + 1):
            for q_max in range(0, max_sigma_power + 1):
                for r_twopi in twopi_power_range:
                    candidates.append(
                        evaluate_uv_candidate(
                            n_inst=n_inst,
                            p_sigma_lat=p_lat,
                            q_sigma_max=q_max,
                            r_twopi=r_twopi,
                        )
                    )
    # Penalize complexity so the report does not reward arbitrary products of
    # many near-1/2 factors over simpler action-scale clues.
    candidates.sort(key=lambda c: (abs(c.log_error) + 0.03 * c.complexity, abs(c.log_error)))
    return candidates[:top_n]


# ---------------------------------------------------------------------------
# 2. Lepton Foot/Koide phase and branch tests
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FootPhasePoint:
    """Diagnostic for one Foot phase."""

    label: str
    delta_pi: float
    ratio_mu_e: float
    ratio_tau_e: float
    signed_koide_q: float
    positive_koide_q: float
    error_mu_pct: float
    error_tau_pct: float
    root_signs: tuple[int, int, int]
    verdict: str


@dataclass(frozen=True)
class RationalFootCandidate:
    """A rational-pi Foot phase candidate."""

    p: int
    q: int
    delta_pi: float
    error_mu_pct: float
    error_tau_pct: float
    combined_error_pct: float


def _foot_roots(delta: float) -> tuple[float, float, float]:
    return tuple(
        1.0 + math.sqrt(2.0) * math.cos(2.0 * math.pi * n / 3.0 + delta)
        for n in range(3)
    )


def _positive_koide_q_from_masses(masses: Sequence[float]) -> float:
    roots = [math.sqrt(m) for m in masses]
    return sum(masses) / (sum(roots) ** 2)


def _signed_koide_q_from_roots(roots: Sequence[float]) -> float:
    masses = [x * x for x in roots]
    return sum(masses) / (sum(roots) ** 2)


def foot_phase_point(label: str, delta_pi: float) -> FootPhasePoint:
    """Evaluate a Foot phase using both signed and positive-root Koide tests."""

    roots = _foot_roots(delta_pi * math.pi)
    masses = sorted(x * x for x in roots)
    if masses[0] <= 1.0e-30:
        ratio_mu = float("inf")
        ratio_tau = float("inf")
    else:
        ratio_mu = masses[1] / masses[0]
        ratio_tau = masses[2] / masses[0]
    signed_q = _signed_koide_q_from_roots(roots)
    positive_q = _positive_koide_q_from_masses(masses)
    signs = tuple(1 if x >= 0.0 else -1 for x in roots)
    err_mu = _pct_error(ratio_mu, RATIO_MU_E) if math.isfinite(ratio_mu) else float("inf")
    err_tau = _pct_error(ratio_tau, RATIO_TAU_E) if math.isfinite(ratio_tau) else float("inf")
    if any(s < 0 for s in signs):
        verdict = "signed Koide branch only; positive masses do not keep Q=2/3"
    elif abs(err_mu) < 1.0 and abs(err_tau) < 1.0:
        verdict = "empirical physical branch"
    else:
        verdict = "physical branch but wrong ratios"
    return FootPhasePoint(
        label=label,
        delta_pi=delta_pi % 2.0,
        ratio_mu_e=ratio_mu,
        ratio_tau_e=ratio_tau,
        signed_koide_q=signed_q,
        positive_koide_q=positive_q,
        error_mu_pct=err_mu,
        error_tau_pct=err_tau,
        root_signs=signs,
        verdict=verdict,
    )


def _foot_loss(delta: float) -> float:
    roots = _foot_roots(delta)
    masses = sorted(x * x for x in roots)
    if masses[0] <= 0.0:
        return float("inf")
    return (
        math.log((masses[1] / masses[0]) / RATIO_MU_E) ** 2
        + math.log((masses[2] / masses[0]) / RATIO_TAU_E) ** 2
    )


def _golden_minimize(func, lo: float, hi: float, tol: float = 1.0e-13) -> float:
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - gr * (hi - lo)
    d = lo + gr * (hi - lo)
    fc = func(c)
    fd = func(d)
    while abs(hi - lo) > tol:
        if fc < fd:
            hi = d
            d = c
            fd = fc
            c = hi - gr * (hi - lo)
            fc = func(c)
        else:
            lo = c
            c = d
            fc = fd
            d = lo + gr * (hi - lo)
            fd = func(d)
    return (lo + hi) / 2.0


def find_empirical_foot_phase(samples: int = 20000) -> float:
    """Find the best Foot phase in units of pi without using scipy."""

    step = 2.0 * math.pi / samples
    best_i = min(range(samples), key=lambda i: _foot_loss(i * step))
    center = best_i * step
    lo = center - 3.0 * step
    hi = center + 3.0 * step
    delta = _golden_minimize(_foot_loss, lo, hi) % (2.0 * math.pi)
    # The sorted-mass Foot map has a Z3 permutation degeneracy.  Return the
    # smallest member of that orbit so reports are stable across grid starts.
    delta_pi = delta / math.pi
    return min(((delta_pi + shift) % 2.0) for shift in (0.0, 2.0 / 3.0, 4.0 / 3.0))


def foot_z3_orbit(delta_pi: float) -> tuple[float, float, float]:
    """Return the Z3-shift orbit delta, delta+2/3, delta+4/3 in units of pi."""

    vals = sorted(((delta_pi + shift) % 2.0) for shift in (0.0, 2.0 / 3.0, 4.0 / 3.0))
    return tuple(vals)


def scan_rational_foot_phases(max_q: int = 30, top_n: int = 8) -> list[RationalFootCandidate]:
    """Scan reduced rational phases p/q*pi and rank by lepton-ratio error."""

    seen: set[tuple[int, int]] = set()
    out: list[RationalFootCandidate] = []
    for q in range(1, max_q + 1):
        for p in range(0, 2 * q + 1):
            g = math.gcd(p, q)
            pr = p // g
            qr = q // g
            if (pr, qr) in seen:
                continue
            seen.add((pr, qr))
            point = foot_phase_point(f"{pr}/{qr}", pr / qr)
            combined = math.hypot(point.error_mu_pct, point.error_tau_pct)
            out.append(
                RationalFootCandidate(
                    p=pr,
                    q=qr,
                    delta_pi=pr / qr,
                    error_mu_pct=point.error_mu_pct,
                    error_tau_pct=point.error_tau_pct,
                    combined_error_pct=combined,
                )
            )
    out.sort(key=lambda r: abs(r.combined_error_pct))
    return out[:top_n]


def lepton_phase_diagnostics() -> dict[str, object]:
    """Return the core lepton-phase diagnostics."""

    empirical = find_empirical_foot_phase()
    conjugate = (2.0 - empirical) % 2.0
    return {
        "empirical_delta_pi": empirical,
        "conjugate_delta_pi": conjugate,
        "z3_orbit": foot_z3_orbit(empirical),
        "empirical_point": foot_phase_point("empirical branch", empirical),
        "conjugate_point": foot_phase_point("conjugate branch", conjugate),
        "pi_over_6_point": foot_phase_point("pi/6 topology", 1.0 / 6.0),
        "seven_pi_over_five_point": foot_phase_point("7pi/5 rational", 7.0 / 5.0),
        "best_rationals": scan_rational_foot_phases(),
    }


# ---------------------------------------------------------------------------
# 3. CKM / Cabibbo selector tests
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CKMRationalCandidate:
    """A rational multiple of pi for the Cabibbo angle."""

    numerator: int
    denominator: int
    theta_deg: float
    error_deg_pct: float
    formula: str


def scan_ckm_rationals(
    max_numerator: int = 10,
    max_denominator: int = 120,
    max_abs_error_pct: float = 1.0,
) -> list[CKMRationalCandidate]:
    """Scan theta_C = n*pi/d and show selector degeneracy."""

    out: list[CKMRationalCandidate] = []
    seen: set[tuple[int, int]] = set()
    for den in range(2, max_denominator + 1):
        for num in range(1, min(max_numerator, den - 1) + 1):
            g = math.gcd(num, den)
            nr = num // g
            dr = den // g
            if (nr, dr) in seen:
                continue
            seen.add((nr, dr))
            deg = 180.0 * nr / dr
            err = _pct_error(deg, THETA_C_DEG)
            if abs(err) <= max_abs_error_pct:
                out.append(
                    CKMRationalCandidate(
                        numerator=nr,
                        denominator=dr,
                        theta_deg=deg,
                        error_deg_pct=err,
                        formula=f"{nr}*pi/{dr}",
                    )
                )
    out.sort(key=lambda c: abs(c.error_deg_pct))
    return out


def ckm_selector_examples_for_55() -> tuple[str, ...]:
    """Examples showing that denominator 55 is not uniquely selected."""

    return (
        "55 = 5 * 11",
        "55 = 2 * 3^3 + 1",
        "55 = 7 * 8 - 1",
        "55 = (4 + 1) * 11",
    )


# ---------------------------------------------------------------------------
# 4. Cosmology transfer-window tests
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TransferWindowCandidate:
    """A toy ratio W_gamma(k) / W_m(k)."""

    family: str
    k_cut_mpc: float
    steepness: int
    f_at_acoustic: float
    f_at_galaxy: float
    required_f_vis: float
    passes_visibility: bool
    keeps_acoustic_visible: bool
    verdict: str


def required_visibility(delta_m: float = 0.0236, cmb_delta_t: float = 1.0e-5) -> float:
    """Required radiation leakage factor f_vis <= deltaT/T / delta_m."""

    return cmb_delta_t / delta_m


def _transfer_ratio(family: str, k: float, k_cut: float, steepness: int) -> float:
    x = k / k_cut
    if family == "exp":
        return math.exp(-(x ** steepness))
    if family == "lorentz":
        return 1.0 / (1.0 + x ** steepness)
    raise ValueError(f"unknown transfer family: {family}")


def scan_transfer_windows(
    k_acoustic: float = 0.05,
    k_galaxy: float = 1.0,
    k_cuts: Sequence[float] = (0.05, 0.075, 0.10, 0.15, 0.20, 0.30),
    steepnesses: Sequence[int] = (1, 2, 4, 6, 8),
) -> list[TransferWindowCandidate]:
    """Find toy transfer windows that hide galaxy-scale seeds from the CMB."""

    req = required_visibility()
    out: list[TransferWindowCandidate] = []
    for family in ("exp", "lorentz"):
        for k_cut in k_cuts:
            for steepness in steepnesses:
                f_acoustic = _transfer_ratio(family, k_acoustic, k_cut, steepness)
                f_galaxy = _transfer_ratio(family, k_galaxy, k_cut, steepness)
                passes = f_galaxy <= req
                keeps = f_acoustic >= 0.5
                if passes and keeps:
                    verdict = "mathematically viable; physical opacity mechanism still missing"
                elif passes:
                    verdict = "hides seeds but also suppresses acoustic-scale radiation"
                else:
                    verdict = "too visible in the CMB"
                out.append(
                    TransferWindowCandidate(
                        family=family,
                        k_cut_mpc=k_cut,
                        steepness=steepness,
                        f_at_acoustic=f_acoustic,
                        f_at_galaxy=f_galaxy,
                        required_f_vis=req,
                        passes_visibility=passes,
                        keeps_acoustic_visible=keeps,
                        verdict=verdict,
                    )
                )
    out.sort(
        key=lambda c: (
            not (c.passes_visibility and c.keeps_acoustic_visible),
            abs(math.log10(max(c.f_at_galaxy, 1.0e-300) / c.required_f_vis)),
            -c.f_at_acoustic,
        )
    )
    return out


# ---------------------------------------------------------------------------
# 5. Dark matter dimer cross-section scale
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DarkMatterCrossSection:
    """Geometric self-interaction estimate for a neutral dimer."""

    mass_gev: float
    radius_fm: float
    sigma_cm2: float
    sigma_over_m_cm2_per_g: float
    verdict: str


def dark_matter_geometric_cross_sections(
    mass_gev: float = 48.6,
    radii_fm: Iterable[float] = (0.2, 1.0, 5.0, 10.0, 50.0, 100.0),
) -> list[DarkMatterCrossSection]:
    """Estimate sigma_self/m from a geometric radius."""

    mass_g = mass_gev * GEV_TO_G
    out: list[DarkMatterCrossSection] = []
    for radius_fm in radii_fm:
        radius_cm = radius_fm * FM_TO_CM
        sigma_cm2 = math.pi * radius_cm * radius_cm
        sigma_over_m = sigma_cm2 / mass_g
        if sigma_over_m < 1.0e-2:
            verdict = "effectively collisionless at halo scale"
        elif sigma_over_m < 1.0:
            verdict = "self-interacting only if composite radius is enlarged"
        else:
            verdict = "large self-interaction; likely constrained"
        out.append(
            DarkMatterCrossSection(
                mass_gev=mass_gev,
                radius_fm=radius_fm,
                sigma_cm2=sigma_cm2,
                sigma_over_m_cm2_per_g=sigma_over_m,
                verdict=verdict,
            )
        )
    return out


def radius_for_self_interaction(
    target_sigma_over_m: float,
    mass_gev: float = 48.6,
) -> float:
    """Radius in fm required for a target sigma/m in cm^2/g."""

    mass_g = mass_gev * GEV_TO_G
    radius_cm = math.sqrt(target_sigma_over_m * mass_g / math.pi)
    return radius_cm / FM_TO_CM


# ---------------------------------------------------------------------------
# 6. Matter orientation / no one-shot baryogenesis
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OrientationBiasRequirement:
    """Required bias or relaxation for suppressing antimatter orientation."""

    target_anti_fraction: float
    delta_e_over_t_eff: float
    tau_over_epoch_max: float
    verdict: str


def orientation_bias_requirements(
    target_fractions: Sequence[float] = (1.0e-9, 1.0e-18, 1.0e-30),
) -> list[OrientationBiasRequirement]:
    """Bias needed if anti-orientation fraction is Boltzmann suppressed.

    If f_anti ~ exp(-DeltaE / T_eff), then DeltaE/T_eff >= ln(1/f).
    If f_anti relaxes as exp(-t/tau), then tau/epoch <= 1/ln(1/f).
    """

    out: list[OrientationBiasRequirement] = []
    for frac in target_fractions:
        req = math.log(1.0 / frac)
        out.append(
            OrientationBiasRequirement(
                target_anti_fraction=frac,
                delta_e_over_t_eff=req,
                tau_over_epoch_max=1.0 / req,
                verdict=(
                    "needs a real orientation-selecting phase transition; "
                    "not solved by saying baryogenesis did not happen"
                ),
            )
        )
    return out


# ---------------------------------------------------------------------------
# High-level report bundle
# ---------------------------------------------------------------------------


def run_missing_piece_hypothesis_tests() -> dict[str, object]:
    """Run all missing-piece hypothesis tests."""

    return {
        "uv_target_ratio": uv_target_ratio(),
        "uv_candidates": search_uv_suppressions(),
        "lepton_phase": lepton_phase_diagnostics(),
        "ckm_rationals": scan_ckm_rationals(),
        "ckm_55_examples": ckm_selector_examples_for_55(),
        "transfer_windows": scan_transfer_windows(),
        "dm_cross_sections": dark_matter_geometric_cross_sections(),
        "dm_radius_for_0p1": radius_for_self_interaction(0.1),
        "dm_radius_for_1": radius_for_self_interaction(1.0),
        "orientation_bias": orientation_bias_requirements(),
    }
