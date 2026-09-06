"""Concrete mechanism trials for the current weak sectors.

These are stronger than free scans but still not final derivations.  Each
trial proposes a physical object, computes the implied number, and marks the
remaining missing step.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .missing_piece_hypotheses import (
    THETA_C_DEG,
    dark_matter_geometric_cross_sections,
    find_empirical_foot_phase,
    foot_phase_point,
    required_visibility,
    uv_target_ratio,
)
from .substrate_planck_scale import SIGMA_LAT, SIGMA_MAX


ALPHA_EM = 1.0 / 137.035999084
XI_QCD_FM = 0.2


@dataclass(frozen=True)
class UVMechanismTrial:
    """A phase-slip action candidate for the UV ratio."""

    name: str
    action: float
    chi: float
    target_action: float
    action_error_pct: float
    chi_error_pct: float
    mechanism: str
    verdict: str


def _uv_trial(name: str, action: float, mechanism: str) -> UVMechanismTrial:
    target_chi = uv_target_ratio()
    target_action = -math.log(target_chi)
    chi = math.exp(-action)
    action_error = (action - target_action) / target_action * 100.0
    chi_error = (chi / target_chi - 1.0) * 100.0
    if abs(action_error) < 0.2:
        verdict = "action-level hit; derive the saddle or reject"
    elif abs(action_error) < 2.0:
        verdict = "same action scale; incomplete"
    else:
        verdict = "miss"
    return UVMechanismTrial(
        name=name,
        action=action,
        chi=chi,
        target_action=target_action,
        action_error_pct=action_error,
        chi_error_pct=chi_error,
        mechanism=mechanism,
        verdict=verdict,
    )


def uv_phase_slip_trials() -> list[UVMechanismTrial]:
    """Try closed phase-slip actions built from spin cycles and determinants."""

    return [
        _uv_trial(
            "4 spin cycles + sigma_lat determinant",
            16.0 * math.pi - 2.0 * math.log(SIGMA_LAT),
            "closed 4 x 4pi spin-cycle slip, with two lattice fluctuation modes",
        ),
        _uv_trial(
            "4 spin cycles + half-flux closure + zero mode",
            17.0 * math.pi - math.log(2.0 * math.pi),
            "4 x 4pi slip plus one pi Mobius closure, divided by a rotational zero mode",
        ),
        _uv_trial(
            "4 spin cycles + sigma_max determinant",
            16.0 * math.pi - 2.0 * math.log(SIGMA_MAX),
            "closed 4 x 4pi spin-cycle slip with ideal saturation determinant",
        ),
        _uv_trial(
            "bare 4 spin cycles",
            16.0 * math.pi,
            "closed 4 x 4pi spin-cycle slip with no determinant correction",
        ),
    ]


@dataclass(frozen=True)
class LeptonMechanismTrial:
    """A positive-root Foot branch mechanism candidate."""

    name: str
    delta_pi: float
    empirical_delta_pi: float
    phase_error_pct: float
    mu_error_pct: float
    tau_error_pct: float
    mechanism: str
    verdict: str


def lepton_boundary_loop_trials() -> list[LeptonMechanismTrial]:
    """Try boundary plus loop repulsion formulas for the Foot phase.

    The physical positive-root branch begins at delta/pi = 7/12, where the
    formerly negative square-root amplitude crosses zero.  A loop-scale
    repulsion from that zero root is modeled by 1/(8*pi^2), with an optional
    spin-cycle correction 1 - 1/(16*pi^2).
    """

    empirical = find_empirical_foot_phase()
    candidates = [
        (
            "positive-root boundary + one-loop repulsion",
            7.0 / 12.0 + 1.0 / (8.0 * math.pi**2),
            "zero-root boundary at 7/12 plus one-loop eigenvalue repulsion",
        ),
        (
            "boundary + one-loop with spin-cycle correction",
            7.0 / 12.0
            + (1.0 / (8.0 * math.pi**2)) * (1.0 - 1.0 / (16.0 * math.pi**2)),
            "zero-root boundary plus one-loop repulsion corrected by a 4pi spin cycle",
        ),
        (
            "clean topology pi/6",
            1.0 / 6.0,
            "signed-root Mobius/Z3 topology without positive-root branch selection",
        ),
    ]
    out: list[LeptonMechanismTrial] = []
    for name, delta_pi, mechanism in candidates:
        point = foot_phase_point(name, delta_pi)
        phase_error = (delta_pi - empirical) / empirical * 100.0
        combined = math.hypot(point.error_mu_pct, point.error_tau_pct)
        if combined < 0.5:
            verdict = "promising; now derive the loop term from O_vertex"
        elif combined < 2.0:
            verdict = "near miss; mechanism shape is plausible"
        else:
            verdict = "fails lepton ratios"
        out.append(
            LeptonMechanismTrial(
                name=name,
                delta_pi=delta_pi,
                empirical_delta_pi=empirical,
                phase_error_pct=phase_error,
                mu_error_pct=point.error_mu_pct,
                tau_error_pct=point.error_tau_pct,
                mechanism=mechanism,
                verdict=verdict,
            )
        )
    return out


@dataclass(frozen=True)
class CKMMechanismTrial:
    """A Cabibbo-angle substrate-overlap candidate."""

    name: str
    sin_theta: float
    theta_deg: float
    sin_error_pct: float
    theta_error_pct: float
    mechanism: str
    verdict: str


def ckm_overlap_trials() -> list[CKMMechanismTrial]:
    """Try simple normalized overlap amplitudes for Cabibbo mixing."""

    candidates = [
        (
            "two-axis half-flux overlap",
            1.0 / (math.pi * math.sqrt(2.0)),
            "overlap of two orthogonal half-flux angular modes normalized by pi*sqrt(2)",
        ),
        (
            "overlap with spin-cycle correction",
            (1.0 / (math.pi * math.sqrt(2.0))) * (1.0 + 1.0 / (16.0 * math.pi**2)),
            "same overlap with a small 4pi spin-cycle correction",
        ),
    ]
    out: list[CKMMechanismTrial] = []
    for name, sin_theta, mechanism in candidates:
        theta_deg = math.degrees(math.asin(sin_theta))
        sin_error = (sin_theta / 0.2255 - 1.0) * 100.0
        theta_error = (theta_deg / THETA_C_DEG - 1.0) * 100.0
        if abs(sin_error) < 0.5:
            verdict = "promising overlap scale; derive H_mix"
        else:
            verdict = "near scale only"
        out.append(
            CKMMechanismTrial(
                name=name,
                sin_theta=sin_theta,
                theta_deg=theta_deg,
                sin_error_pct=sin_error,
                theta_error_pct=theta_error,
                mechanism=mechanism,
                verdict=verdict,
            )
        )
    return out


@dataclass(frozen=True)
class CosmologyMechanismTrial:
    """A concrete opacity/transfer window candidate."""

    name: str
    k_cut_mpc: float
    steepness: int
    f_acoustic: float
    f_galaxy: float
    required_f_vis: float
    mechanism: str
    verdict: str


def cosmology_biharmonic_opacity_trial() -> CosmologyMechanismTrial:
    """Fourth-order elastic opacity window for saturated-phase radiation."""

    k_cut = 0.1
    steepness = 4
    k_acoustic = 0.05
    k_galaxy = 1.0
    f_acoustic = 1.0 / (1.0 + (k_acoustic / k_cut) ** steepness)
    f_galaxy = 1.0 / (1.0 + (k_galaxy / k_cut) ** steepness)
    req = required_visibility()
    if f_galaxy <= req and f_acoustic >= 0.5:
        verdict = "passes transfer numerics; derive k_cut from phase-front physics"
    else:
        verdict = "fails transfer requirements"
    return CosmologyMechanismTrial(
        name="biharmonic saturated-opacity window",
        k_cut_mpc=k_cut,
        steepness=steepness,
        f_acoustic=f_acoustic,
        f_galaxy=f_galaxy,
        required_f_vis=req,
        mechanism="fourth-order elastic opacity W_gamma/W_m = 1/(1+(k/k_c)^4)",
        verdict=verdict,
    )


@dataclass(frozen=True)
class DarkMatterMechanismTrial:
    """A neutral dimer radius/cross-section candidate."""

    name: str
    radius_fm: float
    sigma_over_m: float
    mechanism: str
    verdict: str


def dark_matter_polarization_halo_trial() -> DarkMatterMechanismTrial:
    """Elastic polarization halo radius R = xi_QCD / alpha."""

    radius = XI_QCD_FM / ALPHA_EM
    row = dark_matter_geometric_cross_sections(radii_fm=(radius,))[0]
    if 0.1 <= row.sigma_over_m_cm2_per_g <= 1.0:
        verdict = "self-interacting range; derive halo profile dynamically"
    elif row.sigma_over_m_cm2_per_g < 0.1:
        verdict = "mostly collisionless"
    else:
        verdict = "likely over-interacting"
    return DarkMatterMechanismTrial(
        name="neutral polarization halo",
        radius_fm=radius,
        sigma_over_m=row.sigma_over_m_cm2_per_g,
        mechanism="compact neutral core with elastic halo R = xi_QCD/alpha",
        verdict=verdict,
    )


@dataclass(frozen=True)
class OrientationMechanismTrial:
    """A de-saturation orientation bias candidate."""

    name: str
    action: float
    anti_fraction: float
    mechanism: str
    verdict: str


def orientation_vortex_bias_trials() -> list[OrientationMechanismTrial]:
    """Try closed U(1) orientation vortex actions."""

    candidates = [
        (
            "closed orientation vortex",
            4.0 * math.pi**2,
            "one closed U(1) orientation vortex action",
        ),
        (
            "orientation vortex + zero-mode determinant",
            4.0 * math.pi**2 + math.log(2.0 * math.pi),
            "closed U(1) orientation vortex plus rotational determinant",
        ),
        (
            "orientation vortex + saturation offset",
            4.0 * math.pi**2 + 2.0,
            "closed orientation vortex with saturation-barrier offset",
        ),
    ]
    out: list[OrientationMechanismTrial] = []
    for name, action, mechanism in candidates:
        anti_fraction = math.exp(-action)
        if anti_fraction <= 2.0e-18:
            verdict = "near required suppression; derive orientation action"
        elif anti_fraction <= 1.0e-16:
            verdict = "same suppression scale; not enough alone"
        else:
            verdict = "too weak"
        out.append(
            OrientationMechanismTrial(
                name=name,
                action=action,
                anti_fraction=anti_fraction,
                mechanism=mechanism,
                verdict=verdict,
            )
        )
    return out


def run_mechanism_trials() -> dict[str, object]:
    """Run all concrete mechanism trials."""

    return {
        "uv": uv_phase_slip_trials(),
        "leptons": lepton_boundary_loop_trials(),
        "ckm": ckm_overlap_trials(),
        "cosmology": cosmology_biharmonic_opacity_trial(),
        "dark_matter": dark_matter_polarization_halo_trial(),
        "orientation": orientation_vortex_bias_trials(),
    }
