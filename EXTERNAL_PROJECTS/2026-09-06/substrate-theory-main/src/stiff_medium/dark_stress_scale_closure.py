"""Scale-closure candidates for dark-stress coherence length and speed.

The memory-clock trial narrowed the dimensional clock to

    tau_clock = ell_pol / v_dark.

This module tests whether both factors can be expressed from existing
dimensionless structure and the cosmological Hubble scale.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .dark_stress_memory_clock import relaxation_count
from .dark_stress_hybrid import hybrid_cluster_check
from .substrate_polarization_dm import C_SI, KPC_M, MPC_M, YEAR_S


ALPHA_EM = 1.0 / 137.035999084
H0_KM_S_MPC = 67.4
MYR_S = 1.0e6 * YEAR_S


def h0_si(h0_km_s_mpc: float = H0_KM_S_MPC) -> float:
    """Return H0 in s^-1."""

    return h0_km_s_mpc * 1000.0 / MPC_M


def hubble_length_kpc(h0_km_s_mpc: float = H0_KM_S_MPC) -> float:
    """Return c/H0 in kpc."""

    return C_SI / h0_si(h0_km_s_mpc) / KPC_M


@dataclass(frozen=True)
class PolarizationCoherenceClosure:
    """Candidate for ell_pol."""

    formula: str
    ell_pol_kpc: float
    mechanism: str
    verdict: str


def polarization_coherence_length_closure() -> PolarizationCoherenceClosure:
    """ell_pol = alpha^3 (c/H0) / sqrt(3)."""

    ell = ALPHA_EM**3 * hubble_length_kpc() / math.sqrt(3.0)
    return PolarizationCoherenceClosure(
        formula="alpha^3 * (c/H0) / sqrt(3)",
        ell_pol_kpc=ell,
        mechanism=(
            "three alpha-suppressed charge-symmetric filters applied to the "
            "Hubble coherence length, projected across three spatial axes"
        ),
        verdict=(
            "near 1 kpc coherence scale; derive alpha^3 filtering and 3D projection"
            if 0.8 <= ell <= 1.2
            else "wrong coherence scale"
        ),
    )


@dataclass(frozen=True)
class DarkSpeedClosure:
    """Candidate for v_dark."""

    formula: str
    v_dark_km_s: float
    mechanism: str
    verdict: str


def dark_stress_speed_closure() -> DarkSpeedClosure:
    """v_dark = alpha c / sqrt(5)."""

    v_km_s = ALPHA_EM * C_SI / math.sqrt(5.0) / 1000.0
    return DarkSpeedClosure(
        formula="alpha * c / sqrt(5)",
        v_dark_km_s=v_km_s,
        mechanism=(
            "alpha-suppressed substrate propagation distributed over the five "
            "symmetric-traceless shear components of neutral stress"
        ),
        verdict=(
            "cluster-scale dark-stress speed; derive the five-mode projection"
            if 800.0 <= v_km_s <= 1200.0
            else "wrong speed scale"
        ),
    )


@dataclass(frozen=True)
class ScaleClosureAssessment:
    """Combined coherence-length/speed closure."""

    coherence: PolarizationCoherenceClosure
    speed: DarkSpeedClosure
    tau_clock_myr: float
    tau_pol_myr: float
    required_tau_myr: float
    tau_error_pct: float
    verdict: str


def assess_dark_stress_scale_closure() -> ScaleClosureAssessment:
    """Assess ell_pol and v_dark closure together."""

    coherence = polarization_coherence_length_closure()
    speed = dark_stress_speed_closure()
    tau_clock_myr = (
        coherence.ell_pol_kpc * KPC_M / (speed.v_dark_km_s * 1000.0) / MYR_S
    )
    tau_pol = relaxation_count() * tau_clock_myr
    required = hybrid_cluster_check().required_memory_myr
    err = (tau_pol / required - 1.0) * 100.0
    verdict = (
        "scale closure is strong; derive alpha^3/sqrt(3) and alpha^2 neutral stiffness from dynamics"
        if abs(err) < 1.0
        else "scale closure misses memory target"
    )
    return ScaleClosureAssessment(
        coherence=coherence,
        speed=speed,
        tau_clock_myr=tau_clock_myr,
        tau_pol_myr=tau_pol,
        required_tau_myr=required,
        tau_error_pct=err,
        verdict=verdict,
    )
