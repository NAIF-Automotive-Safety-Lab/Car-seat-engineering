"""Neutral-coupling suppression audit for the dark-stress speed.

The tensor-mode audit gives the denominator in

    v_dark = alpha*c/sqrt(5).

This module asks what effective stiffness ratio is required for the remaining
alpha factor.  For a wave mode,

    v/c = sqrt(K_eff/K) / sqrt(N_modes).

Therefore the observed closure requires K_eff/K = alpha^2 when N_modes = 5.
That is exactly the scaling expected if neutral stress is a second-order
charge-symmetric response.  This does not derive alpha itself; it makes the
remaining derivation target precise.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .dark_stress_scale_closure import ALPHA_EM, dark_stress_speed_closure
from .neutral_stress_tensor_modes import assess_neutral_stress_modes
from .substrate_polarization_dm import C_SI


def speed_from_stiffness_ratio(
    stiffness_ratio: float,
    mode_count: int = 5,
) -> float:
    """Return mode speed in km/s for K_eff/K and mode count."""

    if stiffness_ratio < 0.0:
        raise ValueError("stiffness_ratio must be non-negative")
    if mode_count <= 0:
        raise ValueError("mode_count must be positive")
    return C_SI * math.sqrt(stiffness_ratio / mode_count) / 1000.0


def required_stiffness_ratio_for_speed(
    speed_km_s: float,
    mode_count: int = 5,
) -> float:
    """Return K_eff/K required for a given mode speed."""

    if speed_km_s < 0.0:
        raise ValueError("speed_km_s must be non-negative")
    if mode_count <= 0:
        raise ValueError("mode_count must be positive")
    return mode_count * (speed_km_s * 1000.0 / C_SI) ** 2


@dataclass(frozen=True)
class StiffnessCandidate:
    """One candidate neutral-stress stiffness scaling."""

    name: str
    stiffness_ratio: float
    speed_km_s: float
    speed_error_pct: float
    verdict: str


@dataclass(frozen=True)
class NeutralCouplingSuppressionAssessment:
    """Assessment of the alpha speed-suppression mechanism."""

    mode_count: int
    target_speed_km_s: float
    required_stiffness_ratio: float
    alpha_squared: float
    stiffness_error_pct: float
    candidates: tuple[StiffnessCandidate, ...]
    verdict: str


def _candidate(name: str, stiffness_ratio: float, target_speed: float) -> StiffnessCandidate:
    speed = speed_from_stiffness_ratio(stiffness_ratio)
    err = (speed / target_speed - 1.0) * 100.0
    if abs(err) < 1.0:
        verdict = "matches dark-stress speed"
    elif speed > target_speed:
        verdict = "too fast"
    else:
        verdict = "too slow"
    return StiffnessCandidate(
        name=name,
        stiffness_ratio=stiffness_ratio,
        speed_km_s=speed,
        speed_error_pct=err,
        verdict=verdict,
    )


def assess_neutral_coupling_suppression() -> NeutralCouplingSuppressionAssessment:
    """Assess whether alpha speed suppression equals alpha^2 stiffness."""

    modes = assess_neutral_stress_modes()
    target = dark_stress_speed_closure().v_dark_km_s
    required = required_stiffness_ratio_for_speed(target, modes.projector_rank)
    alpha_squared = ALPHA_EM**2
    stiffness_error = (required / alpha_squared - 1.0) * 100.0
    candidates = (
        _candidate("unsuppressed neutral stiffness: K_eff/K = 1", 1.0, target),
        _candidate("first-order neutral stiffness: K_eff/K = alpha", ALPHA_EM, target),
        _candidate(
            "second-order neutral stiffness: K_eff/K = alpha^2",
            alpha_squared,
            target,
        ),
    )

    if abs(stiffness_error) < 1.0e-9 and candidates[-1].verdict.startswith("matches"):
        verdict = (
            "dark-stress speed is equivalent to a second-order neutral stiffness "
            "K_eff/K=alpha^2 over five symmetric-traceless modes; derive this "
            "quadratic neutral response from the substrate action"
        )
    else:
        verdict = "alpha suppression is not explained by second-order neutral stiffness"

    return NeutralCouplingSuppressionAssessment(
        mode_count=modes.projector_rank,
        target_speed_km_s=target,
        required_stiffness_ratio=required,
        alpha_squared=alpha_squared,
        stiffness_error_pct=stiffness_error,
        candidates=candidates,
        verdict=verdict,
    )
