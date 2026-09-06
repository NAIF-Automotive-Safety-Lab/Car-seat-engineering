"""Candidate self-dual exchange mechanism for the 45 degree cone.

The lattice-invariant audit showed that the cone quartic needs a symmetry that
removes the lower-order mismatch m = p - q, where

    p = |grad_parallel phi|^2
    q = |grad_perp phi|^2

This module tests a minimal local mechanism: a paired exchange cell with two
dual branches.  One branch is frustrated by longitudinal loading, the other by
transverse loading.  If the two branches have equal weight, their linear
anisotropies cancel and the remaining orientation selector is positive m^2.

This is not yet a microscopic proof.  It moves the gap to a sharper condition:
derive exact equal branch weight / local detailed balance for the dual pair.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


def _validate_branch_parameters(branch_weight: float, branch_stiffness: float) -> None:
    if not 0.0 <= branch_weight <= 1.0:
        raise ValueError("branch_weight must be in [0, 1]")
    if branch_stiffness <= 0.0:
        raise ValueError("branch_stiffness must be positive")


def mismatch_from_angle(theta_rad: float) -> float:
    """Return m = p - q at fixed p + q = 1."""

    return math.cos(2.0 * theta_rad)


def dual_branch_energy(
    theta_rad: float,
    *,
    branch_weight: float = 0.5,
    branch_stiffness: float = 1.0,
) -> float:
    """Return the paired-branch exchange energy at fixed |grad phi|.

    branch_weight is the population/weight of the branch frustrated by +m.  The
    dual branch has weight 1 - branch_weight and is frustrated by -m.
    """

    _validate_branch_parameters(branch_weight, branch_stiffness)
    m = mismatch_from_angle(theta_rad)
    return 0.25 * branch_stiffness * (
        branch_weight * (1.0 + m) ** 2
        + (1.0 - branch_weight) * (1.0 - m) ** 2
    )


def effective_linear_bias(
    *,
    branch_weight: float = 0.5,
    branch_stiffness: float = 1.0,
) -> float:
    """Return the m coefficient after expanding the paired-branch energy."""

    _validate_branch_parameters(branch_weight, branch_stiffness)
    return 0.5 * branch_stiffness * (2.0 * branch_weight - 1.0)


def effective_beta(
    *,
    branch_weight: float = 0.5,
    branch_stiffness: float = 1.0,
) -> float:
    """Return beta in E_orient = const + linear*m + beta/4*m^2."""

    _validate_branch_parameters(branch_weight, branch_stiffness)
    return branch_stiffness


def analytic_minimum_angle(
    *,
    branch_weight: float = 0.5,
    branch_stiffness: float = 1.0,
) -> float:
    """Return the analytic minimum angle in radians."""

    _validate_branch_parameters(branch_weight, branch_stiffness)
    delta = 2.0 * branch_weight - 1.0
    minimum_m = max(-1.0, min(1.0, -delta))
    return 0.5 * math.acos(minimum_m)


@dataclass(frozen=True)
class SelfDualExchangeAssessment:
    """Assessment of the paired-branch self-dual mechanism."""

    balanced_linear_bias: float
    balanced_beta: float
    balanced_minimum_angle_deg: float
    balanced_energy_at_minimum: float
    imbalanced_weight: float
    imbalanced_linear_bias: float
    imbalanced_minimum_angle_deg: float
    imbalanced_shift_deg: float
    single_branch_minimum_angle_deg: float
    beta_positive_from_branch_stability: bool
    dual_pair_cancels_quadratic_bias: bool
    conditional_cone_closure: bool
    fully_derived: bool
    verdict: str


def assess_self_dual_exchange_mechanism(
    *,
    imbalanced_weight: float = 0.55,
    branch_stiffness: float = 1.0,
) -> SelfDualExchangeAssessment:
    """Assess whether a dual exchange pair can produce the cone quartic."""

    balanced_angle = analytic_minimum_angle(
        branch_weight=0.5,
        branch_stiffness=branch_stiffness,
    )
    imbalanced_angle = analytic_minimum_angle(
        branch_weight=imbalanced_weight,
        branch_stiffness=branch_stiffness,
    )
    single_branch_angle = analytic_minimum_angle(
        branch_weight=1.0,
        branch_stiffness=branch_stiffness,
    )
    balanced_bias = effective_linear_bias(
        branch_weight=0.5,
        branch_stiffness=branch_stiffness,
    )
    imbalanced_bias = effective_linear_bias(
        branch_weight=imbalanced_weight,
        branch_stiffness=branch_stiffness,
    )
    beta = effective_beta(
        branch_weight=0.5,
        branch_stiffness=branch_stiffness,
    )
    balanced_energy = dual_branch_energy(
        balanced_angle,
        branch_weight=0.5,
        branch_stiffness=branch_stiffness,
    )
    shift = math.degrees(imbalanced_angle - balanced_angle)

    beta_positive = beta > 0.0
    bias_cancelled = abs(balanced_bias) < 1.0e-12
    conditional_closure = (
        bias_cancelled
        and beta_positive
        and abs(math.degrees(balanced_angle) - 45.0) < 1.0e-12
        and abs(shift) > 1.0
    )

    if conditional_closure:
        verdict = (
            "paired self-dual branches produce zero linear bias and positive "
            "beta when branch weights are exactly equal; deriving exact equal "
            "weight remains the open microphysics"
        )
    else:
        verdict = "paired self-dual exchange does not close the cone condition"

    return SelfDualExchangeAssessment(
        balanced_linear_bias=balanced_bias,
        balanced_beta=beta,
        balanced_minimum_angle_deg=math.degrees(balanced_angle),
        balanced_energy_at_minimum=balanced_energy,
        imbalanced_weight=imbalanced_weight,
        imbalanced_linear_bias=imbalanced_bias,
        imbalanced_minimum_angle_deg=math.degrees(imbalanced_angle),
        imbalanced_shift_deg=shift,
        single_branch_minimum_angle_deg=math.degrees(single_branch_angle),
        beta_positive_from_branch_stability=beta_positive,
        dual_pair_cancels_quadratic_bias=bias_cancelled,
        conditional_cone_closure=conditional_closure,
        fully_derived=False,
        verdict=verdict,
    )
