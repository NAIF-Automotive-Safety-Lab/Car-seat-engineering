"""Local detailed-balance candidate for exact dual-branch weights.

The paired self-dual cone mechanism needs exact 50/50 weights for two dual
branches.  This module tests the smallest dynamical route to that condition:
a two-state local exchange generator whose states are related by a dual swap
involution.

If the local generator commutes with the swap, the transition rates are equal,
the stationary branch weights are exactly 1/2, and the cone quartic has no
linear anisotropic bias.  If a rate or energy splitting is introduced, detailed
balance shifts the stationary weight and the cone angle moves away from 45 deg.

This closes the weight algebra conditionally.  The remaining microphysics is to
derive the swap-degenerate local exchange generator from the substrate cell.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .cone_self_dual_exchange import analytic_minimum_angle, effective_linear_bias


def exchange_generator(
    *,
    rate_plus_to_minus: float,
    rate_minus_to_plus: float,
) -> np.ndarray:
    """Return a two-state continuous-time generator for column probabilities."""

    if rate_plus_to_minus <= 0.0 or rate_minus_to_plus <= 0.0:
        raise ValueError("rates must be positive")
    return np.array(
        [
            [-rate_plus_to_minus, rate_minus_to_plus],
            [rate_plus_to_minus, -rate_minus_to_plus],
        ],
        dtype=float,
    )


def swap_involution() -> np.ndarray:
    """Return the branch-dual swap matrix J with J^2 = I."""

    return np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float)


def involution_commutator_norm(generator: np.ndarray) -> float:
    """Return ||GJ - JG||_F for the branch-dual swap."""

    g = np.asarray(generator, dtype=float)
    if g.shape != (2, 2):
        raise ValueError("generator must have shape (2, 2)")
    j = swap_involution()
    return float(np.linalg.norm(g @ j - j @ g))


def stationary_branch_weight(
    *,
    rate_plus_to_minus: float,
    rate_minus_to_plus: float,
) -> float:
    """Return stationary probability of the plus branch."""

    if rate_plus_to_minus <= 0.0 or rate_minus_to_plus <= 0.0:
        raise ValueError("rates must be positive")
    return rate_minus_to_plus / (rate_plus_to_minus + rate_minus_to_plus)


def boltzmann_branch_weight(*, energy_splitting_over_temp: float) -> float:
    """Return plus-branch weight for delta=(E_plus-E_minus)/T."""

    delta = energy_splitting_over_temp
    if delta >= 0.0:
        exp_neg = math.exp(-delta)
        return exp_neg / (1.0 + exp_neg)
    exp_pos = math.exp(delta)
    return 1.0 / (1.0 + exp_pos)


def cone_angle_from_detailed_balance(
    *,
    rate_plus_to_minus: float,
    rate_minus_to_plus: float,
) -> float:
    """Return the cone minimum angle in degrees implied by stationary weights."""

    weight = stationary_branch_weight(
        rate_plus_to_minus=rate_plus_to_minus,
        rate_minus_to_plus=rate_minus_to_plus,
    )
    return math.degrees(analytic_minimum_angle(branch_weight=weight))


@dataclass(frozen=True)
class DetailedBalanceAssessment:
    """Assessment of the local detailed-balance cone closure."""

    symmetric_commutator_norm: float
    symmetric_stationary_weight: float
    symmetric_linear_bias: float
    symmetric_minimum_angle_deg: float
    split_energy_over_temp: float
    split_stationary_weight: float
    split_linear_bias: float
    split_minimum_angle_deg: float
    split_angle_shift_deg: float
    rate_imbalance: float
    imbalanced_stationary_weight: float
    imbalanced_commutator_norm: float
    imbalanced_minimum_angle_deg: float
    detailed_balance_closes_equal_weight: bool
    fully_derived: bool
    verdict: str


def assess_detailed_balance_closure(
    *,
    split_energy_over_temp: float = 0.2,
    rate_imbalance: float = 0.1,
) -> DetailedBalanceAssessment:
    """Assess exact 50/50 branch weight from a swap-symmetric generator."""

    if rate_imbalance <= -1.0:
        raise ValueError("rate_imbalance must be > -1")

    symmetric_g = exchange_generator(
        rate_plus_to_minus=1.0,
        rate_minus_to_plus=1.0,
    )
    symmetric_weight = stationary_branch_weight(
        rate_plus_to_minus=1.0,
        rate_minus_to_plus=1.0,
    )
    symmetric_angle = cone_angle_from_detailed_balance(
        rate_plus_to_minus=1.0,
        rate_minus_to_plus=1.0,
    )
    symmetric_bias = effective_linear_bias(branch_weight=symmetric_weight)

    split_weight = boltzmann_branch_weight(
        energy_splitting_over_temp=split_energy_over_temp,
    )
    split_angle = math.degrees(analytic_minimum_angle(branch_weight=split_weight))
    split_bias = effective_linear_bias(branch_weight=split_weight)

    imbalanced_g = exchange_generator(
        rate_plus_to_minus=1.0,
        rate_minus_to_plus=1.0 + rate_imbalance,
    )
    imbalanced_weight = stationary_branch_weight(
        rate_plus_to_minus=1.0,
        rate_minus_to_plus=1.0 + rate_imbalance,
    )
    imbalanced_angle = cone_angle_from_detailed_balance(
        rate_plus_to_minus=1.0,
        rate_minus_to_plus=1.0 + rate_imbalance,
    )

    closes_equal_weight = (
        involution_commutator_norm(symmetric_g) < 1.0e-12
        and abs(symmetric_weight - 0.5) < 1.0e-12
        and abs(symmetric_bias) < 1.0e-12
        and abs(symmetric_angle - 45.0) < 1.0e-12
        and abs(split_angle - 45.0) > 1.0
        and abs(imbalanced_angle - 45.0) > 1.0
    )

    if closes_equal_weight:
        verdict = (
            "swap-symmetric local detailed balance gives exact 50/50 branch "
            "weights; rate or energy splitting reintroduces cone drift, so the "
            "remaining task is deriving the swap-degenerate generator"
        )
    else:
        verdict = "local detailed balance does not close the branch-weight gate"

    return DetailedBalanceAssessment(
        symmetric_commutator_norm=involution_commutator_norm(symmetric_g),
        symmetric_stationary_weight=symmetric_weight,
        symmetric_linear_bias=symmetric_bias,
        symmetric_minimum_angle_deg=symmetric_angle,
        split_energy_over_temp=split_energy_over_temp,
        split_stationary_weight=split_weight,
        split_linear_bias=split_bias,
        split_minimum_angle_deg=split_angle,
        split_angle_shift_deg=split_angle - symmetric_angle,
        rate_imbalance=rate_imbalance,
        imbalanced_stationary_weight=imbalanced_weight,
        imbalanced_commutator_norm=involution_commutator_norm(imbalanced_g),
        imbalanced_minimum_angle_deg=imbalanced_angle,
        detailed_balance_closes_equal_weight=closes_equal_weight,
        fully_derived=False,
        verdict=verdict,
    )
