"""Anchor-mediated origin for the cone branch-exchange spring.

The diamond-cell selection audit left two local-cell assumptions:

    1. two saturated anchors,
    2. a direct L-T branch-exchange spring.

This module tests whether the second assumption can be derived from the first.
If the anchors are not infinitely rigid but are saturated clamps with finite
compliance, eliminating their local variables by a Schur complement induces an
effective negative off-diagonal term in the L/T branch Hessian.  That term is
exactly the branch-exchange spring used in the diamond cell.

This narrows the remaining cone microgeometry problem: derive the shared
finite-compliance saturated anchors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cone_detailed_balance import (
    cone_angle_from_detailed_balance,
    involution_commutator_norm,
    stationary_branch_weight,
)
from .cone_diamond_cell_geometry import SpringEdge, graph_laplacian
from .cone_self_dual_exchange import effective_linear_bias
from .cone_swap_generator_origin import (
    branch_energy_splitting,
    cell_automorphism_residual,
    generator_from_cell_hessian,
    rates_from_cell_hessian,
)


NODES: tuple[str, ...] = ("L", "T", "A", "B")
BRANCH_INDICES: tuple[int, int] = (0, 1)
ANCHOR_INDICES: tuple[int, int] = (2, 3)


def shared_anchor_edges(*, branch_anchor_stiffness: float = 1.0) -> tuple[SpringEdge, ...]:
    """Return the no-direct-LT graph where both branches share both anchors."""

    if branch_anchor_stiffness <= 0.0:
        raise ValueError("branch_anchor_stiffness must be positive")
    k = branch_anchor_stiffness
    return (
        SpringEdge("L", "A", k),
        SpringEdge("L", "B", k),
        SpringEdge("T", "A", k),
        SpringEdge("T", "B", k),
    )


def pinned_anchor_hessian(
    *,
    branch_anchor_stiffness: float = 1.0,
    anchor_pin_stiffness: float = 1.0,
) -> np.ndarray:
    """Return the full L,T,A,B Hessian with finite anchor pinning to ground."""

    if anchor_pin_stiffness < 0.0:
        raise ValueError("anchor_pin_stiffness must be non-negative")
    laplacian = graph_laplacian(
        shared_anchor_edges(branch_anchor_stiffness=branch_anchor_stiffness),
        nodes=NODES,
    )
    for index in ANCHOR_INDICES:
        laplacian[index, index] += anchor_pin_stiffness
    return laplacian


def schur_branch_hessian(full_hessian: np.ndarray) -> np.ndarray:
    """Eliminate anchors and return the effective L/T branch Hessian."""

    h = np.asarray(full_hessian, dtype=float)
    if h.shape != (4, 4):
        raise ValueError("full_hessian must have shape (4, 4)")
    h_bb = h[np.ix_(BRANCH_INDICES, BRANCH_INDICES)]
    h_ba = h[np.ix_(BRANCH_INDICES, ANCHOR_INDICES)]
    h_aa = h[np.ix_(ANCHOR_INDICES, ANCHOR_INDICES)]
    return h_bb - h_ba @ np.linalg.inv(h_aa) @ h_ba.T


def induced_exchange_strength(branch_hessian: np.ndarray) -> float:
    """Return g_eff for H_eff = [[k_eff, -g_eff], [-g_eff, k_eff]]."""

    h = np.asarray(branch_hessian, dtype=float)
    if h.shape != (2, 2):
        raise ValueError("branch_hessian must have shape (2, 2)")
    return float(-h[0, 1])


def analytic_induced_exchange(
    *,
    branch_anchor_stiffness: float = 1.0,
    anchor_pin_stiffness: float = 1.0,
) -> float:
    """Return the expected induced exchange for two shared anchors."""

    if branch_anchor_stiffness <= 0.0:
        raise ValueError("branch_anchor_stiffness must be positive")
    if anchor_pin_stiffness < 0.0:
        raise ValueError("anchor_pin_stiffness must be non-negative")
    k = branch_anchor_stiffness
    s = anchor_pin_stiffness
    return 2.0 * k * k / (2.0 * k + s)


def fixed_anchor_limit_hessian(*, branch_anchor_stiffness: float = 1.0) -> np.ndarray:
    """Return branch Hessian for the same square cell with anchors fixed."""

    if branch_anchor_stiffness <= 0.0:
        raise ValueError("branch_anchor_stiffness must be positive")
    k = branch_anchor_stiffness
    return np.array([[2.0 * k, 0.0], [0.0, 2.0 * k]], dtype=float)


@dataclass(frozen=True)
class AnchorInducedExchangeAssessment:
    """Assessment of anchor-mediated branch exchange."""

    branch_anchor_stiffness: float
    anchor_pin_stiffness: float
    induced_exchange: float
    analytic_exchange: float
    fixed_anchor_exchange: float
    branch_automorphism_residual: float
    stationary_weight: float
    linear_bias: float
    minimum_angle_deg: float
    generator_commutator_norm: float
    soft_anchor_limit_exchange: float
    rigid_anchor_limit_exchange: float
    exchange_induced_by_finite_anchors: bool
    fully_derived: bool
    verdict: str


def assess_anchor_induced_exchange(
    *,
    branch_anchor_stiffness: float = 1.0,
    anchor_pin_stiffness: float = 1.0,
) -> AnchorInducedExchangeAssessment:
    """Assess whether shared finite anchors induce the branch-exchange spring."""

    full_hessian = pinned_anchor_hessian(
        branch_anchor_stiffness=branch_anchor_stiffness,
        anchor_pin_stiffness=anchor_pin_stiffness,
    )
    branch_hessian = schur_branch_hessian(full_hessian)
    induced = induced_exchange_strength(branch_hessian)
    analytic = analytic_induced_exchange(
        branch_anchor_stiffness=branch_anchor_stiffness,
        anchor_pin_stiffness=anchor_pin_stiffness,
    )
    fixed_exchange = induced_exchange_strength(
        fixed_anchor_limit_hessian(branch_anchor_stiffness=branch_anchor_stiffness)
    )
    rate_plus, rate_minus = rates_from_cell_hessian(branch_hessian)
    weight = stationary_branch_weight(
        rate_plus_to_minus=rate_plus,
        rate_minus_to_plus=rate_minus,
    )
    angle = cone_angle_from_detailed_balance(
        rate_plus_to_minus=rate_plus,
        rate_minus_to_plus=rate_minus,
    )
    generator = generator_from_cell_hessian(branch_hessian)
    soft_exchange = analytic_induced_exchange(
        branch_anchor_stiffness=branch_anchor_stiffness,
        anchor_pin_stiffness=1.0e-9,
    )
    rigid_exchange = analytic_induced_exchange(
        branch_anchor_stiffness=branch_anchor_stiffness,
        anchor_pin_stiffness=1.0e9,
    )

    induced_ok = (
        induced > 0.0
        and abs(induced - analytic) < 1.0e-12
        and fixed_exchange == 0.0
        and cell_automorphism_residual(branch_hessian) < 1.0e-12
        and abs(weight - 0.5) < 1.0e-12
        and abs(angle - 45.0) < 1.0e-12
        and soft_exchange > induced
        and rigid_exchange < 1.0e-8
    )

    if induced_ok:
        verdict = (
            "finite-compliance shared saturated anchors induce the effective "
            "L-T exchange spring; the direct exchange no longer needs to be a "
            "separate primitive"
        )
    else:
        verdict = "shared anchors do not induce the needed branch exchange"

    return AnchorInducedExchangeAssessment(
        branch_anchor_stiffness=branch_anchor_stiffness,
        anchor_pin_stiffness=anchor_pin_stiffness,
        induced_exchange=induced,
        analytic_exchange=analytic,
        fixed_anchor_exchange=fixed_exchange,
        branch_automorphism_residual=cell_automorphism_residual(branch_hessian),
        stationary_weight=weight,
        linear_bias=effective_linear_bias(branch_weight=weight),
        minimum_angle_deg=angle,
        generator_commutator_norm=involution_commutator_norm(generator),
        soft_anchor_limit_exchange=soft_exchange,
        rigid_anchor_limit_exchange=rigid_exchange,
        exchange_induced_by_finite_anchors=induced_ok,
        fully_derived=False,
        verdict=verdict,
    )
