"""Diamond-cell geometry candidate for the cone branch-swap automorphism.

The previous step showed that a two-branch elastic Hessian with

    J^T H J = H

is sufficient to force the swap-degenerate exchange generator.  This module
constructs a concrete local spring graph that has that automorphism:

    saturated anchor A
       /          \\
    branch L -- branch T
       \\          /
    saturated anchor B

The branch reservoirs L and T are coupled to the same saturated anchors with
equal spring weights.  Swapping L and T leaves the graph Laplacian invariant,
and the fixed-anchor branch Hessian inherits J^T H J = H.

This is still a candidate geometry, not a proof that this is the substrate's
actual microscopic cell.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cone_detailed_balance import (
    cone_angle_from_detailed_balance,
    involution_commutator_norm,
    stationary_branch_weight,
)
from .cone_self_dual_exchange import effective_linear_bias
from .cone_swap_generator_origin import (
    branch_energy_splitting,
    cell_automorphism_residual,
    generator_from_cell_hessian,
    rates_from_cell_hessian,
)


NODES: tuple[str, ...] = ("L", "T", "A", "B")
BRANCH_NODES: tuple[str, ...] = ("L", "T")


@dataclass(frozen=True)
class SpringEdge:
    """One scalar spring edge in the local cell graph."""

    node_a: str
    node_b: str
    stiffness: float


def diamond_cell_edges(
    *,
    anchor_stiffness: float = 0.5,
    exchange_stiffness: float = 0.25,
    anchor_split: float = 0.0,
) -> tuple[SpringEdge, ...]:
    """Return the diamond cell's spring graph.

    anchor_split is a failure-control perturbation.  It makes the L-anchor
    springs stronger and T-anchor springs weaker by the same amount.
    """

    if anchor_stiffness <= 0.0:
        raise ValueError("anchor_stiffness must be positive")
    if exchange_stiffness < 0.0:
        raise ValueError("exchange_stiffness must be non-negative")
    if abs(anchor_split) >= anchor_stiffness:
        raise ValueError("anchor_split must be smaller than anchor_stiffness")
    l_anchor = anchor_stiffness + anchor_split
    t_anchor = anchor_stiffness - anchor_split
    return (
        SpringEdge("L", "A", l_anchor),
        SpringEdge("L", "B", l_anchor),
        SpringEdge("T", "A", t_anchor),
        SpringEdge("T", "B", t_anchor),
        SpringEdge("L", "T", exchange_stiffness),
    )


def graph_laplacian(edges: tuple[SpringEdge, ...], nodes: tuple[str, ...] = NODES) -> np.ndarray:
    """Return the scalar spring graph Laplacian."""

    index = {node: i for i, node in enumerate(nodes)}
    laplacian = np.zeros((len(nodes), len(nodes)), dtype=float)
    for edge in edges:
        if edge.stiffness < 0.0:
            raise ValueError("edge stiffness must be non-negative")
        if edge.node_a not in index or edge.node_b not in index:
            raise ValueError("edge uses a node outside the cell")
        i = index[edge.node_a]
        j = index[edge.node_b]
        laplacian[i, i] += edge.stiffness
        laplacian[j, j] += edge.stiffness
        laplacian[i, j] -= edge.stiffness
        laplacian[j, i] -= edge.stiffness
    return laplacian


def branch_swap_permutation(nodes: tuple[str, ...] = NODES) -> np.ndarray:
    """Return the full-cell permutation that swaps L and T."""

    swapped = {"L": "T", "T": "L", "A": "A", "B": "B"}
    index = {node: i for i, node in enumerate(nodes)}
    permutation = np.zeros((len(nodes), len(nodes)), dtype=float)
    for node in nodes:
        permutation[index[swapped[node]], index[node]] = 1.0
    return permutation


def graph_automorphism_residual(laplacian: np.ndarray) -> float:
    """Return ||P^T L P - L|| for the full branch-swap graph automorphism."""

    lap = np.asarray(laplacian, dtype=float)
    if lap.shape != (len(NODES), len(NODES)):
        raise ValueError("laplacian has wrong shape")
    p = branch_swap_permutation()
    return float(np.linalg.norm(p.T @ lap @ p - lap))


def fixed_anchor_branch_hessian(laplacian: np.ndarray) -> np.ndarray:
    """Return the L/T branch Hessian with saturated anchors fixed."""

    lap = np.asarray(laplacian, dtype=float)
    if lap.shape != (len(NODES), len(NODES)):
        raise ValueError("laplacian has wrong shape")
    branch_indices = [NODES.index(node) for node in BRANCH_NODES]
    return lap[np.ix_(branch_indices, branch_indices)]


@dataclass(frozen=True)
class DiamondCellGeometryAssessment:
    """Assessment of the diamond-cell branch-swap candidate."""

    symmetric_graph_residual: float
    symmetric_branch_residual: float
    symmetric_stationary_weight: float
    symmetric_linear_bias: float
    symmetric_minimum_angle_deg: float
    broken_anchor_split: float
    broken_graph_residual: float
    broken_branch_residual: float
    broken_branch_energy_over_temp: float
    broken_commutator_norm: float
    broken_stationary_weight: float
    broken_linear_bias: float
    broken_minimum_angle_deg: float
    broken_angle_shift_deg: float
    diamond_cell_forces_automorphism: bool
    fully_derived: bool
    verdict: str


def assess_diamond_cell_geometry(
    *,
    anchor_split: float = 0.05,
    temperature: float = 1.0,
) -> DiamondCellGeometryAssessment:
    """Assess whether the diamond spring cell gives the needed automorphism."""

    symmetric_lap = graph_laplacian(diamond_cell_edges(anchor_split=0.0))
    symmetric_h = fixed_anchor_branch_hessian(symmetric_lap)
    symmetric_rate_plus, symmetric_rate_minus = rates_from_cell_hessian(
        symmetric_h,
        temperature=temperature,
    )
    symmetric_weight = stationary_branch_weight(
        rate_plus_to_minus=symmetric_rate_plus,
        rate_minus_to_plus=symmetric_rate_minus,
    )
    symmetric_angle = cone_angle_from_detailed_balance(
        rate_plus_to_minus=symmetric_rate_plus,
        rate_minus_to_plus=symmetric_rate_minus,
    )

    broken_lap = graph_laplacian(diamond_cell_edges(anchor_split=anchor_split))
    broken_h = fixed_anchor_branch_hessian(broken_lap)
    broken_g = generator_from_cell_hessian(broken_h, temperature=temperature)
    broken_rate_plus, broken_rate_minus = rates_from_cell_hessian(
        broken_h,
        temperature=temperature,
    )
    broken_weight = stationary_branch_weight(
        rate_plus_to_minus=broken_rate_plus,
        rate_minus_to_plus=broken_rate_minus,
    )
    broken_angle = cone_angle_from_detailed_balance(
        rate_plus_to_minus=broken_rate_plus,
        rate_minus_to_plus=broken_rate_minus,
    )
    broken_delta = branch_energy_splitting(broken_h) / temperature

    closes = (
        graph_automorphism_residual(symmetric_lap) < 1.0e-12
        and cell_automorphism_residual(symmetric_h) < 1.0e-12
        and abs(symmetric_weight - 0.5) < 1.0e-12
        and abs(symmetric_angle - 45.0) < 1.0e-12
        and graph_automorphism_residual(broken_lap) > 1.0e-6
        and cell_automorphism_residual(broken_h) > 1.0e-6
        and involution_commutator_norm(broken_g) > 1.0e-6
        and abs(broken_angle - symmetric_angle) > 1.0
    )

    if closes:
        verdict = (
            "a symmetric saturated diamond cell supplies the branch-swap "
            "elastic automorphism; deriving this cell from deeper substrate "
            "microgeometry remains open"
        )
    else:
        verdict = "diamond spring cell does not supply the cone automorphism"

    return DiamondCellGeometryAssessment(
        symmetric_graph_residual=graph_automorphism_residual(symmetric_lap),
        symmetric_branch_residual=cell_automorphism_residual(symmetric_h),
        symmetric_stationary_weight=symmetric_weight,
        symmetric_linear_bias=effective_linear_bias(branch_weight=symmetric_weight),
        symmetric_minimum_angle_deg=symmetric_angle,
        broken_anchor_split=anchor_split,
        broken_graph_residual=graph_automorphism_residual(broken_lap),
        broken_branch_residual=cell_automorphism_residual(broken_h),
        broken_branch_energy_over_temp=broken_delta,
        broken_commutator_norm=involution_commutator_norm(broken_g),
        broken_stationary_weight=broken_weight,
        broken_linear_bias=effective_linear_bias(branch_weight=broken_weight),
        broken_minimum_angle_deg=broken_angle,
        broken_angle_shift_deg=broken_angle - symmetric_angle,
        diamond_cell_forces_automorphism=closes,
        fully_derived=False,
        verdict=verdict,
    )
