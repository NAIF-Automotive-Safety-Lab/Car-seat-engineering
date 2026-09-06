"""Minimal-cell selection audit for the saturated diamond cone geometry.

The diamond-cell candidate supplies the branch-swap automorphism needed for the
45 degree cone.  This module asks a sharper question: among small spring graphs
with two branch reservoirs (L,T) and two saturated anchors (A,B), is the diamond
the minimal graph selected by explicit local constraints?

Result: yes, conditionally.  If the local cell must

    1. preserve the L <-> T branch swap,
    2. load both saturated anchors symmetrically, and
    3. include a direct branch-exchange spring,

then the unique minimal graph is the saturated diamond:

    L-A, T-A, L-B, T-B, L-T

Dropping condition 2 selects a one-anchor wedge; dropping condition 3 selects
the square without direct branch exchange.  This keeps the remaining assumption
honest: derive two saturated anchors and direct branch exchange from the
substrate microgeometry.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from .cone_diamond_cell_geometry import (
    NODES,
    SpringEdge,
    fixed_anchor_branch_hessian,
    graph_automorphism_residual,
    graph_laplacian,
)
from .cone_swap_generator_origin import cell_automorphism_residual


Edge = tuple[str, str]


ALL_EDGES: tuple[Edge, ...] = (
    ("L", "T"),
    ("L", "A"),
    ("L", "B"),
    ("T", "A"),
    ("T", "B"),
    ("A", "B"),
)

DIAMOND_EDGES: frozenset[Edge] = frozenset(
    {
        ("L", "T"),
        ("L", "A"),
        ("L", "B"),
        ("T", "A"),
        ("T", "B"),
    }
)


def edge_label(edge: Edge) -> str:
    """Return a stable human-readable edge label."""

    return f"{edge[0]}-{edge[1]}"


def graph_edges_to_springs(edges: frozenset[Edge]) -> tuple[SpringEdge, ...]:
    """Return unit-stiffness spring edges for an abstract graph."""

    return tuple(SpringEdge(a, b, 1.0) for a, b in sorted(edges))


def powerset_edges() -> tuple[frozenset[Edge], ...]:
    """Return every undirected graph on the four local cell nodes."""

    graphs: list[frozenset[Edge]] = []
    for size in range(len(ALL_EDGES) + 1):
        for subset in combinations(ALL_EDGES, size):
            graphs.append(frozenset(subset))
    return tuple(graphs)


def is_branch_swap_invariant(edges: frozenset[Edge]) -> bool:
    """Return True if the graph Laplacian is invariant under L <-> T."""

    laplacian = graph_laplacian(graph_edges_to_springs(edges))
    return graph_automorphism_residual(laplacian) < 1.0e-12


def branch_hessian_has_swap_automorphism(edges: frozenset[Edge]) -> bool:
    """Return True if the fixed-anchor branch Hessian has J^T H J = H."""

    laplacian = graph_laplacian(graph_edges_to_springs(edges))
    hessian = fixed_anchor_branch_hessian(laplacian)
    return cell_automorphism_residual(hessian) < 1.0e-12


def uses_both_saturated_anchors(edges: frozenset[Edge]) -> bool:
    """Return True if both A and B directly load both branch reservoirs."""

    return {
        ("L", "A"),
        ("T", "A"),
        ("L", "B"),
        ("T", "B"),
    }.issubset(edges)


def has_direct_branch_exchange(edges: frozenset[Edge]) -> bool:
    """Return True if the L-T branch-exchange spring is present."""

    return ("L", "T") in edges


def branch_hessian_is_positive(edges: frozenset[Edge]) -> bool:
    """Return True if the branch Hessian is positive definite."""

    laplacian = graph_laplacian(graph_edges_to_springs(edges))
    hessian = fixed_anchor_branch_hessian(laplacian)
    eigenvalues = np.linalg.eigvalsh(hessian)
    return bool(np.min(eigenvalues) > 1.0e-12)


def admissible_graphs(
    *,
    require_two_anchors: bool = True,
    require_direct_exchange: bool = True,
) -> tuple[frozenset[Edge], ...]:
    """Return graphs satisfying the stated local-cell constraints."""

    accepted: list[frozenset[Edge]] = []
    for edges in powerset_edges():
        if not edges:
            continue
        if not is_branch_swap_invariant(edges):
            continue
        if not branch_hessian_has_swap_automorphism(edges):
            continue
        if not branch_hessian_is_positive(edges):
            continue
        if require_two_anchors and not uses_both_saturated_anchors(edges):
            continue
        if require_direct_exchange and not has_direct_branch_exchange(edges):
            continue
        accepted.append(edges)
    return tuple(accepted)


def minimal_graphs(graphs: tuple[frozenset[Edge], ...]) -> tuple[frozenset[Edge], ...]:
    """Return all graphs with the smallest edge count."""

    if not graphs:
        return ()
    min_size = min(len(edges) for edges in graphs)
    return tuple(edges for edges in graphs if len(edges) == min_size)


def graph_signature(edges: frozenset[Edge]) -> tuple[str, ...]:
    """Return stable edge labels for reports and tests."""

    return tuple(sorted(edge_label(edge) for edge in edges))


@dataclass(frozen=True)
class DiamondCellSelectionAssessment:
    """Assessment of the minimal-cell selection audit."""

    total_graphs_scanned: int
    selected_min_edge_count: int
    selected_min_graph_count: int
    selected_signature: tuple[str, ...]
    selected_is_diamond: bool
    without_direct_min_edge_count: int
    without_direct_signature: tuple[str, ...]
    without_two_anchors_min_edge_count: int
    without_two_anchors_signature: tuple[str, ...]
    diamond_unique_under_constraints: bool
    fully_derived: bool
    verdict: str


def assess_diamond_cell_selection() -> DiamondCellSelectionAssessment:
    """Assess whether explicit constraints select the saturated diamond cell."""

    selected = minimal_graphs(
        admissible_graphs(
            require_two_anchors=True,
            require_direct_exchange=True,
        )
    )
    without_direct = minimal_graphs(
        admissible_graphs(
            require_two_anchors=True,
            require_direct_exchange=False,
        )
    )
    without_two_anchors = minimal_graphs(
        admissible_graphs(
            require_two_anchors=False,
            require_direct_exchange=True,
        )
    )

    selected_graph = selected[0] if selected else frozenset()
    without_direct_graph = without_direct[0] if without_direct else frozenset()
    without_two_anchors_graph = (
        without_two_anchors[0] if without_two_anchors else frozenset()
    )
    selected_is_diamond = selected_graph == DIAMOND_EDGES
    unique = len(selected) == 1 and selected_is_diamond

    if unique:
        verdict = (
            "two saturated anchors plus direct branch exchange uniquely select "
            "the saturated diamond among four-node spring cells; these two "
            "constraints still need substrate derivation"
        )
    else:
        verdict = "minimal-cell constraints do not uniquely select the diamond"

    return DiamondCellSelectionAssessment(
        total_graphs_scanned=len(powerset_edges()),
        selected_min_edge_count=len(selected_graph),
        selected_min_graph_count=len(selected),
        selected_signature=graph_signature(selected_graph),
        selected_is_diamond=selected_is_diamond,
        without_direct_min_edge_count=len(without_direct_graph),
        without_direct_signature=graph_signature(without_direct_graph),
        without_two_anchors_min_edge_count=len(without_two_anchors_graph),
        without_two_anchors_signature=graph_signature(without_two_anchors_graph),
        diamond_unique_under_constraints=unique,
        fully_derived=False,
        verdict=verdict,
    )
