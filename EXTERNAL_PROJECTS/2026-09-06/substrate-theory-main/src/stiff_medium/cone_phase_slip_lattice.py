"""Discrete-lattice phase-slip origin for the cone anchor cell.

The previous anchor-origin pass used endpoint neutrality abstractly.  This
module makes that endpoint structure a concrete lattice object: an oriented
saturated 1-chain.  On a lattice, the boundary of the minimal nonzero 1-chain
is a pair of opposite endpoint charges.  That supplies the two anchor sites.

This still does not derive the anchor/branch stiffness ratio.  Topology fixes
the endpoint count; stiffness has to come from the substrate stiffness tensor
or a solved saturated-bond saddle.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from .cone_anchor_induced_exchange import assess_anchor_induced_exchange


Node = tuple[int, int, int]
OrientedEdge = tuple[Node, Node]


def unit_phase_slip_segment() -> tuple[OrientedEdge, ...]:
    """Return the minimal oriented saturated lattice segment."""

    return (((0, 0, 0), (1, 0, 0)),)


def unit_plaquette_loop() -> tuple[OrientedEdge, ...]:
    """Return the minimal closed square phase-slip loop."""

    return (
        ((0, 0, 0), (1, 0, 0)),
        ((1, 0, 0), (1, 1, 0)),
        ((1, 1, 0), (0, 1, 0)),
        ((0, 1, 0), (0, 0, 0)),
    )


def boundary_charges(edges: Iterable[OrientedEdge]) -> dict[Node, float]:
    """Return the discrete boundary charge of an oriented 1-chain."""

    charges: dict[Node, float] = {}
    for start, end in edges:
        charges[start] = charges.get(start, 0.0) - 1.0
        charges[end] = charges.get(end, 0.0) + 1.0
    return {
        node: charge
        for node, charge in charges.items()
        if abs(charge) > 1.0e-12
    }


def net_boundary_charge(charges: Mapping[Node, float]) -> float:
    """Return the total signed boundary charge."""

    return float(sum(charges.values()))


def endpoint_count(charges: Mapping[Node, float]) -> int:
    """Return the number of nonzero endpoint charges."""

    return len(charges)


def is_single_anchor_boundary(charges: Mapping[Node, float]) -> bool:
    """Return True only if the charge map is a forbidden one-anchor boundary."""

    return endpoint_count(charges) == 1 and abs(net_boundary_charge(charges)) > 0.0


def exchange_scan_for_anchor_ratios(
    ratios: Iterable[float],
) -> tuple[tuple[float, float, float], ...]:
    """Return (ratio, induced_exchange, cone_angle) for finite anchor ratios."""

    rows: list[tuple[float, float, float]] = []
    for ratio in ratios:
        if ratio <= 0.0:
            raise ValueError("anchor stiffness ratios must be positive")
        result = assess_anchor_induced_exchange(
            branch_anchor_stiffness=1.0,
            anchor_pin_stiffness=ratio,
        )
        rows.append((float(ratio), result.induced_exchange, result.minimum_angle_deg))
    return tuple(rows)


@dataclass(frozen=True)
class PhaseSlipLatticeAssessment:
    """Assessment of the lattice phase-slip anchor origin."""

    segment_edge_count: int
    segment_endpoint_count: int
    segment_net_charge: float
    loop_edge_count: int
    loop_endpoint_count: int
    single_anchor_is_boundary: bool
    topology_selects_open_segment: bool
    ratio_scan: tuple[tuple[float, float, float], ...]
    min_induced_exchange: float
    max_induced_exchange: float
    max_angle_error_deg: float
    stiffness_ratio_fixed_by_topology: bool
    fully_derived: bool
    verdict: str


def assess_phase_slip_lattice_origin(
    *,
    ratios: tuple[float, ...] = (1.0e-6, 1.0e-3, 1.0, 1.0e3, 1.0e6),
) -> PhaseSlipLatticeAssessment:
    """Assess whether a discrete phase-slip bond supplies the anchor topology."""

    segment = unit_phase_slip_segment()
    segment_boundary = boundary_charges(segment)
    loop = unit_plaquette_loop()
    loop_boundary = boundary_charges(loop)
    single_anchor = {(0, 0, 0): 1.0}
    scan = exchange_scan_for_anchor_ratios(ratios)
    exchanges = tuple(row[1] for row in scan)
    angle_errors = tuple(abs(row[2] - 45.0) for row in scan)

    topology_ok = (
        len(segment) == 1
        and endpoint_count(segment_boundary) == 2
        and abs(net_boundary_charge(segment_boundary)) < 1.0e-12
        and len(loop) == 4
        and endpoint_count(loop_boundary) == 0
        and is_single_anchor_boundary(single_anchor)
    )

    max_angle_error = max(angle_errors) if angle_errors else float("inf")
    min_exchange = min(exchanges) if exchanges else 0.0
    max_exchange = max(exchanges) if exchanges else 0.0
    ratio_fixed = False
    conditional_closure = (
        topology_ok
        and min_exchange > 0.0
        and max_angle_error < 1.0e-12
        and not ratio_fixed
    )

    if conditional_closure:
        verdict = (
            "a discrete saturated phase-slip bond supplies the paired endpoint "
            "anchors and preserves a 45 degree cone for any finite symmetric "
            "anchor stiffness ratio; topology does not fix that ratio"
        )
    else:
        verdict = "discrete phase-slip lattice origin fails the anchor topology gate"

    return PhaseSlipLatticeAssessment(
        segment_edge_count=len(segment),
        segment_endpoint_count=endpoint_count(segment_boundary),
        segment_net_charge=net_boundary_charge(segment_boundary),
        loop_edge_count=len(loop),
        loop_endpoint_count=endpoint_count(loop_boundary),
        single_anchor_is_boundary=is_single_anchor_boundary(single_anchor),
        topology_selects_open_segment=topology_ok,
        ratio_scan=scan,
        min_induced_exchange=min_exchange,
        max_induced_exchange=max_exchange,
        max_angle_error_deg=max_angle_error,
        stiffness_ratio_fixed_by_topology=ratio_fixed,
        fully_derived=False,
        verdict=verdict,
    )
