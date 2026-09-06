"""Candidate origin for the two finite saturated anchors.

The anchor-induced exchange result removes the direct L-T spring as a primitive.
The remaining local-cell assumption is now the shared pair of finite saturated
anchors.  This module tests one possible origin:

    - a local saturated phase-slip segment has two endpoint charges;
    - endpoint neutrality forbids an unpaired single anchor;
    - the saturation cap has finite curvature below the exact cap, so anchors
      can be stiff but not infinitely rigid.

This is still conditional.  It derives why two finite anchors are natural if
the local saturated object is a finite phase-slip segment, but it does not yet
derive that segment from a 3D substrate lattice or stiffness tensor.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .cone_anchor_induced_exchange import assess_anchor_induced_exchange


def phase_slip_segment_boundary() -> np.ndarray:
    """Return signed endpoint charges for an oriented saturated segment."""

    return np.array([-1.0, 1.0], dtype=float)


def total_endpoint_charge(boundary: np.ndarray) -> float:
    """Return the net signed endpoint charge."""

    charges = np.asarray(boundary, dtype=float)
    if charges.ndim != 1:
        raise ValueError("boundary must be one-dimensional")
    return float(np.sum(charges))


def minimal_neutral_endpoint_count() -> int:
    """Return the minimal nonzero endpoint count with zero signed charge."""

    return 2


def single_anchor_boundary() -> np.ndarray:
    """Return the forbidden one-anchor boundary-control case."""

    return np.array([1.0], dtype=float)


def saturation_barrier_curvature(
    sigma_fraction: float,
    *,
    sigma_max: float = 0.5,
) -> float:
    """Return d2/dsigma2 of (1 - (sigma/sigma_max)^2)^-1/2.

    The value is finite below the cap and diverges only at the exact cap.
    """

    if sigma_max <= 0.0:
        raise ValueError("sigma_max must be positive")
    if not 0.0 <= sigma_fraction < 1.0:
        raise ValueError("sigma_fraction must be in [0, 1)")

    x = sigma_fraction
    numerator = 1.0 + 2.0 * x * x
    denominator = (1.0 - x * x) ** 2.5
    return numerator / (sigma_max * sigma_max * denominator)


@dataclass(frozen=True)
class TwoAnchorOriginAssessment:
    """Assessment for the two-anchor origin mechanism."""

    endpoint_count: int
    signed_endpoint_charge: float
    single_anchor_charge: float
    reference_barrier_curvature: float
    near_cap_barrier_curvature: float
    finite_anchor_compliance: bool
    two_anchor_topology_selected: bool
    shared_anchor_exchange_strength: float
    cone_angle_deg: float
    fully_derived: bool
    verdict: str


def assess_two_anchor_origin() -> TwoAnchorOriginAssessment:
    """Assess whether phase-slip neutrality supplies the two-anchor cell."""

    boundary = phase_slip_segment_boundary()
    single = single_anchor_boundary()
    endpoint_count = len(boundary)
    signed_charge = total_endpoint_charge(boundary)
    single_charge = total_endpoint_charge(single)
    reference_curvature = saturation_barrier_curvature(0.5)
    near_cap_curvature = saturation_barrier_curvature(0.99)
    anchor_exchange = assess_anchor_induced_exchange()

    finite_compliance = (
        math.isfinite(reference_curvature)
        and math.isfinite(near_cap_curvature)
        and near_cap_curvature > reference_curvature
    )
    topology_selected = (
        endpoint_count == minimal_neutral_endpoint_count()
        and abs(signed_charge) < 1.0e-12
        and abs(single_charge) > 0.0
    )
    conditional_closure = (
        topology_selected
        and finite_compliance
        and anchor_exchange.exchange_induced_by_finite_anchors
    )

    if conditional_closure:
        verdict = (
            "neutral phase-slip endpoints select a two-anchor pair and the "
            "saturation barrier supplies finite anchor compliance below the "
            "exact cap; the remaining gap is deriving the phase-slip segment "
            "and stiffness ratio from the substrate lattice"
        )
    else:
        verdict = "phase-slip endpoint mechanism does not select the anchor cell"

    return TwoAnchorOriginAssessment(
        endpoint_count=endpoint_count,
        signed_endpoint_charge=signed_charge,
        single_anchor_charge=single_charge,
        reference_barrier_curvature=reference_curvature,
        near_cap_barrier_curvature=near_cap_curvature,
        finite_anchor_compliance=finite_compliance,
        two_anchor_topology_selected=topology_selected,
        shared_anchor_exchange_strength=anchor_exchange.induced_exchange,
        cone_angle_deg=anchor_exchange.minimum_angle_deg,
        fully_derived=False,
        verdict=verdict,
    )
