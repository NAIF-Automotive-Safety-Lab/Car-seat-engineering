"""Symmetry audit for the cone equal-partition quartic.

The variational cone candidate uses

    (|grad_parallel phi|^2 - |grad_perp phi|^2)^2

to select a 45 degree propagation cone.  This module asks whether that term is
the first orientation-selecting invariant allowed by local substrate symmetry.

Result: ordinary axial symmetry is not enough.  It allows a lower-order
anisotropic bias, |grad_parallel phi|^2 - |grad_perp phi|^2, which would select
parallel or transverse propagation.  The quartic becomes the first selector only
after imposing a self-dual exchange symmetry between the longitudinal and
transverse strain reservoirs.

This tightens the gap: derive the self-dual exchange symmetry and beta > 0, or
the 45 degree cone is not forced.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def parallel_perpendicular_squares(
    gradient: np.ndarray | list[float] | tuple[float, float, float],
    axis: np.ndarray | list[float] | tuple[float, float, float] | None = None,
) -> tuple[float, float]:
    """Return (parallel^2, perpendicular^2) relative to an internal axis."""

    grad = np.asarray(gradient, dtype=float)
    if grad.shape != (3,):
        raise ValueError("gradient must have shape (3,)")
    if axis is None:
        axis = np.array([0.0, 0.0, 1.0])
    axis_arr = np.asarray(axis, dtype=float)
    if axis_arr.shape != (3,):
        raise ValueError("axis must have shape (3,)")
    norm = np.linalg.norm(axis_arr)
    if norm <= 0.0:
        raise ValueError("axis must be nonzero")
    n = axis_arr / norm
    parallel = float(np.dot(grad, n))
    grad_sq = float(np.dot(grad, grad))
    parallel_sq = parallel * parallel
    perpendicular_sq = grad_sq - parallel_sq
    return parallel_sq, perpendicular_sq


def normalized_mismatch(theta_rad: float) -> float:
    """Return (parallel^2 - perpendicular^2) at fixed |grad phi| = 1."""

    return math.cos(2.0 * theta_rad)


def orientation_energy(
    theta_rad: float,
    *,
    linear_bias: float = 0.0,
    beta: float = 1.0,
) -> float:
    """Return fixed-magnitude orientation energy through quartic order.

    The linear-bias term is quadratic in the gradient.  The beta term is
    quartic in the gradient.  A nonzero linear bias is allowed unless the
    substrate has the self-dual reservoir exchange p <-> q.
    """

    mismatch = normalized_mismatch(theta_rad)
    return linear_bias * mismatch + 0.25 * beta * mismatch * mismatch


def scan_orientation_minimum(
    *,
    linear_bias: float = 0.0,
    beta: float = 1.0,
    num_samples: int = 40001,
) -> tuple[float, float]:
    """Return the minimum angle and energy on theta in [0, pi/2]."""

    if num_samples < 3:
        raise ValueError("num_samples must be >= 3")
    angles = np.linspace(0.0, 0.5 * math.pi, num_samples)
    energies = np.array(
        [
            orientation_energy(float(theta), linear_bias=linear_bias, beta=beta)
            for theta in angles
        ]
    )
    index = int(np.argmin(energies))
    return float(angles[index]), float(energies[index])


def numerical_second_derivative(fn, x: float, dx: float = 1.0e-5) -> float:
    """Return centered numerical second derivative."""

    return (fn(x + dx) - 2.0 * fn(x) + fn(x - dx)) / dx**2


def allowed_orientation_selectors(*, self_dual_exchange: bool) -> tuple[str, ...]:
    """Return lowest local orientation selectors through quartic order.

    Let p = |grad_parallel phi|^2 and q = |grad_perp phi|^2.  Axial symmetry and
    gradient reversal permit functions of p and q.  Self-dual exchange adds
    p <-> q symmetry, removing odd powers of the mismatch m = p - q.
    """

    if self_dual_exchange:
        return (
            "s = p + q (constant at fixed gradient magnitude)",
            "m^2 = (p - q)^2 (quartic, first orientation selector)",
        )
    return (
        "s = p + q (constant at fixed gradient magnitude)",
        "m = p - q (quadratic anisotropic selector)",
        "m^2 = (p - q)^2 (quartic selector)",
    )


@dataclass(frozen=True)
class ConeLatticeMicrogeometryAssessment:
    """Symmetry assessment for the cone quartic origin."""

    quadratic_bias_allowed_without_dual: bool
    self_dual_exchange_required: bool
    lowest_selector_without_dual: str
    lowest_selector_with_dual: str
    quartic_minimum_angle_deg: float
    quartic_curvature_at_minimum: float
    biased_minimum_angle_deg: float
    bias_shift_deg: float
    negative_beta_minimum_angle_deg: float
    beta_positive_required: bool
    cone_forced_by_current_symmetry: bool
    verdict: str


def assess_cone_lattice_microgeometry() -> ConeLatticeMicrogeometryAssessment:
    """Assess the symmetry status of the cone quartic."""

    no_dual_terms = allowed_orientation_selectors(self_dual_exchange=False)
    dual_terms = allowed_orientation_selectors(self_dual_exchange=True)

    theta_quartic, _ = scan_orientation_minimum(linear_bias=0.0, beta=1.0)
    theta_biased, _ = scan_orientation_minimum(linear_bias=0.05, beta=1.0)
    theta_negative_beta, _ = scan_orientation_minimum(linear_bias=0.0, beta=-1.0)
    theta_45 = 0.25 * math.pi
    curvature = numerical_second_derivative(
        lambda theta: orientation_energy(theta, linear_bias=0.0, beta=1.0),
        theta_45,
    )

    bias_shift = math.degrees(theta_biased - theta_quartic)
    quadratic_bias_allowed = any("m = p - q" in term for term in no_dual_terms)
    self_dual_required = "m = p - q" not in " ".join(dual_terms)
    beta_positive_required = abs(math.degrees(theta_negative_beta) - 45.0) > 1.0
    cone_forced = (
        not quadratic_bias_allowed
        and self_dual_required
        and curvature > 0.0
    )

    if quadratic_bias_allowed and self_dual_required and beta_positive_required:
        verdict = (
            "ordinary axial lattice symmetry does not force the cone; the "
            "quartic is the first selector only with self-dual longitudinal/"
            "transverse exchange and positive beta"
        )
    else:
        verdict = "symmetry audit did not isolate the cone-origin gap"

    return ConeLatticeMicrogeometryAssessment(
        quadratic_bias_allowed_without_dual=quadratic_bias_allowed,
        self_dual_exchange_required=self_dual_required,
        lowest_selector_without_dual=no_dual_terms[1],
        lowest_selector_with_dual=dual_terms[1],
        quartic_minimum_angle_deg=math.degrees(theta_quartic),
        quartic_curvature_at_minimum=curvature,
        biased_minimum_angle_deg=math.degrees(theta_biased),
        bias_shift_deg=bias_shift,
        negative_beta_minimum_angle_deg=math.degrees(theta_negative_beta),
        beta_positive_required=beta_positive_required,
        cone_forced_by_current_symmetry=cone_forced,
        verdict=verdict,
    )
