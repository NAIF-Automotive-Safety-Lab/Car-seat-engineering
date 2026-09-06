"""Variational origin candidate for the 45 degree cone geometry.

The compact geometric action replaces the explicit Lagrange multiplier with a
cone metric.  This module asks what local elastic term would make that cone
geometry a stable variational minimum.

Candidate:

    E_balance = beta/4 * (|grad_parallel phi|^2 - |grad_perp phi|^2)^2
                / |grad phi|^2

At fixed gradient magnitude, this is minimized when the longitudinal and
transverse strain magnitudes are equal.  That condition is exactly the
45-degree cone/null condition.

This is a candidate local elastic origin, not a microscopic proof.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def partition_mismatch(theta_rad: float) -> float:
    """Return cos^2(theta) - sin^2(theta) = cos(2 theta)."""

    return math.cos(theta_rad) ** 2 - math.sin(theta_rad) ** 2


def balanced_elastic_penalty(theta_rad: float, beta: float = 1.0) -> float:
    """Return normalized equal-partition penalty at fixed |grad phi|."""

    if beta < 0.0:
        raise ValueError("beta must be non-negative")
    mismatch = partition_mismatch(theta_rad)
    return 0.25 * beta * mismatch**2


def numerical_second_derivative(
    fn,
    x: float,
    dx: float = 1.0e-5,
) -> float:
    """Return centered numerical second derivative."""

    return (fn(x + dx) - 2.0 * fn(x) + fn(x - dx)) / dx**2


def cone_null_residual_local(gradient: np.ndarray) -> float:
    """Return |parallel|^2 - |perpendicular|^2 for the z-axis cone."""

    grad = np.asarray(gradient, dtype=float)
    if grad.shape != (3,):
        raise ValueError("gradient must have shape (3,)")
    parallel = grad[2] ** 2
    perpendicular = grad[0] ** 2 + grad[1] ** 2
    return float(parallel - perpendicular)


@dataclass(frozen=True)
class ConeVariationalAssessment:
    """Assessment of the equal-partition variational candidate."""

    minimum_angle_deg: float
    penalty_at_minimum: float
    penalty_at_0_deg: float
    penalty_at_90_deg: float
    curvature_at_minimum: float
    cone_residual_at_minimum: float
    selected_without_balance_term: bool
    verdict: str


def scan_minimum_angle(num_samples: int = 20001) -> tuple[float, float]:
    """Scan theta in [0, pi/2] and return the minimum angle and penalty."""

    if num_samples < 3:
        raise ValueError("num_samples must be >= 3")
    angles = np.linspace(0.0, 0.5 * math.pi, num_samples)
    penalties = np.array([balanced_elastic_penalty(float(theta)) for theta in angles])
    index = int(np.argmin(penalties))
    return float(angles[index]), float(penalties[index])


def assess_cone_variational_origin() -> ConeVariationalAssessment:
    """Assess whether equal-partition elastic energy selects 45 degrees."""

    theta_min, penalty_min = scan_minimum_angle()
    theta_45 = 0.25 * math.pi
    curvature = numerical_second_derivative(balanced_elastic_penalty, theta_45)
    grad_45 = np.array([math.sin(theta_45), 0.0, math.cos(theta_45)])
    residual = cone_null_residual_local(grad_45)

    # Isotropic quadratic elasticity alone is independent of orientation at
    # fixed |grad phi|, so it does not select a cone angle.
    selected_without_balance = False

    if (
        abs(math.degrees(theta_min) - 45.0) < 1.0e-3
        and penalty_min < 1.0e-12
        and curvature > 0.0
        and abs(residual) < 1.0e-12
    ):
        verdict = (
            "equal-partition elastic mismatch gives a stable 45 degree minimum; "
            "derive beta and this quartic term from substrate microgeometry"
        )
    else:
        verdict = "equal-partition elastic mismatch does not select the cone"

    return ConeVariationalAssessment(
        minimum_angle_deg=math.degrees(theta_min),
        penalty_at_minimum=penalty_min,
        penalty_at_0_deg=balanced_elastic_penalty(0.0),
        penalty_at_90_deg=balanced_elastic_penalty(0.5 * math.pi),
        curvature_at_minimum=curvature,
        cone_residual_at_minimum=residual,
        selected_without_balance_term=selected_without_balance,
        verdict=verdict,
    )
