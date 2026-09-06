"""Elastic-cell origin candidate for the swap-degenerate cone generator.

The detailed-balance closure reduces the cone gap to a local exchange generator
G that commutes with the dual-branch swap J.  This module tests the smallest
elastic-cell condition that forces that result.

Model the two dual branch reservoirs with a local stiffness/Hessian matrix H.
If the substrate cell has an exact branch-swap automorphism,

    J^T H J = H,

then the branch energies are degenerate.  Arrhenius/detailed-balance exchange
rates built from that Hessian are equal, so G commutes with J and the stationary
branch weights are exactly 50/50.

This is still conditional: it derives the generator from a cell automorphism,
but the actual substrate lattice/saturation geometry must still supply that
automorphism.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .cone_detailed_balance import (
    cone_angle_from_detailed_balance,
    exchange_generator,
    involution_commutator_norm,
    stationary_branch_weight,
    swap_involution,
)
from .cone_self_dual_exchange import effective_linear_bias


def elastic_cell_hessian(
    *,
    common_stiffness: float = 1.0,
    exchange_coupling: float = 0.25,
    branch_split: float = 0.0,
) -> np.ndarray:
    """Return a two-branch elastic Hessian for the local exchange cell."""

    if common_stiffness <= 0.0:
        raise ValueError("common_stiffness must be positive")
    if exchange_coupling < 0.0:
        raise ValueError("exchange_coupling must be non-negative")
    if common_stiffness <= abs(branch_split) + exchange_coupling:
        raise ValueError("hessian must remain positive definite")
    return np.array(
        [
            [common_stiffness + branch_split, -exchange_coupling],
            [-exchange_coupling, common_stiffness - branch_split],
        ],
        dtype=float,
    )


def cell_automorphism_residual(hessian: np.ndarray) -> float:
    """Return ||J^T H J - H||_F for the branch-swap automorphism."""

    h = np.asarray(hessian, dtype=float)
    if h.shape != (2, 2):
        raise ValueError("hessian must have shape (2, 2)")
    j = swap_involution()
    return float(np.linalg.norm(j.T @ h @ j - h))


def branch_energy_splitting(hessian: np.ndarray) -> float:
    """Return E_plus - E_minus inferred from branch diagonal stiffness."""

    h = np.asarray(hessian, dtype=float)
    if h.shape != (2, 2):
        raise ValueError("hessian must have shape (2, 2)")
    return float(h[0, 0] - h[1, 1])


def rates_from_cell_hessian(
    hessian: np.ndarray,
    *,
    temperature: float = 1.0,
    attempt_rate: float = 1.0,
) -> tuple[float, float]:
    """Return detailed-balance rates (plus->minus, minus->plus) from H."""

    if temperature <= 0.0:
        raise ValueError("temperature must be positive")
    if attempt_rate <= 0.0:
        raise ValueError("attempt_rate must be positive")
    delta = branch_energy_splitting(hessian) / temperature
    return (
        attempt_rate * math.exp(0.5 * delta),
        attempt_rate * math.exp(-0.5 * delta),
    )


def generator_from_cell_hessian(
    hessian: np.ndarray,
    *,
    temperature: float = 1.0,
    attempt_rate: float = 1.0,
) -> np.ndarray:
    """Return the local exchange generator implied by the elastic cell."""

    rate_plus_to_minus, rate_minus_to_plus = rates_from_cell_hessian(
        hessian,
        temperature=temperature,
        attempt_rate=attempt_rate,
    )
    return exchange_generator(
        rate_plus_to_minus=rate_plus_to_minus,
        rate_minus_to_plus=rate_minus_to_plus,
    )


@dataclass(frozen=True)
class SwapGeneratorOriginAssessment:
    """Assessment of the elastic-cell origin candidate."""

    automorphic_residual: float
    automorphic_commutator_norm: float
    automorphic_stationary_weight: float
    automorphic_linear_bias: float
    automorphic_minimum_angle_deg: float
    split_branch_energy_over_temp: float
    split_automorphism_residual: float
    split_commutator_norm: float
    split_stationary_weight: float
    split_linear_bias: float
    split_minimum_angle_deg: float
    split_angle_shift_deg: float
    cell_automorphism_closes_generator: bool
    fully_derived: bool
    verdict: str


def assess_swap_generator_origin(
    *,
    branch_split: float = 0.1,
    temperature: float = 1.0,
) -> SwapGeneratorOriginAssessment:
    """Assess whether cell automorphism forces the swap-degenerate generator."""

    automorphic_h = elastic_cell_hessian(branch_split=0.0)
    automorphic_g = generator_from_cell_hessian(
        automorphic_h,
        temperature=temperature,
    )
    auto_rate_plus, auto_rate_minus = rates_from_cell_hessian(
        automorphic_h,
        temperature=temperature,
    )
    auto_weight = stationary_branch_weight(
        rate_plus_to_minus=auto_rate_plus,
        rate_minus_to_plus=auto_rate_minus,
    )
    auto_angle = cone_angle_from_detailed_balance(
        rate_plus_to_minus=auto_rate_plus,
        rate_minus_to_plus=auto_rate_minus,
    )

    split_h = elastic_cell_hessian(branch_split=branch_split)
    split_g = generator_from_cell_hessian(split_h, temperature=temperature)
    split_rate_plus, split_rate_minus = rates_from_cell_hessian(
        split_h,
        temperature=temperature,
    )
    split_weight = stationary_branch_weight(
        rate_plus_to_minus=split_rate_plus,
        rate_minus_to_plus=split_rate_minus,
    )
    split_angle = cone_angle_from_detailed_balance(
        rate_plus_to_minus=split_rate_plus,
        rate_minus_to_plus=split_rate_minus,
    )
    split_delta = branch_energy_splitting(split_h) / temperature

    closes_generator = (
        cell_automorphism_residual(automorphic_h) < 1.0e-12
        and involution_commutator_norm(automorphic_g) < 1.0e-12
        and abs(auto_weight - 0.5) < 1.0e-12
        and abs(auto_angle - 45.0) < 1.0e-12
        and cell_automorphism_residual(split_h) > 1.0e-6
        and involution_commutator_norm(split_g) > 1.0e-6
        and abs(split_angle - auto_angle) > 1.0
    )

    if closes_generator:
        verdict = (
            "an exact dual-branch cell automorphism forces a swap-degenerate "
            "exchange generator; deriving that automorphism from the actual "
            "substrate cell remains open"
        )
    else:
        verdict = "elastic-cell automorphism does not force the cone generator"

    return SwapGeneratorOriginAssessment(
        automorphic_residual=cell_automorphism_residual(automorphic_h),
        automorphic_commutator_norm=involution_commutator_norm(automorphic_g),
        automorphic_stationary_weight=auto_weight,
        automorphic_linear_bias=effective_linear_bias(branch_weight=auto_weight),
        automorphic_minimum_angle_deg=auto_angle,
        split_branch_energy_over_temp=split_delta,
        split_automorphism_residual=cell_automorphism_residual(split_h),
        split_commutator_norm=involution_commutator_norm(split_g),
        split_stationary_weight=split_weight,
        split_linear_bias=effective_linear_bias(branch_weight=split_weight),
        split_minimum_angle_deg=split_angle,
        split_angle_shift_deg=split_angle - auto_angle,
        cell_automorphism_closes_generator=closes_generator,
        fully_derived=False,
        verdict=verdict,
    )
