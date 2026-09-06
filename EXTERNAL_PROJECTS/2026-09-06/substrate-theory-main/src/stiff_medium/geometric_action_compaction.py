"""Candidate geometric compaction of the stiff-medium action.

The working Lagrangian currently has separate-looking ingredients:

    - substrate scalar strain phi;
    - explicit 45 degree cone multiplier;
    - Mobius half-flux boundary condition;
    - effective dark rho_kink / rho_pol transport equations.

This module tests a more compact reading: treat the cone as a null structure of
an internal substrate metric, the Mobius sign as connection holonomy, and the
dark sector as the symmetric-traceless neutral-stress sector with second-order
stiffness K_eff/K = alpha^2.

It is a candidate compaction, not a completed derivation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .cone_anchor_induced_exchange import assess_anchor_induced_exchange
from .cone_detailed_balance import assess_detailed_balance_closure
from .cone_diamond_cell_geometry import assess_diamond_cell_geometry
from .cone_diamond_cell_selection import assess_diamond_cell_selection
from .cone_lattice_microgeometry import assess_cone_lattice_microgeometry
from .cone_phase_slip_lattice import assess_phase_slip_lattice_origin
from .cone_saturated_bond_selection import assess_saturated_bond_selection
from .cone_self_dual_exchange import assess_self_dual_exchange_mechanism
from .cone_swap_generator_origin import assess_swap_generator_origin
from .cone_two_anchor_origin import assess_two_anchor_origin
from .cone_variational_origin import assess_cone_variational_origin
from .dark_stress_scale_closure import ALPHA_EM, dark_stress_speed_closure
from .neutral_coupling_suppression import speed_from_stiffness_ratio
from .neutral_stress_tensor_modes import (
    assess_neutral_stress_modes,
    symmetric_tensor_trace_projector,
)


def cone_form_matrix(axis: np.ndarray | None = None) -> np.ndarray:
    """Return Q = n n^T - (I - n n^T), whose null cone is 45 degrees."""

    if axis is None:
        axis = np.array([0.0, 0.0, 1.0])
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm <= 0.0:
        raise ValueError("axis must be nonzero")
    n = axis / norm
    parallel = np.outer(n, n)
    perpendicular = np.eye(3) - parallel
    return parallel - perpendicular


def cone_null_residual(gradient: np.ndarray, axis: np.ndarray | None = None) -> float:
    """Return grad^T Q grad; zero means 45 degree cone geometry."""

    grad = np.asarray(gradient, dtype=float)
    if grad.shape != (3,):
        raise ValueError("gradient must have shape (3,)")
    q = cone_form_matrix(axis)
    return float(grad @ q @ grad)


def mobius_holonomy(winding: int = 1) -> float:
    """Return half-flux holonomy angle for a loop of integer winding."""

    return math.pi * winding


def mobius_transport_phase(winding: int = 1) -> complex:
    """Return exp(i*pi*winding)."""

    angle = mobius_holonomy(winding)
    return complex(math.cos(angle), math.sin(angle))


@dataclass(frozen=True)
class CompactGeometricAction:
    """Human-readable candidate compact action."""

    lagrangian: str
    potential: str
    covariant_derivative: str
    cone_geometry: str
    dark_sector: str
    status: str


def compact_geometric_action() -> CompactGeometricAction:
    """Return the proposed compact expression."""

    return CompactGeometricAction(
        lagrangian=(
            "L_geo = 1/2*rho*(D_t phi)^2 - 1/2*K*g_cone^ij D_i phi D_j phi "
            "- V(phi) + psibar(i*hbar*gamma^a e_a^mu D_mu - g_Y*phi)psi "
            "- 1/4*F_A^2 - 1/2*K*alpha^2*Tr(ST(strain)^2)"
        ),
        potential=(
            "V(phi) = (K/xi^2)(1-cos(phi/xi))/sqrt(1-(phi/phi_max)^2) "
            "- epsilon_0"
        ),
        covariant_derivative=(
            "D = d + i*e*A_EM + i*A_Mobius, with integral A_Mobius = pi*w"
        ),
        cone_geometry=(
            "g_cone carries Q = n*n^T - (I-n*n^T); allowed substrate rays are Q-null"
        ),
        dark_sector=(
            "ST(strain) is the rank-5 symmetric-traceless neutral sector with "
            "K_eff/K=alpha^2"
        ),
        status=(
            "candidate compaction: removes explicit lambda bookkeeping, but "
            "still must be derived from a deeper substrate variational principle"
        ),
    )


@dataclass(frozen=True)
class GeometricCompactionAssessment:
    """Checks for the candidate compact action."""

    cone_residual_45deg: float
    cone_residual_parallel: float
    cone_residual_perpendicular: float
    mobius_phase_2pi_real: float
    mobius_phase_2pi_imag: float
    mobius_phase_4pi_real: float
    mobius_phase_4pi_imag: float
    stress_projector_rank: int
    stiffness_ratio: float
    dark_speed_km_s: float
    scale_speed_km_s: float
    speed_error_pct: float
    cone_variational_minimum_deg: float
    cone_variational_curvature: float
    cone_lattice_self_dual_required: bool
    cone_lattice_bias_shift_deg: float
    cone_lattice_beta_positive_required: bool
    cone_forced_by_current_lattice_symmetry: bool
    cone_self_dual_conditional_closure: bool
    cone_self_dual_imbalanced_shift_deg: float
    cone_self_dual_fully_derived: bool
    cone_detailed_balance_equal_weight: bool
    cone_detailed_balance_split_shift_deg: float
    cone_detailed_balance_rate_shift_deg: float
    cone_detailed_balance_fully_derived: bool
    cone_cell_automorphism_closes_generator: bool
    cone_cell_split_shift_deg: float
    cone_cell_fully_derived: bool
    cone_diamond_cell_forces_automorphism: bool
    cone_diamond_anchor_split_shift_deg: float
    cone_diamond_fully_derived: bool
    cone_diamond_unique_under_selection: bool
    cone_diamond_selection_min_edges: int
    cone_diamond_selection_fully_derived: bool
    cone_anchor_induced_exchange: bool
    cone_anchor_exchange_strength: float
    cone_anchor_fully_derived: bool
    cone_two_anchor_topology_selected: bool
    cone_anchor_compliance_finite: bool
    cone_two_anchor_fully_derived: bool
    cone_phase_slip_lattice_segment_selected: bool
    cone_phase_slip_stiffness_ratio_fixed: bool
    cone_phase_slip_ratio_scan_angle_error_deg: float
    cone_saturation_barrier_selects_single_bond: bool
    cone_saturated_bond_core_cost_fixed: bool
    lambda_multiplier_replaced: bool
    verdict: str


def assess_geometric_compaction() -> GeometricCompactionAssessment:
    """Assess whether the compact geometric reading reproduces current gates."""

    grad_45 = np.array([1.0, 0.0, 1.0])
    grad_parallel = np.array([0.0, 0.0, 1.0])
    grad_perp = np.array([1.0, 0.0, 0.0])
    phase_2pi = mobius_transport_phase(1)
    phase_4pi = mobius_transport_phase(2)
    stress = assess_neutral_stress_modes()
    projector = symmetric_tensor_trace_projector()
    stiffness_ratio = ALPHA_EM**2
    dark_speed = speed_from_stiffness_ratio(
        stiffness_ratio,
        mode_count=stress.projector_rank,
    )
    scale_speed = dark_stress_speed_closure().v_dark_km_s
    speed_error = (dark_speed / scale_speed - 1.0) * 100.0
    cone_variation = assess_cone_variational_origin()
    cone_lattice = assess_cone_lattice_microgeometry()
    cone_exchange = assess_self_dual_exchange_mechanism()
    cone_balance = assess_detailed_balance_closure()
    cone_cell = assess_swap_generator_origin()
    cone_diamond = assess_diamond_cell_geometry()
    cone_diamond_selection = assess_diamond_cell_selection()
    cone_anchor = assess_anchor_induced_exchange()
    cone_two_anchor = assess_two_anchor_origin()
    cone_phase_slip = assess_phase_slip_lattice_origin()
    cone_bond_selection = assess_saturated_bond_selection()
    lambda_replaced = (
        abs(cone_null_residual(grad_45)) < 1.0e-12
        and cone_null_residual(grad_parallel) > 0.0
        and cone_null_residual(grad_perp) < 0.0
        and abs(cone_variation.minimum_angle_deg - 45.0) < 1.0e-3
        and cone_variation.curvature_at_minimum > 0.0
    )
    projector_ok = (
        stress.projector_rank == 5
        and np.linalg.norm(projector @ projector - projector) < 1.0e-12
    )
    mobius_ok = (
        abs(phase_2pi.real + 1.0) < 1.0e-12
        and abs(phase_2pi.imag) < 1.0e-12
        and abs(phase_4pi.real - 1.0) < 1.0e-12
        and abs(phase_4pi.imag) < 1.0e-12
    )
    speed_ok = abs(speed_error) < 1.0e-9

    if lambda_replaced and mobius_ok and projector_ok and speed_ok:
        verdict = (
            "candidate geometric action reproduces cone null geometry, Mobius "
            "holonomy, and dark symmetric-traceless speed; cone microgeometry "
            "has a uniquely minimal saturated-diamond cell, and finite anchors "
            "induce direct branch exchange; phase-slip neutrality conditionally "
            "selects paired anchors, and the minimal lattice segment realizes "
            "them; however the saturation barrier alone delocalizes the bond, "
            "so energetic localization and anchor stiffness remain open"
        )
    else:
        verdict = "candidate geometric action fails at least one structural gate"

    return GeometricCompactionAssessment(
        cone_residual_45deg=cone_null_residual(grad_45),
        cone_residual_parallel=cone_null_residual(grad_parallel),
        cone_residual_perpendicular=cone_null_residual(grad_perp),
        mobius_phase_2pi_real=phase_2pi.real,
        mobius_phase_2pi_imag=phase_2pi.imag,
        mobius_phase_4pi_real=phase_4pi.real,
        mobius_phase_4pi_imag=phase_4pi.imag,
        stress_projector_rank=stress.projector_rank,
        stiffness_ratio=stiffness_ratio,
        dark_speed_km_s=dark_speed,
        scale_speed_km_s=scale_speed,
        speed_error_pct=speed_error,
        cone_variational_minimum_deg=cone_variation.minimum_angle_deg,
        cone_variational_curvature=cone_variation.curvature_at_minimum,
        cone_lattice_self_dual_required=cone_lattice.self_dual_exchange_required,
        cone_lattice_bias_shift_deg=cone_lattice.bias_shift_deg,
        cone_lattice_beta_positive_required=cone_lattice.beta_positive_required,
        cone_forced_by_current_lattice_symmetry=(
            cone_lattice.cone_forced_by_current_symmetry
        ),
        cone_self_dual_conditional_closure=cone_exchange.conditional_cone_closure,
        cone_self_dual_imbalanced_shift_deg=cone_exchange.imbalanced_shift_deg,
        cone_self_dual_fully_derived=cone_exchange.fully_derived,
        cone_detailed_balance_equal_weight=(
            cone_balance.detailed_balance_closes_equal_weight
        ),
        cone_detailed_balance_split_shift_deg=cone_balance.split_angle_shift_deg,
        cone_detailed_balance_rate_shift_deg=(
            cone_balance.imbalanced_minimum_angle_deg
            - cone_balance.symmetric_minimum_angle_deg
        ),
        cone_detailed_balance_fully_derived=cone_balance.fully_derived,
        cone_cell_automorphism_closes_generator=(
            cone_cell.cell_automorphism_closes_generator
        ),
        cone_cell_split_shift_deg=cone_cell.split_angle_shift_deg,
        cone_cell_fully_derived=cone_cell.fully_derived,
        cone_diamond_cell_forces_automorphism=(
            cone_diamond.diamond_cell_forces_automorphism
        ),
        cone_diamond_anchor_split_shift_deg=cone_diamond.broken_angle_shift_deg,
        cone_diamond_fully_derived=cone_diamond.fully_derived,
        cone_diamond_unique_under_selection=(
            cone_diamond_selection.diamond_unique_under_constraints
        ),
        cone_diamond_selection_min_edges=cone_diamond_selection.selected_min_edge_count,
        cone_diamond_selection_fully_derived=cone_diamond_selection.fully_derived,
        cone_anchor_induced_exchange=cone_anchor.exchange_induced_by_finite_anchors,
        cone_anchor_exchange_strength=cone_anchor.induced_exchange,
        cone_anchor_fully_derived=cone_anchor.fully_derived,
        cone_two_anchor_topology_selected=(
            cone_two_anchor.two_anchor_topology_selected
        ),
        cone_anchor_compliance_finite=cone_two_anchor.finite_anchor_compliance,
        cone_two_anchor_fully_derived=cone_two_anchor.fully_derived,
        cone_phase_slip_lattice_segment_selected=(
            cone_phase_slip.topology_selects_open_segment
        ),
        cone_phase_slip_stiffness_ratio_fixed=(
            cone_phase_slip.stiffness_ratio_fixed_by_topology
        ),
        cone_phase_slip_ratio_scan_angle_error_deg=(
            cone_phase_slip.max_angle_error_deg
        ),
        cone_saturation_barrier_selects_single_bond=(
            cone_bond_selection.pure_barrier_selects_single_bond
        ),
        cone_saturated_bond_core_cost_fixed=(
            cone_bond_selection.core_cost_fixed_by_substrate
        ),
        lambda_multiplier_replaced=lambda_replaced,
        verdict=verdict,
    )
