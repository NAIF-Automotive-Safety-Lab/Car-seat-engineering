from stiff_medium.geometric_action_compaction import (
    assess_geometric_compaction,
    compact_geometric_action,
    cone_null_residual,
    mobius_transport_phase,
)


def test_cone_geometry_replaces_explicit_multiplier_gate():
    assert abs(cone_null_residual([1.0, 0.0, 1.0])) < 1.0e-12
    assert cone_null_residual([0.0, 0.0, 1.0]) > 0.0
    assert cone_null_residual([1.0, 0.0, 0.0]) < 0.0


def test_mobius_connection_gives_spinor_sign_flip():
    phase_2pi = mobius_transport_phase(1)
    phase_4pi = mobius_transport_phase(2)

    assert abs(phase_2pi.real + 1.0) < 1.0e-12
    assert abs(phase_2pi.imag) < 1.0e-12
    assert abs(phase_4pi.real - 1.0) < 1.0e-12
    assert abs(phase_4pi.imag) < 1.0e-12


def test_compact_action_contains_current_open_structures():
    action = compact_geometric_action()

    assert "g_cone" in action.lagrangian
    assert "A_Mobius" in action.covariant_derivative
    assert "K_eff/K=alpha^2" in action.dark_sector
    assert "candidate" in action.status


def test_geometric_compaction_reproduces_dark_speed_gate():
    result = assess_geometric_compaction()

    assert result.lambda_multiplier_replaced
    assert abs(result.cone_variational_minimum_deg - 45.0) < 1.0e-3
    assert result.cone_variational_curvature > 0.0
    assert result.cone_lattice_self_dual_required
    assert result.cone_lattice_beta_positive_required
    assert result.cone_lattice_bias_shift_deg > 1.0
    assert not result.cone_forced_by_current_lattice_symmetry
    assert result.cone_self_dual_conditional_closure
    assert result.cone_self_dual_imbalanced_shift_deg > 1.0
    assert not result.cone_self_dual_fully_derived
    assert result.cone_detailed_balance_equal_weight
    assert result.cone_detailed_balance_split_shift_deg < -1.0
    assert result.cone_detailed_balance_rate_shift_deg > 1.0
    assert not result.cone_detailed_balance_fully_derived
    assert result.cone_cell_automorphism_closes_generator
    assert result.cone_cell_split_shift_deg < -1.0
    assert not result.cone_cell_fully_derived
    assert result.cone_diamond_cell_forces_automorphism
    assert result.cone_diamond_anchor_split_shift_deg < -1.0
    assert not result.cone_diamond_fully_derived
    assert result.cone_diamond_unique_under_selection
    assert result.cone_diamond_selection_min_edges == 5
    assert not result.cone_diamond_selection_fully_derived
    assert result.cone_anchor_induced_exchange
    assert abs(result.cone_anchor_exchange_strength - 2.0 / 3.0) < 1.0e-12
    assert not result.cone_anchor_fully_derived
    assert result.cone_two_anchor_topology_selected
    assert result.cone_anchor_compliance_finite
    assert not result.cone_two_anchor_fully_derived
    assert result.cone_phase_slip_lattice_segment_selected
    assert not result.cone_phase_slip_stiffness_ratio_fixed
    assert result.cone_phase_slip_ratio_scan_angle_error_deg < 1.0e-12
    assert not result.cone_saturation_barrier_selects_single_bond
    assert not result.cone_saturated_bond_core_cost_fixed
    assert result.stress_projector_rank == 5
    assert abs(result.speed_error_pct) < 1.0e-9
    assert "saturation barrier alone delocalizes" in result.verdict
