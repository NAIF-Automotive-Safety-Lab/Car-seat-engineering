from stiff_medium.cone_phase_slip_lattice import (
    assess_phase_slip_lattice_origin,
    boundary_charges,
    endpoint_count,
    exchange_scan_for_anchor_ratios,
    is_single_anchor_boundary,
    net_boundary_charge,
    unit_phase_slip_segment,
    unit_plaquette_loop,
)


def test_unit_phase_slip_segment_has_two_endpoint_anchors():
    charges = boundary_charges(unit_phase_slip_segment())

    assert charges == {(0, 0, 0): -1.0, (1, 0, 0): 1.0}
    assert endpoint_count(charges) == 2
    assert net_boundary_charge(charges) == 0.0


def test_closed_plaquette_loop_has_no_endpoint_anchors():
    charges = boundary_charges(unit_plaquette_loop())

    assert charges == {}
    assert endpoint_count(charges) == 0
    assert len(unit_plaquette_loop()) == 4


def test_single_anchor_charge_is_not_a_valid_chain_boundary():
    charges = {(0, 0, 0): 1.0}

    assert is_single_anchor_boundary(charges)
    assert net_boundary_charge(charges) == 1.0


def test_exchange_scan_keeps_cone_angle_for_any_finite_symmetric_ratio():
    scan = exchange_scan_for_anchor_ratios((1.0e-6, 1.0, 1.0e6))
    exchanges = [row[1] for row in scan]
    angles = [row[2] for row in scan]

    assert exchanges[0] > exchanges[1] > exchanges[2] > 0.0
    assert all(abs(angle - 45.0) < 1.0e-12 for angle in angles)


def test_phase_slip_lattice_origin_closes_topology_but_not_ratio():
    result = assess_phase_slip_lattice_origin()

    assert result.segment_edge_count == 1
    assert result.segment_endpoint_count == 2
    assert result.segment_net_charge == 0.0
    assert result.loop_edge_count == 4
    assert result.loop_endpoint_count == 0
    assert result.single_anchor_is_boundary
    assert result.topology_selects_open_segment
    assert result.min_induced_exchange > 0.0
    assert result.max_induced_exchange < 1.0
    assert result.max_angle_error_deg < 1.0e-12
    assert not result.stiffness_ratio_fixed_by_topology
    assert not result.fully_derived
    assert "topology does not fix that ratio" in result.verdict


def test_invalid_ratio_scan_rejects_nonpositive_ratios():
    try:
        exchange_scan_for_anchor_ratios((1.0, 0.0))
    except ValueError as exc:
        assert "positive" in str(exc)
    else:
        raise AssertionError("nonpositive anchor ratio should fail")
