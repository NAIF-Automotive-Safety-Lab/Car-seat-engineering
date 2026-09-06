import numpy as np

from stiff_medium.cone_anchor_induced_exchange import (
    analytic_induced_exchange,
    assess_anchor_induced_exchange,
    fixed_anchor_limit_hessian,
    induced_exchange_strength,
    pinned_anchor_hessian,
    schur_branch_hessian,
    shared_anchor_edges,
)
from stiff_medium.cone_swap_generator_origin import cell_automorphism_residual


def test_shared_anchor_graph_has_no_explicit_lt_edge():
    labels = {(edge.node_a, edge.node_b) for edge in shared_anchor_edges()}

    assert ("L", "T") not in labels
    assert labels == {("L", "A"), ("L", "B"), ("T", "A"), ("T", "B")}


def test_schur_complement_induces_effective_branch_exchange():
    full_hessian = pinned_anchor_hessian(
        branch_anchor_stiffness=1.0,
        anchor_pin_stiffness=1.0,
    )
    branch_hessian = schur_branch_hessian(full_hessian)

    assert np.allclose(
        branch_hessian,
        [[4.0 / 3.0, -2.0 / 3.0], [-2.0 / 3.0, 4.0 / 3.0]],
    )
    assert abs(induced_exchange_strength(branch_hessian) - 2.0 / 3.0) < 1.0e-12
    assert abs(analytic_induced_exchange(
        branch_anchor_stiffness=1.0,
        anchor_pin_stiffness=1.0,
    ) - 2.0 / 3.0) < 1.0e-12


def test_fixed_anchor_limit_has_no_induced_exchange():
    branch_hessian = fixed_anchor_limit_hessian(branch_anchor_stiffness=1.0)

    assert induced_exchange_strength(branch_hessian) == 0.0


def test_induced_exchange_preserves_branch_swap_automorphism():
    branch_hessian = schur_branch_hessian(pinned_anchor_hessian())

    assert cell_automorphism_residual(branch_hessian) < 1.0e-12


def test_anchor_induced_exchange_closes_direct_exchange_assumption():
    result = assess_anchor_induced_exchange()

    assert result.exchange_induced_by_finite_anchors
    assert abs(result.induced_exchange - 2.0 / 3.0) < 1.0e-12
    assert result.fixed_anchor_exchange == 0.0
    assert result.stationary_weight == 0.5
    assert abs(result.minimum_angle_deg - 45.0) < 1.0e-12
    assert result.soft_anchor_limit_exchange > result.induced_exchange
    assert result.rigid_anchor_limit_exchange < 1.0e-8
    assert not result.fully_derived
    assert "no longer needs to be a separate primitive" in result.verdict


def test_invalid_anchor_parameters_are_rejected():
    for fn in (shared_anchor_edges, pinned_anchor_hessian, analytic_induced_exchange):
        try:
            fn(branch_anchor_stiffness=0.0)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid branch-anchor stiffness should fail")

    try:
        pinned_anchor_hessian(anchor_pin_stiffness=-1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid anchor pin stiffness should fail")

    try:
        schur_branch_hessian(np.eye(3))
    except ValueError as exc:
        assert "shape" in str(exc)
    else:
        raise AssertionError("invalid Schur input shape should fail")
