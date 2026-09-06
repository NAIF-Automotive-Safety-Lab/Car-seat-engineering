import numpy as np

from stiff_medium.cone_diamond_cell_geometry import (
    assess_diamond_cell_geometry,
    branch_swap_permutation,
    diamond_cell_edges,
    fixed_anchor_branch_hessian,
    graph_automorphism_residual,
    graph_laplacian,
)
from stiff_medium.cone_swap_generator_origin import cell_automorphism_residual


def test_branch_swap_permutation_is_involution():
    permutation = branch_swap_permutation()

    assert np.allclose(permutation @ permutation, np.eye(4))


def test_symmetric_diamond_graph_has_branch_swap_automorphism():
    laplacian = graph_laplacian(diamond_cell_edges(anchor_split=0.0))
    branch_hessian = fixed_anchor_branch_hessian(laplacian)

    assert graph_automorphism_residual(laplacian) < 1.0e-12
    assert cell_automorphism_residual(branch_hessian) < 1.0e-12


def test_symmetric_diamond_branch_hessian_matches_exchange_cell_form():
    laplacian = graph_laplacian(diamond_cell_edges(anchor_split=0.0))
    branch_hessian = fixed_anchor_branch_hessian(laplacian)

    assert np.allclose(branch_hessian, [[1.25, -0.25], [-0.25, 1.25]])


def test_anchor_split_breaks_graph_and_branch_automorphism():
    laplacian = graph_laplacian(diamond_cell_edges(anchor_split=0.05))
    branch_hessian = fixed_anchor_branch_hessian(laplacian)

    assert graph_automorphism_residual(laplacian) > 0.0
    assert cell_automorphism_residual(branch_hessian) > 0.0
    assert abs(branch_hessian[0, 0] - branch_hessian[1, 1] - 0.2) < 1.0e-12


def test_diamond_cell_geometry_closes_conditionally():
    result = assess_diamond_cell_geometry()

    assert result.diamond_cell_forces_automorphism
    assert result.symmetric_stationary_weight == 0.5
    assert abs(result.symmetric_minimum_angle_deg - 45.0) < 1.0e-12
    assert abs(result.broken_branch_energy_over_temp - 0.2) < 1.0e-12
    assert result.broken_minimum_angle_deg < 45.0
    assert result.broken_angle_shift_deg < -1.0
    assert not result.fully_derived
    assert "diamond cell" in result.verdict


def test_invalid_diamond_cell_parameters_are_rejected():
    for kwargs in (
        {"anchor_stiffness": 0.0},
        {"exchange_stiffness": -0.1},
        {"anchor_split": 0.5},
    ):
        try:
            diamond_cell_edges(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid diamond parameters should be rejected")


def test_invalid_graph_inputs_are_rejected():
    bad = np.eye(3)

    for fn in (graph_automorphism_residual, fixed_anchor_branch_hessian):
        try:
            fn(bad)
        except ValueError as exc:
            assert "shape" in str(exc)
        else:
            raise AssertionError("invalid graph shape should be rejected")
