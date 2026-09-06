import numpy as np

from stiff_medium.cone_swap_generator_origin import (
    assess_swap_generator_origin,
    branch_energy_splitting,
    cell_automorphism_residual,
    elastic_cell_hessian,
    generator_from_cell_hessian,
    rates_from_cell_hessian,
)
from stiff_medium.cone_detailed_balance import involution_commutator_norm


def test_automorphic_cell_hessian_commutes_with_branch_swap():
    hessian = elastic_cell_hessian(branch_split=0.0)

    assert cell_automorphism_residual(hessian) < 1.0e-12
    assert abs(branch_energy_splitting(hessian)) < 1.0e-12


def test_automorphic_cell_implies_swap_degenerate_generator():
    hessian = elastic_cell_hessian(branch_split=0.0)
    generator = generator_from_cell_hessian(hessian)
    rate_plus_to_minus, rate_minus_to_plus = rates_from_cell_hessian(hessian)

    assert abs(rate_plus_to_minus - rate_minus_to_plus) < 1.0e-12
    assert involution_commutator_norm(generator) < 1.0e-12


def test_branch_split_breaks_cell_automorphism_and_shifts_generator():
    hessian = elastic_cell_hessian(branch_split=0.1)
    generator = generator_from_cell_hessian(hessian)
    rate_plus_to_minus, rate_minus_to_plus = rates_from_cell_hessian(hessian)

    assert cell_automorphism_residual(hessian) > 0.0
    assert rate_plus_to_minus > rate_minus_to_plus
    assert involution_commutator_norm(generator) > 0.0


def test_swap_generator_origin_closes_conditionally():
    result = assess_swap_generator_origin()

    assert result.cell_automorphism_closes_generator
    assert result.automorphic_stationary_weight == 0.5
    assert abs(result.automorphic_minimum_angle_deg - 45.0) < 1.0e-12
    assert abs(result.split_branch_energy_over_temp - 0.2) < 1.0e-12
    assert result.split_minimum_angle_deg < 45.0
    assert result.split_angle_shift_deg < -1.0
    assert not result.fully_derived
    assert "cell automorphism" in result.verdict


def test_invalid_cell_parameters_are_rejected():
    for kwargs in (
        {"common_stiffness": 0.0},
        {"exchange_coupling": -0.1},
        {"common_stiffness": 0.2, "exchange_coupling": 0.25},
    ):
        try:
            elastic_cell_hessian(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid hessian parameters should be rejected")

    hessian = elastic_cell_hessian()
    for kwargs in ({"temperature": 0.0}, {"attempt_rate": 0.0}):
        try:
            rates_from_cell_hessian(hessian, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid rate parameters should be rejected")


def test_invalid_matrix_shapes_are_rejected():
    bad = np.eye(3)

    for fn in (cell_automorphism_residual, branch_energy_splitting):
        try:
            fn(bad)
        except ValueError as exc:
            assert "shape" in str(exc)
        else:
            raise AssertionError("invalid matrix shape should be rejected")
