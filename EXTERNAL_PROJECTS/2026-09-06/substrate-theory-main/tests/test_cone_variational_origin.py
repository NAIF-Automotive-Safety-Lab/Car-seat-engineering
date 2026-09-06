import math

from stiff_medium.cone_variational_origin import (
    assess_cone_variational_origin,
    balanced_elastic_penalty,
    partition_mismatch,
)


def test_equal_partition_occurs_only_at_45_degrees():
    assert partition_mismatch(0.0) == 1.0
    assert abs(partition_mismatch(math.pi / 4.0)) < 1.0e-12
    assert abs(partition_mismatch(math.pi / 2.0) + 1.0) < 1.0e-12


def test_balanced_elastic_penalty_selects_45_degrees():
    result = assess_cone_variational_origin()

    assert abs(result.minimum_angle_deg - 45.0) < 1.0e-3
    assert result.penalty_at_minimum < 1.0e-12
    assert result.penalty_at_0_deg > 0.0
    assert result.penalty_at_90_deg > 0.0
    assert result.curvature_at_minimum > 0.0


def test_cone_residual_vanishes_at_variational_minimum():
    result = assess_cone_variational_origin()

    assert abs(result.cone_residual_at_minimum) < 1.0e-12


def test_isotropic_quadratic_term_alone_does_not_select_angle():
    result = assess_cone_variational_origin()

    assert not result.selected_without_balance_term
    assert "quartic term" in result.verdict


def test_negative_balance_stiffness_is_rejected():
    try:
        balanced_elastic_penalty(math.pi / 4.0, beta=-1.0)
    except ValueError as exc:
        assert "beta" in str(exc)
    else:
        raise AssertionError("negative beta should be rejected")
