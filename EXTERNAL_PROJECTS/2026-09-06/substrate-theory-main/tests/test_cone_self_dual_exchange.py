import math

from stiff_medium.cone_self_dual_exchange import (
    analytic_minimum_angle,
    assess_self_dual_exchange_mechanism,
    dual_branch_energy,
    effective_beta,
    effective_linear_bias,
)


def test_equal_dual_branches_cancel_linear_bias():
    assert abs(effective_linear_bias(branch_weight=0.5)) < 1.0e-12
    assert effective_beta(branch_weight=0.5) > 0.0


def test_balanced_dual_pair_selects_45_degrees():
    theta_min = analytic_minimum_angle(branch_weight=0.5)

    assert abs(math.degrees(theta_min) - 45.0) < 1.0e-12
    assert dual_branch_energy(theta_min, branch_weight=0.5) < dual_branch_energy(
        0.0,
        branch_weight=0.5,
    )


def test_imbalanced_dual_pair_reintroduces_cone_shift():
    result = assess_self_dual_exchange_mechanism(imbalanced_weight=0.55)

    assert result.imbalanced_linear_bias > 0.0
    assert result.imbalanced_minimum_angle_deg > 45.0
    assert result.imbalanced_shift_deg > 1.0


def test_single_branch_does_not_select_45_degrees():
    result = assess_self_dual_exchange_mechanism()

    assert abs(result.single_branch_minimum_angle_deg - 45.0) > 1.0


def test_conditional_closure_is_not_full_derivation():
    result = assess_self_dual_exchange_mechanism()

    assert result.dual_pair_cancels_quadratic_bias
    assert result.beta_positive_from_branch_stability
    assert result.conditional_cone_closure
    assert not result.fully_derived
    assert "equal weight remains" in result.verdict


def test_invalid_branch_parameters_are_rejected():
    for kwargs in ({"branch_weight": -0.1}, {"branch_weight": 1.1}):
        try:
            effective_linear_bias(**kwargs)
        except ValueError as exc:
            assert "branch_weight" in str(exc)
        else:
            raise AssertionError("invalid branch weight should be rejected")

    try:
        effective_beta(branch_stiffness=0.0)
    except ValueError as exc:
        assert "branch_stiffness" in str(exc)
    else:
        raise AssertionError("invalid branch stiffness should be rejected")
