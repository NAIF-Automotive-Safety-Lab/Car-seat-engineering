import math

from stiff_medium.cone_detailed_balance import (
    assess_detailed_balance_closure,
    boltzmann_branch_weight,
    cone_angle_from_detailed_balance,
    exchange_generator,
    involution_commutator_norm,
    stationary_branch_weight,
    swap_involution,
)


def test_swap_involution_squares_to_identity():
    j = swap_involution()

    assert (j @ j == [[1.0, 0.0], [0.0, 1.0]]).all()


def test_symmetric_exchange_generator_commutes_with_swap():
    generator = exchange_generator(rate_plus_to_minus=1.0, rate_minus_to_plus=1.0)

    assert involution_commutator_norm(generator) < 1.0e-12


def test_symmetric_detailed_balance_gives_equal_branch_weight():
    assert abs(stationary_branch_weight(rate_plus_to_minus=1.0, rate_minus_to_plus=1.0) - 0.5) < 1.0e-12
    assert abs(cone_angle_from_detailed_balance(rate_plus_to_minus=1.0, rate_minus_to_plus=1.0) - 45.0) < 1.0e-12


def test_rate_imbalance_breaks_swap_and_shifts_cone():
    generator = exchange_generator(rate_plus_to_minus=1.0, rate_minus_to_plus=1.1)
    angle = cone_angle_from_detailed_balance(rate_plus_to_minus=1.0, rate_minus_to_plus=1.1)

    assert involution_commutator_norm(generator) > 0.0
    assert angle > 46.0


def test_energy_splitting_shifts_boltzmann_weight_and_cone():
    weight = boltzmann_branch_weight(energy_splitting_over_temp=0.2)
    result = assess_detailed_balance_closure(split_energy_over_temp=0.2)

    assert weight < 0.5
    assert result.split_stationary_weight == weight
    assert result.split_minimum_angle_deg < 45.0
    assert abs(result.split_angle_shift_deg) > 1.0


def test_detailed_balance_closure_is_conditional_not_full_derivation():
    result = assess_detailed_balance_closure()

    assert result.detailed_balance_closes_equal_weight
    assert not result.fully_derived
    assert "swap-degenerate generator" in result.verdict


def test_invalid_rates_are_rejected():
    for kwargs in (
        {"rate_plus_to_minus": 0.0, "rate_minus_to_plus": 1.0},
        {"rate_plus_to_minus": 1.0, "rate_minus_to_plus": -1.0},
    ):
        try:
            exchange_generator(**kwargs)
        except ValueError as exc:
            assert "rates" in str(exc)
        else:
            raise AssertionError("invalid rates should be rejected")

    try:
        assess_detailed_balance_closure(rate_imbalance=-1.0)
    except ValueError as exc:
        assert "rate_imbalance" in str(exc)
    else:
        raise AssertionError("invalid rate imbalance should be rejected")
