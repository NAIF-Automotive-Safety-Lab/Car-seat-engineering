from stiff_medium.dark_stress_scale_closure import ALPHA_EM, dark_stress_speed_closure
from stiff_medium.neutral_coupling_suppression import (
    assess_neutral_coupling_suppression,
    required_stiffness_ratio_for_speed,
    speed_from_stiffness_ratio,
)


def test_required_stiffness_ratio_is_alpha_squared():
    target = dark_stress_speed_closure().v_dark_km_s
    required = required_stiffness_ratio_for_speed(target, mode_count=5)

    assert abs(required / ALPHA_EM**2 - 1.0) < 1.0e-12


def test_second_order_stiffness_reproduces_dark_speed():
    target = dark_stress_speed_closure().v_dark_km_s
    speed = speed_from_stiffness_ratio(ALPHA_EM**2, mode_count=5)

    assert abs(speed - target) < 1.0e-9


def test_unsuppressed_and_linear_stiffness_are_too_fast():
    result = assess_neutral_coupling_suppression()
    verdicts = {candidate.name: candidate.verdict for candidate in result.candidates}

    assert verdicts["unsuppressed neutral stiffness: K_eff/K = 1"] == "too fast"
    assert verdicts["first-order neutral stiffness: K_eff/K = alpha"] == "too fast"
    assert (
        verdicts["second-order neutral stiffness: K_eff/K = alpha^2"]
        == "matches dark-stress speed"
    )


def test_neutral_coupling_verdict_keeps_derivation_gap_visible():
    result = assess_neutral_coupling_suppression()

    assert abs(result.stiffness_error_pct) < 1.0e-9
    assert "second-order neutral stiffness" in result.verdict
    assert "derive" in result.verdict
