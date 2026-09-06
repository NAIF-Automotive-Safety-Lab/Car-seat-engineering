import numpy as np

from stiff_medium.dark_stress_scale_closure import dark_stress_speed_closure
from stiff_medium.neutral_stress_tensor_modes import (
    assess_neutral_stress_modes,
    projector_rank,
    symmetric_tensor_trace_projector,
)


def test_symmetric_traceless_projector_has_rank_five():
    projector = symmetric_tensor_trace_projector()

    assert projector.shape == (6, 6)
    assert projector_rank(projector) == 5
    assert np.linalg.norm(projector @ projector - projector) < 1.0e-12


def test_trace_mode_is_removed():
    projector = symmetric_tensor_trace_projector()
    trace_direction = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

    assert np.linalg.norm(projector @ trace_direction) < 1.0e-12


def test_neutral_stress_speed_matches_scale_closure():
    result = assess_neutral_stress_modes()
    scale_speed = dark_stress_speed_closure()

    assert result.projector_rank == 5
    assert abs(result.v_dark_km_s - scale_speed.v_dark_km_s) < 1.0e-9
    assert 900.0 < result.v_dark_km_s < 1100.0


def test_neutral_stress_verdict_keeps_alpha_gap_visible():
    result = assess_neutral_stress_modes()

    assert "sqrt(5)" in result.verdict
    assert "alpha" in result.verdict
    assert "derivation" in result.verdict
