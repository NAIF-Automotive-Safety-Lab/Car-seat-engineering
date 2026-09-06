import math

import numpy as np

from stiff_medium.cone_two_anchor_origin import (
    assess_two_anchor_origin,
    minimal_neutral_endpoint_count,
    phase_slip_segment_boundary,
    saturation_barrier_curvature,
    single_anchor_boundary,
    total_endpoint_charge,
)


def test_phase_slip_segment_has_neutral_paired_endpoints():
    boundary = phase_slip_segment_boundary()

    assert np.allclose(boundary, [-1.0, 1.0])
    assert len(boundary) == 2
    assert total_endpoint_charge(boundary) == 0.0
    assert minimal_neutral_endpoint_count() == 2


def test_single_anchor_is_forbidden_by_endpoint_neutrality():
    boundary = single_anchor_boundary()

    assert len(boundary) == 1
    assert total_endpoint_charge(boundary) != 0.0


def test_saturation_barrier_gives_finite_but_growing_curvature():
    reference = saturation_barrier_curvature(0.5)
    near_cap = saturation_barrier_curvature(0.99)

    assert math.isfinite(reference)
    assert math.isfinite(near_cap)
    assert near_cap > 1.0e4 * reference


def test_two_anchor_origin_closes_anchor_count_conditionally():
    result = assess_two_anchor_origin()

    assert result.endpoint_count == 2
    assert result.signed_endpoint_charge == 0.0
    assert result.single_anchor_charge != 0.0
    assert result.finite_anchor_compliance
    assert result.two_anchor_topology_selected
    assert abs(result.shared_anchor_exchange_strength - 2.0 / 3.0) < 1.0e-12
    assert abs(result.cone_angle_deg - 45.0) < 1.0e-12
    assert not result.fully_derived
    assert "remaining gap is deriving the phase-slip segment" in result.verdict


def test_invalid_two_anchor_inputs_are_rejected():
    try:
        total_endpoint_charge(np.eye(2))
    except ValueError as exc:
        assert "one-dimensional" in str(exc)
    else:
        raise AssertionError("nonvector boundary should fail")

    for kwargs in ({"sigma_fraction": 1.0}, {"sigma_fraction": -0.1}):
        try:
            saturation_barrier_curvature(**kwargs)
        except ValueError as exc:
            assert "sigma_fraction" in str(exc)
        else:
            raise AssertionError("invalid saturation fraction should fail")

    try:
        saturation_barrier_curvature(0.5, sigma_max=0.0)
    except ValueError as exc:
        assert "sigma_max" in str(exc)
    else:
        raise AssertionError("invalid sigma_max should fail")
