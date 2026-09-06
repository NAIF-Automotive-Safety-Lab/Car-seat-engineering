"""Tests for the 3D extension. Validates Neutrino3D's cone constraint
and that the 3D dynamics functions preserve the spec invariants."""

import numpy as np
import pytest

from stiff_medium.neutrino import C
from stiff_medium.three_d import (
    Neutrino3D,
    detect_overlap,
    displace,
    make_on_cone,
    propagate,
    step,
)


# Helpers ---------------------------------------------------------------

def _on_cone_z(azimuth: float = 0.0) -> np.ndarray:
    """Velocity on the 45° cone around +z."""
    return make_on_cone(np.array([0.0, 0.0, 1.0]), azimuth)


# Neutrino3D validation -------------------------------------------------

def test_construct_valid_3d_neutrino():
    n = Neutrino3D(
        position=np.array([0.0, 0.0, 0.0]),
        velocity=_on_cone_z(0.0),
        axis=np.array([0.0, 0.0, 1.0]),
    )
    assert np.isclose(np.linalg.norm(n.velocity), C)


def test_reject_non_unit_axis():
    with pytest.raises(ValueError, match="unit"):
        Neutrino3D(
            position=np.array([0.0, 0.0, 0.0]),
            velocity=_on_cone_z(0.0),
            axis=np.array([0.0, 0.0, 2.0]),
        )


def test_reject_velocity_off_cone():
    """A velocity along the axis (not at 45°) should be rejected."""
    with pytest.raises(ValueError, match="45"):
        Neutrino3D(
            position=np.array([0.0, 0.0, 0.0]),
            velocity=np.array([0.0, 0.0, C]),  # along axis, 0° not 45°
            axis=np.array([0.0, 0.0, 1.0]),
        )


def test_reject_velocity_perpendicular_to_axis():
    """A velocity perpendicular to the axis (90°, not 45°) should be rejected."""
    with pytest.raises(ValueError, match="45"):
        Neutrino3D(
            position=np.array([0.0, 0.0, 0.0]),
            velocity=np.array([C, 0.0, 0.0]),  # ⟂ to axis
            axis=np.array([0.0, 0.0, 1.0]),
        )


def test_make_on_cone_produces_valid_velocity():
    axis = np.array([0.0, 0.0, 1.0])
    for azimuth in (0.0, np.pi / 4, np.pi / 2, np.pi, 1.5 * np.pi):
        v = make_on_cone(axis, azimuth)
        # Magnitude is C
        assert np.isclose(np.linalg.norm(v), C)
        # Angle to axis is 45°
        cos_theta = float(np.dot(v, axis)) / C
        assert np.isclose(cos_theta, 1.0 / np.sqrt(2.0))


def test_make_on_cone_arbitrary_axis():
    """make_on_cone works for an arbitrary (non-axis-aligned) axis."""
    axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    v = make_on_cone(axis, np.pi / 3)
    n = Neutrino3D(
        position=np.array([0.0, 0.0, 0.0]),
        velocity=v,
        axis=axis,
    )
    assert np.isclose(np.linalg.norm(n.velocity), C)


# Dynamics --------------------------------------------------------------

def test_propagate_advances_3d_position_and_preserves_velocity_and_axis():
    axis = np.array([0.0, 0.0, 1.0])
    n = Neutrino3D(
        position=np.array([0.0, 0.0, 0.0]),
        velocity=_on_cone_z(0.0),
        axis=axis,
    )
    moved = propagate(n, dt=0.1)
    expected = n.position + n.velocity * 0.1
    assert np.allclose(moved.position, expected)
    assert np.allclose(moved.velocity, n.velocity)
    assert np.allclose(moved.axis, n.axis)


def test_detect_overlap_3d():
    axis = np.array([0.0, 0.0, 1.0])
    n1 = Neutrino3D(
        position=np.array([0.0, 0.0, 0.0]),
        velocity=_on_cone_z(0.0),
        axis=axis,
    )
    n2 = Neutrino3D(
        position=np.array([0.05, 0.0, 0.0]),
        velocity=_on_cone_z(np.pi),
        axis=axis,
    )
    assert detect_overlap(n1, n2, r_overlap=0.1) is True
    assert detect_overlap(n1, n2, r_overlap=0.01) is False


def test_displace_preserves_velocity_and_axis_3d():
    axis = np.array([0.0, 0.0, 1.0])
    n1 = Neutrino3D(
        position=np.array([0.0, 0.0, 0.0]),
        velocity=_on_cone_z(0.0),
        axis=axis,
    )
    n2 = Neutrino3D(
        position=np.array([0.05, 0.0, 0.0]),
        velocity=_on_cone_z(np.pi),
        axis=axis,
    )
    m1, m2 = displace(n1, n2, push=0.1)
    assert np.allclose(m1.velocity, n1.velocity)
    assert np.allclose(m2.velocity, n2.velocity)
    assert np.allclose(m1.axis, n1.axis)
    assert np.allclose(m2.axis, n2.axis)
    new_dist = float(np.linalg.norm(m1.position - m2.position))
    old_dist = float(np.linalg.norm(n1.position - n2.position))
    assert new_dist > old_dist


def test_step_propagates_isolated_3d_neutrinos():
    axis = np.array([0.0, 0.0, 1.0])
    n1 = Neutrino3D(
        position=np.array([0.0, 0.0, 0.0]),
        velocity=_on_cone_z(0.0),
        axis=axis,
    )
    n2 = Neutrino3D(
        position=np.array([10.0, 10.0, 10.0]),
        velocity=_on_cone_z(np.pi),
        axis=axis,
    )
    new_state = step([n1, n2], dt=1.0, r_overlap=0.1, push=0.1)
    assert np.allclose(new_state[0].position, n1.position + n1.velocity)
    assert np.allclose(new_state[1].position, n2.position + n2.velocity)
