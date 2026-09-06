"""Tests for medium back-reaction (spec §5.5).

Verifies cone projection preserves the spec's invariants and that the
back-reaction force has the right structure (push close, pull far).
"""

import numpy as np

from stiff_medium.neutrino import C
from stiff_medium.back_reaction import (
    back_reaction_force,
    project_to_cone,
    vverlet_step,
)


# project_to_cone --------------------------------------------------------

def test_project_to_cone_produces_unit_speed_and_45_angle():
    axis = np.array([0.0, 0.0, 1.0])
    v_in = np.array([0.5, 0.0, 0.5])  # arbitrary off-cone vector
    v_out = project_to_cone(v_in, axis)
    assert np.isclose(np.linalg.norm(v_out), C)
    cos_theta = float(np.dot(v_out, axis)) / C
    assert np.isclose(cos_theta, 1.0 / np.sqrt(2.0))


def test_project_to_cone_preserves_azimuthal_direction():
    """If v already has a clear perpendicular direction, projection should
    preserve that direction (up to magnitude rescaling)."""
    axis = np.array([0.0, 0.0, 1.0])
    # v has perpendicular component along +x
    v_in = np.array([3.0, 0.0, 1.0])
    v_out = project_to_cone(v_in, axis)
    # The perpendicular direction of v_out should still be along +x
    v_along = float(np.dot(v_out, axis)) * axis
    v_perp = v_out - v_along
    assert v_perp[0] > 0  # +x direction preserved
    assert np.isclose(v_perp[1], 0.0)


def test_project_to_cone_handles_axis_aligned_input():
    """Velocity exactly along axis (degenerate: no perpendicular).
    Projection must still produce a valid on-cone vector."""
    axis = np.array([0.0, 0.0, 1.0])
    v_in = np.array([0.0, 0.0, 1.0])  # along axis only
    v_out = project_to_cone(v_in, axis)
    assert np.isclose(np.linalg.norm(v_out), C)
    cos_theta = float(np.dot(v_out, axis)) / C
    assert np.isclose(cos_theta, 1.0 / np.sqrt(2.0))


def test_project_to_cone_arbitrary_axis():
    """Works for non-axis-aligned axes too."""
    axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    v_in = np.array([0.7, 0.2, 0.5])
    v_out = project_to_cone(v_in, axis)
    assert np.isclose(np.linalg.norm(v_out), C)
    cos_theta = float(np.dot(v_out, axis)) / C
    assert np.isclose(cos_theta, 1.0 / np.sqrt(2.0))


# back_reaction_force ----------------------------------------------------

def test_force_is_zero_beyond_r_capture():
    pos_a = np.array([0.0, 0.0, 0.0])
    pos_b = np.array([10.0, 0.0, 0.0])
    f = back_reaction_force(pos_a, pos_b, r_eq=0.2, r_capture=1.0, k_push=5.0, k_pull=5.0)
    assert np.allclose(f, [0.0, 0.0, 0.0])


def test_force_is_repulsive_when_closer_than_r_eq():
    """At d < r_eq, force on A points AWAY from B (negative direction along connecting line)."""
    pos_a = np.array([0.0, 0.0, 0.0])
    pos_b = np.array([0.1, 0.0, 0.0])  # d=0.1 < r_eq=0.2
    f = back_reaction_force(pos_a, pos_b, r_eq=0.2, r_capture=1.0, k_push=5.0, k_pull=5.0)
    # Force on A should be in -x direction (away from B which is at +x).
    assert f[0] < 0


def test_force_is_attractive_when_farther_than_r_eq():
    """At r_eq < d < r_capture, force on A points TOWARD B."""
    pos_a = np.array([0.0, 0.0, 0.0])
    pos_b = np.array([0.5, 0.0, 0.0])  # d=0.5, between r_eq=0.2 and r_capture=1.0
    f = back_reaction_force(pos_a, pos_b, r_eq=0.2, r_capture=1.0, k_push=5.0, k_pull=5.0)
    # Force on A should be in +x direction (toward B at +x).
    assert f[0] > 0


def test_force_zero_at_equilibrium_distance():
    """At d = r_eq exactly, force should be zero."""
    pos_a = np.array([0.0, 0.0, 0.0])
    pos_b = np.array([0.2, 0.0, 0.0])  # d=0.2 = r_eq exactly
    f = back_reaction_force(pos_a, pos_b, r_eq=0.2, r_capture=1.0, k_push=5.0, k_pull=5.0)
    # At d=r_eq, the attractive branch fires (d > r_eq is false but d == r_eq so we're at the boundary).
    # Implementation uses `d > r_eq` for attract, so at d == r_eq, repulsive branch fires with magnitude 0.
    assert np.allclose(f, [0.0, 0.0, 0.0])


def test_force_symmetric_about_r_eq():
    """Pull at d = r_eq + δ has same magnitude as push at d = r_eq - δ
    when k_push = k_pull."""
    delta = 0.05
    r_eq = 0.2
    pos_a = np.array([0.0, 0.0, 0.0])
    pos_close = np.array([r_eq - delta, 0.0, 0.0])
    pos_far = np.array([r_eq + delta, 0.0, 0.0])

    f_push = back_reaction_force(pos_a, pos_close, r_eq=r_eq, r_capture=1.0, k_push=5.0, k_pull=5.0)
    f_pull = back_reaction_force(pos_a, pos_far, r_eq=r_eq, r_capture=1.0, k_push=5.0, k_pull=5.0)

    # Magnitudes equal, signs opposite (pull is +x toward far B; push is -x away from close B).
    assert np.isclose(abs(f_push[0]), abs(f_pull[0]))
    assert np.sign(f_push[0]) != np.sign(f_pull[0])


# vverlet_step ----------------------------------------------------------

def test_vverlet_step_preserves_speed_via_cone_projection():
    """After one step, |v_A| and |v_B| should still equal C
    because cone projection enforces the magnitude."""
    axis = np.array([0.0, 0.0, 1.0])
    s = C / np.sqrt(2.0)

    pos_a = np.array([-0.15, 0.0, 0.0])
    vel_a = np.array([0.0, s, s])
    pos_b = np.array([0.15, 0.0, 0.0])
    vel_b = np.array([0.0, -s, s])

    def force_fn(pa, pb):
        return back_reaction_force(pa, pb, r_eq=0.2, r_capture=1.0, k_push=5.0, k_pull=5.0)

    new_pa, new_va, new_pb, new_vb = vverlet_step(
        pos_a, vel_a, axis, pos_b, vel_b, axis, dt=0.005, force_fn=force_fn
    )

    assert np.isclose(np.linalg.norm(new_va), C)
    assert np.isclose(np.linalg.norm(new_vb), C)


def test_vverlet_step_keeps_45_angle_to_axis():
    """After one step, both velocities should still be at 45° to their axes."""
    axis = np.array([0.0, 0.0, 1.0])
    s = C / np.sqrt(2.0)

    pos_a = np.array([-0.15, 0.0, 0.0])
    vel_a = np.array([0.0, s, s])
    pos_b = np.array([0.15, 0.0, 0.0])
    vel_b = np.array([0.0, -s, s])

    def force_fn(pa, pb):
        return back_reaction_force(pa, pb, r_eq=0.2, r_capture=1.0, k_push=5.0, k_pull=5.0)

    new_pa, new_va, new_pb, new_vb = vverlet_step(
        pos_a, vel_a, axis, pos_b, vel_b, axis, dt=0.005, force_fn=force_fn
    )

    cos_a = float(np.dot(new_va, axis)) / C
    cos_b = float(np.dot(new_vb, axis)) / C
    target = 1.0 / np.sqrt(2.0)
    assert np.isclose(cos_a, target)
    assert np.isclose(cos_b, target)
