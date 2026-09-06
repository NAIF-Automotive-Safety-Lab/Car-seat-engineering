"""Medium back-reaction (spec §5.5).

Implements the centripetal companion to the centrifugal displace rule:
- d < r_eq: repulsive push
- r_eq < d < r_capture: attractive pull
- d > r_capture: no interaction

Plus a 45° cone projection that keeps velocities on the spec's
required cone after the back-reaction modifies them.
"""

from typing import Callable
import numpy as np

from stiff_medium.neutrino import C


def project_to_cone(v: np.ndarray, axis: np.ndarray, speed: float = C) -> np.ndarray:
    """Project a 3D velocity onto the 45° cone around `axis`.

    Result has |v_new| = `speed` and angle 45° to axis. The azimuthal
    direction (rotation around the axis) is taken from the input v's
    perpendicular component. If v has no perpendicular component, an
    arbitrary perpendicular is chosen.
    """
    s = speed / np.sqrt(2.0)
    v_along = float(np.dot(v, axis)) * axis
    v_perp = v - v_along
    perp_norm = float(np.linalg.norm(v_perp))
    if perp_norm < 1e-9:
        # Pick an arbitrary perpendicular to axis.
        if abs(axis[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])
        v_perp = ref - float(np.dot(ref, axis)) * axis
        v_perp = v_perp / float(np.linalg.norm(v_perp))
    else:
        v_perp = v_perp / perp_norm
    return s * axis + s * v_perp


def back_reaction_force(
    pos_a: np.ndarray,
    pos_b: np.ndarray,
    r_eq: float,
    r_capture: float,
    k_push: float,
    k_pull: float,
) -> np.ndarray:
    """Return the back-reaction force on A due to B. Force on B is the
    negative of this (Newton's third law).

    Spring-like potential with minimum at r_eq:
    - d < r_eq: repulsive (push apart), magnitude k_push * (r_eq - d).
    - r_eq < d < r_capture: attractive (pull together), magnitude k_pull * (d - r_eq).
    - d > r_capture or d ≈ 0: zero.
    """
    diff = pos_b - pos_a
    d = float(np.linalg.norm(diff))
    if d > r_capture or d < 1e-12:
        return np.zeros(3)
    unit = diff / d
    if d > r_eq:
        return k_pull * (d - r_eq) * unit  # attractive: A pulled toward B
    else:
        return k_push * (r_eq - d) * (-unit)  # repulsive: A pushed away from B


def vverlet_step(
    pos_a: np.ndarray,
    vel_a: np.ndarray,
    axis_a: np.ndarray,
    pos_b: np.ndarray,
    vel_b: np.ndarray,
    axis_b: np.ndarray,
    dt: float,
    force_fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    mass_a: float = 1.0,
    mass_b: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One velocity-Verlet step with cone projection at every velocity update.

    Returns updated (pos_a, vel_a, pos_b, vel_b). Axes are unchanged.
    `mass_a` and `mass_b` allow unequal-mass dynamics (acceleration = force/mass).
    Defaults to equal unit masses (the original behavior).
    """
    f_a = force_fn(pos_a, pos_b)
    f_b = -f_a
    a_a = f_a / mass_a
    a_b = f_b / mass_b

    half_vel_a = project_to_cone(vel_a + 0.5 * a_a * dt, axis_a)
    half_vel_b = project_to_cone(vel_b + 0.5 * a_b * dt, axis_b)

    new_pos_a = pos_a + half_vel_a * dt
    new_pos_b = pos_b + half_vel_b * dt

    new_f_a = force_fn(new_pos_a, new_pos_b)
    new_f_b = -new_f_a
    new_a_a = new_f_a / mass_a
    new_a_b = new_f_b / mass_b

    new_vel_a = project_to_cone(half_vel_a + 0.5 * new_a_a * dt, axis_a)
    new_vel_b = project_to_cone(half_vel_b + 0.5 * new_a_b * dt, axis_b)

    return new_pos_a, new_vel_a, new_pos_b, new_vel_b
