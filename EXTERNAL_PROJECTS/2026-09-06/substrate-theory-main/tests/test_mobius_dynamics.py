"""Tests for Möbius dynamics — Pauli-via-twist behavior."""

import numpy as np

from stiff_medium.neutrino import C
from stiff_medium.back_reaction import back_reaction_force, project_to_cone
from stiff_medium.mobius_dynamics import (
    MobiusState,
    mobius_aware_force,
    mobius_vverlet_step,
    update_accumulated_azimuth,
)


DT = 0.005
S = C / np.sqrt(2.0)
R_EQ = 0.20
R_CAPTURE = 1.0
K_PUSH = 5.0
K_PULL = 5.0


def _base_force(pos_a, pos_b):
    return back_reaction_force(
        pos_a, pos_b, r_eq=R_EQ, r_capture=R_CAPTURE, k_push=K_PUSH, k_pull=K_PULL
    )


def _force(a, b):
    return mobius_aware_force(a, b, base_force_fn=_base_force)


def _make_pair(initial_phase_a: float, initial_phase_b: float):
    z = np.array([0.0, 0.0, 1.0])
    a = MobiusState(
        position=np.array([-1.5 * R_EQ / 2, 0.0, 0.0]),
        velocity=np.array([0.0, S, S]),
        axis=z,
        accumulated_azimuth=0.0,
        initial_phase=initial_phase_a,
    )
    b = MobiusState(
        position=np.array([1.5 * R_EQ / 2, 0.0, 0.0]),
        velocity=np.array([0.0, -S, S]),
        axis=z,
        accumulated_azimuth=0.0,
        initial_phase=initial_phase_b,
    )
    return a, b


def _run_for_steps(a, b, n_steps):
    for _ in range(n_steps):
        a, b = mobius_vverlet_step(a, b, dt=DT, force_fn=_force, project_to_cone_fn=project_to_cone)
    return a, b


def test_opposite_mobius_pair_binds():
    """e⁺e⁻ analog: opposite slope signs → bound state."""
    a, b = _make_pair(0.0, np.pi)
    a, b = _run_for_steps(a, b, 2000)
    final_dist = float(np.linalg.norm(b.position - a.position))
    assert final_dist < 1.0, (
        f"Opposite-Möbius pair should bind (distance < r_capture=1.0), "
        f"got {final_dist}"
    )


def test_same_mobius_pair_does_not_bind():
    """e⁻e⁻ analog: same slope signs → unbound, distance grows."""
    a, b = _make_pair(0.0, 0.0)
    a, b = _run_for_steps(a, b, 2000)
    final_dist = float(np.linalg.norm(b.position - a.position))
    assert final_dist > 5.0, (
        f"Same-Möbius pair should NOT bind (distance > 5), "
        f"got {final_dist}"
    )


def test_slope_sign_flips_under_orbital_evolution():
    """As the orbit progresses, slope signs flip (Möbius half-integer winding).
    Run a bound pair for many steps; verify both signs and unflipped signs occur."""
    a, b = _make_pair(0.0, np.pi)
    seen_signs_a = set()
    for _ in range(3000):
        a, b = mobius_vverlet_step(a, b, dt=DT, force_fn=_force, project_to_cone_fn=project_to_cone)
        seen_signs_a.add(a.slope_sign)
    assert -1 in seen_signs_a and 1 in seen_signs_a, (
        f"Slope sign should flip during orbit; saw only {seen_signs_a}"
    )


def test_mobius_state_initial_slope():
    """initial_phase=0 → slope +1; initial_phase=π → slope −1."""
    z = np.array([0.0, 0.0, 1.0])
    a = MobiusState(
        position=np.zeros(3),
        velocity=np.array([0.0, S, S]),
        axis=z,
        accumulated_azimuth=0.0,
        initial_phase=0.0,
    )
    b = MobiusState(
        position=np.zeros(3),
        velocity=np.array([0.0, S, S]),
        axis=z,
        accumulated_azimuth=0.0,
        initial_phase=np.pi,
    )
    assert a.slope_sign == 1
    assert b.slope_sign == -1


def test_update_accumulated_azimuth_handles_wraparound():
    """If cone azimuth jumps from near 2π to near 0 (wrap), the accumulated
    delta should be small (positive), not nearly −2π."""
    z = np.array([0.0, 0.0, 1.0])
    # Construct state with velocity at azimuth ≈ 2π − ε
    e1 = np.cross(z, np.array([1.0, 0.0, 0.0]))
    e1 = e1 / float(np.linalg.norm(e1))
    e2 = np.cross(z, e1)
    azimuth_old = 2 * np.pi - 0.05
    v_old = S * z + S * (np.cos(azimuth_old) * e1 + np.sin(azimuth_old) * e2)
    state = MobiusState(
        position=np.zeros(3),
        velocity=v_old,
        axis=z,
        accumulated_azimuth=10.0,  # arbitrary starting cumulative value
        initial_phase=0.0,
    )

    # New velocity at azimuth ≈ 0 + ε (wrapped around)
    azimuth_new = 0.05
    v_new = S * z + S * (np.cos(azimuth_new) * e1 + np.sin(azimuth_new) * e2)

    new_acc = update_accumulated_azimuth(state, v_new)
    delta = new_acc - state.accumulated_azimuth
    # Real change is +0.10, not -2π+0.10 ≈ -6.18
    assert abs(delta) < 0.2, f"Expected small delta, got {delta}"
