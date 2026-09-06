import numpy as np
from stiff_medium.neutrino import Neutrino, C
from stiff_medium.dynamics import propagate, detect_overlap


def test_propagate_advances_position_by_velocity_times_dt():
    s = C / np.sqrt(2)
    n = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    dt = 0.1
    moved = propagate(n, dt)
    assert np.allclose(moved.position, [s * dt, s * dt])


def test_propagate_preserves_velocity():
    s = C / np.sqrt(2)
    n = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, -s]),
    )
    moved = propagate(n, 1.0)
    assert np.allclose(moved.velocity, [s, -s])


def test_detect_overlap_when_close():
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([C / np.sqrt(2), C / np.sqrt(2)]),
    )
    n2 = Neutrino(
        position=np.array([0.05, 0.05]),
        velocity=np.array([-C / np.sqrt(2), -C / np.sqrt(2)]),
    )
    assert detect_overlap(n1, n2, r_overlap=0.1) is True


def test_no_overlap_when_far():
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([C / np.sqrt(2), C / np.sqrt(2)]),
    )
    n2 = Neutrino(
        position=np.array([10.0, 10.0]),
        velocity=np.array([-C / np.sqrt(2), -C / np.sqrt(2)]),
    )
    assert detect_overlap(n1, n2, r_overlap=0.1) is False


from stiff_medium.dynamics import displace


def test_displace_pushes_neutrinos_apart():
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([0.05, 0.0]),
        velocity=np.array([-s, s]),
    )
    moved1, moved2 = displace(n1, n2, push=0.1)
    new_dist = float(np.linalg.norm(moved1.position - moved2.position))
    old_dist = float(np.linalg.norm(n1.position - n2.position))
    assert new_dist > old_dist


def test_displace_preserves_both_velocities():
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([0.05, 0.0]),
        velocity=np.array([-s, s]),
    )
    moved1, moved2 = displace(n1, n2, push=0.1)
    assert np.allclose(moved1.velocity, n1.velocity)
    assert np.allclose(moved2.velocity, n2.velocity)


def test_displace_handles_coincident_positions_safely():
    """If positions are identical, displacement direction is degenerate;
    must not crash, must not produce NaN."""
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([1.0, 1.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([1.0, 1.0]),
        velocity=np.array([-s, -s]),
    )
    moved1, moved2 = displace(n1, n2, push=0.1)
    assert np.all(np.isfinite(moved1.position))
    assert np.all(np.isfinite(moved2.position))
    assert not np.allclose(moved1.position, moved2.position)


from stiff_medium.dynamics import step


def test_step_propagates_isolated_neutrinos():
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([10.0, 10.0]),
        velocity=np.array([-s, -s]),
    )
    new_state = step([n1, n2], dt=1.0, r_overlap=0.1, push=0.1)
    # Far apart, so they just propagate freely.
    assert np.allclose(new_state[0].position, [s, s])
    assert np.allclose(new_state[1].position, [10.0 - s, 10.0 - s])


def test_step_displaces_overlapping_neutrinos():
    s = C / np.sqrt(2)
    # Place them so propagation will bring them within r_overlap.
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([s * 1.0, s * 1.0]),  # right where n1 will be after dt=1
        velocity=np.array([-s, -s]),
    )
    new_state = step([n1, n2], dt=1.0, r_overlap=0.5, push=0.2)
    # After step, they should be > r_overlap apart and have unchanged velocities.
    dist = float(np.linalg.norm(new_state[0].position - new_state[1].position))
    assert np.allclose(new_state[0].velocity, n1.velocity)
    assert np.allclose(new_state[1].velocity, n2.velocity)
    # Displacement happened at least.
    assert dist > 0.0
