import numpy as np
from stiff_medium.neutrino import Neutrino, C
from stiff_medium.detector import BoundStateTracker


def _n(x, y, vx, vy):
    return Neutrino(
        position=np.array([x, y], dtype=float),
        velocity=np.array([vx, vy], dtype=float),
    )


def test_tracker_starts_unbound():
    s = C / np.sqrt(2)
    tracker = BoundStateTracker(r_bound=0.5, persistence=10)
    a, b = _n(0, 0, s, s), _n(10, 10, -s, -s)
    assert tracker.update([a, b]) is False


def test_tracker_flags_bound_after_persistence_steps():
    s = C / np.sqrt(2)
    tracker = BoundStateTracker(r_bound=0.5, persistence=3)
    a, b = _n(0, 0, s, s), _n(0.1, 0, -s, -s)
    assert tracker.update([a, b]) is False  # step 1 close
    assert tracker.update([a, b]) is False  # step 2 close
    assert tracker.update([a, b]) is True   # step 3 close — flagged


def test_tracker_resets_on_separation():
    s = C / np.sqrt(2)
    tracker = BoundStateTracker(r_bound=0.5, persistence=3)
    a, b = _n(0, 0, s, s), _n(0.1, 0, -s, -s)
    far_a, far_b = _n(0, 0, s, s), _n(10, 10, -s, -s)
    tracker.update([a, b])
    tracker.update([a, b])
    tracker.update([far_a, far_b])  # separation resets counter
    # Need 3 more close steps to flag.
    assert tracker.update([a, b]) is False
    assert tracker.update([a, b]) is False
    assert tracker.update([a, b]) is True
