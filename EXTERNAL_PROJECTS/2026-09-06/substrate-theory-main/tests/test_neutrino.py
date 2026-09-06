import numpy as np
import pytest
from stiff_medium.neutrino import Neutrino, C


def test_construct_valid_neutrino():
    n = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([C / np.sqrt(2), C / np.sqrt(2)]),
    )
    assert np.allclose(n.position, [0.0, 0.0])
    assert np.allclose(n.velocity, [C / np.sqrt(2), C / np.sqrt(2)])


def test_reject_non_45_velocity():
    with pytest.raises(ValueError, match="45"):
        Neutrino(
            position=np.array([0.0, 0.0]),
            velocity=np.array([1.0, 0.0]),  # along x-axis, not 45°
        )


def test_reject_wrong_speed():
    with pytest.raises(ValueError, match="magnitude"):
        Neutrino(
            position=np.array([0.0, 0.0]),
            velocity=np.array([1.0, 1.0]),  # 45° but magnitude √2 ≠ C
        )


def test_accept_all_four_45_directions():
    s = C / np.sqrt(2)
    for vx, vy in [(s, s), (s, -s), (-s, s), (-s, -s)]:
        n = Neutrino(
            position=np.array([0.0, 0.0]),
            velocity=np.array([vx, vy]),
        )
        assert np.isclose(np.linalg.norm(n.velocity), C)
