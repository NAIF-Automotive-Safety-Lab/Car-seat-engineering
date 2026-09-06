import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for tests

from stiff_medium.neutrino import Neutrino, C
from stiff_medium.visualize import animate


def test_animate_returns_figure_without_error():
    s = C / np.sqrt(2)
    history = [
        [
            Neutrino(np.array([float(i) * 0.1, 0.0]), np.array([s, s])),
            Neutrino(np.array([1.0 - float(i) * 0.1, 0.0]), np.array([-s, s])),
        ]
        for i in range(20)
    ]
    bound_flags = [False] * 15 + [True] * 5

    fig, anim = animate(history, bound_flags, xlim=(-2, 2), ylim=(-2, 2))
    assert fig is not None
    assert anim is not None
