"""Matplotlib animation for the simulation. v1: 2D scatter of positions,
arrows for velocities, color change when bound state is flagged."""

from typing import Sequence
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

from stiff_medium.neutrino import Neutrino


def animate(
    history: Sequence[Sequence[Neutrino]],
    bound_flags: Sequence[bool],
    xlim: tuple[float, float] = (-5, 5),
    ylim: tuple[float, float] = (-5, 5),
    interval_ms: int = 50,
) -> tuple[plt.Figure, FuncAnimation]:
    """Create animation; caller can plt.show() or save."""
    if len(history) != len(bound_flags):
        raise ValueError("history and bound_flags must be same length")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_title("Stiff-Medium Path C: neutrino dynamics")
    ax.grid(True, alpha=0.3)

    scatter = ax.scatter([], [], s=80, c="steelblue")
    quivers: list = []
    label = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top")

    def init():
        scatter.set_offsets(np.zeros((0, 2)))
        return [scatter, label]

    def update(frame_idx: int):
        for q in quivers:
            q.remove()
        quivers.clear()

        state = history[frame_idx]
        positions = np.array([n.position for n in state])
        velocities = np.array([n.velocity for n in state])

        scatter.set_offsets(positions)
        scatter.set_color("crimson" if bound_flags[frame_idx] else "steelblue")

        q = ax.quiver(
            positions[:, 0], positions[:, 1],
            velocities[:, 0], velocities[:, 1],
            angles="xy", scale_units="xy", scale=2, color="gray", alpha=0.6,
        )
        quivers.append(q)

        label.set_text(
            f"step {frame_idx}/{len(history) - 1}  "
            f"{'BOUND' if bound_flags[frame_idx] else 'free'}"
        )
        return [scatter, q, label]

    anim = FuncAnimation(
        fig, update, frames=len(history),
        init_func=init, interval=interval_ms, blit=False, repeat=False,
    )
    return fig, anim
