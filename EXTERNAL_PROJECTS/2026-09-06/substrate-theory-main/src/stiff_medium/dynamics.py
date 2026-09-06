"""Dynamics rules. See spec §5 for the load-bearing displacement rule."""

import numpy as np
from stiff_medium.neutrino import Neutrino


def propagate(n: Neutrino, dt: float) -> Neutrino:
    """Advance position by velocity*dt. Velocity unchanged (spec §5)."""
    return Neutrino(
        position=n.position + n.velocity * dt,
        velocity=n.velocity.copy(),
    )


def detect_overlap(a: Neutrino, b: Neutrino, r_overlap: float) -> bool:
    """Return True if two neutrinos are within r_overlap of each other.

    This is the trigger for the displacement rule (spec §5).
    """
    distance = float(np.linalg.norm(a.position - b.position))
    return distance < r_overlap


def displace(
    a: Neutrino, b: Neutrino, push: float
) -> tuple[Neutrino, Neutrino]:
    """Push two overlapping neutrinos apart along the line connecting them.

    Velocities are NOT changed (spec §5). If positions coincide exactly,
    use the velocity-difference vector as a fallback; if that's also zero,
    fall back to (1, 0).
    """
    diff = b.position - a.position
    norm = float(np.linalg.norm(diff))

    if norm < 1e-12:
        # Coincident: use velocity difference, then a fixed fallback.
        diff = b.velocity - a.velocity
        norm = float(np.linalg.norm(diff))
    if norm < 1e-12:
        diff = np.array([1.0, 0.0])
        norm = 1.0

    unit = diff / norm
    shift = unit * (push / 2.0)

    moved_a = Neutrino(position=a.position - shift, velocity=a.velocity.copy())
    moved_b = Neutrino(position=b.position + shift, velocity=b.velocity.copy())
    return moved_a, moved_b


def step(
    neutrinos: list[Neutrino],
    dt: float,
    r_overlap: float,
    push: float,
) -> list[Neutrino]:
    """One simulation step: propagate all, then resolve pairwise overlaps.

    Overlap resolution iterates pairwise (O(n²) — fine for small n; spec v1
    only requires 2-particle experiments).
    """
    moved = [propagate(n, dt) for n in neutrinos]

    # Resolve pairwise overlaps. Single pass is enough for n=2; for larger
    # n, repeat until no overlaps remain (left for v2).
    for i in range(len(moved)):
        for j in range(i + 1, len(moved)):
            if detect_overlap(moved[i], moved[j], r_overlap):
                moved[i], moved[j] = displace(moved[i], moved[j], push)

    return moved
