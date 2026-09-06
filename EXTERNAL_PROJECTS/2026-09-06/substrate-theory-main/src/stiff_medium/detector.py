"""Bound-state detection. Spec §6: a stable orbital pattern is two
neutrinos persisting near each other long enough that they're not just
passing through. Operationally: within r_bound for `persistence` steps."""

import numpy as np
from stiff_medium.neutrino import Neutrino


class BoundStateTracker:
    """Flags when two neutrinos have been within r_bound for N consecutive steps."""

    def __init__(self, r_bound: float, persistence: int) -> None:
        if r_bound <= 0:
            raise ValueError("r_bound must be positive")
        if persistence < 1:
            raise ValueError("persistence must be >= 1")
        self.r_bound = r_bound
        self.persistence = persistence
        self._consecutive = 0

    def update(self, neutrinos: list[Neutrino]) -> bool:
        """Update tracker with current state. Returns True if a bound
        state has been detected (persistence threshold reached this step).

        v1: only handles the 2-neutrino case; ignores n != 2.
        """
        if len(neutrinos) != 2:
            self._consecutive = 0
            return False

        dist = float(np.linalg.norm(neutrinos[0].position - neutrinos[1].position))
        if dist < self.r_bound:
            self._consecutive += 1
        else:
            self._consecutive = 0

        return self._consecutive >= self.persistence
