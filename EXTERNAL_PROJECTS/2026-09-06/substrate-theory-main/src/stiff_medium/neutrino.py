"""Neutrino: 1D-vector particle in the 2D simulation. See spec §5."""

from dataclasses import dataclass
import numpy as np

C: float = 1.0  # natural units; c is the medium's wave speed


@dataclass(frozen=True)
class Neutrino:
    """A neutrino is a position + velocity vector at 45° to an axis, speed c.

    The velocity vector is never reoriented during simulation (spec §5).
    Only the position is displaced when overlap conflicts arise.
    """

    position: np.ndarray  # shape (2,), float
    velocity: np.ndarray  # shape (2,), float, |vx|=|vy|, magnitude C

    def __post_init__(self) -> None:
        if self.position.shape != (2,):
            raise ValueError(f"position must be shape (2,), got {self.position.shape}")
        if self.velocity.shape != (2,):
            raise ValueError(f"velocity must be shape (2,), got {self.velocity.shape}")
        if not np.isclose(abs(self.velocity[0]), abs(self.velocity[1])):
            raise ValueError(
                f"velocity must be 45° to an axis (|vx|=|vy|), got {self.velocity}"
            )
        speed = float(np.linalg.norm(self.velocity))
        if not np.isclose(speed, C):
            raise ValueError(f"velocity magnitude must be C={C}, got {speed}")
        # Lock arrays read-only so callers can't mutate position/velocity in
        # place — pairs with frozen=True (which only blocks attribute reassign).
        self.position.setflags(write=False)
        self.velocity.setflags(write=False)
