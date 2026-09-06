"""3D extension of the simulation. Neutrino3D, propagate, detect_overlap,
displace, step — all in 3D — plus the 45°-cone validation that the spec
requires.

The 2D modules (neutrino.py, dynamics.py) remain the v1/v2 simulation;
this module is the v3 (3D) work. The two coexist so the v1 results stay
reproducible.
"""

from dataclasses import dataclass
import numpy as np

from stiff_medium.neutrino import C  # share the natural-units C=1.0


@dataclass(frozen=True)
class Neutrino3D:
    """3D neutrino: position, velocity, and per-particle axis.

    The velocity must lie on the cone of 45° to the axis, with magnitude C.
    Per spec §5, 45° is the *uniquely stable* angle, not just permitted.
    """

    position: np.ndarray  # shape (3,), float
    velocity: np.ndarray  # shape (3,), float, on 45° cone of axis, magnitude C
    axis: np.ndarray      # shape (3,), float, unit vector

    def __post_init__(self) -> None:
        for name, arr in (
            ("position", self.position),
            ("velocity", self.velocity),
            ("axis", self.axis),
        ):
            if arr.shape != (3,):
                raise ValueError(f"{name} must be shape (3,), got {arr.shape}")

        axis_norm = float(np.linalg.norm(self.axis))
        if not np.isclose(axis_norm, 1.0):
            raise ValueError(f"axis must be unit vector, got |axis|={axis_norm}")

        speed = float(np.linalg.norm(self.velocity))
        if not np.isclose(speed, C):
            raise ValueError(f"velocity magnitude must be C={C}, got {speed}")

        cos_theta = float(np.dot(self.velocity, self.axis)) / speed
        target = 1.0 / np.sqrt(2.0)  # cos(45°)
        if not np.isclose(cos_theta, target):
            raise ValueError(
                f"velocity must be at 45° to axis "
                f"(cos θ = 1/√2 ≈ {target:.4f}), got cos θ = {cos_theta:.4f}"
            )

        self.position.setflags(write=False)
        self.velocity.setflags(write=False)
        self.axis.setflags(write=False)


def make_on_cone(
    axis: np.ndarray, azimuth: float, speed: float = C
) -> np.ndarray:
    """Construct a 3D velocity on the 45° cone around `axis`, at azimuth φ.

    The along-axis component is speed/√2 (cos 45°). The perpendicular
    component is speed/√2, rotated by `azimuth` around the axis.
    """
    axis = axis / float(np.linalg.norm(axis))
    # Build an orthonormal frame (axis, e1, e2).
    if abs(axis[0]) < 0.9:
        e1 = np.cross(axis, np.array([1.0, 0.0, 0.0]))
    else:
        e1 = np.cross(axis, np.array([0.0, 1.0, 0.0]))
    e1 = e1 / float(np.linalg.norm(e1))
    e2 = np.cross(axis, e1)

    s = speed / np.sqrt(2.0)
    along = s * axis
    perp = s * (np.cos(azimuth) * e1 + np.sin(azimuth) * e2)
    return along + perp


def propagate(n: Neutrino3D, dt: float) -> Neutrino3D:
    """Advance position by velocity*dt; velocity and axis unchanged."""
    return Neutrino3D(
        position=n.position + n.velocity * dt,
        velocity=n.velocity.copy(),
        axis=n.axis.copy(),
    )


def detect_overlap(a: Neutrino3D, b: Neutrino3D, r_overlap: float) -> bool:
    distance = float(np.linalg.norm(a.position - b.position))
    return distance < r_overlap


def displace(
    a: Neutrino3D, b: Neutrino3D, push: float
) -> tuple[Neutrino3D, Neutrino3D]:
    """Push two overlapping neutrinos apart along the connecting line.

    Velocities and axes are NOT changed (spec §5). Same fallback policy
    as the 2D version: if positions coincide, fall back to velocity
    difference, then to (1, 0, 0).
    """
    diff = b.position - a.position
    norm = float(np.linalg.norm(diff))
    if norm < 1e-12:
        diff = b.velocity - a.velocity
        norm = float(np.linalg.norm(diff))
    if norm < 1e-12:
        diff = np.array([1.0, 0.0, 0.0])
        norm = 1.0

    unit = diff / norm
    shift = unit * (push / 2.0)

    moved_a = Neutrino3D(
        position=a.position - shift,
        velocity=a.velocity.copy(),
        axis=a.axis.copy(),
    )
    moved_b = Neutrino3D(
        position=b.position + shift,
        velocity=b.velocity.copy(),
        axis=b.axis.copy(),
    )
    return moved_a, moved_b


def step(
    neutrinos: list[Neutrino3D],
    dt: float,
    r_overlap: float,
    push: float,
) -> list[Neutrino3D]:
    moved = [propagate(n, dt) for n in neutrinos]
    for i in range(len(moved)):
        for j in range(i + 1, len(moved)):
            if detect_overlap(moved[i], moved[j], r_overlap):
                moved[i], moved[j] = displace(moved[i], moved[j], push)
    return moved
