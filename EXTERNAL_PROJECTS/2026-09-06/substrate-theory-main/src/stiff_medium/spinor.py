"""Spin-½ via Möbius internal twist on the 45° cone (spec §13 gap #1).

The neutrino's strain pattern is given a half-period boundary condition:
the slope's sign flips after one full azimuthal rotation around the
particle's intrinsic axis (2π → -1, 4π → +1, returning only after
720°). This is the topological signature of spin-½.

Functions:
- cone_azimuth(velocity, axis): extract the cone's azimuthal angle φ.
- mobius_phase(azimuth): compute the internal twist phase ψ = φ/2.
- slope_sign(azimuth): sign(cos(ψ)) — flips at φ = π and returns at 2π.

These are passive observables for analyzing simulation results; the
underlying dynamics (back-reaction, velocity-Verlet, cone projection)
operate on positions and velocities only. Spin-½ emerges from the
interpretation of those tracked quantities under Möbius topology.
"""

from typing import Sequence
import numpy as np


def cone_azimuth(velocity: np.ndarray, axis: np.ndarray) -> float:
    """Return the azimuthal angle φ ∈ [0, 2π) of `velocity` on the
    45° cone around `axis`.

    The cone's perpendicular component is parameterized by φ around
    the axis. We pick a reference perpendicular direction e1 and
    measure φ as the angle from e1 to the velocity's perpendicular
    component (signed by the right-hand rule about `axis`).
    """
    axis = axis / float(np.linalg.norm(axis))
    if abs(axis[0]) < 0.9:
        e1 = np.cross(axis, np.array([1.0, 0.0, 0.0]))
    else:
        e1 = np.cross(axis, np.array([0.0, 1.0, 0.0]))
    e1 = e1 / float(np.linalg.norm(e1))
    e2 = np.cross(axis, e1)

    v_perp = velocity - float(np.dot(velocity, axis)) * axis
    perp_norm = float(np.linalg.norm(v_perp))
    if perp_norm < 1e-12:
        return 0.0
    cos_phi = float(np.dot(v_perp, e1)) / perp_norm
    sin_phi = float(np.dot(v_perp, e2)) / perp_norm
    phi = float(np.arctan2(sin_phi, cos_phi))
    if phi < 0:
        phi += 2.0 * np.pi
    return phi


def mobius_phase(azimuth: float) -> float:
    """Internal twist phase ψ = azimuth / 2.

    Period of ψ is 2π (returns at 4π of azimuth → 720° rotation).
    This is the Möbius half-integer winding.
    """
    return azimuth / 2.0


def slope_sign(azimuth: float) -> int:
    """Sign of the neutrino's slope under Möbius topology.

    Returns +1 for ψ in [0, π) (i.e., azimuth in [0, 2π)),
    and -1 for ψ in [π, 2π) (i.e., azimuth in [2π, 4π)).

    Note: this requires tracking the *unwrapped* (cumulative) azimuth,
    not just the value mod 2π. Use unwrap_azimuth_history first.
    """
    psi = mobius_phase(azimuth)
    return 1 if int(psi // np.pi) % 2 == 0 else -1


def mobius_state_returned(azimuth_unwrapped: float, atol: float = 1e-6) -> bool:
    """True when the Möbius state has fully returned to its original sign
    (after a multiple of 720° = 4π in the cumulative cone azimuth).

    Returns True if azimuth_unwrapped is approximately a multiple of 4π.
    Use this to detect when a system has completed a fermionic 720° cycle.
    """
    if abs(azimuth_unwrapped) < atol:
        return True
    rem = abs(azimuth_unwrapped) % (4.0 * np.pi)
    return rem < atol or (4.0 * np.pi - rem) < atol


def mobius_state_flipped(azimuth_unwrapped: float, atol: float = 1e-6) -> bool:
    """True when the Möbius state is fully sign-flipped from its original
    (after an *odd* multiple of 360° = 2π in the cumulative cone azimuth).

    A spin-½ system at this point has its 'wavefunction' equal to MINUS
    the original — half-way through the 720° return cycle.
    """
    if abs(azimuth_unwrapped) < atol:
        return False
    cycles = azimuth_unwrapped / (2.0 * np.pi)
    nearest_int = round(cycles)
    if abs(cycles - nearest_int) > atol / (2.0 * np.pi):
        return False
    return nearest_int % 2 != 0


def unwrap_azimuth_history(azimuths: Sequence[float]) -> np.ndarray:
    """Convert a sequence of azimuths in [0, 2π) into cumulative
    (unwrapped) azimuths so that the total winding is preserved.

    Uses np.unwrap, which detects 2π discontinuities and removes them.
    """
    return np.unwrap(np.asarray(azimuths))


def spin_half_check(
    azimuths_unwrapped: np.ndarray,
    orbital_revolutions: float,
) -> dict:
    """Check whether the cone azimuth rotation matches a spin-½ pattern
    over a known number of orbital revolutions.

    Returns a dict with:
    - total_azimuth_rotation: cumulative Δφ in radians over the run
    - azimuth_rotations: same / (2π)
    - azimuth_per_orbit: azimuth_rotations / orbital_revolutions
    - implied_spin: based on azimuth_per_orbit, the integer spin
      (1 = bosonic, 1/2 = fermionic, etc.)

    Spin-½ requires azimuth_per_orbit = 1 (one full azimuth rotation
    per orbital revolution; under Möbius topology, this gives a slope
    sign flip per orbit and a 720° return).
    """
    total = float(azimuths_unwrapped[-1] - azimuths_unwrapped[0])
    rotations = total / (2.0 * np.pi)
    if abs(orbital_revolutions) < 1e-9:
        return {
            "total_azimuth_rotation": total,
            "azimuth_rotations": rotations,
            "azimuth_per_orbit": float("nan"),
            "implied_spin": "undefined (no orbit)",
        }
    per_orbit = rotations / orbital_revolutions
    if abs(abs(per_orbit) - 1.0) < 0.2:
        spin = "1/2 (Möbius half-integer winding)"
    elif abs(abs(per_orbit) - 2.0) < 0.2:
        spin = "1 (bosonic, integer winding)"
    elif abs(per_orbit) < 0.2:
        spin = "0 (no internal rotation)"
    else:
        spin = f"non-standard ({per_orbit:.2f} per orbit)"
    return {
        "total_azimuth_rotation": total,
        "azimuth_rotations": rotations,
        "azimuth_per_orbit": per_orbit,
        "implied_spin": spin,
    }
