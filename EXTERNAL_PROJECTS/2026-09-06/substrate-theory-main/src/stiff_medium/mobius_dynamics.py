"""Möbius internal twist as a dynamical degree of freedom (spec §13 gap #1).

Each neutrino carries an internal phase ψ that advances with the cone
azimuth. The slope sign is determined by ψ. Under Möbius topology
(half-integer winding), ψ has period 2π and the slope sign flips
when ψ crosses π.

This module promotes Möbius from interpretation overlay (passive observable
in spinor.py) to a dynamical state variable: each step of the simulation
updates ψ explicitly, and the back-reaction force depends on the relative
slope signs of the interacting pair.

Test of fermion behavior: two same-Möbius particles should fail to form
stable bound states (Pauli-like exclusion of identical fermions).
Two opposite-Möbius particles should bind normally.
"""

from dataclasses import dataclass
import numpy as np

from stiff_medium.neutrino import C
from stiff_medium.spinor import cone_azimuth


@dataclass
class MobiusState:
    """A neutrino's full state including Möbius internal phase.

    Position, velocity, axis are the physical quantities. accumulated_azimuth
    tracks the cumulative cone-azimuth winding (NOT modulo 2π). The slope
    sign is computed from accumulated_azimuth via Möbius topology.

    Not frozen: this state must be updated each simulation step. The
    underlying neutrino's identity is preserved via initial_axis.
    """
    position: np.ndarray
    velocity: np.ndarray
    axis: np.ndarray
    accumulated_azimuth: float  # cumulative, NOT mod 2π
    initial_phase: float  # ψ at t=0; usually 0 or π for the two species

    @property
    def mobius_phase(self) -> float:
        """ψ = (initial_phase + accumulated_azimuth / 2) mod 2π."""
        psi = (self.initial_phase + self.accumulated_azimuth / 2.0) % (2.0 * np.pi)
        return float(psi)

    @property
    def slope_sign(self) -> int:
        """+1 for ψ in [0, π); −1 for ψ in [π, 2π). Flips per orbital revolution."""
        return 1 if self.mobius_phase < np.pi else -1


def update_accumulated_azimuth(
    state: MobiusState,
    new_velocity: np.ndarray,
) -> float:
    """Compute the new accumulated azimuth given a new velocity.

    Takes the change in cone azimuth from the old velocity to the new,
    correctly handling wrap-around (the change should be the small angle,
    not the long way around).
    """
    old_az = cone_azimuth(state.velocity, state.axis)
    new_az = cone_azimuth(new_velocity, state.axis)
    delta = new_az - old_az
    # Wrap to [-π, π] so we don't accumulate spurious 2π jumps from arctan2 wrap
    if delta > np.pi:
        delta -= 2.0 * np.pi
    elif delta < -np.pi:
        delta += 2.0 * np.pi
    return state.accumulated_azimuth + delta


def mobius_aware_force(
    a: MobiusState,
    b: MobiusState,
    base_force_fn,
    same_sign_repulsion_factor: float = -1.0,
) -> np.ndarray:
    """Force on A from B, modulated by the relative slope signs.

    - opposite slope signs (e.g., one electron-mode + one positron-mode):
      use the base force unchanged. Standard back-reaction binds them.
    - same slope signs (two electron-modes, or two positron-modes):
      multiply the attractive part of the force by `same_sign_repulsion_factor`.
      Default −1 inverts: what would be pull becomes push, so same-sign
      particles never have an attractive zone and cannot bind.
      The repulsive (close-range) part is unchanged so they still can't coincide.
    """
    base = base_force_fn(a.position, b.position)
    if a.slope_sign == b.slope_sign:
        diff = b.position - a.position
        d = float(np.linalg.norm(diff))
        if d > 1e-12:
            unit = diff / d
            radial = float(np.dot(base, unit))
            if radial > 0:  # attractive component (toward B)
                base = base + (same_sign_repulsion_factor - 1.0) * radial * unit
    return base


def mobius_vverlet_step(
    a: MobiusState,
    b: MobiusState,
    dt: float,
    force_fn,
    project_to_cone_fn,
) -> tuple[MobiusState, MobiusState]:
    """One velocity-Verlet step with Möbius phase tracking and cone projection.

    - Computes Möbius-aware force using the slope signs at this step.
    - Updates positions and velocities.
    - Projects velocities back to the 45° cone.
    - Updates accumulated_azimuth based on the change in cone azimuth.
    """
    f_a = force_fn(a, b)
    f_b = -f_a

    half_va = project_to_cone_fn(a.velocity + 0.5 * f_a * dt, a.axis)
    half_vb = project_to_cone_fn(b.velocity + 0.5 * f_b * dt, b.axis)

    new_pa = a.position + half_va * dt
    new_pb = b.position + half_vb * dt

    # Build temporary states for second-half force evaluation; phase still old
    tmp_a = MobiusState(
        position=new_pa,
        velocity=half_va,
        axis=a.axis,
        accumulated_azimuth=update_accumulated_azimuth(a, half_va),
        initial_phase=a.initial_phase,
    )
    tmp_b = MobiusState(
        position=new_pb,
        velocity=half_vb,
        axis=b.axis,
        accumulated_azimuth=update_accumulated_azimuth(b, half_vb),
        initial_phase=b.initial_phase,
    )

    new_f_a = force_fn(tmp_a, tmp_b)
    new_f_b = -new_f_a

    new_va = project_to_cone_fn(half_va + 0.5 * new_f_a * dt, a.axis)
    new_vb = project_to_cone_fn(half_vb + 0.5 * new_f_b * dt, b.axis)

    final_a = MobiusState(
        position=new_pa,
        velocity=new_va,
        axis=a.axis,
        accumulated_azimuth=update_accumulated_azimuth(a, new_va),
        initial_phase=a.initial_phase,
    )
    final_b = MobiusState(
        position=new_pb,
        velocity=new_vb,
        axis=b.axis,
        accumulated_azimuth=update_accumulated_azimuth(b, new_vb),
        initial_phase=b.initial_phase,
    )
    return final_a, final_b
