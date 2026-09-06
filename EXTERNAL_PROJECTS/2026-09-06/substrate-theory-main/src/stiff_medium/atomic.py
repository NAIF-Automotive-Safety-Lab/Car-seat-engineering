"""Atomic-scale dynamics: bound-state COMs, NOT cone-constrained.

The cone constraint of spec §5 applies to the underlying neutrinos.
At the atomic scale, the relevant objects are already-bound multi-neutrino
configurations (electrons, nucleons) whose center-of-mass dynamics is
free — the COM velocity is not constrained to c, even though the
internal neutrinos all move at c on their cones.

This module provides Newton-like dynamics for COM particles, suitable
for modeling atoms (electron + nucleus) and their isotope variants.

The interaction force at this scale is Coulomb-like (1/r² attractive
for opposite slope shapes; per spec §10), reflecting the medium's
back-reaction averaged over the bound states' internal structure.
"""

from typing import Callable
import numpy as np


def coulomb_attraction(
    pos_a: np.ndarray,
    pos_b: np.ndarray,
    coupling: float = 1.0,
) -> np.ndarray:
    """Attractive 1/r² central force on A from B (and minus this on B).

    F = -coupling / r² · r̂  (attractive: pulls A toward B).

    Per spec §10: at the atomic scale, slope-shape complementarity
    (electron's trough fits proton's hill) gives an attractive force
    that scales as 1/r² in the far field. The coupling parameter
    encodes the medium's effective response strength at this scale.
    """
    diff = pos_b - pos_a
    d = float(np.linalg.norm(diff))
    if d < 1e-12:
        return np.zeros(3)
    unit = diff / d
    return (coupling / (d * d)) * unit  # toward B


def newton_step(
    pos_a: np.ndarray,
    vel_a: np.ndarray,
    mass_a: float,
    pos_b: np.ndarray,
    vel_b: np.ndarray,
    mass_b: float,
    dt: float,
    force_fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One velocity-Verlet step for two COM particles. No cone constraint.

    `force_fn(pos_a, pos_b)` returns the force ON a from b. Force on b
    is the negative (Newton's third law).
    """
    f_a = force_fn(pos_a, pos_b)
    f_b = -f_a

    half_vel_a = vel_a + 0.5 * (f_a / mass_a) * dt
    half_vel_b = vel_b + 0.5 * (f_b / mass_b) * dt

    new_pos_a = pos_a + half_vel_a * dt
    new_pos_b = pos_b + half_vel_b * dt

    new_f_a = force_fn(new_pos_a, new_pos_b)
    new_f_b = -new_f_a

    new_vel_a = half_vel_a + 0.5 * (new_f_a / mass_a) * dt
    new_vel_b = half_vel_b + 0.5 * (new_f_b / mass_b) * dt

    return new_pos_a, new_vel_a, new_pos_b, new_vel_b


def reduced_mass(m_light: float, m_heavy: float) -> float:
    """μ = m_light · m_heavy / (m_light + m_heavy)."""
    return m_light * m_heavy / (m_light + m_heavy)


def n_body_force(
    positions: list[np.ndarray],
    charges: list[float],
    coupling: float = 1.0,
) -> list[np.ndarray]:
    """Compute Coulomb-like force on each particle from all others.

    Force on particle i = Σ_{j ≠ i} (charge_i × charge_j / d_ij²) × r̂_{j→i}
    For attractive potential (opposite charges) this pulls particles together;
    for like charges it pushes apart.

    Per spec §18.8: multi-particle dynamics is additive pairwise at leading order.
    """
    n = len(positions)
    forces = [np.zeros(3) for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            diff = positions[j] - positions[i]
            d = float(np.linalg.norm(diff))
            if d < 1e-12:
                continue
            unit = diff / d
            # Force on i from j: attractive if charges opposite, repulsive if same
            # F_on_i = +q_i q_j / d² × unit when opposite signs (attractive: toward j)
            # F_on_i = -q_i q_j / d² × unit when same signs (repulsive: away from j)
            # Equivalent: F_on_i = -(q_i q_j / d²) × unit always (sign self-corrects)
            # No wait: for q_i=+1, q_j=-1, q_i q_j = -1, force should be attractive (+unit toward j),
            #         so F = -q_i q_j/d² × unit = +1/d² × unit ✓
            # For q_i=-1, q_j=-1, q_i q_j = +1, force should be repulsive (-unit away from j),
            #         so F = -q_i q_j/d² × unit = -1/d² × unit ✓
            f_on_i = -coupling * charges[i] * charges[j] / (d * d) * unit
            forces[i] = forces[i] + f_on_i
            forces[j] = forces[j] - f_on_i
    return forces


def n_body_newton_step(
    positions: list[np.ndarray],
    velocities: list[np.ndarray],
    masses: list[float],
    charges: list[float],
    dt: float,
    coupling: float = 1.0,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """One velocity-Verlet step for N COM particles with pairwise Coulomb forces."""
    f = n_body_force(positions, charges, coupling)
    half_v = [v + 0.5 * (f[i] / masses[i]) * dt for i, v in enumerate(velocities)]
    new_pos = [p + half_v[i] * dt for i, p in enumerate(positions)]
    new_f = n_body_force(new_pos, charges, coupling)
    new_v = [half_v[i] + 0.5 * (new_f[i] / masses[i]) * dt for i in range(len(positions))]
    return new_pos, new_v


def n_body_force_with_pauli(
    positions: list[np.ndarray],
    charges: list[float],
    spins: list[int],
    coupling: float = 1.0,
    pauli_strength: float = 1.0,
    pauli_radius: float = 0.1,
) -> list[np.ndarray]:
    """N-body force = Coulomb + Pauli exclusion for same-charge same-spin pairs.

    Per spec §5.5.1 + §18.5: identical fermions (same charge, same Möbius/spin)
    cannot occupy the same coordinate. Implemented here as a short-range Yukawa-
    like repulsion that fires only between same-charge AND same-spin particles.

    F_pauli(i, j) = pauli_strength × (pauli_radius/d)⁴ × r̂ (repulsive)
    fires only when charges[i] == charges[j] AND spins[i] == spins[j].
    """
    n = len(positions)
    forces = [np.zeros(3) for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            diff = positions[j] - positions[i]
            d = float(np.linalg.norm(diff))
            if d < 1e-12:
                continue
            unit = diff / d

            # Coulomb component (always)
            f_coulomb_on_i = -coupling * charges[i] * charges[j] / (d * d) * unit

            # Pauli component (only same charge AND same spin)
            f_pauli_on_i = np.zeros(3)
            if charges[i] == charges[j] and spins[i] == spins[j]:
                # Steeply repulsive at d ~ pauli_radius, fades fast at d > pauli_radius
                f_pauli_on_i = -pauli_strength * (pauli_radius / d) ** 4 * unit

            f_on_i = f_coulomb_on_i + f_pauli_on_i
            forces[i] = forces[i] + f_on_i
            forces[j] = forces[j] - f_on_i
    return forces


def n_body_step_with_pauli(
    positions: list[np.ndarray],
    velocities: list[np.ndarray],
    masses: list[float],
    charges: list[float],
    spins: list[int],
    dt: float,
    coupling: float = 1.0,
    pauli_strength: float = 1.0,
    pauli_radius: float = 0.1,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Velocity-Verlet step with Pauli exclusion for identical fermions."""
    f = n_body_force_with_pauli(
        positions, charges, spins, coupling, pauli_strength, pauli_radius
    )
    half_v = [v + 0.5 * (f[i] / masses[i]) * dt for i, v in enumerate(velocities)]
    new_pos = [p + half_v[i] * dt for i, p in enumerate(positions)]
    new_f = n_body_force_with_pauli(
        new_pos, charges, spins, coupling, pauli_strength, pauli_radius
    )
    new_v = [half_v[i] + 0.5 * (new_f[i] / masses[i]) * dt for i in range(len(positions))]
    return new_pos, new_v


def em_radiation_reaction(
    positions: list[np.ndarray],
    velocities: list[np.ndarray],
    nucleus_idx: int,
    charges: list[float],
    bohr_radii: list[float],
    radiation_strength: float = 0.1,
) -> list[np.ndarray]:
    """Larmor-style radiation reaction force, opposing radial drift.

    Physical motivation (per spec §11): accelerating charges radiate
    EM. The reaction force damps deviation from a preferred orbital
    radius. The preferred radius for each electron is its Bohr-quantized
    n=1, 2, 3, ... orbit; we provide bohr_radii[i] for each particle.

    For nucleus (nucleus_idx), no reaction force.
    For electrons: force opposes the radial velocity component when
    far from the nearest Bohr radius:

        F_em = -radiation_strength × (r - r_bohr_nearest) × v_radial

    This is a phenomenological capture of EM radiation reaction; the
    full Abraham-Lorentz expression involves the third time derivative
    of position, but this simpler form damps drift effectively.
    """
    n = len(positions)
    forces = [np.zeros(3) for _ in range(n)]
    nucleus_pos = positions[nucleus_idx]
    for i in range(n):
        if i == nucleus_idx:
            continue
        r_vec = positions[i] - nucleus_pos
        r = float(np.linalg.norm(r_vec))
        if r < 1e-12:
            continue
        r_hat = r_vec / r

        # Find the nearest Bohr radius for this particle
        r_bohr = bohr_radii[i]

        # Radial velocity component (positive = outward)
        v_rel = velocities[i] - velocities[nucleus_idx]
        v_radial = float(np.dot(v_rel, r_hat))

        # Damping force opposing radial drift, scaled by deviation from r_bohr
        deviation = r - r_bohr
        # Force is along -r_hat × v_radial (damps radial motion)
        # Stronger when far from r_bohr
        damping_factor = radiation_strength * abs(deviation) * np.sign(v_radial)
        forces[i] = -damping_factor * r_hat
    return forces


def n_body_step_with_em_damping(
    positions: list[np.ndarray],
    velocities: list[np.ndarray],
    masses: list[float],
    charges: list[float],
    spins: list[int],
    bohr_radii: list[float],
    nucleus_idx: int,
    dt: float,
    coupling: float = 1.0,
    pauli_strength: float = 1.0,
    pauli_radius: float = 0.1,
    radiation_strength: float = 0.1,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Verlet step with Coulomb + Pauli + EM radiation reaction.

    Adds Larmor-style damping toward Bohr-quantized orbits per spec §11
    (energy loss from non-resonant configurations as EM radiation).
    """
    f_pauli = n_body_force_with_pauli(
        positions, charges, spins, coupling, pauli_strength, pauli_radius
    )
    f_em = em_radiation_reaction(
        positions, velocities, nucleus_idx, charges, bohr_radii, radiation_strength
    )
    f = [f_pauli[i] + f_em[i] for i in range(len(positions))]

    half_v = [v + 0.5 * (f[i] / masses[i]) * dt for i, v in enumerate(velocities)]
    new_pos = [p + half_v[i] * dt for i, p in enumerate(positions)]

    new_f_pauli = n_body_force_with_pauli(
        new_pos, charges, spins, coupling, pauli_strength, pauli_radius
    )
    new_f_em = em_radiation_reaction(
        new_pos, half_v, nucleus_idx, charges, bohr_radii, radiation_strength
    )
    new_f = [new_f_pauli[i] + new_f_em[i] for i in range(len(positions))]

    new_v = [half_v[i] + 0.5 * (new_f[i] / masses[i]) * dt for i in range(len(positions))]
    return new_pos, new_v
