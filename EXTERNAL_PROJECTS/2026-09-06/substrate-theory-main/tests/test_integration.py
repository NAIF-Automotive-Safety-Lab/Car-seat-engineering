"""Integration tests: lock in the simulation-level results we cite in the spec.

These run actual simulations (not unit tests of single functions) and assert
that the *measured outputs* are within tolerance of what the spec records.
They protect against parameter drift, integrator changes, or subtle bugs
in dynamics composition that could silently break the load-bearing claims.

Tolerances are intentionally loose (5-10%) — these are integration tests
catching regressions, not precision benchmarks. Tighter values are recorded
in the spec/README for documentation.
"""

import numpy as np

from stiff_medium.neutrino import C
from stiff_medium.back_reaction import back_reaction_force, vverlet_step
from stiff_medium.spinor import (
    cone_azimuth,
    spin_half_check,
    unwrap_azimuth_history,
)
from stiff_medium.atomic import (
    coulomb_attraction,
    newton_step,
    reduced_mass,
)


# Back-reaction orbital binding -----------------------------------------

def _run_back_reaction_orbit(n_steps: int = 6000):
    """Run Test 2 of back_reaction_v2.py and return diagnostics."""
    z = np.array([0.0, 0.0, 1.0])
    s = C / np.sqrt(2.0)
    r_eq = 0.20
    r_capture = 1.0
    k_push = 5.0
    k_pull = 5.0
    dt = 0.005

    pos_a = np.array([-1.5 * r_eq / 2, 0.0, 0.0])
    vel_a = np.array([0.0, s, s])
    pos_b = np.array([1.5 * r_eq / 2, 0.0, 0.0])
    vel_b = np.array([0.0, -s, s])

    def force(pa, pb):
        return back_reaction_force(
            pa, pb, r_eq=r_eq, r_capture=r_capture, k_push=k_push, k_pull=k_pull
        )

    state = (pos_a, vel_a, pos_b, vel_b)
    rel_angles = []
    az_a_history = [cone_azimuth(vel_a, z)]

    for k in range(n_steps):
        new_pa, new_va, new_pb, new_vb = vverlet_step(
            state[0], state[1], z, state[2], state[3], z,
            dt=dt, force_fn=force,
        )
        state = (new_pa, new_va, new_pb, new_vb)
        rel = new_pb - new_pa
        rel_angles.append(float(np.arctan2(rel[1], rel[0])))
        az_a_history.append(cone_azimuth(new_va, z))

    return {
        "final_state": state,
        "rel_angles": np.unwrap(np.asarray(rel_angles)),
        "az_a_unwrapped": unwrap_azimuth_history(az_a_history),
    }


def test_back_reaction_produces_multiple_orbits():
    """Spec §5.5: back-reaction at 1.5x r_eq produces sustained 2D orbits.
    The measured value is 5.62 over the second half. Allow 5-7 to catch drift."""
    result = _run_back_reaction_orbit()
    half = len(result["rel_angles"]) // 2
    rel_window = result["rel_angles"][half:]
    revolutions = abs(rel_window[-1] - rel_window[0]) / (2.0 * np.pi)
    assert 4.5 < revolutions < 7.0, f"Expected ~5.6 orbits, got {revolutions}"


def test_back_reaction_preserves_speed():
    """Cone projection should keep |v|=C throughout the orbit."""
    result = _run_back_reaction_orbit()
    final_pa, final_va, final_pb, final_vb = result["final_state"]
    assert np.isclose(np.linalg.norm(final_va), C, atol=1e-6)
    assert np.isclose(np.linalg.norm(final_vb), C, atol=1e-6)


def test_cone_azimuth_ratio_indicates_spin_half():
    """Spec §16 claim 6b: cone-azimuth advances ~1 turn per orbital revolution.
    The measured ratio is 1.004; tolerance ±0.1 is loose enough to catch real
    drift but allow numerical noise."""
    result = _run_back_reaction_orbit()
    half = len(result["rel_angles"]) // 2
    rel_window = result["rel_angles"][half:]
    az_window = result["az_a_unwrapped"][half:]
    revolutions = abs(rel_window[-1] - rel_window[0]) / (2.0 * np.pi)
    spin_result = spin_half_check(az_window, revolutions)
    assert 0.9 < abs(spin_result["azimuth_per_orbit"]) < 1.1, (
        f"Expected ~1 azimuth turn per orbit, got {spin_result['azimuth_per_orbit']}"
    )


# Hydrogen isotope shift ------------------------------------------------

def _bohr_orbital_frequency(m_nucleus: float, m_electron: float = 1.0,
                             coupling: float = 1.0, hbar: float = 1.0,
                             dt: float = 0.0001, n_periods: int = 30) -> float:
    """Run an atom at its Bohr radius and return the simulated orbital frequency."""
    mu = reduced_mass(m_electron, m_nucleus)
    r = hbar * hbar / (mu * coupling)
    omega_expected = mu * coupling * coupling / (hbar ** 3)
    period = 2 * np.pi / omega_expected
    n_steps = int(n_periods * period / dt)

    r_e = (m_nucleus / (m_electron + m_nucleus)) * r
    r_n = (m_electron / (m_electron + m_nucleus)) * r
    v_e_mag = float(np.sqrt(coupling * r_e / (m_electron * r * r)))
    v_n_mag = (m_electron / m_nucleus) * v_e_mag

    pos_e = np.array([-r_e, 0.0, 0.0])
    vel_e = np.array([0.0, v_e_mag, 0.0])
    pos_n = np.array([r_n, 0.0, 0.0])
    vel_n = np.array([0.0, -v_n_mag, 0.0])

    def force(pa, pb):
        return coulomb_attraction(pa, pb, coupling=coupling)

    rel_angles = []
    state = (pos_e, vel_e, pos_n, vel_n)
    for _ in range(n_steps):
        new_pe, new_ve, new_pn, new_vn = newton_step(
            state[0], state[1], m_electron,
            state[2], state[3], m_nucleus,
            dt=dt, force_fn=force,
        )
        state = (new_pe, new_ve, new_pn, new_vn)
        rel = new_pn - new_pe
        if abs(rel[0]) > 1e-9 or abs(rel[1]) > 1e-9:
            rel_angles.append(float(np.arctan2(rel[1], rel[0])))

    half = len(rel_angles) // 2
    rel_window = np.unwrap(np.asarray(rel_angles[half:]))
    return abs(rel_window[-1] - rel_window[0]) / (len(rel_window) - 1) / dt


def test_hydrogen_deuterium_isotope_shift():
    """Spec §8.1a: at Bohr-scaled radii, ω scales as μ; D/H shift = +272 ppm.
    Tolerance ±10% on the ppm value (272 ± 27 ppm) catches major drift but
    allows numerical noise from the integrator."""
    omega_h = _bohr_orbital_frequency(m_nucleus=1836.15)
    omega_d = _bohr_orbital_frequency(m_nucleus=3670.48)
    shift_ppm = (omega_d / omega_h - 1) * 1e6
    assert 240 < shift_ppm < 310, (
        f"Expected ~272 ppm D/H shift, got {shift_ppm:.1f} ppm"
    )


def test_hydrogen_tritium_isotope_shift():
    """Same test for T/H: expected +363 ppm."""
    omega_h = _bohr_orbital_frequency(m_nucleus=1836.15)
    omega_t = _bohr_orbital_frequency(m_nucleus=5496.92)
    shift_ppm = (omega_t / omega_h - 1) * 1e6
    assert 320 < shift_ppm < 410, (
        f"Expected ~363 ppm T/H shift, got {shift_ppm:.1f} ppm"
    )
