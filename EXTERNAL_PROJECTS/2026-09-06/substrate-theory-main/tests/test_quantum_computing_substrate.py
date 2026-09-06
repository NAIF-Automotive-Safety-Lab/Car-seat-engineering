"""Tests for stiff_medium.quantum_computing_substrate."""

from __future__ import annotations

import math

import numpy as np
import pytest

from stiff_medium.quantum_computing_substrate import (
    GATE_H,
    GATE_S,
    GATE_T,
    GATE_X,
    GATE_Y,
    GATE_Z,
    QuantumComputingGeometry,
    QuantumComputingSimulator,
    bell_state_vector,
    fidelity_vs_gamma,
    grover_run_4qubits,
)


# ---------------------------------------------------------------------------
# Geometry — Möbius two-sheet picture
# ---------------------------------------------------------------------------


def test_geometry_dim_and_basis() -> None:
    geom = QuantumComputingGeometry(n_qubits=3)
    assert geom.dim() == 8
    assert geom.basis_index([0, 0, 0]) == 0
    assert geom.basis_index([0, 0, 1]) == 1
    assert geom.basis_index([1, 0, 1]) == 5
    psi = geom.basis_state([1, 0, 1])
    assert psi.shape == (8,)
    assert psi[5] == pytest.approx(1.0 + 0.0j)
    assert np.allclose(np.abs(psi).sum(), 1.0)


def test_geometry_amplitudes_to_bloch_basis_states() -> None:
    geom = QuantumComputingGeometry(n_qubits=1)
    # |0> -> +z
    assert geom.amplitudes_to_bloch(1.0, 0.0) == pytest.approx((0.0, 0.0, 1.0))
    # |1> -> -z
    assert geom.amplitudes_to_bloch(0.0, 1.0) == pytest.approx((0.0, 0.0, -1.0))
    # |+> -> +x
    inv2 = 1.0 / math.sqrt(2.0)
    rx, ry, rz = geom.amplitudes_to_bloch(inv2, inv2)
    assert (rx, ry, rz) == pytest.approx((1.0, 0.0, 0.0), abs=1e-12)
    # |+i> -> +y
    rx, ry, rz = geom.amplitudes_to_bloch(inv2, inv2 * 1j)
    assert (rx, ry, rz) == pytest.approx((0.0, 1.0, 0.0), abs=1e-12)


def test_geometry_bloch_roundtrip() -> None:
    geom = QuantumComputingGeometry(n_qubits=1)
    inv2 = 1.0 / math.sqrt(2.0)
    for alpha, beta in [
        (1.0 + 0j, 0.0 + 0j),
        (0.0 + 0j, 1.0 + 0j),
        (inv2 + 0j, inv2 + 0j),
        (inv2 + 0j, inv2 * 1j),
    ]:
        rx, ry, rz = geom.amplitudes_to_bloch(alpha, beta)
        a, b = geom.bloch_to_amplitudes(rx, ry, rz)
        rx2, ry2, rz2 = geom.amplitudes_to_bloch(a, b)
        assert (rx2, ry2, rz2) == pytest.approx((rx, ry, rz), abs=1e-10)


def test_geometry_two_sheet_phase() -> None:
    geom = QuantumComputingGeometry(n_qubits=1)
    # f=0 -> 1
    assert geom.two_sheet_phase(0.0) == pytest.approx(1.0 + 0.0j)
    # f=1/4 -> i  (S phase)
    assert geom.two_sheet_phase(0.25) == pytest.approx(0.0 + 1.0j, abs=1e-12)
    # f=1/2 -> -1  (Z phase)
    assert geom.two_sheet_phase(0.5) == pytest.approx(-1.0 + 0.0j, abs=1e-12)
    # f=1   -> 1   (Möbius full traversal = identity)
    assert geom.two_sheet_phase(1.0) == pytest.approx(1.0 + 0.0j, abs=1e-12)


# ---------------------------------------------------------------------------
# Single-qubit gate unitarity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,gate",
    [
        ("X", GATE_X),
        ("Y", GATE_Y),
        ("Z", GATE_Z),
        ("H", GATE_H),
        ("S", GATE_S),
        ("T", GATE_T),
    ],
)
def test_single_qubit_gates_are_unitary(
    name: str, gate: np.ndarray,
) -> None:
    prod = gate @ gate.conjugate().T
    assert np.allclose(prod, np.eye(2), atol=1e-12), f"{name} not unitary"


def test_two_qubit_gates_are_unitary() -> None:
    cnot = QuantumComputingSimulator.cnot_matrix()
    cz = QuantumComputingSimulator.cz_matrix()
    swap = QuantumComputingSimulator.swap_matrix()
    for name, g in [("CNOT", cnot), ("CZ", cz), ("SWAP", swap)]:
        prod = g @ g.conjugate().T
        assert np.allclose(prod, np.eye(4), atol=1e-12), f"{name} not unitary"


# ---------------------------------------------------------------------------
# Single-qubit gate behaviour
# ---------------------------------------------------------------------------


def test_x_flips_qubit() -> None:
    sim = QuantumComputingSimulator(n_qubits=1)
    assert sim.state[0] == pytest.approx(1.0)
    sim.x(0)
    assert sim.state[0] == pytest.approx(0.0, abs=1e-15)
    assert sim.state[1] == pytest.approx(1.0)


def test_z_phase_on_one() -> None:
    sim = QuantumComputingSimulator(n_qubits=1)
    sim.x(0)  # |1>
    sim.z(0)
    assert sim.state[1] == pytest.approx(-1.0)


def test_h_creates_superposition() -> None:
    sim = QuantumComputingSimulator(n_qubits=1)
    sim.h(0)
    inv2 = 1.0 / math.sqrt(2.0)
    assert np.allclose(sim.state, [inv2, inv2], atol=1e-12)


def test_h_h_is_identity() -> None:
    sim = QuantumComputingSimulator(n_qubits=1)
    sim.x(0)
    sim.h(0)
    sim.h(0)
    assert sim.state[1] == pytest.approx(1.0, abs=1e-12)
    assert sim.state[0] == pytest.approx(0.0, abs=1e-12)


def test_t_squared_is_s() -> None:
    """T^2 = S."""
    t_sq = GATE_T @ GATE_T
    assert np.allclose(t_sq, GATE_S, atol=1e-12)


def test_s_squared_is_z() -> None:
    """S^2 = Z."""
    s_sq = GATE_S @ GATE_S
    assert np.allclose(s_sq, GATE_Z, atol=1e-12)


def test_state_norm_preserved_after_gates() -> None:
    sim = QuantumComputingSimulator(n_qubits=3)
    sim.h(0); sim.h(1); sim.h(2)
    sim.x(1); sim.z(2); sim.s(0); sim.t(2)
    norm = float(np.sqrt(np.real(np.vdot(sim.state, sim.state))))
    assert norm == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Two-qubit gate behaviour
# ---------------------------------------------------------------------------


def test_cnot_classical_truth_table() -> None:
    sim = QuantumComputingSimulator(n_qubits=2)
    # |00> -> |00>
    sim.reset(); sim.cnot(0, 1)
    assert sim.state[0] == pytest.approx(1.0)
    # |01> -> |01>  (control=0)
    sim.reset(); sim.x(1); sim.cnot(0, 1)
    assert sim.state[1] == pytest.approx(1.0)
    # |10> -> |11>  (control=1)
    sim.reset(); sim.x(0); sim.cnot(0, 1)
    assert sim.state[3] == pytest.approx(1.0)
    # |11> -> |10>
    sim.reset(); sim.x(0); sim.x(1); sim.cnot(0, 1)
    assert sim.state[2] == pytest.approx(1.0)


def test_swap_swaps_basis_states() -> None:
    sim = QuantumComputingSimulator(n_qubits=2)
    sim.x(0)  # |10>
    sim.swap(0, 1)
    assert sim.state[1] == pytest.approx(1.0)  # |01>


def test_cz_phase_on_11() -> None:
    sim = QuantumComputingSimulator(n_qubits=2)
    sim.x(0); sim.x(1)  # |11>
    sim.cz(0, 1)
    assert sim.state[3] == pytest.approx(-1.0)


# ---------------------------------------------------------------------------
# Bell-state preparation + entanglement check
# ---------------------------------------------------------------------------


def test_bell_phi_plus_matches_closed_form() -> None:
    sim = QuantumComputingSimulator(n_qubits=2)
    sim.prepare_bell("phi_plus")
    assert np.allclose(sim.state, bell_state_vector("phi_plus"), atol=1e-12)


@pytest.mark.parametrize(
    "kind", ["phi_plus", "phi_minus", "psi_plus", "psi_minus"],
)
def test_bell_states_are_normalised(kind: str) -> None:
    sim = QuantumComputingSimulator(n_qubits=2)
    sim.prepare_bell(kind)
    norm = float(np.real(np.vdot(sim.state, sim.state)))
    assert norm == pytest.approx(1.0, abs=1e-12)


def test_bell_phi_plus_is_entangled_via_schmidt() -> None:
    """For |phi+> = (|00>+|11>)/sqrt2, the reduced rho_A has rank 2.

    Equivalent: the Schmidt coefficients are both 1/sqrt(2), so
    the entanglement entropy is log(2), the maximum for 2 qubits.
    """
    sim = QuantumComputingSimulator(n_qubits=2)
    sim.prepare_bell("phi_plus")
    psi = sim.state.reshape(2, 2)
    # Reduced density matrix on qubit A: rho_A = Tr_B (|psi><psi|)
    rho_a = psi @ psi.conjugate().T
    eigs = np.linalg.eigvalsh(rho_a)
    eigs = np.clip(np.real(eigs), 1e-15, None)
    s = -float(np.sum(eigs * np.log(eigs)))
    assert s == pytest.approx(math.log(2.0), abs=1e-10)
    # Reduced rho_A should be (1/2) I
    assert np.allclose(rho_a, 0.5 * np.eye(2), atol=1e-12)


def test_bell_state_correlation_classical_limit_violated() -> None:
    """All four Bell states give |<XX>| or |<ZZ>| = 1: perfect correlation."""
    sim = QuantumComputingSimulator(n_qubits=2)
    sim.prepare_bell("phi_plus")
    # <ZZ> = sum_i (-1)^(z_a + z_b) |c_i|^2
    probs = sim.probabilities()
    # |00>=0, |01>=1, |10>=2, |11>=3; ZZ eigenvalue: +1, -1, -1, +1
    zz = probs[0] - probs[1] - probs[2] + probs[3]
    assert abs(zz) == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Grover's algorithm — amplitude amplification
# ---------------------------------------------------------------------------


def test_grover_optimal_iterations_4_qubits_is_3() -> None:
    """For N=16, optimal = floor(pi/4 * 4) = floor(pi) = 3."""
    assert QuantumComputingSimulator.grover_optimal_iterations(4) == 3


def test_grover_amplifies_marked_state() -> None:
    """After the optimal number of iterations, marked-state probability > 0.9."""
    sim = QuantumComputingSimulator(n_qubits=4)
    out = sim.run_grover(marked_index=5)
    final_prob = float(out["marked_probability"][-1])
    assert final_prob > 0.9, f"final marked prob {final_prob:.3f} too low"


def test_grover_initial_amplitude_is_uniform() -> None:
    """Before any iteration the marked amplitude is 1/sqrt(N) = 1/4."""
    out = grover_run_4qubits(marked_index=5)
    assert float(out["marked_amplitude"][0]) == pytest.approx(0.25, abs=1e-12)


def test_grover_simulation_matches_theory() -> None:
    """Closed-form |sin((2k+1) theta)| matches the simulated trajectory.

    The simulator stores |<m|psi_k>| (always non-negative); the theoretical
    sin((2k+1) theta) becomes negative past the optimum.  Match magnitudes.
    """
    out = grover_run_4qubits(marked_index=5)
    sim_amp = out["marked_amplitude"]
    th_amp = out["theory_amplitude"]
    assert sim_amp.shape == th_amp.shape
    assert np.allclose(sim_amp, np.abs(th_amp), atol=1e-10)


def test_grover_amplitude_oscillates_above_optimum() -> None:
    """If we keep going past k=3, amplitude eventually decreases."""
    out = grover_run_4qubits(marked_index=5)
    amps = out["marked_amplitude"]
    # k=3 is the optimal; amps[4] should be smaller than amps[3].
    assert float(amps[4]) < float(amps[3])


def test_grover_state_is_normalised_throughout() -> None:
    sim = QuantumComputingSimulator(n_qubits=4)
    sim.reset()
    for q in range(4):
        sim.h(q)
    oracle = sim.grover_oracle(7)
    diff = sim.grover_diffusion()
    for _ in range(5):
        sim.state = oracle @ sim.state
        sim.state = diff @ sim.state
        norm = float(np.real(np.vdot(sim.state, sim.state)))
        assert norm == pytest.approx(1.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Decoherence — substrate drag γ
# ---------------------------------------------------------------------------


def test_decoherence_zero_gamma_preserves_purity() -> None:
    sim = QuantumComputingSimulator(n_qubits=1, gamma=0.0)
    sim.h(0)
    rho = sim.decohere(t=10.0)
    purity = float(np.real(np.trace(rho @ rho)))
    assert purity == pytest.approx(1.0, abs=1e-12)


def test_decoherence_kills_off_diagonals() -> None:
    sim = QuantumComputingSimulator(n_qubits=1, gamma=1.0)
    sim.h(0)  # rho_01 = 1/2 initially
    rho_init = sim.density_matrix()
    assert abs(rho_init[0, 1]) == pytest.approx(0.5, abs=1e-12)
    rho_late = sim.decohere(t=10.0)
    # exp(-1*10) = 4.54e-5
    assert abs(rho_late[0, 1]) < 1e-3
    # populations preserved
    assert np.real(rho_late[0, 0]) == pytest.approx(0.5, abs=1e-12)
    assert np.real(rho_late[1, 1]) == pytest.approx(0.5, abs=1e-12)


def test_fidelity_is_exp_minus_gamma_t() -> None:
    out = fidelity_vs_gamma(t=1.0)
    assert np.allclose(out["fidelity"], np.exp(-out["gammas"] * 1.0))
    # boundary cases
    assert QuantumComputingSimulator.gate_fidelity(0.0, 1.0) == 1.0
    assert QuantumComputingSimulator.gate_fidelity(1.0, 0.0) == 1.0
    assert QuantumComputingSimulator.gate_fidelity(1.0, 1.0) == pytest.approx(
        math.exp(-1.0), abs=1e-12)


# ---------------------------------------------------------------------------
# Sanity: simulator summary
# ---------------------------------------------------------------------------


def test_simulator_summary_reports_unitary_gates() -> None:
    sim = QuantumComputingSimulator(n_qubits=2)
    info = sim.summary()
    for key in ("x_unitary", "y_unitary", "z_unitary", "h_unitary",
                "s_unitary", "t_unitary",
                "cnot_unitary", "cz_unitary", "swap_unitary"):
        assert info[key] is True, f"{key} reports non-unitary"
