"""Quantum computing — paired GEOMETRY + SIM + VIZ.

Substrate ontology of "quantum computation"
-------------------------------------------
In the B3 substrate framework a *qubit* is not a primitive object but a
distinguished sector of the substrate strain field — specifically a
two-sheet Möbius excitation whose two stable strain configurations on
the local Möbius bundle are identified with |0> and |1>.  Because the
Möbius surface is non-orientable, traversing it once flips the local
sheet labelling: this is the topological origin of the relative phase
that supports superposition

        |psi> = alpha |0> + beta |1>,   |alpha|^2 + |beta|^2 = 1.

Quantum *gates* are unitary rotations on this two-sheet bundle: the
Pauli operators X, Y, Z permute / sign-flip the sheet labels, the
Hadamard H is a 90 degree rotation that maps a sheet eigenstate to an
equal-amplitude superposition (the very Möbius half-twist), and S, T
are diagonal phase gates implementing partial Möbius traversals.

Two-qubit gates are *entanglers*: they identify pairs of Möbius
bundles by gauge-coupling their substrate strain fields, so the joint
state can no longer be written as a product.  CNOT, CZ, SWAP all
arise from two-bundle exchange operators on the substrate field;
their unitarity is the same gauge-symmetry that protects substrate
Lagrangian invariance.

Decoherence in this picture is *not* a separate axiom.  It is exactly
the same substrate drag gamma that appears in
L = 1/2 rho (d_t u)^2 - 1/2 K |grad u|^2 - V(u) - gamma u (d_t u);
when gamma > 0 the off-diagonal coherence rho_{ij}(t) decays as
exp(-gamma t), populations stay put, and gates progressively lose
fidelity F(t) = exp(-gamma t).  T_2 = 1/gamma is the dephasing time.

WHAT THIS MODULE DOES
---------------------
1. ``QuantumComputingGeometry`` — qubits as Möbius two-sheet states.
   Holds a register of n qubits, exposes the canonical computational
   basis state index, and converts between (alpha, beta) amplitudes
   and the Bloch vector (r_x, r_y, r_z).  Knows the Möbius half-twist
   period (2 pi) and the Hadamard half-twist angle (pi/2).

2. ``QuantumComputingSimulator`` — full state-vector simulator.
   Single-qubit gates: X, Y, Z, H, S, T.
   Two-qubit gates: CNOT, CZ, SWAP.
   Bell-state preparation: phi+, phi-, psi+, psi-.
   Grover search on n=4 qubits with one marked element.
   Pure-dephasing decoherence channel under substrate drag gamma.

VISUALS PRODUCED BY ``render_qc``
---------------------------------
90_quantum_gates.png      — Bloch sphere trajectories for X/Y/Z/H
                            single-qubit gates + a circuit diagram
                            of the standard Bell-state preparation.
91_quantum_algorithms.png — Grover amplitude vs iteration with
                            optimum at floor(pi/4 * sqrt(2^n)) and
                            the gate-fidelity decay vs substrate drag
                            gamma showing the F(t) = exp(-gamma t)
                            substrate prediction.

REFERENCES
----------
src/stiff_medium/quantum_measurement_paired.py  -- decoherence channel.
src/stiff_medium/bell_inequality.py             -- CHSH calculator.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Common single-qubit gate matrices
# ---------------------------------------------------------------------------

_SQRT2_INV = 1.0 / np.sqrt(2.0)

GATE_I: NDArray[np.complex128] = np.eye(2, dtype=np.complex128)
GATE_X: NDArray[np.complex128] = np.array(
    [[0, 1], [1, 0]], dtype=np.complex128)
GATE_Y: NDArray[np.complex128] = np.array(
    [[0, -1j], [1j, 0]], dtype=np.complex128)
GATE_Z: NDArray[np.complex128] = np.array(
    [[1, 0], [0, -1]], dtype=np.complex128)
GATE_H: NDArray[np.complex128] = _SQRT2_INV * np.array(
    [[1, 1], [1, -1]], dtype=np.complex128)
GATE_S: NDArray[np.complex128] = np.array(
    [[1, 0], [0, 1j]], dtype=np.complex128)
GATE_T: NDArray[np.complex128] = np.array(
    [[1, 0], [0, np.exp(1j * np.pi / 4.0)]], dtype=np.complex128)


# ---------------------------------------------------------------------------
# 1. QuantumComputingGeometry — qubits as Möbius two-sheet states
# ---------------------------------------------------------------------------

@dataclass
class QuantumComputingGeometry:
    """A register of n qubits as a substrate Möbius two-sheet bundle.

    The Möbius bundle has two sheets identified with |0> and |1>; a
    superposition |psi> = alpha |0> + beta |1> corresponds to a strain
    pattern that simultaneously occupies both sheets with weights
    |alpha|^2 and |beta|^2 and a relative Möbius phase arg(beta/alpha).
    The bundle is closed by the Möbius identification: a 2 pi rotation
    of the local sheet returns to the original state, while a pi
    rotation (Hadamard) maps |0> to (|0> + |1>)/sqrt(2).

    Attributes
    ----------
    n_qubits : int
        Number of qubits in the register.  The Hilbert space dimension
        is 2**n_qubits.
    mobius_period : float
        Period of the Möbius traversal that takes |0> back to |0>.
        Set to 2 pi (one full sheet rotation = identity).
    hadamard_angle : float
        Möbius traversal angle that implements Hadamard.  Set to pi/2
        (a quarter rotation = equal superposition).
    """

    n_qubits: int = 1
    mobius_period: float = 2.0 * np.pi
    hadamard_angle: float = np.pi / 2.0

    def __post_init__(self) -> None:
        if self.n_qubits < 1:
            raise ValueError("n_qubits must be >= 1")

    # ------------------------------------------------------------------ basic geometry
    def dim(self) -> int:
        """Hilbert-space dimension D = 2**n_qubits."""
        return 1 << self.n_qubits

    def basis_index(self, bits: Sequence[int]) -> int:
        """Convert bit-string to integer index in the computational basis.

        Order convention: bits[0] is the most-significant qubit.
        """
        if len(bits) != self.n_qubits:
            raise ValueError(
                f"bit string length {len(bits)} != n_qubits {self.n_qubits}")
        idx = 0
        for b in bits:
            if b not in (0, 1):
                raise ValueError(f"non-binary bit value: {b!r}")
            idx = (idx << 1) | int(b)
        return idx

    def basis_state(self, bits: Sequence[int]) -> NDArray[np.complex128]:
        """Computational-basis ket |bits> as a length-D vector."""
        psi = np.zeros(self.dim(), dtype=np.complex128)
        psi[self.basis_index(bits)] = 1.0 + 0.0j
        return psi

    def superposition(
        self,
        amplitudes: Sequence[complex],
    ) -> NDArray[np.complex128]:
        """Construct an arbitrary normalised state from amplitudes."""
        amp = np.asarray(amplitudes, dtype=np.complex128)
        if amp.size != self.dim():
            raise ValueError(
                f"amplitude vector size {amp.size} != dim {self.dim()}")
        norm = float(np.sqrt(np.real(np.vdot(amp, amp))))
        if norm <= 0.0:
            raise ValueError("amplitude vector is the zero vector")
        return amp / norm

    # ------------------------------------------------------------------ Bloch / Möbius
    @staticmethod
    def amplitudes_to_bloch(
        alpha: complex, beta: complex,
    ) -> Tuple[float, float, float]:
        """One-qubit (alpha, beta) -> Bloch vector (r_x, r_y, r_z).

        Standard QM convention: r_x = 2 Re(alpha* beta),
        r_y = 2 Im(alpha* beta), r_z = |alpha|^2 - |beta|^2.
        """
        a = complex(alpha)
        b = complex(beta)
        rx = 2.0 * float(np.real(np.conjugate(a) * b))
        ry = 2.0 * float(np.imag(np.conjugate(a) * b))
        rz = float(np.real(np.conjugate(a) * a) - np.real(np.conjugate(b) * b))
        return (rx, ry, rz)

    @staticmethod
    def bloch_to_amplitudes(
        rx: float, ry: float, rz: float,
    ) -> Tuple[complex, complex]:
        """Bloch vector -> (alpha, beta) up to global phase.

        Uses the standard parameterisation alpha = cos(theta/2),
        beta = exp(i phi) sin(theta/2).
        """
        rmag = float(np.sqrt(rx * rx + ry * ry + rz * rz))
        if rmag < 1e-12:
            return (1.0 + 0.0j, 0.0 + 0.0j)
        # theta in [0, pi]; phi = atan2(ry, rx).
        cos_theta = max(-1.0, min(1.0, rz / max(rmag, 1e-15)))
        theta = float(np.arccos(cos_theta))
        phi = float(np.arctan2(ry, rx))
        alpha = np.cos(theta / 2.0) + 0.0j
        beta = np.exp(1j * phi) * np.sin(theta / 2.0)
        return (alpha, beta)

    def two_sheet_phase(self, traversal_fraction: float) -> complex:
        """Möbius two-sheet phase exp(i 2 pi f) for traversal fraction f.

        f = 0    -> 1   (identity)
        f = 1/4  -> i   (S gate phase)
        f = 1/2  -> -1  (Z gate phase)
        f = 1    -> 1   (full Möbius traversal returns to identity)
        """
        return complex(np.exp(1j * traversal_fraction * self.mobius_period))

    # ------------------------------------------------------------------ summary
    def summary(self) -> Dict[str, object]:
        return {
            "n_qubits": self.n_qubits,
            "dim": self.dim(),
            "mobius_period": self.mobius_period,
            "hadamard_angle": self.hadamard_angle,
            "two_sheet_phase_at_S": self.two_sheet_phase(0.25),
            "two_sheet_phase_at_Z": self.two_sheet_phase(0.5),
        }


# ---------------------------------------------------------------------------
# 2. QuantumComputingSimulator — full state-vector simulator
# ---------------------------------------------------------------------------

@dataclass
class QuantumComputingSimulator:
    """State-vector simulator with substrate-drag decoherence.

    The state vector psi is held as a length-D = 2**n_qubits complex
    array.  Single-qubit and two-qubit gates act by tensor-product
    expansion over the register.  Decoherence is the same pure-dephasing
    channel from the substrate Lagrangian: rho_{ij}(t) = rho_{ij}(0)
    exp(-gamma t) for i != j.

    Attributes
    ----------
    n_qubits : int
        Number of qubits.  Hilbert dim = 2**n_qubits.
    gamma : float
        Substrate decoherence rate (T_2 = 1/gamma).  Used by
        ``decohere`` to compute fidelity loss.
    state : NDArray[complex]
        Current state vector, initialised to |0...0>.
    """

    n_qubits: int = 1
    gamma: float = 0.0
    state: NDArray[np.complex128] = field(init=False)

    def __post_init__(self) -> None:
        if self.n_qubits < 1:
            raise ValueError("n_qubits must be >= 1")
        self.state = np.zeros(self.dim(), dtype=np.complex128)
        self.state[0] = 1.0 + 0.0j  # |00...0>

    # ------------------------------------------------------------------ helpers
    def dim(self) -> int:
        return 1 << self.n_qubits

    def reset(self) -> None:
        """Reinitialise to |00...0>."""
        self.state[:] = 0.0
        self.state[0] = 1.0 + 0.0j

    def set_state(self, psi: NDArray[np.complex128]) -> None:
        """Overwrite the state vector (no normalisation check)."""
        if psi.shape != (self.dim(),):
            raise ValueError(
                f"state shape {psi.shape} != ({self.dim()},)")
        self.state = psi.astype(np.complex128, copy=True)

    def density_matrix(self) -> NDArray[np.complex128]:
        """rho = |psi><psi| as a (D, D) complex matrix."""
        psi = self.state
        return np.outer(psi, psi.conjugate())

    def probabilities(self) -> NDArray[np.float64]:
        """|c_i|^2 over the computational basis."""
        return np.real(self.state * self.state.conjugate()).astype(np.float64)

    # ------------------------------------------------------------------ gate application
    def _apply_single(
        self,
        gate: NDArray[np.complex128],
        qubit: int,
    ) -> None:
        """Apply a 2x2 unitary to ``qubit`` (0 = most significant).

        Builds the full Kronecker-padded operator I (x) ... (x) gate (x) ... (x) I
        and multiplies the state vector.  Trivially correct for any n.
        """
        if not (0 <= qubit < self.n_qubits):
            raise ValueError(
                f"qubit index {qubit} out of range [0, {self.n_qubits})")
        op = None
        for k in range(self.n_qubits):
            piece = gate if k == qubit else GATE_I
            op = piece if op is None else np.kron(op, piece)
        self.state = op @ self.state

    def _apply_two(
        self,
        gate_4x4: NDArray[np.complex128],
        qubit_a: int,
        qubit_b: int,
    ) -> None:
        """Apply a 4x4 unitary to qubits (a, b) (a is the upper index).

        Builds the full 2^n x 2^n matrix by Kronecker-padding identity
        on the other qubits — clean and obviously correct for n <= 8.
        """
        if qubit_a == qubit_b:
            raise ValueError("two-qubit gate needs distinct qubits")
        if not (0 <= qubit_a < self.n_qubits and 0 <= qubit_b < self.n_qubits):
            raise ValueError(
                f"qubit indices ({qubit_a}, {qubit_b}) out of range")
        # Build full matrix by iterating over basis states.
        d = self.dim()
        full = np.zeros((d, d), dtype=np.complex128)
        for i in range(d):
            bits_i = [(i >> (self.n_qubits - 1 - k)) & 1
                      for k in range(self.n_qubits)]
            ia = bits_i[qubit_a]
            ib = bits_i[qubit_b]
            row_in = (ia << 1) | ib  # 0..3
            for row_out in range(4):
                amp = gate_4x4[row_out, row_in]
                if amp == 0:
                    continue
                ja = (row_out >> 1) & 1
                jb = row_out & 1
                bits_j = list(bits_i)
                bits_j[qubit_a] = ja
                bits_j[qubit_b] = jb
                j = 0
                for b in bits_j:
                    j = (j << 1) | b
                full[j, i] += amp
        self.state = full @ self.state

    # ------------------------------------------------------------------ single-qubit gates
    def x(self, qubit: int = 0) -> None:
        """Pauli-X (NOT): |0> <-> |1>."""
        self._apply_single(GATE_X, qubit)

    def y(self, qubit: int = 0) -> None:
        """Pauli-Y."""
        self._apply_single(GATE_Y, qubit)

    def z(self, qubit: int = 0) -> None:
        """Pauli-Z (phase flip): |1> -> -|1>."""
        self._apply_single(GATE_Z, qubit)

    def h(self, qubit: int = 0) -> None:
        """Hadamard: half-Möbius rotation, |0> -> (|0>+|1>)/sqrt(2)."""
        self._apply_single(GATE_H, qubit)

    def s(self, qubit: int = 0) -> None:
        """Phase gate S: diag(1, i)."""
        self._apply_single(GATE_S, qubit)

    def t(self, qubit: int = 0) -> None:
        """T gate: diag(1, exp(i pi/4))."""
        self._apply_single(GATE_T, qubit)

    # ------------------------------------------------------------------ two-qubit gates
    @staticmethod
    def cnot_matrix() -> NDArray[np.complex128]:
        """CNOT (control, target) on two qubits, control is upper index."""
        return np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
        ], dtype=np.complex128)

    @staticmethod
    def cz_matrix() -> NDArray[np.complex128]:
        """Controlled-Z."""
        return np.diag([1.0, 1.0, 1.0, -1.0]).astype(np.complex128)

    @staticmethod
    def swap_matrix() -> NDArray[np.complex128]:
        """SWAP."""
        return np.array([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ], dtype=np.complex128)

    def cnot(self, control: int, target: int) -> None:
        """Apply CNOT(control, target)."""
        self._apply_two(self.cnot_matrix(), control, target)

    def cz(self, qubit_a: int, qubit_b: int) -> None:
        """Apply controlled-Z (symmetric in its two qubits)."""
        self._apply_two(self.cz_matrix(), qubit_a, qubit_b)

    def swap(self, qubit_a: int, qubit_b: int) -> None:
        """Apply SWAP."""
        self._apply_two(self.swap_matrix(), qubit_a, qubit_b)

    # ------------------------------------------------------------------ Bell states
    def prepare_bell(self, kind: str = "phi_plus") -> None:
        """Prepare one of the four Bell states on qubits (0, 1).

        kind in {phi_plus, phi_minus, psi_plus, psi_minus}.

        Standard recipe:
          phi_+   = H(0) . CNOT(0,1) . |00>
          phi_-   = H(0) . CNOT(0,1) . |10>
          psi_+   = H(0) . CNOT(0,1) . |01>
          psi_-   = H(0) . CNOT(0,1) . |11>
        """
        if self.n_qubits < 2:
            raise ValueError("Bell-state prep needs at least 2 qubits")
        self.reset()
        if kind == "phi_plus":
            pass
        elif kind == "phi_minus":
            self.x(0)  # -> |10>
        elif kind == "psi_plus":
            self.x(1)  # -> |01>
        elif kind == "psi_minus":
            self.x(0)
            self.x(1)  # -> |11>
        else:
            raise ValueError(f"unknown Bell state {kind!r}")
        self.h(0)
        self.cnot(0, 1)

    # ------------------------------------------------------------------ Grover's algorithm
    def grover_oracle(
        self,
        marked_index: int,
    ) -> NDArray[np.complex128]:
        """Diagonal oracle that flips the sign of the marked basis state."""
        d = self.dim()
        if not (0 <= marked_index < d):
            raise ValueError(
                f"marked_index {marked_index} out of [0, {d})")
        diag = np.ones(d, dtype=np.complex128)
        diag[marked_index] = -1.0
        return np.diag(diag)

    def grover_diffusion(self) -> NDArray[np.complex128]:
        """Inversion-about-the-mean operator 2|s><s| - I, where |s> = H^n |0>."""
        d = self.dim()
        s = np.full(d, 1.0 / np.sqrt(d), dtype=np.complex128)
        return 2.0 * np.outer(s, s.conjugate()) - np.eye(d, dtype=np.complex128)

    @staticmethod
    def grover_optimal_iterations(n_qubits: int) -> int:
        """floor(pi/4 sqrt(N)) optimal-iteration count for one marked item.

        For n=4: N=16, optimal = floor(pi/4 * 4) = floor(pi) = 3 iterations.
        """
        n = 1 << n_qubits
        return int(np.floor(np.pi / 4.0 * np.sqrt(n)))

    def run_grover(
        self,
        marked_index: int,
        n_iterations: Optional[int] = None,
    ) -> Dict[str, NDArray]:
        """Run Grover search; return amplitude trajectory of the marked element.

        Parameters
        ----------
        marked_index : int
            The basis-state index whose amplitude we want to amplify.
        n_iterations : int, optional
            Number of Grover iterations.  Defaults to the optimal
            floor(pi/4 sqrt(N)).

        Returns
        -------
        dict with:
            iterations : (k+1,) integer array 0..k
            marked_amplitude : (k+1,) |<m|psi>| trajectory
            marked_probability : (k+1,) |<m|psi>|^2 trajectory
            final_state : full state vector after all iterations
            optimal_iters : the optimal iteration count
        """
        if n_iterations is None:
            n_iterations = self.grover_optimal_iterations(self.n_qubits)
        # Initialize uniform superposition |s> = H^n |0>
        self.reset()
        for q in range(self.n_qubits):
            self.h(q)
        oracle = self.grover_oracle(marked_index)
        diff = self.grover_diffusion()
        amps: List[float] = [float(abs(self.state[marked_index]))]
        for _ in range(n_iterations):
            self.state = oracle @ self.state
            self.state = diff @ self.state
            amps.append(float(abs(self.state[marked_index])))
        amp_arr = np.asarray(amps, dtype=float)
        prob_arr = amp_arr ** 2
        iters = np.arange(amp_arr.size, dtype=int)
        return {
            "iterations": iters,
            "marked_amplitude": amp_arr,
            "marked_probability": prob_arr,
            "final_state": self.state.copy(),
            "optimal_iters": np.int64(
                self.grover_optimal_iterations(self.n_qubits)),
        }

    @staticmethod
    def grover_amplitude_theory(
        n_qubits: int, n_iters: int,
    ) -> NDArray[np.float64]:
        """Closed-form Grover marked amplitude after k iterations.

        With N = 2^n_qubits and theta = arcsin(1/sqrt(N)),
        |<m|psi_k>| = sin((2k+1) theta).
        """
        n = 1 << n_qubits
        theta = float(np.arcsin(1.0 / np.sqrt(n)))
        ks = np.arange(n_iters + 1, dtype=float)
        return np.sin((2.0 * ks + 1.0) * theta).astype(np.float64)

    # ------------------------------------------------------------------ decoherence
    def decohere(self, t: float) -> NDArray[np.complex128]:
        """Apply pure dephasing for time t to the *density matrix* and return rho.

        Population diagonals are conserved; off-diagonals decay as
        exp(-gamma t).  The substrate state vector ceases to be pure
        for gamma t > 0; this method returns the resulting rho.
        """
        rho = self.density_matrix()
        if self.gamma <= 0.0 or t <= 0.0:
            return rho
        decay = float(np.exp(-self.gamma * t))
        d = self.dim()
        for i in range(d):
            for j in range(d):
                if i != j:
                    rho[i, j] *= decay
        return rho

    @staticmethod
    def gate_fidelity(gamma: float, t: float) -> float:
        """Fidelity bound F(t) = exp(-gamma t) for one gate of duration t."""
        if gamma <= 0.0 or t <= 0.0:
            return 1.0
        return float(np.exp(-gamma * t))

    @classmethod
    def is_unitary(
        cls, gate: NDArray[np.complex128], tol: float = 1e-10,
    ) -> bool:
        """Check U U^dagger = I to tolerance ``tol``."""
        d = gate.shape[0]
        prod = gate @ gate.conjugate().T
        return bool(np.allclose(prod, np.eye(d), atol=tol))

    # ------------------------------------------------------------------ summary
    def summary(self) -> Dict[str, object]:
        return {
            "n_qubits": self.n_qubits,
            "dim": self.dim(),
            "gamma": self.gamma,
            "T_2": (1.0 / self.gamma) if self.gamma > 0 else float("inf"),
            "x_unitary": self.is_unitary(GATE_X),
            "y_unitary": self.is_unitary(GATE_Y),
            "z_unitary": self.is_unitary(GATE_Z),
            "h_unitary": self.is_unitary(GATE_H),
            "s_unitary": self.is_unitary(GATE_S),
            "t_unitary": self.is_unitary(GATE_T),
            "cnot_unitary": self.is_unitary(self.cnot_matrix()),
            "cz_unitary": self.is_unitary(self.cz_matrix()),
            "swap_unitary": self.is_unitary(self.swap_matrix()),
        }


# ---------------------------------------------------------------------------
# 3. Convenience helpers used by the renderer / tests
# ---------------------------------------------------------------------------

def bell_state_vector(kind: str = "phi_plus") -> NDArray[np.complex128]:
    """Return the closed-form Bell state as a length-4 amplitude vector.

    Order: |00>, |01>, |10>, |11>.
    """
    inv2 = 1.0 / np.sqrt(2.0)
    if kind == "phi_plus":
        return np.array([inv2, 0, 0, inv2], dtype=np.complex128)
    if kind == "phi_minus":
        return np.array([inv2, 0, 0, -inv2], dtype=np.complex128)
    if kind == "psi_plus":
        return np.array([0, inv2, inv2, 0], dtype=np.complex128)
    if kind == "psi_minus":
        return np.array([0, inv2, -inv2, 0], dtype=np.complex128)
    raise ValueError(f"unknown Bell state {kind!r}")


def fidelity_vs_gamma(
    t: float = 1.0,
    gammas: Optional[Sequence[float]] = None,
) -> Dict[str, NDArray[np.float64]]:
    """Sample the gate fidelity F(t) = exp(-gamma t) over a gamma ladder."""
    if gammas is None:
        gammas = np.logspace(-3.0, 1.5, 50)
    gam = np.asarray(gammas, dtype=float)
    fid = np.exp(-gam * t).astype(np.float64)
    return {"gammas": gam, "fidelity": fid}


def grover_run_4qubits(marked_index: int = 5) -> Dict[str, NDArray]:
    """Convenience wrapper: 4-qubit Grover search on a single marked index.

    Returns a dict with the full simulation trajectory plus the
    closed-form theoretical amplitude curve for comparison.
    """
    sim = QuantumComputingSimulator(n_qubits=4)
    out = sim.run_grover(marked_index=marked_index, n_iterations=8)
    out["theory_amplitude"] = QuantumComputingSimulator.grover_amplitude_theory(
        n_qubits=4, n_iters=8,
    )
    return out
