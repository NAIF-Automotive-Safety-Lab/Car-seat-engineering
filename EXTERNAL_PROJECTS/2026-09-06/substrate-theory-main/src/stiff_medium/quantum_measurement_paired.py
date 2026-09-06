"""Quantum measurement & decoherence — paired GEOMETRY + SIM + VIZ.

Substrate ontology of "measurement" and "decoherence"
-----------------------------------------------------
In the B3 substrate framework there is no special "observer" axiom.
Both interference and its destruction are **substrate dynamics**:

  * INTERFERENCE arises because the same coherent strain pattern is
    advected through both slits and superposes on the screen.  The
    relative phase ϕ = (k · Δx) is just the path-length difference of
    the substrate wave operator, and the screen intensity reads
    I(x) = |ψ_1(x) + ψ_2(x)|² with ψ_i(x) = √I₀ exp(i k r_i)/r_i in the
    Fraunhofer limit.

  * DECOHERENCE arises when the slit / environment couples to the
    wave with a finite drag γ > 0 — the same γ that appears in the
    substrate Lagrangian L = ½ρ(∂_t u)² − ½K|∇u|² − V(u) − γ u(∂_t u).
    A drag γ on the off-diagonal of the density matrix sends the
    coherence ρ_{12}(t) = ρ_{12}(0) · exp(−γ t).  The diagonal
    populations are conserved (the strain energy stays put), so the
    mixed-state limit is the *classical* "which-slit" detector
    pattern: I_det(x) = I_1(x) + I_2(x) (no cross-term).

  * BELL VIOLATION is **not** about hidden information being signalled;
    in B3 the singlet correlation E(a, b) = −cos(a − b) is built
    into the substrate strain field's transverse polarisation, so the
    CHSH value saturates the Tsirelson bound 2√2 ≈ 2.828 just as in
    standard QM.  The substrate is non-local because the SAME elastic
    medium is shared between both detectors.

WHAT THIS MODULE DOES
---------------------
1. ``DoubleSlitGeometry``       — Fraunhofer two-slit geometry.  Knows
   the slit separation d, slit width a, screen distance L, wavelength
   λ.  Computes the interference pattern I(x) and the fringe spacing
   Δx = λ L / d.  Optional decoherence parameter γ progressively
   destroys the cross-term, smoothly interpolating from interference
   (γ = 0) to detector pattern (γ → ∞).

2. ``DecoherenceSimulator``     — 2-level density-matrix evolution
   under substrate drag.  Lindblad-like equation:

        dρ/dt = −γ · L[ρ]      with   L[ρ]_{ij} = ρ_{ij}·(1 − δ_{ij})

   so populations stay constant and coherences decay as exp(−γ t).
   Includes the purity Tr(ρ²), the von Neumann entropy −Tr(ρ log ρ),
   and the Bell test simulator (CHSH value as function of γ).

3. ``substrate_chsh_value``     — closed-form CHSH for the singlet
   under ideal substrate strain coupling: |S| = 2√2.

VISUALS PRODUCED BY ``render_quantum_measurement``
--------------------------------------------------
42_double_slit.png      — interference pattern (γ = 0, low decoherence)
                          vs detector pattern (γ ≫ 1, high decoherence),
                          plus the fringe-visibility crossover curve.
43_decoherence.png      — density-matrix coherence trace |ρ_{12}(t)|
                          decay vs γ, plus Bloch-sphere shrink and
                          Tsirelson saturation under no decoherence.

REFERENCES
----------
src/stiff_medium/wave_particle_duality.py  -- earlier double-slit code.
src/stiff_medium/bell_inequality.py        -- CHSH calculator used here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# 1. DoubleSlitGeometry — substrate field through two slits
# ---------------------------------------------------------------------------

@dataclass
class DoubleSlitGeometry:
    """Fraunhofer two-slit geometry with optional substrate decoherence γ.

    The substrate-wave intensity at screen position x is

        I(x; γ) = I_1(x) + I_2(x)
                  + 2 · √(I_1·I_2) · cos(ϕ(x)) · exp(−γ t_path(x))

    where ϕ(x) = (2π / λ) · d · x / L is the geometric phase, and
    t_path(x) ≈ r(x) / c is the substrate-wave transit time.  As γ
    grows the cross-term dies; the pattern smoothly reduces to the
    incoherent sum I_1 + I_2 (the "which-slit" detector pattern).

    The single-slit envelope is a sinc² with width set by the slit
    width a.  In the limit a → 0 the envelope flattens (point-source
    slits) and the pattern is a pure cosine fringe.

    Attributes
    ----------
    slit_separation : float
        Centre-to-centre slit separation d.  Same units as wavelength
        (default arbitrary; the result is dimensionless after dividing
        by λ L / d).
    slit_width : float
        Single-slit width a.  Sets the sinc² envelope.
    screen_distance : float
        Slit-to-screen distance L.  Fraunhofer limit assumes L ≫ d²/λ.
    wavelength : float
        de Broglie wavelength λ of the substrate-wave packet.
    gamma : float
        Substrate decoherence rate.  γ = 0 → full interference,
        γ → ∞ → classical sum of single-slit intensities.
    wave_speed : float
        Substrate wave speed c (used to convert path lengths to
        decoherence "time exposure").  Default 1.0 in natural units.
    """

    slit_separation: float = 1.0e-5      # 10 µm
    slit_width: float = 2.0e-6           # 2 µm
    screen_distance: float = 1.0         # 1 m
    wavelength: float = 5.0e-7           # 500 nm
    gamma: float = 0.0
    wave_speed: float = 1.0

    # ------------------------------------------------------------------ basic geometry
    def fringe_spacing(self) -> float:
        """Δx = λ · L / d.  Distance between adjacent intensity maxima."""
        return float(self.wavelength * self.screen_distance / self.slit_separation)

    def num_fringes_in_envelope(self) -> float:
        """Approximate count of bright fringes inside the sinc² envelope.

        The envelope half-width (first zero of sinc(πa·x/(λL))) is
        Δ_env = λ L / a; the fringe spacing is Δx = λ L / d, so the
        ratio is just d / a.
        """
        return float(self.slit_separation / self.slit_width)

    # ------------------------------------------------------------------ field helpers
    def _phase_from_path_diff(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Geometric phase ϕ(x) = (2π/λ) · d · x / L.

        This is the Fraunhofer-limit phase-difference between the two
        slits at screen position x.
        """
        k = 2.0 * np.pi / self.wavelength
        return k * self.slit_separation * x / self.screen_distance

    def _envelope(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Single-slit sinc² envelope.  E(x) = sinc²(π a x / (λ L))."""
        if self.slit_width <= 0.0:
            return np.ones_like(x)
        u = np.pi * self.slit_width * x / (self.wavelength * self.screen_distance)
        # numpy's np.sinc is sin(πx)/(πx); we want sin(u)/u ⇒ pass u/π.
        return np.asarray(np.sinc(u / np.pi)) ** 2

    def _path_length_diff_time(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Approximate "exposure time" of the substrate wave on the
        coherence-killing channel.

        The decoherence factor exp(−γ τ(x)) needs a time-like quantity.
        We use τ(x) = (path length / wave speed) of the central ray;
        in Fraunhofer this is just L/c plus an x-dependent quadratic
        correction.  For γ-only decoherence the absolute scale is
        absorbed into γ itself, so the units do not matter — what does
        matter is that τ(x) is finite and roughly uniform.
        """
        # central-ray transit time + tiny x²/(2L) correction
        return (self.screen_distance + x * x / (2.0 * self.screen_distance)) / self.wave_speed

    # ------------------------------------------------------------------ intensity
    def intensity(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Normalised screen intensity I(x; γ).

        I(x) = E(x) · [1 + cos(ϕ(x)) · exp(−γ τ(x))]

        which interpolates smoothly between
        * γ = 0  → I(x) = E(x) · (1 + cos ϕ)   (full interference)
        * γ → ∞ → I(x) = E(x)                  (which-slit detector pattern)
        """
        x = np.asarray(x, dtype=np.float64)
        env = self._envelope(x)
        phi = self._phase_from_path_diff(x)
        if self.gamma <= 0.0:
            decoh = 1.0
        else:
            tau = self._path_length_diff_time(x)
            decoh = np.exp(-self.gamma * tau)
        return env * (1.0 + np.cos(phi) * decoh)

    def coherence_factor(self, x_eval: float = 0.0) -> float:
        """C = exp(−γ · τ(0)).  Visibility multiplier of the cross-term.

        Returns the dimensionless cross-term suppression at screen
        centre x_eval.  C = 1 means full coherence (perfect fringes);
        C = 0 means full decoherence (which-slit pattern).
        """
        if self.gamma <= 0.0:
            return 1.0
        tau0 = float(self._path_length_diff_time(np.array([x_eval])))
        return float(np.exp(-self.gamma * tau0))

    def visibility(self, x: NDArray[np.float64]) -> float:
        """Michelson fringe visibility V = (I_max − I_min) / (I_max + I_min).

        Computed over the supplied screen sample x.  V = 1 for a pure
        cosine fringe (full coherence) and V = 0 for a pure envelope
        (full decoherence).
        """
        i_arr = self.intensity(x)
        imax, imin = float(np.max(i_arr)), float(np.min(i_arr))
        if imax + imin <= 0.0:
            return 0.0
        return float((imax - imin) / (imax + imin))

    # ------------------------------------------------------------------ summary
    def summary(self) -> Dict[str, object]:
        x = np.linspace(-5.0 * self.fringe_spacing(),
                        5.0 * self.fringe_spacing(), 401)
        return {
            "slit_separation": self.slit_separation,
            "slit_width": self.slit_width,
            "screen_distance": self.screen_distance,
            "wavelength": self.wavelength,
            "gamma": self.gamma,
            "fringe_spacing": self.fringe_spacing(),
            "num_fringes_in_envelope": self.num_fringes_in_envelope(),
            "coherence_factor_center": self.coherence_factor(0.0),
            "visibility": self.visibility(x),
        }


# ---------------------------------------------------------------------------
# 2. DecoherenceSimulator — density-matrix evolution under drag γ
# ---------------------------------------------------------------------------

@dataclass
class DecoherenceSimulator:
    """Two-level density-matrix evolution under substrate drag γ.

    The model is the simplest pure-dephasing (a.k.a. T₂-only) Lindblad
    channel.  In the σ_z basis,

        dρ/dt = −γ · (ρ − Diag(ρ))

    so the diagonal populations are conserved and the off-diagonal
    coherences decay as ρ_{ij}(t) = ρ_{ij}(0) · exp(−γ t) for i ≠ j.

    This is exactly what you get when the substrate field that supports
    the off-diagonal coherence is coupled to a dissipative bath with
    drag rate γ — the same γ that appears in the substrate Lagrangian.

    Attributes
    ----------
    gamma : float
        Substrate decoherence rate.  Sets the off-diagonal decay
        timescale T_2 = 1/γ.
    n_steps : int
        Number of time-steps in the evolution.
    dt : float
        Time-step.  Stable for any dt because the equation is
        linear and integrated analytically.
    initial_state : str
        Either ``'plus'`` (pure |+⟩ ⟨+|, full coherence ρ_{01} = ½)
        or ``'mixed'`` (½·I, no coherence).
    """

    gamma: float = 0.5
    n_steps: int = 200
    dt: float = 0.05
    initial_state: str = "plus"

    # state held during evolution
    rho: NDArray[np.complex128] = field(init=False)
    history: List[NDArray[np.complex128]] = field(init=False, default_factory=list)
    times: List[float] = field(init=False, default_factory=list)

    def __post_init__(self) -> None:
        self.rho = self._initial_density_matrix(self.initial_state)
        self.history = [self.rho.copy()]
        self.times = [0.0]

    # ------------------------------------------------------------------ initial states
    @staticmethod
    def _initial_density_matrix(name: str) -> NDArray[np.complex128]:
        if name == "plus":
            psi = np.array([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0)
            return np.outer(psi, psi.conjugate())
        if name == "minus":
            psi = np.array([1.0, -1.0], dtype=np.complex128) / np.sqrt(2.0)
            return np.outer(psi, psi.conjugate())
        if name == "zero":
            psi = np.array([1.0, 0.0], dtype=np.complex128)
            return np.outer(psi, psi.conjugate())
        if name == "one":
            psi = np.array([0.0, 1.0], dtype=np.complex128)
            return np.outer(psi, psi.conjugate())
        if name == "mixed":
            return 0.5 * np.eye(2, dtype=np.complex128)
        raise ValueError(f"Unknown initial state: {name!r}")

    # ------------------------------------------------------------------ one step (analytic)
    def step(self) -> None:
        """Advance one dt under pure dephasing.

        Analytic update: diagonals unchanged, off-diagonals × exp(−γ dt).
        """
        decay = float(np.exp(-self.gamma * self.dt))
        self.rho[0, 1] *= decay
        self.rho[1, 0] *= decay
        # Snapshot
        self.times.append(self.times[-1] + self.dt)
        self.history.append(self.rho.copy())

    def run(self) -> Dict[str, NDArray]:
        """Run the full evolution and return the trajectory."""
        for _ in range(self.n_steps):
            self.step()
        return {
            "t": np.array(self.times, dtype=float),
            "rho": np.stack(self.history, axis=0),
            "coherence": np.array(
                [abs(rho[0, 1]) for rho in self.history], dtype=float),
            "purity": np.array(
                [self.purity(rho) for rho in self.history], dtype=float),
            "von_neumann_entropy": np.array(
                [self.von_neumann_entropy(rho) for rho in self.history],
                dtype=float),
        }

    # ------------------------------------------------------------------ derived quantities
    @staticmethod
    def purity(rho: NDArray[np.complex128]) -> float:
        """Tr(ρ²).  Pure state → 1; maximally mixed (2-level) → ½."""
        return float(np.real(np.trace(rho @ rho)))

    @staticmethod
    def von_neumann_entropy(rho: NDArray[np.complex128]) -> float:
        """S(ρ) = −Tr(ρ log ρ).  Zero on pure states; log 2 on ½·I."""
        eigs = np.linalg.eigvalsh(rho)
        eigs = np.clip(np.real(eigs), 1e-15, None)
        return float(-np.sum(eigs * np.log(eigs)))

    @staticmethod
    def bloch_vector(rho: NDArray[np.complex128]) -> Tuple[float, float, float]:
        """(r_x, r_y, r_z) such that ρ = ½(I + r·σ).  |r| = 1 ⇒ pure."""
        r_x = 2.0 * float(np.real(rho[0, 1]))
        r_y = -2.0 * float(np.imag(rho[0, 1]))
        r_z = float(np.real(rho[0, 0] - rho[1, 1]))
        return (r_x, r_y, r_z)

    def coherence_at(self, t: float) -> float:
        """|ρ_{01}(t)| under the analytic decay law."""
        return float(abs(self.history[0][0, 1]) * np.exp(-self.gamma * t))

    # ------------------------------------------------------------------ Bell / CHSH
    @staticmethod
    def singlet_correlation(theta_a: float, theta_b: float,
                            decoherence_factor: float = 1.0) -> float:
        """E(a,b) = −cos(a − b) · decoherence_factor.

        decoherence_factor ∈ [0, 1] is the substrate coherence factor;
        with decoherence_factor = 1 the singlet is unaffected and CHSH
        saturates Tsirelson 2√2.  decoherence_factor → 0 destroys the
        correlation and CHSH → 0.
        """
        return float(-np.cos(theta_a - theta_b) * decoherence_factor)

    @classmethod
    def chsh_value(cls,
                   a1: float = 0.0,
                   a2: float = np.pi / 2.0,
                   b1: float = np.pi / 4.0,
                   b2: float = 3.0 * np.pi / 4.0,
                   decoherence_factor: float = 1.0) -> float:
        """|S| = |E(a1,b1) − E(a1,b2) + E(a2,b1) + E(a2,b2)|.

        With the canonical optimal angles and decoherence_factor = 1
        this returns 2√2 (the Tsirelson bound), which is the substrate
        framework's prediction for a perfect singlet.
        """
        e = lambda a, b: cls.singlet_correlation(a, b, decoherence_factor)
        s = e(a1, b1) - e(a1, b2) + e(a2, b1) + e(a2, b2)
        return float(abs(s))

    @classmethod
    def substrate_chsh_value(cls) -> float:
        """Closed-form Tsirelson bound 2√2 = substrate's CHSH prediction."""
        return float(2.0 * np.sqrt(2.0))

    @staticmethod
    def classical_chsh_bound() -> float:
        """The local-hidden-variable bound: any classical theory ≤ 2."""
        return 2.0

    # ------------------------------------------------------------------ summary
    def summary(self) -> Dict[str, object]:
        run = self.run()
        return {
            "gamma": self.gamma,
            "T_2": (1.0 / self.gamma) if self.gamma > 0 else float("inf"),
            "n_steps": self.n_steps,
            "dt": self.dt,
            "coherence_initial": float(run["coherence"][0]),
            "coherence_final": float(run["coherence"][-1]),
            "purity_initial": float(run["purity"][0]),
            "purity_final": float(run["purity"][-1]),
            "entropy_initial": float(run["von_neumann_entropy"][0]),
            "entropy_final": float(run["von_neumann_entropy"][-1]),
            "tsirelson_bound": float(2.0 * np.sqrt(2.0)),
            "chsh_no_decoherence": self.chsh_value(decoherence_factor=1.0),
            "chsh_with_decoherence": self.chsh_value(
                decoherence_factor=float(run["coherence"][-1] /
                                          max(run["coherence"][0], 1e-15))),
        }


# ---------------------------------------------------------------------------
# 3. Convenience top-level helpers used by the renderer
# ---------------------------------------------------------------------------

def substrate_chsh_value() -> float:
    """Closed-form: substrate predicts the Tsirelson bound 2√2."""
    return DecoherenceSimulator.substrate_chsh_value()


def classical_chsh_bound() -> float:
    """The local-hidden-variable bound 2.0."""
    return DecoherenceSimulator.classical_chsh_bound()


def visibility_vs_gamma(geom: Optional[DoubleSlitGeometry] = None,
                        gammas: Optional[Sequence[float]] = None,
                        n_screen: int = 801) -> Dict[str, NDArray[np.float64]]:
    """Sample fringe visibility V(γ) over a γ ladder.

    Returns dict with ``gammas`` and ``visibility`` arrays, suitable
    for plotting the V(γ) crossover curve.
    """
    if geom is None:
        geom = DoubleSlitGeometry()
    if gammas is None:
        gammas = np.logspace(-3.0, 1.5, 25)
    fs = geom.fringe_spacing()
    x = np.linspace(-6.0 * fs, 6.0 * fs, n_screen)
    base = geom.gamma
    vis = []
    for g in gammas:
        geom.gamma = float(g)
        vis.append(geom.visibility(x))
    geom.gamma = base
    return {
        "gammas": np.asarray(gammas, dtype=float),
        "visibility": np.asarray(vis, dtype=float),
    }
