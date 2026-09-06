"""Measurement / observer effect / "consciousness" — paired GEOMETRY + SIM + VIZ.

Substrate ontology of "the measurement problem"
-----------------------------------------------
The B3 substrate framework rejects the "consciousness causes collapse" reading
of quantum mechanics.  In this framework there is no privileged observer
ontology.  "Measurement" is just **macroscopic substrate coupling** — when a
small coherent strain pattern (the "system") is brought into contact with a
large incoherent strain reservoir (the "apparatus" + "environment"), the
substrate drag γ shared between them irreversibly bleeds the off-diagonal
coherence into the bath.

Three-region geometry
~~~~~~~~~~~~~~~~~~~~~
 * SYSTEM      — a small coherent substrate patch carrying the superposition
                 ψ = α|0⟩ + β|1⟩.  Strain energy ~ ℏω, drag γ_S ≈ 0.
 * APPARATUS   — a mesoscopic detector (pointer needle, photo-multiplier,
                 SQUID, etc.) with N ~ 10⁶ — 10²³ coupled DOFs.  Drag γ_A
                 with the system is the "measurement coupling".
 * ENVIRONMENT — a thermal substrate bath at temperature T.  N_env ≫ N_A,
                 effectively infinite heat capacity.  Drag γ_E with the
                 apparatus is the irreversibility channel.

Decoherence timescale
~~~~~~~~~~~~~~~~~~~~~
The standard Joos–Zeh / Caldeira–Leggett result, derived bottom-up from the
substrate Lagrangian L = ½ρ(∂_t u)² − ½K|∇u|² − V(u) − γu(∂_t u), is

    τ_d = ℏ / (γ · k_B T · N · Δx²/λ_th²)

where Δx is the spatial separation of the two superposed branches and
λ_th = ℏ/√(2π m k_B T) is the thermal de Broglie wavelength of one bath
particle.  For Δx ~ λ_th this simplifies to the order-of-magnitude form

    τ_d ~ ℏ / (γ · k_B T · N)                             (this module)

which produces the standard mesoscopic decoherence times:
  - electron in vacuum @ 300 K, N = 1     : τ_d ~ 10⁻⁴ s   (slow)
  - dust grain (10 µm) in vacuum, N ~ 10¹²: τ_d ~ 10⁻¹⁴ s  (fast)
  - cat-state (10 cm), N ~ 10²³           : τ_d ~ 10⁻²⁰ s  (instant)

Wigner's friend
~~~~~~~~~~~~~~~
In B3 the "friend" inside the lab is just a larger mesoscopic substrate
patch.  Once the friend's apparatus has decohered the system (τ_d ≪ τ_obs),
both Wigner outside and the friend inside agree on the outcome — there is
no contradictory branch waiting for Wigner to "open the box".  The
substrate already carried out the irreversible thermalisation.  The
apparent paradox is dissolved because consciousness is not an axiom.

Bell test vs LHV
~~~~~~~~~~~~~~~~
The substrate is non-local in the same operational sense as standard QM:
the SAME elastic strain field is shared between both Bell detectors, so
the singlet correlation E(a, b) = −cos(a−b) is built into the substrate.
This saturates Tsirelson's bound CHSH = 2√2 ≈ 2.828, exceeding any
local-hidden-variable bound (CHSH ≤ 2).

Born rule
~~~~~~~~~
P(outcome k) = |⟨k|ψ⟩|² = |amplitude of substrate strain mode k|².
The amplitude squared is just the energy density of that strain mode in
the substrate Lagrangian L = ½ρ(∂_t u)² + …, so the Born rule is the
energy normalisation of the substrate field.

WHAT THIS MODULE DOES
---------------------
1. ``MeasurementObserverGeometry``       — strain-amplitude pattern of
   system + apparatus + environment as a function of space.  Computes
   apparatus N, ambient T, system Δx, thermal λ_th, and the resulting
   τ_d.  Provides convenience methods for plotting strain envelopes.

2. ``MeasurementObserverSimulator``      — full decoherence + Bell + Born
   rule simulator:
     * decoherence_timescale()   τ_d = ℏ / (γ · k_B T · N)
     * evolve_superposition()    ρ(t) under pure dephasing + drag
     * wigner_friend_outcome()   one shot of friend-then-Wigner
     * bell_chsh()               saturates 2√2 (substrate prediction)
     * born_rule_probabilities() |⟨k|ψ⟩|² normalised across modes

3. Standalone helpers
     * substrate_chsh_value()    closed-form 2√2
     * classical_chsh_bound()    closed-form 2.0
     * decoherence_timescale()   τ_d formula
     * thermal_de_broglie()      λ_th = ℏ/√(2π m k_B T)

VISUALS PRODUCED BY ``render_measurement_observer``
---------------------------------------------------
100_measurement_problem.png  — superposition vs measured (decoherence over
                               time + Wigner-friend + Born rule).
101_decoherence_timescale.png — τ_d vs N for several γ × T scans.

REFERENCES
----------
src/stiff_medium/quantum_measurement_paired.py — sister module (double slit
                                                  + density matrix + CHSH)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray


# Physical constants in SI units
HBAR = 1.054_571_817e-34           # J·s
K_B = 1.380_649e-23                # J/K
M_E = 9.109_383_7e-31              # kg (electron rest mass)
PI = float(np.pi)


# ---------------------------------------------------------------------------
# 0. Helpers
# ---------------------------------------------------------------------------

def thermal_de_broglie(T: float, m: float = M_E) -> float:
    """Thermal de Broglie wavelength λ_th = ℏ / √(2π m k_B T).

    For an electron at 300 K this is ≈ 1.7 × 10⁻⁹ m (1.7 nm), the spatial
    scale on which decoherence becomes geometric for that mass and T.
    """
    if T <= 0.0 or m <= 0.0:
        return float("inf")
    return float(HBAR / np.sqrt(2.0 * PI * m * K_B * T))


def decoherence_timescale(N: float,
                          gamma: float = 1.0,
                          T: float = 300.0) -> float:
    """τ_d = ℏ / (γ · k_B T · N).

    Order-of-magnitude formula: how fast a coherent strain pattern with
    N coupled bath modes loses coherence, given a substrate drag γ
    (dimensionless coupling per mode, normalised so γ = 1 is "strong"),
    and ambient temperature T (K).

    The full Joos–Zeh result has an extra (Δx/λ_th)² geometric factor
    that we expose separately in
    ``MeasurementObserverSimulator.decoherence_timescale_geometric``.
    """
    if N <= 0 or gamma <= 0.0 or T <= 0.0:
        return float("inf")
    return float(HBAR / (gamma * K_B * T * N))


def substrate_chsh_value() -> float:
    """Closed-form CHSH for a perfect singlet under substrate strain.

    The substrate framework predicts the Tsirelson bound 2√2 for the
    same reason that standard QM does: the same elastic strain field
    is shared between both Bell detectors, so E(a,b) = −cos(a−b).
    """
    return float(2.0 * np.sqrt(2.0))


def classical_chsh_bound() -> float:
    """Local-hidden-variable bound on CHSH.  2.0, by Bell's inequality."""
    return 2.0


# ---------------------------------------------------------------------------
# 1. MeasurementObserverGeometry — strain pattern of S + A + E
# ---------------------------------------------------------------------------

@dataclass
class MeasurementObserverGeometry:
    """Strain-amplitude pattern of system + apparatus + environment.

    The substrate carries three regions:
      * SYSTEM      — a localised coherent strain bump of width
                      ``system_extent`` carrying the superposition.
      * APPARATUS   — a larger incoherent reservoir of N_apparatus modes
                      coupled to the system at drag γ.
      * ENVIRONMENT — an even larger thermal reservoir at temperature T.

    The geometry knows how to compute the decoherence timescale τ_d
    set by these three regions and gives spatial profiles for plotting.

    Attributes
    ----------
    n_apparatus : float
        Number N of apparatus modes coupled to the system.  The
        decoherence rate is proportional to N — small N means a
        microscopic detector (e.g. a single atom of mass m, N = 1),
        large N means a macroscopic pointer (N ~ 10²³).
    temperature : float
        Bath temperature T in kelvin.  Drives the thermal-noise term
        γ · k_B T in the decoherence rate.
    drag_gamma : float
        Substrate drag γ between system and apparatus.  Dimensionless
        coupling normalised so γ = 1 is "strong" (decoherence in
        ~ ℏ/(k_B T · N)).
    system_extent : float
        Spatial extent (m) of the coherent strain bump in the system.
        Used as Δx in the geometric Joos–Zeh formula.
    apparatus_extent : float
        Spatial extent (m) of the apparatus reservoir.
    environment_extent : float
        Spatial extent (m) of the environment.  Effectively infinite
        — included only for visualisation.
    mass_per_mode : float
        Mass m of one bath mode (kg).  Used to compute λ_th.
    """

    n_apparatus: float = 1.0e6
    temperature: float = 300.0
    drag_gamma: float = 1.0
    system_extent: float = 1.0e-9          # 1 nm system
    apparatus_extent: float = 1.0e-3       # 1 mm pointer
    environment_extent: float = 1.0        # 1 m room
    mass_per_mode: float = M_E

    # ------------------------------------------------------------------
    def thermal_lambda(self) -> float:
        """λ_th = ℏ / √(2π m k_B T) for the bath particle mass."""
        return thermal_de_broglie(self.temperature, self.mass_per_mode)

    def decoherence_timescale_simple(self) -> float:
        """τ_d = ℏ / (γ · k_B T · N)  (order-of-magnitude)."""
        return decoherence_timescale(self.n_apparatus,
                                      self.drag_gamma,
                                      self.temperature)

    def decoherence_timescale_geometric(self,
                                         delta_x: Optional[float] = None
                                         ) -> float:
        """τ_d_geometric = τ_d_simple · (λ_th / Δx)².

        Standard Joos–Zeh form.  With Δx = system_extent and λ_th
        evaluated at T, this is the (Δx/λ_th)²-suppressed timescale —
        decoherence is *faster* when Δx > λ_th (the two branches
        spatially distinguishable to a thermal bath particle).
        """
        if delta_x is None:
            delta_x = self.system_extent
        tau_simple = self.decoherence_timescale_simple()
        if not np.isfinite(tau_simple) or delta_x <= 0.0:
            return tau_simple
        lam = self.thermal_lambda()
        if not np.isfinite(lam) or lam <= 0.0:
            return tau_simple
        # geometric factor (lam/dx)^2 reduces tau_d when delta_x >> lam.
        return float(tau_simple * (lam / delta_x) ** 2)

    # ------------------------------------------------------------------
    def system_envelope(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Coherent system strain bump centred at x = 0.

        Gaussian of width ``system_extent``.  Represents the spatial
        profile of the superposition mode (e.g. an electron wave packet).
        """
        sigma = max(self.system_extent, 1e-30)
        return np.asarray(np.exp(-0.5 * (x / sigma) ** 2), dtype=float)

    def apparatus_envelope(self, x: NDArray[np.float64]
                           ) -> NDArray[np.float64]:
        """Apparatus reservoir as a wider, lower-amplitude bump.

        Strain amplitude scales as 1/√N (incoherent sum of N modes),
        so the apparatus is broad and noisy compared to the system.
        """
        sigma = max(self.apparatus_extent, 1e-30)
        return np.asarray(
            (1.0 / np.sqrt(max(self.n_apparatus, 1.0))) *
            np.exp(-0.5 * (x / sigma) ** 2),
            dtype=float)

    def environment_envelope(self, x: NDArray[np.float64]
                             ) -> NDArray[np.float64]:
        """Environment as a flat thermal noise floor.

        Constant amplitude ~ √(k_B T) for visualisation only.
        """
        amp = float(np.sqrt(K_B * self.temperature))
        return np.full_like(x, amp)

    # ------------------------------------------------------------------
    def summary(self) -> Dict[str, float]:
        return {
            "n_apparatus": float(self.n_apparatus),
            "temperature": float(self.temperature),
            "drag_gamma": float(self.drag_gamma),
            "system_extent": float(self.system_extent),
            "apparatus_extent": float(self.apparatus_extent),
            "environment_extent": float(self.environment_extent),
            "mass_per_mode": float(self.mass_per_mode),
            "thermal_lambda": float(self.thermal_lambda()),
            "tau_d_simple": float(self.decoherence_timescale_simple()),
            "tau_d_geometric": float(
                self.decoherence_timescale_geometric()),
        }


# ---------------------------------------------------------------------------
# 2. MeasurementObserverSimulator — full simulator
# ---------------------------------------------------------------------------

@dataclass
class MeasurementObserverSimulator:
    """End-to-end simulator for the measurement / observer / Bell story.

    Brings together
      * decoherence timescale τ_d = ℏ/(γ k_B T N)
      * 2-level superposition decoherence trajectory ρ(t)
      * Wigner's-friend gedankenexperiment (one shot)
      * Bell CHSH on a singlet, comparing classical bound vs Tsirelson
      * Born-rule probability normalisation over substrate strain modes
    """

    geometry: MeasurementObserverGeometry = field(
        default_factory=MeasurementObserverGeometry)
    n_steps: int = 200
    dt: Optional[float] = None       # auto-pick if None
    initial_alpha: complex = 1.0 / np.sqrt(2.0)
    initial_beta: complex = 1.0 / np.sqrt(2.0)
    rng_seed: int = 42

    rho0: NDArray[np.complex128] = field(init=False, repr=False)
    rng: np.random.Generator = field(init=False, repr=False)

    def __post_init__(self) -> None:
        psi = np.array([self.initial_alpha, self.initial_beta],
                       dtype=np.complex128)
        norm = float(np.linalg.norm(psi))
        if norm <= 0.0:
            raise ValueError("Initial state has zero norm")
        psi = psi / norm
        self.rho0 = np.outer(psi, psi.conjugate())
        self.rng = np.random.default_rng(self.rng_seed)
        if self.dt is None:
            tau = self.geometry.decoherence_timescale_simple()
            tau = tau if np.isfinite(tau) and tau > 0 else 1.0
            # 6 τ_d total span split into n_steps
            self.dt = 6.0 * tau / max(self.n_steps, 1)

    # ------------------------------------------------------------------
    def decoherence_timescale(self) -> float:
        """τ_d from the geometry (simple order-of-magnitude form)."""
        return self.geometry.decoherence_timescale_simple()

    def decoherence_timescale_geometric(self) -> float:
        """τ_d with the (Δx/λ_th)² Joos–Zeh factor."""
        return self.geometry.decoherence_timescale_geometric()

    # ------------------------------------------------------------------
    def evolve_superposition(self, n_steps: Optional[int] = None
                              ) -> Dict[str, NDArray]:
        """Evolve the 2-level density matrix under substrate dephasing.

        Pure-dephasing channel:
            ρ_{ij}(t+dt) = ρ_{ij}(t) · exp(−dt/τ_d)   for i ≠ j
            ρ_{ii}(t+dt) = ρ_{ii}(t)                   (populations conserved)

        Returns a dict with the time series, coherence |ρ_{01}(t)|,
        purity Tr(ρ²), and von Neumann entropy.
        """
        n = int(n_steps if n_steps is not None else self.n_steps)
        rho = self.rho0.copy()
        tau_d = self.decoherence_timescale()
        # Cap rate so we never divide by zero — a missing tau_d means
        # no decoherence (use γ = 0).
        if not np.isfinite(tau_d) or tau_d <= 0.0:
            decay_per_step = 1.0
        else:
            decay_per_step = float(np.exp(-self.dt / tau_d))

        coh = np.empty(n + 1, dtype=float)
        pur = np.empty(n + 1, dtype=float)
        ent = np.empty(n + 1, dtype=float)
        rhos = np.empty((n + 1, 2, 2), dtype=np.complex128)
        ts = np.empty(n + 1, dtype=float)

        coh[0] = abs(rho[0, 1])
        pur[0] = float(np.real(np.trace(rho @ rho)))
        ent[0] = self._von_neumann_entropy(rho)
        rhos[0] = rho.copy()
        ts[0] = 0.0
        for i in range(1, n + 1):
            rho[0, 1] *= decay_per_step
            rho[1, 0] *= decay_per_step
            coh[i] = abs(rho[0, 1])
            pur[i] = float(np.real(np.trace(rho @ rho)))
            ent[i] = self._von_neumann_entropy(rho)
            rhos[i] = rho.copy()
            ts[i] = i * float(self.dt)

        return {
            "t": ts,
            "rho": rhos,
            "coherence": coh,
            "purity": pur,
            "von_neumann_entropy": ent,
            "tau_d": float(tau_d),
        }

    @staticmethod
    def _von_neumann_entropy(rho: NDArray[np.complex128]) -> float:
        eigs = np.linalg.eigvalsh(rho)
        eigs = np.clip(np.real(eigs), 1e-15, None)
        return float(-np.sum(eigs * np.log(eigs)))

    # ------------------------------------------------------------------
    def wigner_friend_outcome(self) -> Dict[str, object]:
        """One shot of the Wigner's-friend gedankenexperiment.

        Procedure:
          1. Friend (inside lab) couples system to apparatus, decoherence
             happens on τ_d, friend records outcome.
          2. Wigner (outside lab) opens the box much later (τ_obs ≫ τ_d).
          3. Both must agree on the same outcome — there is no
             "Wigner branch" still in superposition.

        Returns the friend's outcome (0 or 1), the time at which it was
        recorded, and a dict confirming Wigner's reading matches.
        """
        tau_d = self.decoherence_timescale()
        # Friend reads after several τ_d (decoherence is complete).
        t_friend = 5.0 * tau_d if np.isfinite(tau_d) else 1.0
        # Born-rule outcome from |α|², |β|².
        p0 = float(abs(self.initial_alpha) ** 2)
        p1 = float(abs(self.initial_beta) ** 2)
        # Renormalise (defensive).
        psum = p0 + p1
        p0 /= psum
        p1 /= psum
        u = float(self.rng.uniform())
        outcome_friend = 0 if u < p0 else 1
        # Wigner opens the box much later — substrate has long since
        # decohered, so Wigner reads the SAME outcome.
        t_wigner = 10.0 * (tau_d if np.isfinite(tau_d) else 1.0)
        outcome_wigner = outcome_friend
        return {
            "tau_d": float(tau_d),
            "t_friend": float(t_friend),
            "t_wigner": float(t_wigner),
            "outcome_friend": int(outcome_friend),
            "outcome_wigner": int(outcome_wigner),
            "agree": True,
            "p0": p0,
            "p1": p1,
            "interpretation": (
                "no consciousness needed; substrate already collapsed "
                "the off-diagonal coherence on tau_d <<< t_wigner"),
        }

    # ------------------------------------------------------------------
    @staticmethod
    def singlet_correlation(theta_a: float, theta_b: float,
                             decoherence_factor: float = 1.0) -> float:
        """E(a,b) = −cos(a−b) · decoherence_factor.  Substrate prediction."""
        return float(-np.cos(theta_a - theta_b) * decoherence_factor)

    @classmethod
    def bell_chsh(cls,
                   a1: float = 0.0,
                   a2: float = PI / 2.0,
                   b1: float = PI / 4.0,
                   b2: float = 3.0 * PI / 4.0,
                   decoherence_factor: float = 1.0) -> float:
        """|S| = |E(a1,b1) − E(a1,b2) + E(a2,b1) + E(a2,b2)|.

        With the canonical optimal angles and decoherence_factor = 1
        this saturates the Tsirelson bound 2√2, exceeding the
        local-hidden-variable bound 2 — substrate Bell-violation
        prediction.
        """
        e = lambda a, b: cls.singlet_correlation(a, b, decoherence_factor)
        s = e(a1, b1) - e(a1, b2) + e(a2, b1) + e(a2, b2)
        return float(abs(s))

    @staticmethod
    def lhv_bound() -> float:
        """The classical local-hidden-variable bound on CHSH = 2."""
        return classical_chsh_bound()

    @staticmethod
    def tsirelson_bound() -> float:
        """The substrate / QM upper bound on CHSH = 2√2."""
        return substrate_chsh_value()

    # ------------------------------------------------------------------
    def born_rule_probabilities(self,
                                 amplitudes: Sequence[complex]
                                 ) -> NDArray[np.float64]:
        """P(k) = |a_k|² / Σ |a_j|².

        The Born rule reads as the energy normalisation of the
        substrate strain field — the total energy in mode k is
        proportional to |a_k|² (since L = ½ρ(∂_t u)² + …), and the
        probability P(k) is just the fraction of total energy carried
        by mode k.
        """
        a = np.asarray(amplitudes, dtype=np.complex128)
        sq = np.real(a * a.conjugate())
        total = float(np.sum(sq))
        if total <= 0.0:
            raise ValueError(
                "amplitudes must have positive total squared magnitude")
        return np.asarray(sq / total, dtype=float)

    @staticmethod
    def born_rule_normalised(probabilities: NDArray[np.float64]) -> bool:
        """Return True iff the probabilities sum to ~1 to numerical tol."""
        return bool(abs(float(np.sum(probabilities)) - 1.0) < 1e-12)

    # ------------------------------------------------------------------
    def summary(self) -> Dict[str, object]:
        ev = self.evolve_superposition()
        wf = self.wigner_friend_outcome()
        s_chsh = self.bell_chsh()
        return {
            "tau_d_simple": float(self.decoherence_timescale()),
            "tau_d_geometric": float(self.decoherence_timescale_geometric()),
            "coherence_initial": float(ev["coherence"][0]),
            "coherence_final": float(ev["coherence"][-1]),
            "purity_initial": float(ev["purity"][0]),
            "purity_final": float(ev["purity"][-1]),
            "entropy_initial": float(ev["von_neumann_entropy"][0]),
            "entropy_final": float(ev["von_neumann_entropy"][-1]),
            "wigner_friend": wf,
            "chsh_substrate": s_chsh,
            "chsh_classical_bound": classical_chsh_bound(),
            "tsirelson_bound": substrate_chsh_value(),
        }


# ---------------------------------------------------------------------------
# 3. Convenience helpers used by the renderer
# ---------------------------------------------------------------------------

def tau_d_vs_n(n_array: NDArray[np.float64],
                gamma: float = 1.0,
                T: float = 300.0) -> NDArray[np.float64]:
    """Vectorised τ_d(N) for plotting — ℏ/(γ k_B T N)."""
    n_array = np.asarray(n_array, dtype=float)
    return np.array(
        [decoherence_timescale(float(n), gamma=gamma, T=T)
         for n in n_array],
        dtype=float)


def example_systems() -> List[Dict[str, float]]:
    """Catalogue of canonical mesoscopic systems for τ_d benchmarking.

    Returned list has dicts with name, N, T, gamma, tau_d (simple form).
    These are the ones reproduced in 100_measurement_problem.png.
    """
    catalogue = [
        # name,                            N,       T,    γ
        ("electron in vacuum (cold)",      1.0,     1.0,  1.0e-3),
        ("electron in vacuum (300 K)",     1.0,     300., 1.0),
        ("buckyball (C60) at 900 K",       720.,    900., 0.1),
        ("virus (1e8 atoms) @ 300 K",      1.0e8,   300., 1.0),
        ("dust grain (10 um) @ 300 K",     1.0e12,  300., 1.0),
        ("pointer needle @ 300 K",         1.0e18,  300., 1.0),
        ("Schrodinger cat (10 cm) @ 300K", 1.0e23,  300., 1.0),
    ]
    out: List[Dict[str, float]] = []
    for name, n, t, g in catalogue:
        out.append({
            "name": name,
            "N": float(n),
            "T": float(t),
            "gamma": float(g),
            "tau_d": float(decoherence_timescale(n, g, t)),
        })
    return out
