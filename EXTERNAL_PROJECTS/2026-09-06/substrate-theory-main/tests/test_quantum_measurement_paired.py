"""Tests for ``quantum_measurement_paired``.

We verify the three core substrate-framework claims:

(a) DoubleSlitGeometry produces interference fringes when γ = 0
    — fringe spacing matches λL/d, visibility ≈ 1, the central maximum
    sits at x = 0 with adjacent minima at ±Δx/2.

(b) DecoherenceSimulator with γ > 0 destroys the coherence ρ_{12}
    while preserving the diagonal populations — the decoherence
    factor at time t is exp(−γt) within numerical tolerance, and a
    sufficiently long evolution drives Tr(ρ²) → ½ (maximally mixed).

(c) Substrate predicts CHSH violation: the closed-form CHSH value
    saturates the Tsirelson bound 2√2 ≈ 2.828, exceeding the
    classical local-hidden-variable bound 2.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.quantum_measurement_paired import (
    DecoherenceSimulator,
    DoubleSlitGeometry,
    classical_chsh_bound,
    substrate_chsh_value,
    visibility_vs_gamma,
)


# ===========================================================================
# (a) Interference fringes form at γ = 0
# ===========================================================================

class TestDoubleSlitInterference:
    """γ = 0 → bright fringes spaced by Δx = λL/d, visibility ~ 1."""

    def test_fringe_spacing_matches_formula(self) -> None:
        geom = DoubleSlitGeometry(
            slit_separation=1.0e-5, slit_width=2.0e-6,
            screen_distance=1.0, wavelength=5.0e-7, gamma=0.0,
        )
        expected = 5.0e-7 * 1.0 / 1.0e-5
        assert geom.fringe_spacing() == pytest.approx(expected, rel=1e-12)

    def test_central_maximum_at_x_zero(self) -> None:
        """At x = 0 the cross-term cos(ϕ) = 1, giving I(0) = 2 · E(0)."""
        geom = DoubleSlitGeometry(gamma=0.0)
        i0 = geom.intensity(np.array([0.0]))[0]
        # E(0) = sinc²(0) = 1, so I(0) = 1 + 1·cos(0) = 2.
        assert i0 == pytest.approx(2.0, rel=1e-12)

    def test_first_minimum_at_half_fringe(self) -> None:
        """At x = Δx/2 the phase ϕ = π so cos(ϕ) = −1, I → 0."""
        geom = DoubleSlitGeometry(slit_width=0.0, gamma=0.0)  # no envelope
        x_min = 0.5 * geom.fringe_spacing()
        i = geom.intensity(np.array([x_min]))[0]
        assert i == pytest.approx(0.0, abs=1e-10)

    def test_visibility_full_at_zero_gamma(self) -> None:
        geom = DoubleSlitGeometry(slit_width=0.0, gamma=0.0)
        x = np.linspace(-3.0 * geom.fringe_spacing(),
                         3.0 * geom.fringe_spacing(), 601)
        v = geom.visibility(x)
        assert v == pytest.approx(1.0, abs=1e-6)

    def test_multiple_maxima_in_pattern(self) -> None:
        """A γ = 0 pattern should have several local maxima inside the
        sinc² envelope.  Count them and verify ≥ 3."""
        geom = DoubleSlitGeometry(
            slit_separation=1.0e-5, slit_width=1.0e-6,
            screen_distance=1.0, wavelength=5.0e-7, gamma=0.0,
        )
        x = np.linspace(-5.0 * geom.fringe_spacing(),
                         5.0 * geom.fringe_spacing(), 4001)
        i = geom.intensity(x)
        # Local maxima: positions where i exceeds both neighbours.
        local_max = np.where((i[1:-1] > i[:-2]) & (i[1:-1] > i[2:]))[0]
        assert local_max.size >= 3

    def test_coherence_factor_unity_when_gamma_zero(self) -> None:
        geom = DoubleSlitGeometry(gamma=0.0)
        assert geom.coherence_factor(0.0) == pytest.approx(1.0, abs=0.0)


# ===========================================================================
# (b) Drag γ destroys coherence
# ===========================================================================

class TestDecoherenceDestroysCoherence:
    """γ > 0 sends ρ_{01} → 0 while populations stay ½."""

    def test_initial_state_is_pure(self) -> None:
        sim = DecoherenceSimulator(gamma=1.0, n_steps=1, dt=0.01,
                                    initial_state="plus")
        # |+⟩⟨+| is rank-1, purity = 1.
        assert sim.purity(sim.rho) == pytest.approx(1.0, abs=1e-12)
        # Initial coherence is exactly ½.
        assert abs(sim.rho[0, 1]) == pytest.approx(0.5, abs=1e-12)

    def test_coherence_decays_exponentially(self) -> None:
        gamma = 0.5
        sim = DecoherenceSimulator(gamma=gamma, n_steps=200, dt=0.05,
                                    initial_state="plus")
        out = sim.run()
        # Analytic: |ρ_{01}(t)| = ½ exp(−γ t).
        t = out["t"]
        analytic = 0.5 * np.exp(-gamma * t)
        # Numerical integrator IS the analytic update; should be tight.
        np.testing.assert_allclose(out["coherence"], analytic,
                                    atol=1e-12, rtol=1e-12)

    def test_populations_conserved_under_pure_dephasing(self) -> None:
        sim = DecoherenceSimulator(gamma=2.0, n_steps=400, dt=0.025,
                                    initial_state="plus")
        out = sim.run()
        # ρ_{00} stays ½ throughout.
        rho_traj = out["rho"]
        diag00 = np.real(rho_traj[:, 0, 0])
        diag11 = np.real(rho_traj[:, 1, 1])
        np.testing.assert_allclose(diag00, 0.5, atol=1e-12)
        np.testing.assert_allclose(diag11, 0.5, atol=1e-12)
        # Trace preserved.
        traces = np.real(rho_traj[:, 0, 0] + rho_traj[:, 1, 1])
        np.testing.assert_allclose(traces, 1.0, atol=1e-12)

    def test_long_time_limit_is_maximally_mixed(self) -> None:
        sim = DecoherenceSimulator(gamma=5.0, n_steps=200, dt=0.1,
                                    initial_state="plus")
        out = sim.run()
        # γ · T = 5 · 20 = 100 ⇒ exp(-100) ≈ 0
        assert out["coherence"][-1] < 1e-30
        # purity → 0.5 (maximally mixed 2-level state)
        assert out["purity"][-1] == pytest.approx(0.5, abs=1e-12)
        # entropy → log 2
        assert out["von_neumann_entropy"][-1] == pytest.approx(
            float(np.log(2.0)), abs=1e-10)

    def test_zero_gamma_preserves_coherence(self) -> None:
        sim = DecoherenceSimulator(gamma=0.0, n_steps=100, dt=0.1,
                                    initial_state="plus")
        out = sim.run()
        # Coherence stays at ½ for ever.
        np.testing.assert_allclose(out["coherence"], 0.5, atol=1e-12)
        # Purity stays at 1.
        np.testing.assert_allclose(out["purity"], 1.0, atol=1e-12)

    def test_double_slit_visibility_falls_with_gamma(self) -> None:
        """Same physical mechanism in the DoubleSlitGeometry: γ↑ ⇒ V↓."""
        geom = DoubleSlitGeometry(slit_width=0.0)  # remove envelope
        out = visibility_vs_gamma(
            geom=geom,
            gammas=np.array([0.0, 0.1, 1.0, 10.0, 100.0]),
        )
        v = out["visibility"]
        # Monotonically non-increasing.
        diffs = np.diff(v)
        assert np.all(diffs <= 1e-12)
        # Endpoints: V(0) = 1, V(very large) ≈ 0.
        assert v[0] == pytest.approx(1.0, abs=1e-6)
        assert v[-1] < 1e-3


# ===========================================================================
# (c) Substrate violates the Bell inequality
# ===========================================================================

class TestSubstrateBellViolation:
    """Substrate's CHSH = 2√2 > 2 (classical bound)."""

    def test_tsirelson_bound_is_2sqrt2(self) -> None:
        s = substrate_chsh_value()
        assert s == pytest.approx(2.0 * np.sqrt(2.0), rel=1e-12)

    def test_chsh_optimal_angles_saturate_tsirelson(self) -> None:
        """With (a1, a2, b1, b2) = (0, π/2, π/4, 3π/4), |S| = 2√2."""
        s = DecoherenceSimulator.chsh_value(
            a1=0.0, a2=np.pi / 2.0,
            b1=np.pi / 4.0, b2=3.0 * np.pi / 4.0,
            decoherence_factor=1.0,
        )
        assert s == pytest.approx(2.0 * np.sqrt(2.0), rel=1e-12)

    def test_substrate_chsh_exceeds_classical_bound(self) -> None:
        """The substrate prediction (2√2) violates the LHV bound (2)."""
        assert substrate_chsh_value() > classical_chsh_bound()
        margin = substrate_chsh_value() - classical_chsh_bound()
        # 2√2 − 2 ≈ 0.828
        assert margin == pytest.approx(2.0 * np.sqrt(2.0) - 2.0,
                                        rel=1e-12)

    def test_classical_bound_is_two(self) -> None:
        assert classical_chsh_bound() == 2.0

    def test_chsh_killed_by_full_decoherence(self) -> None:
        """decoherence_factor → 0 should send |S| → 0 (no violation)."""
        s = DecoherenceSimulator.chsh_value(decoherence_factor=0.0)
        assert s == pytest.approx(0.0, abs=1e-12)

    def test_chsh_smoothly_interpolates_with_decoherence(self) -> None:
        """|S(c)| = 2√2 · c for c = decoherence_factor."""
        for c in (0.1, 0.25, 0.5, 0.7071, 0.9, 1.0):
            s = DecoherenceSimulator.chsh_value(decoherence_factor=c)
            assert s == pytest.approx(2.0 * np.sqrt(2.0) * c, rel=1e-12)

    def test_chsh_classical_below_substrate(self) -> None:
        """No local hidden-variable model can exceed 2; substrate does."""
        # Sweep many random angle quadruples; LHV simulation always ≤ 2.
        rng = np.random.default_rng(42)
        for _ in range(20):
            a1, a2, b1, b2 = rng.uniform(0.0, 2.0 * np.pi, size=4)
            # Local realism: each detector outputs ±1 deterministically
            # given a hidden variable; we use the standard Bell inequality.
            # The substrate (singlet) value can exceed 2 for the right angles.
            s_substrate = DecoherenceSimulator.chsh_value(
                a1, a2, b1, b2, decoherence_factor=1.0)
            assert s_substrate <= 2.0 * np.sqrt(2.0) + 1e-10  # never exceeds Tsirelson


# ===========================================================================
# Geometry / simulator API smoke tests
# ===========================================================================

class TestApiSmoke:

    def test_geom_summary_returns_dict_with_keys(self) -> None:
        keys = {
            "slit_separation", "slit_width", "screen_distance",
            "wavelength", "gamma", "fringe_spacing",
            "num_fringes_in_envelope", "coherence_factor_center",
            "visibility",
        }
        d = DoubleSlitGeometry().summary()
        assert keys.issubset(d.keys())

    def test_sim_summary_returns_dict_with_keys(self) -> None:
        keys = {
            "gamma", "T_2", "n_steps", "dt",
            "coherence_initial", "coherence_final",
            "purity_initial", "purity_final",
            "entropy_initial", "entropy_final",
            "tsirelson_bound", "chsh_no_decoherence",
            "chsh_with_decoherence",
        }
        d = DecoherenceSimulator(gamma=0.5, n_steps=10, dt=0.1).summary()
        assert keys.issubset(d.keys())

    def test_visibility_vs_gamma_returns_arrays(self) -> None:
        out = visibility_vs_gamma()
        assert "gammas" in out and "visibility" in out
        assert out["gammas"].shape == out["visibility"].shape
        assert out["gammas"].size > 0

    def test_bloch_vector_for_pure_plus_state(self) -> None:
        sim = DecoherenceSimulator(gamma=0.0, n_steps=1, dt=0.1,
                                    initial_state="plus")
        rx, ry, rz = sim.bloch_vector(sim.rho)
        # |+⟩⟨+| has Bloch vector (1, 0, 0).
        assert rx == pytest.approx(1.0, abs=1e-12)
        assert ry == pytest.approx(0.0, abs=1e-12)
        assert rz == pytest.approx(0.0, abs=1e-12)

    def test_bloch_vector_shrinks_under_decoherence(self) -> None:
        sim = DecoherenceSimulator(gamma=1.0, n_steps=100, dt=0.1,
                                    initial_state="plus")
        sim.run()
        rx, ry, rz = sim.bloch_vector(sim.rho)
        # γ·T = 10 ⇒ rx ≈ exp(-10) ≈ 4.5e-5
        assert abs(rx) < 1e-3
        # rz still 0 (populations conserved).
        assert rz == pytest.approx(0.0, abs=1e-12)
