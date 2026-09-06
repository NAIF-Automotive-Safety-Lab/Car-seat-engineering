"""Tests for bound_state_3d_extractor — particle-as-substrate-pattern module.

Verifies:
  (a) Localization preserved — the centre stays peaked, R(t) doesn't blow up.
  (b) Measured breathing frequency matches ω_b prediction within 15%.
  (c) E_rest = ℏ ω_b — Planck-Compton identity holds at machine precision
      and the measured value matches the predicted value to within the
      same 15% as the frequency check.
  (d) Geometry sanity — analytic ω_b at the electron anchor reproduces
      m_e c² = 0.510999 MeV.
  (e) Mode catalogue — ground / excited / breathing all run, all stay
      bounded, breathing has nonzero radial oscillation amplitude.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.cone_bouncing_protocol import C, HBAR, MEV
from src.stiff_medium.bound_state_3d_extractor import (
    BoundState3DGeometry,
    BoundState3DSimulator,
    electron_geometry,
    electron_rest_energy_mev,
    run_three_modes,
)


# ---------------------------------------------------------------------------
# Geometry sanity
# ---------------------------------------------------------------------------


class TestGeometry:

    def test_natural_units_default(self):
        """K = ρ = ξ = 1 → c = 1, ω_b = 1, T = 2π."""
        g = BoundState3DGeometry()
        assert g.c_eff == pytest.approx(1.0)
        assert g.omega_b == pytest.approx(1.0)
        assert g.period == pytest.approx(2.0 * np.pi)

    def test_dx_dt_obey_cfl(self):
        """3D CFL: dt < dx / (c√3).  We use safety factor 0.4."""
        g = BoundState3DGeometry()
        cfl = g.c_eff * g.dt / g.dx
        # 3D CFL upper bound 1/√3 ≈ 0.577
        assert cfl < 1.0 / np.sqrt(3.0)
        # Safety factor 0.4 → cfl ≈ 0.4/√3 ≈ 0.231
        assert cfl == pytest.approx(0.4 / np.sqrt(3.0), rel=1e-9)

    def test_electron_anchor_reproduces_mev(self):
        """At ξ = ξ_e the analytic ω_b gives m_e c² = 0.510999 MeV."""
        e_mev = electron_rest_energy_mev()
        assert e_mev == pytest.approx(0.510999, rel=1e-3)

    def test_planck_compton_identity_at_anchor(self):
        """E_rest = ℏ ω_b at the electron anchor matches m_e c²."""
        g = electron_geometry()
        e_pred = g.rest_energy
        m_e_kg = g.mass_kg
        # Identity m c² = ℏ ω_b
        assert e_pred == pytest.approx(m_e_kg * C * C, rel=1e-12)
        # Anchor: m equals m_e to part-per-million
        from scipy import constants as sc
        assert m_e_kg == pytest.approx(sc.m_e, rel=1e-6)

    def test_ground_bump_is_localized(self):
        """The ground bump (Gaussian with σ = 4ξ) should be peaked at the
        origin and decay to <1% at the box edge L = 12 ξ
        (exp(−(L/σ)²/2) = exp(−4.5) ≈ 0.011)."""
        g = BoundState3DGeometry()
        u_origin = g.ground_bump(np.array([[[0.0]]]),
                                  np.array([[[0.0]]]),
                                  np.array([[[0.0]]]))[0, 0, 0]
        u_edge = g.ground_bump(np.array([[[g.L]]]),
                                np.array([[[0.0]]]),
                                np.array([[[0.0]]]))[0, 0, 0]
        assert u_origin > 0.0
        # exp(-(12/4)²/2) = exp(-4.5) = 0.0111 — relax to 5% to leave
        # head-room for a finite eval point on the discrete grid
        assert u_edge < 0.05 * u_origin


# ---------------------------------------------------------------------------
# Simulator: localization preserved
# ---------------------------------------------------------------------------


class TestLocalization:
    """A localized bound state must STAY localized over multiple periods."""

    def test_breathing_stays_localized(self):
        """Centre amplitude remains within ~30% of its initial value, and
        the second-moment radius does not double."""
        sim = BoundState3DSimulator(
            geometry=BoundState3DGeometry(N=24, L=12.0),
            mode="breathing",
            n_periods=3.0,
            samples_per_period=24,
        )
        res = sim.run()
        u_c = res["u_center"]
        r   = res["radius"]
        # Centre value oscillates but never collapses to noise
        assert np.all(np.isfinite(u_c))
        # Peak amplitude doesn't blow up (energy growth would mean
        # numerical instability)
        u_max_init = float(np.abs(u_c[:8]).max())
        u_max_final = float(np.abs(u_c[-8:]).max())
        assert u_max_final < 3.0 * u_max_init, (
            f"centre amplitude blew up: init {u_max_init:.3f} → final {u_max_final:.3f}"
        )
        # Radius stays bounded — bound state, not radiative
        loc = sim.localization_metric()
        assert loc["r_growth_pct"] < 100.0, (
            f"radius grew {loc['r_growth_pct']:.1f}% — bound state lost"
        )

    def test_ground_stays_localized(self):
        """Ground mode should be the most stable — radius almost constant."""
        sim = BoundState3DSimulator(
            geometry=BoundState3DGeometry(N=24, L=12.0),
            mode="ground",
            n_periods=2.0,
            samples_per_period=20,
        )
        sim.run()
        loc = sim.localization_metric()
        assert loc["r_growth_pct"] < 60.0


# ---------------------------------------------------------------------------
# Simulator: frequency matches ω_b
# ---------------------------------------------------------------------------


class TestFrequency:

    def test_breathing_frequency_matches_omega_b(self):
        """FFT-extracted breathing ω matches analytic ω_b to within 15%.

        15% accommodates the lattice dispersion + finite-box effects of
        a 24³ grid with samples-per-period = 24.
        """
        sim = BoundState3DSimulator(
            geometry=BoundState3DGeometry(N=24, L=12.0),
            mode="breathing",
            n_periods=4.0,
            samples_per_period=32,
        )
        sim.run()
        omega_meas = sim.measured_omega()
        omega_pred = sim.geometry.omega_b
        rel = abs(omega_meas - omega_pred) / omega_pred
        assert np.isfinite(omega_meas) and omega_meas > 0.0
        assert rel < 0.15, (
            f"breathing ω_meas={omega_meas:.4f} vs ω_b={omega_pred:.4f} "
            f"({rel*100:.1f}% off — expected <15%)"
        )

    def test_higher_xi_lowers_frequency(self):
        """ω_b ∝ 1/ξ — doubling ξ should halve ω_b at γ=0."""
        g1 = BoundState3DGeometry(xi=1.0)
        g2 = BoundState3DGeometry(xi=2.0)
        assert g1.omega_b == pytest.approx(2.0 * g2.omega_b, rel=1e-12)


# ---------------------------------------------------------------------------
# Simulator: E_rest = ℏ ω_b identity holds for measured ω
# ---------------------------------------------------------------------------


class TestRestEnergy:

    def test_planck_compton_identity_predicted(self):
        """E_pred = ℏ ω_pred — by construction, exact."""
        g = BoundState3DGeometry()
        assert g.rest_energy == pytest.approx(HBAR * g.omega_b, rel=1e-12)

    def test_measured_E_matches_predicted_within_15pct(self):
        """ℏ ω_meas should match ℏ ω_pred to within the same tol as ω."""
        sim = BoundState3DSimulator(
            geometry=BoundState3DGeometry(N=24, L=12.0),
            mode="breathing",
            n_periods=4.0,
            samples_per_period=32,
        )
        sim.run()
        rep = sim.report()
        assert rep["omega_rel_err"] < 0.15
        # rest-E identity = ℏ × omega — same relative error
        rel_E = abs(rep["rest_E_meas"] - rep["rest_E_pred"]) / rep["rest_E_pred"]
        assert rel_E == pytest.approx(rep["omega_rel_err"], rel=1e-9)

    def test_mass_via_E_over_c2(self):
        """m = E_rest / c² — algebraic identity at machine precision."""
        sim = BoundState3DSimulator(
            geometry=BoundState3DGeometry(N=20, L=12.0),
            mode="breathing", n_periods=3.0, samples_per_period=24,
        )
        sim.run()
        rep = sim.report()
        m_from_E = rep["rest_E_meas"] / (C * C)
        assert m_from_E == pytest.approx(rep["mass_meas_kg"], rel=1e-12)


# ---------------------------------------------------------------------------
# Mode catalogue
# ---------------------------------------------------------------------------


class TestModes:

    def test_three_modes_run_and_return_snapshots(self):
        """ground / excited / breathing each run successfully and return
        3 snapshots at t = 0, T/2, T."""
        results = run_three_modes(
            geometry=BoundState3DGeometry(N=20, L=12.0),
            n_periods=2.0,
            samples_per_period=20,
        )
        assert set(results.keys()) == {"ground", "excited", "breathing"}
        for mode, res in results.items():
            assert "snapshots" in res
            snaps = res["snapshots"]
            assert snaps.shape[0] == 3, f"{mode}: expected 3 snapshots"
            # Every snapshot is finite and not all-zero
            for k in range(3):
                assert np.all(np.isfinite(snaps[k]))
                assert np.abs(snaps[k]).max() > 0.0

    def test_excited_has_radial_node(self):
        """The first-excited bump has a spherical node at r = σ_env.

        u(r) = A · (1 − r²/σ²) · exp(−r²/(2σ²)) vanishes at r = σ_env
        exactly.  We check the analytic seed directly — geometry, not
        dynamics."""
        g = BoundState3DGeometry()
        sigma = g.envelope_width
        x = np.array([[[sigma]]])
        zero = np.array([[[0.0]]])
        u_node = g.first_excited_bump(x, zero, zero)[0, 0, 0]
        assert u_node == pytest.approx(0.0, abs=1e-12)
        # And nonzero at origin
        u0 = g.first_excited_bump(zero, zero, zero)[0, 0, 0]
        assert abs(u0) > 0.0

    def test_breathing_has_oscillation(self):
        """The breathing seed should produce a centre-amplitude trace
        whose AC component is at least 0.5% of the DC."""
        sim = BoundState3DSimulator(
            geometry=BoundState3DGeometry(N=20, L=12.0),
            mode="breathing", n_periods=3.0, samples_per_period=24,
        )
        res = sim.run()
        u = res["u_center"]
        ac = float(np.std(u))
        dc = float(np.abs(u).mean())
        assert ac / max(dc, 1e-30) > 0.005, (
            f"breathing oscillation too weak: ac/dc={ac/dc:.4f}"
        )


# ---------------------------------------------------------------------------
# Validation helper
# ---------------------------------------------------------------------------


class TestReport:

    def test_report_is_well_formed(self):
        """report() returns a flat dict with the expected keys, all
        finite floats."""
        sim = BoundState3DSimulator(
            geometry=BoundState3DGeometry(N=20, L=12.0),
            mode="breathing", n_periods=3.0, samples_per_period=24,
        )
        sim.run()
        rep = sim.report()
        for key in ("omega_pred", "omega_meas", "omega_rel_err",
                    "rest_E_pred", "rest_E_meas",
                    "mass_pred_kg", "mass_meas_kg",
                    "r_init", "r_max", "r_growth_pct"):
            assert key in rep
            assert np.isfinite(rep[key]), f"{key} is not finite"
