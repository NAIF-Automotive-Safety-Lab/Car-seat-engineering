"""Tests for substrate_3d_movie — 32³ full-evolution rendering pipeline.

Verifies the three required guarantees of a movie-renderable bound state:

  (a) Energy conserved — at γ = 0 the §18.11 Lagrangian is conservative.
      Total energy must drift by less than 5% over one breathing period
      (lattice + Verlet floor).

  (b) Field stays localized — the iso-surface sits inside the box: the
      bound-state seed does NOT radiate into the boundary, so the
      iso-volume cell count remains a tiny fraction of the full lattice.

  (c) Iso-surface volume oscillates — the breathing mode means the cell
      count V_iso(t) at the canonical iso-level changes during the
      period.  Equivalently the centre-cell amplitude has nonzero
      standard deviation.

Plus standard sanity (geometry defaults, snapshot shape, iso-point
extraction returns coordinate triples bounded by the box).
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.substrate_3d_movie import (
    Substrate3DGeometry,
    Substrate3DSimulator,
)


# ---------------------------------------------------------------------------
# Geometry sanity
# ---------------------------------------------------------------------------


class TestGeometry:

    def test_default_lattice_is_32_cubed(self):
        g = Substrate3DGeometry()
        assert g.N == 32
        assert g.L == 12.0

    def test_natural_units_default(self):
        """K = ρ = ξ = 1 → c = 1, ω_b = 1, T = 2π."""
        g = Substrate3DGeometry()
        assert g.c_eff == pytest.approx(1.0)
        assert g.omega_b == pytest.approx(1.0)
        assert g.period == pytest.approx(2.0 * np.pi)

    def test_dt_obeys_3d_cfl(self):
        """3D CFL with safety factor 0.4: cfl = c·dt/dx = 0.4/√3."""
        g = Substrate3DGeometry()
        cfl = g.c_eff * g.dt / g.dx
        assert cfl < 1.0 / np.sqrt(3.0)
        assert cfl == pytest.approx(0.4 / np.sqrt(3.0), rel=1e-9)

    def test_seed_is_localized_at_origin(self):
        """The Gaussian seed peaks at the origin and decays to <5% at
        the box edge L = 12 ξ."""
        g = Substrate3DGeometry()
        zero = np.array([[[0.0]]])
        edge = np.array([[[g.L]]])
        u_origin = g.bound_state_seed(zero, zero, zero)[0, 0, 0]
        u_edge = g.bound_state_seed(edge, zero, zero)[0, 0, 0]
        assert u_origin > 0.0
        assert abs(u_edge) < 0.05 * abs(u_origin)


# ---------------------------------------------------------------------------
# (a) Energy conserved
# ---------------------------------------------------------------------------


class TestEnergyConservation:
    """The full §18.11 Lagrangian without drag (γ = 0) must conserve E."""

    def test_total_energy_drifts_less_than_5pct(self):
        sim = Substrate3DSimulator(
            geometry=Substrate3DGeometry(N=32, L=12.0, gamma=0.0),
            n_frames=8, n_periods=1.0,
        )
        sim.run()
        diag = sim.diagnostics()
        assert diag["energy_init"] > 0.0
        assert diag["energy_drift_pct"] < 5.0, (
            f"energy drifted {diag['energy_drift_pct']:.2f}% — "
            "expected <5% with drag-free Verlet over 1 period"
        )


# ---------------------------------------------------------------------------
# (b) Field stays localized
# ---------------------------------------------------------------------------


class TestLocalization:
    """The iso-surface should remain a tiny fraction of the full lattice
    — i.e. the bound state did not radiate into the boundary."""

    def test_iso_volume_is_tiny_fraction_of_lattice(self):
        sim = Substrate3DSimulator(
            geometry=Substrate3DGeometry(N=32, L=12.0),
            n_frames=8, n_periods=1.0,
        )
        sim.run()
        diag = sim.diagnostics()
        N3 = 32 ** 3
        # Iso shell at 35% should be well under 10% of all cells
        assert diag["iso_vol_mean"] / N3 < 0.10, (
            f"iso volume {diag['iso_vol_mean']:.0f} too large vs "
            f"lattice size {N3}"
        )
        assert diag["iso_vol_mean"] > 0, "iso volume vanished"

    def test_iso_points_lie_inside_box(self):
        """All iso-surface points must have |x|, |y|, |z| ≤ L."""
        g = Substrate3DGeometry(N=32, L=12.0)
        sim = Substrate3DSimulator(geometry=g, n_frames=4, n_periods=1.0)
        sim.run()
        xs, ys, zs, vals = sim.iso_points(frame_index=0, level_frac=0.35)
        if len(xs) > 0:
            assert np.all(np.abs(xs) <= g.L + 1e-9)
            assert np.all(np.abs(ys) <= g.L + 1e-9)
            assert np.all(np.abs(zs) <= g.L + 1e-9)
            assert np.all(vals > 0.0)


# ---------------------------------------------------------------------------
# (c) Iso-surface volume oscillates with the breathing mode
# ---------------------------------------------------------------------------


class TestBreathingOscillation:
    """The bound state breathes — its centre amplitude swings, and the
    iso-volume V_iso(t) changes between frames."""

    def test_centre_amplitude_oscillates(self):
        sim = Substrate3DSimulator(
            geometry=Substrate3DGeometry(N=32, L=12.0),
            n_frames=12, n_periods=1.0,
        )
        sim.run()
        diag = sim.diagnostics()
        # AC component at least 0.1% of the initial DC level.  The seed
        # is intentionally small (5% of σ_max), so even small breathing
        # passes this gate without false negatives.
        assert diag["u_center_amp"] / max(abs(diag["u_center_init"]), 1e-30) > 1e-3

    def test_iso_volume_oscillates(self):
        """V_iso(t) at the canonical iso-level should change across
        frames — it is constant only if u(x,y,z) does not evolve."""
        sim = Substrate3DSimulator(
            geometry=Substrate3DGeometry(N=32, L=12.0),
            n_frames=12, n_periods=1.0,
        )
        sim.run()
        diag = sim.diagnostics()
        # Either absolute or relative oscillation must be nonzero
        assert diag["iso_vol_amp"] > 0, (
            "iso volume did not change between any two frames — "
            "evolution did not occur"
        )


# ---------------------------------------------------------------------------
# Snapshot machinery
# ---------------------------------------------------------------------------


class TestSnapshots:

    def test_snapshot_shape(self):
        """``u_frames`` has shape (n_frames, N, N, N) and is finite."""
        sim = Substrate3DSimulator(
            geometry=Substrate3DGeometry(N=32, L=12.0),
            n_frames=4, n_periods=1.0,
        )
        res = sim.run()
        assert res["u_frames"].shape == (4, 32, 32, 32)
        assert np.all(np.isfinite(res["u_frames"]))

    def test_diagnostics_well_formed(self):
        sim = Substrate3DSimulator(
            geometry=Substrate3DGeometry(N=32, L=12.0),
            n_frames=4, n_periods=1.0,
        )
        sim.run()
        diag = sim.diagnostics()
        for key in (
            "n_frames", "T_period", "iso_level",
            "energy_init", "energy_max_drift", "energy_drift_pct",
            "u_center_init", "u_center_amp",
            "iso_vol_mean", "iso_vol_amp", "iso_vol_amp_frac",
            "omega_pred", "rest_E_pred", "mass_pred_kg", "abs_max",
        ):
            assert key in diag
            assert np.isfinite(diag[key]), f"{key} not finite"

    def test_iso_points_subsamples_when_large(self):
        """``max_points`` cap is honoured."""
        sim = Substrate3DSimulator(
            geometry=Substrate3DGeometry(N=32, L=12.0),
            n_frames=2, n_periods=1.0,
        )
        sim.run()
        xs, ys, zs, vals = sim.iso_points(
            frame_index=0, level_frac=0.05,  # very low → many points
            max_points=500,
        )
        assert len(xs) <= 500
        assert len(xs) == len(ys) == len(zs) == len(vals)
