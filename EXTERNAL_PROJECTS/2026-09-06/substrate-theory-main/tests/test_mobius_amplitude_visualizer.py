"""Tests for mobius_amplitude_visualizer -- 11/12 amplitude on K_4 + alpha.

Verifies:
  (a) Combinatorial route gives EXACTLY 11/12.
  (b) Monte Carlo agrees within 1%.
  (c) Quadrature equals 11/12 to machine precision.
  (d) Geometry inventory: 12 cells = 4 faces * 3 dihedrals; 1 pinched.
  (e) Moebius Z_2 twist identity s(theta + 2 pi) = -s(theta) holds.
  (f) alpha = 11 / (48 pi^3) * exp(-3 pi / 737) matches CODATA to <1%.
  (g) Factor breakdown reproduces the full alpha formula exactly.
  (h) Figure builders return matplotlib Figure objects without errors.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.mobius_amplitude_visualizer import (
    ALPHA_B3,
    ALPHA_CODATA,
    ANALYTIC_RATIO,
    INT_DRAG_RATIO,
    INT_K4_FACES,
    INT_SURVIVING,
    INT_TOTAL,
    AmplitudeReport,
    AmplitudeSimulator,
    MobiusAmplitudeGeometry,
    alpha_factor_breakdown,
    make_alpha_derivation_figure,
    make_mobius_k4_figure,
)


# ---------------------------------------------------------------------------
# Geometry inventory
# ---------------------------------------------------------------------------

class TestGeometry:
    def test_constants_match_alpha_formula(self):
        """The integers 11, 12, 4, 737 are exactly the formula constants."""
        assert INT_SURVIVING == 11
        assert INT_TOTAL == 12
        assert INT_K4_FACES == 4
        assert INT_DRAG_RATIO == 737

    def test_geometry_has_12_cells(self):
        """K_4 face-dihedral inventory: 4 faces * 3 dihedrals = 12 cells."""
        geom = MobiusAmplitudeGeometry()
        assert geom.n_cells == INT_TOTAL
        assert len(geom.faces) == INT_K4_FACES
        assert len(geom.edges) == 6   # C(4, 2)
        assert geom.vertices.shape == (4, 3)

    def test_one_cell_is_pinched(self):
        """Exactly one cell is removed by the Moebius identification."""
        geom = MobiusAmplitudeGeometry()
        assert geom.n_pinched == 1
        assert geom.n_surviving == INT_SURVIVING

    def test_pinch_index_default_is_zero(self):
        geom = MobiusAmplitudeGeometry()
        assert geom.pinch_index == 0
        assert geom.cells[0].is_pinched
        assert all(not c.is_pinched for c in geom.cells[1:])

    def test_pinch_index_can_be_changed(self):
        geom = MobiusAmplitudeGeometry(pinch_index=5)
        assert geom.cells[5].is_pinched
        assert geom.n_pinched == 1
        assert geom.surviving_ratio() == pytest.approx(11.0 / 12.0)

    def test_invalid_pinch_index_raises(self):
        with pytest.raises(ValueError):
            MobiusAmplitudeGeometry(pinch_index=-1)
        with pytest.raises(ValueError):
            MobiusAmplitudeGeometry(pinch_index=12)

    def test_cell_centroids_all_inside_tetrahedron(self):
        """Each cell centroid is interior to the tetrahedron (positive
        barycentric coords)."""
        geom = MobiusAmplitudeGeometry()
        V = geom.vertices
        # Build affine basis (V[1:]-V[0])
        basis = (V[1:] - V[0]).T
        for c in geom.cells:
            local = c.centroid - V[0]
            bary_rest = np.linalg.solve(basis, local)
            bary0 = 1.0 - bary_rest.sum()
            bary = np.concatenate([[bary0], bary_rest])
            # All non-negative within tolerance, sum to 1
            assert np.all(bary > -1e-9), f"cell {c.cell_index} outside: {bary}"
            assert bary.sum() == pytest.approx(1.0, abs=1e-9)

    def test_cell_base_angles_span_full_circle(self):
        """The 12 cells partition [0, 2 pi) at equal spacing."""
        geom = MobiusAmplitudeGeometry()
        angles = geom.cell_base_angles()
        expected = np.linspace(0.0, 2.0 * math.pi, geom.n_cells, endpoint=False)
        assert np.allclose(angles, expected)


# ---------------------------------------------------------------------------
# Moebius Z_2 transition function
# ---------------------------------------------------------------------------

class TestMobiusTwist:
    def test_z2_twist_identity_holds(self):
        """s(theta + 2 pi) = -s(theta) by construction."""
        geom = MobiusAmplitudeGeometry()
        assert geom.verify_z2_twist(n_samples=51)

    def test_z2_twist_explicit_samples(self):
        """Explicit numerical check on hand-picked theta values."""
        geom = MobiusAmplitudeGeometry()
        thetas = np.array([0.0, math.pi / 3, math.pi, 1.7 * math.pi])
        s = geom.mobius_section(thetas)
        s_shifted = geom.mobius_section(thetas + 2.0 * math.pi)
        assert np.allclose(s_shifted, -s, atol=1e-13)

    def test_trivial_section_is_periodic(self):
        """The reference section returns to itself after 2 pi."""
        geom = MobiusAmplitudeGeometry()
        thetas = np.linspace(0.0, math.pi, 17)
        s = geom.trivial_section(thetas)
        s_shifted = geom.trivial_section(thetas + 2.0 * math.pi)
        assert np.allclose(s_shifted, s, atol=1e-13)


# ---------------------------------------------------------------------------
# (a) Combinatorial 11/12
# ---------------------------------------------------------------------------

class TestCombinatorial:
    def test_combinatorial_ratio_is_exactly_11_over_12(self):
        sim = AmplitudeSimulator()
        ratio, info = sim.combinatorial_ratio()
        assert ratio == pytest.approx(11.0 / 12.0, rel=0, abs=1e-15)
        assert info["n_total"] == 12
        assert info["n_pinched"] == 1
        assert info["n_surviving"] == 11
        assert info["n_faces"] == 4
        assert info["dihedrals_per_face"] == 3

    def test_combinatorial_matches_analytic_constant(self):
        sim = AmplitudeSimulator()
        ratio, _ = sim.combinatorial_ratio()
        assert ratio == ANALYTIC_RATIO


# ---------------------------------------------------------------------------
# (b) Monte Carlo agrees within 1%
# ---------------------------------------------------------------------------

class TestMonteCarlo:
    def test_monte_carlo_agrees_with_analytic_within_1pct(self):
        """MC ratio must be within 1% of the analytic 11/12."""
        sim = AmplitudeSimulator()
        ratio, info = sim.monte_carlo_ratio(n_samples=200_000, seed=42)
        rel = abs(ratio - ANALYTIC_RATIO) / ANALYTIC_RATIO
        assert rel < 0.01, (
            f"Monte Carlo ratio {ratio} differs from 11/12 = "
            f"{ANALYTIC_RATIO} by {rel * 100:.3f}% (>1%)"
        )

    def test_monte_carlo_reports_std_err(self):
        sim = AmplitudeSimulator()
        _, info = sim.monte_carlo_ratio(n_samples=50_000, seed=7)
        assert "std_err" in info and info["std_err"] >= 0.0

    def test_monte_carlo_independent_of_seed_within_tolerance(self):
        sim = AmplitudeSimulator()
        r1, _ = sim.monte_carlo_ratio(n_samples=100_000, seed=1)
        r2, _ = sim.monte_carlo_ratio(n_samples=100_000, seed=2)
        # Both within 1% of 11/12, hence within 2% of each other
        assert abs(r1 - r2) / ANALYTIC_RATIO < 0.02


# ---------------------------------------------------------------------------
# Deterministic quadrature route
# ---------------------------------------------------------------------------

class TestQuadrature:
    def test_quadrature_is_exactly_11_over_12(self):
        sim = AmplitudeSimulator()
        ratio, info = sim.quadrature_ratio()
        assert ratio == pytest.approx(11.0 / 12.0, rel=1e-13)
        # Per-cell integral of cos^2 over [0, 2 pi] is pi
        assert info["cell_integral"] == pytest.approx(math.pi, rel=1e-12)


# ---------------------------------------------------------------------------
# Combined report
# ---------------------------------------------------------------------------

class TestReport:
    def test_run_all_assembles_three_routes(self):
        sim = AmplitudeSimulator()
        report = sim.run_all(n_samples=100_000)
        assert isinstance(report, AmplitudeReport)
        assert report.analytic == ANALYTIC_RATIO
        re = report.relative_errors()
        assert re["combinatorial"] < 1e-12
        assert re["monte_carlo"] < 0.01
        assert re["quadrature"] < 1e-12


# ---------------------------------------------------------------------------
# (c) alpha = 11 / (48 pi^3) * exp(-3 pi / 737) <1% of CODATA
# ---------------------------------------------------------------------------

class TestAlphaDerivation:
    def test_alpha_b3_matches_codata_within_1pct(self):
        """The alpha formula must match CODATA to better than 1%."""
        rel = abs(ALPHA_B3 - ALPHA_CODATA) / ALPHA_CODATA
        assert rel < 0.01, (
            f"alpha_B3 = {ALPHA_B3} differs from alpha_CODATA = "
            f"{ALPHA_CODATA} by {rel * 100:.4f}% (>1%)"
        )

    def test_alpha_b3_is_within_5e_minus_4(self):
        """Tighter sanity check: framework anchors at ~4e-5 relative."""
        rel = abs(ALPHA_B3 - ALPHA_CODATA) / ALPHA_CODATA
        assert rel < 5e-4

    def test_alpha_factor_breakdown_is_consistent(self):
        """Multiplying all factors reproduces ALPHA_B3 exactly."""
        factors, alpha_pred, rel_err = alpha_factor_breakdown()
        assert len(factors) == 4
        product = 1.0
        for f in factors:
            product *= f.value
        assert alpha_pred == pytest.approx(product, rel=1e-15)
        assert alpha_pred == pytest.approx(ALPHA_B3, rel=1e-15)

    def test_alpha_factor_includes_11_over_12(self):
        """The Moebius factor 11/12 must appear explicitly."""
        factors, _, _ = alpha_factor_breakdown()
        labels = [f.label for f in factors]
        assert "11 / 12" in labels
        # Find it and check value
        f_mobius = next(f for f in factors if f.label == "11 / 12")
        assert f_mobius.value == pytest.approx(11.0 / 12.0, rel=1e-15)

    def test_alpha_factor_includes_exponential(self):
        factors, _, _ = alpha_factor_breakdown()
        labels = [f.label for f in factors]
        assert any("exp" in lbl for lbl in labels)
        f_exp = next(f for f in factors if "exp" in f.label)
        assert f_exp.value == pytest.approx(
            math.exp(-3.0 * math.pi / INT_DRAG_RATIO), rel=1e-15
        )


# ---------------------------------------------------------------------------
# Visualization builders
# ---------------------------------------------------------------------------

class TestFigures:
    def test_make_mobius_k4_figure_returns_figure(self):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.figure import Figure as MplFigure
        fig = make_mobius_k4_figure()
        assert isinstance(fig, MplFigure)
        # Two axes: the 3D K_4 panel + the 2D theta-section panel
        assert len(fig.axes) == 2
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_make_alpha_derivation_figure_returns_figure(self):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.figure import Figure as MplFigure
        fig = make_alpha_derivation_figure()
        assert isinstance(fig, MplFigure)
        # 3 axes: bar chart + running product + caption panel
        assert len(fig.axes) == 3
        import matplotlib.pyplot as plt
        plt.close(fig)
