"""Tests for generation_tower_visualizer -- 3-generation lepton/quark tower.

Verifies:
  (a) 3 generations are FORCED from D = 3 spatial dimensions.
  (b) m_mu / m_e prediction matches PDG to 0.01%.
  (c) m_tau / m_mu prediction matches PDG within 1%.
  (d) Geometry: three orthogonal K_4 cells, one per substrate axis.
  (e) Renderer entry points produce non-empty matplotlib figures.
"""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from src.stiff_medium.b3_constants import K_pair, R, n_M, n_R
from src.stiff_medium.generation_map import LEPTON_MASSES_MEV
from src.stiff_medium.generation_tower_visualizer import (
    AXIS_NAMES,
    D_SPATIAL,
    GenerationTowerGeometry,
    GenerationTowerSimulator,
    LEPTON_LABELS,
    make_tower_figures,
)
from src.stiff_medium.tau_mass_unified import M_MUON_MEV, M_TAU_MEV


# ---------------------------------------------------------------------------
# (a) 3 generations forced from D = 3
# ---------------------------------------------------------------------------

class TestThreeGenerationsForced:

    def test_default_spatial_dim_is_three(self):
        geom = GenerationTowerGeometry()
        assert geom.spatial_dim == D_SPATIAL == 3

    def test_n_generations_equals_spatial_dim(self):
        """One generation per orthogonal substrate axis."""
        geom = GenerationTowerGeometry()
        assert geom.n_generations == geom.spatial_dim == 3

    def test_three_orthogonal_axis_directions(self):
        """The three axis directions are pairwise orthogonal unit vectors."""
        geom = GenerationTowerGeometry()
        ax = geom.axis_directions
        assert ax.shape == (3, 3)
        # Each axis is unit-length.
        for i in range(3):
            assert ax[i] @ ax[i] == pytest.approx(1.0, rel=1e-12)
        # Pairwise orthogonal.
        assert ax[0] @ ax[1] == pytest.approx(0.0, abs=1e-12)
        assert ax[0] @ ax[2] == pytest.approx(0.0, abs=1e-12)
        assert ax[1] @ ax[2] == pytest.approx(0.0, abs=1e-12)

    def test_one_K4_cell_per_generation(self):
        """Three K_4 cells, one per axis."""
        geom = GenerationTowerGeometry()
        assert len(geom.cells) == 3
        for cell in geom.cells:
            assert cell.n_vertices == 4
            assert cell.n_edges == 6
            assert cell.n_faces == 4

    def test_is_3gen_forced_returns_true(self):
        geom = GenerationTowerGeometry()
        assert geom.is_3gen_forced() is True

    def test_fourth_generation_blocked_by_dim(self):
        """A spatial_dim != 3 should be rejected (no fourth orthogonal axis)."""
        with pytest.raises(ValueError):
            GenerationTowerGeometry(spatial_dim=4)
        with pytest.raises(ValueError):
            GenerationTowerGeometry(spatial_dim=2)

    def test_axis_for_generation_returns_xyz(self):
        geom = GenerationTowerGeometry()
        assert geom.axis_for_generation(1) == "x"
        assert geom.axis_for_generation(2) == "y"
        assert geom.axis_for_generation(3) == "z"
        with pytest.raises(ValueError):
            geom.axis_for_generation(4)

    def test_axis_names_consistent(self):
        assert AXIS_NAMES == ("x", "y", "z")
        assert LEPTON_LABELS == ("electron", "muon", "tau")


# ---------------------------------------------------------------------------
# (b) m_mu / m_e to 0.01%
# ---------------------------------------------------------------------------

class TestMuonElectronRatio:

    def test_mu_over_e_within_0_01_percent(self):
        geom = GenerationTowerGeometry()
        pred = geom.ratio_12()
        obs = LEPTON_MASSES_MEV["mu"] / LEPTON_MASSES_MEV["e"]
        err_pct = 100.0 * (pred - obs) / obs
        assert abs(err_pct) < 0.01, (
            f"m_mu/m_e prediction {pred:.6f} vs PDG {obs:.6f} -> err = {err_pct:.4f}% "
            f"(target: <0.01%)"
        )

    def test_log_ratio_12_matches_268_over_16pi(self):
        geom = GenerationTowerGeometry()
        assert geom.log_ratio_12() == pytest.approx(
            268.0 / (16.0 * math.pi), rel=1e-12
        )

    def test_predicted_muon_mass_close_to_pdg(self):
        geom = GenerationTowerGeometry()
        m_mu_pred = geom.predicted_lepton_masses_mev()["mu"]
        m_mu_obs = LEPTON_MASSES_MEV["mu"]
        err = 100.0 * (m_mu_pred - m_mu_obs) / m_mu_obs
        assert abs(err) < 0.01


# ---------------------------------------------------------------------------
# (c) m_tau / m_mu within 1%
# ---------------------------------------------------------------------------

class TestTauMuonRatio:

    def test_tau_over_mu_within_1_percent(self):
        geom = GenerationTowerGeometry()
        pred = geom.ratio_23()
        obs = M_TAU_MEV / M_MUON_MEV
        err_pct = 100.0 * (pred - obs) / obs
        assert abs(err_pct) < 1.0, (
            f"m_tau/m_mu prediction {pred:.6f} vs PDG {obs:.6f} -> err = {err_pct:.4f}% "
            f"(target: <1%)"
        )

    def test_log_ratio_23_uses_canonical_formula(self):
        """log(m_tau/m_mu) = n_M/(K_pair**4 pi) - (n_R - R)/(K_pair*(K_pair+1))."""
        geom = GenerationTowerGeometry()
        expected = (n_M / (K_pair ** 4 * math.pi)
                    - (n_R - R) / (K_pair * (K_pair + 1)))
        assert geom.log_ratio_23() == pytest.approx(expected, rel=1e-12)

    def test_ratio_13_equals_ratio_12_times_ratio_23(self):
        geom = GenerationTowerGeometry()
        assert geom.ratio_13() == pytest.approx(
            geom.ratio_12() * geom.ratio_23(), rel=1e-12
        )


# ---------------------------------------------------------------------------
# (d) Geometric tower (orthogonal K_4 cells)
# ---------------------------------------------------------------------------

class TestGeometricTower:

    def test_cells_translated_along_orthogonal_axes(self):
        """Each cell's centre lies on a distinct orthogonal axis."""
        geom = GenerationTowerGeometry()
        spacing = geom.axis_spacing
        for gen, cell in enumerate(geom.cells):
            expected = spacing * geom.axis_directions[gen]
            assert np.allclose(cell.centre, expected, atol=1e-12), (
                f"cell {gen} centre {cell.centre} != expected {expected}"
            )

    def test_cell_centres_are_orthogonal_from_origin(self):
        geom = GenerationTowerGeometry()
        c0, c1, c2 = (cell.centre for cell in geom.cells)
        assert c0 @ c1 == pytest.approx(0.0, abs=1e-12)
        assert c0 @ c2 == pytest.approx(0.0, abs=1e-12)
        assert c1 @ c2 == pytest.approx(0.0, abs=1e-12)

    def test_summary_contains_required_keys(self):
        sim = GenerationTowerSimulator()
        s = sim.geometry.summary()
        for key in ("spatial_dim", "n_generations", "axes",
                    "log_ratio_12", "log_ratio_23",
                    "ratio_12_pred", "ratio_23_pred", "ratio_13_pred",
                    "predicted_masses_mev"):
            assert key in s


# ---------------------------------------------------------------------------
# (e) Renderer entry points
# ---------------------------------------------------------------------------

class TestRenderers:

    def test_draw_geometric_tower_returns_axis(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        sim = GenerationTowerSimulator()
        ret = sim.draw_geometric_tower(ax=ax)
        assert ret is ax
        # Title is non-empty.
        assert len(ax.get_title()) > 0
        plt.close(fig)

    def test_draw_full_ladder_returns_axis(self):
        fig, ax = plt.subplots()
        sim = GenerationTowerSimulator()
        ret = sim.draw_full_ladder(ax=ax)
        assert ret is ax
        # The y-axis should be log-scale for the mass ladder.
        assert ax.get_yscale() == "log"
        # 3 ticks on x for 3 generations.
        assert len(ax.get_xticks()) == 3
        plt.close(fig)

    def test_make_tower_figures_returns_two_figures(self):
        fig_tower, fig_ladder = make_tower_figures()
        assert fig_tower is not None
        assert fig_ladder is not None
        plt.close(fig_tower)
        plt.close(fig_ladder)

    def test_lepton_errors_pct_within_targets(self):
        sim = GenerationTowerSimulator()
        errs = sim.lepton_errors_pct()
        assert abs(errs["mu/e"]) < 0.01
        assert abs(errs["tau/mu"]) < 1.0
        # tau/e error compounds 1->2 and 2->3.
        assert abs(errs["tau/e"]) < 1.0
