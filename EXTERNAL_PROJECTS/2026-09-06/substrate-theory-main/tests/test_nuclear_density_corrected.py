"""
test_nuclear_density_corrected.py

Tests for the corrected nuclear matter saturation density derivation
in src/stiff_medium/nuclear_density_corrected.py.
"""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.nuclear_density_corrected import (
    HBARC_MEV_FM,
    LAMBDA_QCD_MEV,
    M_PI_PLUS_MEV,
    RHO_0_OBS_FM3,
    corrected_density,
    fcc_density,
    fcc_volume_per_atom,
    k4_edge_internal,
    original_geom_07_density,
    pion_yukawa_density,
    report,
    tetrahedron_circumradius,
    tetrahedron_height,
    tetrahedron_inradius,
    tetrahedron_volume,
)


# ---------------------------------------------------------------------------
# Geometry sanity checks for the regular tetrahedron
# ---------------------------------------------------------------------------
class TestTetrahedronGeometry:

    def test_volume_formula(self):
        """V = a^3/(6 sqrt(2)) for a regular tetrahedron."""
        a = 1.0
        v = tetrahedron_volume(a)
        assert v == pytest.approx(1.0 / (6.0 * math.sqrt(2.0)), rel=1e-12)

    def test_height_formula(self):
        """h = a sqrt(2/3)."""
        a = 1.0
        h = tetrahedron_height(a)
        assert h == pytest.approx(math.sqrt(2.0 / 3.0), rel=1e-12)

    def test_circumradius_formula(self):
        """R_out = a sqrt(6)/4."""
        a = 1.0
        R = tetrahedron_circumradius(a)
        assert R == pytest.approx(math.sqrt(6.0) / 4.0, rel=1e-12)

    def test_inradius_formula(self):
        """r_in = a/(2 sqrt(6))."""
        a = 1.0
        r = tetrahedron_inradius(a)
        assert r == pytest.approx(1.0 / (2.0 * math.sqrt(6.0)), rel=1e-12)

    def test_inradius_circumradius_ratio(self):
        """For a regular tetrahedron R_out / r_in = 3."""
        a = 1.234
        ratio = tetrahedron_circumradius(a) / tetrahedron_inradius(a)
        assert ratio == pytest.approx(3.0, rel=1e-12)

    def test_circumradius_is_3_times_centroid_to_face(self):
        """Centroid-to-face distance = h/4; R_out = 3 (h/4)."""
        a = 2.5
        h = tetrahedron_height(a)
        R = tetrahedron_circumradius(a)
        # Centroid sits at h/4 from base, R = h - h/4 = 3h/4.
        assert R == pytest.approx(3.0 * h / 4.0, rel=1e-12)


# ---------------------------------------------------------------------------
# K_4 cell at the QCD scale
# ---------------------------------------------------------------------------
class TestK4Cell:

    def test_internal_edge_value(self):
        """a_q = ℏc/Λ_QCD ≈ 0.987 fm with the B3 anchor 200 MeV."""
        a = k4_edge_internal()
        assert a == pytest.approx(HBARC_MEV_FM / LAMBDA_QCD_MEV, rel=1e-12)
        assert a == pytest.approx(0.9866, abs=0.001)

    def test_circumradius_close_to_proton_charge_radius(self):
        """R_out for the K_4 ≈ 0.6 fm; same OOM as proton charge radius 0.84 fm."""
        a = k4_edge_internal()
        R = tetrahedron_circumradius(a)
        assert 0.5 < R < 0.7
        assert R == pytest.approx(0.604, abs=0.005)


# ---------------------------------------------------------------------------
# FCC packing arithmetic
# ---------------------------------------------------------------------------
class TestFCC:

    def test_fcc_volume_per_atom_formula(self):
        """V/atom = d_NN^3 / sqrt(2)."""
        d = 2.0
        v = fcc_volume_per_atom(d)
        assert v == pytest.approx(8.0 / math.sqrt(2.0), rel=1e-12)

    def test_fcc_density_formula(self):
        """rho = sqrt(2) / d_NN^3."""
        d = 2.0
        rho = fcc_density(d)
        assert rho == pytest.approx(math.sqrt(2.0) / 8.0, rel=1e-12)

    def test_density_volume_consistency(self):
        """rho * V/atom == 1 always."""
        d = 1.7
        assert fcc_density(d) * fcc_volume_per_atom(d) == pytest.approx(1.0, rel=1e-12)

    def test_fcc_relation_to_cubic_edge(self):
        """In FCC the conventional cubic edge a_FCC = sqrt(2) * d_NN
        and 4 atoms per cubic cell → V/atom = a_FCC^3/4."""
        d = 1.6
        a_fcc = math.sqrt(2.0) * d
        v_per_atom_cubic = a_fcc ** 3 / 4.0
        v_per_atom_fcc   = fcc_volume_per_atom(d)
        assert v_per_atom_cubic == pytest.approx(v_per_atom_fcc, rel=1e-12)


# ---------------------------------------------------------------------------
# The original (broken) prediction is reproduced and then beaten
# ---------------------------------------------------------------------------
class TestOriginalVsCorrected:

    def test_original_is_too_dense(self):
        """Original gives ~1.47 fm^-3, ~9x too dense."""
        rho = original_geom_07_density()
        assert rho == pytest.approx(1.4725, abs=0.001)
        assert rho > 5.0 * RHO_0_OBS_FM3   # confirms the gap

    def test_original_relative_error_is_huge(self):
        """Confirms the originally-reported 820% error."""
        rho = original_geom_07_density()
        err = abs(rho - RHO_0_OBS_FM3) / RHO_0_OBS_FM3 * 100.0
        assert err > 800.0
        assert err < 850.0

    def test_corrected_is_much_better(self):
        """Corrected derivation reduces error by >50x."""
        rho_orig = original_geom_07_density()
        rho_new  = corrected_density()
        err_orig = abs(rho_orig - RHO_0_OBS_FM3) / RHO_0_OBS_FM3
        err_new  = abs(rho_new  - RHO_0_OBS_FM3) / RHO_0_OBS_FM3
        assert err_new < err_orig / 50.0

    def test_corrected_within_50_percent(self):
        """Substrate prediction is within the 50% requested tolerance."""
        rho = corrected_density()
        err_pct = abs(rho - RHO_0_OBS_FM3) / RHO_0_OBS_FM3 * 100.0
        assert err_pct < 50.0

    def test_corrected_within_20_percent(self):
        """Tighter check: error is actually well under 20%."""
        rho = corrected_density()
        err_pct = abs(rho - RHO_0_OBS_FM3) / RHO_0_OBS_FM3 * 100.0
        assert err_pct < 20.0

    def test_corrected_specific_value(self):
        """rho_corrected = sqrt(2) / (2 a_q)^3 with a_q = ℏc/Λ_QCD."""
        a_q  = HBARC_MEV_FM / LAMBDA_QCD_MEV
        d_nn = 2.0 * a_q
        expected = math.sqrt(2.0) / d_nn ** 3
        assert corrected_density() == pytest.approx(expected, rel=1e-12)
        assert corrected_density() == pytest.approx(0.184, abs=0.005)


# ---------------------------------------------------------------------------
# Pion-Yukawa cross-check
# ---------------------------------------------------------------------------
class TestPionYukawa:

    def test_pion_density_within_15_pct(self):
        """Cross-check with d_NN = sqrt(2) λ_π gives rho ≈ 0.177."""
        rho = pion_yukawa_density()
        err_pct = abs(rho - RHO_0_OBS_FM3) / RHO_0_OBS_FM3 * 100.0
        assert err_pct < 15.0

    def test_pion_density_specific(self):
        lam_pi = HBARC_MEV_FM / M_PI_PLUS_MEV
        d_nn   = math.sqrt(2.0) * lam_pi
        expected = math.sqrt(2.0) / d_nn ** 3
        assert pion_yukawa_density() == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------
# Full report runs cleanly and returns sane structure
# ---------------------------------------------------------------------------
class TestReport:

    def test_report_runs_and_returns_three_results(self):
        results = report()
        assert set(results.keys()) == {
            "original_step10",
            "corrected_2aq",
            "pion_yukawa",
        }

    def test_report_corrected_better_than_original(self):
        results = report()
        assert results["corrected_2aq"].rel_error_pct < \
            results["original_step10"].rel_error_pct
