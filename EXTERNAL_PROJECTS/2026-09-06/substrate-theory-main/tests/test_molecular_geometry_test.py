"""Tests for stiff_medium.molecular_geometry_test.

Verifies:
  (a) K_4 forced angle = arccos(-1/3) = 109.4712°
  (b) Q_3 forced angle = 90° exactly
  (c) Tetrahedral family (7 species): pred error ≤ 0.05% (essentially exact)
  (d) Octahedral family (3 species): pred error ≤ 0.01%
  (e) Linear / trigonal-planar (no lone pair): pred error = 0
  (f) Trigonal-pyramidal & bent: substrate adopts literature VSEPR offsets
       (NOT a free fit — same Δ_lp any orbital theory uses)
  (g) ≥ 22 angle entries + ≥ 9 bond-length entries
  (h) Bond-length mean abs error < 1%
  (i) Visual file is created on disk
"""

from __future__ import annotations

import os
import math
import pytest

from stiff_medium.molecular_geometry_test import (
    REFERENCE_ANGLES,
    REFERENCE_LENGTHS,
    THETA_K4_DEG,
    THETA_Q3_DEG,
    THETA_SP_DEG,
    THETA_SP2_DEG,
    angle_residuals,
    length_residuals,
    family_summary,
    overall_stats,
    make_visual,
)


# ---------------------------------------------------------------------------
# (a, b) Substrate-FORCED angles
# ---------------------------------------------------------------------------

def test_k4_angle_is_arccos_minus_one_third() -> None:
    assert THETA_K4_DEG == pytest.approx(math.degrees(math.acos(-1.0 / 3.0)),
                                          abs=1e-12)
    assert THETA_K4_DEG == pytest.approx(109.4712, abs=1e-3)


def test_q3_angle_is_ninety() -> None:
    assert THETA_Q3_DEG == 90.0


def test_sp_and_sp2_angles_forced() -> None:
    assert THETA_SP_DEG  == 180.0
    assert THETA_SP2_DEG == 120.0


# ---------------------------------------------------------------------------
# (c) Tetrahedral family — substrate is essentially exact
# ---------------------------------------------------------------------------

def test_tetrahedral_family_essentially_exact() -> None:
    fam = [(r, e, p) for r, e, p in angle_residuals() if r.family == "tetrahedral"]
    assert len(fam) >= 7  # CH4, CCl4, SiH4, GeH4, NH4+, BF4-, ClO4-
    for ref, abs_err, rel_err in fam:
        # Experimental data is rounded to 0.01°; any disagreement reflects
        # rounding in the reference, not substrate inaccuracy.
        assert rel_err < 0.05, f"{ref.species}: {rel_err:.4f}%"


# ---------------------------------------------------------------------------
# (d) Octahedral family — Q_3 forces 90°
# ---------------------------------------------------------------------------

def test_octahedral_family_exact() -> None:
    fam = [(r, e, p) for r, e, p in angle_residuals() if r.family == "octahedral"]
    assert len(fam) >= 3
    for ref, abs_err, rel_err in fam:
        assert rel_err < 0.01, f"{ref.species}: {rel_err:.4f}%"


# ---------------------------------------------------------------------------
# (e) Linear / trigonal-planar — no lone pairs, predictions exact
# ---------------------------------------------------------------------------

def test_linear_and_trigonal_planar_exact() -> None:
    for ref, _, rel_err in angle_residuals():
        if ref.family in ("linear", "trigonal_planar"):
            assert rel_err == pytest.approx(0.0, abs=1e-9), \
                f"{ref.species} ({ref.family}): {rel_err:.4f}%"


# ---------------------------------------------------------------------------
# (f) Trigonal-pyramidal & bent — substrate adopts literature VSEPR offset
# ---------------------------------------------------------------------------

def test_lone_pair_families_match_vsepr() -> None:
    for ref, abs_err, rel_err in angle_residuals():
        if ref.family in ("trigonal_pyramid", "bent"):
            # Substrate uses the LITERATURE Δ_lp, so prediction matches
            # experiment to better than 1° (the reported precision).
            assert abs_err < 1.0, f"{ref.species}: {abs_err:.3f}°"


# ---------------------------------------------------------------------------
# (g) Coverage — at least 22 + 9
# ---------------------------------------------------------------------------

def test_coverage_min_counts() -> None:
    assert len(REFERENCE_ANGLES)  >= 22
    assert len(REFERENCE_LENGTHS) >= 9
    fams = {r.family for r in REFERENCE_ANGLES}
    assert fams == {"tetrahedral", "trigonal_pyramid", "bent",
                    "linear", "trigonal_planar", "octahedral"}


# ---------------------------------------------------------------------------
# (h) Bond lengths — mean abs error < 1%
# ---------------------------------------------------------------------------

def test_bond_length_mean_error_under_one_percent() -> None:
    s = overall_stats()
    assert s["length_mean_pct"] < 1.0
    assert s["length_max_pct"]  < 2.0


def test_bond_length_individual_within_two_percent() -> None:
    for ref, abs_err, rel_err in length_residuals():
        assert rel_err < 2.0, f"{ref.bond}: {rel_err:.3f}%"


# ---------------------------------------------------------------------------
# (i) Visual is written
# ---------------------------------------------------------------------------

def test_visual_is_created(tmp_path) -> None:
    out = make_visual(str(tmp_path / "131_molecular_geometry.png"))
    assert os.path.exists(out)
    assert os.path.getsize(out) > 5_000  # not an empty placeholder


# ---------------------------------------------------------------------------
# Aggregate sanity
# ---------------------------------------------------------------------------

def test_overall_stats_sane() -> None:
    s = overall_stats()
    assert s["n_total"] >= 31  # 22 angles + 9 lengths
    # Substrate's geometric predictions should be sub-percent overall
    assert s["angle_mean_pct"]  < 0.5
    assert s["length_mean_pct"] < 1.0
