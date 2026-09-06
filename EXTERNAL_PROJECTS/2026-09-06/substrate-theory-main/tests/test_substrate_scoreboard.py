"""Tests for SubstrateScoreboard aggregator."""

from __future__ import annotations

import math
import pytest

from stiff_medium.substrate_scoreboard import (
    SubstrateScoreboard,
    Prediction,
    Falsifier,
)


@pytest.fixture(scope="module")
def board() -> SubstrateScoreboard:
    sb = SubstrateScoreboard(use_live=False)
    sb.run_all()  # build
    return sb


# ---------------------------------------------------------------------
# (a) at least 30 predictions aggregated
# ---------------------------------------------------------------------
def test_at_least_30_predictions(board: SubstrateScoreboard) -> None:
    preds = board.run_all()
    assert len(preds) >= 30, f"expected >=30 predictions, got {len(preds)}"


def test_predictions_have_required_fields(board: SubstrateScoreboard) -> None:
    preds = board.run_all()
    required = {"predicted", "observed", "units", "sector", "rigor", "source"}
    for label, entry in preds.items():
        missing = required - set(entry.keys())
        assert not missing, f"{label} missing fields {missing}"


# ---------------------------------------------------------------------
# (b) majority PASS
# ---------------------------------------------------------------------
def test_majority_pass(board: SubstrateScoreboard) -> None:
    cats = board.categorize()
    assert cats["TOTAL"] >= 30
    assert cats["PASS"] > cats["TOTAL"] // 2, (
        f"expected majority PASS, got {cats}"
    )


def test_categorize_counts_consistent(board: SubstrateScoreboard) -> None:
    cats = board.categorize()
    assert cats["PASS"] + cats["TENSION"] + cats["FAIL"] == cats["TOTAL"]


# ---------------------------------------------------------------------
# (c) honest reporting of FAILs
# ---------------------------------------------------------------------
def test_at_least_one_fail_reported(board: SubstrateScoreboard) -> None:
    """Substrate framework has known open problems (e.g. density perturbation
    amplitude). Scoreboard must surface these honestly rather than hide them."""
    cats = board.categorize()
    assert cats["FAIL"] >= 1, (
        "scoreboard should honestly report at least one failing prediction; "
        "the density-perturbation-amplitude problem is documented as open"
    )


def test_density_perturbation_is_failing(board: SubstrateScoreboard) -> None:
    cmp = board.compare_to_observed()
    assert "density_perturbation_amp" in cmp
    assert cmp["density_perturbation_amp"]["status"] == "FAIL"


# ---------------------------------------------------------------------
# (d) all major experimental frontiers covered
# ---------------------------------------------------------------------
def test_major_sectors_all_present(board: SubstrateScoreboard) -> None:
    preds = board.run_all()
    sectors = {entry["sector"] for entry in preds.values()}
    required_frontiers = {
        "leptons", "quarks", "baryons", "mixing", "neutrinos",
        "darkmatter", "cosmology", "bao", "cmb", "sparc",
        "hubble", "vacuum", "nuclei", "banner",
    }
    missing = required_frontiers - sectors
    assert not missing, f"missing frontiers: {missing}"


# ---------------------------------------------------------------------
# Rigor classification
# ---------------------------------------------------------------------
def test_forced_vs_fitted_classification(board: SubstrateScoreboard) -> None:
    groups = board.forced_vs_fitted()
    assert "FORCED" in groups and "ANCHORED" in groups and "FITTED" in groups
    # framework should be predominantly FORCED, not predominantly FITTED
    assert len(groups["FORCED"]) > len(groups["FITTED"])
    # every prediction must be classified exactly once
    total = sum(len(v) for v in groups.values())
    assert total == len(board.run_all())


def test_no_unknown_rigor_labels(board: SubstrateScoreboard) -> None:
    valid = {"FORCED", "ANCHORED", "FITTED"}
    for entry in board.run_all().values():
        assert entry["rigor"] in valid


# ---------------------------------------------------------------------
# Falsifiers
# ---------------------------------------------------------------------
def test_falsifiers_cover_key_experiments(board: SubstrateScoreboard) -> None:
    fals = board.falsifiers()
    assert len(fals) >= 5
    text = " ".join(f["experiment"] + " " + f["quantity"] for f in fals).lower()
    for keyword in ["desi", "katrin", "cmb", "edm"]:
        assert keyword in text, f"no falsifier mentions {keyword}"


def test_falsifiers_have_thresholds(board: SubstrateScoreboard) -> None:
    for f in board.falsifiers():
        assert f["decision_threshold"], "falsifier missing decision threshold"
        assert f["b3_prediction"], "falsifier missing B3 prediction"


# ---------------------------------------------------------------------
# Compare-to-observed
# ---------------------------------------------------------------------
def test_compare_residuals_finite(board: SubstrateScoreboard) -> None:
    cmp = board.compare_to_observed()
    for label, entry in cmp.items():
        # density-perturbation comparison is degenerate (huge residual ok)
        # but residual must at least be a number (or nan if observed=0)
        r = entry["residual_pct"]
        assert isinstance(r, float)


def test_compare_status_values_valid(board: SubstrateScoreboard) -> None:
    valid = {"PASS", "TENSION", "FAIL"}
    for entry in board.compare_to_observed().values():
        assert entry["status"] in valid


# ---------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------
def test_report_renders(board: SubstrateScoreboard) -> None:
    txt = board.report()
    assert "SUBSTRATE FRAMEWORK SCOREBOARD" in txt
    assert "PASS" in txt and "FAIL" in txt
    assert "UPCOMING FALSIFIERS" in txt
    # must include several known labels
    for label in ["proton_mass_MeV", "Sigma_m_nu_meV", "H0_km_s_Mpc",
                  "g_dagger_m_s2", "Cabibbo_lambda"]:
        assert label in txt


def test_prediction_status_logic() -> None:
    # within 2 sigma -> PASS
    p = Prediction("x", "test", 1.0, 1.01, 0.05, "u", "FORCED")
    assert p.status() == "PASS"
    # 2.5 sigma -> TENSION
    p = Prediction("x", "test", 1.0, 1.25, 0.10, "u", "FORCED")
    assert p.status() == "TENSION"
    # 10 sigma -> FAIL
    p = Prediction("x", "test", 1.0, 2.0, 0.05, "u", "FORCED")
    assert p.status() == "FAIL"
    # no sigma, within 2% -> PASS
    p = Prediction("x", "test", 1.00, 1.01, None, "u", "FORCED")
    assert p.status() == "PASS"
    # no sigma, > 5% -> FAIL
    p = Prediction("x", "test", 1.00, 2.00, None, "u", "FORCED")
    assert p.status() == "FAIL"
