from __future__ import annotations

from pathlib import Path

from engineering_recovery.features.deterministic import feature_report
from engineering_recovery.pmi.extract import pmi_report

ROOT = Path(__file__).resolve().parents[2]
STEP = ROOT / "R4.1" / "R4.1.step"
EXPECTED_SHA = "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"


def test_feature_inventory_is_deterministic_and_nonsemantic() -> None:
    first = feature_report(STEP)
    second = feature_report(STEP)
    assert first["status"] == "VERIFIED"
    assert first["source_sha256"] == EXPECTED_SHA
    assert first["surface_type_counts"] == second["surface_type_counts"]
    assert first["candidate_features"]["holes"] == "UNDEFINED"


def test_pmi_absence_is_reported_without_inference() -> None:
    report = pmi_report(STEP)
    assert report["source_sha256"] == EXPECTED_SHA
    assert report["pmi_entity_count"] == 0
    assert report["status"] == "NOT_PRESENT_OR_NOT_RECOVERED"
    assert report["dimensions"] == []
