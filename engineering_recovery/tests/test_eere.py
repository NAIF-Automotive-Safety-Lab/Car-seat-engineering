from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engineering_recovery.acquisition.service import ArtifactAcquisitionService
from engineering_recovery.brep.validate import brep_report
from engineering_recovery.step.forensic import forensic_report


STEP = ROOT / "R4.1" / "R4.1.step"
EXPECTED_SHA = "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"


def test_acquisition_manifest_and_replay(tmp_path: Path) -> None:
    service = ArtifactAcquisitionService(ROOT, "test-revision")
    acquired, manifest = service.acquire_file(STEP, tmp_path / "evidence")
    assert manifest.sha256 == EXPECTED_SHA
    assert manifest.file_size_bytes == STEP.stat().st_size
    assert service.verify_identity(acquired, manifest)
    assert service.replay_artifact(acquired, manifest)["status"] == "PASS"


def test_step_forensic_is_reproducible() -> None:
    first = forensic_report(STEP)
    second = forensic_report(STEP)
    assert first == second
    assert first["source_sha256"] == EXPECTED_SHA
    assert first["header_present"] is True
    assert first["footer_present"] is True
    assert first["topology_entity_count"] > 0


def test_brep_is_reproducible_and_real() -> None:
    first = brep_report(STEP)
    second = brep_report(STEP)
    assert first["source_sha256"] == EXPECTED_SHA
    assert first["solid_count"] == 62
    assert first["shape_validity"] is True
    assert first["bounding_box"] == second["bounding_box"]
    assert first["volume"] == pytest.approx(second["volume"])


def test_evidence_artifact_sha_matches_source() -> None:
    digest = hashlib.sha256(STEP.read_bytes()).hexdigest()
    assert digest == EXPECTED_SHA
