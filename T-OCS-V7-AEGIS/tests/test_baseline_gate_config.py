from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import pytest

AEGIS_ROOT = Path(__file__).resolve().parents[1]
if str(AEGIS_ROOT) not in sys.path:
    sys.path.insert(0, str(AEGIS_ROOT))

from aegis import gate_status


def _valid_lock(expected: str) -> dict:
    return {
        "file": "R4.1/R4.1.step",
        "expected_sha256": expected,
        "sha256_before": expected,
        "sha256_after": expected,
        "sha_match": True,
        "solid_count": 62,
        "solid_count_expected": 62,
        "solid_count_pass": True,
        "step_read_status": "IFSelect_ReturnStatus.IFSelect_RetDone",
        "brep_valid": True,
        "mutation_detected": False,
        "status": "PASS",
    }


def test_gate_baseline_config_001_reads_lock_field(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """GATE-BASELINE-CONFIG-001: fail if baseline SHA is hard-coded."""
    artifact = tmp_path / "R4.1" / "R4.1.step"
    artifact.parent.mkdir(parents=True)
    artifact_bytes = b"baseline-gate-config-001-runtime-artifact"
    artifact.write_bytes(artifact_bytes)
    actual = hashlib.sha256(artifact_bytes).hexdigest()

    lock = tmp_path / "R4_1_BASELINE_LOCK.json"
    lock.write_text(json.dumps(_valid_lock(actual)), encoding="utf-8")

    monkeypatch.setattr(gate_status, "sha256_file", lambda _path: actual)
    monkeypatch.setattr(gate_status, "_native_step_check", lambda _path: ("PASS", None, 62))

    status, reasons, info = gate_status.baseline_gate(tmp_path)
    assert status == "PASS", reasons
    assert info["actual_sha256"] == actual
    assert info["reference_sha256"] == actual

    replacement = hashlib.sha256(b"different-governing-value").hexdigest()
    lock.write_text(json.dumps(_valid_lock(replacement)), encoding="utf-8")
    status2, reasons2, info2 = gate_status.baseline_gate(tmp_path)
    assert status2 == "BLOCKED"
    assert any("SHA mismatch" in reason for reason in reasons2)
    assert info2["reference_sha256"] == replacement


def test_gate_baseline_config_001_has_no_embedded_sha256_literal() -> None:
    source = (AEGIS_ROOT / "aegis" / "gate_status.py").read_text(encoding="utf-8")
    assert "EXPECTED_R41_SHA256" not in source
    embedded = re.findall(r"(?<![0-9a-fA-F])[0-9a-fA-F]{64}(?![0-9a-fA-F])", source)
    assert embedded == []
