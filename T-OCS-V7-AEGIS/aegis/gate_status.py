from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

from .provenance import sha256_file

PROJECT_ROOT = Path(__file__).resolve().parents[2]
AEGIS_ROOT = PROJECT_ROOT / "T-OCS-V7-AEGIS"
BASELINE_ARTIFACT = PROJECT_ROOT / "R4.1" / "R4.1.step"
BASELINE_LOCK = PROJECT_ROOT / "R4_1_BASELINE_LOCK.json"
VALID_STATES = {"PASS", "FAIL", "BLOCKED", "NOT_PROVEN"}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


def load_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return value if isinstance(value, dict) else {}


def read_gate_record(path: Path) -> tuple[str, dict]:
    obj = load_json(path)
    value = obj.get("gate_status", obj.get("status", "NOT_PROVEN"))
    return (value if value in VALID_STATES else "NOT_PROVEN"), obj


def _native_step_check(path: Path) -> tuple[str, str | None, int | None]:
    try:
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.STEPControl import STEPControl_Reader
        from OCP.TopAbs import TopAbs_SOLID
        from OCP.TopExp import TopExp_Explorer
    except Exception as exc:
        return "BLOCKED", f"OCP unavailable: {exc}", None

    try:
        reader = STEPControl_Reader()
        status = reader.ReadFile(str(path))
        if status != IFSelect_RetDone:
            return "FAIL", f"STEP read status={status}", None
        transferred = reader.TransferRoots()
        if transferred < 1:
            return "FAIL", f"STEP transfer roots={transferred}", None
        shape = reader.OneShape()
        if shape.IsNull():
            return "FAIL", "B-Rep OneShape is null", None
        explorer = TopExp_Explorer(shape, TopAbs_SOLID)
        solids = 0
        while explorer.More():
            solids += 1
            explorer.Next()
        return "PASS", None, solids
    except Exception as exc:
        return "FAIL", f"native STEP/B-Rep check error: {exc}", None


def baseline_gate(project_root: Path | None = None) -> tuple[str, list[str], dict]:
    root = Path(project_root).resolve() if project_root is not None else PROJECT_ROOT
    artifact = root / "R4.1" / "R4.1.step"
    lock_path = root / "R4_1_BASELINE_LOCK.json"
    blockers: list[str] = []

    if not artifact.exists():
        return "BLOCKED", ["R4.1/R4.1.step is absent from execution checkout"], {}
    if artifact.is_symlink() or not artifact.is_file():
        return "BLOCKED", ["R4.1 artifact is not a regular non-symlink file"], {}

    try:
        lock_raw = lock_path.read_text(encoding="utf-8")
        lock = json.loads(lock_raw)
    except Exception as exc:
        return "BLOCKED", [f"R4_1_BASELINE_LOCK.json unreadable/invalid JSON: {exc}"], {}
    if not isinstance(lock, dict):
        return "BLOCKED", ["R4_1_BASELINE_LOCK.json root is not an object"], {}

    if "expected_sha256" not in lock:
        return "BLOCKED", ["R4_1_BASELINE_LOCK.json missing expected_sha256"], {}
    expected = lock["expected_sha256"]
    if not isinstance(expected, str):
        return "BLOCKED", ["expected_sha256 is not a string"], {"reference": expected}
    if len(expected) != 64:
        return "BLOCKED", [f"expected_sha256 length={len(expected)}; required=64"], {
            "reference_sha256": expected,
            "reference_length": len(expected),
        }
    if SHA256_RE.fullmatch(expected) is None:
        return "BLOCKED", ["expected_sha256 is not 64 hexadecimal characters"], {
            "reference_sha256": expected,
            "reference_length": len(expected),
        }

    actual = sha256_file(artifact)
    if len(actual) != 64 or SHA256_RE.fullmatch(actual) is None:
        return "FAIL", ["computed artifact SHA-256 is not a valid 64-character hex digest"], {
            "actual_sha256": actual,
            "actual_length": len(actual),
        }
    if actual != expected:
        return "BLOCKED", [f"R4.1 SHA mismatch: expected {expected}, actual {actual}"], {
            "actual_sha256": actual,
            "reference_sha256": expected,
            "artifact_size_bytes": artifact.stat().st_size,
        }

    required_lock_evidence = {
        "sha256_before": actual,
        "sha256_after": actual,
        "sha_match": True,
        "solid_count": 62,
        "solid_count_expected": 62,
        "solid_count_pass": True,
        "step_read_status": "IFSelect_ReturnStatus.IFSelect_RetDone",
        "brep_valid": True,
        "mutation_detected": False,
    }
    missing = [key for key in required_lock_evidence if key not in lock]
    if missing:
        return "BLOCKED", [f"baseline-lock evidence missing required fields: {missing}"], {
            "actual_sha256": actual
        }
    for key, expected_value in required_lock_evidence.items():
        if lock.get(key) != expected_value:
            blockers.append(f"baseline-lock evidence mismatch for {key}: {lock.get(key)!r}")

    step_status, step_reason, solid_count = _native_step_check(artifact)
    if step_status != "PASS":
        blockers.append(step_reason or f"native STEP check={step_status}")
    if solid_count != 62:
        blockers.append(f"native solid count={solid_count}; required=62")

    info = {
        "artifact": str(artifact),
        "artifact_size_bytes": artifact.stat().st_size,
        "actual_sha256": actual,
        "actual_length": len(actual),
        "reference_sha256": expected,
        "reference_length": len(expected),
        "step_read": step_status,
        "brep": step_status,
        "solid_count": solid_count,
        "lock_mutation_detected": lock.get("mutation_detected"),
    }
    if blockers:
        return "FAIL", blockers, info
    return "PASS", [], info


def main() -> int:
    blockers: list[str] = []
    g0, g0_reasons, g0_info = baseline_gate()
    blockers.extend(f"G0 {g0}: {reason}" for reason in g0_reasons)

    environment, env_obj = read_gate_record(AEGIS_ROOT / "ENVIRONMENT" / "environment_manifest.json")
    if environment != "PASS":
        blockers.append(
            f"ENVIRONMENT {environment}: "
            f"{env_obj.get('stop_condition', 'environment qualification not closed')}"
        )

    g1 = "NOT_PROVEN"
    blockers.append("G1 NOT_PROVEN: no separate G1 acceptance contract was supplied by the work order.")

    g2, g2_obj = read_gate_record(AEGIS_ROOT / "G2" / "RUNS" / "LATEST_G2_STATUS.json")
    g3, g3_obj = read_gate_record(AEGIS_ROOT / "G3" / "RUNS" / "LATEST_G3_STATUS.json")
    g4, g4_obj = read_gate_record(AEGIS_ROOT / "G4" / "EVIDENCE_ENGINE" / "self_check_result.json")
    if g4_obj.get("self_check") == "PASS":
        g4 = g4_obj.get("gate_status", "NOT_PROVEN")
        if g4 not in VALID_STATES:
            g4 = "NOT_PROVEN"

    for name, state, obj in [("G2", g2, g2_obj), ("G3", g3, g3_obj), ("G4", g4, g4_obj)]:
        if state != "PASS":
            blockers.append(f"{name} {state}: {obj.get('reason', 'gate not proven')}")

    states = (environment, g0, g1, g2, g3, g4)
    if "FAIL" in states:
        overall = "FAIL"
    elif "BLOCKED" in states:
        overall = "BLOCKED"
    elif "NOT_PROVEN" in states:
        overall = "NOT_PROVEN"
    else:
        overall = "PASS"

    result = {
        "schema_version": "1.0",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "fail_closed": True,
        "ENVIRONMENT": environment,
        "G0": g0,
        "G1": g1,
        "G2": g2,
        "G3": g3,
        "G4": g4,
        "OVERALL": overall,
        "claim_boundary": {
            "CRASH_VALIDATION": "NOT_PROVEN",
            "SAFETY_VALIDATION": "NOT_PROVEN",
            "MANUFACTURING_VALIDATION": "NOT_PROVEN",
        },
        "baseline": g0_info,
        "blockers": blockers,
    }
    out = AEGIS_ROOT / "VALIDATION" / "FINAL_GATE_RESULT.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"ENVIRONMENT = {environment}")
    print(f"G0 = {g0}")
    print(f"G1 = {g1}")
    print(f"G2 = {g2}")
    print(f"G3 = {g3}")
    print(f"G4 = {g4}")
    print(f"OVERALL = {overall}")
    for reason in blockers:
        print(f"BLOCKER = {reason}")
    return 0 if overall == "PASS" else 1
