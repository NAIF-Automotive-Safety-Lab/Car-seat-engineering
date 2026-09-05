from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

EXPECTED_SHA256 = "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"
EXPECTED_SOLIDS = 62
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SR11_ROOT = PROJECT_ROOT / "SR11"
STEP_PATH = PROJECT_ROOT / "R4.1" / "R4.1.step"
BODY_MAP = SR11_ROOT / "manifests" / "V1_SR11_BODY_MAPPING.json"
JOINT_MAP = SR11_ROOT / "joints" / "V1_SR11_JOINT_MAPPING.json"
CONTACT_MAP = SR11_ROOT / "contacts" / "V1_SR11_CONTACT_MAPPING.json"
ASSUMPTIONS = SR11_ROOT / "parameters" / "V1_SR11_ASSUMPTIONS.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(str(path))
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON_ROOT_NOT_OBJECT:{path}")
    return value


def _validate_mapping(manifest: dict) -> list[dict]:
    checks: list[dict] = []
    bodies = manifest.get("bodies")
    if not isinstance(bodies, list) or not bodies:
        raise ValueError("BODY_MAPPING_BODIES_NOT_NONEMPTY_LIST")
    if manifest.get("solids") != EXPECTED_SOLIDS:
        raise ValueError("BODY_MAPPING_SOLIDS_DECLARATION_MISMATCH")
    ids: list[str] = []
    solid_indices: list[int] = []
    for index, body in enumerate(bodies):
        if not isinstance(body, dict):
            raise ValueError(f"BODY_RECORD_NOT_OBJECT:{index}")
        body_id = body.get("mbd_body_id")
        if not isinstance(body_id, str) or not body_id.strip():
            raise ValueError(f"EMPTY_MBD_BODY_ID:{index}")
        if "step_solid_indices" not in body or not isinstance(body["step_solid_indices"], list):
            raise ValueError(f"MISSING_STEP_SOLID_INDICES:{body_id}")
        ids.append(body_id.strip())
        for solid in body["step_solid_indices"]:
            if not isinstance(solid, int) or isinstance(solid, bool):
                raise ValueError(f"NON_INTEGER_STEP_SOLID_INDEX:{body_id}")
            solid_indices.append(solid)
    if len(ids) != len(set(ids)):
        raise ValueError("DUPLICATE_MBD_BODY_ID")
    expected_indices = list(range(1, EXPECTED_SOLIDS + 1))
    if sorted(solid_indices) != expected_indices or len(set(solid_indices)) != EXPECTED_SOLIDS:
        raise ValueError("STEP_SOLID_MAPPING_NOT_UNIQUE_1_TO_62")
    checks.append({"check": "body_records", "actual": len(bodies), "pass": True})
    checks.append({"check": "solid_mapping", "expected": "62/62", "actual": "62/62", "pass": True})
    return checks


def _validate_references(joint_map: dict, contact_map: dict, assumptions: dict) -> list[dict]:
    checks: list[dict] = []
    joints = joint_map.get("joints")
    contacts = contact_map.get("contacts")
    parameters = assumptions.get("parameters")
    if not isinstance(joints, list) or not joints:
        raise ValueError("JOINT_MAPPING_NOT_NONEMPTY_LIST")
    if not isinstance(contacts, list) or not contacts:
        raise ValueError("CONTACT_MAPPING_NOT_NONEMPTY_LIST")
    if not isinstance(parameters, list) or not parameters:
        raise ValueError("ASSUMPTIONS_NOT_NONEMPTY_LIST")
    if any(not isinstance(j, dict) or not j.get("body_A") or not j.get("body_B") for j in joints):
        raise ValueError("JOINT_REFERENCE_MISSING")
    if any(not isinstance(c, dict) or not c.get("body_A") or not c.get("body_B") for c in contacts):
        raise ValueError("CONTACT_REFERENCE_MISSING")
    allowed = {"UNKNOWN", "ASSUMED", "SOURCE_VERIFIED", "DERIVED"}
    if any(not isinstance(p, dict) or p.get("status") not in allowed for p in parameters):
        raise ValueError("PARAMETER_STATUS_INVALID")
    checks.extend([
        {"check": "joint_references", "count": len(joints), "pass": True},
        {"check": "contact_references", "count": len(contacts), "pass": True},
        {"check": "no_silent_defaults", "count": len(parameters), "pass": True},
    ])
    return checks


def static_validate() -> dict:
    checks: list[dict] = []
    if not STEP_PATH.is_file():
        return {"status": "FAIL", "failure": "STEP_MISSING", "path": str(STEP_PATH), "checks": checks}
    actual = sha256(STEP_PATH)
    checks.append({"check": "sha256", "expected": EXPECTED_SHA256, "actual": actual, "pass": actual == EXPECTED_SHA256})
    if actual != EXPECTED_SHA256:
        return {"status": "FAIL", "failure": "R4.1_SHA_MISMATCH", "checks": checks}
    try:
        from OCP.STEPControl import STEPControl_Reader
        from OCP.TopAbs import TopAbs_SOLID
        from OCP.TopExp import TopExp_Explorer
    except Exception as exc:
        return {"status": "BLOCKED", "failure": "OCP_UNAVAILABLE", "error": repr(exc), "checks": checks}
    reader = STEPControl_Reader()
    read_status = reader.ReadFile(str(STEP_PATH))
    read_pass = "RetDone" in str(read_status)
    if not read_pass:
        return {"status": "FAIL", "failure": "STEP_READ_FAILED", "read_status": str(read_status), "checks": checks}
    reader.TransferRoots()
    shape = reader.OneShape()
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    solid_count = 0
    while explorer.More():
        solid_count += 1
        explorer.Next()
    checks.append({"check": "step_read", "actual_status": str(read_status), "pass": True})
    checks.append({"check": "solid_count", "expected": EXPECTED_SOLIDS, "actual": solid_count, "pass": solid_count == EXPECTED_SOLIDS})
    if solid_count != EXPECTED_SOLIDS:
        return {"status": "FAIL", "failure": "SOLID_COUNT_MISMATCH", "checks": checks}
    try:
        body_map = _load_json(BODY_MAP)
        joint_map = _load_json(JOINT_MAP)
        contact_map = _load_json(CONTACT_MAP)
        assumptions = _load_json(ASSUMPTIONS)
        checks.extend(_validate_mapping(body_map))
        checks.extend(_validate_references(joint_map, contact_map, assumptions))
    except Exception as exc:
        return {"status": "FAIL", "failure": "STATIC_MAPPING_VALIDATION_FAILED", "error": repr(exc), "checks": checks}
    return {
        "status": "PASS",
        "validation_class": "STATIC_ADAPTER_VALIDATION",
        "solver_execution_performed": False,
        "step_path": str(STEP_PATH),
        "body_record_count": len(body_map["bodies"]),
        "solid_coverage": "62/62",
        "checks": checks,
    }


def build_chrono_system():
    try:
        import pychrono.core as chrono
        import pychrono.cascade as cascade
    except Exception as exc:
        raise RuntimeError(f"PYCHRONO_REQUIRED_BUT_UNAVAILABLE: {type(exc).__name__}: {exc}") from exc
    static = static_validate()
    if static["status"] != "PASS":
        raise RuntimeError(f"STATIC_VALIDATION_{static['status']}:{static.get('failure')}")
    assumptions = _load_json(ASSUMPTIONS)
    unknown = [p["parameter"] for p in assumptions["parameters"] if p["status"] == "UNKNOWN"]
    if unknown:
        raise RuntimeError("MBD_BUILD_BLOCKED_BY_UNRESOLVED_PARAMETERS:" + ",".join(unknown))
    document = cascade.ChCascadeDoc()
    if not document.LoadSTEP(str(STEP_PATH)):
        raise RuntimeError("CHRONO_CASCADE_STEP_LOAD_FAILED")
    system = chrono.ChSystemNSC()
    return system, document


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--static-validate", action="store_true")
    parser.add_argument("--build", action="store_true")
    args = parser.parse_args()
    if args.static_validate:
        result = static_validate()
        print(json.dumps(result, indent=2))
        return 0 if result["status"] == "PASS" else 2
    if args.build:
        try:
            system, document = build_chrono_system()
            print(json.dumps({"status": "CHRONO_SYSTEM_CREATED", "system_type": type(system).__name__, "step_loaded": bool(document)}, indent=2))
            return 0
        except Exception as exc:
            print(json.dumps({"status": "BLOCKED", "error": str(exc)}, indent=2))
            return 2
    parser.error("one of --static-validate or --build is required")
    return 2


if __name__ == "__main__":
    sys.exit(main())
