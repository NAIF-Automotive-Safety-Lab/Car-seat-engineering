from __future__ import annotations
import json
from datetime import datetime, timezone
from pathlib import Path
from .provenance import sha256_file

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT.parent / "R4.1" / "R4.1.step"
EXPECTED_R41_SHA256 = "fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"
VALID_STATES = {"PASS","FAIL","BLOCKED","NOT_PROVEN"}

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

def baseline_gate() -> tuple[str, list[str], dict]:
    if not BASELINE.is_file():
        return "BLOCKED", ["R4.1/R4.1.step is absent from execution checkout"], {}
    actual = sha256_file(BASELINE)
    if len(EXPECTED_R41_SHA256) != 64:
        return "BLOCKED", [
            "R4.1 reference SHA is malformed: SHA-256 must contain exactly 64 hexadecimal characters",
            f"supplied_reference_length={len(EXPECTED_R41_SHA256)}",
            f"computed_actual_sha256={actual}"
        ], {"actual_sha256":actual,"actual_length":len(actual),
            "reference_sha256":EXPECTED_R41_SHA256,"reference_length":len(EXPECTED_R41_SHA256)}
    try:
        int(EXPECTED_R41_SHA256, 16)
    except ValueError:
        return "BLOCKED", ["R4.1 reference SHA contains non-hexadecimal characters"], {"actual_sha256":actual}
    if actual != EXPECTED_R41_SHA256:
        return "BLOCKED", [f"R4.1 SHA mismatch: expected {EXPECTED_R41_SHA256}, actual {actual}"], {"actual_sha256":actual}
    return "PASS", [], {"actual_sha256":actual}

def main() -> int:
    blockers=[]
    g0,g0_reasons,g0_info=baseline_gate()
    for reason in g0_reasons:
        blockers.append(f"G0 {g0}: {reason}")

    environment,env_obj=read_gate_record(ROOT/"ENVIRONMENT"/"environment_manifest.json")
    if environment != "PASS":
        blockers.append(f"ENVIRONMENT {environment}: {env_obj.get('stop_condition','environment qualification not closed')}")

    g1="NOT_PROVEN"
    blockers.append("G1 NOT_PROVEN: no separate G1 acceptance contract was supplied by the work order.")

    g2,g2_obj=read_gate_record(ROOT/"G2"/"RUNS"/"LATEST_G2_STATUS.json")
    g3,g3_obj=read_gate_record(ROOT/"G3"/"RUNS"/"LATEST_G3_STATUS.json")
    g4,g4_obj=read_gate_record(ROOT/"G4"/"EVIDENCE_ENGINE"/"self_check_result.json")
    if g4_obj.get("self_check")=="PASS":
        g4=g4_obj.get("gate_status","NOT_PROVEN")
        if g4 not in VALID_STATES:
            g4="NOT_PROVEN"

    for name,state,obj in [("G2",g2,g2_obj),("G3",g3,g3_obj),("G4",g4,g4_obj)]:
        if state != "PASS":
            blockers.append(f"{name} {state}: {obj.get('reason','gate not proven')}")

    states=(environment,g0,g1,g2,g3,g4)
    if "FAIL" in states:
        overall="FAIL"
    elif "BLOCKED" in states:
        overall="BLOCKED"
    elif "NOT_PROVEN" in states:
        overall="NOT_PROVEN"
    else:
        overall="PASS"

    result={
        "schema_version":"1.0",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "fail_closed":True,
        "ENVIRONMENT":environment,
        "G0":g0,"G1":g1,"G2":g2,"G3":g3,"G4":g4,
        "OVERALL":overall,
        "claim_boundary":{
            "CRASH_VALIDATION":"NOT_PROVEN",
            "SAFETY_VALIDATION":"NOT_PROVEN",
            "MANUFACTURING_VALIDATION":"NOT_PROVEN"
        },
        "baseline":g0_info,
        "blockers":blockers
    }
    out=ROOT/"VALIDATION"/"FINAL_GATE_RESULT.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+"
",encoding="utf-8")
    print(f"ENVIRONMENT = {environment}")
    print(f"G0 = {g0}")
    print(f"G1 = {g1}")
    print(f"G2 = {g2}")
    print(f"G3 = {g3}")
    print(f"G4 = {g4}")
    print(f"OVERALL = {overall}")
    for reason in blockers:
        print(f"BLOCKER = {reason}")
    return 0 if overall=="PASS" else 1
