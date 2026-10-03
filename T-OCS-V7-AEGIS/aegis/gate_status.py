from __future__ import annotations
import json
from datetime import datetime, timezone
from pathlib import Path
from .provenance import sha256_file

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT.parent / "R4.1" / "R4.1.step"
EXPECTED_R41_SHA256 = "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"
VALID_STATES = {"PASS","FAIL","BLOCKED","NOT_PROVEN"}

def load_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return value if isinstance(value, dict) else {}

def read_state(path: Path) -> str:
    value = load_json(path).get("status", "NOT_PROVEN")
    return value if value in VALID_STATES else "NOT_PROVEN"

def main() -> int:
    blockers = []
    if not BASELINE.is_file():
        g0 = "BLOCKED"
        blockers.append("R4.1/R4.1.step is absent from the execution checkout")
    else:
        actual = sha256_file(BASELINE)
        if actual == EXPECTED_R41_SHA256:
            g0 = "PASS"
        else:
            g0 = "FAIL"
            blockers.append(f"R4.1 SHA-256 mismatch: expected {EXPECTED_R41_SHA256}, actual {actual}")

    g1 = "NOT_PROVEN"
    blockers.append("No separate G1 acceptance contract was supplied by the work order.")
    g2 = read_state(ROOT/"G2"/"RUNS"/"LATEST_G2_STATUS.json")
    g3 = read_state(ROOT/"G3"/"RUNS"/"LATEST_G3_STATUS.json")
    g4 = read_state(ROOT/"G4"/"EVIDENCE_ENGINE"/"self_check_result.json")

    for name,state in [("G2",g2),("G3",g3),("G4",g4)]:
        if state == "NOT_PROVEN":
            blockers.append(f"{name} status is absent or not proven")

    states = (g0,g2,g3,g4)
    if any(s == "FAIL" for s in states):
        overall = "FAIL"
    elif any(s == "BLOCKED" for s in states):
        overall = "BLOCKED"
    elif any(s == "NOT_PROVEN" for s in (g0,g1,g2,g3,g4)):
        overall = "NOT_PROVEN"
    else:
        overall = "PASS"

    result = {
        "schema_version":"1.0",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "fail_closed":True,
        "G0":g0,"G1":g1,"G2":g2,"G3":g3,"G4":g4,
        "OVERALL":overall,
        "claim_boundary":{
            "CRASH_VALIDATION":"NOT_PROVEN",
            "SAFETY_VALIDATION":"NOT_PROVEN",
            "MANUFACTURING_VALIDATION":"NOT_PROVEN"
        },
        "blockers":blockers
    }
    out = ROOT/"VALIDATION"/"FINAL_GATE_RESULT.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(f"G0 = {g0}")
    print(f"G1 = {g1}")
    print(f"G2 = {g2}")
    print(f"G3 = {g3}")
    print(f"G4 = {g4}")
    print(f"OVERALL = {overall}")
    for b in blockers:
        print(f"BLOCKER = {b}")
    return 0 if overall == "PASS" else 1
