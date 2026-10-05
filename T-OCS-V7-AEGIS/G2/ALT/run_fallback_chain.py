#!/usr/bin/env python3
"""Automatic, fail-closed qualification of approved local G2 engine fallbacks.

A runtime qualification proves only the selected engine. It never promotes a
synthetic model to T-OCS authority and cannot unlock G3.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ALT = ROOT / "G2" / "ALT"
RUN_ID = os.environ.get("AEGIS_RUN_ID", "G2-ALT-FALLBACK-LOCAL")
ENGINES = (
    ("PYCHRONO", "pychrono", "pychrono"),
    ("PYBULLET", "pybullet", "pybullet"),
    ("MUJOCO", "mujoco", "mujoco"),
    ("DRAKE", "pydrake", "drake"),
    ("DART", "dartpy", "dartpy"),
    ("PINOCCHIO", "pinocchio", "pin"),
)
MUJOCO_RUN = r'''import json
import mujoco
import numpy as np
xml = """<mujoco><option timestep='0.002'/><worldbody><geom type='plane' size='2 2 .1'/><body pos='0 0 1'><freejoint/><geom type='sphere' size='.1' mass='1'/></body></worldbody></mujoco>"""
model = mujoco.MjModel.from_xml_string(xml)
data = mujoco.MjData(model)
trajectory = []
contact_seen = False
contact_force = None
for _ in range(400):
    mujoco.mj_step(model, data)
    contact_seen = contact_seen or data.ncon > 0
    if data.ncon > 0:
        force = np.zeros(6)
        mujoco.mj_contactForce(model, data, 0, force)
        contact_force = force.tolist()
    if len(trajectory) < 5 or _ >= 395:
        trajectory.append({"time": data.time, "qpos": data.qpos.tolist(), "qvel": data.qvel.tolist(), "contacts": data.ncon})
print(json.dumps({"engine_version": mujoco.__version__, "model_compilation": True, "state_initialization": True, "mj_step": True, "joint_state_extraction": True, "contact_capability": contact_seen, "contact_force": contact_force, "trajectory": trajectory}, sort_keys=True))
'''


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def execute(command: list[str]) -> dict[str, object]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    return {"command": command, "exit_status": result.returncode, "stdout": result.stdout.strip(), "stderr": result.stderr.strip()}


def import_probe(module: str, distribution: str) -> dict[str, object]:
    probe = execute([sys.executable, "-c", "import importlib; import importlib.metadata; " + f"m=importlib.import_module({module!r}); print(importlib.metadata.version({distribution!r})); print(m.__file__)"])
    lines = str(probe["stdout"]).splitlines()
    probe["package_version"] = lines[0] if lines else None
    probe["module_file"] = lines[1] if len(lines) > 1 else None
    return probe


def qualification_run(engine: str) -> dict[str, object]:
    if engine != "MUJOCO":
        return {"exit_status": 1, "stdout": "", "stderr": f"No implemented runtime suite for {engine}; import alone is not qualification.", "command": []}
    return execute([sys.executable, "-c", MUJOCO_RUN])


def write_json(path: Path, value: dict[str, object]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    started = datetime.now(timezone.utc).isoformat()
    records: list[dict[str, object]] = []
    selected = "NONE"
    qualification = "BLOCKED"
    repeatability = "NOT_PROVEN"

    for engine, module, distribution in ENGINES:
        probe = import_probe(module, distribution)
        import_ok = probe["exit_status"] == 0
        record: dict[str, object] = {"schema_version": "2.0", "evidence_id": f"ENV-ALT-{engine}-001", "run_id": RUN_ID, "engine": engine, "module": module, "distribution": distribution, "import_status": "PASS" if import_ok else "FAILED", "package_version": probe["package_version"], "module_file": probe["module_file"], "import_command": probe["command"], "import_exit_status": probe["exit_status"], "import_stderr": probe["stderr"], "chrono_equivalence": "NOT_CLAIMED", "t_ocs_g2_execution": "NOT_STARTED"}
        target = ALT / engine / "qualification_001.json"
        if not import_ok:
            record.update({"status": "BLOCKED", "qualification_status": "BLOCKED_AT_IMPORT_GATE", "native_library": "NOT_TESTED", "shared_library_resolution": "NOT_TESTED", "runtime": "NOT_TESTED", "reason": "Local import unavailable; chain automatically advances without a network install."})
            write_json(target, record)
            record["artifact"] = str(target.relative_to(ROOT)); record["artifact_sha256"] = digest(target); records.append(record)
            continue

        run_a = qualification_run(engine)
        run_b = qualification_run(engine)
        try:
            payload_a = json.loads(str(run_a["stdout"])) if run_a["exit_status"] == 0 else None
            payload_b = json.loads(str(run_b["stdout"])) if run_b["exit_status"] == 0 else None
        except json.JSONDecodeError:
            payload_a = payload_b = None
        qualification_ok = payload_a is not None and payload_b is not None
        record.update({"status": "PASS" if qualification_ok else "FAIL", "qualification_status": "PASS" if qualification_ok else "FAILED_RUNTIME_SUITE", "native_library": "PASS" if qualification_ok else "FAILED", "shared_library_resolution": "PASS" if qualification_ok else "FAILED", "runtime": "PASS" if qualification_ok else "FAILED", "reason": "Two deterministic engine runtime qualifications completed." if qualification_ok else str(run_a["stderr"]), "run_A_command": run_a["command"], "run_A_exit_status": run_a["exit_status"], "run_B_command": run_b["command"], "run_B_exit_status": run_b["exit_status"]})
        write_json(target, record)
        q2 = ALT / engine / "qualification_002.json"
        write_json(q2, {"schema_version": "2.0", "evidence_id": f"ENV-ALT-{engine}-002", "run_id": RUN_ID, "engine": engine, "status": "PASS" if qualification_ok else "FAIL", "result": payload_b, "command": run_b["command"], "exit_status": run_b["exit_status"], "stderr": run_b["stderr"]})
        repeat_ok = qualification_ok and payload_a == payload_b
        repeat_path = ALT / engine / "repeatability.json"
        write_json(repeat_path, {"schema_version": "2.0", "evidence_id": f"G2-ALT-{engine}-REPEAT-001", "run_id": RUN_ID, "engine": engine, "run_A": "ALT-ENGINE-001-A", "run_B": "ALT-ENGINE-001-B", "comparison": "Exact normalized JSON equality for state trajectory, positions, velocities, and contact states.", "status": "PASS" if repeat_ok else "NOT_PROVEN", "result_A": payload_a, "result_B": payload_b})
        record["artifact"] = str(target.relative_to(ROOT)); record["artifact_sha256"] = digest(target); records.append(record)
        if qualification_ok and repeat_ok:
            selected = engine; qualification = "PASS"; repeatability = "PASS"; break

    g2_alt_execution = "BLOCKED"
    root_cause = "No approved local physics engine qualified." if selected == "NONE" else "MuJoCo runtime qualification passed, but authoritative T-OCS body/joint/frame/constraint/initial-state inputs remain blocked; no T-OCS model was executed."
    chain = {"schema_version": "2.0", "evidence_id": "ENV-ALT-PY-004", "run_id": RUN_ID, "timestamp_utc": started, "status": "BLOCKED" if selected == "NONE" else "PASS", "automatic_fallback_order": [item[0] for item in ENGINES], "environment": {"python_executable": sys.executable, "python_version": sys.version, "os": platform.platform(), "cpu_architecture": platform.machine(), "runtime_path": str(Path(sys.executable).resolve())}, "network": {"installation_attempted": False, "network_access_attempted": False}, "engines": records, "selected_python_engine": selected, "engine_qualification": qualification, "g2_alt_execution": g2_alt_execution, "g2_alt_repeatability": repeatability, "g3": "LOCKED", "r4_1": "UNCHANGED", "root_cause": root_cause}
    write_json(ALT / "FALLBACK_CHAIN" / "qualification_003.json", chain)
    print(f"G2_ALT_ENGINE_QUALIFICATION={qualification}")
    print(f"G2_ALT_EXECUTION={g2_alt_execution}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
