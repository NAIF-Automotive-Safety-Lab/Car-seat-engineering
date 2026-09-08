#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import platform
import sys
import time
import traceback
from pathlib import Path


def vector_value(obj):
    try:
        return obj.GetXYZ() if hasattr(obj, "GetXYZ") else [obj.x, obj.y, obj.z]
    except Exception:
        return str(obj)


def main() -> int:
    out = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "engine": "Project Chrono / PyChrono",
        "host": platform.platform(),
        "arch": platform.machine(),
        "python": platform.python_version(),
        "expected_source_commit": os.environ.get("PYCHRONO_EXPECTED_SOURCE_COMMIT"),
    }
    try:
        import pychrono as chrono
        import pychrono.core as core

        out["package_file"] = str(Path(chrono.__file__).resolve())
        out["core_module_file"] = str(Path(core.__file__).resolve())
        out["compiled_core_present"] = Path(core.__file__).with_name("_core.so").is_file()
        out["import"] = "PASS"
        out["version"] = getattr(chrono, "__version__", "UNKNOWN")
        out["core_class_present"] = hasattr(core, "ChSystemNSC")
        if not out["compiled_core_present"] or not out["core_class_present"]:
            raise RuntimeError("REAL_PYCHRONO_CORE_BINDING_NOT_VERIFIED")

        system = core.ChSystemNSC()
        out["basic_system_creation"] = True
        fixed = core.ChBody()
        fixed.SetFixed(True)
        fixed.SetMass(1)
        system.AddBody(fixed)
        body = core.ChBody()
        body.SetMass(1)
        body.SetPos(core.ChVector3d(1, 0, 0))
        system.AddBody(body)
        out["basic_body_creation"] = True

        joint = core.ChLinkLockRevolute()
        joint.Initialize(fixed, body, core.ChFramed(core.ChVector3d(0, 0, 0)))
        system.AddLink(joint)
        out["basic_joint_creation"] = True

        for _ in range(10):
            system.DoStepDynamics(1e-3)
        out["basic_solver_execution"] = system.GetChTime() > 0
        out["simulation_time"] = system.GetChTime()
        out["position"] = vector_value(body.GetPos())
        out["velocity"] = vector_value(body.GetPosDt())
        reaction = joint.GetReaction1()
        out["reaction_extracted"] = reaction is not None
        out["reaction_type"] = type(reaction).__name__
        out["status"] = "PASS"
    except Exception as exc:
        out["status"] = "FAIL"
        out["error"] = f"{type(exc).__name__}: {exc}"
        out["traceback"] = traceback.format_exc()

    artifact = Path("artifacts/pychrono_smoke.json")
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print("PYCHRONO_IMPORT=" + out.get("import", "FAIL"))
    print("PYCHRONO_COMPILED_CORE=" + ("PASS" if out.get("compiled_core_present") else "FAIL"))
    print("PYCHRONO_SYSTEM=" + ("PASS" if out.get("basic_system_creation") else "FAIL"))
    print("PYCHRONO_SOLVER=" + ("PASS" if out.get("basic_solver_execution") else "FAIL"))
    print("PYCHRONO_SMOKE=" + out["status"])
    return 0 if out["status"] == "PASS" else 7


if __name__ == "__main__":
    sys.exit(main())
