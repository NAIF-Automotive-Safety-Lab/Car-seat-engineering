from __future__ import annotations

import importlib.util
import platform
import sys
import traceback
from typing import Any


def dependency_status() -> dict[str, Any]:
    deps = {
        "occt": {"module": "OCP", "role": "STEP/B-Rep/XDE"},
        "analysis_situs": {"module": None, "role": "B-Rep analysis/feature recognition"},
        "stepcode": {"module": "stepcode", "role": "STEP semantic parser"},
        "step_p21": {"module": "step_p21", "role": "ISO-10303-21 parser"},
        "chrono": {"module": "pychrono", "role": "dynamics"},
        "freecad": {"module": "FreeCAD", "role": "optional CAD interoperability"},
        "cadquery": {"module": "cadquery", "role": "Python CAD automation"},
        "brepmfr": {"module": None, "role": "optional candidate-only AI feature recognition"},
    }
    for item in deps.values():
        module = item["module"]
        item["status"] = "NOT_APPLICABLE" if module is None and item["role"].startswith("optional") else ("VERIFIED" if module and importlib.util.find_spec(module) else "NOT_AVAILABLE")
    deps["occt"]["version"] = "OCP-installed; exact OCCT runtime version not exposed by binding"
    return {"python": sys.version, "platform": platform.platform(), "dependencies": deps}


def chrono_smoke() -> dict[str, Any]:
    result: dict[str, Any] = {"runtime_available": False, "status": "BLOCKED", "platform": platform.platform(), "python": sys.version}
    try:
        import pychrono.core as chrono
        result["runtime_available"] = True
        result["version"] = getattr(chrono, "__version__", "UNKNOWN")
        system = chrono.ChSystemNSC()
        result["basic_system_creation"] = True
        body = chrono.ChBody(); body.SetFixed(True); body.SetMass(1.0); system.AddBody(body)
        result["basic_body_creation"] = True
        result["basic_joint_creation"] = False
        result["basic_solver_execution"] = False
        result["status"] = "VERIFIED"
    except Exception as exc:
        result.update({"error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc(), "status": "BLOCKED"})
    return result
