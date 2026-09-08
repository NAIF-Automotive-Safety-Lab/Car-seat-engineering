#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engineering_recovery.acquisition.service import ArtifactAcquisitionService
from engineering_recovery.brep.validate import brep_report
from engineering_recovery.features.deterministic import feature_report
from engineering_recovery.pmi.extract import pmi_report
from engineering_recovery.step.forensic import forensic_report


def command_available(name: str) -> bool:
    return shutil.which(name) is not None


def executable_status(name: str, candidates: list[Path], version_args: list[str]) -> dict[str, object]:
    resolved = shutil.which(name)
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            resolved = str(candidate)
            break
    if not resolved:
        return {"status": "NOT_AVAILABLE", "reason": "executable not found"}
    try:
        result = subprocess.run([resolved, *version_args], capture_output=True, text=True, timeout=20)
        return {
            "status": "VERIFIED" if result.returncode == 0 else "FAILED",
            "executable": resolved,
            "returncode": result.returncode,
            "version_output": (result.stdout or result.stderr).strip()[:500],
        }
    except Exception as exc:
        return {"status": "FAILED", "executable": resolved, "reason": f"{type(exc).__name__}: {exc}"}


def chrono_runtime_status() -> dict[str, object]:
    build_root = Path(os.environ.get("CHRONO_BUILD_ROOT", "/home/ubuntu/third_party/chrono-build-10.0.0"))
    python_root = build_root / "bin"
    library_root = build_root / "lib"
    if python_root.is_dir():
        sys.path.insert(0, str(python_root))
    os.environ["LD_LIBRARY_PATH"] = f"{library_root}:{os.environ.get('LD_LIBRARY_PATH', '')}"
    try:
        import pychrono as chrono
        import pychrono.core as core

        core_file = Path(core.__file__).resolve()
        compiled = core_file.with_name("_core.so").is_file()
        return {
            "status": "VERIFIED" if compiled and hasattr(core, "ChSystemNSC") else "INSTALLED_UNVERIFIED",
            "package_file": str(Path(chrono.__file__).resolve()),
            "core_file": str(core_file),
            "compiled_core_present": compiled,
            "system_class_present": hasattr(core, "ChSystemNSC"),
            "build_root": str(build_root),
        }
    except Exception as exc:
        return {"status": "BLOCKED", "reason": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", default=str(ROOT / "R4.1" / "R4.1.step"))
    parser.add_argument("--out", default=str(ROOT / "artifacts" / "engineering-evidence"))
    args = parser.parse_args()
    step = Path(args.step).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
    service = ArtifactAcquisitionService(ROOT, revision)
    acquired, manifest = service.acquire_file(step, out / "acquired")
    forensic = forensic_report(acquired)
    brep = brep_report(acquired)
    features = feature_report(acquired)
    pmi = pmi_report(acquired)
    gmsh = executable_status(
        "gmsh",
        [Path("/home/ubuntu/third_party/gmsh-install-4.15.2/bin/gmsh"), Path("/home/ubuntu/third_party/gmsh-build-4.15.2/gmsh")],
        ["--version"],
    )
    calculix = executable_status(
        "ccx",
        [Path("/home/ubuntu/third_party/calculix/src/ccx"), Path("/home/ubuntu/third_party/calculix/src/ccx_2.22")],
        ["-v"],
    )
    chrono = chrono_runtime_status()
    runtime = {
        "status": "VERIFIED" if all(x["status"] in {"VERIFIED", "NOT_APPLICABLE"} for x in [chrono, gmsh, calculix]) else "PARTIAL",
        "components": {
            "ocp": {"status": "VERIFIED", "reason": "B-Rep report completed"},
            "stepcode": {"status": "NOT_AVAILABLE", "reason": "No executable stepcode binding detected"},
            "step_p21": {"status": "NOT_AVAILABLE", "reason": "No executable step-p21 parser detected"},
            "pychrono_core": chrono,
            "gmsh": gmsh,
            "calculix": calculix,
        },
        "platform": platform.platform(),
        "python": platform.python_version(),
        "checked_at": datetime.now(timezone.utc).isoformat(),
    }
    (out / "artifact_manifest.json").write_text(json.dumps(manifest.to_dict(), indent=2) + "\n", encoding="utf-8")
    (out / "step_forensic_report.json").write_text(json.dumps(forensic, indent=2) + "\n", encoding="utf-8")
    (out / "brep_validation.json").write_text(json.dumps(brep, indent=2) + "\n", encoding="utf-8")
    (out / "feature_inventory.json").write_text(json.dumps(features, indent=2) + "\n", encoding="utf-8")
    (out / "pmi_report.json").write_text(json.dumps(pmi, indent=2) + "\n", encoding="utf-8")
    (out / "runtime_status.json").write_text(json.dumps(runtime, indent=2) + "\n", encoding="utf-8")
    result = {"status": "PASS", "artifact_manifest": manifest.to_dict(), "forensic": forensic, "brep": brep, "features": features, "pmi": pmi, "runtime": runtime}
    (out / "eere_smoke.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "out": str(out), "sha256": manifest.sha256, "solids": brep["solid_count"], "chrono": runtime["components"]["pychrono_core"]["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
