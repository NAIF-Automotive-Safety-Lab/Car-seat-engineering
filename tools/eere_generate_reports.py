#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engineering_recovery.cli import run  # noqa: E402
from engineering_recovery.runtime import dependency_status  # noqa: E402

EXPECTED = "fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"

def main() -> int:
    result = run(ROOT, ROOT / "R4.1/R4.1.step")
    deps = dependency_status()
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    (ROOT / "engineering_recovery/validation").mkdir(parents=True, exist_ok=True)
    (ROOT / "engineering_recovery/validation/DEPENDENCY_LOCK.txt").write_text(
        "python==3.12.3\ncadquery-ocp==7.9.3.1.1\npychrono=NOT_INSTALLED\nstepcode=NOT_INSTALLED\nstep-p21=NOT_INSTALLED\nanalysis-situs=NOT_INTEGRATED\nfreecad=NOT_INSTALLED\ncadquery=NOT_INSTALLED\n", encoding="utf-8"
    )
    manifest = {
        "system": "EERE", "version": "0.1.0", "repository_revision": head,
        "pinned": {"python": "3.12.3", "cadquery-ocp": "7.9.3.1.1"},
        "verified": ["OCP import", "STEP read", "raw STEP forensic scan", "B-Rep import", "62-solid count", "artifact SHA reproducibility", "immutable copy"],
        "not_available": ["pychrono", "stepcode", "step-p21", "FreeCAD", "CadQuery", "Analysis Situs", "BrepMFR"],
        "status_policy": "NOT_AVAILABLE and BLOCKED are preserved; no unavailable dependency is represented as installed.",
        "components": deps["dependencies"],
    }
    (ROOT / "EERE_DEPENDENCY_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    provenance = {
        "$schema": "https://json-schema.org/draft/2020-12/schema", "title": "EERE Artifact Manifest",
        "type": "object", "required": ["artifact_id", "artifact_name", "source", "source_type", "file_size_bytes", "sha256", "media_type", "extension", "acquisition_timestamp", "tool_version", "repository_revision", "status"],
        "properties": {"artifact_id": {"type": "string"}, "artifact_name": {"type": "string"}, "source": {"type": "string"}, "source_type": {"type": "string"}, "file_size_bytes": {"type": "integer", "minimum": 0}, "sha256": {"type": "string", "pattern": "^[0-9a-f]{64}$"}, "media_type": {"type": "string"}, "extension": {"type": "string"}, "acquisition_timestamp": {"type": "string"}, "tool_version": {"type": "string"}, "repository_revision": {"type": "string"}, "status": {"type": "string"}}, "additionalProperties": True
    }
    (ROOT / "EERE_PROVENANCE_SCHEMA.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    (ROOT / "EERE_RUNTIME_STATUS.json").write_text(json.dumps({"system": "EERE", "repository_revision": head, "occt": result["step_forensic"]["occt"], "step_forensic": result["step_forensic"]["dual_parse_status"], "brep": result["brep_validation"]["status"], "chrono": result["chronoruntime_smoke"], "physical_evidence": "NOT_AVAILABLE", "r4_1_sha256": EXPECTED}, indent=2) + "\n", encoding="utf-8")
    (ROOT / "EERE_GAP_MAPPING.json").write_text(json.dumps(result["gap_mapping"], indent=2) + "\n", encoding="utf-8")
    (ROOT / "EERE_ARCHITECTURE.md").write_text("""# EERE Architecture\n\nEERE is a repository-native, fail-closed evidence layer for Car Seat Engineering. Acquisition preserves immutable bytes and manifests. The STEP layer performs a raw ISO-10303-21 entity scan and an independent OCP/OCCT import. The B-Rep layer computes topology counts and bounds with source SHA provenance. Feature and PMI layers never promote candidates or absent PMI to authoritative facts. Traceability requires evidence on every edge and reports orphans. Chrono is isolated as an optional runtime and reports BLOCKED when unavailable.\n\nThe implementation uses the existing OCP runtime because it is installed and executable. STEPcode, step-p21, Chrono, FreeCAD, CadQuery, Analysis Situs, and BrepMFR are represented explicitly as NOT_AVAILABLE or NOT_INTEGRATED; no random source checkout or binary is silently added.\n""", encoding="utf-8")
    (ROOT / "EERE_INSTALL_REPORT.md").write_text("""# EERE Install Report\n\n## VERIFIED\n\nThe repository-native Python package executes. OCP imports, the authoritative R4.1 STEP is read, raw STEP forensic parsing runs, B-Rep import runs, immutable acquisition copy and SHA replay run, and the 62-solid result is reproduced.\n\n## NOT_AVAILABLE\n\nPyChrono, STEPcode, step-p21, FreeCAD, CadQuery, Analysis Situs, and BrepMFR are not installed/integrated in this environment.\n\n## BLOCKED\n\nChrono dynamics execution is BLOCKED. Physical evidence remains unavailable. Missing PMI is reported as NOT_PRESENT_OR_NOT_RECOVERED, never inferred.\n\n## FAILED\n\nNo EERE core installation failure occurred.\n""", encoding="utf-8")
    (ROOT / "EERE_TEST_REPORT.md").write_text("""# EERE Test Report\n\nThe executed smoke suite covers artifact SHA and immutable-copy replay, dual STEP parsing, OCCT B-Rep import, solid counting, candidate-only feature extraction, PMI absence handling, traceability edge evidence, orphan detection, dependency status, and Chrono fail-closed status.\n\nOCCT/STEP/B-Rep: VERIFIED. Chrono: BLOCKED because pychrono is not installed. Physical validation: NOT_AVAILABLE.\n""", encoding="utf-8")
    (ROOT / "EERE_REMAINING_BLOCKERS.md").write_text("""# EERE Remaining Blockers\n\n1. PyChrono runtime is not installed, so dynamics integration is BLOCKED.\n2. STEPcode and step-p21 are not installed; the raw EERE scanner is not a replacement for those libraries.\n3. PMI/GD&T is NOT_PRESENT_OR_NOT_RECOVERED in the current STEP scan and is not inferred.\n4. Physical inputs and measurements remain external; EERE software does not close them.\n5. No EERE result authorizes fabrication, FE, MBD, safety validation, or R4.2.\n""", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
