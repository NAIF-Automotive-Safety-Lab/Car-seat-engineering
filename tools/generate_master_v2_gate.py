#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "master_v2"
R41 = ROOT / "R4.1" / "R4.1.step"
R41_SHA = hashlib.sha256(R41.read_bytes()).hexdigest()
NOW = datetime.now(timezone.utc).isoformat()


def write(name: str, value: object) -> None:
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def sha(path: Path) -> str | None:
    if not path.is_file():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def version(command: list[str]) -> str | None:
    try:
        return subprocess.run(command, capture_output=True, text=True, check=True).stdout.strip().splitlines()[0]
    except Exception:
        return None


def status_for(path: str) -> str:
    return "VERIFIED" if Path(path).exists() else "BLOCKED"


parts = [
    ("100", "vehicle interface/base frame"), ("110L", "left longitudinal load path"),
    ("110R", "right longitudinal load path"), ("120", "seat carriage"),
    ("130", "ride-down module"), ("140", "pelvic-control pan"),
    ("150L", "left rotation-control link"), ("150R", "right rotation-control link"),
    ("160", "seatback frame"), ("170", "multi-state lock"),
    ("180", "rebound-control"), ("190", "restraint interface"),
    ("200", "torso/head guidance"), ("210", "sensor/trigger interface"),
]

write("repository_audit.json", {
    "status": "VERIFIED_PARTIAL",
    "repository": str(ROOT),
    "git_revision": version(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
    "platform": platform.platform(), "python": platform.python_version(),
    "detected": {
        "ocp": status_for("/usr/local/lib/python3.12/dist-packages/OCP"),
        "r41_step": {"status": "VERIFIED", "path": str(R41), "sha256": R41_SHA, "size_bytes": R41.stat().st_size},
        "project_chrono_core": {"status": "VERIFIED", "source_commit": "9faf13dd8f1128dd75ed233a9627027b0422c3f7"},
        "gmsh": {"status": status_for("/home/ubuntu/third_party/gmsh-install-4.15.2/bin/gmsh"), "version": "4.15.2", "binary": "/home/ubuntu/third_party/gmsh-install-4.15.2/bin/gmsh"},
        "calculix": {"status": status_for("/home/ubuntu/third_party/calculix/src/CalculiX"), "binary": "/home/ubuntu/third_party/calculix/src/CalculiX"},
        "stepcode": {"status": status_for("/home/ubuntu/third_party/stepcode-install-0.8.2/bin/exp2cxx"), "version": "0.8.2"},
        "occt_7_9_2": {"status": status_for("/home/ubuntu/third_party/occt-install-7.9.2/lib"), "version": "7.9.2"},
        "pychrono_cascade": {"status": "BLOCKED", "reason": "OCCT build was stopped before installation"},
    },
    "missing_or_unverified": ["V7 native CAD payload", "released V7 BOM", "released materials", "released PMI/GD&T", "tolerances", "physical measurements", "V7-specific CAE inputs"],
    "compatibility_risks": ["R4.1 is not evidence that V7 is geometrically released", "synthetic engine smoke outputs are not project CAE evidence", "no external CAD-AI Engineering OS repository was present locally; its public code was inspected as a reference only"],
})

write("dependency_manifest.json", {
    "status": "PARTIAL",
    "generated_at": NOW,
    "dependencies": [
        {"name": "OpenCascade/OCP", "version": "8.0.1.0.0 Python binding", "status": "VERIFIED", "use": "R4.1 STEP/B-Rep"},
        {"name": "Project Chrono", "version": "10.0.0", "commit": "9faf13dd8f1128dd75ed233a9627027b0422c3f7", "status": "VERIFIED", "use": "core dynamics smoke"},
        {"name": "Gmsh", "version": "4.15.2", "commit": "657c8e915f60405e6cad0c8ec7faf812bfff1a60", "status": "VERIFIED", "use": "synthetic mesh smoke only"},
        {"name": "CalculiX", "commit": "3593817a4a7658eb8af9e9827fbe98f9552fbd44", "status": "VERIFIED", "use": "synthetic solver smoke only"},
        {"name": "STEPcode", "version": "0.8.2", "commit": "5cf8fc1a8907983365b233280f128566956713fa", "status": "VERIFIED_BUILD", "use": "independent STEP semantic layer; project run pending"},
        {"name": "step-p21", "commit": "2baab1148970a4a328429a0a5f1e459c979e77ba", "status": "BLOCKED", "use": "Part-21 parser", "reason": "build stopped by explicit user stop"},
        {"name": "Chrono Cascade", "status": "BLOCKED", "reason": "OCCT 7.9.2 build stopped before install"},
    ],
})

write("v7_forensic_audit.json", {
    "status": "BLOCKED_RELEASE",
    "design_candidate": "V7",
    "baseline": "V5",
    "source_payload": {"status": "BLOCKED", "reason": "No native V7 CAD payload in repository; only R4.1 is available", "r41_sha256": R41_SHA},
    "components": [{"reference": ref, "part": part, "geometry": "BLOCKED", "interfaces": "BLOCKED", "dimensions": "BLOCKED", "thicknesses": "BLOCKED", "materials": "BLOCKED", "fasteners": "BLOCKED", "joints_welds": "BLOCKED", "mounting_points": "BLOCKED", "clearances": "BLOCKED", "mechanisms": "BLOCKED", "assembly_feasibility": "BLOCKED", "service_access": "BLOCKED", "critical_load_paths": "BLOCKED", "next_action": "Provide released V7 native CAD and associated engineering records"} for ref, part in parts],
    "critical_gap_count": len(parts),
    "release_decision": "STOP_RELEASE",
})

write("manufacturability_proof.json", {
    "status": "BLOCKED_RELEASE",
    "classification_policy": ["PROVEN", "ASSUMED", "UNVERIFIED", "BLOCKED"],
    "components": [{"part_reference": ref, "part": part, "material": {"status": "BLOCKED", "value": None}, "manufacturing_process": {"status": "BLOCKED", "value": None}, "joining_method": {"status": "BLOCKED", "value": None}, "critical_dimensions": {"status": "BLOCKED", "value": None}, "tolerances": {"status": "BLOCKED", "value": None}, "assembly_method": {"status": "BLOCKED", "value": None}, "manufacturing_risk": {"status": "BLOCKED", "value": "Cannot assess without released V7 definition"}, "evidence": [], "release_status": "BLOCKED"} for ref, part in parts],
    "critical_manufacturability_gap_count": len(parts),
})

write("step_forensic_report.json", {"status": "PARTIAL", "source_sha256": R41_SHA, "source": str(R41), "available_evidence": "artifacts/engineering-evidence/step_forensic_report.json", "independent_semantic_parser": "BLOCKED", "parser_errors": [], "parser_warnings": ["V7 payload absent; R4.1 evidence is not promoted to V7"], "pmi": "NOT_PRESENT_OR_NOT_RECOVERED"})
write("brep_validation.json", {"status": "VERIFIED_FOR_R4_1_ONLY", "source_sha256": R41_SHA, "solid_count": 62, "available_evidence": "artifacts/engineering-evidence/brep_validation.json", "v7_status": "BLOCKED"})
write("feature_extraction_report.json", {"status": "BLOCKED", "source_sha256": R41_SHA, "features": [], "reason": "No released V7 payload and no deterministic V7 feature extraction result", "ai_authority": False})
write("pmi_report.json", {"status": "NOT_PRESENT_OR_NOT_RECOVERED", "source_sha256": R41_SHA, "dimensions": [], "tolerances": [], "gd_t": [], "datums": [], "material_callouts": [], "reason": "No released V7 PMI/drawing/native PMI evidence"})
write("chrono_runtime_status.json", {"status": "VERIFIED_CORE_PARTIAL", "source": "artifacts/pychrono_smoke.json", "core": "VERIFIED", "cascade": "BLOCKED", "v7_analysis": "BLOCKED", "reason": "No V7 physics inputs or Cascade STEP path"})

write("cae_evidence_ledger.json", {"status": "PARTIAL", "analyses": [
    {"analysis": "Gmsh synthetic mesh smoke", "input": "synthetic_engine_smoke/cube.geo", "solver": "Gmsh 4.15.2", "command": "gmsh -3 cube.geo -format msh2", "output": "synthetic_engine_smoke/run/cube.msh", "result": "VERIFIED", "evidence_class": "CALCULATED", "project_status": "NOT_V7_EVIDENCE"},
    {"analysis": "CalculiX synthetic static smoke", "input": "synthetic_engine_smoke/cube.inp", "solver": "CalculiX source build", "command": "CalculiX cube", "output": "synthetic_engine_smoke/run/cube.frd", "result": "VERIFIED", "evidence_class": "CALCULATED", "project_status": "NOT_V7_EVIDENCE"},
    {"analysis": "R4.1 OCP B-Rep validation", "input": "R4.1/R4.1.step", "solver": "OCP/OpenCascade", "command": "validate_r41_step.py", "output": "r41_step_validation.json", "result": "VERIFIED_FOR_R4_1", "evidence_class": "DERIVED_FROM_CAD", "project_status": "NOT_V7_EVIDENCE"},
    {"analysis": "V7 CAE", "input": None, "solver": None, "command": None, "output": None, "result": "BLOCKED", "evidence_class": "UNVERIFIED", "project_status": "BLOCKED"},
]})

write("engineering_evidence_graph.json", {"status": "PARTIAL", "nodes": [{"id": f"PART-{ref}", "type": "PART", "label": part, "evidence_class": "DERIVED_FROM_DOCUMENT", "verification_state": "UNVERIFIED"} for ref, part in parts] + [{"id": "R4.1", "type": "CAD_ARTIFACT", "label": "R4.1 STEP", "evidence_class": "DERIVED_FROM_CAD", "source_sha256": R41_SHA, "verification_state": "VERIFIED"}], "edges": [], "orphan_audit": "BLOCKED_PENDING_V7_PAYLOAD"})
write("traceability_audit.json", {"status": "PARTIAL", "source_sha256": R41_SHA, "orphan_nodes": [f"PART-{ref}" for ref, _ in parts], "reason": "V7 CAD-to-feature-to-CAE links cannot be established without V7 payload and released engineering records"})

physical = []
for i in range(1, 15):
    physical.append({"gap_id": f"O{i:02d}", "parameter": "physical input not released", "current_value": None, "unit": None, "evidence_class": "BLOCKED", "source": None, "source_hash": None, "acquisition_method": "authoritative document or calibrated physical test", "required_test": "define and execute after release", "uncertainty": None, "status": "BLOCKED", "owner": "Jon / engineering authority", "next_action": "provide authoritative input and provenance"})
write("physical_input_matrix.json", {"status": "BLOCKED", "physical_measurements_count": 0, "records": physical})

categories = [("A", 12, "Native CAD/Payload"), ("B", 20, "Deterministic Extraction"), ("C", 17, "Geometry/Interfaces"), ("D", 10, "PMI/Drawing/GD&T"), ("E", 11, "Joints/Hinges"), ("F", 9, "Materials"), ("G", 11, "Fasteners/Preload"), ("H", 12, "Lock-170"), ("I", 11, "Rebound-180"), ("J", 11, "Absorber"), ("K", 12, "Instrumentation"), ("L", 10, "Safety/Containment"), ("M", 11, "Configuration Management"), ("N", 10, "Traceability"), ("O", 14, "Physical Inputs")]
gaps = []
for letter, count, category in categories:
    for n in range(1, count + 1):
        gid = f"{letter}{n:02d}"
        gaps.append({"gap_id": gid, "description": f"{category} evidence for frozen V7 release item {gid}", "current_status": "BLOCKED", "blocker_type": "MISSING_AUTHORITATIVE_EVIDENCE", "automatable": letter in "ABCDMN", "required_evidence": ["released V7 CAD/document/test evidence"], "required_engine": ["OCP/STEPcode/Chrono/Gmsh/CalculiX as applicable"], "verification_method": "deterministic audit plus source-hash-bound evidence", "current_evidence": [], "next_action": "obtain authoritative V7 evidence and rerun audit", "owner_class": "engineering authority", "dependencies": ["V7 native CAD payload"], "closure_criteria": "source artifact, SHA, deterministic result, and review record", "last_verified": NOW})
write("comprehensive_gap_register.json", {"status": "BLOCKED_RELEASE", "design_candidate": "V7", "baseline": "V5", "gap_count": len(gaps), "gaps": gaps})
write("gap_automation_matrix.json", {"status": "PARTIAL", "records": [{"gap_id": g["gap_id"], "closure_class": "SOFTWARE-CLOSABLE" if g["automatable"] else "EVIDENCE-CLOSABLE", "status": "BLOCKED", "reason": "V7 source evidence absent"} for g in gaps]})
write("reproducibility_report.json", {"status": "PARTIAL", "verified_replays": ["R4.1 STEP validator", "OCP B-Rep inventory", "repository regression suite", "Gmsh synthetic smoke", "CalculiX synthetic smoke"], "blocked_replays": ["V7 CAD audit", "V7 CAE", "step-p21 project parse", "Chrono Cascade SR11 adapter"], "source_sha256": R41_SHA})
write("remaining_blockers.json", {"status": "BLOCKED_RELEASE", "critical": ["V7 native CAD payload absent", "released V7 BOM/materials/PMI/tolerances absent", "V7 interfaces and critical load paths not deterministically auditable", "V7 CAE inputs and acceptance thresholds absent"], "noncritical_or_external": ["OCCT 7.9.2 installation stopped", "step-p21 build stopped", "feature recognition not implemented", "physical measurements count is zero"], "release_decision": "STOP_RELEASE"})

write("master_validation_report.json", {"status": "STOP_RELEASE", "gates": {"A_repository": "VERIFIED_PARTIAL", "B_dependencies": "PARTIAL", "C_evidence": "VERIFIED_FOR_R4_1", "D_CAD": "BLOCKED_FOR_V7", "E_features": "BLOCKED", "F_PMI": "NOT_PRESENT_OR_NOT_RECOVERED", "G_traceability": "PARTIAL", "H_physics": "PARTIAL", "I_physical_evidence": "BLOCKED", "J_gap_closure": "BLOCKED", "K_reproducibility": "PARTIAL", "L_zero_bypass": "VERIFIED_FOR_PIPELINE_NOT_V7_RELEASE"}, "critical_gap_count": len(gaps), "manufacturability_critical_gap_count": len(parts), "prototype_release": "NOT_AUTHORIZED", "next_action": "Provide released V7 source package and engineering records; then rerun audit"})

# A candidate directory is created only as a blocked handoff container, never as a release claim.
candidate = OUT / "V7-CAD-RELEASE-CANDIDATE"
candidate.mkdir(parents=True, exist_ok=True)
(candidate / "STATUS.json").write_text(json.dumps({"status": "BLOCKED_RELEASE_CANDIDATE", "authorization": "NOT_AUTHORIZED", "reason": "V7 native CAD and critical evidence absent", "revision": "V7-CANDIDATE-BLOCKED-001", "source_sha256": R41_SHA}, indent=2) + "\n")
(candidate / "README.md").write_text("# V7-CAD-RELEASE-CANDIDATE\n\nThis is a blocked evidence container, not a release package. V7 native CAD, BOM, materials, PMI/GD&T, tolerances, interfaces, and V7 CAE evidence are absent. No manufacturer handoff is authorized.\n", encoding="utf-8")

print(json.dumps({"status": "STOP_RELEASE", "out": str(OUT), "gap_count": len(gaps), "critical_components": len(parts), "r41_sha256": R41_SHA}, indent=2))
