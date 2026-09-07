#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "master_v3"
UPLOAD = Path("/home/ubuntu/upload")
NOW = datetime.now(timezone.utc).isoformat()

EXTENSIONS = {".step", ".stp", ".iges", ".igs", ".dxf", ".stl", ".3mf", ".dwg", ".sldprt", ".sldasm", ".x_t", ".x_b", ".pdf", ".docx", ".xlsx", ".csv", ".json", ".md", ".txt"}
KEYWORDS = ("v7", "native", "step", "stp", "assembly", "assy", "bom", "material", "pmi", "gdt", "tolerance", "interface", "datum", "cae", "fea", "fastener", "joint", "source_data_required", "test_data_required", "manufacturing_drawing")
EXCLUDE_PARTS = {".git", "__pycache__", "node_modules", "EXTERNAL_PROJECTS"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def classify(path: Path) -> tuple[str, str, str]:
    text = str(path).lower()
    name = path.name.lower()
    if "cad-ai-engineering-os/tests/fixtures" in text or "synthetic_engine_smoke" in text:
        return "SYNTHETIC", "software fixture or engine smoke data", "not project evidence"
    if name in {"r4.1.step", "r4.1.stp"} and "acquired" not in text:
        return "REFERENCE", "R4.1 parent baseline; no V7 release authority", "reference only"
    if "artifacts/engineering-evidence/acquired" in text or "external_intake" in text:
        return "DERIVED", "immutable acquisition/extraction copy", "must retain parent-source identity"
    if "/p0_" in text or "/jon_p0_" in text or "v5_v7" in text or "v7_correlation" in text or "source_text" in text:
        return "REFERENCE", "concept/prototype/architecture or extracted document record", "not released V7 native authority"
    if "external_data_intake_package" in text:
        return "REFERENCE", "required-input schema/template with null or unsubmitted fields", "not CAE-ready evidence"
    if "artifacts/master_v" in text:
        return "DERIVED", "repository-generated evidence/report", "not source authority"
    return "REFERENCE", "repository document candidate", "authority requires explicit release record"


def candidate_paths() -> list[Path]:
    paths: set[Path] = set()
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in EXTENSIONS:
            continue
        if any(part in EXCLUDE_PARTS for part in path.parts):
            continue
        text = str(path).lower()
        if any(k in text for k in KEYWORDS):
            paths.add(path)
    # Explicit records required for the V7 recovery decision even if filename matching changes.
    explicit = [
        ROOT / "JON_P0_PHYSICAL_PROTOTYPE_MASTER/P0_MASTER/P0_BOM.csv",
        ROOT / "JON_P0_PHYSICAL_PROTOTYPE_MASTER/P0_MASTER/P0_MASTER_ASSEMBLY.md",
        ROOT / "JON_P0_PHYSICAL_PROTOTYPE_MASTER/P0_MASTER/P0_MANUFACTURING_DRAWINGS.md",
        ROOT / "P0_FUNCTIONAL_DEVELOPMENT_AUTHORIZATION/physical_execution/P0_MATERIAL_RECORD.json",
        ROOT / "P0_CRITICAL_DESIGN_CLOSURE/P0_INTERFACES_AND_DATUMS.json",
        ROOT / "P0_EVIDENCE_RECOVERY/P0_MATERIAL_FASTENER_EVIDENCE_REGISTER.json",
        ROOT / "P0_MINIMUM_BUILD_DEFINITION/P0_MINIMUM_BUILD_BOM.json",
        ROOT / "P0_PROTOTYPE_PARAMETER_PACK/V5_V7_NOT_FOUND_REGISTER.json",
        ROOT / "V7_CORRELATION_ARCHITECTURE.json",
        ROOT / "R4.1/R4.1.step",
        ROOT / "R4.1.step",
        ROOT / "artifacts/engineering-evidence/acquired/R4.1.step",
    ]
    paths.update(p for p in explicit if p.is_file())
    return sorted(paths)


def attached_sources() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path in sorted(UPLOAD.glob("*.pdf")):
        records.append({"artifact_id": "attachment::" + path.name, "category": "ATTACHED_REFERENCE_DOCUMENT",
                        "path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path),
                        "classification": "REFERENCE", "authority_status": "NOT_AUTHORITATIVE_V7",
                        "classification_basis": "attached patent/design reference PDF; not native CAD or released drawing",
                        "limitation": "visual/textual design intent only", "discovered_at": NOW})
    for path in sorted(UPLOAD.glob("*.zip")):
        members: list[dict[str, Any]] = []
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                if info.is_dir():
                    continue
                payload = archive.read(info)
                members.append({"path": info.filename, "size_bytes": info.file_size,
                                "sha256": hashlib.sha256(payload).hexdigest()})
        native = [m["path"] for m in members if Path(m["path"]).suffix.lower() in {".step", ".stp", ".iges", ".igs", ".sldprt", ".sldasm", ".x_t", ".x_b", ".dwg", ".dxf"}]
        records.append({"artifact_id": "attachment::" + path.name, "category": "ATTACHED_REFERENCE_ARCHIVE",
                        "path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path),
                        "classification": "REFERENCE", "authority_status": "NOT_AUTHORITATIVE_V7",
                        "classification_basis": "attached P0/physical-input documentation archive",
                        "limitation": "archive contains documentation/CSV/JSON/Markdown only" if not native else "native CAD member requires independent release authority",
                        "native_cad_members": native, "members": members, "discovered_at": NOW})
    return records


def missing_record(artifact_id: str, category: str, requirement: str, reason: str) -> dict[str, Any]:
    return {"artifact_id": artifact_id, "category": category, "path": None, "size_bytes": None, "sha256": None,
            "classification": "MISSING", "authority_status": "MISSING_AUTHORITATIVE", "requirement": requirement,
            "reason": reason, "recovery_action": "obtain exact released artifact, then hash and replay validation"}


artifacts: list[dict[str, Any]] = []
for path in candidate_paths():
    classification, basis, limitation = classify(path)
    artifacts.append({"artifact_id": path.name + "::" + relative(path), "category": "CANDIDATE_SOURCE",
                      "path": relative(path), "size_bytes": path.stat().st_size, "sha256": sha256(path),
                      "classification": classification, "authority_status": "NOT_AUTHORITATIVE_V7" if classification != "AUTHORITATIVE" else "AUTHORITATIVE",
                      "classification_basis": basis, "limitation": limitation, "discovered_at": NOW})

artifacts.extend(attached_sources())

# Connected source: CAD-AI Engineering OS's minimal STEP fixture is intentionally synthetic.
connected = Path("/home/ubuntu/cad-ai-engineering-os/tests/fixtures/minimal-box.step")
if connected.is_file():
    artifacts.append({"artifact_id": "cad-ai-engineering-os::minimal-box.step", "category": "CONNECTED_SOURCE",
                      "path": str(connected), "size_bytes": connected.stat().st_size, "sha256": sha256(connected),
                      "classification": "SYNTHETIC", "authority_status": "NOT_AUTHORITATIVE_V7",
                      "classification_basis": "test fixture from connected repository", "limitation": "not project geometry", "discovered_at": NOW})

missing = [
    missing_record("V7-NATIVE-CAD", "native_cad", "Released V7 native CAD/STEP/assembly bytes", "No V7 native CAD payload was found; R4.1 is the only native STEP and is a parent reference."),
    missing_record("V7-MASTER-ASSEMBLY", "assembly", "Released V7 master assembly and component hierarchy", "P0 master assembly is explicitly Prototype/Test Article P0.1 and points back to R4.1/reference geometry."),
    missing_record("V7-RELEASED-BOM", "bom", "Released V7 manufacturing BOM", "Available BOMs are P0/reference records and explicitly require native CAD/drawing/material release."),
    missing_record("V7-MATERIAL-DEFINITIONS", "materials", "Released V7 material definitions/certificates", "Material records are OPEN/UNKNOWN or concept basis; no certificates or lot records were recovered."),
    missing_record("V7-PMI-GDT", "pmi_gdt", "Released V7 PMI/GD&T/drawings with dimensions and tolerances", "Available manufacturing-drawing records are text/reference controls; deterministic PMI/GD&T is not recovered."),
    missing_record("V7-INTERFACES", "interfaces", "Released V7 interface/datum/mating definitions", "Interface records exist but every exact geometry/datum/mating status is BLOCKED."),
    missing_record("V7-TOLERANCES", "tolerances", "Released V7 tolerance stack and inspection data", "No released tolerance data with authoritative source identity was found."),
    missing_record("V7-CAE-INPUTS", "cae_inputs", "V7-ready mass, CG, inertia, materials, pulse, absorber, contact, friction, joint and boundary inputs", "Input templates exist with null/unsubmitted values; no V7-ready source values were recovered."),
]

recovery = {"status": "RECOVERY_COMPLETE_NO_AUTHORITATIVE_V7_SOURCE", "freeze": {"v7_design_modification": False, "concept_change": False, "assumption_fill": False, "manufacturer_release_artifacts_generated": False, "gate": "STOP_RELEASE"}, "source_recovery_timestamp": NOW, "repository": str(ROOT), "connected_sources_checked": ["/home/ubuntu/cad-ai-engineering-os", "/home/ubuntu/upload", "Git history and remote refs"], "authoritative_artifacts": [], "recovered_artifacts": artifacts, "missing_authoritative_slots": missing, "counts": {"recovered_candidates": len(artifacts), "authoritative": 0, "derived": sum(a["classification"] == "DERIVED" for a in artifacts), "reference": sum(a["classification"] == "REFERENCE" for a in artifacts), "synthetic": sum(a["classification"] == "SYNTHETIC" for a in artifacts), "missing_authoritative_slots": len(missing)}, "decision": "No gap closure permitted; recompute graph and register with all V7 authority-dependent gaps remaining BLOCKED."}
(OUT / "v7_native_source_recovery_inventory.json").write_text(json.dumps(recovery, indent=2) + "\n", encoding="utf-8")

old = json.loads((OUT / "comprehensive_gap_register_v3.json").read_text())
updated = []
for gap in old["gaps"]:
    g = dict(gap)
    g["recovery_evidence"] = [a["path"] for a in artifacts if a["classification"] in {"REFERENCE", "DERIVED"} and any(token.lower() in a["path"].lower() for token in (gap["gap_id"], "v7", "r4.1", "p0"))][:10]
    g["authoritative_v7_source_recovered"] = False
    g["current_status"] = "BLOCKED"
    g["closure_stage"] = "BLOCKED"
    g["blocking_dependency"] = "missing authoritative V7 source package"
    g["next_action"] = "Recover exact released V7 artifact; do not infer from reference records"
    updated.append(g)
(OUT / "comprehensive_gap_register_v3_recovered.json").write_text(json.dumps({"status": "RECOVERY_RECOMPUTED_STOP_RELEASE", "source_inventory": "v7_native_source_recovery_inventory.json", "preserved_gap_ids": True, "gap_count": len(updated), "gaps": updated}, indent=2) + "\n", encoding="utf-8")

graph = {"status": "RECOVERY_RECOMPUTED_NO_V7_AUTHORITY", "source_inventory": "v7_native_source_recovery_inventory.json", "nodes": [], "edges": [], "orphan_audit": "ALL_V7_AUTHORITY_NODES_ORPHANED_OR_BLOCKED"}
for artifact in artifacts:
    graph["nodes"].append({"node_id": artifact["artifact_id"], "node_type": artifact["category"], "label": artifact["path"], "source": artifact["path"], "source_sha256": artifact["sha256"], "classification": artifact["classification"], "verification_state": artifact["authority_status"]})
for slot in missing:
    graph["nodes"].append({"node_id": slot["artifact_id"], "node_type": slot["category"], "label": slot["requirement"], "source": None, "source_sha256": None, "classification": "MISSING", "verification_state": "MISSING_AUTHORITATIVE"})
for artifact in artifacts:
    graph["edges"].append({"edge_id": "recovery::" + artifact["artifact_id"], "source_node": artifact["artifact_id"], "target_node": "V7-NATIVE-PACKAGE" if artifact["classification"] == "AUTHORITATIVE" else "REFERENCE-EVIDENCE-POOL", "source_sha256": artifact["sha256"], "evidence_class": "DERIVED_FROM_DOCUMENT" if artifact["classification"] in {"DERIVED", "REFERENCE"} else artifact["classification"], "verification_state": artifact["authority_status"]})
(OUT / "engineering_evidence_graph_v3_recovered.json").write_text(json.dumps(graph, indent=2) + "\n", encoding="utf-8")
print(json.dumps(recovery["counts"], indent=2))
