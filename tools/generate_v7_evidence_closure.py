#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "master_v4"
OUT.mkdir(parents=True, exist_ok=True)
MATRIX = json.loads((OUT / "pdf_functional_standards_matrix.json").read_text(encoding="utf-8"))
INVENTORY = json.loads((ROOT / "artifacts" / "master_v3" / "v7_native_source_recovery_inventory.json").read_text(encoding="utf-8"))
REVISION = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
NOW = datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def evidence_for(requirement: dict) -> dict[str, object]:
    rid = requirement["id"]
    ref_map = {
        "FR-01": {"cad_entity": "110L/110R load-path entities", "interface": "left/right structural interfaces", "load_path": "bilateral longitudinal reaction path", "test_method": "left/right seat-strength and reaction-channel test", "acceptance": "released criterion for bilateral load sharing"},
        "FR-02": {"cad_entity": "130 energy-dissipation cartridge entities", "interface": "cartridge-to-carriage mount", "load_path": "longitudinal force-stroke path", "test_method": "force-stroke, hysteresis, temperature, cycling, sled test", "acceptance": "released force/stroke and repeatability limits"},
        "FR-03": {"cad_entity": "150L/150R rotation-control entities", "interface": "seatback-to-carriage rotational interfaces", "load_path": "distinct rotational reaction path", "test_method": "isolated rotation and link-load test", "acceptance": "released rotation/load and interaction limits"},
        "FR-04": {"cad_entity": "170 state-control entities", "interface": "state trigger and restraint interfaces", "load_path": "armed/ride-down/rebound/secure transition path", "test_method": "state-transition, latency, unintended-release and fail-safe tests", "acceptance": "released transition and fault criteria"},
        "FR-05": {"cad_entity": "180 rebound-control entities", "interface": "rebound stop/capture interface", "load_path": "reverse-travel reaction path", "test_method": "reverse-travel and rebound engagement test", "acceptance": "released rebound travel/load limits"},
        "FR-06": {"cad_entity": "190 restraint-path entities", "interface": "released anchorage and belt-routing interfaces", "load_path": "restraint load path during ride-down", "test_method": "anchorage strength and dynamic restraint test", "acceptance": "released anchorage geometry and strength criteria"},
        "FR-07": {"cad_entity": "140 pelvic-control entities", "interface": "occupant/ATD restraint interfaces", "load_path": "pelvic restraint and forward-migration path", "test_method": "ATD/sled kinematics and load-distribution test", "acceptance": "released migration/load criteria"},
        "FR-08": {"cad_entity": "replaceable cartridge entities", "interface": "cartridge replacement and mounting interfaces", "load_path": "cartridge force-stroke path", "test_method": "cartridge characterization and serial traceability test", "acceptance": "released region, interchangeability and traceability criteria"},
        "FR-09": {"cad_entity": "left/right reaction measurement features", "interface": "sensor mounting and channel interfaces", "load_path": "separate left/right reaction channels", "test_method": "calibrated ISO 6487/SAE J211 measurement test", "acceptance": "released calibration, CFC, timing and uncertainty limits"},
        "FR-10": {"cad_entity": "all real parts, joints and travel-limit entities", "interface": "assembly interfaces and joint definitions", "load_path": "all declared structural and dynamic paths", "test_method": "assembly inspection, joint/travel measurement and material verification", "acceptance": "released CAD/BOM/material/measurement criteria"},
        "FR-11": {"cad_entity": "stable reference numerals across assembly", "interface": "CAD-test-CAE identifier mapping", "load_path": "traceability chain", "test_method": "hashed mapping-graph audit", "acceptance": "every released reference numeral resolves to source and evidence"},
        "FR-12": {"cad_entity": "V0-V6 gate-specific entities", "interface": "gate handoff interfaces", "load_path": "gate-specific load and evidence paths", "test_method": "gate-specific physical, CAE and inspection methods", "acceptance": "released gate criteria and signed evidence package"},
    }
    extra = ref_map[rid]
    return {
        "requirement_id": rid,
        "requirement": requirement["requirement"],
        "source": requirement["source"],
        "standards": requirement["standards"],
        "mapping_baseline": requirement["mapping"],
        "cad_entity": extra["cad_entity"],
        "interface": extra["interface"],
        "load_path": extra["load_path"],
        "material": "MISSING_AUTHORITATIVE_V7_MATERIAL_DEFINITION",
        "tolerance": "MISSING_AUTHORITATIVE_V7_PMI_GDT_TOLERANCE",
        "cae_input": "MISSING_AUTHORITATIVE_V7_CAE_INPUT",
        "test_method": extra["test_method"],
        "acceptance_criterion": extra["acceptance"],
        "available_evidence": [],
        "evidence_class": "MISSING",
        "verification_state": "BLOCKED",
        "closure_status": "NOT_CLOSED",
        "reason": "No released V7 authoritative native source, material, PMI/GD&T, interface, tolerance, CAE input, or physical test evidence was recovered.",
    }


requirements = [evidence_for(req) for req in MATRIX["requirements"]]
missing = [
    {"slot_id": "V7-NATIVE-CAD", "required_for": [r["requirement_id"] for r in requirements], "description": "Released V7 native CAD/STEP/assembly bytes", "classification": "MISSING"},
    {"slot_id": "V7-ASSEMBLY", "required_for": ["FR-01", "FR-03", "FR-04", "FR-05", "FR-06", "FR-10", "FR-11", "FR-12"], "description": "Released V7 assembly hierarchy, joints, travel limits, clearances and interfaces", "classification": "MISSING"},
    {"slot_id": "V7-BOM", "required_for": ["FR-10", "FR-11", "FR-12"], "description": "Released V7 manufacturing BOM and stable item identifiers", "classification": "MISSING"},
    {"slot_id": "V7-MATERIALS", "required_for": ["FR-01", "FR-02", "FR-03", "FR-05", "FR-08", "FR-10", "FR-12"], "description": "Released material definitions, certificates, lot mapping and CAE material cards", "classification": "MISSING"},
    {"slot_id": "V7-PMI-GDT", "required_for": ["FR-01", "FR-03", "FR-05", "FR-06", "FR-08", "FR-10", "FR-11", "FR-12"], "description": "Released PMI/GD&T, drawings, dimensions, datums and tolerances", "classification": "MISSING"},
    {"slot_id": "V7-CAE-INPUTS", "required_for": ["FR-01", "FR-02", "FR-03", "FR-04", "FR-05", "FR-06", "FR-07", "FR-08", "FR-09", "FR-12"], "description": "Released mass, CG, inertia, pulse, absorber, contact, friction, joint and boundary inputs", "classification": "MISSING"},
    {"slot_id": "V7-PHYSICAL-TESTS", "required_for": [r["requirement_id"] for r in requirements], "description": "V7 physical prototype test records, calibrated channels, raw data and acceptance results", "classification": "MISSING"},
]

inventory_summary = {
    "source_inventory": "artifacts/master_v3/v7_native_source_recovery_inventory.json",
    "inventory_sha256": sha256(ROOT / "artifacts/master_v3/v7_native_source_recovery_inventory.json"),
    "authoritative_artifacts": INVENTORY.get("authoritative_artifacts", []),
    "counts": INVENTORY.get("counts", {}),
    "classification_policy": ["AUTHORITATIVE", "DERIVED", "REFERENCE", "SYNTHETIC", "MISSING"],
    "recovery_result": INVENTORY.get("status"),
    "sha_binding": "Every recovered file record in the source inventory carries a SHA-256; no authoritative V7 artifact was found.",
}

trace = {
    "status": "TRACEABILITY_COMPLETE_WITH_ALL_V7_CLOSURES_BLOCKED",
    "generated_at": NOW,
    "repository_revision": REVISION,
    "baseline_matrix": "artifacts/master_v4/pdf_functional_standards_matrix.json",
    "baseline_matrix_sha256": sha256(OUT / "pdf_functional_standards_matrix.json"),
    "freeze": {"v7_design_modification": False, "concept_change": False, "assumption_fill": False, "manufacturer_release_artifacts_generated": False},
    "requirements": requirements,
    "summary": {"requirement_count": len(requirements), "closed": 0, "blocked": len(requirements), "missing_authoritative_evidence_slots": len(missing)},
    "release_gate": "STOP_RELEASE",
}
(OUT / "v7_12_requirement_traceability_matrix.json").write_text(json.dumps(trace, indent=2) + "\n", encoding="utf-8")
(OUT / "v7_missing_authoritative_evidence.json").write_text(json.dumps({"status": "MISSING_AUTHORITATIVE_EVIDENCE", "requirements": [r["requirement_id"] for r in requirements], "missing": missing, "release_gate": "STOP_RELEASE"}, indent=2) + "\n", encoding="utf-8")
(OUT / "v7_package_inventory_closure.json").write_text(json.dumps(inventory_summary, indent=2) + "\n", encoding="utf-8")

nodes = [{"node_id": r["requirement_id"], "node_type": "V7_REQUIREMENT", "label": r["requirement"], "classification": "MISSING", "verification_state": "BLOCKED"} for r in requirements]
nodes += [{"node_id": s["slot_id"], "node_type": "MISSING_EVIDENCE", "label": s["description"], "classification": "MISSING", "verification_state": "MISSING_AUTHORITATIVE"} for s in missing]
for artifact in INVENTORY.get("recovered_artifacts", []):
    nodes.append({"node_id": artifact["artifact_id"], "node_type": artifact["category"], "label": artifact["path"], "source": artifact["path"], "source_sha256": artifact["sha256"], "classification": artifact["classification"], "verification_state": artifact.get("authority_status")})
edges = [{"edge_id": f"{r['requirement_id']}::missing", "source_node": r["requirement_id"], "target_node": "V7-NATIVE-CAD", "evidence_class": "MISSING", "verification_state": "BLOCKED"} for r in requirements]
(OUT / "v7_evidence_graph_closure.json").write_text(json.dumps({"status": "RECOMPUTED_STOP_RELEASE", "generated_at": NOW, "source_inventory": inventory_summary["inventory_sha256"], "nodes": nodes, "edges": edges, "orphan_audit": "ALL_12_REQUIREMENTS_REMAIN_BLOCKED_BY_MISSING_AUTHORITATIVE_V7_EVIDENCE"}, indent=2) + "\n", encoding="utf-8")

old = json.loads((ROOT / "artifacts/master_v3/comprehensive_gap_register_v3_recovered.json").read_text(encoding="utf-8"))
for gap in old.get("gaps", []):
    gap["authoritative_v7_source_recovered"] = False
    gap["current_status"] = "BLOCKED"
    gap["closure_stage"] = "BLOCKED"
    gap["blocking_dependency"] = "missing authoritative V7 source package"
    gap["next_action"] = "Recover and SHA-bind exact released V7 evidence; do not infer from reference or synthetic records"
    gap["traceability_requirements"] = [r["requirement_id"] for r in requirements if gap.get("gap_id", "").startswith(r["requirement_id"])]
old.update({"status": "V7_EVIDENCE_CLOSURE_RECOMPUTED_STOP_RELEASE", "closure_matrix": "v7_12_requirement_traceability_matrix.json", "all_v7_gaps_closed": False})
(ROOT / "artifacts/master_v3/comprehensive_gap_register_v3_recovered.json").write_text(json.dumps(old, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"requirements": len(requirements), "closed": 0, "blocked": len(requirements), "missing_slots": len(missing), "gate": "STOP_RELEASE"}, indent=2))
