#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "master_v4"
NOW = datetime.now(timezone.utc).isoformat()
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(name: str, data: object) -> None:
    (OUT / name).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

bom_path = ROOT / "JON_P0_PHYSICAL_PROTOTYPE_MASTER/P0_MASTER/P0_BOM.csv"
bom_sha = sha256(bom_path)
components = []
with bom_path.open(encoding="utf-8-sig", newline="") as f:
    for row in csv.DictReader(f):
        components.append({
            "part_id": row["P0_ID"],
            "reference": row["REF"],
            "component": row["NAME"],
            "part_function": row["FUNCTION"],
            "material": row["MATERIAL_BASIS"],
            "manufacturing_process": row["MANUFACTURING"],
            "joining_method": "NOT_RELEASED; derive from authoritative native CAD/drawing",
            "critical_dimensions": "NOT_RELEASED; native R4.1/V7 CAD required",
            "tolerances": "NOT_RELEASED",
            "assembly_method": "P0 prototype/test article only; production assembly not authorized",
            "manufacturing_risk": "HIGH — source row requires native CAD, material certificate, drawing release, or characterization",
            "evidence": [{"path": str(bom_path.relative_to(ROOT)), "sha256": bom_sha, "classification": "REFERENCE"}],
            "status": "BLOCKED",
            "release_authority": "NOT_PRESENT",
        })

inventory = json.loads((ROOT / "artifacts/master_v3/v7_native_source_recovery_inventory.json").read_text())
attachment_refs = [{"path": a["path"], "sha256": a["sha256"], "classification": a["classification"]} for a in inventory["recovered_artifacts"] if a["category"].startswith("ATTACHED")]

forensic = {
    "status": "BLOCKED",
    "scope": "FROZEN_V7_FORENSIC_AUDIT",
    "design_modification": False,
    "concept_change": False,
    "assumption_fill": False,
    "source_inventory": "../master_v3/v7_native_source_recovery_inventory.json",
    "native_v7_authoritative_sources": 0,
    "r4_1_reference_sha256": "fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68",
    "audited_domains": ["geometry", "interfaces", "dimensions", "thicknesses", "materials", "fasteners", "joints_welds", "mounting_points", "clearances", "mechanisms", "assembly_feasibility", "service_accessibility", "critical_load_paths"],
    "findings": [
        {"domain": "geometry", "status": "BLOCKED", "finding": "V7 native CAD not recovered; R4.1 is reference only"},
        {"domain": "interfaces", "status": "BLOCKED", "finding": "exact datums/mating/clearances not released"},
        {"domain": "materials", "status": "UNVERIFIED", "finding": "concept material basis exists; certificates and lot mapping absent"},
        {"domain": "fasteners", "status": "UNVERIFIED", "finding": "concept grades exist; released identity/preload evidence absent"},
        {"domain": "mechanisms", "status": "REFERENCE", "finding": "stable numerals and functional architecture documented, exact implementation not released"},
        {"domain": "assembly_feasibility", "status": "BLOCKED", "finding": "P0 assembly is explicitly a prototype/test article and not production hardware"},
        {"domain": "critical_load_paths", "status": "REFERENCE", "finding": "architecture and intended paths documented; no V7 structural validation"},
    ],
    "component_count": len(components),
    "timestamp": NOW,
}
dump("v4_forensic_audit.json", forensic)

dump("manufacturability_proof.json", {"status": "BLOCKED", "scope": "P0 reference records applied to frozen V7 audit; not a manufacturer release", "source": {"path": str(bom_path.relative_to(ROOT)), "sha256": bom_sha, "classification": "REFERENCE"}, "components": components, "critical_manufacturability_gaps": ["released native CAD/drawings", "material certificates and lot mapping", "released BOM", "PMI/GD&T and tolerance stack", "joining/weld specifications", "inspection and acceptance criteria"]})

dump("cae_evidence_ledger.json", {"status": "BLOCKED", "entries": [
    {"analysis": "R4.1 STEP/B-Rep forensic reference audit", "input": "R4.1/R4.1.step", "input_sha256": "fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68", "solver": "OpenCascade/OCP", "version": "repository-verified", "command": "python3 validate_r41_step.py R4.1/R4.1.step", "output_artifact": "artifacts/r41_step_validation.json", "result": "reference geometry parsed; not V7", "hash": None, "status": "PARTIAL"},
    {"analysis": "V7 structural FEA", "input": None, "input_sha256": None, "solver": "CalculiX", "version": None, "command": None, "output_artifact": None, "result": "V7 native geometry/material/loads unavailable", "hash": None, "status": "BLOCKED"},
    {"analysis": "V7 multibody dynamics", "input": None, "input_sha256": None, "solver": "Project Chrono", "version": "10.0.0 core verified", "command": None, "output_artifact": None, "result": "V7 configuration and calibrated physical inputs unavailable", "hash": None, "status": "BLOCKED"},
    {"analysis": "Synthetic engine smokes", "input": "artifacts/master_v2/synthetic_engine_smoke", "input_sha256": None, "solver": "Gmsh/CalculiX", "version": "repository-recorded", "command": "synthetic fixture smoke commands", "output_artifact": "synthetic outputs", "result": "software availability only", "hash": None, "status": "VERIFIED_NOT_PROJECT_EVIDENCE"}
]})

dump("v4_release_gate.json", {"status": "STOP_RELEASE", "v7_status": "BLOCKED", "critical_gap_count": "NONZERO", "manufacturability_critical_gap_count": "NONZERO", "available_cae_evidence_collected": False, "all_blocked_items_explicitly_identified": True, "prototype_release_package_generated": False, "manufacturer_handoff_authorized": False, "reasons": ["no authoritative V7 native CAD/assembly", "no released V7 BOM/material/PMI/tolerance package", "no V7-ready CAE input set", "physical validation and calibration absent", "Peter reauthorization not present"]})
print(json.dumps({"status": "STOP_RELEASE", "components": len(components), "attachment_refs": len(attachment_refs)}, indent=2))
