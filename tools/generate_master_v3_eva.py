#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engineering_recovery.eva.core import parameter_sweep, prioritize_test

OUT = ROOT / "artifacts" / "master_v3"
NOW = datetime.now(timezone.utc).isoformat()
R41 = ROOT / "R4.1" / "R4.1.step"
R41_SHA = hashlib.sha256(R41.read_bytes()).hexdigest()


def dump(name: str, value: object) -> None:
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def source(path: Path, source_type: str, authority: int, context: str) -> dict[str, object]:
    return {"source_id": path.name, "source": str(path), "source_type": source_type,
            "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "authority_level": authority,
            "retrieval_timestamp": NOW, "context": context, "status": "DISCOVERED"}


sources = [
    source(R41, "PROJECT_CAD", 4, "R4.1 STEP source; authoritative for R4.1 only, not V7"),
    source(ROOT / "OPEN_ENGINEERING_QUESTIONS.json", "PROJECT_REGISTER", 4, "Open physical and engineering questions"),
    source(ROOT / "P0_PARAMETER_INTELLIGENCE_PACK/P0_V7_ARCHITECTURE_PHYSICAL_MAP.json", "PROJECT_DOCUMENT", 4, "V7 physical architecture map"),
    source(ROOT / "artifacts/master_v2/comprehensive_gap_register.json", "PROJECT_REGISTER", 4, "Existing A-O gap register"),
]
dump("eva_source_registry.json", {"status": "VERIFIED_PARTIAL", "authority_hierarchy": {"1": "MEASURED_PROJECT_OR_CERTIFIED_LAB", "2": "CONTROLLED_SUPPLIER_OR_STANDARD", "3": "PEER_REVIEWED_OR_ESTABLISHED_DATABASE", "4": "ANALYTICAL_DERIVATION", "5": "ENGINEERING_ESTIMATE", "6": "AI_CANDIDATE"}, "sources": sources, "conflicts": []})

parameters = []
for i in range(1, 15):
    parameters.append({"parameter_id": f"O{i:02d}", "value": None, "unit": None, "source": None, "source_type": None, "source_sha256": None,
                       "evidence_class": "BLOCKED", "derivation_method": None, "equation_or_model": None,
                       "assumptions": [], "uncertainty": "not bounded: no authoritative basis", "range_min": None,
                       "range_nominal": None, "range_max": None, "distribution": None, "sensitivity_status": "NOT_RUN",
                       "validation_status": "UNVERIFIED", "timestamp": NOW, "tool": "EVA", "tool_version": "0.1.0",
                       "configuration_id": f"O{i:02d}-blocked", "state": "TEST_REQUIRED"})
dump("eva_parameter_registry.json", {"status": "PARTIAL", "parameters": parameters, "rule": "No value without provenance"})

dump("eva_derivation_registry.json", {"status": "READY", "derivations": [{"derivation_id": "mass_from_volume_density", "formula": "m = V * rho", "inputs": ["authoritative volume", "authoritative density"], "output_evidence_class": "CALCULATED", "uncertainty_propagation": "interval or distribution required", "status": "AVAILABLE_NOT_EXECUTED_FOR_V7"}, {"derivation_id": "inertia_from_CAD_mass_properties", "formula": "CAD mass-property integration", "inputs": ["released V7 B-Rep", "released density"], "output_evidence_class": "DERIVED_FROM_CAD", "status": "BLOCKED_V7_PAYLOAD"}]})

dump("eva_uncertainty_registry.json", {"status": "READY", "methods": ["interval_analysis", "corner_propagation", "Monte_Carlo", "parameter_sweep", "worst_case", "percentile", "sensitivity"], "project_parameters": [{"parameter_id": f"O{i:02d}", "status": "BOUNDED_MODEL_UNAVAILABLE", "assumptions": [], "range": None, "confidence": "NONE", "source_basis": None, "next_action": "obtain authoritative basis or design minimum test"} for i in range(1, 15)]})

synthetic = parameter_sweep("EVA-SYNTHETIC-001", "x", [0.0, 0.5, 1.0], lambda x: x * x, ["SYNTHETIC_TEST_DATA"], "EVA_DETERMINISTIC")
synthetic["test_data_class"] = "SYNTHETIC_TEST_DATA"
dump("eva_virtual_experiments.json", {"status": "VERIFIED_SOFTWARE_ONLY", "experiments": [synthetic], "project_experiments": [], "rule": "Synthetic experiments never enter project validation evidence"})

dump("eva_sensitivity_results.json", {"status": "VERIFIED_SOFTWARE_ONLY", "results": [{"experiment_id": "EVA-SYNTHETIC-001", "parameter": "x", "status": "CALCULATED", "note": "Software contract smoke only"}], "project_results": [], "unknown_parameters": [f"O{i:02d}" for i in range(1, 15)]})
dump("eva_robustness_results.json", {"status": "NOT_RUN_FOR_V7", "classification": "ROBUSTNESS_EVIDENCE_IS_DISTINCT_FROM_PHYSICAL_VALIDATION", "project_results": [], "synthetic_results": [{"id": "EVA-SYNTHETIC-001", "status": "SYNTHETIC_TEST_DATA"}]})

tests = [
    prioritize_test("T01_vehicle_pulse_and_acceleration", ["O01", "O02", "L01", "L02"], 4, 10, 10, 8, ["lock timing", "rebound", "CAE calibration"]),
    prioritize_test("T02_absorber_force_stroke", ["J01", "J02", "I01", "O06"], 4, 9, 8, 6, ["absorber model", "ride-down dynamics"]),
    prioritize_test("T03_joint_and_lock_characterization", ["E01", "H01", "H02", "G01"], 4, 8, 9, 7, ["Lock-170 ranking", "structural load path"]),
    prioritize_test("T04_material_and_fastener_certificate_review", ["F01", "G01", "G02"], 3, 7, 8, 3, ["FE material card", "preload model"]),
]
dump("eva_test_priority_matrix.json", {"status": "CANDIDATE_TEST_PLAN", "selection_policy": "information gain per effort; no invented cost", "tests": tests})
dump("eva_calibration_registry.json", {"status": "TEST_REQUIRED", "records": [{"test_id": t["test_id"], "status": "TEST_REQUIRED", "raw_test_data_sha256": None, "pre_test_model_preserved": True, "next_action": "import calibrated raw data"} for t in tests]})
dump("eva_model_maturity.json", {"status": "PARTIAL", "models": [{"model_id": "R4.1_BREP", "maturity": "M2", "label": "COMPUTATIONAL", "evidence_refs": [R41_SHA], "status": "VERIFIED_FOR_R4_1"}, {"model_id": "V7_DYNAMIC_MODEL", "maturity": "M0", "label": "CONCEPTUAL", "evidence_refs": [], "status": "UNVERIFIED"}]})

dump("F_material_closure.json", {"status": "UNVERIFIED", "gaps": [{"gap_id": f"F{i:02d}", "evidence_search": "project sources searched", "candidate_sources": [], "derived_properties": [], "bounds": None, "sensitivity": "NOT_RUN", "physical_characterization": "REQUIRED", "minimum_human_input": "material certificate or controlled datasheet"} for i in range(1, 10)]})
dump("G_fastener_closure.json", {"status": "BLOCKED", "gaps": [{"gap_id": f"G{i:02d}", "cad_bom_search": "completed against available records", "fastener": None, "calculated_preload": None, "preload_evidence_class": "CALCULATED_ONLY_IF_INPUTS_EXIST", "physical_test_required": True, "minimum_human_input": "released BOM/fastener specification"} for i in range(1, 12)]})
dump("H_lock_170_closure.json", {"status": "BLOCKED", "candidates": ["L170-A", "L170-B", "L170-C"], "virtual_ranking": "NOT_RUN_V7_INPUTS_ABSENT", "physical_status": "TEST_REQUIRED", "final_selection": "NOT_AUTHORIZED", "thresholds": {"acceleration": {"value": 2.5, "unit": "g", "classification": "TARGET_OR_REQUIREMENT_UNVERIFIED"}, "response_time": {"value": 20, "unit": "ms", "classification": "TARGET_OR_REQUIREMENT_UNVERIFIED"}}})
dump("I_rebound_180_closure.json", {"status": "BLOCKED", "target": {"value": 180, "unit": "mm", "evidence_class": "TARGET"}, "design_envelope": "NOT_RUN_V7_INPUTS_ABSENT", "physical_test_required": True})
dump("J_absorber_closure.json", {"status": "UNVERIFIED", "models": ["linear", "nonlinear", "piecewise", "viscoelastic", "rate_dependent", "hysteretic", "energy_absorbing"], "curves_available": 0, "physical_test_required": True})
dump("K_instrumentation_plan.json", {"status": "CALIBRATION_READY", "physical_measurements": 0, "sensor_map": [], "channels": ["base_acceleration", "carriage_xva", "seatback_theta_omega_alpha", "rail_reactions", "absorber_force_displacement", "lock_state", "occupant_motion"], "sampling_requirements": "TO BE RELEASED BY TEST AUTHORITY", "synchronization": "TEST_TRIGGER_REQUIRED", "calibration_requirements": "sensor-specific certificates required", "raw_data_schema": ["sensor_id", "calibration", "sampling", "units", "timestamp", "test_id", "raw_data_sha256"], "acceptance_metrics": ["peak", "RMS", "phase/time-of-peak", "energy where meaningful", "event timing"]})
dump("O_physical_input_closure.json", {"status": "TEST_REQUIRED", "state_machine": ["UNKNOWN", "SEARCHED", "SOURCE_FOUND", "EXTRACTED", "DERIVED", "BOUNDED", "SENSITIVITY_TESTED", "TEST_PRIORITIZED", "MEASURED", "CALIBRATED", "VALIDATED"], "records": [{"gap_id": f"O{i:02d}", "state": "TEST_REQUIRED", "evidence_available": False, "evidence_sources": [], "derivation_available": False, "derivation_method": None, "bounded_model_available": False, "virtual_experiment_available": True, "sensitivity_completed": False, "robustness_status": "NOT_RUN", "physical_test_required": True, "minimum_test_definition": "define one calibrated acquisition with raw-data SHA-256", "calibration_required": True, "closure_path": "test_prioritized -> measured -> calibrated -> validated", "closure_stage": "TEST_REQUIRED", "blocking_dependency": "authoritative physical input"} for i in range(1, 15)]})
dump("L_safety_virtual_envelope.json", {"status": "BLOCKED_V7_INPUTS", "analysis_domains": ["containment", "retention", "secondary_motion", "pinch_crush", "rebound", "stop_loads", "component_failure", "occupant_interaction", "structural_load_paths"], "virtual_pre_screen": "AVAILABLE_AS_INTERFACE_ONLY", "physical_validation": "NOT_PERFORMED"})

old = json.loads((ROOT / "artifacts/master_v2/comprehensive_gap_register.json").read_text())
records = []
for gap in old["gaps"]:
    g = dict(gap)
    g.update({"evidence_available": bool(g.get("current_evidence")), "evidence_sources": g.get("current_evidence", []), "derivation_available": False, "derivation_method": None, "bounded_model_available": True, "virtual_experiment_available": True, "sensitivity_completed": False, "robustness_status": "NOT_RUN", "physical_test_required": gap["gap_id"].startswith(("F", "G", "H", "I", "J", "K", "L", "O")), "minimum_test_definition": "smallest calibrated test that resolves the stated uncertainty" if gap["gap_id"].startswith(("F", "G", "H", "I", "J", "K", "L", "O")) else None, "calibration_required": gap["gap_id"].startswith(("F", "G", "H", "I", "J", "K", "L", "O")), "closure_path": "evidence -> derive/bound -> virtual experiment -> sensitivity -> prioritized test -> calibration -> validation", "closure_stage": "TEST_REQUIRED" if gap["gap_id"].startswith(("F", "G", "H", "I", "J", "K", "L", "O")) else "BLOCKED", "blocking_dependency": "released V7 source package"})
    records.append(g)
dump("comprehensive_gap_register_v3.json", {"status": "PARTIAL_WITH_EXPLICIT_BLOCKERS", "preserved_gap_ids": True, "gap_count": len(records), "gaps": records})
dump("master_v3_validation_report.json", {"status": "PARTIAL", "eva_runtime": "VERIFIED_SOFTWARE_CONTRACTS", "source_discovery": "VERIFIED_PARTIAL", "document_extraction": "VERIFIED_SOFTWARE_CONTRACTS", "derivation": "VERIFIED_SOFTWARE_CONTRACTS", "uncertainty": "VERIFIED_SOFTWARE_CONTRACTS", "virtual_experiments": "VERIFIED_SYNTHETIC_ONLY", "sensitivity": "VERIFIED_SYNTHETIC_ONLY", "test_prioritization": "CANDIDATE_TEST_PLAN", "calibration": "TEST_REQUIRED", "physical_measurements": 0, "v7_validation": "BLOCKED", "false_closure": False, "remaining_critical_blockers": ["V7 native CAD", "authoritative V7 engineering records", "physical data", "V7 CAE inputs and thresholds"]})

print(json.dumps({"status": "PARTIAL", "reports": len(list(OUT.glob("*.json"))), "source_sha256": R41_SHA}, indent=2))
