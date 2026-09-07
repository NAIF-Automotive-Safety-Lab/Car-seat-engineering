"""Fail-closed engineering record contracts for non-CAD evidence domains."""
from __future__ import annotations

from typing import Any

from engineering_recovery.evidence.contracts import EvidenceClass, GapStatus, make_blocked_record


def joint_record(joint_id: str, body_a: str, body_b: str, joint_type: str) -> dict[str, Any]:
    return {"joint_id": joint_id, "body_a": body_a, "body_b": body_b, "joint_type": joint_type,
            "location": None, "axis": None, "limits": None, "clearance": None,
            "stop_condition": None, "fastener_reference": None, "bushing_reference": None,
            "compliance_reference": None, "source_evidence": EvidenceClass.BLOCKED.value,
            "confidence": "NONE", "status": GapStatus.BLOCKED.value,
            "reason": "No authoritative joint geometry or test evidence"}


def material_record(part_reference: str) -> dict[str, Any]:
    fields = ["material", "grade", "heat_treatment", "density", "youngs_modulus", "poisson_ratio",
              "yield_strength", "ultimate_strength", "plasticity", "strain_rate", "temperature_dependence", "fatigue"]
    return {"part_reference": part_reference, "status": GapStatus.BLOCKED.value,
            "fields": [make_blocked_record(field, "Provide certificate, datasheet, standard, or test", "No authoritative material evidence") for field in fields]}


def fastener_record(part_reference: str) -> dict[str, Any]:
    fields = ["type", "diameter", "thread", "grade", "quantity", "location", "washer", "nut", "torque", "preload", "friction", "joint_stiffness"]
    return {"part_reference": part_reference, "status": GapStatus.BLOCKED.value,
            "fields": [make_blocked_record(field, "Provide released BOM or fastener record", "No authoritative fastener evidence") for field in fields]}


def pmi_record() -> dict[str, Any]:
    return {"status": "NOT_PRESENT_OR_NOT_RECOVERED", "dimensions": [], "tolerances": [], "gd_t": [], "datums": [], "surface_finish": [], "material_callouts": [], "assembly_notes": [], "manufacturing_notes": []}


def feature_candidate(feature_type: str, entity_ref: str) -> dict[str, Any]:
    return {"feature_type": feature_type, "entity_ref": entity_ref, "status": "CANDIDATE", "evidence_class": EvidenceClass.CANDIDATE.value, "deterministic_verification": "REQUIRED", "authoritative": False}
