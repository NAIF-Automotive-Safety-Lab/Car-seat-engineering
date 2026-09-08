from __future__ import annotations

from pathlib import Path
from typing import Any

from ..acquisition.service import sha256_file, utc_now
from ..step.forensics import parse_step_text


def extract_features(path: Path) -> dict[str, Any]:
    report = parse_step_text(path)
    counts = report["entity_counts"]
    return {
        "source_artifact_sha256": sha256_file(path),
        "extraction_timestamp": utc_now(),
        "status": "CANDIDATE_ONLY",
        "deterministic": True,
        "features": {
            "planes": counts.get("PLANE", 0),
            "cylinders": counts.get("CYLINDRICAL_SURFACE", 0),
            "cones": counts.get("CONICAL_SURFACE", 0),
            "spheres": counts.get("SPHERICAL_SURFACE", 0),
            "holes": "UNDEFINED_REQUIRES_TOPOLOGY_RULES",
            "slots": "UNDEFINED_REQUIRES_TOPOLOGY_RULES",
            "fillets": "UNDEFINED_REQUIRES_TOPOLOGY_RULES",
            "chamfers": "UNDEFINED_REQUIRES_TOPOLOGY_RULES",
            "candidate_pivot_axes": "UNDEFINED_UNTIL_ASSEMBLY_CORRELATION",
            "candidate_hinge_axes": "UNDEFINED_UNTIL_ASSEMBLY_CORRELATION",
            "candidate_stop_faces": "UNDEFINED_UNTIL_FEATURE_RULES",
            "mating_surfaces": "UNDEFINED_UNTIL_INTERFACE_RULES",
            "clearance_candidates": "UNDEFINED_UNTIL_INTERFACE_RULES",
            "datum_candidates": "UNDEFINED_UNTIL_PMI_OR_RULES",
        },
        "authority_rule": "AI_OR_HEURISTIC_FEATURES_CANNOT_BE_AUTHORITATIVE_WITHOUT_INDEPENDENT_VERIFICATION",
    }


def pmi_report(path: Path) -> dict[str, Any]:
    report = parse_step_text(path)
    status = "PRESENT" if report["pmi_entity_count"] else "NOT_PRESENT_OR_NOT_RECOVERED"
    return {"source_artifact_sha256": sha256_file(path), "extraction_timestamp": utc_now(), "PMI_STATUS": status, "items": [], "note": "No PMI is inferred."}
