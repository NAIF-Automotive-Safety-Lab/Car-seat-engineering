from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from ..acquisition.service import sha256_file, utc_now

ENTITY_RE = re.compile(r"#(\d+)\s*=\s*([A-Z0-9_]+)", re.MULTILINE)
HEADER_RE = re.compile(r"HEADER;(?P<body>.*?)ENDSEC;", re.DOTALL | re.IGNORECASE)


def _header_fields(text: str) -> dict[str, Any]:
    match = HEADER_RE.search(text)
    if not match:
        return {"present": False, "raw": ""}
    body = match.group("body")
    return {
        "present": True,
        "file_description": re.findall(r"FILE_DESCRIPTION\s*\((.*?)\);", body, re.DOTALL | re.IGNORECASE),
        "file_name": re.findall(r"FILE_NAME\s*\((.*?)\);", body, re.DOTALL | re.IGNORECASE),
        "file_schema": re.findall(r"FILE_SCHEMA\s*\((.*?)\);", body, re.DOTALL | re.IGNORECASE),
    }


def parse_step_text(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    text = raw.decode("utf-8", errors="replace")
    entities = ENTITY_RE.findall(text)
    counts: dict[str, int] = {}
    for _, name in entities:
        counts[name] = counts.get(name, 0) + 1
    required_markers = {
        "iso_10303_21": text.startswith("ISO-10303-21;"),
        "end_iso_10303_21": text.rstrip().endswith("END-ISO-10303-21;"),
    }
    geometry_names = ("CARTESIAN_POINT", "DIRECTION", "AXIS2_PLACEMENT_3D", "CYLINDRICAL_SURFACE", "PLANE", "CIRCLE", "LINE")
    topology_names = ("MANIFOLD_SOLID_BREP", "CLOSED_SHELL", "ADVANCED_FACE", "EDGE_CURVE", "VERTEX_POINT")
    pmi_names = tuple(name for name in counts if any(token in name for token in ("TOLERANCE", "DATUM", "DIMENSION", "GEOMETRIC")))
    return {
        "source_sha256": sha256_file(path),
        "source_size_bytes": len(raw),
        "extraction_timestamp": utc_now(),
        "parser": "EERE_RAW_STEP_ENTITY_SCANNER/0.1.0",
        "status": "VERIFIED" if all(required_markers.values()) else "BLOCKED",
        "header": _header_fields(text),
        "schema": counts.get("FILE_SCHEMA", 0),
        "application_protocol": "AUTOMOTIVE_DESIGN" if "AUTOMOTIVE_DESIGN" in text else "NOT_DETECTED",
        "entity_count": len(entities),
        "entity_counts": counts,
        "product_definitions": counts.get("PRODUCT_DEFINITION", 0),
        "product_definition_relationships": counts.get("PRODUCT_DEFINITION_RELATIONSHIP", 0),
        "representation_relationships": counts.get("SHAPE_REPRESENTATION_RELATIONSHIP", 0),
        "geometry_entity_count": sum(counts.get(n, 0) for n in geometry_names),
        "topology_entity_count": sum(counts.get(n, 0) for n in topology_names),
        "pmi_entity_count": sum(counts.get(n, 0) for n in pmi_names),
        "pmi_entities": {n: counts[n] for n in pmi_names},
        "assembly_entities": {n: counts[n] for n in counts if "ASSEMBLY" in n or "NEXT_ASSEMBLY" in n},
        "unknown_entities": [],
        "parse_errors": [] if all(required_markers.values()) else ["INVALID_ISO_10303_21_BOUNDARIES"],
        "warnings": ["Raw scanner is semantic-light; compare with OCCT report."],
    }


def import_with_occt(path: Path) -> tuple[Any, dict[str, Any]]:
    try:
        from OCP.STEPControl import STEPControl_Reader
    except Exception as exc:
        return None, {"status": "BLOCKED", "reason": f"OCCT_UNAVAILABLE:{type(exc).__name__}:{exc}"}
    reader = STEPControl_Reader()
    status = reader.ReadFile(str(path))
    if "RetDone" not in str(status):
        return None, {"status": "BLOCKED", "reason": f"STEP_READ_FAILED:{status}"}
    reader.TransferRoots()
    return reader.OneShape(), {"status": "VERIFIED", "reader_status": str(status), "kernel": "OCP/OCCT"}


def build_step_forensic_report(path: Path) -> dict[str, Any]:
    report = parse_step_text(path)
    _, occt = import_with_occt(path)
    report["occt"] = occt
    report["dual_parse_status"] = "VERIFIED" if report["status"] == "VERIFIED" and occt.get("status") == "VERIFIED" else "BLOCKED"
    return report
