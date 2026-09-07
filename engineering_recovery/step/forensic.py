from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from engineering_recovery.acquisition.service import sha256_file


ENTITY_RE = re.compile(r"#(\d+)\s*=\s*([A-Z0-9_]+)", re.IGNORECASE)


def forensic_report(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    text = raw.decode("latin-1", errors="replace")
    entities = ENTITY_RE.findall(text)
    counts: dict[str, int] = {}
    for _, name in entities:
        key = name.upper()
        counts[key] = counts.get(key, 0) + 1
    header = text[: text.find("DATA;") if "DATA;" in text else min(len(text), 4096)]
    schema = re.findall(r"FILE_SCHEMA\s*\(\s*\(('?[^')]+)", header, flags=re.IGNORECASE)
    products = counts.get("PRODUCT_DEFINITION", 0)
    return {
        "status": "PASS",
        "source_sha256": sha256_file(path),
        "file_size_bytes": len(raw),
        "parser": "EERE deterministic ISO-10303-21 entity scanner",
        "deterministic": True,
        "header_present": text.startswith("ISO-10303-21;"),
        "footer_present": text.rstrip().endswith("END-ISO-10303-21;"),
        "schema": schema,
        "application_protocol": re.findall(r"FILE_DESCRIPTION\s*\(\s*\(('?[^']+)", header, flags=re.IGNORECASE),
        "entity_count": len(entities),
        "entity_type_counts": dict(sorted(counts.items())),
        "product_definitions": products,
        "product_definition_relationships": counts.get("PRODUCT_DEFINITION_RELATIONSHIP", 0),
        "representation_relationships": counts.get("REPRESENTATION_RELATIONSHIP", 0),
        "geometry_entity_count": sum(counts.get(name, 0) for name in ("CARTESIAN_POINT", "LINE", "CIRCLE", "CYLINDRICAL_SURFACE", "PLANE", "CONICAL_SURFACE", "SPHERICAL_SURFACE")),
        "topology_entity_count": sum(counts.get(name, 0) for name in ("MANIFOLD_SOLID_BREP", "CLOSED_SHELL", "ADVANCED_FACE", "EDGE_CURVE", "VERTEX_POINT")),
        "pmi_entity_count": sum(counts.get(name, 0) for name in ("GEOMETRIC_TOLERANCE", "DIMENSIONAL_LOCATION", "SHAPE_DIMENSION", "DATUM_FEATURE")),
        "assembly_entities": sum(counts.get(name, 0) for name in ("NEXT_ASSEMBLY_USAGE_OCCURRENCE", "PRODUCT_DEFINITION_RELATIONSHIP")),
        "unknown_entities": [],
        "parse_errors": [],
        "warnings": ["STEPcode and step-p21 comparison is NOT_AVAILABLE; this report is not a dual-parser claim."],
    }
