from __future__ import annotations

from pathlib import Path
from typing import Any

from ..acquisition.service import sha256_file, utc_now
from ..step.forensics import import_with_occt


def _count(shape: Any, kind: Any) -> int:
    from OCP.TopExp import TopExp_Explorer
    explorer = TopExp_Explorer(shape, kind)
    total = 0
    while explorer.More():
        total += 1
        explorer.Next()
    return total


def validate_brep(path: Path, extractor_version: str = "eere-brep/0.1.0") -> dict[str, Any]:
    shape, readback = import_with_occt(path)
    base = {
        "source_artifact_sha256": sha256_file(path),
        "extraction_timestamp": utc_now(),
        "cad_kernel_version": "OCP/OCCT",
        "extraction_tool_version": extractor_version,
        "deterministic": True,
        "status": "BLOCKED",
        "readback": readback,
    }
    if shape is None:
        base["reason"] = readback.get("reason", "OCCT_IMPORT_BLOCKED")
        return base
    from OCP.TopAbs import TopAbs_COMPSOLID, TopAbs_COMPOUND, TopAbs_EDGE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box(); BRepBndLib.Add_s(shape, box)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    result = dict(base)
    result.update({
        "solid_count": _count(shape, TopAbs_SOLID),
        "shell_count": _count(shape, TopAbs_SHELL),
        "face_count": _count(shape, TopAbs_FACE),
        "edge_count": _count(shape, TopAbs_EDGE),
        "vertex_count": _count(shape, TopAbs_VERTEX),
        "compound_count": _count(shape, TopAbs_COMPOUND),
        "compsolid_count": _count(shape, TopAbs_COMPSOLID),
        "bounding_box": {"xmin": xmin, "ymin": ymin, "zmin": zmin, "xmax": xmax, "ymax": ymax, "zmax": zmax},
        "volume": "UNDEFINED_UNTIL_DETERMINISTIC_MASS_PROPERTIES_POLICY",
        "surface_area": "UNDEFINED_UNTIL_DETERMINISTIC_MASS_PROPERTIES_POLICY",
        "center_of_mass": "UNDEFINED_UNTIL_DETERMINISTIC_MASS_PROPERTIES_POLICY",
        "shape_validity": "IMPORTED; DETAILED_CHECK_NOT_RUN",
        "self_intersection_checks": "NOT_RUN",
        "status": "VERIFIED",
    })
    return result
