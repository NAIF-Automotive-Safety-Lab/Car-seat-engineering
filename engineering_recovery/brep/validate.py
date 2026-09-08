from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from engineering_recovery.acquisition.service import sha256_file


def _box(shape: Any) -> dict[str, float]:
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    xmin, ymin, zmin = box.GetXMin(), box.GetYMin(), box.GetZMin()
    xmax, ymax, zmax = box.GetXMax(), box.GetYMax(), box.GetZMax()
    return {"xmin": xmin, "ymin": ymin, "zmin": zmin, "xmax": xmax, "ymax": ymax, "zmax": zmax}


def brep_report(path: Path) -> dict[str, Any]:
    from OCP.BRep import BRep_Tool
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_COMPOUND, TopAbs_EDGE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX
    from OCP.TopExp import TopExp_Explorer
    from OCP.GProp import GProp_GProps
    from OCP.BRepGProp import BRepGProp

    reader = STEPControl_Reader()
    status = reader.ReadFile(str(path))
    if "RetDone" not in str(status):
        raise RuntimeError(f"STEP_READ_FAILED:{status}")
    reader.TransferRoots()
    shape = reader.OneShape()

    def count(kind: Any) -> int:
        explorer = TopExp_Explorer(shape, kind)
        total = 0
        while explorer.More():
            total += 1
            explorer.Next()
        return total

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    center = props.CentreOfMass()
    return {
        "status": "PASS",
        "source_sha256": sha256_file(path),
        "extraction_timestamp": datetime.now(timezone.utc).isoformat(),
        "kernel": "OCP/OpenCascade",
        "kernel_version": "cadquery-ocp 8.0.1.0.0 (environment evidence)",
        "deterministic": True,
        "read_status": str(status),
        "solid_count": count(TopAbs_SOLID),
        "shell_count": count(TopAbs_SHELL),
        "face_count": count(TopAbs_FACE),
        "edge_count": count(TopAbs_EDGE),
        "vertex_count": count(TopAbs_VERTEX),
        "compound_count": count(TopAbs_COMPOUND),
        "bounding_box": _box(shape),
        "volume": float(props.Mass()),
        "surface_area": None,
        "center_of_mass": {"x": float(center.X()), "y": float(center.Y()), "z": float(center.Z())},
        "shape_validity": bool(BRepCheck_Analyzer(shape).IsValid()),
        "self_intersection_check": "NOT_PERFORMED",
        "notes": ["Surface-area and surface/curve type inventories are not claimed until a dedicated deterministic extractor is added."],
    }
