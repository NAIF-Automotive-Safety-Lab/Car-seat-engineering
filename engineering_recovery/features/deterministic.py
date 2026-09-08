from __future__ import annotations

from pathlib import Path
from typing import Any

from engineering_recovery.acquisition.service import sha256_file


def feature_report(path: Path) -> dict[str, Any]:
    """Inventory deterministic face/topology candidates; never labels design intent as fact."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Plane, GeomAbs_Sphere

    reader = STEPControl_Reader()
    status = reader.ReadFile(str(path))
    if "RetDone" not in str(status):
        raise RuntimeError(f"STEP_READ_FAILED:{status}")
    reader.TransferRoots()
    shape = reader.OneShape()
    counts = {"plane": 0, "cylinder": 0, "cone": 0, "sphere": 0, "other_surface": 0}
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        kind = BRepAdaptor_Surface(TopoDS.Face(explorer.Current())).GetType()
        if kind == GeomAbs_Plane:
            counts["plane"] += 1
        elif kind == GeomAbs_Cylinder:
            counts["cylinder"] += 1
        elif kind == GeomAbs_Cone:
            counts["cone"] += 1
        elif kind == GeomAbs_Sphere:
            counts["sphere"] += 1
        else:
            counts["other_surface"] += 1
        explorer.Next()
    solids = TopExp_Explorer(shape, TopAbs_SOLID)
    solid_count = 0
    while solids.More():
        solid_count += 1
        solids.Next()
    return {
        "status": "VERIFIED",
        "source_sha256": sha256_file(path),
        "kernel": "OCP/OpenCascade",
        "deterministic": True,
        "solid_count": solid_count,
        "surface_type_counts": counts,
        "candidate_features": {
            "planar_faces": counts["plane"],
            "cylindrical_faces": counts["cylinder"],
            "conical_faces": counts["cone"],
            "spherical_faces": counts["sphere"],
            "holes": "UNDEFINED",
            "slots": "UNDEFINED",
            "fillets": "UNDEFINED",
            "chamfers": "UNDEFINED",
            "interfaces": "UNDEFINED",
            "pivot_axes": "UNDEFINED",
            "hinge_axes": "UNDEFINED",
        },
        "authority_rule": "surface counts are deterministic; semantic feature labels remain UNDEFINED unless independently verified",
    }


__all__ = ["feature_report"]
