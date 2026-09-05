from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bridge.r41_paths import (
    BODY_MAPPING,
    CONTACT_MAPPING,
    EXPECTED_R41_SHA256,
    JOINT_MAPPING,
    PARAMETERS,
    R41_STEP,
    load_json,
    require_authoritative_inputs,
)


@dataclass(frozen=True)
class R41StaticModel:
    step_path: Path
    step_sha256: str
    body_records: int
    solid_coverage: tuple[int, ...]
    joints: tuple[dict[str, Any], ...]
    contacts: tuple[dict[str, Any], ...]
    parameters: tuple[dict[str, Any], ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def import_static_model() -> R41StaticModel:
    require_authoritative_inputs()
    actual = _sha256(R41_STEP)
    if actual != EXPECTED_R41_SHA256:
        raise ValueError(f"R4.1_SHA_MISMATCH:{actual}")
    body_map = load_json(BODY_MAPPING)
    joint_map = load_json(JOINT_MAPPING)
    contact_map = load_json(CONTACT_MAPPING)
    parameters = load_json(PARAMETERS)
    bodies = body_map.get("bodies")
    joints = joint_map.get("joints")
    contacts = contact_map.get("contacts")
    parameter_rows = parameters.get("parameters")
    if not isinstance(bodies, list) or not isinstance(joints, list) or not isinstance(contacts, list) or not isinstance(parameter_rows, list):
        raise ValueError("STATIC_MODEL_MAPPING_SCHEMA_INVALID")
    solids = [solid for body in bodies for solid in body.get("step_solid_indices", [])]
    if sorted(solids) != list(range(1, 63)) or len(set(solids)) != 62:
        raise ValueError("STATIC_MODEL_SOLID_COVERAGE_INVALID")
    if any(not row.get("body_A") or not row.get("body_B") for row in joints + contacts):
        raise ValueError("STATIC_MODEL_REFERENCE_INVALID")
    return R41StaticModel(
        step_path=R41_STEP,
        step_sha256=actual,
        body_records=len(bodies),
        solid_coverage=tuple(sorted(solids)),
        joints=tuple(joints),
        contacts=tuple(contacts),
        parameters=tuple(parameter_rows),
    )


def import_step_with_ocp() -> Any:
    model = import_static_model()
    try:
        from OCP.STEPControl import STEPControl_Reader
    except Exception as exc:
        raise RuntimeError(f"OCP_UNAVAILABLE:{type(exc).__name__}:{exc}") from exc
    reader = STEPControl_Reader()
    status = reader.ReadFile(str(model.step_path))
    if "RetDone" not in str(status):
        raise RuntimeError(f"STEP_READ_FAILED:{status}")
    reader.TransferRoots()
    return reader.OneShape()
