from __future__ import annotations

import json
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
R41_STEP = PROJECT_ROOT / "R4.1" / "R4.1.step"
SR11_ROOT = PROJECT_ROOT / "SR11"
BODY_MAPPING = SR11_ROOT / "manifests" / "V1_SR11_BODY_MAPPING.json"
JOINT_MAPPING = SR11_ROOT / "joints" / "V1_SR11_JOINT_MAPPING.json"
CONTACT_MAPPING = SR11_ROOT / "contacts" / "V1_SR11_CONTACT_MAPPING.json"
PARAMETERS = SR11_ROOT / "parameters" / "V1_SR11_ASSUMPTIONS.json"
EXPECTED_R41_SHA256 = "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(str(path))
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON_ROOT_NOT_OBJECT:{path}")
    return value


def require_authoritative_inputs() -> None:
    required = (R41_STEP, BODY_MAPPING, JOINT_MAPPING, CONTACT_MAPPING, PARAMETERS)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("MISSING_AUTHORITATIVE_INPUTS:" + ",".join(missing))
