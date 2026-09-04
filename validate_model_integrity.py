#!/usr/bin/env python3
"""Fail-closed validation for the R4.1 SR11 model manifest."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "SR11/manifests/V1_SR11_BODY_MAPPING.json"

REQUIRED_BODY_IDS = {
    "BASE",
    "RAIL_L",
    "RAIL_R",
    "CARRIAGE",
    "ABSORBER_L",
    "ABSORBER_R",
    "SEATBACK",
    "ROT_CTRL_L",
    "ROT_CTRL_R",
    "LOCK_L",
    "LOCK_R",
    "ENDSTOPS",
    "HINGE",
    "ANTI_RACK",
    "SHOULDER_ANCHOR",
    "SHOULDER_GUSSET",
    "SHOULDER_ROOT",
    "SHOULDER_BRACE",
    "SHOULDER_FRAME_BRACKET",
}
EXPECTED_SOLID_INDICES = set(range(1, 63))


def fail(reason: str, **extra: object) -> int:
    print("MODEL_INTEGRITY=FAIL")
    print(f"REASON={reason}")
    for key, value in extra.items():
        print(f"{key}={value}")
    return 1


def main() -> int:
    if not MANIFEST.is_file():
        return fail("MISSING_BODY_MAPPING", path=str(MANIFEST))

    try:
        data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except Exception as exc:
        return fail("BODY_MAPPING_JSON_INVALID", error=repr(exc))

    bodies = data.get("bodies")
    if not isinstance(bodies, list) or not bodies:
        return fail("NO_BODIES")

    body_ids = [b.get("mbd_body_id") for b in bodies]
    if any(not isinstance(x, str) or not x for x in body_ids):
        return fail("EMPTY_MBD_BODY_ID")

    if len(set(body_ids)) != len(body_ids):
        return fail("DUPLICATE_MBD_BODY_ID")

    missing_required = sorted(REQUIRED_BODY_IDS - set(body_ids))
    if missing_required:
        return fail("MISSING_REQUIRED_BODY_IDS", missing=",".join(missing_required))

    step_indices: list[int] = []
    for body in bodies:
        indices = body.get("step_solid_indices", [])
        if not isinstance(indices, list):
            return fail("INVALID_STEP_SOLID_INDICES", body=body.get("mbd_body_id"))
        for idx in indices:
            if not isinstance(idx, int):
                return fail("NON_INTEGER_STEP_SOLID_INDEX", body=body.get("mbd_body_id"))
            step_indices.append(idx)

    if set(step_indices) != EXPECTED_SOLID_INDICES:
        missing = sorted(EXPECTED_SOLID_INDICES - set(step_indices))
        extra = sorted(set(step_indices) - EXPECTED_SOLID_INDICES)
        return fail(
            "STEP_SOLID_COVERAGE_NOT_1_TO_62",
            missing=missing,
            extra=extra,
            count=len(step_indices),
            unique=len(set(step_indices)),
        )

    if len(step_indices) != 62 or len(set(step_indices)) != 62:
        return fail(
            "STEP_SOLID_MAPPING_NOT_UNIQUE",
            count=len(step_indices),
            unique=len(set(step_indices)),
        )

    print("SOLID_MAPPING_COVERAGE=62/62")
    print("UNIQUE_MAPPING=PASS")
    print("REQUIRED_MODEL_ELEMENTS=PASS")
    print("MODEL_INTEGRITY=PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
