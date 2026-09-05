from __future__ import annotations

import hashlib
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bridge.r41_paths import (  # noqa: E402
    BODY_MAPPING,
    CONTACT_MAPPING,
    EXPECTED_R41_SHA256,
    JOINT_MAPPING,
    PARAMETERS,
    R41_STEP,
)
from model.r41_model import import_static_model  # noqa: E402


class R41StaticIntegrityTests(unittest.TestCase):
    def test_r41_sha_is_immutable(self) -> None:
        digest = hashlib.sha256(R41_STEP.read_bytes()).hexdigest()
        self.assertEqual(digest, EXPECTED_R41_SHA256)

    def test_mapping_is_31_records_with_62_solid_coverage(self) -> None:
        data = json.loads(BODY_MAPPING.read_text(encoding="utf-8"))
        bodies = data["bodies"]
        self.assertEqual(len(bodies), 31)
        solids = [solid for body in bodies for solid in body["step_solid_indices"]]
        self.assertEqual(sorted(solids), list(range(1, 63)))
        self.assertEqual(len(set(solids)), 62)

    def test_joint_and_contact_references_are_explicit(self) -> None:
        joints = json.loads(JOINT_MAPPING.read_text(encoding="utf-8"))["joints"]
        contacts = json.loads(CONTACT_MAPPING.read_text(encoding="utf-8"))["contacts"]
        self.assertTrue(all(row.get("body_A") and row.get("body_B") for row in joints))
        self.assertTrue(all(row.get("body_A") and row.get("body_B") for row in contacts))

    def test_static_model_import(self) -> None:
        model = import_static_model()
        self.assertEqual(model.step_sha256, EXPECTED_R41_SHA256)
        self.assertEqual(model.body_records, 31)
        self.assertEqual(model.solid_coverage, tuple(range(1, 63)))

    def test_physical_unknowns_are_not_silent_defaults(self) -> None:
        parameters = json.loads(PARAMETERS.read_text(encoding="utf-8"))["parameters"]
        allowed = {"UNKNOWN", "ASSUMED", "SOURCE_VERIFIED", "DERIVED"}
        self.assertTrue(all(row.get("status") in allowed for row in parameters))


if __name__ == "__main__":
    unittest.main()
