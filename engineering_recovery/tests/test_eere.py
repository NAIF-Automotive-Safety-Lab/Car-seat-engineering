from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engineering_recovery.acquisition import ArtifactAcquisitionService  # noqa: E402
from engineering_recovery.brep import validate_brep  # noqa: E402
from engineering_recovery.features import extract_features, pmi_report  # noqa: E402
from engineering_recovery.runtime import chrono_smoke  # noqa: E402
from engineering_recovery.step import build_step_forensic_report  # noqa: E402
from engineering_recovery.traceability import EvidenceGraph  # noqa: E402


class EereTests(unittest.TestCase):
    def test_artifact_sha_reproducibility_and_immutable_copy(self) -> None:
        source = ROOT / "R4.1/R4.1.step"
        expected = "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"
        with tempfile.TemporaryDirectory() as tmp:
            service = ArtifactAcquisitionService(ROOT, Path(tmp))
            destination, manifest = service.immutable_copy(source, "test-r41")
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), expected)
            self.assertEqual(manifest.sha256, expected)
            self.assertEqual(hashlib.sha256(destination.read_bytes()).hexdigest(), expected)
            self.assertTrue(service.replay_artifact(destination.parent / "manifest.json").sha256 == expected)

    def test_dual_step_forensics_is_verified(self) -> None:
        report = build_step_forensic_report(ROOT / "R4.1/R4.1.step")
        self.assertEqual(report["dual_parse_status"], "VERIFIED")
        self.assertEqual(report["source_sha256"], "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68")
        self.assertGreater(report["entity_count"], 0)

    def test_brep_solid_count_is_reproducible(self) -> None:
        report = validate_brep(ROOT / "R4.1/R4.1.step")
        self.assertEqual(report["status"], "VERIFIED")
        self.assertEqual(report["solid_count"], 62)
        self.assertTrue(report["deterministic"])

    def test_feature_and_pmi_firewall(self) -> None:
        path = ROOT / "R4.1/R4.1.step"
        features = extract_features(path)
        pmi = pmi_report(path)
        self.assertEqual(features["status"], "CANDIDATE_ONLY")
        self.assertIn(pmi["PMI_STATUS"], {"PRESENT", "NOT_PRESENT_OR_NOT_RECOVERED"})
        self.assertEqual(features["features"]["holes"], "UNDEFINED_REQUIRES_TOPOLOGY_RULES")

    def test_traceability_requires_evidence_and_reports_orphans(self) -> None:
        graph = EvidenceGraph()
        graph.add_node("a", "PARAMETER")
        graph.add_node("b", "CAD_ENTITY")
        graph.add_node("orphan", "JOINT")
        graph.add_edge("a", "b", {"sha256": "x"})
        self.assertEqual(graph.orphan_audit()["orphan_nodes"], ["orphan"])
        with self.assertRaises(ValueError):
            graph.add_edge("a", "missing", {"sha256": "x"})

    def test_chrono_is_honest(self) -> None:
        result = chrono_smoke()
        self.assertIn(result["status"], {"VERIFIED", "BLOCKED"})
        if result["status"] == "BLOCKED":
            self.assertIn("error", result)


if __name__ == "__main__":
    unittest.main()
