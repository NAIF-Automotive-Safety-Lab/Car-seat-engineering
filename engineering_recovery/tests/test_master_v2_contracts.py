from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from engineering_recovery.evidence.contracts import EvidenceClass, EvidenceRef, reject_physical_claim
from engineering_recovery.instrumentation.interfaces import MeasurementRecord
from engineering_recovery.physics.interfaces import absorber_record, lock_170_record, rebound_180_record
from engineering_recovery.gap_closure.workflow import GapStatus, close_gap, new_blocked_gap
from engineering_recovery.interfaces.contracts import feature_candidate
from engineering_recovery.traceability.graph import EvidenceEdge, EvidenceGraph, EvidenceNode


def test_blocked_physics_never_selects_lock() -> None:
    record = lock_170_record("L170-A")
    assert record["final_selection"] == "NOT_AUTHORIZED"
    assert record["status"] == "BLOCKED"


def test_rebound_target_is_not_validation() -> None:
    record = rebound_180_record()
    assert record["target_180_mm"]["evidence_class"] == "TARGET"
    assert record["target_180_mm"]["status"] == "CANDIDATE"


def test_absorber_requires_measured_curves() -> None:
    assert absorber_record()["absorber_status"] == "UNVERIFIED"


def test_measurement_requires_raw_hash() -> None:
    with pytest.raises(ValueError):
        MeasurementRecord("ACC-1", "CAL-1", 1000, "m/s2", "now", "TEST-1", "bad")


def test_physical_claim_rejects_missing_provenance() -> None:
    with pytest.raises(ValueError):
        reject_physical_claim(EvidenceRef(None, None, EvidenceClass.MEASURED_PHYSICALLY, "tool", "1", "VERIFIED", "now"))


def test_gap_cannot_close_without_evidence() -> None:
    gap = new_blocked_gap("L01", "containment", "provide test evidence")
    with pytest.raises(ValueError):
        close_gap(gap, (), GapStatus.VERIFIED)


def test_candidate_is_not_authoritative() -> None:
    assert feature_candidate("plane", "entity-1")["authoritative"] is False


def test_graph_requires_source_hash_for_derived_edge() -> None:
    graph = EvidenceGraph()
    graph.add_node(EvidenceNode("a", "CAD", "source", "x.step", "a" * 64, EvidenceClass.DERIVED_FROM_CAD, "VERIFIED"))
    graph.add_node(EvidenceNode("b", "FEATURE", "feature", None, None, EvidenceClass.CANDIDATE, "CANDIDATE"))
    with pytest.raises(ValueError):
        graph.add_edge(EvidenceEdge("e", "a", "b", "x.step", None, EvidenceClass.DERIVED_FROM_CAD, "now", "tool", "1", "LOW", "CANDIDATE"))
