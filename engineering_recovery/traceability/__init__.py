from .graph import EvidenceEdge, EvidenceGraph, EvidenceNode
from .gap_mapping import build_gap_mapping, gap_ids


def normalized_joint(joint_id: str, body_a: str, body_b: str, source_evidence: dict | None = None) -> dict:
    return {
        "joint_id": joint_id,
        "body_a": body_a,
        "body_b": body_b,
        "joint_type": "UNDEFINED",
        "location": "UNDEFINED",
        "axis": "UNDEFINED",
        "limits": "UNDEFINED",
        "source_evidence": source_evidence or {},
        "status": "UNVERIFIED",
    }


__all__ = ["EvidenceEdge", "EvidenceGraph", "EvidenceNode", "normalized_joint", "build_gap_mapping", "gap_ids"]
