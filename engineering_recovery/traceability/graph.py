from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class EvidenceGraph:
    nodes: dict[str, dict[str, Any]] = field(default_factory=dict)
    edges: list[dict[str, Any]] = field(default_factory=list)

    def add_node(self, node_id: str, node_type: str, **data: Any) -> None:
        self.nodes[node_id] = {"id": node_id, "type": node_type, **data}

    def add_edge(self, source: str, target: str, evidence: dict[str, Any]) -> None:
        if source not in self.nodes or target not in self.nodes:
            raise ValueError("TRACEABILITY_EDGE_ENDPOINT_MISSING")
        if not evidence:
            raise ValueError("TRACEABILITY_EDGE_EVIDENCE_REQUIRED")
        self.edges.append({"source": source, "target": target, "evidence": evidence})

    def orphan_audit(self) -> dict[str, list[str]]:
        connected = {x for edge in self.edges for x in (edge["source"], edge["target"])}
        return {"orphan_nodes": sorted(set(self.nodes) - connected), "orphan_edges": []}

    def as_dict(self) -> dict[str, Any]:
        return {"nodes": list(self.nodes.values()), "edges": self.edges, "orphan_audit": self.orphan_audit()}


def normalized_joint(joint_id: str, body_a: str, body_b: str, source_evidence: dict[str, Any] | None = None) -> dict[str, Any]:
    return {
        "joint_id": joint_id, "body_a": body_a, "body_b": body_b,
        "joint_type": "UNDEFINED", "location": "UNDEFINED", "axis": "UNDEFINED",
        "limits": "UNDEFINED", "clearance": "UNDEFINED", "stop_condition": "UNDEFINED",
        "fastener_reference": "UNDEFINED", "bushing_reference": "UNDEFINED",
        "compliance_reference": "UNDEFINED", "source_evidence": source_evidence or {},
        "confidence": "UNDEFINED", "status": "UNVERIFIED",
    }
