"""Evidence graph primitives with explicit provenance on nodes and edges."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from typing import Any

from engineering_recovery.evidence.contracts import EvidenceClass


@dataclass(frozen=True)
class EvidenceNode:
    node_id: str
    node_type: str
    label: str
    source: str | None
    source_sha256: str | None
    evidence_class: EvidenceClass
    verification_state: str


@dataclass(frozen=True)
class EvidenceEdge:
    edge_id: str
    source_node: str
    target_node: str
    source: str | None
    source_sha256: str | None
    evidence_class: EvidenceClass
    timestamp: str
    tool: str | None
    tool_version: str | None
    confidence: str
    verification_state: str


class EvidenceGraph:
    def __init__(self) -> None:
        self.nodes: dict[str, EvidenceNode] = {}
        self.edges: list[EvidenceEdge] = []

    def add_node(self, node: EvidenceNode) -> None:
        self.nodes[node.node_id] = node

    def add_edge(self, edge: EvidenceEdge) -> None:
        if edge.source_node not in self.nodes or edge.target_node not in self.nodes:
            raise ValueError("evidence edge references an unknown node")
        if edge.evidence_class in {EvidenceClass.DERIVED_FROM_CAD, EvidenceClass.CALCULATED, EvidenceClass.SIMULATED} and not edge.source_sha256:
            raise ValueError("derived/calculated/simulated edges require source SHA-256")
        self.edges.append(edge)

    def orphan_nodes(self) -> list[str]:
        linked = {e.source_node for e in self.edges} | {e.target_node for e in self.edges}
        return sorted(set(self.nodes) - linked)

    def as_dict(self) -> dict[str, Any]:
        return {"nodes": [asdict(n) | {"evidence_class": n.evidence_class.value} for n in self.nodes.values()],
                "edges": [asdict(e) | {"evidence_class": e.evidence_class.value} for e in self.edges],
                "orphan_nodes": self.orphan_nodes(), "generated_at": datetime.now(timezone.utc).isoformat()}
