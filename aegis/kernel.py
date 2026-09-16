from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import hashlib
import json
from typing import Any, Iterable


class Decision(str, Enum):
    PASS = "PASS"
    BLOCK = "BLOCK"
    UNKNOWN = "UNKNOWN"


class Phase(str, Enum):
    REQUIREMENTS = "REQUIREMENTS"
    ARCHITECTURE = "ARCHITECTURE"
    CAD = "CAD"
    CAE = "CAE"
    DFM = "DFM"
    TEST = "TEST"
    REVIEW = "REVIEW"
    COMPLIANCE = "COMPLIANCE"
    MANUFACTURING = "MANUFACTURING"
    RELEASE = "RELEASE"


@dataclass(frozen=True)
class Artifact:
    artifact_id: str
    sha256: str
    kind: str
    uri: str
    producer: str
    signed_by: str | None = None

    def __post_init__(self) -> None:
        if not self.artifact_id or not self.kind or not self.uri or not self.producer:
            raise ValueError("artifact identity, type, URI, and producer are required")
        if len(self.sha256) != 64 or any(c not in "0123456789abcdef" for c in self.sha256.lower()):
            raise ValueError("artifact sha256 must be 64 hexadecimal characters")


@dataclass(frozen=True)
class Event:
    sequence: int
    phase: Phase
    decision: Decision
    reason: str
    artifact_ids: tuple[str, ...]
    previous_hash: str
    event_hash: str
    graph_hash: str


class EvidenceGraph:
    """Append-only, content-addressed evidence graph."""
    def __init__(self) -> None:
        self._artifacts: dict[str, Artifact] = {}
        self._links: set[tuple[str, str, str]] = set()

    def add_artifact(self, artifact: Artifact) -> None:
        prior = self._artifacts.get(artifact.artifact_id)
        if prior and prior != artifact:
            raise ValueError("artifact ID is immutable and cannot be overwritten")
        self._artifacts[artifact.artifact_id] = artifact

    def link(self, source: str, target: str, relation: str) -> None:
        if source not in self._artifacts or target not in self._artifacts:
            raise ValueError("evidence link endpoints must exist")
        if not relation:
            raise ValueError("evidence link relation is required")
        self._links.add((source, target, relation))

    def has_kind(self, kind: str) -> bool:
        return any(a.kind == kind for a in self._artifacts.values())

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "artifacts": [asdict(self._artifacts[k]) for k in sorted(self._artifacts)],
            "links": [list(x) for x in sorted(self._links)],
        }

    def digest(self) -> str:
        blob = json.dumps(self.canonical_payload(), sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(blob).hexdigest()


class ReleasePolicy:
    """Fail-closed release requirements; PASS is impossible without all evidence."""
    REQUIRED_KINDS = frozenset({
        "requirements-baseline", "native-cad", "bom", "material-certificate",
        "cae-result", "test-raw-data", "test-report", "regulatory-evidence",
        "manufacturing-plan", "independent-audit",
    })

    def evaluate(self, graph: EvidenceGraph, reviews: Iterable["ReviewFinding"]) -> tuple[Decision, str]:
        missing = sorted(kind for kind in self.REQUIRED_KINDS if not graph.has_kind(kind))
        critical = [r for r in reviews if r.severity == "CRITICAL" and not r.resolved]
        if critical:
            return Decision.BLOCK, "unresolved critical independent-review findings"
        if missing:
            return Decision.BLOCK, "missing required evidence: " + ", ".join(missing)
        return Decision.PASS, "all release evidence and independent audit are present"


@dataclass(frozen=True)
class ReviewFinding:
    finding_id: str
    reviewer: str
    evidence_ids: tuple[str, ...]
    severity: str
    statement: str
    resolved: bool = False

    def __post_init__(self) -> None:
        if not self.finding_id or not self.reviewer or not self.statement or not self.evidence_ids:
            raise ValueError("review findings require identity, statement, and evidence citations")
        if self.severity not in {"INFO", "MINOR", "MAJOR", "CRITICAL"}:
            raise ValueError("invalid finding severity")


class AEOK:
    """Deterministic orchestration boundary for engineering-worker outputs."""
    _ORDER = tuple(Phase)

    def __init__(self, graph: EvidenceGraph | None = None, policy: ReleasePolicy | None = None) -> None:
        self.graph = graph or EvidenceGraph()
        self.policy = policy or ReleasePolicy()
        self.phase = Phase.REQUIREMENTS
        self.events: list[Event] = []
        self.findings: list[ReviewFinding] = []

    def record_artifact(self, artifact: Artifact) -> None:
        self.graph.add_artifact(artifact)

    def add_finding(self, finding: ReviewFinding) -> None:
        if any(item.finding_id == finding.finding_id for item in self.findings):
            raise ValueError("review finding ID is immutable")
        if finding.reviewer == "originating-worker":
            raise ValueError("originating worker cannot self-approve a finding")
        self.findings.append(finding)

    def transition(self, target: Phase, decision: Decision, reason: str, artifact_ids: Iterable[str] = ()) -> Event:
        if not reason:
            raise ValueError("every transition requires a reason")
        if self._ORDER.index(target) < self._ORDER.index(self.phase):
            raise ValueError("state machine cannot move backwards")
        ids = tuple(sorted(set(artifact_ids)))
        if any(item not in self.graph._artifacts for item in ids):
            raise ValueError("transition references unknown artifact")
        if target is Phase.RELEASE:
            decision, reason = self.policy.evaluate(self.graph, self.findings)
        previous = self.events[-1].event_hash if self.events else "0" * 64
        body = {"sequence": len(self.events) + 1, "phase": target.value, "decision": decision.value,
                "reason": reason, "artifact_ids": ids, "previous_hash": previous, "graph_hash": self.graph.digest()}
        event_hash = hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        event = Event(sequence=len(self.events) + 1, phase=target, decision=decision, reason=reason, artifact_ids=ids, previous_hash=previous, event_hash=event_hash, graph_hash=self.graph.digest())
        self.events.append(event)
        self.phase = target
        return event
