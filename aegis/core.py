from __future__ import annotations

import hashlib
import json
import sqlite3
import uuid
import threading
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from enum import StrEnum
from pathlib import Path
from typing import Any

class GateVerdict(StrEnum):
    PASS = "PASS"
    BLOCK = "BLOCK"
    UNKNOWN = "UNKNOWN"

class Lifecycle(StrEnum):
    DRAFT = "DRAFT"
    ARCHITECTED = "ARCHITECTED"
    JOB_CREATED = "JOB_CREATED"
    EVIDENCE_REGISTERED = "EVIDENCE_REGISTERED"
    REVIEWED = "REVIEWED"
    RELEASE_BLOCKED = "RELEASE_BLOCKED"
    RELEASED = "RELEASED"

TRANSITIONS = {
    Lifecycle.DRAFT: {Lifecycle.ARCHITECTED},
    Lifecycle.ARCHITECTED: {Lifecycle.JOB_CREATED},
    Lifecycle.JOB_CREATED: {Lifecycle.EVIDENCE_REGISTERED},
    Lifecycle.EVIDENCE_REGISTERED: {Lifecycle.REVIEWED, Lifecycle.RELEASE_BLOCKED},
    Lifecycle.REVIEWED: {Lifecycle.RELEASED, Lifecycle.RELEASE_BLOCKED},
    Lifecycle.RELEASE_BLOCKED: set(), Lifecycle.RELEASED: set(),
}

ENTITY_TYPES = frozenset({"Project", "Requirement", "Claim", "ArchitectureDecision", "CADArtifact", "STEPArtifact", "Drawing", "BOM", "Material", "LoadCase", "FEModel", "Mesh", "SolverRun", "CAEResult", "OptimizationRun", "Prototype", "TestProtocol", "TestRun", "TestMeasurement", "DVRecord", "PVRecord", "RegulatoryRequirement", "RegulatoryEvidence", "SupplierEvidence", "OEMEvidence", "Audit", "AdversarialReview", "Contradiction", "ChangeSet", "ReleaseCandidate", "ReleaseGate", "WorkerJob", "FederationEndpoint", "FederationHandshake", "HashRecord", "SignatureRecord", "ProvenanceRecord"})

@dataclass(frozen=True)
class Provenance:
    created_at: str
    created_by: str
    input_refs: list[str]
    output_refs: list[str]


def now() -> str: return datetime.now(timezone.utc).isoformat()
def canonical_json(value: Any) -> bytes: return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
def sha256(value: bytes) -> str: return hashlib.sha256(value).hexdigest()

class AegisStore:
    """SQLite ledger; critical writes are append-only and transactionally hashed."""
    def __init__(self, path: str | Path = "aegis.db"):
        self.db = sqlite3.connect(str(path), check_same_thread=False)
        self.lock = threading.RLock()
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA foreign_keys = ON")
        self.db.execute("CREATE TABLE IF NOT EXISTS entities (id TEXT PRIMARY KEY, type TEXT NOT NULL, version INTEGER NOT NULL, parent_id TEXT, status TEXT NOT NULL, content_json TEXT NOT NULL, content_hash TEXT NOT NULL, created_at TEXT NOT NULL, created_by TEXT NOT NULL, provenance_json TEXT NOT NULL)")
        self.db.execute("CREATE TABLE IF NOT EXISTS events (id TEXT PRIMARY KEY, entity_id TEXT NOT NULL, action TEXT NOT NULL, actor TEXT NOT NULL, timestamp TEXT NOT NULL, input_hash TEXT NOT NULL, output_hash TEXT NOT NULL, detail_json TEXT NOT NULL)")
        self.db.execute("CREATE TABLE IF NOT EXISTS nonces (nonce TEXT PRIMARY KEY, seen_at TEXT NOT NULL)")
        self.db.commit()

    def create(self, entity_type: str, content: dict[str, Any], actor: str, *, status: str = "UNKNOWN", parent_id: str | None = None, input_refs: list[str] | None = None) -> dict[str, Any]:
        if entity_type not in ENTITY_TYPES: raise ValueError("unsupported entity type")
        entity_id = str(uuid.uuid4()); timestamp = now(); provenance = Provenance(timestamp, actor, input_refs or [], [])
        body = {**content, "id": entity_id, "type": entity_type, "version": 1, "parent_id": parent_id, "status": status, "created_at": timestamp, "created_by": actor, "provenance": asdict(provenance)}
        digest = sha256(canonical_json(body))
        with self.lock, self.db:
            self.db.execute("INSERT INTO entities VALUES (?, ?, 1, ?, ?, ?, ?, ?, ?, ?)", (entity_id, entity_type, parent_id, status, json.dumps(body), digest, timestamp, actor, json.dumps(asdict(provenance))))
            self._event(entity_id, "CREATE", actor, "", digest, {"type": entity_type})
        return {**body, "content_hash": digest}

    def get(self, entity_id: str) -> dict[str, Any] | None:
        row = self.db.execute("SELECT content_json, content_hash FROM entities WHERE id = ?", (entity_id,)).fetchone()
        return ({**json.loads(row["content_json"]), "content_hash": row["content_hash"]} if row else None)

    def list(self, entity_type: str | None = None) -> list[dict[str, Any]]:
        sql = "SELECT content_json, content_hash FROM entities" + (" WHERE type = ?" if entity_type else "") + " ORDER BY created_at"
        return [{**json.loads(r["content_json"]), "content_hash": r["content_hash"]} for r in self.db.execute(sql, (() if not entity_type else (entity_type,)))]

    def _event(self, entity_id: str, action: str, actor: str, input_hash: str, output_hash: str, detail: dict[str, Any]) -> None:
        self.db.execute("INSERT INTO events VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (str(uuid.uuid4()), entity_id, action, actor, now(), input_hash, output_hash, json.dumps(detail, sort_keys=True)))

    def consume_nonce(self, nonce: str) -> bool:
        try:
            with self.lock, self.db: self.db.execute("INSERT INTO nonces VALUES (?, ?)", (nonce, now()))
            return True
        except sqlite3.IntegrityError: return False

class AEOK:
    def __init__(self, store: AegisStore): self.store = store
    def transition(self, entity_id: str, target: Lifecycle, actor: str) -> dict[str, Any]:
        entity = self.store.get(entity_id)
        if not entity: raise KeyError(entity_id)
        current = Lifecycle(entity["status"])
        if target not in TRANSITIONS[current]: raise ValueError(f"AEOK denied transition {current}->{target}")
        body = {k:v for k,v in entity.items() if k != "content_hash"}; body["status"] = target.value; body["version"] += 1; body["parent_id"] = entity_id; body["created_at"] = now(); body["created_by"] = actor
        return self.store.create(entity["type"], {k:v for k,v in body.items() if k not in {"id", "type", "version", "parent_id", "status", "created_at", "created_by", "provenance"}}, actor, status=target.value, parent_id=entity_id, input_refs=[entity["content_hash"]])
    def release_gate(self, candidate_id: str, checks: dict[str, str], actor: str) -> dict[str, Any]:
        mandatory = {"build", "tests", "security", "aeok", "evidence_registry", "database", "federation", "routes", "frontend_backend_boundary", "mandatory_evidence", "contradictions"}
        zeroes = {"critical_errors", "broken_paths"}
        passed = mandatory <= checks.keys() and zeroes <= checks.keys() and all(checks[k] == "PASS" for k in mandatory) and all(checks[k] in {0, "0"} for k in zeroes)
        verdict = GateVerdict.PASS if passed else GateVerdict.BLOCK
        return self.store.create("ReleaseGate", {"candidate_id": candidate_id, "checks": checks, "decision": verdict.value, "fail_closed": True}, actor, status=verdict.value, input_refs=[candidate_id])
