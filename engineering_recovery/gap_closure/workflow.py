"""Deterministic gap closure helpers; no implicit closure states."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from typing import Any

from engineering_recovery.evidence.contracts import GapStatus


@dataclass(frozen=True)
class GapRecord:
    gap_id: str
    description: str
    current_status: GapStatus
    blocker_type: str
    automatable: bool
    required_evidence: tuple[str, ...]
    required_engine: tuple[str, ...]
    verification_method: str
    current_evidence: tuple[str, ...]
    next_action: str
    owner_class: str
    dependencies: tuple[str, ...]
    closure_criteria: str
    last_verified: str

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["current_status"] = self.current_status.value
        return result


def new_blocked_gap(gap_id: str, description: str, next_action: str, *, automatable: bool = False) -> GapRecord:
    return GapRecord(gap_id, description, GapStatus.BLOCKED, "MISSING_AUTHORITATIVE_EVIDENCE", automatable,
                     ("authoritative source artifact",), ("appropriate deterministic engine",),
                     "source-hash-bound deterministic verification", (), next_action, "engineering authority",
                     ("source artifact",), "source, SHA-256, deterministic result, and review record", datetime.now(timezone.utc).isoformat())


def close_gap(record: GapRecord, evidence: tuple[str, ...], status: GapStatus) -> GapRecord:
    if status not in {GapStatus.DERIVED, GapStatus.VERIFIED, GapStatus.VALIDATED}:
        raise ValueError("close_gap requires a positive evidence status")
    if not evidence:
        raise ValueError("a gap cannot close without evidence references")
    return GapRecord(record.gap_id, record.description, status, record.blocker_type, record.automatable,
                     record.required_evidence, record.required_engine, record.verification_method,
                     evidence, record.next_action, record.owner_class, record.dependencies,
                     record.closure_criteria, datetime.now(timezone.utc).isoformat())
