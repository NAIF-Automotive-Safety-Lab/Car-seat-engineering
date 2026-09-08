"""Typed, fail-closed evidence contracts for the Master v2 engineering layers."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from enum import Enum
from typing import Any


class EvidenceClass(str, Enum):
    DERIVED_FROM_CAD = "DERIVED_FROM_CAD"
    DERIVED_FROM_DOCUMENT = "DERIVED_FROM_DOCUMENT"
    MEASURED_PHYSICALLY = "MEASURED_PHYSICALLY"
    CALCULATED = "CALCULATED"
    SIMULATED = "SIMULATED"
    ASSUMED = "ASSUMED"
    TARGET = "TARGET"
    CANDIDATE = "CANDIDATE"
    UNVERIFIED = "UNVERIFIED"
    BLOCKED = "BLOCKED"
    VALIDATED = "VALIDATED"


class GapStatus(str, Enum):
    OPEN = "OPEN"
    IN_PROGRESS = "IN_PROGRESS"
    CANDIDATE = "CANDIDATE"
    DERIVED = "DERIVED"
    VERIFIED = "VERIFIED"
    VALIDATED = "VALIDATED"
    BLOCKED = "BLOCKED"
    NOT_PRESENT = "NOT_PRESENT"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class EvidenceRef:
    source: str | None
    source_sha256: str | None
    evidence_class: EvidenceClass
    tool: str | None
    tool_version: str | None
    verification_state: str
    timestamp: str
    confidence: str = "UNVERIFIED"

    @classmethod
    def blocked(cls, reason: str) -> "EvidenceRef":
        return cls(source=None, source_sha256=None, evidence_class=EvidenceClass.BLOCKED,
                   tool=None, tool_version=None, verification_state=reason,
                   timestamp=datetime.now(timezone.utc).isoformat(), confidence="NONE")

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["evidence_class"] = self.evidence_class.value
        return value


def reject_physical_claim(evidence: EvidenceRef) -> None:
    """Prevent software-derived records from being labeled physical measurements."""
    if evidence.evidence_class is EvidenceClass.MEASURED_PHYSICALLY:
        if not evidence.source or not evidence.source_sha256:
            raise ValueError("MEASURED_PHYSICALLY requires a hashed raw measurement source")
        if evidence.verification_state not in {"VERIFIED", "VALIDATED"}:
            raise ValueError("MEASURED_PHYSICALLY requires verified or validated evidence")


def make_blocked_record(item: str, next_action: str, reason: str) -> dict[str, Any]:
    return {"item": item, "status": GapStatus.BLOCKED.value, "evidence": EvidenceRef.blocked(reason).to_dict(), "next_action": next_action}
