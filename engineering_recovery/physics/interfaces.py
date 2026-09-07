"""Physics interfaces: contracts only; no unverified parameter defaults."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from engineering_recovery.evidence.contracts import EvidenceClass, GapStatus, make_blocked_record


@dataclass(frozen=True)
class PhysicsParameter:
    name: str
    value: float | str | None
    unit: str | None
    evidence_class: EvidenceClass
    source: str | None
    source_sha256: str | None
    status: GapStatus

    def as_dict(self) -> dict[str, Any]:
        return {"name": self.name, "value": self.value, "unit": self.unit,
                "evidence_class": self.evidence_class.value, "source": self.source,
                "source_sha256": self.source_sha256, "status": self.status.value}


def lock_170_record(candidate: str) -> dict[str, Any]:
    if candidate not in {"L170-A", "L170-B", "L170-C"}:
        raise ValueError("unsupported Lock-170 candidate")
    return {"module": "lock_170", "candidate": candidate, "final_selection": "NOT_AUTHORIZED",
            "status": GapStatus.BLOCKED.value,
            "parameters": [make_blocked_record(p, "Provide authoritative model/test evidence", "No released lock mechanism evidence") for p in
                           ["acceleration_threshold", "response_time", "joint_dynamics", "inertia", "contact", "actuator", "loads", "energy", "structural_response"]]}


def rebound_180_record() -> dict[str, Any]:
    return {"module": "rebound_180", "status": GapStatus.BLOCKED.value,
            "target_180_mm": PhysicsParameter("rebound_target", 180, "mm", EvidenceClass.TARGET, "existing design record", None, GapStatus.CANDIDATE).as_dict(),
            "parameters": [make_blocked_record(p, "Provide measured or authoritative rebound evidence", "No calibrated rebound test or validated model") for p in
                           ["displacement", "velocity", "acceleration", "rebound", "damping", "contact_force", "stop_condition", "energy_dissipation"]]}


def absorber_record() -> dict[str, Any]:
    return {"module": "absorbers", "status": "UNVERIFIED", "absorber_status": "UNVERIFIED",
            "curves": [make_blocked_record(p, "Import hashed measured curve", "No measured absorber characterization") for p in
                       ["force_displacement", "force_velocity", "nonlinear_response", "hysteresis", "damping", "energy_absorption", "temperature_effects", "rate_effects"]]}
