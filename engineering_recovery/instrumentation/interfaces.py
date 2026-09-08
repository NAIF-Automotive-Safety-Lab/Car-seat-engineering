"""Physical evidence firewall and instrumentation metadata contracts."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import hashlib

from engineering_recovery.evidence.contracts import EvidenceClass


@dataclass(frozen=True)
class MeasurementRecord:
    sensor_id: str
    calibration_ref: str
    sampling_hz: float
    units: str
    timestamp: str
    test_id: str
    raw_data_sha256: str
    evidence_class: EvidenceClass = EvidenceClass.MEASURED_PHYSICALLY

    def __post_init__(self) -> None:
        if self.evidence_class is not EvidenceClass.MEASURED_PHYSICALLY:
            raise ValueError("measurement records must remain MEASURED_PHYSICALLY")
        if not self.sensor_id or not self.calibration_ref or not self.test_id:
            raise ValueError("sensor, calibration, and test identity are required")
        if self.sampling_hz <= 0 or len(self.raw_data_sha256) != 64:
            raise ValueError("sampling and raw-data SHA-256 are required")

    def as_dict(self) -> dict[str, object]:
        return {"sensor_id": self.sensor_id, "calibration": self.calibration_ref,
                "sampling_hz": self.sampling_hz, "units": self.units,
                "timestamp": self.timestamp, "test_id": self.test_id,
                "raw_data_sha256": self.raw_data_sha256,
                "evidence_class": self.evidence_class.value}


def hash_raw_data(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def blocked_physical_input(parameter: str, next_action: str) -> dict[str, object]:
    return {"parameter": parameter, "status": "BLOCKED", "evidence_class": "BLOCKED",
            "measured_value": None, "next_action": next_action,
            "reason": "No raw calibrated measurement was supplied"}
