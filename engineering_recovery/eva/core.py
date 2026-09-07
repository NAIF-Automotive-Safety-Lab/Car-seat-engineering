"""Engineering Virtual & Evidence Acquisition (EVA) core.

EVA can discover, derive, bound, simulate, and prioritize work, but never
promotes calculated or simulated values to physical measurements.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from enum import IntEnum
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

from engineering_recovery.evidence.contracts import EvidenceClass, GapStatus


class Authority(IntEnum):
    MEASURED_PROJECT = 1
    CONTROLLED_SUPPLIER = 2
    PEER_REVIEWED = 3
    ANALYTICAL_DERIVATION = 4
    ENGINEERING_ESTIMATE = 5
    AI_CANDIDATE = 6


@dataclass(frozen=True)
class SourceRecord:
    source_id: str
    source: str
    source_type: str
    source_sha256: str | None
    authority_level: int
    retrieval_timestamp: str
    context: str
    status: str = "DISCOVERED"

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ParameterRecord:
    parameter_id: str
    value: float | str | None
    unit: str | None
    source: str | None
    source_type: str | None
    source_sha256: str | None
    evidence_class: EvidenceClass
    derivation_method: str | None
    equation_or_model: str | None
    assumptions: tuple[str, ...]
    uncertainty: str | None
    range_min: float | None
    range_nominal: float | None
    range_max: float | None
    distribution: str | None
    sensitivity_status: str
    validation_status: str
    timestamp: str
    tool: str
    tool_version: str
    configuration_id: str

    def __post_init__(self) -> None:
        if self.evidence_class in {EvidenceClass.DERIVED_FROM_CAD, EvidenceClass.DERIVED_FROM_DOCUMENT,
                                   EvidenceClass.CALCULATED, EvidenceClass.SIMULATED} and not self.source_sha256:
            raise ValueError("derived/calculated/simulated parameter requires source SHA-256")
        if self.range_min is not None and self.range_max is not None and self.range_min > self.range_max:
            raise ValueError("range_min cannot exceed range_max")
        if self.evidence_class is EvidenceClass.MEASURED_PHYSICALLY and self.validation_status not in {"VERIFIED", "VALIDATED"}:
            raise ValueError("physical measurements require verified/validated status")

    def as_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["evidence_class"] = self.evidence_class.value
        value["assumptions"] = list(self.assumptions)
        return value


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_file(path: str | Path, source_type: str = "PROJECT_FILE", context: str = "") -> SourceRecord:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(p)
    return SourceRecord(p.name, str(p), source_type, sha256_file(p), 1 if source_type == "MEASURED_PROJECT" else 4,
                        datetime.now(timezone.utc).isoformat(), context)


def bounded_parameter(parameter_id: str, unit: str, minimum: float, nominal: float | None, maximum: float,
                      source: str, source_sha256: str | None, assumptions: Iterable[str],
                      evidence_class: EvidenceClass = EvidenceClass.ASSUMED) -> ParameterRecord:
    return ParameterRecord(parameter_id, nominal, unit, source, "BOUNDED_ASSUMPTION", source_sha256, evidence_class,
                           "interval representation", "min <= parameter <= max", tuple(assumptions),
                           "bounded interval", minimum, nominal, maximum, "uniform_or_unspecified",
                           "NOT_RUN", "UNVERIFIED", datetime.now(timezone.utc).isoformat(), "EVA", "0.1.0",
                           f"{parameter_id}-bounded")


def derive_mass(volume: ParameterRecord, density: ParameterRecord) -> ParameterRecord:
    if volume.value is None or density.value is None:
        raise ValueError("volume and density values are required")
    if volume.unit not in {"m3", "mm3"}:
        raise ValueError("volume must be in m3 or mm3")
    conversion = 1e-9 if volume.unit == "mm3" else 1.0
    value = float(volume.value) * conversion * float(density.value)
    source_sha = volume.source_sha256 or density.source_sha256
    if not source_sha:
        raise ValueError("derived mass requires a source SHA-256")
    return ParameterRecord("mass_from_volume_density", value, "kg", ";".join(filter(None, [volume.source, density.source])),
                           "DERIVATION_INPUTS", source_sha, EvidenceClass.CALCULATED, "mass = volume * density",
                           "m = V*rho", ("uniform density",), "input uncertainty not propagated", None, value, None,
                           "INPUT_DEPENDENT", "NOT_RUN", "CALCULATED", datetime.now(timezone.utc).isoformat(),
                           "EVA", "0.1.0", "mass-derivation-v1")


def interval_propagate(function: Callable[..., float], inputs: Mapping[str, tuple[float, float]],
                       source_sha256: str, parameter_id: str, unit: str) -> ParameterRecord:
    names = list(inputs)
    corners = [dict(zip(names, values)) for values in itertools.product(*[(lo, hi) for lo, hi in inputs.values()])]
    values = [float(function(**corner)) for corner in corners]
    return ParameterRecord(parameter_id, float(function(**{k: (v[0] + v[1]) / 2 for k, v in inputs.items()})), unit,
                           "interval inputs", "UNCERTAINTY_MODEL", source_sha256, EvidenceClass.CALCULATED,
                           "interval corner propagation", "f(x) over Cartesian interval corners", ("corner bound",),
                           "corner-derived bounds", min(values), float(function(**{k: (v[0] + v[1]) / 2 for k, v in inputs.items()})), max(values),
                           "interval", "NOT_RUN", "CALCULATED", datetime.now(timezone.utc).isoformat(), "EVA", "0.1.0", f"{parameter_id}-interval")


def configuration_id(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def parameter_sweep(experiment_id: str, parameter: str, values: Iterable[float], evaluator: Callable[[float], Any],
                    input_provenance: Iterable[str], solver: str = "EVA_DETERMINISTIC") -> dict[str, Any]:
    rows = [{"parameter": parameter, "value": v, "result": evaluator(v)} for v in values]
    config = {"experiment_id": experiment_id, "parameter": parameter, "values": list(values)}
    return {"experiment_id": experiment_id, "configuration_id": configuration_id(config), "inputs": rows,
            "input_provenance": list(input_provenance), "solver": solver, "solver_version": "0.1.0",
            "results": rows, "failure_state": None, "runtime_status": "VERIFIED", "evidence_class": "SIMULATED"}


def local_sensitivity(evaluator: Callable[[float], float], nominal: float, delta: float) -> dict[str, Any]:
    if delta <= 0:
        raise ValueError("delta must be positive")
    low, high = evaluator(nominal - delta), evaluator(nominal + delta)
    return {"nominal": nominal, "delta": delta, "low_result": low, "high_result": high,
            "local_slope": (high - low) / (2 * delta), "status": "CALCULATED"}


def prioritize_test(test_id: str, affected_gaps: Iterable[str], uncertainty_count: int,
                    sensitivity: float, safety_importance: float, effort: float, unlocks: Iterable[str]) -> dict[str, Any]:
    if effort <= 0:
        raise ValueError("effort must be positive")
    score = (len(set(affected_gaps)) + uncertainty_count + sensitivity + safety_importance + len(set(unlocks))) / effort
    return {"test_id": test_id, "affected_gaps": list(affected_gaps), "uncertainties_reduced": uncertainty_count,
            "sensitivity_importance": sensitivity, "safety_importance": safety_importance,
            "effort_estimate": effort, "dependency_unlocks": list(unlocks), "information_gain_score": score,
            "status": "CANDIDATE", "physical_test_required": True}


def calibration_record(test_id: str, raw_data_sha256: str | None, model_revision: str) -> dict[str, Any]:
    if not raw_data_sha256 or len(raw_data_sha256) != 64:
        return {"test_id": test_id, "status": "TEST_REQUIRED", "raw_test_data_sha256": None,
                "model_revision": model_revision, "reason": "raw measured data unavailable"}
    return {"test_id": test_id, "status": "CALIBRATION_READY", "raw_test_data_sha256": raw_data_sha256,
            "model_revision": model_revision, "pre_test_model_preserved": True}


def model_maturity(model_id: str, level: str, evidence_refs: Iterable[str]) -> dict[str, Any]:
    allowed = {"M0": "CONCEPTUAL", "M1": "ANALYTICAL", "M2": "COMPUTATIONAL", "M3": "COMPONENT_CORRELATED", "M4": "SUBSYSTEM_CORRELATED", "M5": "SYSTEM_VALIDATED"}
    if level not in allowed:
        raise ValueError("invalid maturity level")
    refs = list(evidence_refs)
    if level in {"M3", "M4", "M5"} and not refs:
        raise ValueError("correlated/validated maturity requires evidence")
    return {"model_id": model_id, "maturity": level, "label": allowed[level], "evidence_refs": refs,
            "status": "VERIFIED" if refs else "UNVERIFIED"}
