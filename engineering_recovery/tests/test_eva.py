from __future__ import annotations

from pathlib import Path

import pytest

from engineering_recovery.eva.core import (
    EvidenceClass,
    ParameterRecord,
    bounded_parameter,
    calibration_record,
    derive_mass,
    local_sensitivity,
    model_maturity,
    parameter_sweep,
    prioritize_test,
)
from engineering_recovery.eva.extraction.service import extract_numeric_context, normalize_unit


def _source(name: str) -> str:
    return "a" * 64


def test_bounded_parameter_preserves_assumption_class() -> None:
    p = bounded_parameter("friction", "1", 0.1, 0.2, 0.3, "test", _source("x"), ["synthetic bound"])
    assert p.evidence_class is EvidenceClass.ASSUMED
    assert p.range_min == 0.1 and p.range_max == 0.3


def test_mass_derivation_is_calculated_not_measured() -> None:
    v = bounded_parameter("volume", "m3", 2.0, 2.0, 2.0, "cad", _source("v"), [])
    d = bounded_parameter("density", "kg/m3", 1000, 1000, 1000, "datasheet", _source("d"), [])
    m = derive_mass(v, d)
    assert m.value == 2000
    assert m.evidence_class is EvidenceClass.CALCULATED


def test_derived_parameter_requires_source_hash() -> None:
    with pytest.raises(ValueError):
        ParameterRecord("x", 1, "m", "source", "x", None, EvidenceClass.CALCULATED, "x", "x", (), None, None, 1, None, "", "NOT_RUN", "CALCULATED", "now", "EVA", "1", "c")


def test_synthetic_sweep_has_configuration_and_simulated_class() -> None:
    result = parameter_sweep("SYN-1", "x", [0, 1], lambda x: x + 1, ["SYNTHETIC_TEST_DATA"])
    assert result["configuration_id"]
    assert result["evidence_class"] == "SIMULATED"


def test_sensitivity_and_priority_are_calculated() -> None:
    sens = local_sensitivity(lambda x: x * x, 2.0, 0.1)
    priority = prioritize_test("T1", ["F01"], 1, 2, 3, 1, ["material"])
    assert sens["status"] == "CALCULATED"
    assert priority["information_gain_score"] > 0


def test_calibration_stays_test_required_without_raw_data() -> None:
    assert calibration_record("T1", None, "M0")["status"] == "TEST_REQUIRED"


def test_validation_maturity_requires_evidence() -> None:
    with pytest.raises(ValueError):
        model_maturity("model", "M5", [])


def test_document_extraction_retains_context_and_hash() -> None:
    records = extract_numeric_context(Path("README.md"), page_or_section="README")
    assert records
    assert all(record["source_sha256"] and record["exact_context"] for record in records)


def test_unit_normalization_is_explicit() -> None:
    assert normalize_unit(180, "mm", "m") == pytest.approx(0.18)
    with pytest.raises(ValueError):
        normalize_unit(1, "g", "m/s2")
