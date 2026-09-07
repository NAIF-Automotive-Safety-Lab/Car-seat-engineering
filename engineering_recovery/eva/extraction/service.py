"""Conservative document extraction: every number remains bound to context and source hash."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from engineering_recovery.eva.core import sha256_file

_NUMBER = re.compile(r"(?P<value>[-+]?\d+(?:\.\d+)?)\s*(?P<unit>mm|m|kg|N|Hz|g|ms|MPa|GPa|°|deg)?", re.I)

_UNIT_FACTORS = {
    ("mm", "m"): 1e-3,
    ("m", "m"): 1.0,
    ("kg", "kg"): 1.0,
    ("N", "N"): 1.0,
    ("ms", "s"): 1e-3,
    ("Hz", "Hz"): 1.0,
}


def normalize_unit(value: float, unit: str, target_unit: str) -> float:
    key = (unit, target_unit)
    if key not in _UNIT_FACTORS:
        raise ValueError(f"unsupported unit conversion: {unit} -> {target_unit}")
    return value * _UNIT_FACTORS[key]


def extract_numeric_context(path: str | Path, *, page_or_section: str, context_window: int = 80) -> list[dict[str, Any]]:
    source = Path(path)
    text = source.read_text(encoding="utf-8", errors="replace")
    source_hash = sha256_file(source)
    records: list[dict[str, Any]] = []
    for match in _NUMBER.finditer(text):
        start = max(0, match.start() - context_window)
        end = min(len(text), match.end() + context_window)
        context = text[start:end].replace("\n", " ").strip()
        if not context:
            continue
        records.append({"value": float(match.group("value")), "unit": match.group("unit"),
                        "source": str(source), "source_sha256": source_hash,
                        "page_or_section": page_or_section, "exact_context": context,
                        "evidence_class": "DERIVED_FROM_DOCUMENT", "status": "EXTRACTED"})
    return records
