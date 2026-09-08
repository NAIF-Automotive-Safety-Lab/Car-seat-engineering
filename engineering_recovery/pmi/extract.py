from __future__ import annotations

from pathlib import Path
from typing import Any

from engineering_recovery.acquisition.service import sha256_file
from engineering_recovery.step.forensic import forensic_report


def pmi_report(path: Path) -> dict[str, Any]:
    forensic = forensic_report(path)
    count = int(forensic.get("pmi_entity_count", 0))
    return {
        "status": "VERIFIED" if count > 0 else "NOT_PRESENT_OR_NOT_RECOVERED",
        "source_sha256": sha256_file(path),
        "pmi_entity_count": count,
        "dimensions": [],
        "datum_features": [],
        "geometric_tolerances": [],
        "units": [],
        "authority_rule": "PMI is reported only when present in the inspected STEP entity inventory; no PMI is inferred",
    }


__all__ = ["pmi_report"]
