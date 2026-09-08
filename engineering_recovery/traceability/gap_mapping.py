from __future__ import annotations

from typing import Any


def gap_ids() -> list[str]:
    groups = {"A": 12, "B": 20, "C": 17, "D": 10, "E": 11, "F": 9, "G": 11, "H": 12, "I": 11, "J": 11, "K": 12, "L": 10, "M": 11, "N": 10, "O": 14}
    return [f"{letter}{index:02d}" for letter, total in groups.items() for index in range(1, total + 1)]


def build_gap_mapping(existing_register: dict[str, Any] | None = None) -> dict[str, Any]:
    rows = []
    for gap_id in gap_ids():
        existing = (existing_register or {}).get(gap_id, {})
        rows.append({
            "GAP_ID": gap_id,
            "CURRENT_STATUS": existing.get("CURRENT_STATUS", "UNVERIFIED"),
            "BLOCKER_TYPE": existing.get("BLOCKER_TYPE", "EVIDENCE_REQUIRED"),
            "AUTOMATABLE": existing.get("AUTOMATABLE", "PARTIAL"),
            "REQUIRED_EVIDENCE": existing.get("REQUIRED_EVIDENCE", "PROJECT_REGISTER_OR_AUTHORITATIVE_ARTIFACT"),
            "EERE_COMPONENT": existing.get("EERE_COMPONENT", "engineering_recovery"),
            "VERIFICATION_METHOD": existing.get("VERIFICATION_METHOD", "EVIDENCE_MANIFEST_AND_DETERMINISTIC_TEST"),
            "REMAINING_HUMAN_INPUT": existing.get("REMAINING_HUMAN_INPUT", "OWNER_AND_ENGINEERING_EVIDENCE"),
        })
    return {"system": "EERE", "rule": "Software installation never closes a gap", "gaps": rows}
