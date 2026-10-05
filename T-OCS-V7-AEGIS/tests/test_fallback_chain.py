from __future__ import annotations

import json
import sys
from pathlib import Path

AEGIS_ROOT = Path(__file__).resolve().parents[1]
FALLBACK = AEGIS_ROOT / "G2" / "ALT" / "FALLBACK_CHAIN" / "qualification_003.json"
RUNNER = AEGIS_ROOT / "G2" / "ALT" / "run_fallback_chain.py"


def test_fallback_chain_is_complete_and_fail_closed() -> None:
    """All approved alternatives must be represented without a false PASS."""
    record = json.loads(FALLBACK.read_text(encoding="utf-8"))
    assert record["automatic_fallback_order"] == [
        "PYCHRONO",
        "PYBULLET",
        "MUJOCO",
        "DRAKE",
        "DART",
        "PINOCCHIO",
    ]
    assert record["selected_python_engine"] == "MUJOCO"
    assert record["engine_qualification"] == "PASS"
    assert record["g2_alt_execution"] == "BLOCKED"
    assert record["g2_alt_repeatability"] == "PASS"
    assert record["g3"] == "LOCKED"
    assert record["r4_1"] == "UNCHANGED"
    assert len(record["engines"]) == 3
    assert [item["import_status"] for item in record["engines"]] == ["FAILED", "FAILED", "PASS"]
    assert [item["status"] for item in record["engines"]] == ["BLOCKED", "BLOCKED", "PASS"]


def test_fallback_runner_uses_the_required_order_and_no_installs() -> None:
    source = RUNNER.read_text(encoding="utf-8")
    positions = [source.index(name) for name in ("PYCHRONO", "PYBULLET", "MUJOCO", "DRAKE", "DART", "PINOCCHIO")]
    assert positions == sorted(positions)
    assert "pip install" not in source
    assert "network_access_attempted\": False" in source
    assert sys.executable
