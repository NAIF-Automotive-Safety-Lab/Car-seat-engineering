#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "artifacts/master_v3/comprehensive_gap_register_v3.json"
out_json = ROOT / "artifacts/master_v3/physical_release_requirements.json"
out_md = ROOT / "docs/master_v3/PHYSICAL_RELEASE_REQUIREMENTS.md"
data = json.loads(source.read_text())
gaps = data["gaps"]

physical_categories = {
    "E": "joints/hinges and compliance",
    "F": "materials and characterization",
    "G": "fasteners, preload, and joint stiffness",
    "H": "Lock-170 dynamics and physical response",
    "I": "Rebound-180 dynamics and stop behavior",
    "J": "absorber force/displacement/rate behavior",
    "K": "instrumentation and calibration",
    "L": "safety, containment, retention, and failure modes",
    "O": "physical input acquisition",
}

def prefix(gap_id: str) -> str:
    return gap_id[0]

counts = Counter(prefix(g["gap_id"]) for g in gaps)
selected = [g for g in gaps if prefix(g["gap_id"]) in physical_categories]
by_category: dict[str, list[dict]] = defaultdict(list)
for g in selected:
    by_category[prefix(g["gap_id"])].append(g)

requirements = []
for letter, label in physical_categories.items():
    rows = by_category.get(letter, [])
    requirements.append({
        "category": letter,
        "label": label,
        "gap_count": len(rows),
        "gap_ids": [g["gap_id"] for g in rows],
        "current_statuses": sorted({g["current_status"] for g in rows}),
        "physical_test_required_count": sum(bool(g["physical_test_required"]) for g in rows),
        "minimum_test_definitions": sorted({g["minimum_test_definition"] for g in rows if g["minimum_test_definition"]}),
        "closure_paths": sorted({g["closure_path"] for g in rows}),
        "blocking_dependencies": sorted({g["blocking_dependency"] for g in rows}),
        "required_evidence": sorted({e for g in rows for e in g["required_evidence"]}),
        "next_actions": sorted({g["next_action"] for g in rows}),
    })

summary = {
    "source": str(source),
    "source_status": data["status"],
    "total_gap_count": len(gaps),
    "category_counts": dict(sorted(counts.items())),
    "physical_release_gap_count": len(selected),
    "physical_measurements_current": 0,
    "release_status": "BLOCKED",
    "release_blockers": [
        "released V7 native CAD and authoritative interfaces",
        "released BOM and material/fastener evidence",
        "calibrated physical measurements with raw-data hashes",
        "V7-specific CAE inputs and acceptance thresholds",
    ],
    "requirements": requirements,
    "minimum_release_evidence": [
        "V7 native CAD/STEP with SHA-256 and deterministic STEP/B-Rep audit",
        "released BOM, materials, fastener, joining, PMI/GD&T, dimensions, and tolerances",
        "joint/hinge/lock/rebound/absorber characterization where sensitivity requires it",
        "calibrated instrumentation plan and imported raw measurements",
        "V7 CAE model inputs, solver configuration, output artifact hashes, and limitations",
        "evidence graph links from source to parameter/CAD entity/feature/part/interface/joint/material/physics/test/measurement/result",
    ],
}
out_json.write_text(json.dumps(summary, indent=2) + "\n")

lines = [
    "# Physical Release Requirements Extracted from A–O Gap Register v3",
    "",
    f"**Release decision: {summary['release_status']}**. The register contains **{len(gaps)} gaps**; **{len(selected)}** are directly tied to physical evidence, physical characterization, instrumentation, safety, joints, materials, fasteners, dynamics, or physical inputs. Current verified physical measurements: **0**.",
    "",
    "## Category summary",
    "",
    "| Category | Requirement domain | Gaps | Physical-test-required gaps | Current state |",
    "|---|---|---:|---:|---|",
]
for item in requirements:
    lines.append(f"| {item['category']} | {item['label']} | {item['gap_count']} | {item['physical_test_required_count']} | {', '.join(item['current_statuses'])} |")
lines += ["", "## Minimum evidence required before release", ""]
for item in summary["minimum_release_evidence"]:
    lines.append(f"- {item}")
lines += ["", "## Category-specific requirements", ""]
for item in requirements:
    lines += [f"### {item['category']} — {item['label']}", "", f"Gap IDs: `{', '.join(item['gap_ids'])}`.", ""]
    lines.append(f"Physical-test-required records: **{item['physical_test_required_count']}**.")
    if item["minimum_test_definitions"]:
        lines.append("Minimum test definition: " + "; ".join(item["minimum_test_definitions"]) + ".")
    lines.append("Blocking dependencies: " + "; ".join(item["blocking_dependencies"]) + ".")
    lines.append("Next action: " + "; ".join(item["next_actions"]) + ".")
    lines.append("")
lines += ["## Interpretation", "", "The register distinguishes computational enablement from physical closure. EVA may derive, bound, simulate, rank sensitivity, and prioritize tests, but none of those actions creates a `MEASURED_PHYSICALLY` or `VALIDATED` record. Synthetic experiments remain software tests only. Manufacturer release remains blocked until the listed authoritative and physical evidence is imported, hashed, and replayed.", ""]
out_md.write_text("\n".join(lines))
print(json.dumps({"total_gaps": len(gaps), "physical_release_gaps": len(selected), "output_json": str(out_json), "output_md": str(out_md)}, indent=2))
