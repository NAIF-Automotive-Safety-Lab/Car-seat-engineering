from __future__ import annotations

import argparse
import json
from pathlib import Path

from .acquisition import ArtifactAcquisitionService
from .brep import validate_brep
from .features import extract_features, pmi_report
from .runtime import chrono_smoke, dependency_status
from .step import build_step_forensic_report
from .traceability import EvidenceGraph, build_gap_mapping, normalized_joint


def dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def run(project_root: Path, step_path: Path) -> dict[str, object]:
    project_root = project_root.resolve(); step_path = step_path.resolve()
    evidence = project_root / "artifacts/engineering-evidence"
    service = ArtifactAcquisitionService(project_root, evidence)
    artifact = service.acquire_file(step_path, artifact_id="r4_1_authoritative")
    forensic = build_step_forensic_report(step_path)
    brep = validate_brep(step_path)
    features = extract_features(step_path)
    pmi = pmi_report(step_path)
    runtimes = dependency_status()
    chrono = chrono_smoke()
    graph = EvidenceGraph()
    graph.add_node("R4.1", "CAD_ARTIFACT", sha256=artifact["manifest"]["sha256"])
    graph.add_node("P0", "ENGINEERING_CONFIGURATION", status="UNVERIFIED")
    graph.add_edge("R4.1", "P0", {"source": "current_artifact_manifest", "status": "VERIFIED"})
    graph.add_node("joint:unknown", "JOINT", value=normalized_joint("joint:unknown", "UNDEFINED", "UNDEFINED"))
    result = {
        "artifact": artifact,
        "step_forensic": forensic,
        "brep_validation": brep,
        "features": features,
        "pmi": pmi,
        "runtime": runtimes,
        "chronoruntime_smoke": chrono,
        "traceability": graph.as_dict(),
        "gap_mapping": build_gap_mapping(),
    }
    out = project_root / "artifacts/engineering-evidence/r4_1_eere"
    for key, value in result.items():
        dump(out / f"{key}.json", value)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="eere")
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--step", type=Path, default=None)
    args = parser.parse_args(argv)
    root = args.project_root.resolve(); step = (args.step or root / "R4.1/R4.1.step").resolve()
    if not step.is_file():
        print("EERE_STATUS=BLOCKED"); print(f"REASON=MISSING_STEP:{step}"); return 2
    result = run(root, step)
    print("EERE_ACQUISITION=VERIFIED")
    print(f"EERE_STEP_FORENSIC={result['step_forensic']['dual_parse_status']}")
    print(f"EERE_BREP={result['brep_validation']['status']}")
    print(f"EERE_CHRONO={result['chronoruntime_smoke']['status']}")
    print("EERE_PHYSICAL_EVIDENCE=NOT_AVAILABLE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
