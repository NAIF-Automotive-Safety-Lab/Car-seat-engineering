#!/usr/bin/env python3
import argparse
import hashlib
import json
import sys
from pathlib import Path

EXPECTED_SHA256 = "fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"
EXPECTED_SOLIDS = 62


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_step_solids(path: Path) -> tuple[str, int]:
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    reader = STEPControl_Reader()
    status = reader.ReadFile(str(path))
    if "RetDone" not in str(status):
        raise RuntimeError(f"STEP read failed: {status}")

    reader.TransferRoots()
    shape = reader.OneShape()
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    count = 0
    while explorer.More():
        count += 1
        explorer.Next()
    return str(status), count


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("step_path")
    ap.add_argument("output_path")
    ap.add_argument("--readability-only", action="store_true")
    ap.add_argument("--solid-count-only", action="store_true")
    args = ap.parse_args()

    step = Path(args.step_path)
    out = Path(args.output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    result = {
        "file": str(step),
        "expected_sha256": EXPECTED_SHA256,
        "expected_solids": EXPECTED_SOLIDS,
        "checks": [],
        "geometry_mutated": False,
    }

    if not step.is_file():
        print("R4.1_PRESENT=FAIL")
        result["status"] = "FAIL"
        result["failure"] = "FILE_NOT_FOUND"
        out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return 2

    digest = sha256(step)
    result["checks"].append({
        "name": "sha256",
        "pass": digest == EXPECTED_SHA256,
        "actual": digest,
        "expected": EXPECTED_SHA256,
    })

    try:
        status, solids = read_step_solids(step)
        result["checks"].append({
            "name": "step_read",
            "pass": True,
            "status": status,
        })
        result["checks"].append({
            "name": "solid_count",
            "pass": solids == EXPECTED_SOLIDS,
            "actual": solids,
            "expected": EXPECTED_SOLIDS,
        })
    except Exception as exc:
        result["status"] = "FAIL"
        result["failure"] = f"{type(exc).__name__}: {exc}"
        out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print("STEP_READ=FAIL")
        return 3

    if args.readability_only and not result["checks"][-2]["pass"]:
        out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return 4

    if args.solid_count_only and solids != EXPECTED_SOLIDS:
        print("SOLID_COUNT=FAIL")
        out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return 5

    if digest != EXPECTED_SHA256:
        print("SHA_MATCH=FAIL")
        result["status"] = "FAIL"
        out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return 6

    if solids != EXPECTED_SOLIDS:
        print("SOLID_COUNT=FAIL")
        result["status"] = "FAIL"
        out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return 7

    print("STEP_READ=PASS")
    print("SOLID_COUNT=PASS")
    result["status"] = "PASS"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
