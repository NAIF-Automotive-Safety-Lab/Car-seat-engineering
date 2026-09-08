# EERE Test Report

## VERIFIED

| Command | Result |
|---|---|
| `python3 tools/eere_run.py` | PASS; SHA-bound acquisition, STEP forensic scan, OCP B-Rep, 62 solids |
| `python3 -m pytest -q engineering_recovery/tests tests` | **10 passed** |
| Acquisition replay test | PASS |
| STEP forensic reproducibility test | PASS |
| B-Rep reproducibility and validity test | PASS |
| Existing repository static integrity tests | PASS |
| `tools/run_pychrono_smoke.sh` with Chrono 10.0.0 build tree | PASS; compiled `pychrono.core`, system, body, revolute joint, solver step, and reaction extraction |

Generated evidence is under `artifacts/engineering-evidence/`, including `artifact_manifest.json`, `step_forensic_report.json`, `brep_validation.json`, `runtime_status.json`, and `eere_smoke.json`.

## VERIFIED

The real PyChrono binding was built from Project Chrono tag `10.0.0`, source commit `9faf13dd8f1128dd75ed233a9627027b0422c3f7`. The smoke artifact records `import=PASS`, `compiled_core_present=true`, system/body/joint creation, a 10 ms solver advance, and reaction extraction. This is a runtime smoke result, not a car-seat dynamics validation result.

## BLOCKED

Gmsh and CalculiX remain unavailable. No FE or manufacturing solver result is claimed.

## NOT_AVAILABLE

Dual-parser comparison through STEPcode and step-p21, PMI extraction, deterministic feature extraction, FE, and CAM are not yet executable.

## TEST SCOPE NOTE

An unscoped repository-level pytest run collects unrelated vendored projects under `EXTERNAL_PROJECTS`. The authoritative owned-suite command is the command shown above, which explicitly runs EERE tests and the existing repository tests.
