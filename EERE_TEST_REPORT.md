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

Generated evidence is under `artifacts/engineering-evidence/`, including `artifact_manifest.json`, `step_forensic_report.json`, `brep_validation.json`, `runtime_status.json`, and `eere_smoke.json`.

## BLOCKED

Chrono runtime smoke is blocked because `pychrono.core` cannot be imported. No dynamic system, body, joint, or solver result is claimed.

## NOT_AVAILABLE

Dual-parser comparison through STEPcode and step-p21, PMI extraction, deterministic feature extraction, FE, and CAM are not yet executable.

## TEST SCOPE NOTE

An unscoped repository-level pytest run collects unrelated vendored projects under `EXTERNAL_PROJECTS`. The authoritative owned-suite command is the command shown above, which explicitly runs EERE tests and the existing repository tests.
