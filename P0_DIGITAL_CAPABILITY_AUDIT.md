# P0 Digital Capability Audit

**Date:** 2026-09-06  
**Scope:** CAD/STEP/B-Rep, P0 digital definition, DFAM, FE/MBD readiness, evidence integrity  
**R4.1:** immutable; SHA256 `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

## Relevant skills and actual status

| Skill | Version/runtime | Smoke result | Classification | P0/V4 usefulness |
|---|---|---|---|---|
| `cad` | OCP 7.9.3.1 | `tests/test_step_brep.py`: PASS; STEP read and 62 solids | VERIFIED | R4.1 readback, B-Rep, extraction |
| `cad-viewer` | `cadgen` executable absent | Cannot launch viewer | UNAVAILABLE | Would provide visual CAD review only |
| `dfam-check` | Script present; requires mesh input | No STL/OBJ/3MF input and CAD viewer/export runtime absent | EXTERNAL_REQUIRED | Manufacturing/printability only, not structural proof |
| `ansys-structural-workbench` | No ANSYS MCP/server or executable | No live Workbench capability discovered | UNAVAILABLE | FE orchestration cannot run |
| `forgecad-build-model` | `forgecad` executable absent | Not executable | UNAVAILABLE | No P0 reconstruction or authoring |
| `forgecad-design-spec` | Documentation only | Planning capability; no execution runtime | UNVERIFIED | Not used; P0 architecture frozen |
| `forgecad-grade-model` | `forgecad` executable absent | Not executable | UNAVAILABLE | No grade generated |
| `forgecad-inspect-model` | `forgecad` executable absent | Not executable | UNAVAILABLE | No ForgeCAD inspection bundle |
| `forgecad-reconstruct-cad-file` | `forgecad` executable absent | Not executable | UNAVAILABLE | Reconstruction forbidden for R4.1 anyway |
| `step-parts` | Hosted catalog workflow | Not needed; no part procurement requested | EXTERNAL_REQUIRED | Cannot substitute P0 geometry |
| `pinned-external-cad-runtime-restoration` | verifier present; no runtime manifest/source pin | No pinned runtime to verify | BLOCKED | No external CAD adapter restoration |
| `sourceless-authoritative-cad-execution` | Lifecycle guidance only | No source-less CAD service in this repository | UNAVAILABLE | Not applicable to immutable R4.1 |
| `legacy-cad-adapter-executor-migration` | Migration guidance only | No adapter migration requested or executed | UNAVAILABLE | Not applicable to current read-only evidence |

## Executed digital checks

| Check | Command | Result |
|---|---|---|
| STEP read/B-Rep smoke | `python3 tests/test_step_brep.py` | PASS, return code 0 |
| SHA256 | `sha256sum R4.1.step R4.1/R4.1.step` | Both exact required SHA |
| JSON parsing | Python strict parser over canonical project JSON | 61 valid, 0 invalid outside quarantine |
| Secure intake validator | `python3 EXTERNAL_DATA_INTAKE_PACKAGE/validate_intake.py` | Correct fail-closed rejection; no records CLOSED |
| PyChrono smoke | `python3 pychrono_smoke_test.py` | FAIL before artifact creation: missing `artifacts/pychrono_smoke.json`; not accepted as runtime proof |
| FE executables | `calculix`, `ccx`, `ansys`, `abaqus` | Not present |
| CAD/ForgeCAD executables | `cadgen`, `forgecad` | Not present |
| Physical devices | `/dev/video*`, `/dev/ttyUSB*`, `/dev/ttyACM*` | Not present |

## P0 digital results

- **R4.1/P0 CAD extraction:** PASS for digital STEP readback, 62 solids, mapping 62/62, and B01–B06 second-pass evidence already committed.
- **P0 build-definition completeness:** PARTIAL. The P0 package and work package define architecture, interfaces, sequence, and gates. Build-critical tolerances, travel, clearances, calibration, and physical article evidence remain unavailable where the extraction explicitly marks them unavailable.
- **DFAM:** BLOCKED/NOT APPLICABLE to raw STEP in the current environment. The DFAM tool accepts mesh input and the required CAD-to-mesh runtime is unavailable. No printability claim is made.
- **FE preparation:** BLOCKED. No solver runtime/MCP, no complete 14-input physical closure, and no authorized FE execution.
- **MBD readiness:** BLOCKED/UNVERIFIED. PyChrono smoke failed before artifact creation; no validated MBD runtime or physical-input contract is available.
- **Evidence/package integrity:** PASS for the new repaired matrix and digital evidence hashes; original corrupted bytes remain unchanged in quarantine.

## Physical boundary

Fabrication, assembly, metrology, calibration, DAQ, physical photographs/video, and physical measurements remain **EXTERNAL_EXECUTION_REQUIRED**. No synthetic or simulated values are promoted as physical evidence.

## Frozen constraints

`R4.1` was not modified. `R4.2` was not created. P0 architecture was not redesigned. No physical fabrication or inspection was claimed.
