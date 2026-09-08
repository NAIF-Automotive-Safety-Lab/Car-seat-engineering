# EERE Compatibility Report

## Existing repository

The repository is a Python-based Car Seat Engineering evidence package with an existing OCP-based STEP importer in `model/r41_model.py`, immutable R4.1 path definitions in `bridge/r41_paths.py`, unittest-based integrity tests, and a pre-existing PyChrono smoke harness. There is no root Python packaging manifest or Docker build system. The existing CI is a contract/gate workflow rather than a general dependency installer.

## Selected integration

EERE uses the already-installed `cadquery-ocp==7.9.3.1.1` binding as the executable OCCT layer. It does not replace the existing importer or mutate R4.1. Acquisition and reports are repository-native Python modules using the standard library.

## Compatibility states

| Component | State | Evidence | Decision |
|---|---|---|---|
| OCP/OCCT | VERIFIED | Import, STEP read, B-Rep import, 62 solids | Use existing runtime |
| STEPcode | NOT_AVAILABLE | Module not installed | Do not claim dual-library semantic parse |
| step-p21 | NOT_AVAILABLE | Module not installed | Raw EERE scanner remains clearly labeled, not a substitute |
| Project Chrono/PyChrono | BLOCKED | `ModuleNotFoundError` in smoke test | Preserve blocked status |
| FreeCAD | NOT_AVAILABLE | Module not installed | Optional, not integrated |
| CadQuery | NOT_AVAILABLE | Module not installed; CadQuery-OCP is distinct | Do not conflate bindings |
| Analysis Situs | NOT_INTEGRATED | No compatible local package detected | Do not clone into root |
| BrepMFR | NOT_APPLICABLE | Optional candidate-only AI component | No authoritative use |
| SLSA | NOT_INTEGRATED | No source checkout or binary added | Provenance concepts implemented locally |

No floating dependency is silently used. Missing components remain `NOT_AVAILABLE`, `BLOCKED`, or `NOT_INTEGRATED` as appropriate.
