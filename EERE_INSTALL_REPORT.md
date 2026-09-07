# EERE Installation Report

## INSTALLED

The first EERE layer is installed in the repository under `engineering_recovery/` and `tools/eere_run.py`. The environment contains `cadquery-ocp 8.0.1.0.0` and `pytest 9.1.1`. The implementation uses the existing OCP binding rather than replacing the project architecture or vendoring OCCT. Project Chrono 10.0.0 was built from source commit `9faf13dd8f1128dd75ed233a9627027b0422c3f7` with GCC, CMake, and SWIG.

## VERIFIED

`python3 tools/eere_run.py` completed successfully. It created a read-only immutable copy of `R4.1/R4.1.step`, generated a byte-identity manifest, scanned the STEP exchange structure, imported the model through OCP, and counted 62 solids. With the pinned Chrono build environment, the same runner reports `pychrono_core=VERIFIED`.

## FAILED

The initial EERE smoke execution exposed three OCP 8 binding API differences. These were corrected using the supported `_s` methods and explicit bounding-box getters. The corrected run passed and the fixes are covered by regression tests.

## BLOCKED

Gmsh and CalculiX are unavailable. The earlier similarly named PyPI package did not provide the required Chrono binding; it was not accepted. The real binding is now built and verified separately.

## NOT_AVAILABLE

STEPcode, step-p21, Analysis Situs, FreeCAD, CadQuery, feature recognition, PMI/GD&T extraction, and the full engineering traceability graph are not executable in this phase.

## NOT_APPLICABLE

No AI-assisted geometry recognition was used. No physical value was promoted from UNKNOWN, and no FE/MBD/fabrication result was generated.
