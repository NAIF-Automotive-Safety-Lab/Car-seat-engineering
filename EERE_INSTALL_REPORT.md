# EERE Installation Report

## INSTALLED

EERE is installed as a repository-native layer under `engineering_recovery/` with `tools/eere_run.py`. The implementation preserves the existing Car Seat Engineering paths and uses the installed OCP/OpenCascade binding for STEP and B-Rep execution. Project Chrono 10.0.0 was built from source commit `9faf13dd8f1128dd75ed233a9627027b0422c3f7`.

## VERIFIED

The executed EERE smoke path acquired the exact `R4.1/R4.1.step` payload and created a read-only evidence copy with:

- size: `1,220,000` bytes;
- SHA-256: `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`;
- STEP entity count: `22,119`;
- OCP solids: `62`;
- faces: `402`;
- edges: `1,550`;
- vertices: `3,100`;
- valid B-Rep: `true`;
- deterministic planar/cylindrical surface inventory: `355 / 47`;
- Project Chrono core system/body/joint/solver/reaction smoke: `PASS`;
- Gmsh executable `4.15.2-git-657c8e9`: `VERIFIED`.

The smoke output now includes `feature_inventory.json` and `pmi_report.json`. PMI is reported as `NOT_PRESENT_OR_NOT_RECOVERED` because the inspected STEP payload contains zero PMI entities. No semantic holes, slots, pivots, or hinges are inferred from surface counts.

## FAILED

The first feature smoke attempt failed because the OCP binding requires an explicit `TopoDS.Face` cast for `BRepAdaptor_Surface`. The implementation was corrected and the full suite subsequently passed with `29 passed`.

## BLOCKED

The existing SR11 adapter still requires `pychrono.cascade`, which is not built. CalculiX was not available as an executable in the inspected environment. The EERE runtime aggregate therefore remains `PARTIAL`; unavailable components are not converted to success.

## NOT_AVAILABLE

STEPcode and step-p21 executable bindings are not available. Independent dual-parser comparison, PMI/GD&T semantic extraction, semantic feature labels, and Analysis Situs/FreeCAD/CadQuery integrations are not claimed.

## NOT_APPLICABLE

No AI-generated geometry was used. No physical value was promoted from target, assumed, or unverified status. No manufacturer-release artifact or compliance result was generated.
