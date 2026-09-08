# EERE Architecture

EERE is a repository-native, evidence-first layer around the existing Car Seat Engineering assets. It does not replace the existing SR11, R4.1, or validation paths.

| Layer | Current implementation | Status |
|---|---|---|
| Acquisition/provenance | `engineering_recovery/acquisition/service.py` | VERIFIED |
| Immutable artifact copy | SHA-checked copy with read-only evidence copy | VERIFIED |
| STEP forensic scan | `engineering_recovery/step/forensic.py` | VERIFIED, deterministic single scanner |
| OCCT/OCP B-Rep | `engineering_recovery/brep/validate.py` | VERIFIED |
| Deterministic surface inventory | `engineering_recovery/features/deterministic.py` | VERIFIED; semantic features remain UNDEFINED |
| PMI/GD&T | `engineering_recovery/pmi/extract.py` | VERIFIED absence result; PMI status NOT_PRESENT_OR_NOT_RECOVERED |
| STEPcode / step-p21 comparison | No executable bindings detected | NOT_AVAILABLE |
| Chrono dynamics | Project Chrono 10.0.0 compiled core binding; system/body/joint/solver/reaction smoke | VERIFIED |
| Gmsh | Pinned executable 4.15.2-git-657c8e9 | VERIFIED executable availability |
| CalculiX | No executable detected | NOT_AVAILABLE |
| Chrono Cascade/SR11 adapter | `pychrono.cascade` not built | BLOCKED |
| Traceability graph | Existing contract layer with provenance checks; full CAD-to-test graph pending | PARTIAL |

The executable entry point is `tools/eere_run.py`. It acquires the exact R4.1 STEP payload, verifies immutable identity, performs deterministic ISO-10303-21 scanning, imports the B-Rep through OCP, emits surface and explicit PMI status evidence, and records external runtime status. It returns a smoke `PASS` only for executed layers; aggregate runtime remains `PARTIAL` when required components are unavailable.

The current EERE boundary is intentionally narrow. Future layers must consume the manifest and source SHA rather than bypassing acquisition or substituting historical claims. Feature surface counts are deterministic evidence, while holes, slots, fillets, chamfers, interfaces, pivot axes, hinge axes, materials, fasteners, joints, and PMI remain undefined unless independently recovered and verified.
