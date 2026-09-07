# EERE Architecture

EERE is a repository-native, evidence-first layer around the existing Car Seat Engineering assets. It does not replace the existing SR11, R4.1, or validation paths.

| Layer | Current implementation | Status |
|---|---|---|
| Acquisition/provenance | `engineering_recovery/acquisition/service.py` | VERIFIED |
| Immutable artifact copy | SHA-checked copy with read-only evidence copy | VERIFIED |
| STEP forensic scan | `engineering_recovery/step/forensic.py` | VERIFIED, single deterministic scanner |
| OCCT/OCP B-Rep | `engineering_recovery/brep/validate.py` | VERIFIED |
| Feature extraction | Not implemented | NOT_AVAILABLE |
| PMI/GD&T | Not implemented; no absence inference | NOT_AVAILABLE |
| STEPcode / step-p21 comparison | No executable bindings detected | NOT_AVAILABLE |
| Chrono dynamics | Existing adapter retained; `pychrono.core` unavailable | BLOCKED |
| FE/CAM | No Gmsh or CalculiX executable detected | BLOCKED |
| Traceability graph | Existing project manifests remain authoritative; EERE graph layer pending | NOT_AVAILABLE |

The executable entry point is `tools/eere_run.py`. It acquires the exact R4.1 STEP payload, emits a manifest, performs deterministic ISO-10303-21 scanning, imports the B-Rep through OCP, emits geometry evidence, and records external runtime status. It returns PASS only for the layers it actually executes; blocked runtimes are not converted into success.

The current EERE boundary is intentionally narrow. Future layers must consume the manifest and source SHA rather than bypassing acquisition or substituting historical claims.
