# EERE Architecture

EERE is a repository-native, fail-closed evidence layer for Car Seat Engineering. Acquisition preserves immutable bytes and manifests. The STEP layer performs a raw ISO-10303-21 entity scan and an independent OCP/OCCT import. The B-Rep layer computes topology counts and bounds with source SHA provenance. Feature and PMI layers never promote candidates or absent PMI to authoritative facts. Traceability requires evidence on every edge and reports orphans. Chrono is isolated as an optional runtime and reports BLOCKED when unavailable.

The implementation uses the existing OCP runtime because it is installed and executable. STEPcode, step-p21, Chrono, FreeCAD, CadQuery, Analysis Situs, and BrepMFR are represented explicitly as NOT_AVAILABLE or NOT_INTEGRATED; no random source checkout or binary is silently added.
