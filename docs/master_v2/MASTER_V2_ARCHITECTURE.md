# Master v2 Architecture

The repository-native Master v2 architecture preserves the existing EERE and OpenCascade paths and adds typed, fail-closed contracts for evidence, physics, instrumentation, traceability, and deterministic gap closure. Computational evidence and physical evidence are separate classes. Missing V7 native CAD or engineering records remain BLOCKED.

The authoritative flow is acquisition → provenance → STEP/B-Rep → deterministic features/PMI → interfaces/joints/materials/fasteners → physics → instrumentation → evidence graph → gap closure → validation. No later layer can upgrade an earlier evidence class without source identity and verification.
