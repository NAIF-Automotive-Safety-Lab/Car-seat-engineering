# AEGIS_REPOSITORY_BASELINE

**Inspection timestamp:** 2026-09-16 UTC. **Scope:** repository root, excluding vendored external projects for implementation decisions.

## Current architecture and dependency graph

The repository is a Python engineering-evidence package, not a web application. `engineering_recovery` contains deterministic recovery, B-Rep, STEP, evidence, traceability and runtime probes. Root JSON/CSV/STEP/MD files are existing engineering records; they are not changed by AEGIS. CI consists of EERE tests and an R4.1 MBD workflow. There is no package manifest, frontend entrypoint, database, HTTP server, Firebase configuration, hosting configuration, authentication, federation, or importer service.

`engineering_recovery/runtime.py` probes optional OCP/STEP/Chrono/FreeCAD/CadQuery runtimes. Its MBD CI deliberately reports that real solver work is not started unless requested by validation gates. R4.1 and V7 source records are canonical/review records and remain outside the new `aegis/` implementation path.

## Runtime boundaries and reusable modules

Browser: static `aegis/static/index.html` only. Server: `aegis.server`, SQLite ledger, upload validation, AEOK, and federation verification. Native CAD, mesh and solver work remains external-runtime-only; existing `engineering_recovery` modules remain reusable adapters/probes and are not imported into browser code.

## Baseline blockers and risks

No deployable host or Firebase project is configured; deployment and live smoke tests are therefore **BLOCKED**. No physical test evidence, external CAD/CAE runtime, GPT/DeepSeek credentials, or federation counterparty is present; each is **UNKNOWN** or **EXTERNAL_RUNTIME_REQUIRED**, never promoted. The prior broad evidence corpus contains sensitive engineering material, so this implementation stores only references/hashes and never embeds secrets in the static bundle.
