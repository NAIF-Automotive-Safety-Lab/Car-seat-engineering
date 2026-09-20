# AEOK specification

`aegis.core.AEOK` is the sole implementation path for lifecycle transitions and release-gate creation. The finite state machine is `DRAFT → ARCHITECTED → JOB_CREATED → EVIDENCE_REGISTERED → REVIEWED → RELEASED`; any state can only take explicitly listed forward transitions, while `RELEASE_BLOCKED` and `RELEASED` are terminal. Invalid transitions raise an error.

The SQLite ledger assigns UUID IDs, version 1 records, parent references for successor versions, UTC timestamp, actor, input references, canonical JSON SHA-256, and append-only event entries. Release gates require persisted `PASS` for build, tests, security, AEOK, evidence registry, database, federation, and routes; every other condition is `BLOCK`.

Release evaluation additionally requires `frontend_backend_boundary`, `mandatory_evidence`, and `contradictions` to be `PASS`, plus zero `critical_errors` and `broken_paths`; absent checks block release. It does not treat an externally unavailable runtime, unresolved contradiction, or unknown evidence as a passing check.
