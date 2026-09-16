# V7-AEGIS engineering orchestration architecture

AEGIS is a server-side engineering organization OS layered beside—not inside—the immutable V7/R4.1 evidence corpus. Its deterministic core owns lifecycle and release gates. The persistent SQLite ledger is the evidence graph: requirements, artifacts, jobs, reviews, contradictions, compliance and release records join through immutable IDs, SHA-256 values and provenance input references.

The worker fabric boundary is represented by persistent `WorkerJob`, `SolverRun`, `Mesh`, `FEModel`, and importer records. It has no in-process solver substitute. CAD/mesh/FE execution must be invoked by an external sandboxed runtime and return hashed output records. The evidence firewall admits only validated uploads and ledger records; release evaluation is performed only by AEOK. Change impact is derived from transitive `input_refs` query processing and must be recorded as a `ChangeSet`; it cannot silently approve a release.

Adversarial reviews are persistent `AdversarialReview` records and contradictions are persistent `Contradiction` records. Provider execution is intentionally unavailable until authenticated server-only GPT and DeepSeek adapters are deployed; their status remains `UNKNOWN`, and unresolved reviews block release. Federation is server-only signed HMAC verification with nonce replay prevention.

The browser is a static dashboard. It never imports server code and it never contains secrets. The Python standard-library HTTP service exposes the authenticated API. The API records real database state only—there are no timer-based progress states or fake execution results.

## Closure implementation

`aegis.workflows` now executes deterministic change-impact traversal over persisted provenance references and emits a successor `ChangeSet` with grouped requirement, CAD, CAE, test, compliance, and release-gate impacts. It compiles a requirement with explicit acceptance criteria into a machine-readable `TestProtocol`; physical intake requires test identity, source, capture time, measurement references, safe file validation, and a hash. CAE-to-test correlation accepts only a `CAEResult`/`SolverRun` and a `TestMeasurement` marked `MEASURED_PHYSICALLY`; it remains `UNKNOWN` without matching calibrated numeric channels.

Worker jobs persist `QUEUED` and explicit `BLOCKED` external-runtime outcomes. The GPT and DeepSeek adapters use separate server-only URL/key/model environment variables and persist request/context hashes and outputs. Missing credentials or provider errors persist `EXTERNAL_RUNTIME_REQUIRED`; no response is invented. The contradiction engine compares persisted reviews and records `AGREE`, `DISAGREE`, `INSUFFICIENT_EVIDENCE`, or `UNRESOLVED`; every classification other than `AGREE` blocks the recorded contradiction.
