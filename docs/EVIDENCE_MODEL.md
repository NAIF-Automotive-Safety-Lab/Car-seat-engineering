# Evidence model

The persistent entity registry accepts the required project, engineering, CAE, test, compliance, audit, worker, federation, hash, signature, provenance and release entity types. Each inserted entity has immutable UUID identity, version, status, creator/time, provenance, and SHA-256 content hash. `input_refs` supply the evidence-graph edges; graph rendering is intentionally a query concern over persisted records rather than a second mutable truth store.

Upload records retain filename, byte size, claimed media type, metadata-parser status, and byte hash. Native STEP/STL/OBJ parsing is marked `EXTERNAL_RUNTIME_REQUIRED`, not represented as validated geometry.

Physical records are only created by the physical-evidence intake route and carry `evidence_class = MEASURED_PHYSICALLY`; correlation is a separate `CAEResult` record with `kind = CORRELATION`. A `ChangeSet` analysis result retains the precise affected IDs/hashes by engineering category. Adversarial provider records retain provider/model/request ID, prompt-context SHA-256, artifact references, raw output or external-runtime result, verdict, confidence, and evidence references.
