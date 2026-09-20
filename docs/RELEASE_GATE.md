# Release gate

Release decisions use only `PASS`, `BLOCK`, and `UNKNOWN`. The implemented gate is fail-closed: absent, unknown, malformed, or non-PASS mandatory checks produce `BLOCK`. It explicitly does not equate a claim, simulation, mock, or unknown value to evidence or physical test. A persisted gate is an auditable `ReleaseGate` record; no HTTP endpoint performs a release transition automatically.

Required release checks also include `frontend_backend_boundary`, `mandatory_evidence`, `contradictions`, `critical_errors`, and `broken_paths`. `contradictions` must be `PASS`; an `INSUFFICIENT_EVIDENCE`, `DISAGREE`, or `UNRESOLVED` contradiction therefore cannot be released.
