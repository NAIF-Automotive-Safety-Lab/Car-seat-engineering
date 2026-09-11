# JON Authorization Model

JON is an external request principal, not an executor. The gateway maps an authenticated service token to a configured service user and then calls protected AEGIS-X services.

JON may request tests, query runs, read results/evidence/audit/gates, request reproducibility, and request reruns through approved APIs. JON may not mutate databases, execute shell commands, delete baselines, modify R4.1 or V7-R3, bypass gates, or forge results. Missing service configuration fails closed with `BLOCKED`.
