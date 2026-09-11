# JON API

The JON gateway is a restricted third-party boundary around AEGIS-X. It never writes results directly, never executes shell commands, and cannot modify protected baselines.

## Authentication

Use `Authorization: Bearer <JON_API_KEY>` or `Authorization: Bearer <JON_SERVICE_TOKEN>`. Secrets are read only from server environment variables and are never returned or logged.

## Endpoints

`GET /api/jon/v1/engines` lists registered engines and adapter descriptors. `GET /api/jon/v1/artifacts` returns `BLOCKED` until a trusted artifact resolver is configured. `POST /api/jon/v1/submit-test` validates a request and forwards it to AEGIS-X planning; it returns `QUEUED` or `BLOCKED`, never a direct solver result. `GET /api/jon/v1/jobs/:runId` returns an owned run and its persisted evidence, audit, and gate records. `POST /api/jon/webhook` accepts only signed callbacks.

A real artifact is required for an engineering PASS. Without it, R4.1 remains `BLOCKED / NOT_PROVEN`.
