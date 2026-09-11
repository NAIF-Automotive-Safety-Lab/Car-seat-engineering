# Credential Policy

Supported server-only credentials are `JON_API_KEY`, `JON_SERVICE_TOKEN`, `RUNNER_TOKEN`, and `WEBHOOK_SECRET`. They must be supplied by an environment or secret manager, never hardcoded, committed, returned after creation, or exposed to the frontend.

Rotation is performed by replacing the configured secret and restarting the service. Revocation is performed by removing it or setting it to an unavailable value; the gateway then fails closed. OpenAI credentials are not used for CAD/CAE execution.
