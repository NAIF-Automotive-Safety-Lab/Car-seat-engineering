# Deployment

Run `AEGIS_AUTH_TOKEN=... AEGIS_FEDERATION_SHARED_SECRET=... python -m aegis.server` behind a TLS reverse proxy, with a writable non-public SQLite volume (`AEGIS_DB`) and a public bind explicitly selected with `AEGIS_HOST`. Firebase Hosting is **NOT_AVAILABLE** for this backend-only repository because no Firebase project/configuration or serverless adapter exists. No deployment URL or deployment evidence has been recorded.
