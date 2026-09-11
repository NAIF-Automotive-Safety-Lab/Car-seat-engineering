# Webhook Specification

`POST /api/jon/webhook` uses `WEBHOOK_SECRET` from the server environment. Required headers are `x-jon-signature`, `x-jon-timestamp`, `x-jon-nonce`, and optional `x-jon-request-id`.

Signature: `HMAC-SHA256(WEBHOOK_SECRET, timestamp + "." + nonce + "." + raw JSON body)`. Timestamps outside five minutes are `REJECTED`; reused timestamp/nonce pairs are `REPLAYED`; invalid or missing credentials are `UNAUTHORIZED`; valid requests return `ACCEPTED`.
