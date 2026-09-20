# Federation protocol

The server accepts an envelope containing `source_id`, `destination_endpoint`, `protocol_version`, `request_id`, `timestamp`, `nonce`, `payload_hash`, and `signature`. It verifies a SHA-256 HMAC over canonical JSON using `AEGIS_FEDERATION_SHARED_SECRET`, records a nonce in SQLite before accepting it, and returns `ESTABLISHED` only on cryptographic success. Missing secret, malformed envelope, bad signature, and replays are `BLOCKED`. TLS endpoint transport, key rotation, and a real remote counterparty require external deployment/runtime.

The API persists the cryptographically verified result as a `FederationHandshake` audit entity. The contract covers CONNECT request intake, AUTHENTICATE using the server-only shared secret, CHALLENGE/nonce issuance by the caller or peer, canonical signing, server verification, and `ESTABLISHED` only after signature comparison. This repository does not operate a remote federation peer, so no deployment-level handshake is claimed.
