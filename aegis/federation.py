from __future__ import annotations
import hashlib, hmac
from datetime import datetime, timezone
from typing import Any
from .core import AegisStore, canonical_json
REQUIRED = frozenset({"source_id", "destination_endpoint", "protocol_version", "request_id", "timestamp", "nonce", "payload_hash", "signature"})
def sign(envelope: dict[str, Any], secret: str) -> str:
    payload = {k:v for k,v in envelope.items() if k != "signature"}
    return hmac.new(secret.encode(), canonical_json(payload), hashlib.sha256).hexdigest()
def verify_handshake(envelope: dict[str, Any], secret: str, store: AegisStore) -> str:
    if not REQUIRED <= envelope.keys() or not secret: return "BLOCKED"
    try: datetime.fromisoformat(envelope["timestamp"].replace("Z", "+00:00"))
    except (ValueError, TypeError): return "BLOCKED"
    if not store.consume_nonce(str(envelope["nonce"])): return "BLOCKED"
    return "ESTABLISHED" if hmac.compare_digest(str(envelope["signature"]), sign(envelope, secret)) else "BLOCKED"
