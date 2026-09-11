import { describe, expect, it } from "vitest";
import { buildWebhookSignature, verifyWebhookSignature } from "../server/jon";

describe("JON connectivity fabric security", () => {
  it("accepts a correctly signed webhook once and rejects replay", () => {
    const secret = "test-webhook-secret";
    const body = JSON.stringify({ runId: "run-1", status: "BLOCKED" });
    const timestamp = String(1_700_000_000_000);
    const nonce = "nonce-1";
    const signature = buildWebhookSignature(secret, timestamp, nonce, body);
    const accepted = verifyWebhookSignature({ secret, signature, timestamp, nonce, body, now: 1_700_000_001_000 });
    expect(accepted).toBe("ACCEPTED");
    const replayed = verifyWebhookSignature({ secret, signature, timestamp, nonce, body, now: 1_700_000_001_500 });
    expect(replayed).toBe("REPLAYED");
  });

  it("rejects forged, expired, and secretless webhook requests", () => {
    const body = "{}";
    expect(verifyWebhookSignature({ secret: "secret", signature: "forged", timestamp: "1700000000000", nonce: "nonce-forged", body, now: 1700000001000 })).toBe("UNAUTHORIZED");
    expect(verifyWebhookSignature({ secret: "secret", signature: "forged", timestamp: "1700000000000", nonce: "nonce-expired", body, now: 1700001000000 })).toBe("REJECTED");
    expect(verifyWebhookSignature({ secret: "", signature: undefined, timestamp: undefined, nonce: undefined, body })).toBe("UNAUTHORIZED");
  });
});
