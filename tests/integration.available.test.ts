import { describe, expect, it } from "vitest";
import { deleteUserByOpenId, getDb, getUserByOpenId, upsertUser } from "../server/db";

const apiBaseUrl = process.env.AEGIS_API_URL ?? "http://127.0.0.1:3000";

describe("available V7-AEGIS integration paths", () => {
  it("responds from the live health endpoint", async () => {
    const response = await fetch(`${apiBaseUrl}/api/health`);
    expect(response.ok).toBe(true);
    const body = (await response.json()) as { ok?: boolean; timestamp?: number };
    expect(body.ok).toBe(true);
    expect(typeof body.timestamp).toBe("number");
  });

  it("persists and reads a user record through the configured database", async () => {
    const db = await getDb();
    expect(db).not.toBeNull();

    const openId = `aegis-test-${Date.now()}-${Math.random().toString(36).slice(2)}`;
    try {
      await upsertUser({ openId, name: "AEGIS Integration Test", email: `${openId}@example.invalid`, loginMethod: "integration-test" });
      const saved = await getUserByOpenId(openId);
      expect(saved?.openId).toBe(openId);
      expect(saved?.name).toBe("AEGIS Integration Test");
    } finally {
      await deleteUserByOpenId(openId);
    }
  });
});
