import express from "express";
import { afterEach, describe, expect, it, vi } from "vitest";
import { registerOAuthRoutes } from "../server/_core/oauth";
import { sdk } from "../server/_core/sdk";
import { deleteUserByOpenId } from "../server/db";

function startTestServer() {
  const app = express();
  app.use(express.json());
  registerOAuthRoutes(app);
  return new Promise<{ server: ReturnType<typeof app.listen>; baseUrl: string }>((resolve) => {
    const server = app.listen(0, "127.0.0.1", () => {
      const address = server.address();
      if (!address || typeof address === "string") throw new Error("test server did not bind");
      resolve({ server, baseUrl: `http://127.0.0.1:${address.port}` });
    });
  });
}

describe("OAuth callback security", () => {
  afterEach(() => vi.restoreAllMocks());

  it("rejects missing code", async () => {
    const { server, baseUrl } = await startTestServer();
    try {
      const response = await fetch(`${baseUrl}/api/oauth/mobile?state=${encodeURIComponent(Buffer.from("https://example.test/callback").toString("base64"))}`);
      expect(response.status).toBe(400);
      expect(await response.json()).toEqual({ error: "code and state are required" });
    } finally {
      server.close();
    }
  });

  it("rejects missing state and sessionToken URL injection", async () => {
    const { server, baseUrl } = await startTestServer();
    try {
      const missingState = await fetch(`${baseUrl}/api/oauth/mobile?code=legitimate-code&sessionToken=forged-token`);
      expect(missingState.status).toBe(400);
      expect(await missingState.json()).toEqual({ error: "code and state are required" });
    } finally {
      server.close();
    }
  });

  it("rejects invalid state before token exchange", async () => {
    const exchange = vi.spyOn(sdk, "exchangeCodeForToken");
    const { server, baseUrl } = await startTestServer();
    try {
      const response = await fetch(`${baseUrl}/api/oauth/mobile?code=legitimate-code&state=not-a-valid-state`);
      expect(response.status).toBe(400);
      expect(await response.json()).toEqual({ error: "invalid OAuth state" });
      expect(exchange).not.toHaveBeenCalled();
    } finally {
      server.close();
    }
  });

  it("accepts a valid code/state exchange and returns a session cookie", async () => {
    const openId = `oauth-security-${Date.now()}-${Math.random().toString(36).slice(2)}`;
    const state = Buffer.from("https://example.test/oauth/callback").toString("base64");
    vi.spyOn(sdk, "exchangeCodeForToken").mockResolvedValue({ accessToken: "provider-access-token" } as never);
    vi.spyOn(sdk, "getUserInfo").mockResolvedValue({ openId, name: "OAuth Test", email: `${openId}@example.invalid`, loginMethod: "test" } as never);
    vi.spyOn(sdk, "createSessionToken").mockResolvedValue("server-created-session-token");
    const { server, baseUrl } = await startTestServer();
    try {
      const response = await fetch(`${baseUrl}/api/oauth/mobile?code=legitimate-code&state=${encodeURIComponent(state)}`);
      expect(response.status).toBe(200);
      const body = await response.json();
      expect(body.app_session_id).toBe("server-created-session-token");
      expect(body.user.openId).toBe(openId);
      expect(response.headers.get("set-cookie")).toContain("app_session_id=");
    } finally {
      server.close();
      await deleteUserByOpenId(openId);
    }
  });

  it("returns 401 without a session", async () => {
    const { server, baseUrl } = await startTestServer();
    try {
      const response = await fetch(`${baseUrl}/api/auth/me`);
      expect(response.status).toBe(401);
      expect(await response.json()).toEqual({ error: "Not authenticated", user: null });
    } finally {
      server.close();
    }
  });
});
