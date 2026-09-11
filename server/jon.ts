import crypto from "node:crypto";
import type { Express, NextFunction, Request, Response } from "express";
import { z } from "zod";
import { getKernelRun, listKernelStatus, planKernelRun } from "./kernel";

const replayCache = new Map<string, number>();
const WINDOW_MS = 5 * 60 * 1000;

function configuredSecret(name: "JON_API_KEY" | "JON_SERVICE_TOKEN" | "WEBHOOK_SECRET") {
  return process.env[name] ?? "";
}

function safeEqual(a: string, b: string) {
  const left = Buffer.from(a);
  const right = Buffer.from(b);
  return left.length === right.length && crypto.timingSafeEqual(left, right);
}

export function buildWebhookSignature(secret: string, timestamp: string, nonce: string, body: string) {
  return crypto.createHmac("sha256", secret).update(`${timestamp}.${nonce}.${body}`).digest("hex");
}

export function verifyWebhookSignature(input: { secret: string; signature?: string; timestamp?: string; nonce?: string; body: string; now?: number }) {
  const timestamp = Number(input.timestamp);
  if (!input.secret || !input.signature || !input.timestamp || !input.nonce || !Number.isFinite(timestamp)) return "UNAUTHORIZED" as const;
  if (Math.abs((input.now ?? Date.now()) - timestamp) > WINDOW_MS) return "REJECTED" as const;
  const replayKey = `${input.timestamp}:${input.nonce}`;
  const cachedUntil = replayCache.get(replayKey);
  if (cachedUntil && cachedUntil > (input.now ?? Date.now())) return "REPLAYED" as const;
  const expected = buildWebhookSignature(input.secret, input.timestamp, input.nonce, input.body);
  if (!safeEqual(expected, input.signature)) return "UNAUTHORIZED" as const;
  replayCache.set(replayKey, (input.now ?? Date.now()) + WINDOW_MS);
  for (const [key, expiresAt] of replayCache) if (expiresAt <= (input.now ?? Date.now())) replayCache.delete(key);
  return "ACCEPTED" as const;
}

function requireJonApi(req: Request, res: Response, next: NextFunction) {
  const configured = [configuredSecret("JON_API_KEY"), configuredSecret("JON_SERVICE_TOKEN")].filter(Boolean);
  if (configured.length === 0) {
    res.status(503).json({ status: "BLOCKED", reason: "JON credentials are not configured" });
    return;
  }
  const token = req.header("authorization")?.replace(/^Bearer\s+/i, "") ?? "";
  if (!configured.some((secret) => safeEqual(secret, token))) {
    res.status(401).json({ status: "UNAUTHORIZED" });
    return;
  }
  next();
}

const submitSchema = z.object({
  testId: z.string().min(1),
  testVersion: z.string().min(1),
  repository: z.string().optional(),
  commit: z.string().optional(),
  artifactSha256: z.string().optional(),
  model: z.string().optional(),
  inputHashes: z.record(z.string(), z.string()).optional(),
  engineId: z.string().min(1),
  engineVersion: z.string().min(1),
});

function serviceUserId() {
  const value = Number(process.env.JON_SERVICE_USER_ID ?? "");
  return Number.isInteger(value) && value > 0 ? value : null;
}

export function registerJonRoutes(app: Express) {
  app.get("/api/jon/v1/engines", requireJonApi, async (_req, res) => {
    try {
      const status = await listKernelStatus();
      res.json({ status: "PASS", engines: status.engines, adapters: status.adapters });
    } catch {
      res.status(503).json({ status: "BLOCKED", reason: "kernel status unavailable" });
    }
  });

  app.get("/api/jon/v1/artifacts", requireJonApi, (_req, res) => {
    res.status(200).json({ status: "BLOCKED", artifacts: [], reason: "No trusted artifact resolver is configured" });
  });

  app.post("/api/jon/v1/submit-test", requireJonApi, async (req, res) => {
    const userId = serviceUserId();
    if (!userId) {
      res.status(503).json({ status: "BLOCKED", reason: "JON_SERVICE_USER_ID is not configured" });
      return;
    }
    const parsed = submitSchema.safeParse(req.body);
    if (!parsed.success) {
      res.status(400).json({ status: "REJECTED", reason: "invalid test request" });
      return;
    }
    try {
      const planned = await planKernelRun({ userId, ...parsed.data });
      res.status(202).json({ jobId: planned.runId, status: planned.state === "QUEUED" ? "QUEUED" : "BLOCKED", gateId: planned.gateId, blockReason: planned.blockReason });
    } catch {
      res.status(503).json({ status: "BLOCKED", reason: "AEGIS-X unavailable" });
    }
  });

  app.get("/api/jon/v1/jobs/:runId", requireJonApi, async (req, res) => {
    const userId = serviceUserId();
    if (!userId) {
      res.status(503).json({ status: "BLOCKED", reason: "JON_SERVICE_USER_ID is not configured" });
      return;
    }
    const run = await getKernelRun(userId, req.params.runId);
    if (!run) {
      res.status(404).json({ status: "REJECTED", reason: "run not found" });
      return;
    }
    res.json({ status: run.run.state, run: run.run, evidence: run.evidence, audits: run.audits, gates: run.gates });
  });

  app.post("/api/jon/webhook", (req, res) => {
    const result = verifyWebhookSignature({ secret: configuredSecret("WEBHOOK_SECRET"), signature: req.header("x-jon-signature"), timestamp: req.header("x-jon-timestamp"), nonce: req.header("x-jon-nonce"), body: JSON.stringify(req.body) });
    if (result !== "ACCEPTED") {
      res.status(result === "REPLAYED" ? 409 : result === "UNAUTHORIZED" ? 401 : 400).json({ status: result });
      return;
    }
    res.status(202).json({ status: "ACCEPTED", requestId: req.header("x-jon-request-id") ?? null });
  });
}
