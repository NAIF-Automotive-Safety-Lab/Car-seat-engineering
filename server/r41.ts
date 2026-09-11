import { createHash, randomUUID } from "node:crypto";
import { getDb } from "./db";
import { kernelAudits, kernelEvidence, kernelGates, kernelRuns } from "../drizzle/schema";

export const R41_FIXTURE_TESTS = [
  { testId: "R4.1-001", purpose: "Artifact Resolution" },
  { testId: "R4.1-002", purpose: "STEP Parse" },
  { testId: "R4.1-003", purpose: "B-Rep" },
  { testId: "R4.1-004", purpose: "Immutability" },
  { testId: "R4.1-005", purpose: "Protected Mutation Guard" },
  { testId: "R4.1-006", purpose: "Run→Result→Evidence→Audit→Gate" },
] as const;

function sha256(value: string) {
  return createHash("sha256").update(value, "utf8").digest("hex");
}

async function requireDb() {
  const db = await getDb();
  if (!db) throw new Error("Database is not available");
  return db;
}

export async function qualifyR41IsolatedFixture(input: { userId: number; fixturePayload: string }) {
  const db = await requireDb();
  const fixtureSha = sha256(input.fixturePayload);
  const commit = "fixture-isolated";
  const inputHashes = JSON.stringify({ fixture: fixtureSha });
  const records: Array<{ runId: string; testId: string; version: string; artifactSha: string; commit: string; engineVersion: string; inputHashes: string; outputHash: string; result: "NOT_PROVEN"; evidenceId: string; auditId: string; gate: "NOT_PROVEN" }> = [];
  for (const test of R41_FIXTURE_TESTS) {
    const runId = `r41_fixture_${randomUUID()}`;
    const evidenceId = `r41_ev_${randomUUID()}`;
    const auditId = `r41_audit_${randomUUID()}`;
    const outputHash = sha256(JSON.stringify({ testId: test.testId, fixtureSha, scope: "ISOLATED_FIXTURE" }));
    await db.insert(kernelRuns).values({ runId, requestedBy: input.userId, testId: test.testId, testVersion: "1.0.0", repository: "isolated-fixture://r4.1", commit, artifactSha256: fixtureSha, model: "R4.1-fixture-boundary", inputHashes, engineId: "fixture-boundary", engineVersion: "0.0.0-fixture", state: "NOT_PROVEN", blockReason: "R4.1 REAL ARTIFACT = BLOCKED / NOT_PROVEN", outputHash, result: JSON.stringify({ scope: "ISOLATED_FIXTURE", purpose: test.purpose, realArtifactStatus: "BLOCKED", engineeringResult: "NOT_PROVEN" }), evidenceIds: JSON.stringify([evidenceId]), gate: "NOT_PROVEN", attempts: 0, maxAttempts: 1, resourceClass: "isolated", createdAt: new Date(), completedAt: new Date() });
    await db.insert(kernelEvidence).values({ evidenceId, runId, evidenceClass: "ISOLATED_FIXTURE_PIPELINE", sourceHash: fixtureSha, payload: JSON.stringify({ scope: "ISOLATED_FIXTURE", testId: test.testId, payloadLength: input.fixturePayload.length }), provenance: JSON.stringify({ source: "test-fixture", realArtifactStatus: "BLOCKED", noEngineeringPromotion: true }), createdAt: new Date() });
    await db.insert(kernelGates).values({ gateId: `r41_gate_${randomUUID()}`, runId, policy: "NEVER_PROMOTE_FIXTURE_TO_ENGINEERING_RESULT", decision: "NOT_PROVEN", reason: "Fixture pipeline executed without a real R4.1 artifact.", dependencies: JSON.stringify(["real artifact", "real CAD parser", "real engineering evidence"]), createdAt: new Date() });
    await db.insert(kernelAudits).values({ auditId, runId, eventType: "R41_ISOLATED_FIXTURE_EXECUTED", actorUserId: input.userId, status: "NOT_PROVEN", payload: JSON.stringify({ testId: test.testId, fixtureSha, outputHash }), createdAt: new Date() });
    records.push({ runId, testId: test.testId, version: "1.0.0", artifactSha: fixtureSha, commit, engineVersion: "0.0.0-fixture", inputHashes, outputHash, result: "NOT_PROVEN", evidenceId, auditId, gate: "NOT_PROVEN" });
  }
  return { realArtifactStatus: "BLOCKED / NOT_PROVEN" as const, records };
}
