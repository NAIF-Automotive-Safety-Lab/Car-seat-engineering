import { and, eq } from "drizzle-orm";
import { drizzle } from "drizzle-orm/mysql2";
import { describe, expect, it } from "vitest";
import { appRouter } from "../server/routers";
import { deleteUserByOpenId, getDb, getUserByOpenId, upsertUser } from "../server/db";
import { deleteKernelRunForTest } from "../server/kernel";
import { kernelAudits, kernelEngines, kernelRuns } from "../drizzle/schema";
import type { TrpcContext } from "../server/_core/context";

type TestUser = NonNullable<TrpcContext["user"]>;

function context(user: TestUser | null): TrpcContext {
  return {
    user,
    req: { protocol: "https", hostname: "localhost", headers: {} } as TrpcContext["req"],
    res: {} as TrpcContext["res"],
  };
}

async function testUser() {
  const openId = `aegis-kernel-${Date.now()}-${Math.random().toString(36).slice(2)}`;
  await upsertUser({ openId, name: "AEGIS-X Kernel Test", email: `${openId}@example.invalid`, loginMethod: "kernel-test" });
  const user = await getUserByOpenId(openId);
  if (!user) throw new Error("kernel test user was not persisted");
  return { openId, user };
}

const validPlan = {
  testId: "R4.1-QUALIFICATION",
  testVersion: "1.0.0",
  repository: "https://github.com/NAIF-Automotive-Safety-Lab/Car-seat-engineering.git",
  commit: "0123456789abcdef0123456789abcdef01234567",
  artifactSha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
  model: "R4.1",
  inputHashes: { artifact: "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb" },
  engineId: "opencascade",
  engineVersion: "1.0.0",
};

describe("AEGIS-X kernel trust and execution gates", () => {
  it("blocks R4.1 before execution when upstream trust is missing", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    let runId = "";
    try {
      const planned = await caller.kernel.plan({ ...validPlan, artifactSha256: undefined });
      runId = planned.runId;
      expect(planned.state).toBe("BLOCKED");
      expect(planned.blockReason).toContain("missing artifact SHA");
    } finally {
      if (runId) await deleteKernelRunForTest(user.id, runId);
      await deleteUserByOpenId(openId);
    }
  });

  it("requires a registered engine, evidence, and audit before a PASS can be recorded", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    let runId = "";
    try {
      await caller.kernel.registerEngine({ engineId: "opencascade", name: "OpenCascade", version: "1.0.0", adapterKind: "CAD_KERNEL", capabilities: ["STEP", "B-REP", "TOPOLOGY"] });
      const planned = await caller.kernel.plan(validPlan);
      runId = planned.runId;
      expect(planned.state).toBe("QUEUED");
      await expect(caller.kernel.complete({ runId, outputHash: "output-hash", result: { solids: "not_proven" }, resultStatus: "PASS", evidenceIds: [] })).rejects.toThrow("PASS requires evidence");

      const evidence = await caller.kernel.recordEvidence({ runId, evidenceClass: "ARTIFACT_PROOF", sourceHash: validPlan.artifactSha256, payload: { source: "test-fixture" }, provenance: { repository: validPlan.repository, commit: validPlan.commit } });
      expect(evidence.evidenceId).toMatch(/^k_ev_/);
      const completed = await caller.kernel.complete({ runId, outputHash: "output-hash", result: { solids: "not_proven" }, resultStatus: "PASS", evidenceIds: [evidence.evidenceId] });
      expect(completed.state).toBe("COMPLETED");
      await expect(caller.kernel.complete({ runId, outputHash: "output-hash-2", result: { duplicate: true }, resultStatus: "PASS", evidenceIds: [evidence.evidenceId] })).rejects.toThrow("run is not executable");

      const readBack = await caller.kernel.getRun({ runId });
      expect(readBack?.run.state).toBe("COMPLETED");
      expect(readBack?.run.gate).toBe("PASS");
      expect(readBack?.evidence).toHaveLength(1);
      expect(readBack?.audits.map((row) => row.eventType)).toEqual(expect.arrayContaining(["EXECUTION_QUEUED", "EVIDENCE_RECORDED", "RESULT_RECORDED"]));
    } finally {
      if (runId) await deleteKernelRunForTest(user.id, runId);
      const db = await getDb();
      if (db) {
        await db.delete(kernelAudits).where(and(eq(kernelAudits.actorUserId, user.id), eq(kernelAudits.eventType, "ENGINE_REGISTERED")));
        await db.delete(kernelEngines).where(eq(kernelEngines.engineId, "opencascade"));
      }
      await deleteUserByOpenId(openId);
    }
  });

  it("rejects unauthenticated planning and keeps adapters separate from solver success", async () => {
    const anonymous = appRouter.createCaller(context(null));
    await expect(anonymous.kernel.plan(validPlan)).rejects.toMatchObject({ code: "UNAUTHORIZED" });
    expect(validPlan.engineId).toBe("opencascade");
  });

  it("blocks engine execution when no real solver runtime is available", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    let runId = "";
    try {
      await caller.kernel.registerEngine({ engineId: "opencascade", name: "OpenCascade", version: "1.0.0", adapterKind: "CAD_KERNEL", capabilities: ["STEP"] });
      const planned = await caller.kernel.plan(validPlan);
      runId = planned.runId;
      const executed = await caller.kernel.execute({ runId });
      expect(executed.state).toBe("BLOCKED");
      expect(executed.reason).toContain("ENGINE_RUNTIME_NOT_AVAILABLE");
      await expect(caller.kernel.cancel({ runId })).rejects.toThrow("terminal run cannot be cancelled");
    } finally {
      if (runId) await deleteKernelRunForTest(user.id, runId);
      const db = await getDb();
      if (db) {
        await db.delete(kernelAudits).where(and(eq(kernelAudits.actorUserId, user.id), eq(kernelAudits.eventType, "ENGINE_REGISTERED")));
        await db.delete(kernelEngines).where(eq(kernelEngines.engineId, "opencascade"));
      }
      await deleteUserByOpenId(openId);
    }
  });

  it("detects wrong SHA and compares reproducibility without inventing a result", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    let runId = "";
    try {
      await caller.kernel.registerEngine({ engineId: "opencascade", name: "OpenCascade", version: "1.0.0", adapterKind: "CAD_KERNEL", capabilities: ["STEP"] });
      const wrongSha = await caller.kernel.plan({ ...validPlan, artifactSha256: "not-a-sha" });
      runId = wrongSha.runId;
      expect(wrongSha.state).toBe("BLOCKED");
      expect(wrongSha.blockReason).toContain("artifact SHA format is invalid");
      const comparison = await caller.kernel.reproducibility({ first: { repository: validPlan.repository, commit: validPlan.commit, artifactSha: validPlan.artifactSha256 }, second: { repository: validPlan.repository, commit: validPlan.commit, artifactSha: validPlan.artifactSha256 } });
      expect(comparison.reproducible).toBe(true);
      expect(comparison.differences).toHaveLength(0);
      const gaps = await caller.kernel.gaps();
      expect(gaps.some((gap) => gap.id === "K-GAP-BLOCKED-RUNS")).toBe(true);
    } finally {
      if (runId) await deleteKernelRunForTest(user.id, runId);
      const db = await getDb();
      if (db) {
        await db.delete(kernelAudits).where(and(eq(kernelAudits.actorUserId, user.id), eq(kernelAudits.eventType, "ENGINE_REGISTERED")));
        await db.delete(kernelEngines).where(eq(kernelEngines.engineId, "opencascade"));
      }
      await deleteUserByOpenId(openId);
    }
  });

  it("denies forged ownership and protected baseline mutation", async () => {
    const owner = await testUser();
    const forger = await testUser();
    const ownerCaller = appRouter.createCaller(context(owner.user));
    const forgerCaller = appRouter.createCaller(context(forger.user));
    let runId = "";
    try {
      const blocked = await ownerCaller.kernel.plan({ ...validPlan, artifactSha256: undefined });
      runId = blocked.runId;
      await expect(forgerCaller.kernel.getRun({ runId })).resolves.toBeUndefined();
      await expect(forgerCaller.kernel.execute({ runId })).rejects.toThrow("run not found or not owned");
      await expect(ownerCaller.kernel.registerTest({ testId: "MUTATION-ATTEMPT", version: "1.0.0", purpose: "negative", preconditions: [], inputs: [], inputHashes: [], engineId: "opencascade", engineVersion: "1.0.0", command: "deny", checks: [], acceptance: [], outputs: [], evidenceClass: "NONE", gatePolicy: "BLOCK", protectedArtifact: "R4.1" })).rejects.toThrow("protected baseline mutation denied");
    } finally {
      if (runId) await deleteKernelRunForTest(owner.user.id, runId);
      await deleteUserByOpenId(owner.openId);
      await deleteUserByOpenId(forger.openId);
    }
  });

  it("claims a queued run once with a lease and isolated resource class", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    let runId = "";
    try {
      await caller.kernel.registerEngine({ engineId: "opencascade", name: "OpenCascade", version: "1.0.0", adapterKind: "CAD_KERNEL", capabilities: ["STEP"] });
      const planned = await caller.kernel.plan(validPlan);
      runId = planned.runId;
      const claim = await caller.kernel.claim({ workerId: "worker-test-1", leaseSeconds: 15 });
      expect(claim.claimed).toBe(true);
      if (!claim.claimed) throw new Error("queue claim did not return a run");
      expect(claim.runId).toBe(runId);
      const readBack = await caller.kernel.getRun({ runId });
      expect(readBack?.run.state).toBe("RUNNING");
      expect(readBack?.run.attempts).toBe(1);
      expect(readBack?.run.resourceClass).toBe("isolated");
      expect(readBack?.run.leaseOwner).toBe("worker-test-1");
      const secondClaim = await caller.kernel.claim({ workerId: "worker-test-2", leaseSeconds: 15 });
      expect(secondClaim.claimed).toBe(false);
    } finally {
      if (runId) await deleteKernelRunForTest(user.id, runId);
      const db = await getDb();
      if (db) {
        await db.delete(kernelAudits).where(and(eq(kernelAudits.actorUserId, user.id), eq(kernelAudits.eventType, "ENGINE_REGISTERED")));
        await db.delete(kernelEngines).where(eq(kernelEngines.engineId, "opencascade"));
      }
      await deleteUserByOpenId(openId);
    }
  });

  it("executes all six R4.1 pipeline steps only as isolated NOT_PROVEN fixture runs", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    const runIds: string[] = [];
    try {
      const result = await caller.kernel.r41Fixture({ fixturePayload: "isolated-r4.1-step-fixture" });
      expect(result.realArtifactStatus).toBe("BLOCKED / NOT_PROVEN");
      expect(result.records).toHaveLength(6);
      for (const record of result.records) {
        runIds.push(record.runId);
        expect(record.testId).toMatch(/^R4\.1-00[1-6]$/);
        expect(record.version).toBe("1.0.0");
        expect(record.artifactSha).toHaveLength(64);
        expect(record.commit).toBe("fixture-isolated");
        expect(record.engineVersion).toBe("0.0.0-fixture");
        expect(record.outputHash).toHaveLength(64);
        expect(record.result).toBe("NOT_PROVEN");
        expect(record.evidenceId).toMatch(/^r41_ev_/);
        expect(record.auditId).toMatch(/^r41_audit_/);
        expect(record.gate).toBe("NOT_PROVEN");
      }
    } finally {
      for (const runId of runIds) await deleteKernelRunForTest(user.id, runId);
      await deleteUserByOpenId(openId);
    }
  });

  it("reads a persisted run through a fresh database connection", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    let runId = "";
    try {
      await caller.kernel.registerEngine({ engineId: "opencascade", name: "OpenCascade", version: "1.0.0", adapterKind: "CAD_KERNEL", capabilities: ["STEP"] });
      const planned = await caller.kernel.plan(validPlan);
      runId = planned.runId;
      await caller.kernel.claim({ workerId: "restart-test-worker", leaseSeconds: 15 });
      if (!process.env.DATABASE_URL) throw new Error("DATABASE_URL is required for restart persistence coverage");
      const freshDb = drizzle(process.env.DATABASE_URL);
      const persisted = await freshDb.select().from(kernelRuns).where(eq(kernelRuns.runId, runId)).limit(1);
      expect(persisted[0]?.state).toBe("RUNNING");
      expect(persisted[0]?.attempts).toBe(1);
      expect(persisted[0]?.leaseOwner).toBe("restart-test-worker");
      await freshDb.$client.end();
    } finally {
      if (runId) await deleteKernelRunForTest(user.id, runId);
      const db = await getDb();
      if (db) {
        await db.delete(kernelAudits).where(and(eq(kernelAudits.actorUserId, user.id), eq(kernelAudits.eventType, "ENGINE_REGISTERED")));
        await db.delete(kernelEngines).where(eq(kernelEngines.engineId, "opencascade"));
      }
      await deleteUserByOpenId(openId);
    }
  });

  it("runs one isolated worker tick and marks dependent runs stale", async () => {
    const { openId, user } = await testUser();
    const caller = appRouter.createCaller(context(user));
    let runId = "";
    try {
      await caller.kernel.registerEngine({ engineId: "opencascade", name: "OpenCascade", version: "1.0.0", adapterKind: "CAD_KERNEL", capabilities: ["STEP"] });
      const planned = await caller.kernel.plan(validPlan);
      runId = planned.runId;
      const tick = await caller.kernel.workerTick({ workerId: "worker-tick-test", leaseSeconds: 15 });
      expect(tick.claimed).toBe(true);
      if (!tick.claimed) throw new Error("worker did not claim a queued run");
      expect(tick.execution.state).toBe("BLOCKED");
      const impact = await caller.kernel.dependencyImpact({ changedDependencies: ["R4.1-artifact"], tests: [{ testId: "R4.1-QUALIFICATION", dependencies: ["R4.1-artifact", "repository"] }, { testId: "OTHER", dependencies: ["unrelated"] }], runs: [{ runId, testId: "R4.1-QUALIFICATION", state: "BLOCKED" }] });
      expect(impact.affectedTests).toEqual(["R4.1-QUALIFICATION"]);
      expect(impact.staleRuns[0]).toMatchObject({ runId, status: "STALE" });
      expect(impact.rerunRequired).toBe(true);
    } finally {
      if (runId) await deleteKernelRunForTest(user.id, runId);
      const db = await getDb();
      if (db) {
        await db.delete(kernelAudits).where(and(eq(kernelAudits.actorUserId, user.id), eq(kernelAudits.eventType, "ENGINE_REGISTERED")));
        await db.delete(kernelEngines).where(eq(kernelEngines.engineId, "opencascade"));
      }
      await deleteUserByOpenId(openId);
    }
  });
});
