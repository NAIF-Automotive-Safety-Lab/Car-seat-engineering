import { and, eq } from "drizzle-orm";
import { describe, expect, it } from "vitest";
import { appRouter } from "../server/routers";
import { deleteUserByOpenId, getDb, getUserByOpenId, upsertUser } from "../server/db";
import { deleteKernelRunForTest } from "../server/kernel";
import { kernelAudits, kernelEngines } from "../drizzle/schema";
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
});
