import { describe, expect, it } from "vitest";
import { appRouter } from "../server/routers";
import { deleteUserByOpenId, getUserByOpenId, upsertUser } from "../server/db";
import { deleteQualificationRunForTest } from "../server/qualification";
import type { TrpcContext } from "../server/_core/context";

type TestUser = NonNullable<TrpcContext["user"]>;

function createContext(user: TestUser | null): TrpcContext {
  return {
    user,
    req: { protocol: "https", hostname: "localhost", headers: {} } as TrpcContext["req"],
    res: {} as TrpcContext["res"],
  };
}

async function createTestUser() {
  const openId = `aegis-run-test-${Date.now()}-${Math.random().toString(36).slice(2)}`;
  await upsertUser({ openId, name: "Qualification Test User", email: `${openId}@example.invalid`, loginMethod: "integration-test" });
  const user = await getUserByOpenId(openId);
  if (!user) throw new Error("test user was not persisted");
  return { openId, user };
}

const runInput = {
  executor: "USER_ID",
  testId: "R4.1-IDENTITY",
  testVersion: "1.0.0",
  repository: "https://github.com/NAIF-Automotive-Safety-Lab/Car-seat-engineering.git",
  branch: "main",
  commit: "0123456789abcdef0123456789abcdef01234567",
  modelId: "R4.1",
  modelVersion: "baseline",
  modelSha256: "unverified",
  inputs: { source: "integration-test" },
  engine: "Artifact Hash / Evidence",
  engineVersion: "0.1.0",
  command: "verify artifact identity",
};

describe("protected qualification run lifecycle", () => {
  it("creates a run, stores a result, and reads the persisted result and audit trail", async () => {
    const { openId, user } = await createTestUser();
    const caller = appRouter.createCaller(createContext(user));
    let runId = "";
    try {
      const created = await caller.qualification.create({ ...runInput, executor: String(user.id) });
      runId = created.runId;
      expect(created.state).toBe("CREATED");

      const completed = await caller.qualification.complete({ runId, result: { identity: "PASS" }, output: { sha256: "unverified" }, outputHash: "output-hash-from-test", resultStatus: "PASS", evidenceIds: ["EVID-TEST-001"] });
      expect(completed.state).toBe("COMPLETED");

      const reread = await caller.qualification.get({ runId });
      expect(reread?.run.runId).toBe(runId);
      expect(reread?.run.state).toBe("COMPLETED");
      expect(reread?.result?.resultStatus).toBe("PASS");
      expect(reread?.result?.outputHash).toBe("output-hash-from-test");

      const audits = await caller.qualification.audits({ runId });
      expect(audits.map((audit) => audit.eventType)).toEqual(expect.arrayContaining(["RUN_CREATED", "RESULT_STORED", "RUN_COMPLETED"]));
    } finally {
      if (runId) await deleteQualificationRunForTest(user.id, runId);
      await deleteUserByOpenId(openId);
    }
  });

  it("rejects unauthenticated creation and forged executor identity", async () => {
    const anonymousCaller = appRouter.createCaller(createContext(null));
    await expect(anonymousCaller.qualification.create({ ...runInput, executor: "1" })).rejects.toMatchObject({ code: "UNAUTHORIZED" });

    const { openId, user } = await createTestUser();
    try {
      const caller = appRouter.createCaller(createContext(user));
      await expect(caller.qualification.create({ ...runInput, executor: "forged-user" })).rejects.toThrow("executor does not match authenticated user");
      await expect(caller.qualification.create({ ...runInput, repository: "https://example.invalid/not-canonical", executor: String(user.id) })).rejects.toThrow("repository binding is invalid");
    } finally {
      await deleteUserByOpenId(openId);
    }
  });

  it("rejects duplicate completion of a terminal run", async () => {
    const { openId, user } = await createTestUser();
    const caller = appRouter.createCaller(createContext(user));
    let runId = "";
    try {
      const created = await caller.qualification.create({ ...runInput, executor: String(user.id) });
      runId = created.runId;
      const result = { result: { identity: "BLOCKED" }, output: { reason: "blocked" }, outputHash: "blocked-hash", resultStatus: "BLOCKED" as const, evidenceIds: [] };
      await caller.qualification.complete({ runId, ...result });
      await expect(caller.qualification.complete({ runId, ...result })).rejects.toThrow("run is already terminal");
    } finally {
      if (runId) await deleteQualificationRunForTest(user.id, runId);
      await deleteUserByOpenId(openId);
    }
  });
});
