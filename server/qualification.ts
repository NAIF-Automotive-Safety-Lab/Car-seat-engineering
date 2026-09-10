import { randomUUID } from "node:crypto";
import { and, desc, eq } from "drizzle-orm";
import { auditRecords, qualificationResults, qualificationRuns } from "../drizzle/schema";
import { getDb } from "./db";

export const CANONICAL_REPOSITORY = "https://github.com/NAIF-Automotive-Safety-Lab/Car-seat-engineering.git";

type CreateRunInput = {
  userId: number;
  executor: string;
  testId: string;
  testVersion: string;
  repository: string;
  branch: string;
  commit: string;
  modelId?: string;
  modelVersion?: string;
  modelSha256?: string;
  inputs: unknown;
  engine: string;
  engineVersion: string;
  command: string;
};

function requireNonEmpty(value: string, field: string) {
  if (!value.trim()) throw new Error(`${field} is required`);
}

function validateCreateInput(input: CreateRunInput) {
  for (const [field, value] of Object.entries(input)) {
    if (["modelId", "modelVersion", "modelSha256"].includes(field)) continue;
    if (typeof value === "string") requireNonEmpty(value, field);
  }
  if (input.repository !== CANONICAL_REPOSITORY) throw new Error("repository binding is invalid");
  if (input.commit.length < 7) throw new Error("exact commit is required");
  if (input.executor !== String(input.userId)) throw new Error("executor does not match authenticated user");
}

async function requireDb() {
  const db = await getDb();
  if (!db) throw new Error("Database is not available");
  return db;
}

export async function createQualificationRun(input: CreateRunInput) {
  validateCreateInput(input);
  const db = await requireDb();
  const runId = `run_${randomUUID()}`;
  const now = new Date();
  await db.insert(qualificationRuns).values({
    runId,
    requestedBy: input.userId,
    executor: input.executor,
    testId: input.testId,
    testVersion: input.testVersion,
    repository: input.repository,
    branch: input.branch,
    commit: input.commit,
    modelId: input.modelId,
    modelVersion: input.modelVersion,
    modelSha256: input.modelSha256,
    inputs: JSON.stringify(input.inputs),
    engine: input.engine,
    engineVersion: input.engineVersion,
    command: input.command,
    state: "CREATED",
    startTime: now,
    createdAt: now,
  });
  const auditId = await writeAudit({ eventType: "RUN_CREATED", actorUserId: input.userId, runId, status: "PASS", payload: { testId: input.testId, repository: input.repository, branch: input.branch, commit: input.commit } });
  return { runId, auditId, state: "CREATED" as const };
}

export async function completeQualificationRun(input: { userId: number; runId: string; result: unknown; output: unknown; outputHash: string; resultStatus: "PASS" | "FAIL" | "BLOCKED" | "NOT_PROVEN"; evidenceIds: string[] }) {
  requireNonEmpty(input.runId, "runId");
  requireNonEmpty(input.outputHash, "outputHash");
  const db = await requireDb();
  const rows = await db.select().from(qualificationRuns).where(and(eq(qualificationRuns.runId, input.runId), eq(qualificationRuns.requestedBy, input.userId))).limit(1);
  const run = rows[0];
  if (!run) throw new Error("run not found or not owned by authenticated user");
  if (run.state === "COMPLETED" || run.state === "FAILED" || run.state === "BLOCKED") throw new Error("run is already terminal");
  const resultId = `result_${randomUUID()}`;
  const auditId = await writeAudit({ eventType: "RESULT_STORED", actorUserId: input.userId, runId: input.runId, status: input.resultStatus, payload: { outputHash: input.outputHash, evidenceIds: input.evidenceIds } });
  await db.insert(qualificationResults).values({ resultId, runId: input.runId, result: JSON.stringify(input.result), output: JSON.stringify(input.output), outputHash: input.outputHash, resultStatus: input.resultStatus, evidenceIds: JSON.stringify(input.evidenceIds), auditId, createdAt: new Date() });
  const nextState = input.resultStatus === "PASS" ? "COMPLETED" : input.resultStatus === "FAIL" ? "FAILED" : "BLOCKED";
  await db.update(qualificationRuns).set({ state: nextState, endTime: new Date() }).where(eq(qualificationRuns.runId, input.runId));
  await writeAudit({ eventType: "RUN_COMPLETED", actorUserId: input.userId, runId: input.runId, status: nextState, payload: { resultId, resultStatus: input.resultStatus } });
  return { resultId, auditId, runId: input.runId, state: nextState };
}

export async function getQualificationRun(userId: number, runId: string) {
  const db = await requireDb();
  const runs = await db.select().from(qualificationRuns).where(and(eq(qualificationRuns.runId, runId), eq(qualificationRuns.requestedBy, userId))).limit(1);
  const run = runs[0];
  if (!run) return undefined;
  const results = await db.select().from(qualificationResults).where(eq(qualificationResults.runId, runId)).limit(1);
  return { run, result: results[0] };
}

export async function listAudits(userId: number, runId?: string) {
  const db = await requireDb();
  if (runId) {
    const owned = await db.select({ runId: qualificationRuns.runId }).from(qualificationRuns).where(and(eq(qualificationRuns.runId, runId), eq(qualificationRuns.requestedBy, userId))).limit(1);
    if (!owned[0]) return [];
  }
  const conditions = runId ? eq(auditRecords.runId, runId) : eq(auditRecords.actorUserId, userId);
  return db.select().from(auditRecords).where(conditions).orderBy(desc(auditRecords.createdAt));
}

async function writeAudit(input: { eventType: string; actorUserId: number; runId?: string; status: string; payload: unknown }) {
  const db = await requireDb();
  const auditId = `audit_${randomUUID()}`;
  await db.insert(auditRecords).values({ auditId, eventType: input.eventType, actorUserId: input.actorUserId, runId: input.runId, status: input.status, payload: JSON.stringify(input.payload), createdAt: new Date() });
  return auditId;
}

export async function deleteQualificationRunForTest(userId: number, runId: string) {
  const db = await requireDb();
  await db.delete(qualificationResults).where(eq(qualificationResults.runId, runId));
  await db.delete(auditRecords).where(eq(auditRecords.runId, runId));
  await db.delete(qualificationRuns).where(and(eq(qualificationRuns.runId, runId), eq(qualificationRuns.requestedBy, userId)));
}
