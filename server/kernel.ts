import { randomUUID } from "node:crypto";
import { and, desc, eq, inArray, isNull, lt, or } from "drizzle-orm";
import { getDb } from "./db";
import { kernelAudits, kernelDrifts, kernelEngines, kernelEvidence, kernelGates, kernelRuns, kernelTests } from "../drizzle/schema";

export const ENGINE_ADAPTERS = [
  { engineId: "opencascade", name: "OpenCascade", adapterKind: "CAD_KERNEL", capabilities: ["STEP", "B-REP", "TOPOLOGY"] },
  { engineId: "openradioss", name: "OpenRadioss", adapterKind: "CAE_SOLVER", capabilities: ["STRUCTURAL"] },
  { engineId: "calculix", name: "CalculiX", adapterKind: "CAE_SOLVER", capabilities: ["FEA"] },
  { engineId: "mbd-future", name: "Future MBD Adapter", adapterKind: "MBD_SOLVER", capabilities: [] },
] as const;
export const CANONICAL_REPOSITORY = "https://github.com/NAIF-Automotive-Safety-Lab/Car-seat-engineering.git";

export const R41_TEST_CONTRACT = {
  testId: "R4.1-QUALIFICATION",
  version: "1.0.0",
  purpose: "Qualify artifact identity, SHA, STEP schema, geometry topology, coverage, immutability, and evidence gates without inventing engineering data.",
  preconditions: ["canonical repository", "exact commit", "exact artifact", "artifact SHA", "model", "inputs", "engine", "authorization"],
  inputs: ["repository", "commit", "artifactSha256", "model", "inputHashes"],
  checks: ["Identity", "SHA", "STEP", "Schema", "Body Records", "Solids", "B-Rep", "Topology", "Coverage", "Immutability", "Mutation Guard", "Evidence", "Audit", "Gate"],
  evidenceClass: "ARTIFACT_PROOF",
  gatePolicy: "BLOCK_ON_UNTRUSTED_UPSTREAM",
} as const;

export type KernelState = "BLOCKED" | "QUEUED" | "RUNNING" | "COMPLETED" | "FAILED" | "NOT_PROVEN";

type PlanInput = {
  userId: number;
  testId: string;
  testVersion: string;
  repository?: string;
  commit?: string;
  artifactSha256?: string;
  model?: string;
  inputHashes?: Record<string, string>;
  engineId: string;
  engineVersion: string;
  maxAttempts?: number;
};

async function requireDb() {
  const db = await getDb();
  if (!db) throw new Error("Database is not available");
  return db;
}

function validationFailures(input: PlanInput, contract?: typeof R41_TEST_CONTRACT) {
  const failures: string[] = [];
  if (!input.testId.trim()) failures.push("missing test id");
  if (!input.testVersion.trim()) failures.push("missing test version");
  if (!input.repository?.trim()) failures.push("missing repository trust");
  else if (input.repository !== CANONICAL_REPOSITORY) failures.push("repository trust mismatch");
  if (!input.commit?.trim()) failures.push("missing exact commit");
  else if (!/^[0-9a-f]{7,64}$/i.test(input.commit)) failures.push("exact commit format is invalid");
  if (!input.artifactSha256?.trim()) failures.push("missing artifact SHA");
  else if (!/^[0-9a-f]{64}$/i.test(input.artifactSha256)) failures.push("artifact SHA format is invalid");
  if (!input.model?.trim()) failures.push("missing model");
  if (!input.inputHashes || Object.keys(input.inputHashes).length === 0) failures.push("missing input hashes");
  if (!input.engineId.trim()) failures.push("missing engine");
  if (!input.engineVersion.trim()) failures.push("missing engine version");
  if (contract && (input.testId !== contract.testId || input.testVersion !== contract.version)) failures.push("test contract version mismatch");
  return failures;
}

async function audit(eventType: string, userId: number, runId: string | undefined, status: string, payload: unknown) {
  const db = await requireDb();
  const auditId = `k_audit_${randomUUID()}`;
  await db.insert(kernelAudits).values({ auditId, runId, eventType, actorUserId: userId, status, payload: JSON.stringify(payload), createdAt: new Date() });
  return auditId;
}

export async function listKernelStatus() {
  const db = await requireDb();
  const [engines, tests, runs, evidence, gates, audits, drifts] = await Promise.all([
    db.select().from(kernelEngines).orderBy(desc(kernelEngines.createdAt)),
    db.select().from(kernelTests).orderBy(desc(kernelTests.createdAt)),
    db.select().from(kernelRuns).orderBy(desc(kernelRuns.createdAt)).limit(20),
    db.select().from(kernelEvidence).orderBy(desc(kernelEvidence.createdAt)).limit(20),
    db.select().from(kernelGates).orderBy(desc(kernelGates.createdAt)).limit(20),
    db.select().from(kernelAudits).orderBy(desc(kernelAudits.createdAt)).limit(20),
    db.select().from(kernelDrifts).orderBy(desc(kernelDrifts.createdAt)).limit(20),
  ]);
  return { engines, tests, runs, evidence, gates, audits, drifts, adapters: ENGINE_ADAPTERS };
}

export async function registerKernelEngine(input: { userId: number; engineId: string; name: string; version: string; adapterKind: string; capabilities: string[] }) {
  const db = await requireDb();
  await db.insert(kernelEngines).values({ engineId: input.engineId, name: input.name, version: input.version, adapterKind: input.adapterKind, health: "UNKNOWN", capabilities: JSON.stringify(input.capabilities), createdAt: new Date() }).onDuplicateKeyUpdate({ set: { version: input.version, name: input.name, adapterKind: input.adapterKind, capabilities: JSON.stringify(input.capabilities) } });
  await audit("ENGINE_REGISTERED", input.userId, undefined, "RECORDED", input);
  return { engineId: input.engineId, health: "UNKNOWN" as const };
}

export async function registerKernelTest(input: { userId: number; testId: string; version: string; purpose: string; preconditions: string[]; inputs: string[]; inputHashes: string[]; engineId: string; engineVersion: string; command: string; checks: string[]; acceptance: string[]; outputs: string[]; evidenceClass: string; gatePolicy: string; protectedArtifact?: string }) {
  const db = await requireDb();
  if (input.protectedArtifact && /(?:V7-R3|R4\.1)/i.test(input.protectedArtifact)) {
    await audit("PROTECTED_MUTATION_DENIED", input.userId, undefined, "BLOCKED", { protectedArtifact: input.protectedArtifact, testId: input.testId });
    throw new Error("protected baseline mutation denied");
  }
  await db.insert(kernelTests).values({ testId: input.testId, version: input.version, purpose: input.purpose, preconditions: JSON.stringify(input.preconditions), inputs: JSON.stringify(input.inputs), inputHashes: JSON.stringify(input.inputHashes), engineId: input.engineId, engineVersion: input.engineVersion, command: input.command, checks: JSON.stringify(input.checks), acceptance: JSON.stringify(input.acceptance), outputs: JSON.stringify(input.outputs), evidenceClass: input.evidenceClass, gatePolicy: input.gatePolicy, protectedArtifact: input.protectedArtifact, createdAt: new Date() }).onDuplicateKeyUpdate({ set: { version: input.version, purpose: input.purpose, command: input.command, engineVersion: input.engineVersion, gatePolicy: input.gatePolicy } });
  await audit("TEST_REGISTERED", input.userId, undefined, "RECORDED", { testId: input.testId, version: input.version });
  return { testId: input.testId, version: input.version };
}

export async function planKernelRun(input: PlanInput) {
  const db = await requireDb();
  const contract = input.testId === R41_TEST_CONTRACT.testId ? R41_TEST_CONTRACT : undefined;
  const failures = validationFailures(input, contract);
  const engine = await db.select().from(kernelEngines).where(and(eq(kernelEngines.engineId, input.engineId), eq(kernelEngines.version, input.engineVersion))).limit(1);
  if (!engine[0]) failures.push("engine adapter is not registered");
  const duplicate = await db.select().from(kernelRuns).where(and(eq(kernelRuns.requestedBy, input.userId), eq(kernelRuns.testId, input.testId), eq(kernelRuns.testVersion, input.testVersion), eq(kernelRuns.artifactSha256, input.artifactSha256 ?? ""), inArray(kernelRuns.state, ["QUEUED", "RUNNING"]))).limit(1);
  if (duplicate[0]) return { runId: duplicate[0].runId, state: duplicate[0].state, blockReason: duplicate[0].blockReason, gateId: null, deduplicated: true as const };
  const runId = `k_run_${randomUUID()}`;
  const state: KernelState = failures.length > 0 ? "BLOCKED" : "QUEUED";
  const blockReason = failures.length > 0 ? `BLOCK BEFORE EXECUTION: ${failures.join("; ")}` : null;
  await db.insert(kernelRuns).values({ runId, requestedBy: input.userId, testId: input.testId, testVersion: input.testVersion, repository: input.repository ?? "", commit: input.commit ?? "", artifactSha256: input.artifactSha256 ?? "", model: input.model ?? "", inputHashes: JSON.stringify(input.inputHashes ?? {}), engineId: input.engineId, engineVersion: input.engineVersion, state, blockReason, evidenceIds: JSON.stringify([]), gate: state === "BLOCKED" ? "BLOCKED" : "PENDING", maxAttempts: Math.max(1, Math.min(input.maxAttempts ?? 1, 3)), resourceClass: "isolated", createdAt: new Date() });
  const gateId = `k_gate_${randomUUID()}`;
  await db.insert(kernelGates).values({ gateId, runId, policy: contract?.gatePolicy ?? "BLOCK_ON_UNTRUSTED_UPSTREAM", decision: state === "BLOCKED" ? "BLOCKED" : "NOT_PROVEN", reason: blockReason ?? "Execution is queued; no solver result exists.", dependencies: JSON.stringify(failures), createdAt: new Date() });
  await audit(state === "BLOCKED" ? "EXECUTION_BLOCKED" : "EXECUTION_QUEUED", input.userId, runId, state, { failures, testId: input.testId, engineId: input.engineId });
  return { runId, state, blockReason, gateId, deduplicated: false as const };
}

export async function completeKernelRun(input: { userId: number; runId: string; outputHash: string; result: unknown; resultStatus: "PASS" | "FAIL" | "BLOCKED" | "NOT_PROVEN"; evidenceIds: string[] }) {
  const db = await requireDb();
  const rows = await db.select().from(kernelRuns).where(and(eq(kernelRuns.runId, input.runId), eq(kernelRuns.requestedBy, input.userId))).limit(1);
  const run = rows[0];
  if (!run) throw new Error("run not found or not owned by authenticated user");
  if (run.state !== "QUEUED" && run.state !== "RUNNING") throw new Error("run is not executable");
  if (!input.outputHash.trim()) throw new Error("output hash is required");
  if (input.resultStatus === "PASS" && input.evidenceIds.length === 0) throw new Error("PASS requires evidence");
  const nextState: KernelState = input.resultStatus === "PASS" ? "COMPLETED" : input.resultStatus === "FAIL" ? "FAILED" : input.resultStatus;
  await db.update(kernelRuns).set({ state: nextState, outputHash: input.outputHash, result: JSON.stringify(input.result), evidenceIds: JSON.stringify(input.evidenceIds), gate: nextState === "COMPLETED" ? "PASS" : nextState, leaseOwner: null, leaseExpiresAt: null, completedAt: new Date() }).where(eq(kernelRuns.runId, input.runId));
  await audit("RESULT_RECORDED", input.userId, input.runId, nextState, { outputHash: input.outputHash, evidenceIds: input.evidenceIds });
  return { runId: input.runId, state: nextState };
}

export async function executeKernelRun(input: { userId: number; runId: string }) {
  const db = await requireDb();
  const rows = await db.select().from(kernelRuns).where(and(eq(kernelRuns.runId, input.runId), eq(kernelRuns.requestedBy, input.userId))).limit(1);
  const run = rows[0];
  if (!run) throw new Error("run not found or not owned by authenticated user");
  if (run.state !== "QUEUED" && run.state !== "RUNNING") throw new Error("only queued or running runs can enter execution");
  const reason = "ENGINE_RUNTIME_NOT_AVAILABLE: adapter registration does not prove solver execution";
  await db.update(kernelRuns).set({ state: "BLOCKED", blockReason: reason, gate: "BLOCKED", leaseOwner: null, leaseExpiresAt: null, completedAt: new Date() }).where(eq(kernelRuns.runId, input.runId));
  await db.update(kernelGates).set({ decision: "BLOCKED", reason }).where(eq(kernelGates.runId, input.runId));
  await audit("EXECUTION_BLOCKED", input.userId, input.runId, "BLOCKED", { reason });
  return { runId: input.runId, state: "BLOCKED" as const, reason };
}

export async function cancelKernelRun(input: { userId: number; runId: string }) {
  const db = await requireDb();
  const rows = await db.select().from(kernelRuns).where(and(eq(kernelRuns.runId, input.runId), eq(kernelRuns.requestedBy, input.userId))).limit(1);
  const run = rows[0];
  if (!run) throw new Error("run not found or not owned by authenticated user");
  if (run.state !== "QUEUED" && run.state !== "RUNNING") throw new Error("terminal run cannot be cancelled");
  const reason = "CANCELLED_BY_OPERATOR";
  await db.update(kernelRuns).set({ state: "BLOCKED", blockReason: reason, gate: "BLOCKED", leaseOwner: null, leaseExpiresAt: null, completedAt: new Date() }).where(eq(kernelRuns.runId, input.runId));
  await audit("EXECUTION_CANCELLED", input.userId, input.runId, "BLOCKED", { reason });
  return { runId: input.runId, state: "BLOCKED" as const, reason };
}

export async function claimKernelRun(input: { userId: number; workerId: string; leaseSeconds?: number }) {
  const db = await requireDb();
  const now = new Date();
  const candidate = await db.select().from(kernelRuns).where(and(eq(kernelRuns.requestedBy, input.userId), eq(kernelRuns.state, "QUEUED"), or(isNull(kernelRuns.leaseExpiresAt), lt(kernelRuns.leaseExpiresAt, now)))).orderBy(kernelRuns.createdAt).limit(1);
  const run = candidate[0];
  if (!run) return { claimed: false as const, reason: "NO_QUEUED_RUN" };
  if (run.attempts >= run.maxAttempts) {
    await db.update(kernelRuns).set({ state: "BLOCKED", blockReason: "MAX_ATTEMPTS_EXCEEDED", gate: "BLOCKED", completedAt: now }).where(eq(kernelRuns.runId, run.runId));
    await audit("QUEUE_BLOCKED_MAX_ATTEMPTS", input.userId, run.runId, "BLOCKED", { attempts: run.attempts, maxAttempts: run.maxAttempts });
    return { claimed: false as const, reason: "MAX_ATTEMPTS_EXCEEDED", runId: run.runId };
  }
  const leaseSeconds = Math.max(1, Math.min(input.leaseSeconds ?? 30, 300));
  const leaseExpiresAt = new Date(now.getTime() + leaseSeconds * 1000);
  await db.update(kernelRuns).set({ state: "RUNNING", attempts: run.attempts + 1, leaseOwner: input.workerId, leaseExpiresAt, resourceClass: "isolated" }).where(and(eq(kernelRuns.runId, run.runId), eq(kernelRuns.state, "QUEUED"), or(isNull(kernelRuns.leaseExpiresAt), lt(kernelRuns.leaseExpiresAt, now))));
  await audit("QUEUE_CLAIMED", input.userId, run.runId, "RUNNING", { workerId: input.workerId, leaseExpiresAt, resourceClass: "isolated" });
  return { claimed: true as const, runId: run.runId, leaseOwner: input.workerId, leaseExpiresAt };
}

export async function workerTick(input: { userId: number; workerId: string; leaseSeconds?: number }) {
  const claim = await claimKernelRun(input);
  if (!claim.claimed) return claim;
  const execution = await executeKernelRun({ userId: input.userId, runId: claim.runId });
  return { ...claim, execution };
}

export async function recoverExpiredRuns(input: { userId: number }) {
  const db = await requireDb();
  const now = new Date();
  const expired = await db.select().from(kernelRuns).where(and(eq(kernelRuns.requestedBy, input.userId), eq(kernelRuns.state, "RUNNING"), lt(kernelRuns.leaseExpiresAt, now)));
  let requeued = 0;
  let blocked = 0;
  for (const run of expired) {
    if (run.attempts < run.maxAttempts) {
      await db.update(kernelRuns).set({ state: "QUEUED", leaseOwner: null, leaseExpiresAt: null, blockReason: "WORKER_LEASE_EXPIRED_REQUEUED" }).where(eq(kernelRuns.runId, run.runId));
      await audit("WORKER_CRASH_RECOVERED", input.userId, run.runId, "QUEUED", { attempts: run.attempts, maxAttempts: run.maxAttempts });
      requeued += 1;
    } else {
      await db.update(kernelRuns).set({ state: "BLOCKED", gate: "BLOCKED", leaseOwner: null, leaseExpiresAt: null, blockReason: "WORKER_LEASE_EXPIRED_MAX_ATTEMPTS", completedAt: now }).where(eq(kernelRuns.runId, run.runId));
      await audit("WORKER_CRASH_BLOCKED", input.userId, run.runId, "BLOCKED", { attempts: run.attempts, maxAttempts: run.maxAttempts });
      blocked += 1;
    }
  }
  return { requeued, blocked };
}

export async function detectKernelGaps() {
  const db = await requireDb();
  const [engines, tests, runs] = await Promise.all([db.select().from(kernelEngines), db.select().from(kernelTests), db.select().from(kernelRuns).orderBy(desc(kernelRuns.createdAt)).limit(50)]);
  const gaps: Array<{ id: string; status: "OPEN" | "BLOCKED" | "NOT_PROVEN"; impact: string; remediation: string }> = [];
  if (engines.length === 0) gaps.push({ id: "K-GAP-ENGINE-REGISTRY", status: "OPEN", impact: "No adapter is registered.", remediation: "Register an adapter and separately verify runtime availability." });
  if (tests.length === 0) gaps.push({ id: "K-GAP-TEST-REGISTRY", status: "OPEN", impact: "No versioned contract is persisted.", remediation: "Register a complete test contract." });
  if (runs.some((run) => run.state === "BLOCKED" || run.state === "NOT_PROVEN")) gaps.push({ id: "K-GAP-BLOCKED-RUNS", status: "BLOCKED", impact: "At least one run is blocked or unverified.", remediation: "Resolve upstream trust/evidence gaps; do not promote downstream results." });
  if (runs.length > 0 && runs.every((run) => run.state !== "COMPLETED")) gaps.push({ id: "K-GAP-NO-VERIFIED-RUN", status: "NOT_PROVEN", impact: "No completed kernel run exists.", remediation: "Provide trusted artifact inputs and a real engine runtime." });
  return gaps;
}

export function compareReproducibility(input: { first: Record<string, unknown>; second: Record<string, unknown> }) {
  const keys = ["repository", "commit", "artifactSha", "model", "testVersion", "inputs", "engineVersion", "outputs", "result"];
  const differences = keys.filter((key) => JSON.stringify(input.first[key]) !== JSON.stringify(input.second[key]));
  return { reproducible: differences.length === 0, differences };
}

export function calculateDependencyImpact(input: { changedDependencies: string[]; tests: Array<{ testId: string; dependencies: string[] }>; runs: Array<{ runId: string; testId: string; state: string }> }) {
  const changed = new Set(input.changedDependencies);
  const affectedTests = input.tests.filter((test) => test.dependencies.some((dependency) => changed.has(dependency))).map((test) => test.testId);
  const affectedSet = new Set(affectedTests);
  const staleRuns = input.runs.filter((run) => affectedSet.has(run.testId)).map((run) => ({ runId: run.runId, priorState: run.state, status: "STALE" as const }));
  return { affectedTests, staleRuns, rerunRequired: affectedTests.length > 0 };
}

export async function getKernelRun(userId: number, runId: string) {
  const db = await requireDb();
  const rows = await db.select().from(kernelRuns).where(and(eq(kernelRuns.runId, runId), eq(kernelRuns.requestedBy, userId))).limit(1);
  if (!rows[0]) return undefined;
  const [evidence, gates, audits] = await Promise.all([
    db.select().from(kernelEvidence).where(eq(kernelEvidence.runId, runId)),
    db.select().from(kernelGates).where(eq(kernelGates.runId, runId)),
    db.select().from(kernelAudits).where(eq(kernelAudits.runId, runId)).orderBy(desc(kernelAudits.createdAt)),
  ]);
  return { run: rows[0], evidence, gates, audits };
}

export async function recordKernelEvidence(input: { userId: number; runId: string; evidenceClass: string; sourceHash: string; payload: unknown; provenance: unknown }) {
  const db = await requireDb();
  const owned = await db.select().from(kernelRuns).where(and(eq(kernelRuns.runId, input.runId), eq(kernelRuns.requestedBy, input.userId))).limit(1);
  if (!owned[0]) throw new Error("run not found or not owned by authenticated user");
  const evidenceId = `k_ev_${randomUUID()}`;
  await db.insert(kernelEvidence).values({ evidenceId, runId: input.runId, evidenceClass: input.evidenceClass, sourceHash: input.sourceHash, payload: JSON.stringify(input.payload), provenance: JSON.stringify(input.provenance), createdAt: new Date() });
  await audit("EVIDENCE_RECORDED", input.userId, input.runId, "RECORDED", { evidenceId, sourceHash: input.sourceHash });
  return { evidenceId };
}

export async function recordKernelDrift(input: { userId: number; runId?: string; artifactSha256: string; baselineSha256: string; impact: string[] }) {
  const db = await requireDb();
  const driftId = `k_drift_${randomUUID()}`;
  const decision = input.artifactSha256 === input.baselineSha256 ? "MATCH" : "DRIFT";
  await db.insert(kernelDrifts).values({ driftId, runId: input.runId, artifactSha256: input.artifactSha256, baselineSha256: input.baselineSha256, decision, impact: JSON.stringify(input.impact), createdAt: new Date() });
  await audit("DRIFT_RECORDED", input.userId, input.runId, decision, { driftId, impact: input.impact });
  return { driftId, decision };
}

export async function deleteKernelRunForTest(userId: number, runId: string) {
  const db = await requireDb();
  await db.delete(kernelEvidence).where(eq(kernelEvidence.runId, runId));
  await db.delete(kernelGates).where(eq(kernelGates.runId, runId));
  await db.delete(kernelAudits).where(eq(kernelAudits.runId, runId));
  await db.delete(kernelRuns).where(and(eq(kernelRuns.runId, runId), eq(kernelRuns.requestedBy, userId)));
}
