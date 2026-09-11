import { int, mysqlEnum, mysqlTable, text, timestamp, varchar } from "drizzle-orm/mysql-core";

export const users = mysqlTable("users", {
  id: int("id").autoincrement().primaryKey(),
  openId: varchar("openId", { length: 64 }).notNull().unique(),
  name: text("name"),
  email: varchar("email", { length: 320 }),
  loginMethod: varchar("loginMethod", { length: 64 }),
  role: mysqlEnum("role", ["user", "admin"]).default("user").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
  updatedAt: timestamp("updatedAt").defaultNow().onUpdateNow().notNull(),
  lastSignedIn: timestamp("lastSignedIn").defaultNow().notNull(),
});

export const qualificationRuns = mysqlTable("qualification_runs", {
  runId: varchar("runId", { length: 64 }).primaryKey(),
  requestedBy: int("requestedBy").notNull(),
  executor: varchar("executor", { length: 128 }).notNull(),
  testId: varchar("testId", { length: 128 }).notNull(),
  testVersion: varchar("testVersion", { length: 64 }).notNull(),
  repository: varchar("repository", { length: 512 }).notNull(),
  branch: varchar("branch", { length: 255 }).notNull(),
  commit: varchar("commit", { length: 128 }).notNull(),
  modelId: varchar("modelId", { length: 128 }),
  modelVersion: varchar("modelVersion", { length: 64 }),
  modelSha256: varchar("modelSha256", { length: 128 }),
  inputs: text("inputs").notNull(),
  engine: varchar("engine", { length: 128 }).notNull(),
  engineVersion: varchar("engineVersion", { length: 64 }).notNull(),
  command: text("command").notNull(),
  state: mysqlEnum("state", ["CREATED", "RUNNING", "COMPLETED", "FAILED", "BLOCKED"]).default("CREATED").notNull(),
  startTime: timestamp("startTime").defaultNow().notNull(),
  endTime: timestamp("endTime"),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const qualificationResults = mysqlTable("qualification_results", {
  resultId: varchar("resultId", { length: 64 }).primaryKey(),
  runId: varchar("runId", { length: 64 }).notNull().unique(),
  result: text("result").notNull(),
  output: text("output").notNull(),
  outputHash: varchar("outputHash", { length: 128 }).notNull(),
  resultStatus: mysqlEnum("resultStatus", ["PASS", "FAIL", "BLOCKED", "NOT_PROVEN"]).notNull(),
  evidenceIds: text("evidenceIds").notNull(),
  auditId: varchar("auditId", { length: 64 }).notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const auditRecords = mysqlTable("audit_records", {
  auditId: varchar("auditId", { length: 64 }).primaryKey(),
  eventType: varchar("eventType", { length: 128 }).notNull(),
  actorUserId: int("actorUserId"),
  runId: varchar("runId", { length: 64 }),
  status: varchar("status", { length: 32 }).notNull(),
  payload: text("payload").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const kernelEngines = mysqlTable("kernel_engines", {
  engineId: varchar("engineId", { length: 128 }).primaryKey(),
  name: varchar("name", { length: 128 }).notNull(),
  version: varchar("version", { length: 64 }).notNull(),
  adapterKind: varchar("adapterKind", { length: 64 }).notNull(),
  health: mysqlEnum("health", ["UNKNOWN", "READY", "BLOCKED", "FAILED"]).default("UNKNOWN").notNull(),
  capabilities: text("capabilities").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const kernelTests = mysqlTable("kernel_tests", {
  testId: varchar("testId", { length: 128 }).primaryKey(),
  version: varchar("version", { length: 64 }).notNull(),
  purpose: text("purpose").notNull(),
  preconditions: text("preconditions").notNull(),
  inputs: text("inputs").notNull(),
  inputHashes: text("inputHashes").notNull(),
  engineId: varchar("engineId", { length: 128 }).notNull(),
  engineVersion: varchar("engineVersion", { length: 64 }).notNull(),
  command: text("command").notNull(),
  checks: text("checks").notNull(),
  acceptance: text("acceptance").notNull(),
  outputs: text("outputs").notNull(),
  evidenceClass: varchar("evidenceClass", { length: 64 }).notNull(),
  gatePolicy: text("gatePolicy").notNull(),
  protectedArtifact: varchar("protectedArtifact", { length: 128 }),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const kernelRuns = mysqlTable("kernel_runs", {
  runId: varchar("runId", { length: 64 }).primaryKey(),
  requestedBy: int("requestedBy").notNull(),
  testId: varchar("testId", { length: 128 }).notNull(),
  testVersion: varchar("testVersion", { length: 64 }).notNull(),
  repository: varchar("repository", { length: 512 }).notNull(),
  commit: varchar("commit", { length: 128 }).notNull(),
  artifactSha256: varchar("artifactSha256", { length: 128 }).notNull(),
  model: varchar("model", { length: 256 }).notNull(),
  inputHashes: text("inputHashes").notNull(),
  engineId: varchar("engineId", { length: 128 }).notNull(),
  engineVersion: varchar("engineVersion", { length: 64 }).notNull(),
  state: mysqlEnum("state", ["BLOCKED", "QUEUED", "RUNNING", "COMPLETED", "FAILED", "NOT_PROVEN"]).default("BLOCKED").notNull(),
  blockReason: text("blockReason"),
  outputHash: varchar("outputHash", { length: 128 }),
  result: text("result"),
  evidenceIds: text("evidenceIds").notNull(),
  gate: varchar("gate", { length: 32 }).notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
  completedAt: timestamp("completedAt"),
});

export const kernelEvidence = mysqlTable("kernel_evidence", {
  evidenceId: varchar("evidenceId", { length: 64 }).primaryKey(),
  runId: varchar("runId", { length: 64 }).notNull(),
  evidenceClass: varchar("evidenceClass", { length: 64 }).notNull(),
  sourceHash: varchar("sourceHash", { length: 128 }).notNull(),
  payload: text("payload").notNull(),
  provenance: text("provenance").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const kernelAudits = mysqlTable("kernel_audits", {
  auditId: varchar("auditId", { length: 64 }).primaryKey(),
  runId: varchar("runId", { length: 64 }),
  eventType: varchar("eventType", { length: 128 }).notNull(),
  actorUserId: int("actorUserId"),
  status: varchar("status", { length: 32 }).notNull(),
  payload: text("payload").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const kernelGates = mysqlTable("kernel_gates", {
  gateId: varchar("gateId", { length: 64 }).primaryKey(),
  runId: varchar("runId", { length: 64 }).notNull(),
  policy: varchar("policy", { length: 128 }).notNull(),
  decision: mysqlEnum("decision", ["BLOCKED", "PASS", "FAIL", "NOT_PROVEN"]).notNull(),
  reason: text("reason").notNull(),
  dependencies: text("dependencies").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export const kernelDrifts = mysqlTable("kernel_drifts", {
  driftId: varchar("driftId", { length: 64 }).primaryKey(),
  runId: varchar("runId", { length: 64 }),
  artifactSha256: varchar("artifactSha256", { length: 128 }).notNull(),
  baselineSha256: varchar("baselineSha256", { length: 128 }).notNull(),
  decision: mysqlEnum("decision", ["UNKNOWN", "MATCH", "DRIFT"]).notNull(),
  impact: text("impact").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
});

export type User = typeof users.$inferSelect;
export type InsertUser = typeof users.$inferInsert;
export type QualificationRun = typeof qualificationRuns.$inferSelect;
export type QualificationResult = typeof qualificationResults.$inferSelect;
export type AuditRecord = typeof auditRecords.$inferSelect;
export type KernelEngine = typeof kernelEngines.$inferSelect;
export type KernelTest = typeof kernelTests.$inferSelect;
export type KernelRun = typeof kernelRuns.$inferSelect;
export type KernelEvidence = typeof kernelEvidence.$inferSelect;
export type KernelAudit = typeof kernelAudits.$inferSelect;
export type KernelGate = typeof kernelGates.$inferSelect;
export type KernelDrift = typeof kernelDrifts.$inferSelect;
