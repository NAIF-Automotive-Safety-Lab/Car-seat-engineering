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

export type User = typeof users.$inferSelect;
export type InsertUser = typeof users.$inferInsert;
export type QualificationRun = typeof qualificationRuns.$inferSelect;
export type QualificationResult = typeof qualificationResults.$inferSelect;
export type AuditRecord = typeof auditRecords.$inferSelect;
