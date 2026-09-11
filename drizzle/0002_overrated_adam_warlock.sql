CREATE TABLE `kernel_audits` (
	`auditId` varchar(64) NOT NULL,
	`runId` varchar(64),
	`eventType` varchar(128) NOT NULL,
	`actorUserId` int,
	`status` varchar(32) NOT NULL,
	`payload` text NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `kernel_audits_auditId` PRIMARY KEY(`auditId`)
);
--> statement-breakpoint
CREATE TABLE `kernel_drifts` (
	`driftId` varchar(64) NOT NULL,
	`runId` varchar(64),
	`artifactSha256` varchar(128) NOT NULL,
	`baselineSha256` varchar(128) NOT NULL,
	`decision` enum('UNKNOWN','MATCH','DRIFT') NOT NULL,
	`impact` text NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `kernel_drifts_driftId` PRIMARY KEY(`driftId`)
);
--> statement-breakpoint
CREATE TABLE `kernel_engines` (
	`engineId` varchar(128) NOT NULL,
	`name` varchar(128) NOT NULL,
	`version` varchar(64) NOT NULL,
	`adapterKind` varchar(64) NOT NULL,
	`health` enum('UNKNOWN','READY','BLOCKED','FAILED') NOT NULL DEFAULT 'UNKNOWN',
	`capabilities` text NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `kernel_engines_engineId` PRIMARY KEY(`engineId`)
);
--> statement-breakpoint
CREATE TABLE `kernel_evidence` (
	`evidenceId` varchar(64) NOT NULL,
	`runId` varchar(64) NOT NULL,
	`evidenceClass` varchar(64) NOT NULL,
	`sourceHash` varchar(128) NOT NULL,
	`payload` text NOT NULL,
	`provenance` text NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `kernel_evidence_evidenceId` PRIMARY KEY(`evidenceId`)
);
--> statement-breakpoint
CREATE TABLE `kernel_gates` (
	`gateId` varchar(64) NOT NULL,
	`runId` varchar(64) NOT NULL,
	`policy` varchar(128) NOT NULL,
	`decision` enum('BLOCKED','PASS','FAIL','NOT_PROVEN') NOT NULL,
	`reason` text NOT NULL,
	`dependencies` text NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `kernel_gates_gateId` PRIMARY KEY(`gateId`)
);
--> statement-breakpoint
CREATE TABLE `kernel_runs` (
	`runId` varchar(64) NOT NULL,
	`requestedBy` int NOT NULL,
	`testId` varchar(128) NOT NULL,
	`testVersion` varchar(64) NOT NULL,
	`repository` varchar(512) NOT NULL,
	`commit` varchar(128) NOT NULL,
	`artifactSha256` varchar(128) NOT NULL,
	`model` varchar(256) NOT NULL,
	`inputHashes` text NOT NULL,
	`engineId` varchar(128) NOT NULL,
	`engineVersion` varchar(64) NOT NULL,
	`state` enum('BLOCKED','QUEUED','RUNNING','COMPLETED','FAILED','NOT_PROVEN') NOT NULL DEFAULT 'BLOCKED',
	`blockReason` text,
	`outputHash` varchar(128),
	`result` text,
	`evidenceIds` text NOT NULL,
	`gate` varchar(32) NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	`completedAt` timestamp,
	CONSTRAINT `kernel_runs_runId` PRIMARY KEY(`runId`)
);
--> statement-breakpoint
CREATE TABLE `kernel_tests` (
	`testId` varchar(128) NOT NULL,
	`version` varchar(64) NOT NULL,
	`purpose` text NOT NULL,
	`preconditions` text NOT NULL,
	`inputs` text NOT NULL,
	`inputHashes` text NOT NULL,
	`engineId` varchar(128) NOT NULL,
	`engineVersion` varchar(64) NOT NULL,
	`command` text NOT NULL,
	`checks` text NOT NULL,
	`acceptance` text NOT NULL,
	`outputs` text NOT NULL,
	`evidenceClass` varchar(64) NOT NULL,
	`gatePolicy` text NOT NULL,
	`protectedArtifact` varchar(128),
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `kernel_tests_testId` PRIMARY KEY(`testId`)
);
