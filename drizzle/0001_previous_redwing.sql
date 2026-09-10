CREATE TABLE `audit_records` (
	`auditId` varchar(64) NOT NULL,
	`eventType` varchar(128) NOT NULL,
	`actorUserId` int,
	`runId` varchar(64),
	`status` varchar(32) NOT NULL,
	`payload` text NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `audit_records_auditId` PRIMARY KEY(`auditId`)
);
--> statement-breakpoint
CREATE TABLE `qualification_results` (
	`resultId` varchar(64) NOT NULL,
	`runId` varchar(64) NOT NULL,
	`result` text NOT NULL,
	`output` text NOT NULL,
	`outputHash` varchar(128) NOT NULL,
	`resultStatus` enum('PASS','FAIL','BLOCKED','NOT_PROVEN') NOT NULL,
	`evidenceIds` text NOT NULL,
	`auditId` varchar(64) NOT NULL,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `qualification_results_resultId` PRIMARY KEY(`resultId`),
	CONSTRAINT `qualification_results_runId_unique` UNIQUE(`runId`)
);
--> statement-breakpoint
CREATE TABLE `qualification_runs` (
	`runId` varchar(64) NOT NULL,
	`requestedBy` int NOT NULL,
	`executor` varchar(128) NOT NULL,
	`testId` varchar(128) NOT NULL,
	`testVersion` varchar(64) NOT NULL,
	`repository` varchar(512) NOT NULL,
	`branch` varchar(255) NOT NULL,
	`commit` varchar(128) NOT NULL,
	`modelId` varchar(128),
	`modelVersion` varchar(64),
	`modelSha256` varchar(128),
	`inputs` text NOT NULL,
	`engine` varchar(128) NOT NULL,
	`engineVersion` varchar(64) NOT NULL,
	`command` text NOT NULL,
	`state` enum('CREATED','RUNNING','COMPLETED','FAILED','BLOCKED') NOT NULL DEFAULT 'CREATED',
	`startTime` timestamp NOT NULL DEFAULT (now()),
	`endTime` timestamp,
	`createdAt` timestamp NOT NULL DEFAULT (now()),
	CONSTRAINT `qualification_runs_runId` PRIMARY KEY(`runId`)
);
