ALTER TABLE `kernel_runs` ADD `attempts` int DEFAULT 0 NOT NULL;--> statement-breakpoint
ALTER TABLE `kernel_runs` ADD `maxAttempts` int DEFAULT 1 NOT NULL;--> statement-breakpoint
ALTER TABLE `kernel_runs` ADD `leaseOwner` varchar(128);--> statement-breakpoint
ALTER TABLE `kernel_runs` ADD `leaseExpiresAt` timestamp;--> statement-breakpoint
ALTER TABLE `kernel_runs` ADD `resourceClass` varchar(64) DEFAULT 'isolated' NOT NULL;