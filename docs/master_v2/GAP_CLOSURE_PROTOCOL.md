# Gap Closure Protocol

Allowed states are `OPEN`, `IN_PROGRESS`, `CANDIDATE`, `DERIVED`, `VERIFIED`, `VALIDATED`, `BLOCKED`, `NOT_PRESENT`, and `NOT_APPLICABLE`. `INSTALLED`, `TEST_CREATED`, and `SIMULATION_COMPLETED` are not closure states.

A gap can become `VERIFIED` only after deterministic criteria pass. It can become `VALIDATED` only after the required engineering or physical evidence passes. Every gap has a blocker, evidence requirement, verification method, next action, owner class, dependencies, closure criteria, and last-verified timestamp.
