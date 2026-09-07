# Engineering Decision Map

## Lock-170 candidate selection

**Decision:** rank L170-A/B/C but do not authorize final selection. **Required parameters:** latch dynamics, acceleration, response time, inertia, contact, friction, compliance, actuator response, loads, and energy. **Available evidence:** candidate definitions only; no released V7 mechanism or calibrated physical response. **Derivation/virtual path:** EVA parameterization, Chrono dynamics, uncertainty, sweep, and sensitivity are available as interfaces. **Minimum test:** one calibrated lock timing/load characterization linked to raw data. **Current status:** BLOCKED; `FINAL_SELECTION = NOT_AUTHORIZED`.

## Rebound-180

**Decision:** preserve 180 mm as a target. **Required parameters:** displacement, velocity, acceleration, damping, contact, stop force, and energy dissipation. **Available evidence:** target only. **Derivation/virtual path:** parameterized rebound envelope and sensitivity are available. **Minimum test:** calibrated displacement/velocity/force acquisition for a representative subsystem. **Current status:** TEST_REQUIRED/BLOCKED.

## Materials and fasteners

**Decision:** do not release a manufacturing or FE card from appearance. **Required parameters:** material grade, treatment, density, stiffness, strength, rate dependence, fastener identity, torque, preload, grip, and joint stiffness. **Available evidence:** no released V7 certificates/BOM. **Derivation/virtual path:** mass, stiffness, preload, and bounds can be calculated after authoritative inputs arrive. **Minimum input:** controlled material and fastener records; physical characterization only where sensitivity shows it is decision-critical. **Current status:** UNVERIFIED/BLOCKED.

## Physical test strategy

**Decision:** minimize tests by information gain. **Required parameters:** affected gaps, uncertainty reduction, sensitivity, safety importance, effort, and dependency unlocks. **Available evidence:** candidate priority matrix, not measured results. **Minimum test set:** vehicle pulse/acceleration, absorber force-stroke, joint/lock characterization, and certificate review, subject to test-authority approval. **Current status:** CANDIDATE_TEST_PLAN.
