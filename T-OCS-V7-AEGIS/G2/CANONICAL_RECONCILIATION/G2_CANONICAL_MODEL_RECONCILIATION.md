# T-OCS V7.0 REVISION-2 — G2 CANONICAL MODEL RECONCILIATION

## CANONICAL MODEL STATUS

**BLOCKED**

## Conflict A — four-joint vs six-joint representation

**SOURCE A:** `artifacts/master_v4/v7_baseline_independent_verification.json`, audit commit `02b7919cddd54a116ba74d6fba13e8983d001231`, `$.kinematic_audit.declared_joints`.

Finding: V7 declares four joint families: `J_RAIL`, `J_SEATBACK`, `J_LINK_L`, `J_LINK_R`. The same source states `executable_joint_graph=false` and DOF/limits are metadata only.

**SOURCE B:** `SR11/joints/V1_SR11_JOINT_MAPPING.json`, V1, blob `10f4c658c77037147f788deeacbfa853122e3636`.

Finding: six records are present: two rail prismatic joints, one seatback revolute, two model-defined linkages, and one state-dependent `LOCK170` constraint.

**Classification:** REPRESENTATION_GRAIN_DIFFERENCE.

**Resolution:** PARTIALLY RESOLVED. The six-record SR11 representation is historical/model-specification decomposition; the four-family representation is the later V7 declaration. No controlled executable crosswalk exists, so no canonical numerical joint inventory may be released.

## Conflict B — 16 vs 31 vs 62

**16** = V7 reference/BOM record grain.

**31** = logical MBD body-group grain in `SR11/manifests/V1_SR11_BODY_MAPPING.json`.

**62** = R4.1 geometric B-Rep solid count.

The SR11 mapping explicitly states deterministic solid-index mapping and covers all 62 solids. The recovery report records 31 grouped body records and 62/62 coverage.

**Classification:** REPRESENTATION_GRAIN_DIFFERENCE.

**Resolution:** RESOLVED AS GRAIN DIFFERENCE. No evidence establishes a single authoritative physical rigid-body count. Therefore no physical-body count is promoted.

## Conflict C — historical R4.1 hash assertion

The current `R4_1_BASELINE_LOCK.json` is internally consistent: the governing SHA matches the artifact, solid count is 62, and mutation is false. An older V4 record reported `BASELINE_HASH_MATCH=false`.

**Classification:** HISTORICAL_RECORD_CONFLICT.

**Resolution:** RESOLVED BY CURRENT AUTHORITY. The old statement remains provenance only; the current immutable lock governs.

## Canonical decision

| Field | Status |
|---|---|
| canonical body set | CONFLICTED |
| canonical joint set | CONFLICTED |
| canonical joint graph | DEFINED_NOT_EXECUTABLE |
| canonical DOF | DEFINED_NOT_EXECUTABLE |
| canonical frames | MISSING |
| canonical axes | PARTIALLY_VERIFIED |
| coordinate convention | VERIFIED |
| canonical constraints | MISSING |
| limits/end-stops | DEFINED_NOT_EXECUTABLE |
| lock/rebound states | DEFINED_NOT_EXECUTABLE |
| canonical initial state | MISSING |
| dynamic fields | NOT_REQUIRED_FOR_DECLARED_SCOPE |

## Release gate

Required authoritative model fields: **NO**.

Representation conflicts resolved for execution: **NO**.

Numerically executable without reinterpretation: **NO**.

Unsupported engineering assumptions introduced: **NO**.

**G2 CANONICAL MODEL = BLOCKED**

## Protection

R4.1 SHA-256 = `fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

Size = 1,220,000 bytes.

Solids = 62.

Mutation = FALSE.

G2 execution = NOT_STARTED.

G2 repeatability = NOT_STARTED.

G3 = LOCKED.
