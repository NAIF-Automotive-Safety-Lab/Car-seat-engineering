# T-OCS V7.0 REVISION-2 — G2 CANONICAL MODEL RECONCILIATION

## Decision
**G2 CANONICAL MODEL = BLOCKED**

### Conflict A — 4-joint vs 6-joint
- V7 audit: `J_RAIL`, `J_SEATBACK`, `J_LINK_L`, `J_LINK_R`; executable graph = false.
- SR11 V1 mapping: two rail prismatics + seatback revolute + two linkages + LOCK170 state constraint.
- Classification: representation-grain difference.
- Resolution: partially resolved. SR11 is historical/model-specification decomposition; V7 is the later declared joint-family representation. No authoritative executable crosswalk exists.

### Conflict B — 16 vs 31 vs 62
- 16 = V7 reference/BOM record grain.
- 31 = SR11 logical MBD body grouping grain.
- 62 = R4.1 B-Rep solid count.
- These are not the same quantity. No evidence establishes a single physical rigid-body count.

### Conflict C — historical baseline hash assertion
- Current R4.1 lock is internally consistent and matches the governing digest.
- The older V4 assertion of `BASELINE_HASH_MATCH=false` is retained only as historical provenance.
- Current authority resolves the conflict without changing the lock.

## Canonical model decision
| Field | Status |
|---|---|
| canonical_body_set | CONFLICTED |
| canonical_joint_set | CONFLICTED |
| canonical_joint_graph | DEFINED_NOT_EXECUTABLE |
| canonical_DOF | DEFINED_NOT_EXECUTABLE |
| canonical_frames | MISSING |
| canonical_axes | PARTIALLY_VERIFIED |
| canonical_coordinate_convention | VERIFIED |
| canonical_constraints | MISSING |
| canonical_limits_end_stops | DEFINED_NOT_EXECUTABLE |
| canonical_lock_rebound_states | DEFINED_NOT_EXECUTABLE |
| canonical_initial_state | MISSING |
| canonical_dynamic_fields | NOT_REQUIRED_FOR_DECLARED_SCOPE |

## Release gate
- All required model fields authoritative: **NO**
- All representation conflicts resolved for execution: **NO**
- Numerically executable without reinterpretation: **NO**
- Unsupported assumptions introduced: **NO**
- **RESULT = BLOCKED**

## R4.1
`fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68` / 1,220,000 bytes / 62 solids / mutation FALSE.

## G2/G3
G2 execution = NOT_STARTED. G2 repeatability = NOT_STARTED. G3 = LOCKED.
