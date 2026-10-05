# T-OCS V7.0 REVISION-2 — G2 CANONICAL MODEL RECONCILIATION

## CANONICAL MODEL STATUS

**BLOCKED**

This is evidence-first and fail-closed. G2 was **not executed**, CasADi was **not run**, R4.1 was **not changed**, and G3 remains **locked**.

## Conflict register

| Conflict | Source A | Source B | Classification | Resolution status |
|---|---|---|---|---|
| A: 4 vs 6 joints | V7 verification: 4 declared families; graph explicitly not executable | SR11 V1: 6 MODEL_SPECIFICATION_ONLY mixed records | Conflicted representation grain; source B mixes physical joints, logical/model-defined linkages, and LOCK170 constraint | **BLOCKED** — no controlled crosswalk proves containment or a current graph |
| B: 16 vs recorded 32 bodies | V7: 16 BOM/reference records | G2 closure register claims 32 SR11 bodies; direct SR11 mapping and scaffold each contain **31** logical groups | Unsupported source-record count plus representation-grain difference | 32 claim **rejected**; canonical physical body set remains **BLOCKED** |
| C-001: R4.1 identity | Current immutable lock | Historical P0 record says hash mismatch | Historical provenance conflict | **RESOLVED** by current immutable lock |

### Conflict A evidence

The V7 source declares `J_RAIL`, `J_SEATBACK`, `J_LINK_L`, and `J_LINK_R`, but expressly says `executable_joint_graph=false`. SR11 V1 contains two rail prismatics, one seatback revolute, two model-defined linkages, and a state-dependent `LOCK170` constraint. The two rail records may be a bilateral decomposition, but that is not a controlled crosswalk. It must not be silently chosen. `LOCK170` is a constraint/state record, not proof of a sixth physical joint.

### Conflict B evidence

16 is V7 BOM/reference grain. The recorded 32-body assertion fails direct count: `SR11/manifests/V1_SR11_BODY_MAPPING.json` has 31 logical MBD groups and `G2/MODEL/bodies.json` also has 31. SR11 maps those groups to 62 R4.1 solids; the V7 reconstruction has 21 B-Rep solids. No controlled source says whether duplicate, reference, visual, or construction bodies are included in a current physical rigid-body model. The counts therefore do **not** establish the same model revision or a supersession relation.

## Canonical-model decision matrix

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
| required_dynamic_fields | NOT_REQUIRED_FOR_DECLARED_SCOPE |

## Separate registers

- **Register A — T-OCS authoritative model inputs:** `G2_MODEL_INPUT_REGISTER.json`.
- **Register B — G2 execution protocol:** `G2_EXECUTION_PROTOCOL_REGISTER.json`. Protocol controls are explicitly not design facts.

## Authoritative release test

Every required field authoritative: **NO**. Conflicts resolved for execution: **NO**. Numerically expressible without reinterpretation: **NO**. Declared scope executable: **NO**. No values were guessed: **YES**.

**G2 CANONICAL MODEL = BLOCKED**

## Decision required

Engineering authority must resolve DEC-001 through DEC-006 in `G2_DECISION_REQUIRED_REGISTER.json`. This is not optional paperwork: without those decisions, any executable model would be an unsupported synthetic substitute.

## R4.1 protection

- SHA-256: `fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`
- Size: `1,220,000 bytes`
- Solids: `62`
- Mutation: `FALSE`

## Final status

```text
CANONICAL MODEL STATUS
BLOCKED

G2 INPUT CLOSURE
BLOCKED

G2 EXECUTION
NOT_STARTED

G2 REPEATABILITY
NOT_STARTED

G3
LOCKED

R4.1
FAIL

DECISION REQUIRED
LIST

FINAL STATUS
BLOCKED
```
