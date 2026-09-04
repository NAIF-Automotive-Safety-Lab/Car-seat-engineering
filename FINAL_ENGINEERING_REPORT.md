# T-OCS R4.1 Evidence Closure — Final Engineering Report

## BASELINE STATUS
R4.1 = IMMUTABLE
SHA-256 = fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68
SOLIDS = 62
STEP_READ = PASS
B-REP_VALIDITY = PASS
MUTATION = FALSE

## EVIDENCE STATUS
ARTIFACT_RECOVERY = PARTIAL
EVIDENCE_INTEGRITY = PASS

## INTERSECTION STATUS
Historical R4 scope:
- R1↔R1 = 69
- R1↔R4 = 88
- R4↔R4 = 27
- R4 TOTAL = 184

R4.1:
- R1↔R4.1 = 90
- R4↔R4.1 = 30
- R4.1 TOTAL = 189

Independent raw B-rep audit:
- positive-volume overlaps = 189
- cross-body = 170
- intra-body = 19

The raw 189-overlap audit is NOT accepted as 189 defects and is NOT interchangeable with the historical scope. Full semantic ownership for all raw overlaps remains incomplete.

The historical 27-target scope has resolved=27 and unresolved=0 in the existing evidence record.

## PHYSICAL-PARAMETER STATUS
Mass / CG / inertia / density / friction / contact law / absorber law / joint compliance data remain UNKNOWN unless source-verified.

## MECHANISM STATUS
R4.1 model DOF = 2 (MODEL EXECUTED).
Hardware verification = NOT DONE.

## LOAD PATH
PROVEN_TOPOLOGY:
190_Shoulder_Anchor
→ R4_190_SHOULDER_GUSSET
→ SHOULDER_ROOT_MEMBER
→ OUTBOARD_BRACE
→ FRAME_OUTBOARD_BRACKET
→ 160_Seatback_Structural_Frame

LOAD_CAPACITY = UNKNOWN.

## MBD STATUS
REAL_MBD = BLOCKED
PYCHRONO = BLOCKED
REACTIONS = NOT_AVAILABLE
FL/FR = NOT_AVAILABLE
CONVERGENCE = BLOCKED

No fake solver results were produced.

## DYNAMIC BILATERAL STATUS
Geometry symmetry = PASS for 19 audited L/R pairs.
Dynamic synchronization = NOT_EXECUTED.
Rack/binding = UNKNOWN.

## ROOT-CAUSE STATUS
No confirmed mechanical design failure was proven.
Current blockers are evidence/toolchain gaps.

## CAD-REPAIR AUTHORIZATION
CAD_REPAIR = NOT_AUTHORIZED
R4.2 = NOT_CREATED

Creating repair geometry now would be unjustified.

## PROMOTION DECISION
PROMOTION = BLOCKED

### VERIFIED FACTS
- R4.1 SHA matches.
- 62 solids read.
- B-Rep validity passes.
- R4.1 remains unchanged.
- R5.1 shoulder topology evidence exists.
- 27-target R4 overlap ownership record is resolved.
- SR11 body mapping is 62/62.
- No silent defaults are present in the current assumptions register.

### UNKNOWN FACTS
- Full semantic ownership of the raw 189-overlap audit.
- Released mass/CG/inertia/material/contact/friction properties.
- Dynamic bilateral response.
- Real MBD reaction histories.
- Dynamic end-stop reactions.
- Dynamic Lock-170 behavior.
- Physical absorber behavior.

### BLOCKERS
1. REAL_MBD runtime unavailable.
2. Physical-property provenance incomplete.
3. Full intersection semantic ownership incomplete.

### NEXT REQUIRED ENGINEERING ACTIONS
1. Provide a real MBD-capable host/runtime and execute the frozen R4.1 model.
2. Provide released physical-property/joint/contact/absorber data.
3. Complete mechanical semantic ownership for the full raw overlap set without changing CAD.
4. Only after verified failure, authorize a minimal repair candidate.

## FINAL DECISION
NO EVIDENCE → NO CLAIM
NO ROOT CAUSE → NO REPAIR
NO REAL SOLVER → NO MBD PASS
NO VALIDATION → NO PROMOTION
