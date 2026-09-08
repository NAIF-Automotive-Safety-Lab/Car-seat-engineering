# V7 Interface Semantics Closure

**Revision:** AUDIT-2026-09-09  
**Source:** V7-RC-004 / V7-RC-005 inputs  
**Mode:** Read-only, evidence-first, zero-bypass

## Final status

| Metric | Result |
|---|---:|
| CLOSED | 0 |
| PARTIAL | 2 |
| UNKNOWN | 7 |
| GEOMETRIC_ERROR | 0 |
| UNSUPPORTED_CLAIMS | 0/FAIL |
| V7_GEOMETRY_CHANGED | NO |
| PHYSICAL_VALIDATION | NOT_DONE |
| Final state | MODEL/ENGINEERING DEFINITION ONLY |

No CAE, dynamics, FE, crash-pulse, absorber simulation, friction calibration, or occupant simulation was started. No baseline geometry was modified.

## Interface register

| ID | PART-A | PART-B | OVERLAP mm³ | CLASS | STATUS | EVIDENCE | JOINT | DOF | CONTACT | CAPTURE | STRUCTURAL ROLE | UNKNOWN REASON |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|
| IF-R4-01 | 120 | 130L | 30250.000000 | UNKNOWN | UNKNOWN | R4 register + geometry evidence + forensic evidence | NONE | NOT_PROVEN | UNKNOWN | UNKNOWN | UNKNOWN | Geometry proves a positive-volume intersection, but the authoritative interface register and forensic record both state mechanical purpose, relative motion, axis datum, DOF, contact rule, clearance rule, load-path role, and fastening are NOT_PROVEN. |
| IF-R4-02 | 120 | 130R | 30250.000000 | UNKNOWN | UNKNOWN | R4 register + geometry evidence + forensic evidence | NONE | NOT_PROVEN | UNKNOWN | UNKNOWN | UNKNOWN | Geometry proves a positive-volume intersection, but the authoritative interface register and forensic record both state mechanical purpose, relative motion, axis datum, DOF, contact rule, clearance rule, load-path role, and fastening are NOT_PROVEN. |
| IF-R4-03 | 120 | 140 | 113575.000000 | UNKNOWN | UNKNOWN | R4 register + geometry evidence + forensic evidence | NONE | NOT_PROVEN | UNKNOWN | UNKNOWN | UNKNOWN | Geometry proves a positive-volume intersection, but the authoritative interface register and forensic record both state mechanical purpose, relative motion, axis datum, DOF, contact rule, clearance rule, load-path role, and fastening are NOT_PROVEN. |
| IF-R4-04 | 120 | 200 | 7671.875043 | UNKNOWN | UNKNOWN | R4 register + geometry evidence + forensic evidence | NONE | NOT_PROVEN | UNKNOWN | UNKNOWN | UNKNOWN | Geometry proves a positive-volume intersection, but the authoritative interface register and forensic record both state mechanical purpose, relative motion, axis datum, DOF, contact rule, clearance rule, load-path role, and fastening are NOT_PROVEN. |
| IF-R4-05 | 130L | 140 | 303345.231682 | UNKNOWN | UNKNOWN | R4 register + geometry evidence + forensic evidence | NONE | NOT_PROVEN | UNKNOWN | UNKNOWN | UNKNOWN | Geometry proves a positive-volume intersection, but the authoritative interface register and forensic record both state mechanical purpose, relative motion, axis datum, DOF, contact rule, clearance rule, load-path role, and fastening are NOT_PROVEN. |
| IF-R4-06 | 130R | 140 | 303345.231682 | UNKNOWN | UNKNOWN | R4 register + geometry evidence + forensic evidence | NONE | NOT_PROVEN | UNKNOWN | UNKNOWN | UNKNOWN | Geometry proves a positive-volume intersection, but the authoritative interface register and forensic record both state mechanical purpose, relative motion, axis datum, DOF, contact rule, clearance rule, load-path role, and fastening are NOT_PROVEN. |
| IF-R4-07 | 140 | 160 | 11200.000000 | UNKNOWN | UNKNOWN | R4 register + geometry evidence + forensic evidence | NONE | NOT_PROVEN | UNKNOWN | UNKNOWN | UNKNOWN | Geometry proves a positive-volume intersection, but the authoritative interface register and forensic record both state mechanical purpose, relative motion, axis datum, DOF, contact rule, clearance rule, load-path role, and fastening are NOT_PROVEN. |
| IF-R4-08 | 150L | 160 | 106752.606055 | JOINT | PARTIAL | R4 register + geometry evidence + forensic evidence | J_LINK_L | q2-coupled | UNKNOWN | UNKNOWN | UNKNOWN | Direct pair matches the model-defined J_LINK_L parent/child relation 150L -> 160; the source marks this as MODEL_DEFINED and physical validation is not done. |
| IF-R4-09 | 160 | J | 67326.535731 | JOINT | PARTIAL | R4 register + geometry evidence + forensic evidence | J_SEATBACK | Ry | UNKNOWN | UNKNOWN | UNKNOWN | Direct pair matches the model-defined J_SEATBACK parent/child relation J -> 160; the source marks this as MODEL_DEFINED and physical validation is not done. |

## Classification decisions

- **IF-R4-08** and **IF-R4-09** are `JOINT / PARTIAL` only because the model explicitly maps their body pairs to `J_LINK_L` and `J_SEATBACK`. This is model-defined evidence, not physical validation or release evidence.
- The other seven interfaces remain `UNKNOWN`; their source records explicitly mark mechanical purpose, relative motion, axis datum, DOF, contact, clearance, load path, and fastening as `NOT_PROVEN`.
- No interface is `CLOSED`.
- No interface is classified `GEOMETRIC_ERROR`; positive volume alone does not prove geometry is wrong.
- No unsupported claims were introduced.

## Evidence hashes

| Path | SHA-256 | Bytes |
|---|---|---:|
| `/home/ubuntu/interface-audit-work/base-r4/INTERFACES/V7_R4_INTERFACE_REGISTER.json` | `e6be7d2326f73e4e9d51f1f2f766e59f66f2edc428dc3748e2b907a1dfc18502` | 7773 |
| `/home/ubuntu/interface-audit-work/base-r4/INTERFACES/V7_R4_INTERFACE_GEOMETRY_EVIDENCE.json` | `ad47a76f28c1c6fcb1b61b3fd100c993bae368ac5a30312bec80c966c9a372c2` | 11454 |
| `/home/ubuntu/interface-audit-work/base-r4/INTERFACES/V7_R4_INTERFACE_FORENSICS.json` | `089e5b2533cd48528e45bd02fb5e84769e625b28687865881b6f6196ad823492` | 11678 |
| `/home/ubuntu/interface-audit-work/r4/V7_R4_KINEMATIC_MODEL.json` | `a02cc7a48eabb9e61fa71b697098bd977362758163389f6f4e986ead43441b73` | 7116 |
| `/home/ubuntu/upload/V7_ENGINEERING_BASELINE_PACKAGE_R4.zip` | `68620362150fd604dfdc3089e12ce346bb5acb1fbb40b02daaa6e60880d08d61` | 681538 |
| `/home/ubuntu/upload/V7_R4_KINEMATIC_CLOSURE_CORRECTED.zip` | `4b0118697b42c15f7b476890bea60f97e7fd16878ccd6f41ab7cb73878fbf274` | 10555 |

## What was verified

The nine interface IDs and body pairs were cross-checked against the uploaded V7-R4 interface register and geometry-forensics records. All nine have positive-volume intersection evidence and face-ID evidence. The sampled motion audit reports 361 samples and nine nominal positive intersections, while explicitly preserving interface semantics as pending and disclaiming continuous CCD.

## What remains open

Physical datums, contact normals, clearance requirements, fastener/interface evidence, load-transfer intent, capture/retention behavior, inspection evidence, physical validation, and production release remain open.
