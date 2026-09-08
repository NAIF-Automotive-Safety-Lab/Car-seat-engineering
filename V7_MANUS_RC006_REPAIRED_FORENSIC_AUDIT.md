# V7 Manus RC-006 Repaired-Package Independent Audit

**Mode:** Strict read-only / independent / zero-bypass  
**Final decision:** **BLOCKED**  
**No repair performed:** This audit did not modify the package or any baseline.

## Package integrity

| Field | Result |
|---|---|
| Package | `RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_REPAIRED.zip` |
| Expected package SHA-256 | `1f7fe60e05982faab269a9d62d374147e4916619c1aa202cc37ac0428327f705` |
| Actual package SHA-256 | `1f7fe60e05982faab269a9d62d374147e4916619c1aa202cc37ac0428327f705` |
| Expected external manifest SHA-256 | `a7f25ac8d07e68ee03f5882245a6e701e078ea56926c5a0ac505d841b0bdc365` |
| Actual external manifest SHA-256 | `a7f25ac8d07e68ee03f5882245a6e701e078ea56926c5a0ac505d841b0bdc365` |
| ZIP member count | 28 |
| External manifest count | 28 |
| Exact paths / sizes / SHA | PASS — 0 external-manifest mismatches |
| ZIP integrity | PASS |
| Missing / extra / duplicate ZIP artifacts | 0 / 0 / 0 |

The external repaired manifest matches the ZIP exactly and the expected package/manifest hashes match. However, the ZIP contains a stale embedded `12_MANIFESTS/RC006_CANONICAL_MANIFEST.json`, which is itself named canonical but declares the pre-repair package state.

## Blocking package/provenance defect

| Embedded manifest field | Declared | Required/current | Result |
|---|---|---|---|
| package identity | `RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE` | `RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_REPAIRED` | FAIL |
| revision | `V7-RC-006` | `RC-006-TRACEABILITY-001` | FAIL |
| member count | 27 | 28 | FAIL |
| canonical ZIP SHA | `366519c3226a4c3793f890e8c012d676cc139eacae9fb7a44784def501cf8682` | `1f7fe60e05982faab269a9d62d374147e4916619c1aa202cc37ac0428327f705` | FAIL |
| coverage | 27 declared paths | 28 actual paths | FAIL |

**Decision impact:** Under zero-bypass, the stale embedded canonical manifest is a package/provenance defect. Therefore the final verdict is `BLOCKED`, even though the external repaired manifest and vehicle-ID repair independently pass.

## Master evidence register

| Check | Result |
|---|---|
| Total records | 331 / 331 |
| Unique IDs | 331 / 331 |
| Duplicate IDs | 0 |
| Required fields | PASS |

## VEH traceability repair

The eight repaired IDs are present and distinct: `VEH-DATUM`, `VEH-SEAT-HARDPOINTS`, `VEH-RESTRAINT-ANCHORS`, `VEH-H-POINT`, `VEH-CLEARANCE-ENVELOPE`, `VEH-MOUNTING-LOADS`, `VEH-CRASH-PULSE`, and `VEH-INITIAL-CONDITIONS`. Each has `CURRENT_VALUE = null`, `STATUS = EXTERNAL_SOURCE_REQUIRED`, `VALUE_CLASS = EXTERNAL_REQUIRED`, and `BLOCKING = true`.

Before/after comparison confirms that engineering values, source status, blocking status, and next action were preserved. Only ID, semantic parameter identity, and traceability metadata changed.

## Zero-bypass

No unauthorized promotion was found: no null became numeric, no model-derived value became physical, no OEM input became known, no test-required item became tested, and no reference became validated. The `180 mm` absorber value remains reference-only; no `18–22 kN` validated T-OCS result exists.

## Engineering immutability

- CAD/STEP/geometry changed: **false**  
- Design intent changed: **false**  
- Engineering values changed: **false**  
- R4.1 changed: **false**  
- R4.2 created: **false**  
- CAE executed: **false**  
- Physical tests executed: **false**  
- Manufacturing release: **false**  

## Remaining engineering blockers

Material/density evidence, mass/CG/inertia, joint/fastener evidence, absorber characterization, lock characterization, rebound characterization, vehicle/OEM inputs, CAE physical inputs, physical measurements and physical test results remain open.

## Final verdict

**BLOCKED.** The repair itself is traceability-valid, but the stale embedded canonical manifest prevents reliable zero-bypass package/provenance verification. No repair was applied by this audit.
