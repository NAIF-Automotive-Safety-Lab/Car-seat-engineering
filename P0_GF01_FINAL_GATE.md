# P0-GF01 FINAL OWNER DECISION GATE

## Owner decision

| Field | Status |
|---|---|
| Decision owner | `Naif` |
| Owner decision | `CONDITIONAL_FUNCTIONAL_DEVELOPMENT` |
| Authorization | `FUNCTIONAL_DEVELOPMENT_AUTHORIZED_WITH_GATES` |
| Physical build | `AUTHORIZED_FOR_FUNCTIONAL_DEVELOPMENT` |
| Functional P0 | `DEVELOPMENT_ONLY` |
| R4.1 | `FROZEN` |
| R4.1 SHA256 | `fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68` |
| Physical-input closure | `0/14` |

**Owner command record:** Explicit conditional functional-development command received 2026-09-07. This record does not constitute safety approval, qualification, certification, homologation, or unrestricted FE authorization.

## Authorized scope

`P0-GF01` may proceed only as controlled engineering-development work to generate physical evidence and characterize the design. Developmental results must not be labeled qualification, certification, homologation, or safety validation.

## Mandatory gates

1. Assign and preserve a unique configuration ID and serial.
2. Record actual fabrication process, actual materials, and every deviation from R4.1.
3. Establish temporary/local datums where released datums do not exist.
4. Do not invent tolerances, loads, friction, stiffness, damping, absorber curves, bolt preload, torque, or vehicle hardpoints.
5. Bind every measured value to source, method, instrument, calibration status, raw data, units, uncertainty/tolerance, timestamp, operator/lab, artifact path, SHA256, and R4.1 mapping.
6. Keep every one of the 14 inputs `OPEN` until the secure validator accepts authoritative evidence.
7. B02 mirror relation remains `UNPROVEN` and is not an acceptance criterion.
8. R4.1 remains immutable; R4.2 remains unauthorized.

## Test permissions

| Activity | Status |
|---|---|
| Geometry/fit development | Development only |
| Assembly/alignment development | Development only |
| Static function | Development only |
| Joint/lock function | Development only |
| Absorber characterization | Development only; not executed |
| Friction characterization | Development only; not executed |
| Contact/damping characterization | Development only; not executed |
| Structural development test | Development only; not executed |
| Material characterization | Development only; not executed |
| Crash test | Development only, **blocked until external qualified-facility gate** |
| Safety validation | Not authorized |
| Formal qualification | Not authorized |
| Certification/homologation | Not authorized |
| Vehicle safety release | Not authorized |
| FE authorization | Not authorized until 14/14 input gate closure |

## Explicit crash hard gate

No crash test may begin merely because functional development is authorized. Before any crash test, provide and independently record:

- approved test plan;
- exact test-article configuration;
- instrumentation and DAQ;
- qualified facility capability;
- applicable test protocol;
- safety and emergency controls;
- raw-data retention;
- independent test record;
- external test authority approval.

If any item is missing: `CRASH_TEST = BLOCKED`.

## Current gate status

```text
FUNCTIONAL_TEST_STATUS = AUTHORIZED_NOT_EXECUTED
CRASH_TEST_STATUS = BLOCKED_BY_EXTERNAL_FACILITY_GATE
STRUCTURAL_TEST_STATUS = DEVELOPMENT_ONLY_NOT_EXECUTED
ABSORBER_TEST_STATUS = DEVELOPMENT_ONLY_NOT_EXECUTED
14_INPUT_CLOSURE = 0/14 BLOCKED
FE_AUTHORIZATION = NOT_AUTHORIZED_UNTIL_INPUT_GATE_CLOSED
SAFETY_VALIDATION = NOT_AUTHORIZED
```
