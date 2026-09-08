# MANUS Real Evidence Forensic Audit

**Mode:** Independent, read-only, zero-bypass audit against JON outputs and project baseline  
**Final verdict:** **EVIDENCE_GATE_FAILED**  
**CAE readiness:** **CAE_BLOCKED**  

## Summary

| Metric | Count / result |
|---|---:|
| `TOTAL_CLAIMS_AUDITED` | 14 |
| `SUPPORTED` | 5 |
| `UNSUPPORTED` | 0 |
| `PARTIALLY_SUPPORTED` | 2 |
| `MODEL_DERIVED` | 3 |
| `SOURCE_REQUIRED` | 9 |
| `TEST_REQUIRED` | 6 |
| `OEM_INPUT_REQUIRED` | 2 |
| `ZERO_BYPASS_VIOLATIONS` | 0 |
| `BASELINE_MODIFICATIONS` | 0 |

## Independent conclusion

The audit found no unauthorized promotion in the inspected controlled records: null values remain null, model-derived values are not presented as physical results, the 180 mm value remains design-reference-only, and no validated 18–22 kN T-OCS result is present. The evidence gate nevertheless fails because the required physical/OEM evidence needed to close the engineering domains and CAE gate is absent.

## Claim audit matrix

| ID | Domain | Classification | Result | Promotion check | Evidence basis |
|---|---|---|---|---|---|
| `MAT-01` | materials | `SOURCE_REQUIRED` | `SOURCE_REQUIRED` | `PASS_NO_PROMOTION` | All 16 records have MATERIAL_FAMILY and GRADE null; no authoritative mapping recovered. |
| `MAT-02` | materials | `SOURCE_REQUIRED` | `SOURCE_REQUIRED` | `PASS_NO_PROMOTION` | Density, yield, ultimate, modulus and Poisson fields are null for 16/16 records. |
| `MAS-01` | mass_cg_inertia | `MODEL_DERIVED` | `SUPPORTED` | `PASS_NO_PROMOTION` | Volumes are present and labeled geometry-derived. |
| `MAS-02` | mass_cg_inertia | `SOURCE_REQUIRED` | `SOURCE_REQUIRED` | `PASS_NO_PROMOTION` | Assembly and part MASS/CG/INERTIA are null; derivation blocked by density/material. |
| `JNT-01` | joints_fasteners | `MODEL_DERIVED` | `PARTIALLY_SUPPORTED` | `PASS_NO_PROMOTION` | Nine interface IDs are defined at model/design level; no released fastener numerical evidence. |
| `JNT-02` | joints_fasteners | `SOURCE_REQUIRED` | `SOURCE_REQUIRED` | `PASS_NO_PROMOTION` | All fastener fields are null; status includes supplier-controlled/test-required. |
| `ABS-01` | absorber | `MODEL_DERIVED` | `PARTIALLY_SUPPORTED` | `PASS_NO_PROMOTION` | 180 mm is explicitly DESIGN_REFERENCE_ONLY; physical_validation false. |
| `ABS-02` | absorber | `TEST_REQUIRED` | `TEST_REQUIRED` | `PASS_NO_PROMOTION` | Characterization is specified; all results are TEST_REQUIRED and raw physical curves are absent. |
| `LCK-01` | lock | `TEST_REQUIRED` | `TEST_REQUIRED` | `PASS_NO_PROMOTION` | Ten lock conditions are specified; every result is TEST_REQUIRED. |
| `REB-01` | rebound | `TEST_REQUIRED` | `TEST_REQUIRED` | `PASS_NO_PROMOTION` | Seven rebound tests are specified; every result is TEST_REQUIRED; no assumed damping/friction. |
| `OEM-01` | oem_vehicle | `OEM_INPUT_REQUIRED` | `OEM_INPUT_REQUIRED` | `PASS_NO_PROMOTION` | All eight values are null and EXTERNAL_SOURCE_REQUIRED with NO_GUESSING true. |
| `CAE-01` | cae | `SOURCE_REQUIRED` | `CAE_BLOCKED` | `PASS_NO_PROMOTION` | Geometry is available; material/density external-required; mass/CG/inertia blocked; contact semantics pending; friction/absorber physical-test-required; load/pulse/restraint external-required; acceptance undefined. |
| `PROM-01` | zero_bypass | `SUPPORTED` | `SUPPORTED` | `NO_VIOLATION_FOUND` | Controlled records preserve nulls and explicit SOURCE_REQUIRED/TEST_REQUIRED/OEM_REQUIRED states. |
| `BASE-01` | immutability | `SUPPORTED` | `SUPPORTED` | `NO_VIOLATION_FOUND` | Working tree clean; baseline expected hash equals recorded before/after hash; R4.2 not authorized/absent; RC-006 and historical artifacts unchanged in this audit. |

## Domain decisions

| Domain | Decision | Reason |
|---|---|---|
| Materials | `SOURCE_REQUIRED` | 16/16 records lack authoritative identity, grade, density and properties. |
| Mass/CG/inertia | `SOURCE_REQUIRED` | Volume is model-derived; mass, CG and inertia remain null and blocked by density/material. |
| Joints/fasteners | `SOURCE_REQUIRED` | Nine interfaces are model-defined; fastener identity, grade, preload, torque and load-path evidence are absent. |
| Absorber | `TEST_REQUIRED` | F-x, F-v, hysteresis, rate, temperature, cycling and energy results are absent. |
| Lock | `TEST_REQUIRED` | Engagement, retention, partial engagement, release, jam and fault-containment results are absent. |
| Rebound | `TEST_REQUIRED` | Velocity, acceleration, energy return, damping/friction and reset results are absent. |
| OEM/vehicle | `OEM_INPUT_REQUIRED` | All eight required IDs have null values and explicit no-guessing status. |
| CAE | `CAE_BLOCKED` | Missing material/density/mass/CG/inertia, physical behavior, load/pulse/restraint and acceptance inputs. |

## Zero-bypass checks

- Assumption promoted to fact: **not found**.
- Model-derived promoted to physical: **not found**.
- Reference promoted to validated: **not found**.
- Target promoted to test result: **not found**.
- Literature promoted to T-OCS result: **not found**.
- Unverified source promoted to authoritative: **not found**.
- `180 mm`: design reference only; not physical validation.
- `18–22 kN`: no validated T-OCS physical result found.

## Repository immutability

- V7-R3: unchanged.
- R4.1: unchanged; baseline hash lock passes.
- R4.2: absent/not authorized.
- CAD/STEP/Geometry/Design Intent/Engineering Values: unchanged.
- RC-006 and historical artifacts: unchanged.
- No CAE execution, physical validation, manufacturing release, safety claim, R4.2 or RC-007 performed.

## Final verdict

**EVIDENCE_GATE_FAILED** — the failure is due to missing required authoritative/physical/OEM evidence and blocked CAE dependencies, not an observed zero-bypass promotion violation.
