# Physical Release Requirements Extracted from A–O Gap Register v3

**Release decision: BLOCKED**. The register contains **181 gaps**; **101** are directly tied to physical evidence, physical characterization, instrumentation, safety, joints, materials, fasteners, dynamics, or physical inputs. Current verified physical measurements: **0**.

## Category summary

| Category | Requirement domain | Gaps | Physical-test-required gaps | Current state |
|---|---|---:|---:|---|
| E | joints/hinges and compliance | 11 | 0 | BLOCKED |
| F | materials and characterization | 9 | 9 | BLOCKED |
| G | fasteners, preload, and joint stiffness | 11 | 11 | BLOCKED |
| H | Lock-170 dynamics and physical response | 12 | 12 | BLOCKED |
| I | Rebound-180 dynamics and stop behavior | 11 | 11 | BLOCKED |
| J | absorber force/displacement/rate behavior | 11 | 11 | BLOCKED |
| K | instrumentation and calibration | 12 | 12 | BLOCKED |
| L | safety, containment, retention, and failure modes | 10 | 10 | BLOCKED |
| O | physical input acquisition | 14 | 14 | BLOCKED |

## Minimum evidence required before release

- V7 native CAD/STEP with SHA-256 and deterministic STEP/B-Rep audit
- released BOM, materials, fastener, joining, PMI/GD&T, dimensions, and tolerances
- joint/hinge/lock/rebound/absorber characterization where sensitivity requires it
- calibrated instrumentation plan and imported raw measurements
- V7 CAE model inputs, solver configuration, output artifact hashes, and limitations
- evidence graph links from source to parameter/CAD entity/feature/part/interface/joint/material/physics/test/measurement/result

## Category-specific requirements

### E — joints/hinges and compliance

Gap IDs: `E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11`.

Physical-test-required records: **0**.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### F — materials and characterization

Gap IDs: `F01, F02, F03, F04, F05, F06, F07, F08, F09`.

Physical-test-required records: **9**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### G — fasteners, preload, and joint stiffness

Gap IDs: `G01, G02, G03, G04, G05, G06, G07, G08, G09, G10, G11`.

Physical-test-required records: **11**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### H — Lock-170 dynamics and physical response

Gap IDs: `H01, H02, H03, H04, H05, H06, H07, H08, H09, H10, H11, H12`.

Physical-test-required records: **12**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### I — Rebound-180 dynamics and stop behavior

Gap IDs: `I01, I02, I03, I04, I05, I06, I07, I08, I09, I10, I11`.

Physical-test-required records: **11**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### J — absorber force/displacement/rate behavior

Gap IDs: `J01, J02, J03, J04, J05, J06, J07, J08, J09, J10, J11`.

Physical-test-required records: **11**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### K — instrumentation and calibration

Gap IDs: `K01, K02, K03, K04, K05, K06, K07, K08, K09, K10, K11, K12`.

Physical-test-required records: **12**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### L — safety, containment, retention, and failure modes

Gap IDs: `L01, L02, L03, L04, L05, L06, L07, L08, L09, L10`.

Physical-test-required records: **10**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

### O — physical input acquisition

Gap IDs: `O01, O02, O03, O04, O05, O06, O07, O08, O09, O10, O11, O12, O13, O14`.

Physical-test-required records: **14**.
Minimum test definition: smallest calibrated test that resolves the stated uncertainty.
Blocking dependencies: released V7 source package.
Next action: obtain authoritative V7 evidence and rerun audit.

## Interpretation

The register distinguishes computational enablement from physical closure. EVA may derive, bound, simulate, rank sensitivity, and prioritize tests, but none of those actions creates a `MEASURED_PHYSICALLY` or `VALIDATED` record. Synthetic experiments remain software tests only. Manufacturer release remains blocked until the listed authoritative and physical evidence is imported, hashed, and replayed.
