# JON_FIRST_REAL_PHYSICAL_TEST_PROTOTYPE_MASTER

**Status:** DESIGN / PROTOTYPE ENGINEERING AUTHORIZED — FIRST REAL PHYSICAL TEST ARTICLE

**R4.1:** IMMUTABLE — SHA256 `fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

## 1. What exactly are we building?
A non-occupant, modular, instrumentable **P0-E / V0 integrated exploratory mechanism** derived from the R4.1 architecture. It is a physical test article, not a production seat and not R4.2. Its job is to expose mechanical truth through controlled tests.

The source architecture explicitly defines bilateral longitudinal paths, a distinct ride-down energy path, a separate seatback-rotation-control path, a multi-state lock, rebound control, pelvic-control geometry and a sensor/trigger interface. The stable reference numerals are 100–210 as defined in the V7 addendum.

## 2. Ruthless correction to the previous P0 posture
The P0 must **not** be frozen as though it were production hardware before testing. The correct architecture is modular: test one variable, measure the response, change the component, and retest. A single non-modifiable prototype would be poor experimental engineering because it would confound design change with learning.

## 3. First physical article
**P0-E / V0** contains real hardware for 100, 110L, 110R, 120, 130, 150L, 150R, 160, 170, 180, J and 210. 130/170/180 are explicitly replaceable. 150L/150R and joint stacks are configurable. Early testing is occupant-free and bench/subsystem controlled.

## 4. Source/reference values — not measured results
The concept package contains nominal/reference values such as 180 mm rearward stroke, 210 mm total travel, a -10° to +25° seatback motion envelope, 520/560 mm package references, a 25.6 g example pulse, and absorber force references. These remain design/source values only. The V7 document explicitly warns not to present such values as T-OCS measured crash results.

## 5. Exploratory architecture
- **100 / Vehicle/test-fixture interface base frame** — TEST-CRITICAL, TEST-FIXTURE, MODULAR; Anchors the P0 mechanism to a controlled bench or sled fixture and provides repeatable datum retention.
- **110L / Left longitudinal rail/load path S1** — TEST-CRITICAL, AS-DESIGNED, INSTRUMENTED, REPLACEABLE; Guided longitudinal load path for the carriage.
- **110R / Right longitudinal rail/load path S2** — TEST-CRITICAL, AS-DESIGNED, INSTRUMENTED, REPLACEABLE; Second bilateral longitudinal load path for load sharing/torsional control.
- **120 / Seat carriage** — TEST-CRITICAL, INSTRUMENTED, REPLACEABLE; Moving member translating on 110L/110R and carrying seat/absorber reactions.
- **130 / Ride-down energy-management cartridge interface** — TEST-CRITICAL, VARIABLE, REPLACEABLE, INSTRUMENTED; Allows interchangeable energy absorber configurations without changing the structural surrounding architecture.
- **140 / Pelvic-control / anti-submarining pan** — TEST-CRITICAL, AS-DESIGNED, REPLACEABLE; Carries the pelvic guidance/load-distribution geometry used in integrated mechanism studies.
- **150L / Left seatback rotation-control link** — TEST-CRITICAL, VARIABLE, REPLACEABLE, INSTRUMENTED; Independent mechanical rotation-control path on the left side.
- **150R / Right seatback rotation-control link** — TEST-CRITICAL, VARIABLE, REPLACEABLE, INSTRUMENTED; Independent mechanical rotation-control path on the right side.
- **160 / Seatback frame** — TEST-CRITICAL, INSTRUMENTED, REPLACEABLE; Carries torso/head support geometry and the independent rotation-control interfaces.
- **170 / Multi-state lock/trigger module** — TEST-CRITICAL, VARIABLE, REPLACEABLE, INSTRUMENTED; Tests the mechanical state transition architecture: normal, armed/capture, ride-down, rebound-control, secure/end state as applicable to the embodiment.
- **180 / Rebound-control module** — TEST-CRITICAL, VARIABLE, REPLACEABLE, INSTRUMENTED; Controls reverse carriage motion after primary travel.
- **190 / Occupant-restraint interface test mount** — TEST-CRITICAL, INSTRUMENTED, TEST-FIXTURE; Preserves the restraint load-path interface for non-human bench/sled studies without making a production anchorage claim.
- **200 / Torso/head guidance test geometry** — VARIABLE, TEST-CRITICAL, REPLACEABLE; Provides modular non-occupant surrogate/guide surfaces where trajectory interaction must be studied.
- **210 / Sensor/trigger interface panel** — INSTRUMENTED, REPLACEABLE, TEST-FIXTURE; Carries fixed sensor datums, connector routing, trigger targets and synchronization hardware.
- **J / Hinge/pivot/bushing/fastener set** — TEST-CRITICAL, VARIABLE, REPLACEABLE, INSTRUMENTED; Creates controlled joint interfaces whose compliance, friction and preload effects can be characterized.

## 6. Variable-at-a-time experimental logic
### P0-A — Baseline mechanism
**Hypothesis:** Establish reference function and measurement repeatability.
**Control:** Fixed baseline absorber, baseline lock, baseline rebound and baseline joint stack.
**Single variable:** None; this is the control configuration.
**Primary tests:** P0-T02, P0-T03, P0-T14, P0-T15, P0-T16, P0-T18
**Change rule:** No design change during baseline campaign except safety repair that is logged separately.
### P0-B — Alternative absorber
**Hypothesis:** Absorber characteristic is a causal driver of carriage force/displacement response.
**Control:** All structural geometry, lock, rebound, joints and test command held constant.
**Single variable:** Replace 130 cartridge only.
**Primary tests:** P0-T08, P0-T09, P0-T10, P0-T17, P0-T18
**Change rule:** Only cartridge serial/specification changes; mounting interface remains invariant.
### P0-C — Alternative lock
**Hypothesis:** Lock strategy changes transition timing/state robustness without changing absorber law.
**Control:** Same absorber cartridge, rails, joints, rebound and test command.
**Single variable:** Replace/modify 170 module only.
**Primary tests:** P0-T13, P0-T17, P0-T18
**Change rule:** Only 170 is changed; any trigger hardware change must be recorded as part of the same controlled variable.
### P0-D — Alternative rebound control
**Hypothesis:** Rebound-control strategy changes reverse excursion while primary forward path remains comparable.
**Control:** Same absorber, lock, rail, joint and test command.
**Single variable:** Replace/modify 180 only.
**Primary tests:** P0-T16, P0-T17, P0-T18
**Change rule:** Only 180 changes.
### P0-E — Integrated exploratory configuration
**Hypothesis:** The combined architecture can operate as a coupled mechanism with measurable bilateral translation, independent rotation-control and rebound control.
**Control:** Use the best documented preceding component options, each with provenance.
**Single variable:** None; this is intentionally MULTIVARIABLE and is not used for causal attribution.
**Primary tests:** P0-T14, P0-T15, P0-T16, P0-T18, P0-T19
**Change rule:** Configuration manifest freezes the exact selected components and serials before each test.

## 7. Test progression
| Level | Purpose | Example tests | Escalation rule |
|---|---|---|---|
| 1 | inspection / dimensional verification | P0-T01, P0-T19 | No dynamic testing if critical interference/retention issue remains |
| 2 | manual/mechanical function | P0-T02, P0-T13 | Proceed only after positive retention/state behavior |
| 3 | quasi-static characterization | P0-T03 to P0-T12, P0-T14-T16 | Proof-load/measurement integrity required |
| 4 | controlled dynamic | P0-T09, P0-T17, P0-T18 | Instrumentation valid, containment active, pre-approved load ceiling |
| 5+ | integrated sled / higher energy | outside initial P0 release | Requires separate technical gate and evidence review |

## 8. What counts as failure?
A failure is a test outcome, not a judgment on the whole concept. Record the mode, location, load, displacement, velocity, time, trigger, damage, permanent deformation, lock state, absorber state and evidence files. Use: CONFIRMED BY PHYSICAL TEST / CONFIRMED BY CALCULATION / SUPPORTED BY SIMULATION / ENGINEERING HYPOTHESIS / UNKNOWN.

## 9. Data path
`P0 revision → serial → configuration → calibrated sensors → raw DAQ/video → SHA256 → reduction → uncertainty → evidence artifact → Manus validation → physical-input closure → later MBD/FE correlation`

## 10. DEV-006 split
- **DEV-006-FABRICATION-EXACT:** exact R4.1-derived fabrication remains dependent on native CAD extraction.
- **P0-EXPLORATORY-FABRICATION:** allowed to use documented local adapters, modular mounts and sacrificial parts, provided the tested hypothesis remains valid and all deviations are recorded.

## 11. Hard rule
If a part of P0 is not directly supported by the native R4.1 geometry, it is never silently called an R4.1 dimension. It is tagged as R4.1-SIMPLIFIED, P0-EXPERIMENTAL or P0-TEST-FIXTURE. The article may proceed; the traceability must remain honest.

## 12. Engineering conclusion
P0 is now an experimental platform, not a frozen answer. The only acceptable path to a stronger design is physical evidence followed by controlled iteration. No result in this package is represented as an already-measured value.