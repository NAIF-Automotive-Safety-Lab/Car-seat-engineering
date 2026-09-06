# JON_P0_PHYSICAL_PROTOTYPE_MASTER

**Version:** P0.1  
**Date:** 2026-09-06  
**Status:** PROTOTYPE / TEST ARTICLE - CONDITIONAL MANUFACTURING RELEASE  
**Immutable reference:** R4.1/R4.1.step  
**R4.1 SHA256:** `fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`  

## 0. Executive engineering decision

P0 is a **real, independent physical test article**, not R4.2 and not a production seat. It is intentionally built as a single independently managed module so that the claim-critical architecture can be instrumented and tested without carrying the full packaging burden of a 2+1/3-seat installation. The source disclosure itself treats front/rear and seat-count packaging as embodiments rather than the core inventive center.

The core P0 mechanism is: bilateral longitudinal load paths -> moving carriage -> controlled ride-down energy path; mechanically distinct bilateral seatback rotation-control path; multi-state lock/trigger; explicit rebound-control path; pelvic-control geometry; restraint-interface continuity; instrumentation/trigger interfaces. The source disclosure explicitly identifies these relationships and maps them to physical observables. fileciteturn1file0L63-L79 fileciteturn1file0L89-L95

**Ruthless rule:** the concept boards contain nominal dimensions, targets and illustrative performance. They are not measured P0 data. The source document itself warns that generated CAD is visualization and that example pulse/stroke/force values are not T-OCS test evidence. fileciteturn1file0L71-L77

## 1. Prototype architecture

### 1.1 Primary load path
`TEST FIXTURE / VEHICLE INTERFACE (100)` -> `S1/S2 RAILS (110L/110R)` -> `SEAT CARRIAGE (120)` -> `ENERGY MANAGEMENT (130)` -> controlled longitudinal motion.

### 1.2 Distinct rotation path
`SEATBACK FRAME (160)` -> `150L/150R ROTATION-CONTROL LINKS` -> dedicated hinge/control structure. The source architecture deliberately separates this path from the ride-down absorber path. fileciteturn1file0L85-L91

### 1.3 State path
`170 MULTI-STATE LOCK` coordinates normal / armed / ride-down / rebound-control / secure functions. The source state machine defines S0 through S5 as physical state transitions with observable position, acceleration, travel, angle, link load and rebound variables. fileciteturn1file0L93-L95

### 1.4 Rebound path
`180 REBOUND CONTROL` acts on reverse carriage motion after primary stroke. The prototype will measure reverse velocity, secondary excursion and control force.

### 1.5 Occupant path
`190 RESTRAINT INTERFACE` remains structurally available in the prototype, but **no human occupant** is used in early P0 testing. The source safety section explicitly requires the restraint path to remain available and calls for staged testing before human/vehicle-crash use. fileciteturn1file0L103-L103

## 2. Source-derived nominal design basis (NOT measurements)

The concept/manufacturing boards show a nominal 180 mm ride-down stroke, 210 mm total travel envelope, 1120 mm overall height, rail/package dimensions around 520/560 mm, and a source-study crash-pulse example with 25.6 g peak. The absorber board shows 18-22 kN per-side force and 25-35 kN peak damping, and an operating temperature concept of -40 to +80 C. These are **design-study/source values only** and must never be promoted to P0 results. The V7 document explicitly classifies the 25.6 g pulse, 180 mm stroke and force-range examples as study/target material rather than measured T-OCS evidence. fileciteturn1file0L71-L77

## 3. Component list

See `P0_BOM.csv`. The functional build includes 100, 110L, 110R, 120, 130, 140, 150L, 150R, 160, 170, 180, 190, 200 and 210 plus the replaceable hinge/pivot and instrumentation hardware.

The V7 stable reference list identifies 100-210 and explicitly freezes those numerals for CAD, drawings, claims, tests and FEA. fileciteturn1file0L63-L64

## 4. Parts classified by role

**AS-DESIGNED:** 100, 110L/R, 140, 160 and structural interfaces whose geometry must match the selected R4.1 configuration.

**TEST-CRITICAL:** 100, 110L/R, 120, 130, 140, 150L/R, 160, 170, 180, 190.

**VARIABLE:** 130 cartridge, 170 lock implementation details, 180 rebound element, 200 guidance surface.

**INSTRUMENTED:** 110L/R, 120, 130L/R, 150L/R, 160, 170, 180, 190 and 210.

**REPLACEABLE:** 130, 150L/R pivot/bushing stacks, 170, 180 and instrumentation adapters.

## 5. Manufacturing definition

P0 is manufacturable only through a controlled drawing release process. The uploaded V7 PDF provides the architecture and nominal board dimensions, but it is not safe to infer every hole coordinate/tolerance by raster measurement. The source itself says the design boards are concept/engineering material and calls for a true parametric CAD assembly before formal manufacturing/patent drawing release. fileciteturn1file0L77-L77

Therefore: use the immutable R4.1/native geometry to populate fabrication-critical coordinates. Do **not** create geometry by eyeballing pixels. This is a hard manufacturing gate.

## 6. P0 sub-test articles

### P0-A - Mass properties fixture
Real P0 assembly or rigidly representative structural configuration. Produces mass, CG and inertia.

### P0-B - Contact/interface fixture
Separate specimen/fixture for friction, contact stiffness, damping and restitution, with interface-specific surfaces and normal-load control.

### P0-C - Absorber fixture
Clevis-to-clevis characterization rig for F-x, F-v, hysteresis, temperature response and cycling. The source architecture specifically calls for absorber characterization by force-stroke, time, hysteresis and temperature response. fileciteturn1file0L91-L92

### P0-D - Joint/fastener fixture
Measured bolt/joint preload and compliance fixture for hinges, rail mounts, absorber clevises and lock mounts.

### P0-E - Integrated mechanism
Full structural P0 assembly for controlled bench dynamics and later ATD/sled-level correlation.

## 7. Minimum viable instrumentation

See `P0_INSTRUMENTATION_PLAN.csv`. The required minimum is: independent S1/S2 load channels; absorber force channels; carriage displacement; seat-base and fixture acceleration; seatback angle; rotation-link loads where practical; lock-state sensing; absorber temperature; synchronized high-speed video.

Instrument ranges are sized from the documented design envelope, not from invented test results. For example, the source absorber envelope is 18-22 kN per side with a 25-35 kN peak damping concept, so a preliminary 50 kN force-channel class is a procurement sizing choice with subsequent proof-load recheck. Likewise, a 250 mm displacement class is a preliminary choice because the concept shows 180 mm ride-down and 210 mm total travel. These are **instrument selection estimates, not physical measurements**.

## 8. Test matrix

See `P0_TEST_MATRIX.csv`. The sequence is deliberate: mass properties -> interface/contact -> absorber -> joints/lock -> integrated functionality -> controlled dynamic event. This directly follows the staged validation discipline in the source: V0 CAD, V1 MBD, V2 absorber, V3 lock, V4 sled, then V5 correlation and V6 full vehicle. fileciteturn1file0L75-L75

## 9. Safety and containment

P0 tests are **unoccupied** until the responsible authority explicitly advances maturity. Begin with component/subsystem/bench tests. Hard stops are installed at both travel ends; both rails require positive guide engagement; the fixture is independently retained; dynamic runs are remote-operated and guarded; emergency stop and abort criteria are established. The source safety section specifically requires single-point/common-cause failure analysis, hard-stop analysis, over-travel/jam analysis, asymmetric-load analysis and temperature/contamination checks before occupant/vehicle crash use. fileciteturn1file0L101-L103

No numerical stop load is fabricated. The stop limit for each run is the lower of: verified component/fixture proof envelope and calibrated instrumentation safe operating range, with a defined operating margin approved before the run.

## 10. R4.1 deviation register

See `P0_R4.1_DEVIATION_REGISTER.csv`. Every deviation is explicit. No unlogged substitution is allowed. The most material open item is vehicle hardpoints: the uploaded V7 PDF does not provide enough authoritative native-coordinate data to claim geometric identity to R4.1. P0 therefore uses a dedicated traceable laboratory base fixture.

## 11. Requirements -> evidence -> Manus

`ORIGINAL CONCEPT` -> `R4.1 IMMUTABLE REFERENCE` -> `P0 PART/SUBSYSTEM` -> `TEST CONFIGURATION` -> `RAW DATA` -> `REDUCED DATA` -> `EVIDENCE ARTIFACT` -> `MANUS VALIDATOR` -> `PHYSICAL INPUT CLOSURE` -> `MBD` -> `FE CORRELATION`.

The source validation map already associates controlled displacement with x(t)/absorber F-x, rotation control with theta(t)/link loads, bilateral paths with S1/S2 reactions, rebound with reverse velocity/secondary excursion, and absorber characterization with F-x/hysteresis/temperature/cycling. fileciteturn1file0L95-L95

## 12. Exact data P0 is designed to generate

### Source/CAD-derived or measured physical properties
- Complete P0 mass, CG vector and inertia tensor, with configuration and uncertainty.
- Material-lot density where CAD material assignment is otherwise not authoritative.

### Interface/contact data
- Friction as a function of normal load, direction, velocity and temperature where tested.
- Interface-specific contact stiffness, damping and restitution laws.

### Absorber data
- Force-stroke curve.
- Force-velocity relationship.
- Hysteresis/work loss.
- Temperature sensitivity.
- Cycle-to-cycle drift.

### Joint/lock data
- Bolt/joint stiffness vs preload.
- Joint compliance/clearance/friction.
- Lock state truth table, trigger conditions and measured transition latency.

### Integrated dynamics
- x(t), v(t), a(t).
- seatback theta(t).
- left/right rail reactions.
- absorber force and work.
- rebound velocity/secondary excursion.
- lock state transitions.
- synchronized high-speed kinematics.

## 13. Manus evidence artifact structure

For each input, create:

`<INPUT_ID>__P0__<ARTICLE>__<TEST_ID>__<REVISION>__<YYYYMMDD>.json`

Required top-level fields:
`INPUT_ID`, `VALUE`, `UNIT`, `SOURCE_ID`, `REVISION`, `DATE`, `SHA256`, `PROVENANCE`, `CONFIGURATION_MAPPING`, `VALIDATION_STATUS`.

Recommended sibling files:
- raw data
- calibration certificate(s)
- test setup photo(s)
- test procedure revision
- configuration manifest
- reduction script/version
- reduced data
- uncertainty budget
- inspection report
- SHA256 manifest

## 14. What is still genuinely missing

This package deliberately does **not** pretend that the following are already known:

1. Native R4.1 hardpoint coordinates and machine-readable critical dimensions.
2. Final P0 absorber serial/configuration and its physical law.
3. Material certificates for the actual P0 build lots.
4. Verified fastener preload/torque or instrumented-bolt response for the actual joint stack.
5. External vehicle/sled pulse and initial velocity.

These are acquisition tasks. They are not excuses to stop development; P0 is designed specifically to generate the locally measurable items and to identify the exact external route for the remaining ones.

## 15. Build sequence

1. Freeze P0.1 architecture and deviation register.
2. Extract authoritative R4.1/native hardpoints and critical coordinates.
3. Manufacture 100/110L/110R/120/140/150L/R/160/170/180 and test fixtures.
4. Obtain traceable material and fastener certificates.
5. Bench-fit and inspect all joints/rails.
6. Assemble P0-A/P0-B/P0-C/P0-D.
7. Characterize absorber, joints and interfaces.
8. Assemble P0-E.
9. Run static/low-load functional proof.
10. Run integrated controlled dynamic screening.
11. Package evidence and submit to Manus.
12. Close only the physical inputs validated from accepted artifacts; then parameterize MBD and correlate FE.

## 16. Ruthless mentor verdict

The good news is that the prototype path is now concrete. The bad news is that calling the V7 picture board a manufacturing release would be engineering trash. The board is useful for architecture and nominal envelope; it is not a substitute for authoritative CAD hardpoints, tolerances, material certificates and a selected/serialised absorber. The source document itself says exactly this in substance: the concept is testable, but the true parametric CAD and physical validation still have to be built. fileciteturn1file0L77-L77

So the mandate is: **build real hardware, but build it under configuration control. Measure everything that matters. Never turn a design target into a test result.**

### Deliverables in this package
- `JON_P0_PHYSICAL_PROTOTYPE_MASTER.json`
- `JON_P0_PHYSICAL_PROTOTYPE_MASTER.md`
- `P0_BOM.csv`
- `P0_MANUFACTURING_DRAWING_RELEASE.md`
- `P0_TEST_CONFIGURATION.txt`
- `P0_INSTRUMENTATION_PLAN.csv`
- `P0_TEST_MATRIX.csv`
- `P0_R4.1_DEVIATION_REGISTER.csv`
- `P0_REQUIREMENTS_TRACEABILITY.csv`
- `P0_DATA_ACQUISITION_PLAN.csv`
