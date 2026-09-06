# P0-V0 Manufacturing and Assembly Work Package

## Configuration freeze

| Field | Value |
|---|---|
| P0-V0 serial | `P0V0-20260906-001` |
| Configuration ID | `P0V0-CFG-001` |
| R4.1 artifact | `R4.1/R4.1.step` |
| R4.1 SHA256 | `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68` |
| Geometry evidence | `P0_GEOMETRY_EXTRACTION/B01_B06_extraction_pass1.json` |
| Geometry second pass | `P0_GEOMETRY_EXTRACTION/B01_B06_second_pass_verification.json` |
| Occupant testing | Prohibited |
| Crash testing | Prohibited |
| FE/MBD | Prohibited |
| R4.2 | Prohibited |
| R4.1 mutation | Prohibited |

The serial and configuration identifiers are administrative identifiers assigned to this exploratory article. They do not imply fabrication, validation, production equivalence, or physical-input closure.

## Release boundary

This work package authorizes preparation, controlled fabrication, assembly, inspection, and low-risk motion planning only. It does not authorize a dynamic event, occupant test, crash test, FE/MBD run, or claim of design validation.

No part may be fabricated from an inferred dimension. Class-A geometry must be taken from the SHA-matched R4.1 extraction or a controlled drawing derived from it. The extraction explicitly marks unavailable PMI, tolerances, travel, and clearance as unavailable; these fields must be supplied by controlled drawing/inspection evidence before they are treated as manufacturing requirements.

## Manufacturing bill of work

| Work package | Scope | Required evidence before release |
|---|---|---|
| WP-100 | Fixture/base reaction structure | Controlled drawing, material certificate, hole-pattern inspection, fixture retention review |
| WP-110L/R | Bilateral rails | Rail datum/centerline drawing, material certificate, straightness/parallelism inspection |
| WP-120 | Carriage | Interface drawing, rail-runner fit inspection, positive retention evidence |
| WP-130 | Ride-down/absorber cassette | Selected serialized absorber, axis/alignment drawing, force/stroke evidence reference |
| WP-140 | Pelvic-control geometry | Controlled part drawing and interference review |
| WP-150L/R | Seatback rotation-control links | Pivot-axis drawing, bore/fit inspection, link-load sensor provisions |
| WP-160 | Seatback frame/surrogate interface | Pivot and control-path drawing, rigidity/retention inspection |
| WP-170 | Multi-state lock | Engagement drawing, state sensing, positive retention and release-control inspection |
| WP-180 | Rebound-control module | Stop/rebound drawing, travel-limit inspection, controlled engagement direction |
| WP-190 | Restraint-interface test mount | Qualified test anchor documentation; no human occupancy authorization |
| WP-210 | Sensor/trigger interface | Sensor datum drawing, connector/strain relief inspection, synchronization plan |
| WP-J/F | Hinge, pivot, fastener and instrumentation adapters | Supplier certificates, lot traceability, preload/fit record |

## Assembly sequence

1. Freeze `P0V0-CFG-001` and record every part number, revision, lot, and serial.
2. Inspect the fixture/base and establish the approved reference datum.
3. Install the left and right rails without altering their controlled spacing or axes.
4. Measure rail straightness, parallelism, level, and positive guide engagement.
5. Install the carriage and verify free travel without absorber or lock engagement first.
6. Install the ride-down/absorber cassette and verify the absorber axis and clevis engagement.
7. Install the seatback pivot and left/right rotation-control links.
8. Install the lock and record locked, released, and controlled transition states without energy input.
9. Install rebound-control and hard stops; verify no unintended contact through the allowed manual range.
10. Install the surrogate/instrumentation interfaces and route cables outside moving/contact envelopes.
11. Record torque/preload, shims, clearances, component substitutions, and deviations.
12. Apply the L1 inspection checklist. Stop on any safety or structural anomaly.
13. Execute L2 only after L1 is formally PASS and the responsible authority releases the low-risk motion run.

## L1 and L2 execution boundary

`L1` is dimensional and assembly inspection. It requires a fabricated and assembled article plus inspection evidence. `L2` is low-risk manual or controlled mechanical motion only after L1 PASS. No high-energy actuator, sled, occupant, crash, or dynamic pulse is permitted.

Current status: `L1=NOT_EXECUTED`, `L2=NOT_EXECUTED`.

## Required evidence package

For this serial/configuration, retain:

- controlled drawings and revisions;
- material and fastener certificates;
- part and lot register;
- assembly photos;
- datum setup photos;
- calibrated measurement files;
- torque/preload records;
- L1 checklist and signed result;
- L2 motion log and signed result;
- anomaly records using `OBSERVED → MEASURED → CLASSIFIED → EVIDENCE`;
- SHA256 manifest for every raw file, photo, video, and report.

## Stop rules

Stop immediately for structural cracking, permanent deformation, rail disengagement, uncontrolled release, binding that prevents safe stop, hard-stop damage, unexpected interference, loss of fixture retention, or any measurement outside its calibrated range. Do not label the event a design failure without a defined criterion and evidence.
