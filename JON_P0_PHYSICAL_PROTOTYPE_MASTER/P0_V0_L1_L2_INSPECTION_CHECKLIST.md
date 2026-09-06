# P0-V0 L1/L2 Inspection Checklist

**Serial:** `P0V0-20260906-001`  
**Configuration:** `P0V0-CFG-001`  
**R4.1 SHA256:** `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

## Evidence control

| Field | Required entry |
|---|---|
| Inspector | Name/role/signature |
| Date/time | ISO 8601 with timezone |
| Instrument IDs | ID, calibration certificate, calibration due date |
| Article state | Part/lot/serial/configuration |
| R4.1 mapping | Exact B01–B06 target and source SHA |
| Raw evidence | Filename and SHA256 |
| Photos/video | Filename, view, timestamp, SHA256 |
| Result | PASS / FAIL / BLOCKED |

A blank, estimated, or memory-based entry is not evidence.

## L1 — dimensional and assembly inspection

| ID | Check | Method/equipment | Acceptance record | Status |
|---|---|---|---|---|
| L1-B01 | Rail/carriage interface | Calibrated CMM/height gauge/appropriate gauge; inspect both rails and carriage interfaces against controlled extraction/drawing | Record actuals, units, instrument ID, tolerance source, and photos; no inferred tolerance | NOT_EXECUTED |
| L1-B02 | Seatback pivot/rotation-control | Calibrated arm/angle measurement and datum inspection | Pivot axis, link alignment, retention, and orientation recorded against B02 | NOT_EXECUTED |
| L1-B03 | Ride-down/absorber interface | Calibrated alignment and fit inspection | Absorber axis, mount engagement, clevis retention, and load-path alignment recorded | NOT_EXECUTED |
| L1-B04 | Lock mating interface | Visual plus calibrated engagement/travel measurement | Locked/released states and mating engagement recorded; no unintended release | NOT_EXECUTED |
| L1-B05 | Rebound-control/stop interface | Calibrated travel/clearance measurement through manual safe range | Stop location, direction, engagement, and clearance recorded | NOT_EXECUTED |
| L1-B06 | Fixture reaction boundary | Datum and fastener/retention inspection | Fixture contact, reaction boundary, fastener state, and retention recorded | NOT_EXECUTED |
| L1-07 | Rail straightness/parallelism | Calibrated dimensional equipment | Actuals recorded; acceptance criterion must come from controlled drawing or approved inspection plan | NOT_EXECUTED |
| L1-08 | Carriage free travel | Manual low-force movement with displacement recording | Smooth travel, no binding, both rails engaged | NOT_EXECUTED |
| L1-09 | Hard stops | Manual approach under safe force limit | Stops present, retained, undamaged, and reachable before unsafe over-travel | NOT_EXECUTED |
| L1-10 | Interference/binding | Full permitted manual envelope, visual and tactile inspection | No unexpected interference or jam | NOT_EXECUTED |
| L1-11 | Fasteners and joints | Torque/preload tools with valid calibration | Part/lot/torque/preload recorded; no unapproved substitution | NOT_EXECUTED |
| L1-12 | Sensor/trigger mounts | Datum/clearance inspection | Sensors and cables clear of moving/contact envelopes | NOT_EXECUTED |

**L1 release rule:** L1 is PASS only when every mandatory row is PASS with raw evidence and no safety/structural anomaly. Otherwise L1 is BLOCKED or FAIL and L2 is prohibited.

## L2 — low-risk functional/mechanical motion

L2 begins only after signed L1 PASS. Use manual or controlled low-force motion with physical guards and an emergency stop where powered equipment is used. No occupant, crash pulse, sled event, high-energy actuator, or FE/MBD is permitted.

| ID | Motion/function | Minimum observation | Evidence | Status |
|---|---|---|---|---|
| L2-01 | Rail travel | Carriage completes the approved manual/low-force range without binding | Displacement log, video, configuration ID, SHA256 | NOT_EXECUTED |
| L2-02 | Carriage motion | Smooth bilateral engagement and no rail disengagement | Video and inspection notes | NOT_EXECUTED |
| L2-03 | Seatback pivot | Pivot moves about intended axis within approved safe envelope | Angle/time log and video | NOT_EXECUTED |
| L2-04 | Rotation-control | Left/right control path remains engaged and distinguishable from translation path | Link motion/load observation; no unsupported force claim | NOT_EXECUTED |
| L2-05 | Absorber engagement | Absorber mounts remain aligned and engaged during low-risk motion | Fit/engagement photos and motion record | NOT_EXECUTED |
| L2-06 | Lock engagement/release | Intended states engage/release under controlled manual command only | State log, trigger ID, video | NOT_EXECUTED |
| L2-07 | Rebound-control | Reverse motion is controlled by installed mechanism/stop without uncontrolled return | Low-risk reverse-motion video and travel record | NOT_EXECUTED |
| L2-08 | Interference/binding | No unexpected contact, jam, cable snag, or loss of retention | Video, inspection record, anomaly log | NOT_EXECUTED |
| L2-09 | Hard stops | Stops arrest motion safely before over-travel and remain undamaged | Stop contact record and post-run photos | NOT_EXECUTED |

**L2 release rule:** L2 is PASS only when every mandatory row is PASS, all evidence hashes are recorded, and no safety or structural anomaly exists. L3 is not automatic.

## Anomaly record

For every anomaly, use exactly:

```text
OBSERVED:
MEASURED:
CLASSIFIED: safety / structural / dimensional / functional / instrumentation / unknown
EVIDENCE: filenames + SHA256 + photo/video references
CRITERION:
DISPOSITION: STOP / HOLD / CORRECTIVE ACTION / REVIEW
```

Do not label an anomaly as a design failure without a predefined criterion and evidence.

## Current gate

```text
L1_STATUS = NOT_EXECUTED
L2_STATUS = NOT_EXECUTED
NEXT_ALLOWED_LEVEL = L1 only after fabricated/assembled article and evidence are available
OCCUPANT_TESTING = PROHIBITED
CRASH_TESTING = PROHIBITED
FE/MBD = PROHIBITED
```
