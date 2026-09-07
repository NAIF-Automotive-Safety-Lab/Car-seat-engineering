# P0 Prototype Parameter Pack

**Status:** `PARTIAL` — extraction/classification package only; no manufacturing release, physical test result, safety claim, FE authorization, or R4.2.

**R4.1 SHA256:** `fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

## Integrity firewall
All physical-result counts are zero. Design references, design targets, source-study values, and unresolved numeric tokens are not promoted to physical evidence.

## Architecture

- `100` — Vehicle interface / base frame: vehicle interface
- `110L/110R` — Longitudinal load paths S1/S2: bilateral guidance
- `120` — Seat carriage: seat-to-rail/absorber interface
- `130` — Ride-down energy-management module: controlled longitudinal displacement and energy dissipation
- `140` — Pelvic-control / anti-submarining pan: pelvis trajectory and load distribution
- `150L/150R` — Seatback rotation-control links: independent rotation-limiting path
- `160` — Seatback frame: torso/head structural support
- `170` — Multi-state lock: state transitions
- `180` — Rebound-control module: reverse-motion control
- `190` — Occupant restraint interface: seat-belt/airbag interface
- `200` — Torso/head guidance: distributed trajectory control
- `210` — Sensor/trigger interface: acceleration/position/state recognition and measurement

## Negative register

- `physical mass` — `NOT_FOUND`
- `physical CG` — `NOT_FOUND`
- `physical inertia` — `NOT_FOUND`
- `validated density mapping` — `NOT_FOUND`
- `validated friction` — `NOT_FOUND`
- `validated contact stiffness` — `NOT_FOUND`
- `validated contact damping` — `NOT_FOUND`
- `validated restitution` — `NOT_FOUND`
- `measured absorber Fx` — `NOT_FOUND`
- `measured absorber Fv` — `NOT_FOUND`
- `validated vehicle pulse` — `NOT_FOUND`
- `validated initial velocity` — `NOT_FOUND`
- `validated joint compliance` — `NOT_FOUND`
- `validated bolt stiffness` — `NOT_FOUND`
- `released GD&T` — `NOT_FOUND`
- `released manufacturing BOM` — `NOT_FOUND`

## Counts

```json
{
  "TOTAL_NUMERICAL_VALUES_EXTRACTED": 846,
  "TOTAL_ENGINEERING_RELATIONSHIPS": 5,
  "TOTAL_PARAMETERS": 862,
  "TOTAL_DESIGN_TARGETS": 0,
  "TOTAL_SOURCE_STUDY_VALUES": 158,
  "TOTAL_ASSUMPTIONS": 0,
  "TOTAL_DERIVED_VALUES": 0,
  "TOTAL_PHYSICAL_MEASUREMENTS": 0,
  "TOTAL_TEST_RESULTS": 0,
  "TOTAL_NOT_FOUND_ITEMS": 16,
  "TOTAL_BUILD_CRITICAL_PARAMETERS": 12,
  "TOTAL_AUTHORITY_REQUIRED": 16,
  "TOTAL_PHYSICAL_VALIDATION_REQUIRED": 16
}
```

## Final status

```json
{
  "P0_PROTOTYPE_PARAMETER_PACK": "PARTIAL",
  "V5_EXTRACTION": "PARTIAL",
  "V7_EXTRACTION": "PARTIAL",
  "V5_V7_DELTA": "PASS",
  "TRACEABILITY": "PARTIAL",
  "PHYSICAL_RESULT_INTEGRITY": "PASS",
  "FALSE_MEASUREMENT_PROMOTION": 0
}
```
