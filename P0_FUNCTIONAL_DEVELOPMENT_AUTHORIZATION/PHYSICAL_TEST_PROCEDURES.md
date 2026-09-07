# P0 Functional Development Test Procedures

All procedures are developmental only. No result may be labeled qualification, certification, homologation, or safety validation.

## Universal pre-test gate

Confirm configuration ID, serial, R4.1 SHA, actual materials, deviations, fixture datum, instrument IDs, calibration status, raw-data destination, operator, laboratory, emergency controls, and stop criteria. If any item is missing, status is `BLOCKED`.

## Test sequence

1. Geometry/fit and assembly/alignment.
2. Static function and joint/lock function.
3. Absorber, friction, contact/damping, structural-development, and material characterization only with dedicated fixtures and calibrated instrumentation.
4. Crash testing is not executable under this package until an external qualified-facility gate is complete.

Every measured value must bind source, method, instrument, calibration, raw data, units, uncertainty/tolerance, timestamp, operator/lab, artifact path, SHA256, and R4.1 mapping.

## Crash hard gate

Required before any crash test: approved plan, exact article configuration, instrumentation and DAQ, facility capability, applicable protocol, safety/emergency controls, raw-data retention, and independent test record. Missing any item means `CRASH_TEST = BLOCKED`.

## Stop rules

Stop for cracking, permanent deformation, rail disengagement, uncontrolled release, unsafe binding, hard-stop damage, unexpected interference, fixture-retention loss, or out-of-calibration measurement. Record `OBSERVED -> MEASURED -> CLASSIFIED -> EVIDENCE`.
