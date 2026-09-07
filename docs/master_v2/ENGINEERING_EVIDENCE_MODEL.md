# Engineering Evidence Model

Every record carries an evidence class, source, source SHA-256 where applicable, tool identity, timestamp, and verification state. Allowed classes are `DERIVED_FROM_CAD`, `DERIVED_FROM_DOCUMENT`, `MEASURED_PHYSICALLY`, `CALCULATED`, `SIMULATED`, `ASSUMED`, `TARGET`, `CANDIDATE`, `UNVERIFIED`, `BLOCKED`, and `VALIDATED`.

Software cannot label a record `MEASURED_PHYSICALLY` without a hashed raw measurement source, sensor/calibration identity, sampling metadata, units, timestamp, and test ID. AI candidates cannot become authoritative without deterministic geometry verification.
