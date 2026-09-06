# P0 Data Acquisition Plan

## Raw data
Raw data are immutable/read-only after acquisition. Each raw file receives SHA256 and is referenced by the test record. No interpolation or replacement is performed in raw files.

## Required synchronized channels
- carriage displacement
- S1 rail reaction
- S2 rail reaction
- absorber force per installed cartridge
- seatback angle
- acceleration at the designated moving-mass datum
- lock state
- joint deformation channels when running P0-D/P0-E joint characterization
- absorber temperature where relevant
- high-speed video timecode where used

## Reduction
Velocity is derived only from calibrated displacement or validated inertial data. Energy is calculated only from measured force-displacement data. Every filter is retained by name/version/cutoff and applied to a new reduced-data file; raw data remain untouched.

## Manus physical-input artifact fields
Every candidate physical input output must contain:
`VALUE, UNIT, SOURCE_ID, REVISION, DATE, SHA256, PROVENANCE, CONFIGURATION_MAPPING, VALIDATION_STATUS`.

An input becomes CLOSED only through the existing Manus/Peter evidence gate, not merely because P0 produced a number.
