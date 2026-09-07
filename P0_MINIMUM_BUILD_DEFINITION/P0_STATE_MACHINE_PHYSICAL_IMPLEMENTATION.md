# P0 STATE MACHINE — PHYSICAL IMPLEMENTATION

The source architecture supports a six-state sequence, while numeric transition timing is not consistently released. No timing is invented.

| State | Physical purpose | Minimum mechanism | Observable |
|---|---|---|---|
| S0 Normal | Normal adjustment/retention | carriage + lock in normal state | position, lock state |
| S1 Armed/Capture | Prepare event transition | trigger/lock mechanism | acceleration/trigger/state |
| S2 Pelvis Lock/Capture | Preserve pelvis relationship | 140 + restraint/guides in integrated version | pelvic/guidance state |
| S3 Ride-Down | Controlled longitudinal translation | 120 + 110L/R + 130 | x(t), v(t), a(t), S1/S2, absorber F-x |
| S4 Rotation/Rebound | Control seatback and return | 150L/R + 160 + 180 | theta(t), link loads, rebound velocity |
| S5 Secure/Post-event | Prevent uncontrolled motion | 170 + end-stop/rebound retention | lock state, residual motion |

## Transition rule

The source permits event recognition by inertial/acceleration trigger, sensor fusion, seat position/time logic, ΔV/crash-pulse detection or an optional E-trigger. The exact trigger threshold is **not assumed here**.

`UNRESOLVED_TRIGGER_TIMING = TRUE`

## Development article rule

A bench article may validate the mechanical sequence and observable transitions without claiming that the source trigger thresholds are correct for a vehicle crash.

Source architecture and state definitions are documented in the P0 master and V6/V7 source boards. fileciteturn32file11L1-L9
