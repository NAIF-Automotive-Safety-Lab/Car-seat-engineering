# P0 STATE MACHINE — CONFIGURATION CONTROL

S0 → S1 → S2 → S3 → S4 → S5

| State | Physical intent | Core hardware | Trigger/timing status |
|---|---|---|---|
| S0 Normal | Normal mechanism state | 120, 170, 110L/R | Defined conceptually |
| S1 Armed/Capture | Prepare controlled event transition | 170, 210 | Trigger mechanism candidate; timing UNRESOLVED |
| S2 Ride-down | Controlled longitudinal translation | 110L/R, 120, 130 | Dynamic trigger threshold UNRESOLVED |
| S3 Rotation Control | Control seatback rotation | 150L/R, 160 | Transition condition UNRESOLVED |
| S4 Rebound Control | Control reverse motion | 180, 120 | Rebound criterion UNRESOLVED |
| S5 Secure/Post-event | Stable inspected state | 170, end-stops, structural path | Timing UNRESOLVED |

The values ≥2.5 g and <20 ms remain DESIGN_REFERENCE_ONLY. They are not P0 measured or validated values.

No transition timing has been promoted to a requirement.
