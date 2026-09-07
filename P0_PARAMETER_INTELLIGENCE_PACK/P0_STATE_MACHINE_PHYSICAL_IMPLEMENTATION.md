# P0 State-Machine Physical Implementation

| From | To | Trigger | Mechanical action | Physical component | Observable | Status |
|---|---|---|---|---|---|---|
| S0 Normal | S1 Armed/Capture | Event/arming condition not numerically specified | Engage/capture development state | 170 lock | Lock state | UNRESOLVED
| S1 Armed/Capture | S2 Ride-down | Dynamic event/primary stroke | Permit controlled translation and absorber engagement | 120/130/110L/110R | Travel, force, reactions | DISCLOSED_CONCEPT
| S2 Ride-down | S3 Rotation Control | Rotation-control condition not specified | Limit seatback rotation independently | 150L/150R/160 | Angle/link motion | DISCLOSED_CONCEPT
| S2 Ride-down | S4 Rebound Control | After primary stroke | Directional damping/locking | 180 | Reverse motion/secondary excursion | DISCLOSED_CONCEPT
| S4 Rebound Control | S5 Secure/Post-event | Secure condition not specified | Capture/secure | 170/180 | Secure state | UNRESOLVED

No trigger timing, threshold, damping coefficient, or force is invented.
