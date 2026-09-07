# P0 Lock Implementation Options

| Option | Source basis | Function | Advantage | Risk | Required test | Authority required |
|---|---|---|---|---|---|---|
| Cam/pawl | V5/V7 generic disclosure | State capture/release | Explicitly disclosed alternative | Timing/load capacity unresolved | Static capture/release and transition test | Yes |
| Wedge | V5/V7 generic disclosure | State capture/release | Explicitly disclosed alternative | Interface/preload unresolved | Capture/release and binding test | Yes |
| Dog-clutch/pin | V5/V7 generic disclosure | State capture/release | Explicitly disclosed alternative | Alignment/engagement unresolved | Engagement and misalignment test | Yes |
| Ratchet/equivalent | V5/V7 generic disclosure | State capture/release | Explicitly disclosed alternative | Reverse control unresolved | Transition/rebound test | Yes |

No final mechanism is selected by convenience.
