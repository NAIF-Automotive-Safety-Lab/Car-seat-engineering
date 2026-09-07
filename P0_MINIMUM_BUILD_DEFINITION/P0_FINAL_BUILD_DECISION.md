# P0 FINAL BUILD DECISION

## MINIMUM_PHYSICAL_ARTICLE

**P0-FD-MIN / dual-side non-occupant integrated mechanism demonstrator**

Mandatory:
100 + 110L + 110R + 120 + 130L/130R + 150L/150R + 160 + 170 + 180 + J + 210.

Optional for a later integrated occupant-interface configuration:
140 + 190 + 200.

## WHY THIS IS MINIMUM

The source architecture's core relationship is the combination of:

bilateral longitudinal load paths
+
controlled longitudinal ride-down
+
separate seatback rotation-control path
+
event-dependent state transition
+
rebound-control state.

V7 explicitly defines that combination as the central architecture. fileciteturn29file14L1-L7

Removing one rail destroys the bilateral hypothesis.
Removing 130 destroys the ride-down/energy-management hypothesis.
Removing 150L/150R or 160 destroys the separate rotation-control path.
Removing 170 destroys the state-transition function.
Removing 180 destroys the rebound function.
Removing 210 removes reproducible instrumentation of the mechanism.

## WHAT IT CAN PROVE

At development-rig level only:

1. bilateral guided translation;
2. longitudinal ride-down mechanism function;
3. separate seatback rotation-control path;
4. lock/state transition mechanics;
5. rebound-control mechanics;
6. measurable mechanical kinematics/reactions.

## WHAT IT CANNOT PROVE

It cannot prove vehicle compatibility, occupant safety, injury reduction, regulatory compliance, production durability, crashworthiness, full-vehicle performance, or FE correlation.

## MUST BE KNOWN BEFORE BUILD

- exact critical native-CAD geometry and configuration;
- fixture attachment definition;
- exact joint/hinge interfaces;
- selected/serialized absorber hardware for functional tests;
- lock/rebound hardware configuration;
- material/fastener identity under controlled build release;
- critical inspection dimensions;
- safety containment and hard-stop arrangement;
- instrumentation mounting/configuration.

## CAN BE MEASURED AFTER BUILD

- actual mass;
- CG;
- inertia;
- actual rail/ carriage motion;
- left/right reaction forces;
- absorber F-x/F-v;
- hinge/joint response;
- lock engagement/release behavior;
- rebound response;
- temperature effects of the tested absorber;
- actual as-built clearances and alignment.

## MUST WAIT FOR SPECIALIZED TEST

- interface friction law;
- contact stiffness/damping;
- restitution;
- joint/bolt compliance;
- validated vehicle pulse;
- validated initial velocity/load case;
- FE material/contact inputs;
- occupant/ATD/sled evidence.

## RELEASE STATUS

`R4.1 = FROZEN`
`R4.2 = NOT_AUTHORIZED`
`PHYSICAL_DATA = 0`
`FE = NOT_AUTHORIZED`
`SAFETY_VALIDATION = NOT_AUTHORIZED`

The P0 master explicitly says early testing is non-occupant and that exact fabrication remains blocked until native CAD control of critical geometry and interfaces is available. fileciteturn31file6L1-L7
