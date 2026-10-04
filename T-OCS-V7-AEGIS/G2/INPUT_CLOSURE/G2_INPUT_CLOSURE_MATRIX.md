# G2 INPUT CLOSURE MATRIX

**Decision:** BLOCKED

## Status counts
| Status | Count |
|---|---:|
| VERIFIED | 1 |
| PARTIALLY_VERIFIED | 5 |
| DEFINED_NOT_EXECUTABLE | 11 |
| MISSING | 12 |
| CONFLICTED | 2 |
| NOT_REQUIRED_FOR_DECLARED_SCOPE | 9 |

## Matrix
| ID | Field | Status | Sources | Finding |
|---|---|---|---|---|
| K-001 | body_component_list | PARTIALLY_VERIFIED | SRC-008, SRC-001 | V7 refs: 100,110L,110R,120,130L,130R,140,150L,150R,160,170,180,190,200,210,J. Current historical SR11 G2 scaffold contains 32 model-defined bodies and is NOT_PROVEN. |
| K-002 | declared_joint_list | CONFLICTED | SRC-008, SRC-002 | V7 audit: J_RAIL,J_SEATBACK,J_LINK_L,J_LINK_R. Historical SR11 scaffold: 6 joint records including two rail joints and J-LOCK170. |
| K-003 | parent_child_graph | DEFINED_NOT_EXECUTABLE | SRC-009, SRC-010, SRC-008 | Nine interface relations are model/design-defined; no executable multibody graph exists. |
| K-004 | graph_connectivity | DEFINED_NOT_EXECUTABLE | SRC-009, SRC-010, SRC-003 | Multiple parent relations occur without explicit closed-loop constraint equations. |
| K-005 | DOF_count | DEFINED_NOT_EXECUTABLE | SRC-008, SRC-004 | Declared DOF/limits are metadata only; dof_map.json is empty. |
| K-006 | constrained_DOF | DEFINED_NOT_EXECUTABLE | SRC-009, SRC-003 | Interface constrained/coupled DOF statements exist; no global executable constraint set. |
| K-007 | generalized_coordinates | DEFINED_NOT_EXECUTABLE | SRC-009, SRC-004 | q1/q2 are conceptually defined; executable generalized-coordinate map absent. |
| K-008 | coordinate_convention | VERIFIED | SRC-011 | +X=FWD,+Y=LEFT,+Z=UP; right-hand rotation convention explicitly stated. |
| K-009 | sign_convention | PARTIALLY_VERIFIED | SRC-011, SRC-009 | Right-hand convention explicit; individual joint signs are not released in one executable map. |
| K-010 | units | PARTIALLY_VERIFIED | SRC-010, SRC-011 | Units exist per source record, but no single executable G2 unit manifest. |
| K-011 | native_origin | MISSING | SRC-011, SRC-009 | Native R4.1 origin is mandated, but exact machine-readable origin/instance mapping is not released. |
| K-012 | joint_axes | PARTIALLY_VERIFIED | SRC-009, SRC-010 | q1/q2 axes are conceptually identified; exact local attachment/pin axes are not released. |
| K-013 | reference_frames | MISSING | SRC-005, SRC-011 | frames.json is empty; executable reference-frame definitions absent. |
| K-014 | J_SEATBACK_axis | PARTIALLY_VERIFIED | SRC-009, SRC-010 | IF-R4-09 gives model axis [0,1,0]; physical pin axis remains unreleased. |
| K-015 | travel_limits | DEFINED_NOT_EXECUTABLE | SRC-008 | Limits are metadata only; executable numeric limit set not proven. |
| K-016 | end_stops | DEFINED_NOT_EXECUTABLE | SRC-012, SRC-010 | End-stop requirements are definition-level; executable stop geometry/reaction model not closed. |
| K-017 | lock_states | DEFINED_NOT_EXECUTABLE | SRC-012, SRC-009 | Lock-170 states/fault conditions are definition/test items, not executable state equations. |
| K-018 | rebound_states | DEFINED_NOT_EXECUTABLE | SRC-012, SRC-010 | Rebound states are conceptual; behavior laws not closed. |
| K-019 | contact_interference | DEFINED_NOT_EXECUTABLE | SRC-008, SRC-010 | Positive overlap/interface geometry exists; intended contact versus hard interference unresolved. |
| K-020 | initial_position | MISSING | SRC-007, SRC-006 | No authoritative numerical initial pose. |
| K-021 | initial_orientation | MISSING | SRC-007, SRC-006 | No authoritative numerical initial orientation. |
| K-022 | initial_generalized_coordinates | MISSING | SRC-007, SRC-004 | No executable q1/q2 initial values. |
| K-023 | initial_velocity | MISSING | SRC-008, SRC-007 | Vehicle pulse/initial velocity is explicitly external-data-required and absent. |
| K-024 | initial_angular_velocity | MISSING | SRC-007 | No authoritative angular-velocity initial condition. |
| K-025 | initial_constraint_state | MISSING | SRC-003, SRC-007 | No executable initial constraint-state vector. |
| K-026 | initial_lock_rebound_state | MISSING | SRC-009, SRC-007 | Initial executable lock/rebound assignment absent. |
| D-001 | mass | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-007, SRC-008 | Pure kinematic/constraint gate excludes dynamic force/moment prediction. |
| D-002 | CG | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-007, SRC-008 | Pure kinematic/constraint gate excludes dynamic force/moment prediction. |
| D-003 | inertia_tensor | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-007, SRC-008 | No dynamic rotational equation is accepted at this gate. |
| D-004 | load_definitions | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-007, SRC-008 | No applied-load history is part of declared kinematic-only scope. |
| D-005 | force_laws | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-008 | No force-law prediction in declared kinematic scope. |
| D-006 | torque_laws | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-008 | No torque-law prediction in declared kinematic scope. |
| D-007 | damping | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-008, SRC-012 | Excluded from declared kinematic-only gate. |
| D-008 | friction | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-008, SRC-010 | Excluded from pure kinematic replay; remains unresolved for dynamic/contact extension. |
| D-009 | contact_parameters | NOT_REQUIRED_FOR_DECLARED_SCOPE | SRC-008, SRC-010 | Excluded from pure kinematic replay; not closed for contact dynamics. |
| E-001 | solver_timestep | MISSING | SRC-006, SRC-007 | No authoritative G2 execution timestep. |
| E-002 | integration_method | MISSING | SRC-006, SRC-007 | No authoritative G2 integration method. |
| E-003 | numerical_tolerances | MISSING | SRC-006, SRC-007 | No authoritative G2 numerical tolerance set. |
| E-004 | output_channels | DEFINED_NOT_EXECUTABLE | SRC-008, SRC-009 | Kinematic/state outputs exist conceptually; executable output map absent. |
| E-005 | source_baseline_identity | CONFLICTED | SRC-014, SRC-011 | Current governing lock matches stated R4.1 SHA; historical evidence contains BASELINE_HASH_MATCH=false. |

## Deterministic consistency
- Unit consistency: PARTIALLY_VERIFIED.
- Coordinates: VERIFIED only for the stated +X=FWD,+Y=LEFT,+Z=UP convention; exact machine-readable frames remain missing.
- Connectivity: DEFINED_NOT_EXECUTABLE.
- Joint inventory: CONFLICTED (4 declared V7 joints vs 6 historical SR11 records).
- DOF: DEFINED_NOT_EXECUTABLE; dof_map.json is empty.
- Limits/end-stops/lock/rebound: DEFINED_NOT_EXECUTABLE.
- Initial state: MISSING.
- Dynamic parameters: NOT_REQUIRED_FOR_DECLARED_SCOPE for the kinematic-only gate.
- Source consistency: CONFLICTED due historical R4.1 hash-match assertion.

## Acceptance
G2 INPUT CLOSURE = BLOCKED. No numerical engineering values were invented or promoted.
