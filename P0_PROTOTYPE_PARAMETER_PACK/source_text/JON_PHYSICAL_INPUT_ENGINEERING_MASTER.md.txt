# JON_PHYSICAL_INPUT_ENGINEERING_MASTER

**Scope:** Physical-input engineering definition layer for immutable R4.1. This document defines acquisition/derivation contracts only. It does not close inputs, authorize R4.2, or authorize FE/MBD execution.

## Authority and baseline
- Peter_ChatGPT: architecture / gate / authorization authority.
- JON_ChatGPT: development / innovation / engineering implementation.
- Manus-Tasks: execution / simulation / testing / validation / evidence.
- Immutable artifact: `R4.1/R4.1.step`
- Asserted R4.1 SHA256: `fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`
- Mutation allowed: `FALSE`
- Native design axes: `+X = FWD`, `+Y = LEFT`, `+Z = UP`; use the native STEP origin and never rebase the model.

> **Evidence warning:** the V4 baseline record reports `BASELINE_HASH_MATCH=false` even though it asserts the same R4.1 SHA. This master therefore does not claim independent byte-level verification of the STEP file from the currently uploaded evidence package.

## Current state
- `PHYSICAL_INPUT_CLOSURE = BLOCKED`
- `FINAL_VALIDATION_BY_MANUS = 0/14`
- FE/MBD execution remains blocked by the physical-input gate.
- Current intake audit: 14 rejected because canonical templates were supplied without authoritative external artifacts.

## A. 14-input master matrix

| INPUT_ID | CATEGORY | REQUIRED QUANTITY | UNIT | TARGET TYPE | CAD DERIVABLE | SOURCE/TEST | STATUS |
|---|---|---|---|---|---|---|---|
| mass | SOURCE_DATA_REQUIRED | Total mass of the complete immutable R4.1 seat hardware configuration represented by the physical-input model. | kg | FULL_ASSEMBLY_MASS_PROPERTY | YES | SOURCE_DATA_REQUIRED | CAD_DERIVABLE |
| CG | SOURCE_DATA_REQUIRED | Three-dimensional center-of-mass location of the complete R4.1 assembly, expressed in the native R4.1 coordinate system. | m | FULL_ASSEMBLY_CENTER_OF_GRAVITY | YES | SOURCE_DATA_REQUIRED | CAD_DERIVABLE |
| inertia | SOURCE_DATA_REQUIRED | Full symmetric 3x3 mass-moment-of-inertia tensor of the complete R4.1 assembly about its CG, expressed in the native R4.1 axes. | kg*m^2 | FULL_ASSEMBLY_INERTIA_TENSOR | YES | SOURCE_DATA_REQUIRED | CAD_DERIVABLE |
| density | SOURCE_DATA_REQUIRED | Mass density for every R4.1 FE material region/unique material definition used in the physical model. | kg/m^3 | PER_MATERIAL_REGION | NO | SOURCE_DATA_REQUIRED | EXTERNAL_SOURCE_REQUIRED |
| vehicle_pulse | SOURCE_DATA_REQUIRED | Authoritative longitudinal vehicle/seat-base acceleration time history for the applicable R4.1 rear-impact load case. | m/s^2 versus s | LONGITUDINAL_CRASH_PULSE | NO | SOURCE_DATA_REQUIRED | EXTERNAL_SOURCE_REQUIRED |
| initial_velocity | SOURCE_DATA_REQUIRED | Initial translational velocity vector of the R4.1 seat assembly at the defined load-case start time t0, in the R4.1 native frame. | m/s | INITIAL_TRANSATIONAL_VELOCITY_VECTOR | NO | SOURCE_DATA_REQUIRED | EXTERNAL_SOURCE_REQUIRED |
| bolt_stiffness | SOURCE_DATA_REQUIRED | Axial tangent stiffness of each represented fastener along its bolt axis, with bolt compliance separated from clamped-member/joint compliance. | N/m | PER_FASTENER_AXIAL_STIFFNESS | NO | SOURCE_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |
| friction | TEST_DATA_REQUIRED | Coefficient(s) of friction for every R4.1 contact pair that carries tangential relative motion in the physical model. | dimensionless coefficient | PER_CONTACT_PAIR_FRICTION_LAW | NO | TEST_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |
| contact_stiffness | TEST_DATA_REQUIRED | Normal contact force-displacement law for every R4.1 normal contact pair represented in the physical model. | N/m | PER_CONTACT_PAIR_NORMAL_STIFFNESS_LAW | NO | TEST_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |
| contact_damping | TEST_DATA_REQUIRED | Normal contact damping law for every R4.1 normal contact pair where dissipative normal contact is represented. | N*s/m | PER_CONTACT_PAIR_NORMAL_DAMPING_LAW | NO | TEST_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |
| restitution | TEST_DATA_REQUIRED | Coefficient of restitution for each discrete R4.1 impact interface for which impact separation is modeled. | dimensionless | PER_IMPACT_INTERFACE_RESTITUTION | NO | TEST_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |
| absorber_Fx | TEST_DATA_REQUIRED | Complete installed force-displacement characteristic of each R4.1 energy absorber instance, measured along its working axis through the design stroke. | N versus m | ABSORBER_FORCE_DISPLACEMENT_LAW | NO | TEST_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |
| absorber_Fv | TEST_DATA_REQUIRED | Rate-dependent force-velocity characteristic of each R4.1 energy absorber instance over the authoritative operating velocity range. | N versus m/s | ABSORBER_FORCE_VELOCITY_LAW | NO | TEST_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |
| joint_compliance | TEST_DATA_REQUIRED | Generalized compliance of each structurally relevant R4.1 compliant joint, expressed as relative 6-DOF displacement/rotation per applied generalized force/moment. | m/N, rad/(N*m), with coupled 6-DOF matrix convention explicitly defined | PER_JOINT_6DOF_COMPLIANCE_MATRIX | NO | TEST_DATA_REQUIRED | PHYSICAL_TEST_REQUIRED |

## B. Input-by-input engineering definitions

### 1. `mass` - CAD_DERIVABLE
**REQUIRED_PHYSICAL_QUANTITY:** Total mass of the complete immutable R4.1 seat hardware configuration represented by the physical-input model.
**R4_1_TARGET:** R4.1 assembly mass properties / full assembly configuration
**TARGET_TYPE:** FULL_ASSEMBLY_MASS_PROPERTY
**COORDINATE_SYSTEM:** Scalar; configuration is the immutable R4.1 assembly. Reference pose shall be the exact R4.1 configuration used by the MBD load case.
**SIGN_CONVENTION:** Positive scalar mass.
**LOADING_CONDITION:** Same structural hardware/load configuration as the R4.1 MBD model; no occupant mass or external fixture may be included unless physically represented by R4.1 and explicitly mapped.
**CONFIGURATION:** All 62 R4.1 STEP solids as mapped by the verified geometry set; exact material assignments and assembly state from authoritative CAD/BOM.
**SOURCE_OR_TEST:** SOURCE_DATA_REQUIRED
**CAD_DERIVABLE:** YES
**ORIGINAL_DESIGN_BASIS:** The design package shows the complete seat architecture, dual S1/S2 rail system, base/frame, recline/pan mechanisms, lock/unlock system and energy-absorber module. The design drawings also provide material strategy/BOM information. A deterministic CAD mass-property computation is legitimate only after exact per-solid material densities are authoritative and the exact R4.1 solid set/configuration is used. Design basis: صور.pdf pp.1, 9-15.
**REQUIRED_EVIDENCE:** Either (a) a released CAD mass-properties export computed from the exact R4.1 STEP and released material-density assignment, with geometry/material provenance, or (b) a calibrated full-assembly weigh record tied to the same R4.1 configuration.
**TEST_METHOD_IF_REQUIRED:** Calibrated full-assembly weighing on a level fixture; record tare, fixture configuration, support reactions and assembly state. Physical weighing is verification if CAD-derived mass is the primary value.
**EQUIPMENT:** CAD mass-property engine with exact R4.1 STEP + released material cards, or calibrated load-cell/weigh-scale system for verification.
**CALIBRATION:** For test route: current calibration traceable to the laboratory calibration system; record calibration ID/date/validity. For CAD route: record software/version, calculation settings, STEP hash and material-card hashes.
**RAW_DATA_FORMAT:** JSON summary plus native CAD export/PDF report; if test route, CSV/TXT time-stamped load-cell channels plus calibration metadata.
**IDENTIFICATION_METHOD:** Sum of per-solid m_i from volume*density in the released CAD configuration; cross-check against independent physical mass where available.
**UNCERTAINTY_REQUIREMENT:** Evidence must state a numerical measurement/calculation uncertainty and its basis; reject if uncertainty is absent or configuration sensitivity is undocumented.
**ACCEPTANCE_RULE:** Accept only when value is traceable to the exact R4.1 configuration, all contributing material densities are released, units are kg, provenance/hash/revision/date are present, and independent mass reconciliation is either passed or explicitly documented as not applicable. Never close from a guessed or typical mass.
**CANONICAL_FILENAME:** mass_properties.json
**BLOCKER:** Current evidence contains no released mass value; current intake audit rejected the template because no external artifact was supplied.
**MANUS_MAPPING:** `VALUE=null; UNIT=kg; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 2. `CG` - CAD_DERIVABLE
**REQUIRED_PHYSICAL_QUANTITY:** Three-dimensional center-of-mass location of the complete R4.1 assembly, expressed in the native R4.1 coordinate system.
**R4_1_TARGET:** R4.1 assembly CG and reference frame
**TARGET_TYPE:** FULL_ASSEMBLY_CENTER_OF_GRAVITY
**COORDINATE_SYSTEM:** R4.1 native/global design coordinate system: +X = FWD, +Y = LEFT, +Z = UP; right-hand rotation convention.
**SIGN_CONVENTION:** CG components are signed coordinates from the native R4.1 origin: +X forward, +Y left, +Z up.
**LOADING_CONDITION:** Exact R4.1 assembly configuration used at the MBD event start; if joint states change but mass distribution does not, report the exact pose used for the load case.
**CONFIGURATION:** Same 62-solid R4.1 hardware set used for mass; no inferred occupant/trim/fixture additions.
**SOURCE_OR_TEST:** SOURCE_DATA_REQUIRED
**CAD_DERIVABLE:** YES
**ORIGINAL_DESIGN_BASIS:** The design package explicitly defines the vehicle/seat axes X(FWD), Y(LEFT), Z(UP) and depicts the complete seat architecture. CG is therefore a legitimate CAD mass-property output once authoritative density assignment and exact configuration are known. Design basis: صور.pdf pp.1, 9, 13-15.
**REQUIRED_EVIDENCE:** Released CAD CG report or a calibrated multi-point CG measurement tied to the exact R4.1 configuration, with the native coordinate-frame definition.
**TEST_METHOD_IF_REQUIRED:** Measure support reactions at multiple known support points or equivalent certified CG rig; solve for CG coordinates in the native R4.1 frame.
**EQUIPMENT:** CAD mass-property engine or calibrated multi-point load-cell/CG rig.
**CALIBRATION:** CAD route: exact STEP/material hashes and calculation settings. Test route: current calibration of all load channels and geometric datum measurements.
**RAW_DATA_FORMAT:** JSON/CSV with X,Y,Z; raw support reactions and datum geometry for physical route; PDF/CAD report for CAD route.
**IDENTIFICATION_METHOD:** First mass moments about the native R4.1 origin; physical route solves CG from support-reaction equilibrium.
**UNCERTAINTY_REQUIREMENT:** Report component-wise uncertainty and datum/origin uncertainty; reject if origin definition or coordinate transform is not documented.
**ACCEPTANCE_RULE:** Accept only when X/Y/Z are tied to the immutable R4.1 native frame, exact configuration is identified, mass basis is traceable, and all provenance/revision/date/hash fields are complete. The PDF example geometry is not itself a numerical CG source.
**CANONICAL_FILENAME:** cg_properties.json
**BLOCKER:** No released CG artifact is currently supplied; native STEP origin must be verified against the actual R4.1 file before numerical acquisition.
**MANUS_MAPPING:** `VALUE=null; UNIT=m; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 3. `inertia` - CAD_DERIVABLE
**REQUIRED_PHYSICAL_QUANTITY:** Full symmetric 3x3 mass-moment-of-inertia tensor of the complete R4.1 assembly about its CG, expressed in the native R4.1 axes.
**R4_1_TARGET:** R4.1 assembly inertia tensor, reference point and axes
**TARGET_TYPE:** FULL_ASSEMBLY_INERTIA_TENSOR
**COORDINATE_SYSTEM:** R4.1 native/global design coordinate system: +X = FWD, +Y = LEFT, +Z = UP; right-hand rotation convention. Tensor origin: assembly CG; axes parallel to native R4.1 axes.
**SIGN_CONVENTION:** Return the full Cartesian inertia tensor matrix, not separately signed products of inertia. This removes product-of-inertia sign ambiguity. Rotations follow right-hand convention.
**LOADING_CONDITION:** Exact MBD initial configuration and assembly material state. Tensor is about the CG; if solver requires another point, perform a documented parallel-axis transformation.
**CONFIGURATION:** All 62 R4.1 solids with released densities/material regions.
**SOURCE_OR_TEST:** SOURCE_DATA_REQUIRED
**CAD_DERIVABLE:** YES
**ORIGINAL_DESIGN_BASIS:** The original package defines the complete structural assembly and native X/Y/Z axes. Inertia is deterministically calculable from exact geometry plus authoritative density assignment; no laboratory value is needed solely because the model is multibody. Design basis: صور.pdf pp.1, 9-15.
**REQUIRED_EVIDENCE:** Released CAD inertia export containing the 3x3 tensor, CG, axes, mass, STEP hash, material hashes, configuration and calculation software/version; alternatively a certified inertia measurement linked to the same configuration.
**TEST_METHOD_IF_REQUIRED:** Multi-axis inertia measurement (for example certified torsional/pendulum method) with independent axis/datum verification.
**EQUIPMENT:** CAD mass-property engine or calibrated inertia-measurement rig.
**CALIBRATION:** Physical route requires traceable calibration of torque/angle/time channels and reference inertias; CAD route requires exact geometry/material provenance.
**RAW_DATA_FORMAT:** JSON tensor object plus native CAD report or raw test time histories and reduction report.
**IDENTIFICATION_METHOD:** Direct CAD integration of r^2 and cross terms or validated rigid-body inertia identification from measured oscillation/dynamic data.
**UNCERTAINTY_REQUIREMENT:** Report tensor element uncertainty/correlation and any transformation uncertainty.
**ACCEPTANCE_RULE:** Accept only with a complete symmetric tensor tied to the exact CG and native axes, traceable mass basis, and complete provenance. Do not substitute scalar rotational inertias or literature values.
**CANONICAL_FILENAME:** inertia_tensor.json
**BLOCKER:** No released inertia tensor is currently supplied; exact R4.1 material assignment and native-frame verification are prerequisites for CAD derivation.
**MANUS_MAPPING:** `VALUE=null; UNIT=kg*m^2; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 4. `density` - EXTERNAL_SOURCE_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Mass density for every R4.1 FE material region/unique material definition used in the physical model.
**R4_1_TARGET:** Every R4.1 FE material region
**TARGET_TYPE:** PER_MATERIAL_REGION
**COORDINATE_SYSTEM:** Material scalar; mapping is by R4.1 solid/material-region identity.
**SIGN_CONVENTION:** Positive scalar density.
**LOADING_CONDITION:** Material condition corresponding to the FE/MBD model; any temperature-dependent density must be explicitly tied to the authoritative operating condition.
**CONFIGURATION:** Per-solid material mapping from R4.1 CAD/BOM/material cards; no color/appearance-based inference.
**SOURCE_OR_TEST:** SOURCE_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The design package contains a material strategy/BOM and identifies material families, but a CAD geometric model alone cannot supply authoritative density. Density must come from a released material card/certificate or a traceable material test. Design basis: صور.pdf pp.9-15.
**REQUIRED_EVIDENCE:** Released material card/certificate for each distinct material region, plus exact mapping from each R4.1 solid to the released material identifier.
**TEST_METHOD_IF_REQUIRED:** If no released certificate exists: density measurement on representative material specimens from the production-equivalent material, with specimen identity and lot traceability.
**EQUIPMENT:** Certified analytical balance plus dimensional-volume or pycnometer/density-measurement apparatus, as applicable to material form.
**CALIBRATION:** Current traceable mass and volume calibration; record instrument IDs and certificates.
**RAW_DATA_FORMAT:** JSON per material plus source certificate PDF and specimen/lot CSV where tested.
**IDENTIFICATION_METHOD:** Certificate value or measured mass/volume density at the specified condition; map value to every R4.1 region explicitly.
**UNCERTAINTY_REQUIREMENT:** Report material-property uncertainty or certificate tolerance and specimen/lot traceability.
**ACCEPTANCE_RULE:** Reject any region without a released density source and explicit R4.1 mapping. Never infer density from material name, color, typical handbook value, or nominal grade alone.
**CANONICAL_FILENAME:** material_density.json
**BLOCKER:** The current package provides no released material-density artifact; material names shown in the drawings are not sufficient to close numerical density.
**MANUS_MAPPING:** `VALUE=null; UNIT=kg/m^3; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 5. `vehicle_pulse` - EXTERNAL_SOURCE_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Authoritative longitudinal vehicle/seat-base acceleration time history for the applicable R4.1 rear-impact load case.
**R4_1_TARGET:** R4.1 applicable crash/load case
**TARGET_TYPE:** LONGITUDINAL_CRASH_PULSE
**COORDINATE_SYSTEM:** R4.1 native/global design coordinate system: +X = FWD, +Y = LEFT, +Z = UP; right-hand rotation convention. Pulse channel: longitudinal +X component at the defined vehicle/seat interface.
**SIGN_CONVENTION:** Preserve the authoritative source sign. Positive acceleration is +X; do not convert a signed acceleration trace into positive-magnitude deceleration without recording the transformation.
**LOADING_CONDITION:** The exact crash/event/load case for which R4.1 is being analyzed. Pulse must represent the boundary acceleration applied at the defined vehicle interface, not an arbitrary occupant or seat-back acceleration.
**CONFIGURATION:** R4.1 attached at its vehicle interface/base structure; pulse channel location, filtering and synchronization must match the released load-case definition.
**SOURCE_OR_TEST:** SOURCE_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The design drawings show the vehicle interface/anchor region and a crash load-path timeline. Page 8 also shows a pulse labelled as an example; that example is not acceptable as the authoritative R4.1 vehicle pulse. Design basis: صور.pdf pp.7-8, 14-15.
**REQUIRED_EVIDENCE:** Released crash sled/vehicle test pulse or authoritative load-case report containing channel identity, location, axis convention, sample interval/rate, filter settings, synchronization marker, test/configuration ID and source hash.
**TEST_METHOD_IF_REQUIRED:** Acquire the vehicle/seat-base acceleration channel during the actual released crash/sled event with the specified instrumentation and test setup.
**EQUIPMENT:** Calibrated crash accelerometer(s) and acquisition system appropriate to the released load case.
**CALIBRATION:** Traceable acceleration-channel calibration before/after test; document channel sensitivity, serial number, range, anti-aliasing and digital filtering used.
**RAW_DATA_FORMAT:** Time-stamped CSV/UFF or native acquisition format plus channel metadata, filter metadata and checksum manifest.
**IDENTIFICATION_METHOD:** Use the released time history without synthetic interpolation; any resampling/filtering required for solver ingestion must preserve the released waveform and record the exact transform.
**UNCERTAINTY_REQUIREMENT:** Report acceleration-channel measurement uncertainty and timing/sample-rate uncertainty; do not accept a plotted image without raw data.
**ACCEPTANCE_RULE:** Reject any pulse lacking authoritative load-case identity, raw time history, channel location, sign convention, filter/sample metadata or R4.1 configuration linkage. The example pulse drawn in the design package cannot close this input.
**CANONICAL_FILENAME:** vehicle_pulse.json
**BLOCKER:** EXTERNAL_DATA_REQUIRED: no authoritative vehicle pulse is supplied. The drawing contains an example only.
**MANUS_MAPPING:** `VALUE=null; UNIT=m/s^2 versus s; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 6. `initial_velocity` - EXTERNAL_SOURCE_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Initial translational velocity vector of the R4.1 seat assembly at the defined load-case start time t0, in the R4.1 native frame.
**R4_1_TARGET:** R4.1 initial-condition/load-case definition
**TARGET_TYPE:** INITIAL_TRANSATIONAL_VELOCITY_VECTOR
**COORDINATE_SYSTEM:** R4.1 native/global design coordinate system: +X = FWD, +Y = LEFT, +Z = UP; right-hand rotation convention.
**SIGN_CONVENTION:** Vector components use +X forward, +Y left, +Z up; no magnitude-only replacement is accepted.
**LOADING_CONDITION:** Same event/load-case and same t0 definition as vehicle_pulse. The event-start synchronization definition is mandatory.
**CONFIGURATION:** Exact R4.1 seat/base state at t0, including translation/pitch/lock state from the authoritative event definition.
**SOURCE_OR_TEST:** SOURCE_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The design shows time-phased seat translation/control and a crash load-path timeline, but no released numerical initial velocity is provided. Design basis: صور.pdf pp.8, 13-15.
**REQUIRED_EVIDENCE:** Authoritative initial-condition record from the crash/sled test or load-case definition, explicitly linked to the same pulse and t0.
**TEST_METHOD_IF_REQUIRED:** Velocity measurement from the test event using the certified event instrumentation; alternatively use a released initial-condition record if it is the contractual source.
**EQUIPMENT:** Appropriate calibrated speed/position instrumentation for the event, or authoritative test-data acquisition system.
**CALIBRATION:** Traceable calibration for velocity/position channels; time synchronization with the crash pulse.
**RAW_DATA_FORMAT:** JSON vector plus source CSV/native acquisition data and timing metadata.
**IDENTIFICATION_METHOD:** Use the source-defined t0 velocity vector. Integration of acceleration is permitted only when the authoritative data package explicitly defines that reduction and documents drift correction.
**UNCERTAINTY_REQUIREMENT:** Report component-wise velocity uncertainty and timing uncertainty.
**ACCEPTANCE_RULE:** Reject if t0 is not explicitly synchronized to the vehicle pulse or if only a guessed zero/rest condition is supplied. Do not infer initial velocity from the vehicle narrative.
**CANONICAL_FILENAME:** initial_velocity.json
**BLOCKER:** No authoritative initial-condition value is supplied.
**MANUS_MAPPING:** `VALUE=null; UNIT=m/s; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 7. `bolt_stiffness` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Axial tangent stiffness of each represented fastener along its bolt axis, with bolt compliance separated from clamped-member/joint compliance.
**R4_1_TARGET:** Each R4.1 represented bolt/joint
**TARGET_TYPE:** PER_FASTENER_AXIAL_STIFFNESS
**COORDINATE_SYSTEM:** Per-fastener local frame: local +n follows the bolt axis from the defined reference member toward the mating member; origin at the modeled bolt-joint datum taken from R4.1.
**SIGN_CONVENTION:** Positive axial extension is tensile elongation along +n; stiffness k=dF/d(delta) about the installed operating point.
**LOADING_CONDITION:** Installed fastener geometry, thread engagement, grip length and preload representative of the R4.1 assembly. Joint-member flexibility belongs to joint_compliance, not this scalar.
**CONFIGURATION:** Every R4.1 fastener represented by a bolt/joint element in the crash load path, especially vehicle attachment, rail/base and mechanism fasteners that carry structural load.
**SOURCE_OR_TEST:** SOURCE_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The engineering/manufacturing package explicitly contains fasteners, rail/base interfaces, joint details and lock/recline/pan mechanisms; exact fastener dimensions/grades must be taken from the authoritative R4.1 BOM rather than inferred from drawing pixels. Design basis: صور.pdf pp.6-12.
**REQUIRED_EVIDENCE:** Released fastener/joint calculation or a physical axial stiffness test using the exact production-equivalent fastener, grip and thread engagement.
**TEST_METHOD_IF_REQUIRED:** Axial force-displacement test of the exact fastener configuration; measure elongation over the defined effective gauge/grip length and identify tangent stiffness over the load range used in the model.
**EQUIPMENT:** Servo-hydraulic or electromechanical tensile test frame with calibrated load cell and extensometry suitable for the fastener size.
**CALIBRATION:** Traceable load and displacement calibration; record fastener grade/lot, geometry, grip length, thread engagement and preload setup.
**RAW_DATA_FORMAT:** CSV force-displacement data plus specimen/configuration metadata and reduction report.
**IDENTIFICATION_METHOD:** Local slope of the released force-displacement curve about the model operating point; do not use nominal material modulus alone to substitute for missing fastener configuration.
**UNCERTAINTY_REQUIREMENT:** Report stiffness estimate with fit uncertainty and repeatability across identical specimens/configurations.
**ACCEPTANCE_RULE:** Reject if fastener identity/geometry/preload/configuration cannot be mapped to an R4.1 represented bolt. Reject pure literature/nominal "typical" bolt stiffness without exact fastener traceability.
**CANONICAL_FILENAME:** bolt_stiffness.json
**BLOCKER:** Current design drawings do not expose a complete, machine-readable fastener-to-body mapping in the supplied source package; exact R4.1 fastener identities must be resolved from the authoritative CAD/BOM.
**MANUS_MAPPING:** `VALUE=null; UNIT=N/m; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 8. `friction` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Coefficient(s) of friction for every R4.1 contact pair that carries tangential relative motion in the physical model.
**R4_1_TARGET:** Each R4.1 contact pair
**TARGET_TYPE:** PER_CONTACT_PAIR_FRICTION_LAW
**COORDINATE_SYSTEM:** Each pair uses a local contact frame derived from the R4.1 mating surfaces: +n is separation; tangential axes follow the local R4.1 surface datum.
**SIGN_CONVENTION:** Friction magnitude is positive scalar; tangential force opposes relative slip. Static and kinetic values must be distinguished whenever the contact law requires them.
**LOADING_CONDITION:** Normal-load, sliding-speed, surface condition and temperature state matching the authoritative R4.1 load case. No generic coefficient may be shared across dissimilar pairs.
**CONFIGURATION:** At minimum evaluate the S1/S2 rail-guide sliding interfaces and the lock/pawl-pin/cam interfaces when those surfaces are represented as contact; add every other slip contact explicitly declared by the R4.1 contact map.
**SOURCE_OR_TEST:** TEST_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The drawings explicitly show dual S1/S2 rail mechanisms, multi-state lock/unlock hardware, hinges and absorber interfaces. These create distinct contact mechanisms; one global coefficient is not physically justified. Design basis: صور.pdf pp.1, 6-8, 12-15.
**REQUIRED_EVIDENCE:** Per-contact-pair friction report with exact mating-material/surface identity, normal-load range, speed, temperature, lubrication/finish state, specimen orientation, raw force data and fitted law.
**TEST_METHOD_IF_REQUIRED:** Reciprocating/linear tribology test using production-equivalent mating materials and surface finish; test at the operating load/speed/temperature envelope defined by the load case.
**EQUIPMENT:** Instrumented linear tribometer or representative component friction test rig with calibrated normal/tangential force and displacement channels.
**CALIBRATION:** Traceable force/displacement calibration and documented surface-preparation verification.
**RAW_DATA_FORMAT:** CSV time history for normal load, tangential force, displacement and speed, plus specimen metadata and reduction JSON.
**IDENTIFICATION_METHOD:** Determine static/kinetic or velocity/load-dependent friction law directly from measured tangential-to-normal behavior; map each law to exactly one R4.1 contact pair or clearly identified family proven equivalent.
**UNCERTAINTY_REQUIREMENT:** Report repeatability, scatter and uncertainty versus load/speed/temperature. Reject plots without raw data.
**ACCEPTANCE_RULE:** Reject if two physically different R4.1 interfaces are collapsed into one coefficient without evidence of equivalence. Reject any coefficient without exact pair/configuration traceability.
**CANONICAL_FILENAME:** friction.json
**BLOCKER:** No released friction data exists; exact R4.1 contact-pair registry must be established from the native model before testing.
**MANUS_MAPPING:** `VALUE=null; UNIT=dimensionless coefficient; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 9. `contact_stiffness` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Normal contact force-displacement law for every R4.1 normal contact pair represented in the physical model.
**R4_1_TARGET:** Each R4.1 normal contact pair
**TARGET_TYPE:** PER_CONTACT_PAIR_NORMAL_STIFFNESS_LAW
**COORDINATE_SYSTEM:** Local contact frame per pair; +n is separation/positive gap, while compression/penetration is evaluated along -n. The reported scalar/law uses positive compression force versus positive compression displacement.
**SIGN_CONVENTION:** Positive normal stiffness is dF_n/delta_n under compression; report the complete F-delta curve when nonlinear.
**LOADING_CONDITION:** Same material pair, surface condition, contact area/geometry, load range and temperature relevant to the R4.1 event.
**CONFIGURATION:** At minimum rail/guide contacts and lock/pawl-pin/cam contacts where normal contact is modeled; include rail/end-stop or absorber mechanical stop contacts only when they are actual R4.1 contact entities.
**SOURCE_OR_TEST:** TEST_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The design contains guide/rail interfaces and a multi-state lock/unlock mechanism with physical mating surfaces. Page 6 shows the lock/unlock mechanism states; page 7 shows the rail system and connection details. Design basis: صور.pdf pp.6-7, 11-15.
**REQUIRED_EVIDENCE:** Per-contact-pair normal force versus normal approach/indentation data, contact geometry, specimen identity and reduction method.
**TEST_METHOD_IF_REQUIRED:** Instrumented compression/indentation test of representative mating geometry or extracted production-equivalent coupons/components; characterize the full model-relevant normal range.
**EQUIPMENT:** Calibrated universal test frame or instrumented impact/compression rig with displacement measurement at the actual contact datum.
**CALIBRATION:** Traceable force and displacement calibration; geometry/datum verification of the contact surfaces.
**RAW_DATA_FORMAT:** CSV force-displacement data plus geometry/configuration metadata and reduction report.
**IDENTIFICATION_METHOD:** Fit the measured nonlinear normal law directly; avoid replacing nonlinear contact with one arbitrary linear spring unless the measured law demonstrably supports it over the model range.
**UNCERTAINTY_REQUIREMENT:** Report fit uncertainty, measurement uncertainty and any preload/geometry dependency.
**ACCEPTANCE_RULE:** Reject a contact law if the exact pair/normal direction/geometry cannot be mapped to R4.1. Reject a generic "steel-on-steel" stiffness.
**CANONICAL_FILENAME:** contact_stiffness.json
**BLOCKER:** No released contact stiffness laws are supplied and the exact contact surfaces are not represented in a machine-readable contact registry in the provided package.
**MANUS_MAPPING:** `VALUE=null; UNIT=N/m; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 10. `contact_damping` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Normal contact damping law for every R4.1 normal contact pair where dissipative normal contact is represented.
**R4_1_TARGET:** Each R4.1 normal contact pair
**TARGET_TYPE:** PER_CONTACT_PAIR_NORMAL_DAMPING_LAW
**COORDINATE_SYSTEM:** Local contact normal frame; use relative normal compression velocity at the same contact datum as contact_stiffness.
**SIGN_CONVENTION:** Positive damping coefficient opposes closing/separating relative normal motion according to the defined solver convention; raw evidence must preserve force/velocity signs.
**LOADING_CONDITION:** Same contact pair and geometry as contact_stiffness, with impact/closing velocities spanning the authoritative event range.
**CONFIGURATION:** Each R4.1 normal contact pair that exhibits dissipative impact behavior; do not assign damping to frictional-only contacts unless justified by the physical test.
**SOURCE_OR_TEST:** TEST_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The design includes discrete locking/contact mechanisms and guide interfaces that can produce impact transients. Page 6 explicitly shows locked/prepared/ride-down state transitions, making contact damping a pair-specific dynamic property. Design basis: صور.pdf pp.6-8, 13-15.
**REQUIRED_EVIDENCE:** Per-contact-pair force-time and relative-normal-velocity/indentation history, with identified damping law and exact impact geometry.
**TEST_METHOD_IF_REQUIRED:** Instrumented controlled impact or high-rate compression test on representative contact geometry over the authoritative relative-speed range.
**EQUIPMENT:** High-speed actuator/drop/impact rig with calibrated force and high-bandwidth displacement/velocity measurement.
**CALIBRATION:** Dynamic calibration/traceability of force and displacement channels; sampling/anti-aliasing metadata.
**RAW_DATA_FORMAT:** High-rate CSV/native acquisition file plus synchronized time base and reduction JSON.
**IDENTIFICATION_METHOD:** Identify damping from dissipative force relative to normal compression/separation velocity after removing the elastic contact contribution using the simultaneously measured force-displacement law.
**UNCERTAINTY_REQUIREMENT:** Report dynamic-channel uncertainty, time synchronization uncertainty and fit residuals.
**ACCEPTANCE_RULE:** Reject if damping is inferred from a generic damping ratio or copied from another contact without paired dynamic evidence. Reject plots without raw synchronized force/velocity data.
**CANONICAL_FILENAME:** contact_damping.json
**BLOCKER:** No released normal damping law exists.
**MANUS_MAPPING:** `VALUE=null; UNIT=N*s/m; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 11. `restitution` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Coefficient of restitution for each discrete R4.1 impact interface for which impact separation is modeled.
**R4_1_TARGET:** Each relevant R4.1 impact interface
**TARGET_TYPE:** PER_IMPACT_INTERFACE_RESTITUTION
**COORDINATE_SYSTEM:** Local contact normal defined by the mating geometry at first contact.
**SIGN_CONVENTION:** Scalar 0..1-like parameter only after direct identification; define using pre- and post-impact relative normal velocities with the source sign convention preserved.
**LOADING_CONDITION:** Actual impact configuration, closing speed, normal load and surface/material condition corresponding to the R4.1 mechanism state.
**CONFIGURATION:** Likely candidates include lock pawl/pin/cam engagement and any rail or absorber hard-stop contact actually modeled. Exact interfaces are selected from the R4.1 contact map, not from generic assumptions.
**SOURCE_OR_TEST:** TEST_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The design package depicts discrete lock-state transitions and mechanical stop/end conditions; page 6 shows the EOK mechanism states and page 7 shows the rail/absorber module. Design basis: صور.pdf pp.6-8, 13-15.
**REQUIRED_EVIDENCE:** Impact-test report with pre/post relative normal velocities, contact identity, geometry, surface condition, test-speed range and synchronization.
**TEST_METHOD_IF_REQUIRED:** Instrumented low/mid-speed impact test of production-equivalent interface or representative contact geometry; derive from measured relative velocities immediately before and after impact.
**EQUIPMENT:** High-speed displacement/velocity measurement and calibrated force/acceleration instrumentation with synchronized acquisition.
**CALIBRATION:** Traceable timing and velocity/displacement calibration; verify frame rate/sampling and impact alignment.
**RAW_DATA_FORMAT:** Synchronized high-rate CSV/native acquisition data plus derived pre/post velocity values.
**IDENTIFICATION_METHOD:** e = -v_rel,post/v_rel,pre using the explicitly documented normal direction; report the range and conditions rather than a single unsupported constant.
**UNCERTAINTY_REQUIREMENT:** Report timing and velocity uncertainty and repeatability across identical impacts.
**ACCEPTANCE_RULE:** Reject if the interface is not an actual modeled impact pair or if restitution is inferred from generic material tables.
**CANONICAL_FILENAME:** restitution.json
**BLOCKER:** No released impact-test data or validated impact-pair map is supplied.
**MANUS_MAPPING:** `VALUE=null; UNIT=dimensionless; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 12. `absorber_Fx` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Complete installed force-displacement characteristic of each R4.1 energy absorber instance, measured along its working axis through the design stroke.
**R4_1_TARGET:** Each R4.1 absorber instance
**TARGET_TYPE:** ABSORBER_FORCE_DISPLACEMENT_LAW
**COORDINATE_SYSTEM:** Local absorber axis: +x_a is the force-producing compression direction from the installed extended reference position toward stroke completion; origin at the absorber installed rod-end/clevis load datum.
**SIGN_CONVENTION:** Positive force = compressive load carried through the absorber; positive displacement = absorber stroke/compression from the defined installed reference.
**LOADING_CONDITION:** Exact mounting geometry and loading direction shown by the R4.1 absorber module; characterize the complete installed absorber, including spring/crush/damper elements that contribute to the measured axial law.
**CONFIGURATION:** R4.1 Energy Absorber - Damper + Energy Absorber - Crush Element as physically assembled in the specified mount orientation. Design drawings show a 180 mm total stroke; this is a design-stroke constraint, not a closed force curve.
**SOURCE_OR_TEST:** TEST_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** Pages 3-5 and 7 of the design package explicitly depict the progressive energy absorber, mounting interfaces, force-stroke curve, load path and staged absorption. The drawings show a 180 mm total stroke and separate low/primary/high/end-cushion regions, but the plotted curve is marked typ/target/example and is not accepted as physical evidence. Design basis: صور.pdf pp.3-5, 7.
**REQUIRED_EVIDENCE:** Instrumented installed-absorber force-displacement test over the full design stroke, with raw synchronized force/displacement data, orientation, specimen ID/lot, temperature, preconditioning and cycle history.
**TEST_METHOD_IF_REQUIRED:** Displacement-controlled axial characterization through the full 180 mm design stroke in the actual installed load direction; include loading and unloading/reload segments so hysteresis and end effects are visible. Use the authoritative test protocol for the loading rate; do not invent a rate in the evidence.
**EQUIPMENT:** Servo-hydraulic/electromechanical actuator with calibrated high-capacity load cell and independent stroke/displacement measurement.
**CALIBRATION:** Traceable force and displacement calibration before/after characterization; record fixture alignment and actuator compliance.
**RAW_DATA_FORMAT:** Time-stamped CSV/native acquisition: time, force, stroke, actuator position; JSON reduction + PDF report.
**IDENTIFICATION_METHOD:** Construct F(x) from raw force/stroke data; identify nonlinear regions, preload, hysteresis and end-cushion behavior directly. Do not fit an analytical law unless fit residuals and data coverage are reported.
**UNCERTAINTY_REQUIREMENT:** Report force/stroke uncertainty, repeatability, hysteresis and specimen-to-specimen scatter; include temperature condition.
**ACCEPTANCE_RULE:** Reject if full stroke is not covered, mounting interfaces differ from R4.1, or raw force/stroke data are missing. The drawing curve and design target numbers are not substitutes for measured evidence.
**CANONICAL_FILENAME:** absorber_Fx.json
**BLOCKER:** No released absorber F-x evidence is supplied.
**MANUS_MAPPING:** `VALUE=null; UNIT=N versus m; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 13. `absorber_Fv` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Rate-dependent force-velocity characteristic of each R4.1 energy absorber instance over the authoritative operating velocity range.
**R4_1_TARGET:** Each R4.1 absorber instance
**TARGET_TYPE:** ABSORBER_FORCE_VELOCITY_LAW
**COORDINATE_SYSTEM:** Same local absorber axis and stroke datum as absorber_Fx.
**SIGN_CONVENTION:** Positive force = compression; positive velocity = increasing compression along +x_a.
**LOADING_CONDITION:** Actual absorber mounting, stroke region, temperature and velocity range experienced in the authoritative crash/load case.
**CONFIGURATION:** Complete installed energy absorber instance, including the damping cartridge and crush/spring elements whose combined axial response is seen at the two mounting interfaces.
**SOURCE_OR_TEST:** TEST_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The drawings explicitly describe a progressive energy absorber with a damping cartridge, force-stroke curves, hydraulic/flow details, mounting interfaces and staged energy-management behavior. The load-path drawings establish that absorber force is transmitted between seat structure and vehicle/base. Design basis: صور.pdf pp.3-8, 14-15.
**REQUIRED_EVIDENCE:** Dynamic absorber characterization at the authoritative velocity range, with synchronized force, stroke and velocity data, temperature, preconditioning/cycle count and fixture/mounting configuration.
**TEST_METHOD_IF_REQUIRED:** Dynamic axial actuator characterization using the exact installed absorber; execute multiple controlled velocity levels spanning the load-case operating range, and record full stroke history. The velocity levels must be derived from the authoritative load case, not invented here.
**EQUIPMENT:** High-bandwidth servo-hydraulic test actuator, calibrated load cell, displacement sensor and independent velocity derivation/measurement.
**CALIBRATION:** Dynamic force/displacement/time calibration and documented sampling/filter chain; verify actuator tracking.
**RAW_DATA_FORMAT:** High-rate synchronized CSV/native acquisition with time, force, stroke, velocity and temperature channels.
**IDENTIFICATION_METHOD:** Identify F(v,x,T) behavior or the minimum F(v) representation required by Manus from direct data; capture hysteresis and rate dependence rather than assuming separability.
**UNCERTAINTY_REQUIREMENT:** Report dynamic force/velocity uncertainty, timing uncertainty, repeatability and temperature sensitivity.
**ACCEPTANCE_RULE:** Reject any F-v law not tied to the exact absorber instance and authoritative operating-speed range. Do not use the plotted design curve as physical evidence.
**CANONICAL_FILENAME:** absorber_Fv.json
**BLOCKER:** No released rate-dependent absorber force data is supplied and the authoritative crash-speed range is missing.
**MANUS_MAPPING:** `VALUE=null; UNIT=N versus m/s; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

### 14. `joint_compliance` - PHYSICAL_TEST_REQUIRED
**REQUIRED_PHYSICAL_QUANTITY:** Generalized compliance of each structurally relevant R4.1 compliant joint, expressed as relative 6-DOF displacement/rotation per applied generalized force/moment.
**R4_1_TARGET:** Each R4.1 compliant joint
**TARGET_TYPE:** PER_JOINT_6DOF_COMPLIANCE_MATRIX
**COORDINATE_SYSTEM:** Per-joint local frame: origin at the modeled joint center/datum; +axes from R4.1 geometry; rotations follow right-hand convention. Define generalized displacement q=[dx,dy,dz,rx,ry,rz]^T and generalized load p=[Fx,Fy,Fz,Mx,My,Mz]^T.
**SIGN_CONVENTION:** Compliance matrix C maps q=C*p in the defined local frame. Translation positive along +local axes; rotation positive by right-hand rule. The complete matrix must be reported if coupling terms are physically present.
**LOADING_CONDITION:** Same installed joint configuration, preload, joint angle and load direction as the R4.1 model state being validated.
**CONFIGURATION:** Structurally relevant compliant joints identified by the original design architecture: back recline hinges L/R, pan front hinges L/R, rail/base/seat-frame attachment joints, absorber rod-end/clevis joints and locking-mechanism pivot/joint interfaces where their compliance is represented. Exact instance mapping must come from the R4.1 joint map.
**SOURCE_OR_TEST:** TEST_DATA_REQUIRED
**CAD_DERIVABLE:** NO
**ORIGINAL_DESIGN_BASIS:** The original design explicitly labels back recline hinges, pan front hinges, rail locking mechanism, translation/lateral units, absorber mounting interfaces and compliance/control elements. Page 12 provides dedicated joint details; page 15 lists multi-DOF mechanism axes. Design basis: صور.pdf pp.6-7, 11-15.
**REQUIRED_EVIDENCE:** Per-joint compliance characterization report containing local frame, preload/configuration, generalized load/displacement data, matrix identification method, raw channels and exact R4.1 joint mapping.
**TEST_METHOD_IF_REQUIRED:** Instrument each representative joint assembly in the exact mechanical stack; apply controlled multiaxial or sequential independent loads in the defined local frame and measure relative translation/rotation across the joint datum. Identify the full 6-DOF compliance matrix where coupling matters.
**EQUIPMENT:** Calibrated multi-axis load/torque cell or dedicated single-axis fixtures with optical/DIC/LVDT/encoder displacement and angular measurement.
**CALIBRATION:** Traceable force/torque and displacement/angle calibration; verify joint datum and local-frame geometry.
**RAW_DATA_FORMAT:** Synchronized CSV/native data for applied forces/moments and relative 6-DOF motions, plus JSON matrix and reduction report.
**IDENTIFICATION_METHOD:** Linear/tangent compliance matrix from measured generalized load-displacement response about the model operating point; nonlinear or angle-dependent behavior must be reported as a law or operating-point family rather than collapsed without evidence.
**UNCERTAINTY_REQUIREMENT:** Report each matrix element uncertainty/correlation and repeatability; document any constrained DOF.
**ACCEPTANCE_RULE:** Reject if the local joint frame, preload/state or R4.1 joint identity is missing. Reject scalar spring substitutes where measured coupling or multiple DOF compliance is material.
**CANONICAL_FILENAME:** joint_compliance.json
**BLOCKER:** No released joint-compliance data exists; exact compliant-joint instance map is not present in the current evidence package.
**MANUS_MAPPING:** `VALUE=null; UNIT=m/N, rad/(N*m), with coupled 6-DOF matrix convention explicitly defined; SOURCE_ID=null; REVISION=null; DATE=null; SHA256=null; PROVENANCE=null; CONFIGURATION_MAPPING=null; VALIDATION_STATUS=BLOCKED`

## C. CAD-derivable candidates

- `mass`: YES, provided exact R4.1 geometry and released per-solid densities/material mapping are available.
- `CG`: YES, using the same mass-property basis and native R4.1 frame.
- `inertia`: YES, using the same mass-property basis; return a full symmetric tensor about CG in native axes.

`density` is **not** legitimately derived from geometry alone; it requires released material-property evidence. All other listed physical inputs require source/test evidence.

## D. SOURCE-DATA REQUIRED

`mass`, `CG`, `inertia`, `density`, `vehicle_pulse`, `initial_velocity`, `bolt_stiffness`

## E. TEST-DATA REQUIRED

`friction`, `contact_stiffness`, `contact_damping`, `restitution`, `absorber_Fx`, `absorber_Fv`, `joint_compliance`

## F. Design ambiguities preventing zero-interpretation acquisition

- Native R4.1 STEP origin/body/joint IDs are not included in the currently uploaded evidence package; axes are visible in the design package as X=FWD, Y=LEFT, Z=UP, but exact machine-readable origin and instance IDs must be taken from the authoritative R4.1 STEP/joint map.
- The supplied V4 baseline record reports BASELINE_HASH_MATCH=false; this blocks an independent claim that the currently assessed STEP bytes equal the asserted authority hash.
- The current physical_input_register contains an extra bolt_joint_stiffness identifier outside the mandated 14-input set. The master uses only bolt_stiffness and treats this as register/schema hygiene issue.
- The design package includes example/typical/target curves and performance numbers for the absorber and example crash pulse; these are not physical evidence and are deliberately not used as closure values.
- Exact fastener identities, exact R4.1 contact-pair IDs, and exact compliant-joint IDs require the authoritative CAD/BOM/model map before test specimens can be traced without interpretation.

## G. Exact artifact/file structure for Manus

```text
EXTERNAL_DATA_INTAKE_PACKAGE/
├── source_data_required/
│   ├── mass_properties.json
│   ├── cg_properties.json
│   ├── inertia_tensor.json
│   ├── material_density.json
│   ├── vehicle_pulse.json
│   ├── initial_velocity.json
│   ├── bolt_stiffness.json
├── test_data_required/
│   ├── friction.json
│   ├── contact_stiffness.json
│   ├── contact_damping.json
│   ├── restitution.json
│   ├── absorber_Fx.json
│   ├── absorber_Fv.json
│   ├── joint_compliance.json
└── evidence/<INPUT_ID>/{raw,reduction,report}/ + manifest with SHA256, provenance and R4.1 mapping
```

Every canonical JSON artifact must carry exactly the Manus intake fields: `VALUE`, `UNIT`, `SOURCE_ID`, `REVISION`, `DATE`, `SHA256`, `PROVENANCE`, `CONFIGURATION_MAPPING`, `VALIDATION_STATUS`. This master keeps those handoff fields blocked/null and never writes `FINAL_VALIDATION_BY_MANUS`.

## H. Traceability

`ORIGINAL DESIGN -> R4.1 -> PHYSICAL INPUT -> EVIDENCE -> MANUS VALIDATOR -> FUTURE MBD/FE`

Design-source basis: `صور.pdf`, 15-page visual engineering package; architecture/axes pp.1, 13-15; absorber pp.3-5,7; lock/joint/rail pp.6-7,11-12; crash-load-path pp.8,14-15. PDF SHA256: `9f14f3ca685621e4fb842d18a250cda3f76bc10faa5b96b499906ca57ed8367d`.

## I. Numerical values found in source material

No authoritative numerical value for any of the 14 physical inputs was found in the supplied evidence package. The drawings do contain design/target/example numbers, including an absorber design stroke of 180 mm. That is a configuration/test envelope, not a closed absorber F-x/F-v value. The crash pulse drawn on p.8 is labelled as an example and is not accepted as the authoritative `vehicle_pulse`.

## Governance defect requiring correction

The current `physical_input_register.json` contains an extra `bolt_joint_stiffness` identifier in addition to the mandated `bolt_stiffness`. That is a schema/governance defect. The authoritative 14-input set used here is the command-defined set and matches the Manus intake package expectation: `bolt_stiffness` only.

## Final engineering verdict

The physical-input layer is now **defined**, not **closed**. The definition is sufficient to prevent a later test team from choosing generic coefficients, arbitrary load pulses, nominal densities, or an invented absorber curve. What remains blocked is the acquisition/validation of the actual physical numbers and their evidence hashes.
