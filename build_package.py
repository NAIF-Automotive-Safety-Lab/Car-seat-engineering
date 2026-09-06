import json, hashlib, os
from pathlib import Path

ROOT = Path('/mnt/data/seat_engineering_contracts_r4_1')
CONTRACTS = ROOT/'contracts'
SCHEMAS = ROOT/'schemas'
MANIFESTS = ROOT/'manifests'

BASELINE_HASH='fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68'
STATUSES=['VERIFIED','DEFINED','IMPLEMENTATION_READY','INPUT_REQUIRED','BLOCKED','NOT_TESTED','NOT_VALIDATED']
GATES=['GATE_R4_1_BASELINE','GATE_PHYSICAL_INPUTS','GATE_FE_RUNTIME','GATE_FE_MODEL','GATE_FE_EXECUTION','GATE_FE_VALIDITY','GATE_V5_READINESS','GATE_V6_READINESS','GATE_V7_READINESS','GATE_R4_2_AUTHORIZATION','GATE_PROMOTION']


def dump(name, obj, sub=CONTRACTS):
    p=sub/name
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False, sort_keys=False)+'\n', encoding='utf-8')
    return p

schema_meta={"$schema":"https://json-schema.org/draft/2020-12/schema"}

reqs=[]
for rid, func, param, interface, model, test, evidence, criterion in [
("REQ-V5-001","Represent occupant using an explicitly specified computational representation","P-V5-REP","IF-V5-OCCUPANT","MODEL-V5-OCCUPANT","TEST-V5-READINESS","EV-V5-INPUT-CLOSURE","representation is non-null and source/provenance are present"),
("REQ-V5-002","Supply traceable anthropometry without guessed values","P-V5-ANTHRO","IF-V5-OCCUPANT","MODEL-V5-OCCUPANT","TEST-V5-ANTHRO-INPUT","EV-V5-ANTHRO-SOURCE","source is explicit and validation state is closed"),
("REQ-V5-003","Define occupant mass/CG/inertia from authoritative inputs","P-V5-MASS/P-V5-CG/P-V5-INERTIA","IF-V5-OCCUPANT","MODEL-V5-OC-INERTIA","TEST-V5-READINESS","EV-V5-PHYSICS","each critical field has value or explicit BLOCKED state; no guess"),
("REQ-V5-004","Define occupant joints and DOF","P-V5-JOINTS/P-V5-DOF","IF-V5-JOINTS","MODEL-V5-KINEMATICS","TEST-V5-KINEMATICS","EV-V5-KINEMATICS","joint definitions and DOF are complete for intended model"),
("REQ-V5-005","Define restraint and seat-contact interfaces","P-V5-RESTRAINT/P-V5-SEATCONTACT","IF-V5-REST/CONTACT","MODEL-V5-INTERFACE","TEST-V5-INTERFACE","EV-V5-INTERFACE","all required interfaces are mapped and sourced"),
("REQ-V5-006","Define initial posture/velocity/boundary conditions","P-V5-INIT/P-V5-BC","IF-V5-BC","MODEL-V5-INITIAL","TEST-V5-INITIAL","EV-V5-INITIAL","closed with explicit source and units"),
("REQ-V5-007","Define load and response channels","P-V5-LOAD/P-V5-RESP","IF-V5-CHANNELS","MODEL-V5-CHANNELS","TEST-V5-CHANNELS","EV-V5-CHANNEL-MAP","channel map is complete and unique"),
("REQ-V6-001","Provide executable sled sensor and acquisition specification","P-V6-SENSOR-SPEC","IF-V6-SENSOR","MODEL-V6-INSTRUMENTATION","TEST-V6-SETUP","EV-V6-SETUP","all mandatory sensor fields are closed"),
("REQ-V6-002","Enforce calibration and synchronization traceability","P-V6-CAL/SYNC","IF-V6-CAL","MODEL-V6-DATA-ACQ","TEST-V6-CAL-SYNC","EV-V6-CAL-SYNC","calibration reference and synchronization policy are explicit"),
("REQ-V6-003","Define accepted sled data states and provenance","P-V6-DATA-STATE","IF-V6-DATA","MODEL-V6-SLED-DATA","TEST-V6-DATA-QA","EV-V6-DATA-QA","state transitions and provenance fields are machine-readable"),
("REQ-V7-001","Correlate FE to physical hierarchy with explicit channel mapping","P-V7-MAP","IF-V7-CORRELATION","MODEL-V7-CORRELATION","TEST-V7-CORRELATION","EV-V7-CORRELATION","all compared channels have deterministic mappings"),
("REQ-V7-002","Normalize time/coordinate/sampling/filtering before metric calculation","P-V7-NORM","IF-V7-PREPROCESS","MODEL-V7-PREPROCESS","TEST-V7-PREPROCESS","EV-V7-PREPROCESS","transform sequence is explicit and auditable"),
("REQ-V7-003","Compute generic correlation metrics without fabricated thresholds","P-V7-METRICS","IF-V7-METRICS","MODEL-V7-METRICS","TEST-V7-METRICS","EV-V7-METRICS","metrics are computable; threshold remains REQUIRED when unsourced"),
("REQ-KIN-001","Generate x(t), v(t), a(t), theta(t), FL(t), FR(t) with traceable frames/signs","P-KIN-SIGNFRAME","IF-KIN","MODEL-KIN","TEST-KIN","EV-KIN","each signal declares source/frame/sign/unit/time-base"),
("REQ-KIN-002","Prevent numerical drift and coordinate/channel ambiguity","P-KIN-INTEGRATION","IF-KIN","MODEL-KIN-QA","TEST-KIN-QA","EV-KIN-QA","integration/differentiation and channel identity checks exist"),
("REQ-ENE-001","Account for energy through vehicle-to-occupant path","P-ENE-TERMS","IF-ENE-FLOW","MODEL-ENE","TEST-ENE-CONSERVATION","EV-ENE","all defined terms have equations/source and accounting sign"),
("REQ-ENE-002","Detect energy imbalance without inventing a numerical result","P-ENE-BAL","IF-ENE-BAL","MODEL-ENE-BAL","TEST-ENE-CONSERVATION","EV-ENE-BAL","imbalance is calculated only from actual inputs"),
("REQ-FLT-001","Make fault modes testable and evidence-linked","P-FLT-ATTR","IF-FLT","MODEL-FMEA","TEST-FLT","EV-FLT","required attributes exist; unknown ratings remain explicit"),
("REQ-V4-001","Deliver V4 evidence required to unlock downstream stages","P-V4-EVIDENCE","IF-V4","MODEL-V4-EVIDENCE","TEST-V4-CLOSURE","EV-V4","all minimum evidence classes are present and valid"),
("REQ-GATE-001","Apply fail-closed gates","P-GATE-STATE","IF-GATE","MODEL-GATE","TEST-GATE-NEGATIVE","EV-GATE","unknown/missing/blocked/untested/assumed never evaluate PASS"),
("REQ-NEXT-001","Identify highest-priority blocker automatically","P-NEXT-PRIORITY","IF-NEXT","MODEL-NEXT","TEST-NEXT","EV-NEXT","first unsatisfied blocker in defined priority order is returned"),
("REQ-BASE-001","Protect immutable R4.1 baseline","P-R4-1-SHA","IF-BASELINE","MODEL-BASELINE","TEST-BASELINE-INTEGRITY","EV-BASELINE","actual hash equals authoritative hash and mutation false"),
]:
    reqs.append({"requirement_id":rid,"function":func,"parameter_refs":param.split('/'),"interface_refs":interface.split('/'),"model_refs":model.split('/'),"test_refs":[test],"evidence_refs":[evidence],"acceptance_criterion":criterion,"verification_status":"DEFINED"})

traceability={
 "artifact":"REQUIREMENTS_TRACEABILITY",
 "version":"R4.1-CONTRACT-0.1",
 "authority_state":{"r4_1_baseline":"PASS","r4_1_sha256":BASELINE_HASH,"r4_1_mutation":False,"r4_2":"NOT_AUTHORIZED","promotion":"BLOCKED"},
 "allowed_statuses":STATUSES,
 "requirements":reqs,
 "orphan_rules":[
   "Every requirement_id must resolve to at least one function, parameter, interface, model, test, evidence, and acceptance criterion.",
   "Every parameter/interface/model/test/evidence ID must be referenced by at least one requirement or dependency.",
   "Unreferenced entities are ORPHAN and cause the owning readiness gate to FAIL."
 ]
}
dump('REQUIREMENTS_TRACEABILITY.json', traceability)

# V5
field_names=[
('occupant_representation','required', 'representation type and implementation identity; NULL only when blocked'),
('anthropometry_source','required','authoritative source identifier; no guessed anthropometry'),
('mass','required','numeric mass when supplied'),('CG','required','CG vector in declared frame when supplied'),('inertia','required','mass inertia tensor in declared frame when supplied'),
('joint_definitions','required','joint list with type, parent, child, axes, limits if applicable'),('degrees_of_freedom','required','explicit DOF set'),
('pelvis_reference','required','reference frame/marker definition'),('torso_reference','required','reference frame/marker definition'),('head_reference','required','reference frame/marker definition'),
('restraint_interfaces','required','belt/restraint attachment definitions'),('seat_contact_interfaces','required','contact definitions'),
('initial_posture','required','pose state with frame'),('initial_velocity','required','velocity state with frame and units'),('boundary_conditions','required','BC definitions'),
('load_channels','required','input channel map'),('response_channels','required','response channel map')]
v5_fields=[]
for n, req, desc in field_names:
    v5_fields.append({"field_id":"V5-"+n.upper(),"name":n,"value":None,"state":"INPUT_REQUIRED","unit":None,"source":None,"provenance":None,"required":req=="required","validation_status":"INPUT_REQUIRED","description":desc})
v5={
 "artifact":"V5_EXECUTION_CONTRACT","contract_version":"0.1","status":"IMPLEMENTATION_READY",
 "execution_rule":"No critical field may be silently defaulted. Missing critical inputs keep V5_MODEL_READY=FAIL.",
 "fields":v5_fields,
 "readiness_gate":{"gate_id":"GATE_V5_READINESS","expression":"ALL(required field.state IN {READY,MEASURED,VALIDATED} AND required field has source+provenance+unit where applicable)","pass_state":"PASS","else_state":"FAIL"},
 "manus_required_inputs":[f["field_id"] for f in v5_fields],
 "forbidden_actions":["guess anthropometry","substitute undocumented mass/CG/inertia","treat DEFINE as VALIDATED","consume V4 outputs with status BLOCKED/INVALID/UNVERIFIED"]
}
dump('V5_EXECUTION_CONTRACT.json',v5)

# V6
v6={
 "artifact":"V6_SLED_EXECUTION_CONTRACT","contract_version":"0.1","status":"IMPLEMENTATION_READY",
 "states":["NOT_CONNECTED","NOT_CALIBRATED","READY","MEASURED","REJECTED","VALIDATED"],
 "sensor_schema_fields":["sensor_id","sensor_type","location","axis","sampling_rate","synchronization","calibration_reference","measurement_range","expected_signal","data_format","timestamp_policy","coordinate_system","filtering_policy","channel_name","data_provenance","acceptance_criterion","state"],
 "input_spec":{"scenario_id":None,"test_id":None,"vehicle_or_rig_id":None,"seat_configuration_id":None,"acquisition_system_id":None,"clock_reference":None,"sampling_policy":None,"filter_policy":None,"sensor_manifest_ref":None,"calibration_manifest_ref":None,"state":"INPUT_REQUIRED"},
 "output_spec":{"raw_data_ref":None,"processed_data_ref":None,"channel_manifest_ref":None,"calibration_evidence_ref":None,"sync_evidence_ref":None,"qa_result":None,"state":"NOT_TESTED"},
 "calibration_schema":{"calibration_id":None,"sensor_id":None,"reference_standard":None,"date":None,"operator":None,"pre_test_status":"NOT_CALIBRATED","post_test_status":"NOT_CALIBRATED","certificate_ref":None,"traceability_ref":None},
 "data_quality_rules":["No MEASURED or VALIDATED state without raw-data provenance.","No VALIDATED state without calibration and synchronization evidence.","Reject duplicate channel names unless explicitly versioned.","Coordinate system and sign convention must be declared."],
 "manus_execution_outputs":["sensor manifest","calibration records","raw synchronized data","processed data with transformation log","QA decision","evidence hashes"]
}
dump('V6_SLED_EXECUTION_CONTRACT.json',v6)

# V7
v7={
 "artifact":"V7_CORRELATION_ENGINE_CONTRACT","contract_version":"0.1","status":"IMPLEMENTATION_READY",
 "comparison_levels":["FE","COMPONENT_TEST","SUBSYSTEM_TEST","SLED","VEHICLE_TEST"],
 "pipeline":[
  {"step":"time_alignment","status":"DEFINED","method":None,"threshold_required":False},
  {"step":"coordinate_alignment","status":"DEFINED","method":None,"threshold_required":False},
  {"step":"channel_mapping","status":"DEFINED","method":None,"threshold_required":False},
  {"step":"sampling_normalization","status":"DEFINED","method":None,"threshold_required":False},
  {"step":"filtering","status":"DEFINED","method":None,"threshold_required":False},
  {"step":"metric_calculation","status":"DEFINED","method":"metric registry","threshold_required":True}
 ],
 "signal_comparisons":["displacement","velocity","acceleration","force","energy","rotation","timing"],
 "metrics":[
  {"metric_id":"RMSE","formula":"sqrt(mean((y_model-y_ref)^2))","threshold_required":True,"result":None},
  {"metric_id":"NRMSE","formula":"RMSE / declared normalization scale","threshold_required":True,"result":None},
  {"metric_id":"PEAK_ERROR","formula":"abs(peak_model-peak_ref)/abs(peak_ref) when denominator valid","threshold_required":True,"result":None},
  {"metric_id":"PHASE_ERROR","formula":"declared phase/time lag metric","threshold_required":True,"result":None},
  {"metric_id":"CURVE_SIMILARITY","formula":"declared curve similarity metric","threshold_required":True,"result":None},
  {"metric_id":"ENERGY_DISCREPANCY","formula":"declared energy difference/ratio","threshold_required":True,"result":None},
  {"metric_id":"TIMING_DISCREPANCY","formula":"declared timing difference","threshold_required":True,"result":None}
 ],
 "acceptance_policy":"No threshold may be invented. Unsourced thresholds must remain THRESHOLD_REQUIRED=TRUE and cannot produce PASS.",
 "correlation_record_fields":["correlation_id","source_level","reference_level","channel_map_ref","time_alignment_ref","coordinate_alignment_ref","sampling_ref","filter_ref","metric_results","threshold_source_refs","acceptance_decision","evidence_refs","status"]
}
dump('V7_CORRELATION_ENGINE_CONTRACT.json',v7)

# Kinematics
signals=[]
for name in ['x(t)','v(t)','a(t)','theta(t)','FL(t)','FR(t)']:
    signals.append({"signal_id":"KIN-"+name.replace('(t)','').upper(),"name":name,"source":None,"coordinate_system":None,"sign_convention":None,"reference_frame":None,"unit":None,"sampling_definition":None,"integration_rule":None,"differentiation_rule":None,"measurement_source":None,"model_source":None,"validation_source":None,"state":"INPUT_REQUIRED"})
kin={"artifact":"KINEMATIC_EXECUTION_CONTRACT","contract_version":"0.1","architecture":{"seatback":["150L","150R"],"carriage":[120],"ride_down":["130L","130R"],"dual_rails":["110L","110R"],"vehicle_base":[100]},"signals":signals,"qa_controls":["explicit single authoritative time base","no double integration without drift control","coordinate transform log required","sign convention invariant checks","left/right identity checks","sampling synchronization check"],"derivation_policy":{"position_to_velocity":"declared numerical differentiation only","velocity_to_acceleration":"declared numerical differentiation only","acceleration_to_velocity":"single integration with baseline/drift control if used","velocity_to_position":"single integration with initial condition and drift control if used"}}
dump('KINEMATIC_EXECUTION_CONTRACT.json',kin)

# Energy
terms=["kinetic_energy","absorbed_energy","elastic_energy","plastic_dissipation","contact_work","friction_work","joint_work","residual_energy"]
energy={"artifact":"ENERGY_ACCOUNTING_CONTRACT","contract_version":"0.1","flow":["VEHICLE_INPUT","RAILS","RIDE-DOWN","ABSORBER","CARRIAGE","SEAT_PAN","PELVIS","TORSO","HEAD_SUPPORT","OCCUPANT"],"terms":[{"term_id":t,"value":None,"unit":None,"source":None,"equation":None,"sign_convention":None,"status":"INPUT_REQUIRED"} for t in terms],"checks":[{"check_id":"ENE-BAL-001","definition":"input energy - accounted downstream energy = residual/imbalance","threshold":None,"threshold_required":True,"status":"DEFINED"},{"check_id":"ENE-BAL-002","definition":"component-to-component work/energy transfer continuity","threshold":None,"threshold_required":True,"status":"DEFINED"}],"rule":"No numerical energy result is reported until actual executable/model/test data are supplied."}
dump('ENERGY_ACCOUNTING_CONTRACT.json',energy)

# Fault matrix
components=['LOCK','LEFT_RAIL','RIGHT_RAIL','RIDE-DOWN','ABSORBER','SEATBACK','JOINTS','BOLTS','WELDS','CONTACT','SENSORS','ACTUATORS','CONTROL_LOGIC']
faults=[]
for i,c in enumerate(components,1):
    faults.append({"FAULT_ID":f"FLT-{i:03d}","COMPONENT":c,"INITIATING_EVENT":None,"FAILURE_MECHANISM":None,"PHYSICAL_EFFECT":None,"DETECTABLE_SIGNAL":None,"DETECTION_METHOD":None,"SAFE_STATE":None,"DEGRADED_STATE":None,"SEVERITY":None,"OCCURRENCE":None,"DETECTABILITY":None,"REQUIRED_TEST":f"TEST-FLT-{i:03d}","EVIDENCE_REQUIRED":["failure mechanism evidence","detectable signal evidence","test evidence"],"STATUS":"INPUT_REQUIRED","RATING_POLICY":"NULL/TO_BE_ESTABLISHED until sourced evidence exists"})
fault={"artifact":"FAULT_TEST_CONTRACT","contract_version":"0.1","rating_rule":"No probability/severity/detectability values may be invented.","faults":faults,"test_record_fields":["FAULT_ID","test_id","setup_ref","stimulus_ref","instrumentation_ref","expected_detection","actual_observation","raw_data_ref","analysis_ref","evidence_ref","decision","status"]}
dump('FAULT_TEST_CONTRACT.json',fault)

# V4 dependency
v4={"artifact":"V4_DEPENDENCY_CONTRACT","contract_version":"0.1","status":"BLOCKED","required_inputs":["geometry/mesh definition","material evidence","contact definitions","connection/joint definitions","boundary conditions","solver identity/availability","analysis controls"],"required_outputs":["mesh evidence","solver identity","material evidence","contact evidence","connection evidence","boundary-condition evidence","convergence evidence","stress/strain evidence","deformation evidence","energy evidence","failure criterion evaluation"],"evidence_schema":[{"evidence_id":None,"class":None,"source_file":None,"producer":None,"execution_id":None,"timestamp":None,"hash":None,"status":"INPUT_REQUIRED","validity":None,"review_status":"NOT_VALIDATED"}],"failure_schema":{"failure_id":None,"stage":None,"failure_type":None,"message":None,"input_state":None,"recoverable":None,"root_cause_ref":None,"evidence_ref":None,"status":"INPUT_REQUIRED"},"root_cause_schema":{"root_cause_id":None,"symptom":None,"candidate_causes":[],"tests":[],"disposition":None,"evidence_refs":[],"status":"INPUT_REQUIRED"},"downstream_consumption_rule":"V5/V6/V7 MUST reject V4 results with BLOCKED, INVALID, or UNVERIFIED status."}
dump('V4_DEPENDENCY_CONTRACT.json',v4)

# Gate definition
gate_defs=[]
for gid in GATES:
    gate_defs.append({"gate_id":gid,"default":"FAIL","unknown_is_pass":False,"missing_is_pass":False,"blocked_is_pass":False,"untested_is_pass":False,"assumed_is_verified":False,"required_inputs":[],"expression":None,"decision":"FAIL"})
gate={"artifact":"AUTOMATIC_GATE_DEFINITION","contract_version":"0.1","evaluation_semantics":{"PASS":"all mandatory predicates true","FAIL":"any mandatory predicate false","UNKNOWN":"treated as FAIL"},"gates":gate_defs}
# Fill exact rules
gate_by={g['gate_id']:g for g in gate['gates']}
gate_by['GATE_R4_1_BASELINE'].update({"required_inputs":["actual baseline bytes","expected SHA256","mutation indicator"],"expression":f"sha256(actual)=={BASELINE_HASH} AND mutation==FALSE"})
gate_by['GATE_PHYSICAL_INPUTS'].update({"required_inputs":["all mandatory physical input records"],"expression":"all mandatory physical inputs CLOSED with source/unit/provenance"})
gate_by['GATE_FE_RUNTIME'].update({"required_inputs":["solver identity","runtime availability evidence"],"expression":"solver identity present AND executable runtime evidence present"})
gate_by['GATE_FE_MODEL'].update({"required_inputs":["mesh","material","contacts","connections","BCs"],"expression":"all mandatory FE model artifacts VALIDATED"})
gate_by['GATE_FE_EXECUTION'].update({"required_inputs":["execution record","solver log","result artifacts"],"expression":"execution completed and artifacts are traceable"})
gate_by['GATE_FE_VALIDITY'].update({"required_inputs":["convergence","quality checks","failure evaluation"],"expression":"validity evidence complete and no blocking validity issue"})
gate_by['GATE_V5_READINESS'].update({"required_inputs":["V5 critical fields","V4 valid outputs"],"expression":"all V5 critical inputs closed AND V4 consumed artifacts not BLOCKED/INVALID/UNVERIFIED"})
gate_by['GATE_V6_READINESS'].update({"required_inputs":["sensor manifest","calibration","sync","test setup"],"expression":"all mandatory setup fields closed AND calibration/sync ready"})
gate_by['GATE_V7_READINESS'].update({"required_inputs":["FE results","test data","mapping","preprocessing"],"expression":"all compared datasets valid and mapping/preprocessing complete"})
gate_by['GATE_R4_2_AUTHORIZATION'].update({"required_inputs":["Peter authorization record"],"expression":"explicit authorization exists AND all prerequisite gates PASS"})
gate_by['GATE_PROMOTION'].update({"required_inputs":["R4.2 authorization","verification evidence","promotion decision"],"expression":"explicit promotion authorization AND all required verification evidence PASS"})
dump('AUTOMATIC_GATE_DEFINITION.json',gate)

# Next action
next_action={"artifact":"NEXT_ACTION_ENGINE","contract_version":"0.1","priority_order":["EVIDENCE_INTEGRITY","BASELINE_INTEGRITY","PHYSICAL_INPUTS","SOLVER_AVAILABILITY","FE_MODEL_READINESS","FE_EXECUTION","FE_VALIDATION","V5_READINESS","V6_TESTING","V7_CORRELATION","R4_2_AUTHORIZATION","PROMOTION"],"selection_rule":"return first unsatisfied priority item; do not skip unresolved higher priority blockers","current_snapshot":[{"priority":1,"blocker_id":"BLOCKER-V4-EVIDENCE-INTEGRITY","blocker_description":"V4 execution/input-gap closure evidence is not present in the authoritative state; downstream validation cannot be claimed.","required_input":"V4 execution record + input-gap closure manifest + evidence hashes","responsible_agent":"Manus-Tasks","dependency":"V4 execution environment and complete inputs","unlock_condition":"V4 evidence package passes evidence integrity gate"},{"priority":3,"blocker_id":"BLOCKER-V4-PHYSICAL-INPUTS","blocker_description":"Authoritative state says V4 has input gaps; exact missing inputs must be enumerated and closed.","required_input":"machine-readable V4 input-gap list with source/unit/provenance and closure state","responsible_agent":"Manus-Tasks","dependency":"V4_REQUIRED_INPUTS","unlock_condition":"all mandatory V4 inputs CLOSED"}],"return_schema":["BLOCKER_ID","BLOCKER_DESCRIPTION","REQUIRED_INPUT","RESPONSIBLE_AGENT","DEPENDENCY","UNLOCK_CONDITION"]}
dump('NEXT_ACTION_ENGINE.json',next_action)

# Evidence / provenance
evidence={
 "artifact":"EVIDENCE_SCHEMA",
 "contract_version":"0.1",
 "record_fields":["evidence_id","evidence_type","requirement_refs","test_refs","source_uri","producer_agent","execution_id","timestamp_utc","file_name","sha256","status","validity","reviewer","review_timestamp","parent_evidence_refs"],
 "valid_statuses":["VERIFIED","DEFINED","IMPLEMENTATION_READY","INPUT_REQUIRED","BLOCKED","NOT_TESTED","NOT_VALIDATED","INVALID","UNVERIFIED"],
 "invalidation_rules":["missing hash","missing producer or source","broken requirement/test linkage","result consumed after INVALID/BLOCKED/UNVERIFIED"]
}
dump('EVIDENCE_SCHEMA.json',evidence)
provenance={"artifact":"INPUT_PROVENANCE_SCHEMA","contract_version":"0.1","record_fields":["input_id","name","value","unit","source","source_version","provenance_type","acquisition_method","acquisition_timestamp","operator_or_agent","file_ref","sha256","validation_state","required","status"],"provenance_types":["DESIGN_SOURCE","MEASUREMENT","SIMULATION","LITERATURE","STANDARD","USER_PROVIDED","DERIVED"],"rule":"DERIVED inputs must cite all parent input IDs and derivation method."}
dump('INPUT_PROVENANCE_SCHEMA.json',provenance)

# baseline report
baseline={"artifact":"BASELINE_INTEGRITY_REPORT","contract_version":"0.1","authoritative_expected_sha256":BASELINE_HASH,"authoritative_mutation_flag":False,"independent_verification":"NOT_PERFORMED","verification_reason":"Baseline bytes/file were not supplied in this conversation/package; hash cannot be recomputed here.","required_check":"sha256(actual_R4.1_bytes)==authoritative_expected_sha256 AND mutation_flag==FALSE","on_mismatch":"STOP IMMEDIATELY","status":"ASSERTED_BY_AUTHORITY_NOT_INDEPENDENTLY_VERIFIED"}
dump('BASELINE_INTEGRITY_REPORT.json',baseline)

# machine-readable JSON Schema-ish wrappers
for fn, title in [
('REQUIREMENTS_TRACEABILITY','Requirements Traceability'),('V5_EXECUTION_CONTRACT','V5 Execution Contract'),('V6_SLED_EXECUTION_CONTRACT','V6 Sled Execution Contract'),('V7_CORRELATION_ENGINE_CONTRACT','V7 Correlation Engine Contract'),('KINEMATIC_EXECUTION_CONTRACT','Kinematic Execution Contract'),('ENERGY_ACCOUNTING_CONTRACT','Energy Accounting Contract'),('FAULT_TEST_CONTRACT','Fault Test Contract'),('V4_DEPENDENCY_CONTRACT','V4 Dependency Contract'),('AUTOMATIC_GATE_DEFINITION','Automatic Gate Definition'),('NEXT_ACTION_ENGINE','Next Action Engine'),('EVIDENCE_SCHEMA','Evidence Schema'),('INPUT_PROVENANCE_SCHEMA','Input Provenance Schema'),('BASELINE_INTEGRITY_REPORT','Baseline Integrity Report')]:
    s=dict(schema_meta)
    s.update({"title":title,"type":"object","additionalProperties":True})
    (SCHEMAS/(fn+'.schema.json')).write_text(json.dumps(s,indent=2)+'\n')

# Unified package manifest + all-in-one
files=[p for p in CONTRACTS.glob('*.json')]
manifest={"package":"SEAT_ENGINEERING_EXECUTABLE_CONTRACTS","package_version":"R4.1-CONTRACT-0.1","baseline":{"sha256":BASELINE_HASH,"mutation":False},"authority_state":{"V4":"BLOCKED_BY_EXECUTION_AND_INPUT_GAPS","V5":"DEFINED_NOT_VALIDATED","V6":"DEFINED_NOT_TESTED","V7":"DEFINED_NOT_CORRELATED","R4_2":"NOT_AUTHORIZED","PROMOTION":"BLOCKED"},"artifacts":[]}
for p in sorted(files):
    manifest['artifacts'].append({"name":p.stem,"path":str(p.relative_to(ROOT)),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
(MANIFESTS/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')

unified={"package":manifest,"contracts":{p.stem:json.loads(p.read_text()) for p in sorted(files)}}
(ROOT/'UNIFIED_EXECUTABLE_CONTRACT_PACKAGE.json').write_text(json.dumps(unified,indent=2,ensure_ascii=False)+'\n')

# README engineering handoff
readme='''# Seat Engineering — Executable Contract Package (R4.1)

This package converts the supplied V5/V6/V7 architecture into machine-readable execution contracts. It does not contain fabricated measurements, FE results, occupant data, sled data, validation outcomes, or R4.2 authorization.

## Authority state
- R4.1 baseline hash (authoritative): `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`
- R4.1 mutation: `FALSE`
- V4: blocked by execution + input gaps
- V5: defined / not validated
- V6: defined / not tested
- V7: defined / not correlated
- R4.2: not authorized
- Promotion: blocked

## Agent separation
- Peter_ChatGPT: architecture authority, gates, evidence authority, authorization, revision control, promotion.
- JON_ChatGPT: development, innovation, architecture implementation, parameterization, contracts, design proposals.
- Manus-Tasks: execution, solver/simulation, validation/testing, evidence generation.

## Core rule
No document may be treated as validation evidence. No missing input may silently become a default. No unsourced threshold may produce PASS. All gates fail closed.

## Current highest blocker
The authoritative state itself says V4 is blocked by execution + input gaps. The package therefore makes V4 closure the controlling dependency for downstream consumption.

## Baseline caveat
The expected R4.1 SHA-256 is recorded exactly as supplied by the authority, but the original R4.1 bytes were not supplied in this conversation, so an independent hash recomputation was not possible here. `BASELINE_INTEGRITY_REPORT.json` records that explicitly.
'''
(ROOT/'README.md').write_text(readme,encoding='utf-8')

# zip
import shutil
zip_path=shutil.make_archive(str(ROOT), 'zip', root_dir=ROOT.parent, base_dir=ROOT.name)
print(zip_path)
print('files',len(list(ROOT.rglob('*'))))
