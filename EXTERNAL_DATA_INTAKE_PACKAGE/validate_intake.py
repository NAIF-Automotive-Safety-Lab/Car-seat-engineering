#!/usr/bin/env python3
import copy, hashlib, json, re, sys
from datetime import date
from pathlib import Path

BASELINE = "fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"
REQUIRED = ["VALUE", "UNIT", "SOURCE_ID", "REVISION", "DATE", "SHA256", "PROVENANCE", "CONFIGURATION_MAPPING", "VALIDATION_STATUS"]
EXPECTED = {
    "mass": ("SOURCE_DATA_REQUIRED", "kg"), "CG": ("SOURCE_DATA_REQUIRED", "mm or m in the R4.1 coordinate system"),
    "inertia": ("SOURCE_DATA_REQUIRED", "kg*m^2 or kg*mm^2"), "density": ("SOURCE_DATA_REQUIRED", "kg/m^3"),
    "vehicle_pulse": ("SOURCE_DATA_REQUIRED", "m/s^2 versus s or g versus s"), "initial_velocity": ("SOURCE_DATA_REQUIRED", "m/s"),
    "bolt_stiffness": ("SOURCE_DATA_REQUIRED", "N/m or N/mm"), "friction": ("TEST_DATA_REQUIRED", "dimensionless coefficient"),
    "contact_stiffness": ("TEST_DATA_REQUIRED", "N/m or N/mm"), "contact_damping": ("TEST_DATA_REQUIRED", "N*s/m or explicitly defined equivalent"),
    "restitution": ("TEST_DATA_REQUIRED", "dimensionless"), "absorber_Fx": ("TEST_DATA_REQUIRED", "N versus mm or m"),
    "absorber_Fv": ("TEST_DATA_REQUIRED", "N versus m/s"), "joint_compliance": ("TEST_DATA_REQUIRED", "N/m, N*m/rad, or explicit 6-DOF matrix convention")
}
TARGETS = {
    "mass":"R4.1 assembly mass properties / full assembly configuration", "CG":"R4.1 assembly CG and reference frame",
    "inertia":"R4.1 assembly inertia tensor, reference point and axes", "density":"Every R4.1 FE material region",
    "vehicle_pulse":"R4.1 applicable crash/load case", "initial_velocity":"R4.1 initial-condition/load-case definition",
    "bolt_stiffness":"Each R4.1 represented bolt/joint", "friction":"Each R4.1 contact pair",
    "contact_stiffness":"Each R4.1 normal contact pair", "contact_damping":"Each R4.1 normal contact pair",
    "restitution":"Each relevant R4.1 impact interface", "absorber_Fx":"Each R4.1 absorber instance",
    "absorber_Fv":"Each R4.1 absorber instance", "joint_compliance":"Each R4.1 compliant joint"
}
MARKERS = ("default", "guess", "assum", "fabricat", "synth", "interpolat", "silent", "stale", "tamper", "conflict", "invalid")
HEX64 = re.compile(r"^[0-9a-fA-F]{64}$")

def fail(msg):
    return msg

def validate(root):
    errors=[]
    try: package=json.loads((root/"EXTERNAL_DATA_INTAKE_PACKAGE.json").read_text())
    except Exception as e: return {"PHYSICAL_INPUT_CLOSURE":"BLOCKED","errors":["package JSON/schema parse failure: "+str(e)],"results":[]}
    records=package.get("RECORDS")
    if not isinstance(records,list) or len(records)!=14: errors.append("package must contain exactly 14 records")
    ids=[r.get("INPUT_ID") for r in records] if isinstance(records,list) else []
    if len(ids)!=len(set(ids)): errors.append("duplicate INPUT_ID")
    if set(ids)!=set(EXPECTED): errors.append("input set mismatch")
    results=[]
    for r in records if isinstance(records,list) else []:
        name=r.get("INPUT_ID"); cat=r.get("CATEGORY"); fn=r.get("EXACT_FILENAME")
        local=[]
        if name not in EXPECTED: local.append("unknown INPUT_ID")
        else:
            expcat, expunit=EXPECTED[name]
            if cat != expcat: local.append("category mismatch")
            if r.get("EXPECTED_UNIT") != expunit: local.append("expected unit contract mismatch")
        if not isinstance(fn,str) or Path(fn).name != fn or fn.endswith((".py",".sh")): local.append("unsafe/malformed filename")
        f=root/("source_data_required" if cat=="SOURCE_DATA_REQUIRED" else "test_data_required")/fn if isinstance(fn,str) else root/"missing"
        try: d=json.loads(f.read_text())
        except Exception as e: d={}; local.append("artifact JSON/schema parse failure")
        if not isinstance(d,dict): local.append("artifact must be object")
        else:
            for k in REQUIRED:
                if k not in d or d[k] is None or d[k]=="": local.append("missing "+k)
            if d.get("INPUT_ID") != name: local.append("INPUT_ID mismatch")
            if d.get("CATEGORY") != cat: local.append("CATEGORY mismatch")
            if d.get("VALIDATION_STATUS") != "CLOSED": local.append("not explicitly CLOSED")
            if not isinstance(d.get("VALUE"),(int,float,list,dict)) or isinstance(d.get("VALUE"),bool): local.append("VALUE type/absence invalid")
            if not isinstance(d.get("UNIT"),str) or not d.get("UNIT").strip(): local.append("UNIT invalid")
            if name in EXPECTED and d.get("UNIT") != EXPECTED[name][1]: local.append("UNIT mismatch")
            if not isinstance(d.get("SOURCE_ID"),str) or not d.get("SOURCE_ID").strip(): local.append("SOURCE_ID invalid")
            if not isinstance(d.get("REVISION"),str) or not d.get("REVISION").strip(): local.append("REVISION invalid")
            try: date.fromisoformat(str(d.get("DATE")))
            except Exception: local.append("DATE invalid")
            if not isinstance(d.get("SHA256"),str) or not HEX64.fullmatch(d.get("SHA256","")): local.append("SHA256 invalid")
            if not isinstance(d.get("PROVENANCE"),(dict,list,str)) or not d.get("PROVENANCE"): local.append("PROVENANCE invalid")
            mapping=d.get("CONFIGURATION_MAPPING")
            if not isinstance(mapping,dict) or mapping.get("R4_1_SHA256") != BASELINE or mapping.get("TARGET") != TARGETS.get(name): local.append("wrong/missing R4.1 or part/interface mapping")
            for k,v in d.items():
                if isinstance(v,str) and any(m in v.lower() for m in MARKERS): local.append("forbidden assumption/tamper marker in "+k)
            if "ARTIFACT_PATH" in d and "ARTIFACT_SHA256" in d:
                ap=root/d["ARTIFACT_PATH"]
                if not ap.exists() or hashlib.sha256(ap.read_bytes()).hexdigest()!=d["ARTIFACT_SHA256"]: local.append("artifact hash mismatch")
            if "RAW_DATA_PATH" in d or "REPORT_PATH" in d:
                for pathkey,hashkey in (("RAW_DATA_PATH","RAW_DATA_SHA256"),("REPORT_PATH","REPORT_SHA256")):
                    if pathkey not in d or hashkey not in d: local.append("raw/report hash pair incomplete")
                    elif not (root/d[pathkey]).exists() or hashlib.sha256((root/d[pathkey]).read_bytes()).hexdigest()!=d[hashkey]: local.append("raw/report hash mismatch")
            if "REVISION_STATUS" in d and d["REVISION_STATUS"] != "CURRENT": local.append("stale/non-current revision")
        results.append({"INPUT_ID":name,"STATUS":"REJECTED" if local else "CLOSED","errors":sorted(set(local))})
    if any(x["STATUS"]=="CLOSED" for x in results) and any(x["STATUS"]=="REJECTED" for x in results): errors.append("partial closure is not permitted in package gate")
    if errors: results.append({"INPUT_ID":"PACKAGE","STATUS":"REJECTED","errors":sorted(set(errors))})
    passed=not errors and len(results)==14 and all(x["STATUS"]=="CLOSED" for x in results)
    return {"PHYSICAL_INPUT_CLOSURE":"PASS" if passed else "BLOCKED","results":results,"R4_1_SHA256":BASELINE,"AUDIT_POLICY":"fail-closed; no coercion/default/interpolation/synthesis/automatic approval"}

if __name__ == "__main__":
    out=validate(Path(__file__).parent)
    print(json.dumps(out,indent=2,sort_keys=True))
    sys.exit(0 if out["PHYSICAL_INPUT_CLOSURE"]=="PASS" else 1)
