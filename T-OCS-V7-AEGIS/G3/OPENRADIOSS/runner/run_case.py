from __future__ import annotations
import json, os, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
closure=ROOT/"G3"/"INPUT_CONTROL"/"input_closure.json"
out=ROOT/"G3"/"OPENRADIOSS"/"outputs"
out.mkdir(parents=True,exist_ok=True)

def blocked(reason):
    payload={"status":"BLOCKED","solver":"OpenRadioss","reason":reason,"execution_performed":False}
    (out/"case_status.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("OPENRADIOSS=BLOCKED REASON="+reason)
    return 20

root_value=os.environ.get("OPENRADIOSS_ROOT","")
if not root_value:
    sys.exit(blocked("OPENRADIOSS_RUNTIME_NOT_DECLARED"))
root=Path(root_value)
if not root.exists():
    sys.exit(blocked("OPENRADIOSS_RUNTIME_PATH_MISSING"))
exe=next((p for p in [root/"starter",root/"starter_linux64",root/"engine_linux64_gf"] if p.is_file() and os.access(p,os.X_OK)),None)
if exe is None:
    sys.exit(blocked("OPENRADIOSS_EXECUTABLE_NOT_FOUND_OR_NOT_EXECUTABLE"))
if not closure.is_file():
    sys.exit(blocked("PHYSICAL_INPUT_CLOSURE_MANIFEST_MISSING"))
data=json.loads(closure.read_text(encoding="utf-8"))
statuses={x.get("status") for x in data.get("mandatory_items",[])}
if data.get("status")!="PASS" or statuses & {"MISSING","UNVERIFIED"}:
    sys.exit(blocked("MANDATORY_PHYSICAL_INPUTS_NOT_CLOSED"))
print("OPENRADIOSS_READY=PASS")
