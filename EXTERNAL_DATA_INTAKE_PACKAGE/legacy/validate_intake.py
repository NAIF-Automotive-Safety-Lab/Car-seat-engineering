#!/usr/bin/env python3
import json, hashlib, sys
from pathlib import Path
BASELINE="fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68"
REQUIRED=['VALUE', 'UNIT', 'SOURCE_ID', 'REVISION', 'DATE', 'SHA256', 'PROVENANCE', 'CONFIGURATION_MAPPING', 'VALIDATION_STATUS']
root=Path(__file__).parent
p=json.loads((root/"EXTERNAL_DATA_INTAKE_PACKAGE.json").read_text())
results=[]
for r in p["RECORDS"]:
 f=root/("source_data_required" if r["CATEGORY"]=="SOURCE_DATA_REQUIRED" else "test_data_required")/r["EXACT_FILENAME"]
 d=json.loads(f.read_text()) if f.exists() else {}
 missing=[x for x in REQUIRED if d.get(x) in (None,"")]
 mapping=str(d.get("CONFIGURATION_MAPPING") or "")
 ok=(not missing and BASELINE in mapping and d.get("VALIDATION_STATUS")=="CLOSED")
 results.append((r["INPUT_ID"],"CLOSED" if ok else "OPEN",missing))
print(json.dumps({"PHYSICAL_INPUT_CLOSURE":"PASS" if all(x[1]=="CLOSED" for x in results) else "BLOCKED","results":results},indent=2))
sys.exit(0 if all(x[1]=="CLOSED" for x in results) else 1)
