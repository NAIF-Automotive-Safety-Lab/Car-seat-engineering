from __future__ import annotations
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HASHDIR=ROOT/"G4"/"HASHES"
LEDGERDIR=ROOT/"G4"/"LEDGER"
SELF=ROOT/"G4"/"EVIDENCE_ENGINE"/"self_check_result.json"
HASHFILE=HASHDIR/"SHA256SUMS.json"
LEDGER=LEDGERDIR/"evidence_ledger.jsonl"
HASHDIR.mkdir(parents=True,exist_ok=True); LEDGERDIR.mkdir(parents=True,exist_ok=True)

def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()

rows=[]
for p in sorted(ROOT.rglob("*")):
    if p.is_file() and p not in {SELF,HASHFILE}:
        rows.append({"artifact":str(p.relative_to(ROOT)).replace("\\","/"),"sha256":digest(p),"size_bytes":p.stat().st_size})
HASHFILE.write_text(json.dumps({"schema_version":"1.0","artifacts":rows},indent=2,sort_keys=True)+"\n",encoding="utf-8")
record={"evidence_id":"G4-SELF-001","gate":"G4","claim":"Evidence engine can hash and ledger artifacts","artifact":"G4/EVIDENCE_ENGINE/self_check_result.json","source":"T-OCS-V7-REV2","run_id":os.environ.get("GITHUB_RUN_ID","LOCAL"),"timestamp":datetime.now(timezone.utc).isoformat(),"status":"PASS","hash":digest(HASHFILE),"review_status":"AUTOMATED_SELF_CHECK"}
with LEDGER.open("a",encoding="utf-8") as f: f.write(json.dumps(record,sort_keys=True)+"\n")
SELF.write_text(json.dumps({"schema_version":"1.0","status":"PASS","fail_closed":True,"artifact_count":len(rows),"physical_validation_claimed":False,"engineering_claims_promoted":False},indent=2,sort_keys=True)+"\n",encoding="utf-8")
print("G4_SELF_CHECK=PASS")
