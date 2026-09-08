import hashlib, json, subprocess, zipfile
from pathlib import Path
repo=Path('/home/ubuntu/audit-repo')
head=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
base=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD~2'],text=True).strip()
critical=['RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip','RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json','RC006_REPAIRED_PROVENANCE_SELF_VERIFICATION_FINAL.json','RC006_EMBEDDED_CONTENT_MANIFEST.json','RC006_PROVENANCE_REPAIR_FINAL_REPORT.json','V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.json','V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.md','V7_MANUS_RC006_REPAIRED_FORENSIC_AUDIT.json','V7_MANUS_RC006_REPAIRED_FORENSIC_AUDIT.md','V7_MANUS_RC006_FORENSIC_AUDIT.json','V7_MANUS_RC006_FORENSIC_AUDIT.md','RC006_REPOSITORY_CHECKPOINT.json','RC006_REPOSITORY_CHECKPOINT.md']
def sha(p): return hashlib.sha256((repo/p).read_bytes()).hexdigest()
tracked=[]; missing=[]
for p in critical:
 r=subprocess.run(['git','-C',str(repo),'cat-file','-e',f'HEAD:{p}'])
 (tracked if r.returncode==0 else missing).append(p)
status=subprocess.check_output(['git','-C',str(repo),'status','--porcelain=v1'],text=True).splitlines()
rc_untracked=[x for x in status if ('RC006' in x or 'RC-006' in x or 'rc006' in x)]
with zipfile.ZipFile(repo/critical[0]) as z: zip_test=z.testzip(); members=[n for n in z.namelist() if not n.endswith('/')]
reg=json.loads((repo/'V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.json').read_text())
result={'repository':str(repo),'branch':subprocess.check_output(['git','-C',str(repo),'branch','--show-current'],text=True).strip(),'head_before_checkpoint':base,'checkpoint_commit':'202eb623df893e405468b213738f5332579aa072','record_commit_and_final_head':head,'critical_tracked_and_reachable':len(tracked)==len(critical),'missing_from_head':missing,'rc006_untracked_or_uncommitted':rc_untracked,'other_untracked_outside_scope':status,'critical_sha256':{p:sha(Path(p)) for p in critical},'final_package_sha256':sha(Path(critical[0])),'final_manifest_sha256':sha(Path(critical[1])),'expected_package_sha256':'8ceb3c44b0cc14860ea0bcf5c500e71c6213556bd6d3ecf79665eb6f78f27aef','expected_manifest_sha256':'febffa281b5b7a63c8084c44b7554a627c579ef39f68f76805cc0790417e5e24','package_sha_match':sha(Path(critical[0]))=='8ceb3c44b0cc14860ea0bcf5c500e71c6213556bd6d3ecf79665eb6f78f27aef','manifest_sha_match':sha(Path(critical[1]))=='febffa281b5b7a63c8084c44b7554a627c579ef39f68f76805cc0790417e5e24','zip_test':zip_test is None,'zip_member_count':len(members),'audit_register_total':reg['evidence_register']['total_records'],'audit_register_unique':reg['evidence_register']['unique_ids'],'audit_register_duplicates':reg['evidence_register']['duplicate_ids'],'baseline_immutability':reg['immutability'],'persistence_verdict':'PERSISTENCE_VERIFIED' if len(tracked)==len(critical) and not rc_untracked and not missing and not status else 'PERSISTENCE_VERIFIED_WITH_BLOCKERS'}
Path('/home/ubuntu/rc006-postcommit-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
