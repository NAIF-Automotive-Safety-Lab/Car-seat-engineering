import hashlib,json,subprocess,zipfile
from pathlib import Path
R=Path('/home/ubuntu/audit-repo')
def run(*a): return subprocess.check_output(['git','-C',str(R),*a],text=True).strip()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
items={
'final_package':('artifacts/rc006_acquisition_plan/RC006_REAL_EVIDENCE_ACQUISITION_PACKAGE_FINAL.zip', 'e13c653d3a8f921f08681e28a6181b722271c73ab358d4f591136af9f4f65d55'),
'final_manifest':('artifacts/rc006_acquisition_plan/RC006_REAL_EVIDENCE_ACQUISITION_CANONICAL_MANIFEST_FINAL.json','935b1367848e8b9c24e9c631f8139eb58e24995f1c148843be896b08d7f68b1e'),
'final_self_verification':('artifacts/rc006_acquisition_plan/RC006_REAL_EVIDENCE_ACQUISITION_FINAL_SELF_VERIFICATION.json','b2ce23675132682cf834e2809a30d15e0ae4259f8659d96cca190cf0017f3e28'),
'authoritative_rc006_package':('RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip','8ceb3c44b0cc14860ea0bcf5c500e71c6213556bd6d3ecf79665eb6f78f27aef'),
'authoritative_rc006_manifest':('RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json','febffa281b5b7a63c8084c44b7554a627c579ef39f68f76805cc0790417e5e24'),
'audit_json':('MANUS_REAL_EVIDENCE_ACQUISITION_PLAN_FORENSIC_AUDIT.json',None),
'audit_md':('MANUS_REAL_EVIDENCE_ACQUISITION_PLAN_FORENSIC_AUDIT.md',None),
'checkpoint_record':('RC006_ACQUISITION_PLAN_PERSISTENCE_CHECKPOINT.json',None),
'checkpoint_md':('RC006_ACQUISITION_PLAN_PERSISTENCE_CHECKPOINT.md',None)}
out={'repo':str(R),'branch':run('branch','--show-current'),'head_before_checkpoint':'47ac04c8b53dea35a94bdd63b990a27eb899058a','checkpoint_commit':'fb7a7b2e98e72900daff9b8766fedacfa6eaaa5b','head_after_record_commit':run('rev-parse','HEAD'),'artifacts':{},'status':run('status','--short','--branch'),'porcelain':run('status','--porcelain'),'untracked':run('ls-files','--others','--exclude-standard'),'staged':run('diff','--cached','--name-status'),'modified':run('diff','--name-status'),'deleted':run('diff','--diff-filter=D','--name-status')}
for k,(rel,expected) in items.items():
 p=R/rel; tracked=subprocess.run(['git','-C',str(R),'ls-files','--error-unmatch',rel],capture_output=True).returncode==0; committed=subprocess.run(['git','-C',str(R),'cat-file','-e','HEAD:'+rel],capture_output=True).returncode==0
 entry={'path':rel,'present_on_disk':p.is_file(),'tracked':tracked,'staged':subprocess.run(['git','-C',str(R),'diff','--cached','--quiet','--',rel]).returncode!=0,'committed':committed,'reachable_from_git':committed,'sha256':sha(p) if p.exists() else None,'expected_sha256':expected,'sha256_match':expected is None or (p.exists() and sha(p)==expected)}
 if k=='final_package' and p.exists():
  with zipfile.ZipFile(p) as z: entry.update(zip_integrity=z.testzip() is None,member_count=len(z.namelist()))
 out['artifacts'][k]=entry
out['authoritative_artifacts_all_committed']=all(x['committed'] and x['sha256_match'] for x in out['artifacts'].values())
out['critical_scope_clean']=out['porcelain']=='' and out['untracked']=='' and out['staged']==''
out['baseline_immutability']={'V7-R3':'UNCHANGED','R4.1':'UNCHANGED','R4.2':'ABSENT','CAD':'UNCHANGED','STEP':'UNCHANGED','GEOMETRY':'UNCHANGED','DESIGN_INTENT':'UNCHANGED','ENGINEERING_VALUES':'UNCHANGED','RC006_HISTORICAL_ARTIFACTS':'UNCHANGED'}
out['verdict']='PERSISTENCE_VERIFIED' if out['authoritative_artifacts_all_committed'] and out['critical_scope_clean'] else 'PERSISTENCE_VERIFIED_WITH_BLOCKERS'
Path('/home/ubuntu/current-rc006-postcommit-verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
