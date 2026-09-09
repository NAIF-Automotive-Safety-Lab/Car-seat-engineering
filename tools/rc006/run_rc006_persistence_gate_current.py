import hashlib,json,subprocess,os,zipfile
from pathlib import Path
R=Path('/home/ubuntu/audit-repo')
def run(*a): return subprocess.check_output(['git','-C',str(R),*a],text=True).strip()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tracked_state(path):
 p=R/path
 tracked=(subprocess.run(['git','-C',str(R),'ls-files','--error-unmatch',str(path)],capture_output=True).returncode==0) if p.exists() else False
 staged=bool(subprocess.run(['git','-C',str(R),'diff','--cached','--quiet','--',str(path)]).returncode)
 modified=bool(subprocess.run(['git','-C',str(R),'diff','--quiet','--',str(path)]).returncode)
 committed=subprocess.run(['git','-C',str(R),'cat-file','-e','HEAD:'+str(path)],capture_output=True).returncode==0
 return {'path':str(path),'present_on_disk':p.is_file(),'tracked':tracked,'staged':staged,'modified':modified,'committed':committed,'sha256':sha(p) if p.is_file() else None}
critical={
 'final_package':('RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip','8ceb3c44b0cc14860ea0bcf5c500e71c6213556bd6d3ecf79665eb6f78f27aef'),
 'external_manifest':('RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json','febffa281b5b7a63c8084c44b7554a627c579ef39f68f76805cc0790417e5e24'),
 'embedded_manifest':'RC006_EMBEDDED_CONTENT_MANIFEST.json',
 'self_verification':'RC006_REPAIRED_PROVENANCE_SELF_VERIFICATION_FINAL.json',
 'repair_report':'RC006_PROVENANCE_REPAIR_FINAL_REPORT.json',
 'checkpoint_record':'RC006_REPOSITORY_CHECKPOINT.json',
 'postcommit_verification':'RC006_POSTCOMMIT_VERIFICATION.json',
 'acquisition_plan_audit_json':'MANUS_REAL_EVIDENCE_ACQUISITION_PLAN_FORENSIC_AUDIT.json',
 'acquisition_plan_audit_md':'MANUS_REAL_EVIDENCE_ACQUISITION_PLAN_FORENSIC_AUDIT.md',
 'real_evidence_audit_json':'MANUS_REAL_EVIDENCE_ACQUISITION_FORENSIC_AUDIT.json',
 'real_evidence_audit_md':'MANUS_REAL_EVIDENCE_ACQUISITION_FORENSIC_AUDIT.md'}
out={'repo':str(R),'branch':run('branch','--show-current'),'head_before':run('rev-parse','HEAD'),'status_short':run('status','--short','--branch'),'status_porcelain':run('status','--porcelain'),'staged_files':run('diff','--cached','--name-status'),'unstaged_files':run('diff','--name-status'),'untracked_files':run('ls-files','--others','--exclude-standard'),'ignored_rc006':subprocess.check_output(['git','-C',str(R),'status','--ignored','--short'],text=True).splitlines(),'recent_commits':run('log','--oneline','--decorate','-15')}
for k,v in list(critical.items()):
 if isinstance(v,tuple): critical[k]=tracked_state(v[0]); critical[k]['expected_sha256']=v[1]; critical[k]['sha_match']=critical[k]['sha256']==v[1]
 else: critical[k]=tracked_state(v)
out['critical_artifacts']=critical
# validate frozen package against embedded manifest if present
z=R/'RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip'
if z.exists():
 with zipfile.ZipFile(z) as f: out['final_package_validation']={'zip_integrity':f.testzip() is None,'member_count':len(f.namelist()),'names':f.namelist()}
for p in ['RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json','RC006_EMBEDDED_CONTENT_MANIFEST.json','RC006_REPAIRED_PROVENANCE_SELF_VERIFICATION_FINAL.json','RC006_PROVENANCE_REPAIR_FINAL_REPORT.json','RC006_REPOSITORY_CHECKPOINT.json']:
 q=R/p
 if q.exists():
  try: out.setdefault('json_validity',{})[p]=json.loads(q.read_text())
  except Exception as e: out.setdefault('json_validity',{})[p]={'valid':False,'error':str(e)}
Path('/home/ubuntu/current-rc006-persistence-inventory.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
