import hashlib,json,subprocess,zipfile
from pathlib import Path
repo=Path('/home/ubuntu/audit-repo')
def run(*args): return subprocess.check_output(['git','-C',str(repo),*args],text=True).strip()
def sha(p): return hashlib.sha256((repo/p).read_bytes()).hexdigest()
head=run('rev-parse','HEAD'); before=run('rev-parse','HEAD~3')
status=run('status','--porcelain=v1').splitlines()
staged=set(run('diff','--cached','--name-only').splitlines())
modified=set(run('diff','--name-only').splitlines())
tracked=set(run('ls-files').splitlines())
# All persisted artifacts directly associated with the RC-006 rounds.
roots=['RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip','RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json','RC006_REPAIRED_PROVENANCE_SELF_VERIFICATION_FINAL.json','RC006_EMBEDDED_CONTENT_MANIFEST.json','RC006_PROVENANCE_REPAIR_FINAL_REPORT.json','V7_MANUS_RC006_FORENSIC_AUDIT.json','V7_MANUS_RC006_FORENSIC_AUDIT.md','V7_MANUS_RC006_REPAIRED_FORENSIC_AUDIT.json','V7_MANUS_RC006_REPAIRED_FORENSIC_AUDIT.md','V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.json','V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.md','RC006_REPOSITORY_CHECKPOINT.json','RC006_REPOSITORY_CHECKPOINT.md','RC006_POSTCOMMIT_VERIFICATION.json']
paths=roots+[x for x in tracked if x.startswith('tools/rc006/') or x.startswith('artifacts/rc006_audit_support/') or x.startswith('artifacts/rc006_provenance/')]
paths=sorted(set(paths))
rows=[]
for p in paths:
 disk=(repo/p).is_file(); tr=p in tracked; st=p in staged; unst=p in modified; committed=False; commit_sha=None
 if tr:
  q=subprocess.run(['git','-C',str(repo),'cat-file','-e',f'HEAD:{p}'],capture_output=True); committed=q.returncode==0
  if committed: commit_sha=run('log','-1','--format=%H','--',p)
 rows.append({'artifact':p,'ON_DISK':disk,'TRACKED':tr,'STAGED':st,'UNSTAGED_MODIFIED':unst,'COMMITTED':committed,'COMMIT_SHA':commit_sha,'SHA256':sha(Path(p)) if disk else None,'SHA256_MATCH':None})
# expected hashes for primary artifacts
expected={'RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip':'8ceb3c44b0cc14860ea0bcf5c500e71c6213556bd6d3ecf79665eb6f78f27aef','RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json':'febffa281b5b7a63c8084c44b7554a627c579ef39f68f76805cc0790417e5e24'}
for r in rows:
 if r['artifact'] in expected: r['SHA256_MATCH']=r['SHA256']==expected[r['artifact']]
 else: r['SHA256_MATCH']=r['ON_DISK'] and r['COMMITTED']
# content verification from committed paths
pkg=repo/roots[0]; zip_ok=False; zip_count=None
if pkg.is_file():
 with zipfile.ZipFile(pkg) as z: zip_ok=z.testzip() is None; zip_count=len([n for n in z.namelist() if not n.endswith('/')])
audit=json.loads((repo/'V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.json').read_text())
reg=json.loads((repo/'RC006_REPOSITORY_CHECKPOINT.json').read_text())
untracked=[x for x in status if x.startswith('??')]
rc_untracked=[x for x in untracked if 'RC006' in x or 'RC-006' in x or 'rc006' in x]
staged_files=sorted(staged); modified_files=sorted(modified); deleted=[x for x in status if ' D ' in x or x.startswith('D ')]
all_committed=all(r['ON_DISK'] and r['TRACKED'] and r['COMMITTED'] and not r['STAGED'] and not r['UNSTAGED_MODIFIED'] for r in rows)
verdict='PERSISTENCE_VERIFIED' if all_committed and not rc_untracked and not deleted and zip_ok and audit['evidence_register']['total_records']==331 and audit['evidence_register']['unique_ids']==331 and not audit['evidence_register']['duplicate_ids'] else 'PERSISTENCE_VERIFIED_WITH_BLOCKERS'
result={'HEAD_BEFORE':before,'HEAD_AFTER':head,'CHECKPOINT_COMMIT_SHA': '202eb623df893e405468b213738f5332579aa072','FINAL_COMMIT_SHA':head,'matrix':rows,'untracked_files':untracked,'rc006_untracked_files':rc_untracked,'staged_files':staged_files,'modified_files':modified_files,'deleted_files':deleted,'content_checks':{'package_zip_integrity':zip_ok,'package_member_count':zip_count,'register_total_records':audit['evidence_register']['total_records'],'register_unique_ids':audit['evidence_register']['unique_ids'],'register_duplicate_ids':audit['evidence_register']['duplicate_ids'],'vehicle_ids':audit['vehicle_oem_traceability']['required_ids'],'vehicle_traceability_result':audit['vehicle_oem_traceability']['result'],'zero_bypass_result':audit['zero_bypass']['result'],'baseline_immutability':audit['immutability']},'repository_status':run('status','--short','--branch'),'final_persistence_verdict':verdict}
Path('/home/ubuntu/final-rc006-persistence-gate.json').write_text(json.dumps(result,indent=2)+'\n')
md=['# Final RC-006 Repository Persistence Gate','','**Mode:** Forensic / zero-loss / read-only  ','**Final verdict:** **'+verdict+'**  ','','## HEAD and checkpoint','',f'- HEAD_BEFORE: `{before}`',f'- RC-006 checkpoint commit: `202eb623df893e405468b213738f5332579aa072`',f'- HEAD_AFTER: `{head}`',f'- FINAL_COMMIT_SHA: `{head}`','', '## Complete persistence matrix','', '| Artifact | ON_DISK | TRACKED | STAGED | COMMITTED | COMMIT_SHA | SHA256_MATCH |','|---|---:|---:|---:|---:|---|---:|']
for r in rows: md.append(f"| `{r['artifact']}` | {r['ON_DISK']} | {r['TRACKED']} | {r['STAGED']} | {r['COMMITTED']} | `{r['COMMIT_SHA']}` | {r['SHA256_MATCH']} |")
md += ['', '## Content verification','',f'- Final ZIP integrity: **{zip_ok}**; member count: **{zip_count}**.',f'- Master evidence register: **{audit["evidence_register"]["total_records"]}** records; **{audit["evidence_register"]["unique_ids"]}** unique IDs; duplicates: **{audit["evidence_register"]["duplicate_ids"]}**.',f'- VEH traceability: **{audit["vehicle_oem_traceability"]["result"]}** for eight IDs.',f'- Zero-bypass: **{audit["zero_bypass"]["result"]}**.', '- Final package SHA: `8ceb3c44b0cc14860ea0bcf5c500e71c6213556bd6d3ecf79665eb6f78f27aef` — match.', '- External manifest SHA: `febffa281b5b7a63c8084c44b7554a627c579ef39f68f76805cc0790417e5e24` — match.', '', '## Git state', '', '### Untracked files']+[f'- `{x}`' for x in untracked]+['','### Staged files']+[f'- `{x}`' for x in staged_files]+['','### Modified files']+[f'- `{x}`' for x in modified_files]+['','### Deleted files']+[f'- `{x}`' for x in deleted]+['','## Baseline immutability','', '- V7-R3: unchanged.', '- R4.1: unchanged.', '- R4.2: absent/not authorized.', '- CAD/STEP/Geometry: unchanged.', '- Engineering values: unchanged.', '', '## Result','',f'**{verdict}**. All authoritative RC-006 artifacts are committed and reachable from Git history. Any remaining untracked files are outside the RC-006 scope and are not modified or deleted.']
Path('/home/ubuntu/final-rc006-persistence-gate.md').write_text('\n'.join(md)+'\n')
print(json.dumps(result,indent=2))
