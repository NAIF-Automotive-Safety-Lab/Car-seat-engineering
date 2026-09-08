import hashlib, json, os, subprocess, zipfile
from pathlib import Path
repo=Path('/home/ubuntu/audit-repo')
items=[repo/'V7-RC-004_INTERFACE_DEFINITION_PACKAGE.zip',repo/'V7_MANUS_INTERFACE_DEFINITION_AUDIT.json',repo/'V7_MANUS_INTERFACE_DEFINITION_AUDIT.md',repo/'V7_MANUS_RC005_FORENSIC_AUDIT.json',repo/'V7_MANUS_RC005_FORENSIC_AUDIT.md',repo/'artifacts/interface_semantics_closure']
def sha(p):
 h=hashlib.sha256()
 if p.is_file():
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def files(root): return [p for p in ([root] if root.is_file() else root.rglob('*')) if p.is_file()]
# tracked files and all file hashes for duplicate analysis
tracked=subprocess.check_output(['git','-C',str(repo),'ls-files'],text=True).splitlines()
tracked_hashes={}
for rel in tracked:
 p=repo/rel
 if p.is_file(): tracked_hashes.setdefault(sha(p),[]).append(rel)
report=[]
for item in items:
 fs=files(item)
 rec={'path':str(item.relative_to(repo)),'type':'file' if item.is_file() else 'directory','size_bytes':item.stat().st_size if item.is_file() else sum(p.stat().st_size for p in fs),'mtime':item.stat().st_mtime,'files':[]}
 for p in fs:
  h=sha(p); rel=str(p.relative_to(repo)); dup=tracked_hashes.get(h,[])
  rec['files'].append({'path':rel,'size_bytes':p.stat().st_size,'sha256':h,'mtime':p.stat().st_mtime,'identical_committed_paths':dup})
 if item.suffix=='.zip':
  with zipfile.ZipFile(item) as z:
   rec['zip_test']=z.testzip() is None; rec['member_count']=len([n for n in z.namelist() if not n.endswith('/')]); rec['members']=z.namelist()
 elif item.suffix=='.json':
  try: rec['json_valid']=json.loads(item.read_text()) is not None
  except Exception as e: rec['json_valid']=False; rec['json_error']=str(e)
 report.append(rec)
# keyword/content summaries
for rec in report:
 text=''
 for f in rec['files']:
  p=repo/f['path']
  if p.suffix in {'.json','.md','.txt','.csv'}: text+=p.read_text(errors='replace')+'\n'
 rec['scope_keywords']={k: text.lower().count(k) for k in ['rc-004','rc005','rc-005','interface','closure','forensic','geometry','design intent','physical validation']}
print(json.dumps({'repo':str(repo),'head':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),'items':report},indent=2))
