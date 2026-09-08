import json,zipfile,hashlib
from pathlib import Path
w=Path('/home/ubuntu/rc006-audit-work/package'); zpath=Path('/home/ubuntu/upload/RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE.zip')
r=json.loads((w/'01_MASTER/RC006_MASTER_EVIDENCE_REGISTER.json').read_text())['records']
from collections import defaultdict
g=defaultdict(list)
for x in r:g[x['ID']].append(x)
for k,v in g.items():
 if len(v)>1: print('DUPLICATE',k,json.dumps(v,indent=2))
print('NONNULL RECORDS')
for x in r:
 if x['CURRENT_VALUE'] is not None: print(x['ID'],x['DOMAIN'],x['PARAMETER'],repr(x['CURRENT_VALUE']),x['UNIT'],x['VALUE_CLASS'],x['SOURCE_TYPE'],x['STATUS'],x['BLOCKING'])
print('MANIFEST COMPARE')
ext=json.loads(Path('/home/ubuntu/upload/RC006_CANONICAL_MANIFEST.json').read_text())
with zipfile.ZipFile(zpath) as z:
 names=sorted(n for n in z.namelist() if not n.endswith('/'))
 print('zip_names',len(names),'manifest_names',len(ext['members']),'exact',names==sorted(m['path'] for m in ext['members']))
 for m in ext['members']:
  b=z.read(m['path']); a=(len(b),hashlib.sha256(b).hexdigest()); e=(m['size_bytes'],m['sha256'])
  if a!=e: print('MISMATCH',m['path'],a,e)
 print('internal manifest bytes',len(z.read('12_MANIFESTS/RC006_CANONICAL_MANIFEST.json')))
