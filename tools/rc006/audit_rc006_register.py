import json, hashlib, zipfile, re
from pathlib import Path
work=Path('/home/ubuntu/rc006-audit-work/package')
reg=json.loads((work/'01_MASTER/RC006_MASTER_EVIDENCE_REGISTER.json').read_text())
records=reg['records']
required=['ID','DOMAIN','PARAMETER','CURRENT_VALUE','UNIT','VALUE_CLASS','SOURCE','SOURCE_TYPE','AUTHORITY','TRACEABILITY','APPLICABILITY','STATUS','BLOCKING','NEXT_ACTION']
missing={k:[] for k in required}; ids=[]; nonnull=[]; source_issues=[]; promotions=[]; bad_status=[]
for r in records:
    ids.append(r.get('ID'))
    for k in required:
        if k not in r: missing[k].append(r.get('ID'))
    if r.get('CURRENT_VALUE') is not None:
        nonnull.append(r)
        if not r.get('SOURCE') or not r.get('SOURCE_TYPE') or not r.get('TRACEABILITY'): source_issues.append((r.get('ID'),'nonnull_without_source_traceability'))
    s=str(r.get('SOURCE_TYPE','')).upper(); vc=str(r.get('VALUE_CLASS','')).upper(); st=str(r.get('STATUS','')).upper()
    # promotion patterns: a value claimed physical/measured without test source, or model/geometry derivation promoted to physical.
    if any(t in vc for t in ['PHYSICAL','MEASURED','VALIDATED']) and not any(t in s for t in ['TEST','MEASUREMENT','HARDWARE']): promotions.append((r.get('ID'),vc,s))
    if st in {'PHYSICAL_VALIDATED','VALIDATED','MEASURED'} and not any(t in s for t in ['TEST','MEASUREMENT','HARDWARE']): promotions.append((r.get('ID'),st,s))
    if r.get('BLOCKING') is None: bad_status.append(r.get('ID'))
from collections import Counter
print('record_count',len(records)); print('unique_id_count',len(set(ids))); print('duplicate_ids',[x for x,c in Counter(ids).items() if c>1]); print('missing_fields',{k:v[:20] for k,v in missing.items() if v}); print('nonnull_count',len(nonnull)); print('nonnull_ids',[r['ID'] for r in nonnull][:100]); print('source_issues',source_issues); print('promotion_candidates',promotions); print('missing_blocking_field',bad_status)
print('domain_counts',Counter(r.get('DOMAIN') for r in records)); print('status_counts',Counter(r.get('STATUS') for r in records)); print('value_class_counts',Counter(r.get('VALUE_CLASS') for r in records)); print('source_type_counts',Counter(r.get('SOURCE_TYPE') for r in records)); print('blocking_counts',Counter(r.get('BLOCKING') for r in records))
# Explicit checks for requested numbers and phrases across all records.
for needle in ['180','18','22']:
    hits=[(r['ID'],r['PARAMETER'],r['CURRENT_VALUE'],r['UNIT'],r['VALUE_CLASS'],r['SOURCE_TYPE'],r['STATUS']) for r in records if needle in json.dumps(r)]
    print('needle',needle,'count',len(hits),'hits',hits[:40])
# provenance source availability is reproducible as package file references if source path/hash is present; inspect refs.
prov=[]
for r in records:
    if r.get('SOURCE_TYPE') in ['SOURCE','ENGINEERING_SOURCE','EXISTING_AUTHORITATIVE_SOURCE','REFERENCE','LITERATURE','MODEL','GEOMETRY']:
        prov.append((r['ID'],r.get('SOURCE'),r.get('SOURCE_TYPE'),r.get('TRACEABILITY'),r.get('APPLICABILITY')))
print('provenance_sample',prov[:40])
# domain files concise status
for f in sorted(work.rglob('*.json')):
    if f.name=='RC006_MASTER_EVIDENCE_REGISTER.json': continue
    try: d=json.loads(f.read_text())
    except: continue
    vals=[]
    def walk(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k.upper() in {'STATUS','RESULT_STATUS','EVIDENCE_CLOSURE','FINAL_STATUS','PHYSICAL_VALIDATION','MANUFACTURING_RELEASE','CAE_STATUS','READINESS'}: vals.append((k,v))
                walk(v)
        elif isinstance(x,list):
            for v in x: walk(v)
    walk(d)
    if vals: print('FILE',str(f.relative_to(work)),vals[:80])
