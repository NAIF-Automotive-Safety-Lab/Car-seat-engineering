import json
from pathlib import Path
repo=Path('/home/ubuntu/audit-repo')
j=json.loads((repo/'V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.json').read_text())
md=repo/'V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.md'
s=md.read_text()
old='**Final verdict:** **VERIFIED_WITH_BLOCKERS**  '
new=f'**Final verdict:** **{j["final_decision"]}**  '
s=s.replace(old,new)
old2='**VERIFIED_WITH_BLOCKERS.** Package integrity, manifest architecture, 331-record coverage, VEH traceability, zero-bypass, and immutability all pass. The remaining blockers are genuine unresolved engineering/physical/OEM evidence, not package or provenance defects.'
new2=f'**{j["final_decision"]}.** Package integrity, manifest architecture, 331-record coverage, VEH traceability, zero-bypass, and immutability all pass. The remaining blockers are genuine unresolved engineering/physical/OEM evidence, not package or provenance defects.' if j['final_decision']=='VERIFIED_WITH_BLOCKERS' else '**BLOCKED.** The computed audit decision is blocked; inspect the JSON defect fields for the exact gate failure.'
s=s.replace(old2,new2)
md.write_text(s)
print(j['final_decision'])
