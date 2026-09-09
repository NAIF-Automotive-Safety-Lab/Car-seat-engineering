# GitHub RC-006 Persistence Forensic Report

**Final verdict:** `GITHUB_PERSISTENCE_FAILED`  

## Result

- Local reconciled HEAD: `822cb3bcb9448068e092031315104b113b20cbc8`  
- GitHub main at last fetch: `460fda6b46f27f9670d97e8978e32f7cb3047fd4`  
- Merge base: `460fda6b46f27f9670d97e8978e32f7cb3047fd4`  
- Non-destructive merge commit: `822cb3bcb9448068e092031315104b113b20cbc8`  
- No force-push, reset, rebase, squash, deletion, or history rewrite was used.

## Push result

- Command: `git push origin main`
- Return code: `128`
- Error: `fatal: could not read Username for 'https://github.com': terminal prompts disabled`
- Push was blocked because GitHub credentials were unavailable in the environment.

## Remote verification

- Independent GitHub clone confirmed that the authoritative RC-006 ZIP and external canonical manifest are absent from GitHub main `460fda6b46f27f9670d97e8978e32f7cb3047fd4`.
- `RC006_REPOSITORY_CHECKPOINT.json` is present remotely, but it does not contain the authoritative package bytes.

## Local artifact matrix

| Artifact | On disk | Tracked | Committed | Reachable | SHA match |
|---|---:|---:|---:|---:|---:|
| `RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip` | True | True | True | True | True |
| `RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json` | True | True | True | True | True |
| `RC006_EMBEDDED_CONTENT_MANIFEST.json` | True | True | True | True | True |
| `RC006_REPAIRED_PROVENANCE_SELF_VERIFICATION_FINAL.json` | True | True | True | True | True |
| `RC006_PROVENANCE_REPAIR_FINAL_REPORT.json` | True | True | True | True | True |
| `RC006_REPOSITORY_CHECKPOINT.json` | True | True | True | True | True |

## Baseline protection

- V7-R3 unchanged.
- R4.1 unchanged.
- R4.2 absent.
- CAD/STEP/Geometry/Design Intent/Engineering Values unchanged.
- No CAE or physical validation performed.

## Required next action

Enable authenticated GitHub access for this environment, then rerun the ordinary `git push origin main`. Do not use force-push.
