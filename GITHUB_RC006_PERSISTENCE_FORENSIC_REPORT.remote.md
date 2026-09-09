# GitHub RC-006 Persistence Forensic Report

## Final verdict

# `GITHUB_PERSISTENCE_FAILED`

The requested local source repository `/home/ubuntu/audit-repo` does not exist in the current environment. The two authoritative RC-006 artifacts were also not found in `/home/ubuntu/repo`, `/home/ubuntu/upload`, or `/tmp`. Therefore no safe reconciliation or push of RC-006 content was possible.

No artifact was fabricated, regenerated, substituted, deleted, or silently copied from another repository.

## Required artifact verification

| Artifact | Expected SHA-256 | On disk | Tracked | Committed | Reachable from GitHub main | SHA match |
|---|---|---:|---:|---:|---:|---:|
| `RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip` | `8ceb3c44b0cc14860ea0bcf5c500e71c6213556bd6d3ecf79665eb6f78f27aef` | No | No | No | No | No |
| `RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json` | `febffa281b5b7a63c8084c44b7554a627c579ef39f68f76805cc0790417e5e24` | No | No | No | No | No |

The following additional RC-006 items were also not found: embedded content manifest, provenance self-verification, provenance repair report, RC-006 audit/checkpoint records, the 331-record dataset, and the eight VEH traceability IDs.

## Repository state actually available

The requested `/home/ubuntu/audit-repo` path was absent, so its local HEAD, remote HEAD, merge base, local-only commits, and remote-only commits could not be computed. The available repository `/home/ubuntu/repo` was inspected separately:

| Field | Result |
|---|---|
| Branch | `main` |
| Local HEAD | `460fda6b46f27f9670d97e8978e32f7cb3047fd4` |
| `origin/main` | `460fda6b46f27f9670d97e8978e32f7cb3047fd4` |
| Merge base | `460fda6b46f27f9670d97e8978e32f7cb3047fd4` |
| Divergence | `0 0` |
| Push performed for RC-006 | No |

The available repository contains unrelated untracked V7 verification files. They were not staged or altered.

## Prohibited-operation checks

No force push, reset, rebase, squash, deletion, CAD/STEP/geometry modification, design-intent modification, engineering-value modification, CAE execution, or physical test was performed.

The existing R4.1 STEP hash was recomputed as:

`fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

No R4.2 artifact was found.

## Required recovery action

Make `/home/ubuntu/audit-repo` available in the current environment, or provide the exact RC-006 artifacts and their authoritative Git source. Then repeat byte-level verification, file-level ancestry checks, safe non-fast-forward reconciliation if required, push with ordinary `git push origin main`, and verify the result from a fresh independent clone.

**STOP CONDITION:** GitHub persistence cannot be claimed.
