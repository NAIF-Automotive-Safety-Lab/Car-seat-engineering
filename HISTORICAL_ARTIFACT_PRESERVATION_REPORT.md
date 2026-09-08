# Historical Artifact Preservation Report

**Mode:** Read-only forensic inventory followed by preservation commit  
**Final verdict:** **FULL_REPOSITORY_PERSISTENCE_VERIFIED**  

## Git checkpoint
- HEAD before: `a9203c1a10babfe0b8e3f55917a236db90c783a6`
- Preservation commit SHA: `b44e1b2590dabf27a1ed1021693477606ea9dc7f`
- HEAD after: `b44e1b2590dabf27a1ed1021693477606ea9dc7f`
- Commit message: `preserve: archive historical RC-004/RC-005 evidence`

## Complete artifact inventory

| Artifact | Classification | Size (bytes) | SHA-256 | Tracked | Committed | Reachable |
|---|---|---:|---|---:|---:|---:|
| `V7-RC-004_INTERFACE_DEFINITION_PACKAGE.zip` | `VALID_HISTORICAL_ARTIFACT` | 5899 | `a2545e781afbdfd0a8241bb99a0d2eeb4fd9c3b4576476158e14d6cac4de0e75` | True | True | True |
| `V7_MANUS_INTERFACE_DEFINITION_AUDIT.json` | `VALID_HISTORICAL_ARTIFACT` | 12526 | `a218353818262478225dfb39edfaa321740f0cdb4f927f553d3a51bc8aa0fe91` | True | True | True |
| `V7_MANUS_INTERFACE_DEFINITION_AUDIT.md` | `DERIVED/SUPPORTING_ARTIFACT` | 4823 | `4e963c256d3459ed5c2154b7b12c2f2fe161569eac57962e6e6ff28e346ba009` | True | True | True |
| `V7_MANUS_RC005_FORENSIC_AUDIT.json` | `VALID_HISTORICAL_ARTIFACT` | 10681 | `d2489c241b3a9493cca37151614f3ea7b7c50259d8f311e59fa8bceb1723beec` | True | True | True |
| `V7_MANUS_RC005_FORENSIC_AUDIT.md` | `DERIVED/SUPPORTING_ARTIFACT` | 3749 | `d477b8a89f0fcd2e3ae66143366353a58f0a686856ea15324b7a1070299891e5` | True | True | True |
| `artifacts/interface_semantics_closure/V7_INTERFACE_SEMANTICS_CLOSURE.json` | `VALID_HISTORICAL_ARTIFACT` | 24101 | `29f25f59cbe9ca1aa2fee6b86e39f99967861dd3aedbc739692f1b3052964022` | True | True | True |
| `artifacts/interface_semantics_closure/V7_INTERFACE_SEMANTICS_CLOSURE.md` | `DERIVED/SUPPORTING_ARTIFACT` | 6625 | `d75c824e2ac0ea9479dd00445ef88424b969799061b8949dac2e257b8993d4e2` | True | True | True |

## Duplicate analysis

- No inventory artifact was byte-identical to a pre-existing committed file.
- The RC-004 ZIP has independent package provenance and was preserved as its original historical package.
- The interface closure JSON/MD are supporting artifacts, not replacements for RC-006 authoritative artifacts.
- RC-004 ZIP integrity: PASS; member count: 4.

## Classification summary

- `VALID_HISTORICAL_ARTIFACT`: RC-004 package, RC-004 JSON audit, RC-005 JSON audit, interface closure JSON.
- `DERIVED/SUPPORTING_ARTIFACT`: RC-004 MD audit, RC-005 MD audit, interface closure MD report.
- `TEMPORARY/GENERATED`: none identified.
- `UNKNOWN`: none identified.

## Final repository state

All legitimate historical/project artifacts in the requested scope are tracked, committed, and reachable from Git history. No untracked items remain in the requested preservation scope.

## Engineering protection

- RC-006 was not modified.
- V7-R3 and R4.1 were not modified.
- CAD/STEP/Geometry and engineering values were not modified.
- No R4.2 or RC-007 was created.
- No historical artifact was deleted or replaced.
