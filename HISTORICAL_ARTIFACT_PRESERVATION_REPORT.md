# Historical Artifact Preservation Report

**Mode:** Strict read-only audit → preservation → commit → re-verification  
**Final verdict:** **FULL_REPOSITORY_PERSISTENCE_VERIFIED**  

## Required summary

| Field | Value |
|---|---:|
| `TOTAL_TARGETS` | 6 |
| `TOTAL_FILES_DISCOVERED` | 7 |
| `VALID_HISTORICAL_ARTIFACTS` | 4 |
| `DUPLICATES` | 0 |
| `DERIVED_SUPPORTING` | 3 |
| `TEMPORARY_GENERATED` | 0 |
| `UNKNOWN` | 0 |
| `PRESERVED` | 7 |
| `NOT_PRESERVED` | 0 |
| `DELETED` | 0 |
| `MODIFIED` | 0 |
| `UNAUTHORIZED_MODIFICATIONS` | 0 |

## Git checkpoint
- HEAD before preservation: `a9203c1a10babfe0b8e3f55917a236db90c783a6`
- Preservation commit SHA: `b44e1b2590dabf27a1ed1021693477606ea9dc7f`
- HEAD at report generation: `4923f182c21c3f0343d67686c7579bfb33757867`
- Final report commit will be a documentation-only follow-up; historical artifact preservation remains anchored at `b44e1b2590dabf27a1ed1021693477606ea9dc7f`.

## Artifact matrix

| Artifact | Classification | Size | SHA-256 | Tracked/Committed | Content changed | Commit SHA |
|---|---|---:|---|---|---:|---|
| `V7-RC-004_INTERFACE_DEFINITION_PACKAGE.zip` | `VALID_HISTORICAL_ARTIFACT` | 5899 | `a2545e781afbdfd0a8241bb99a0d2eeb4fd9c3b4576476158e14d6cac4de0e75` | `TRACKED_COMMITTED_REACHABLE` | `False` | `b44e1b2590dabf27a1ed1021693477606ea9dc7f` |
| `V7_MANUS_INTERFACE_DEFINITION_AUDIT.json` | `VALID_HISTORICAL_ARTIFACT` | 12526 | `a218353818262478225dfb39edfaa321740f0cdb4f927f553d3a51bc8aa0fe91` | `TRACKED_COMMITTED_REACHABLE` | `False` | `b44e1b2590dabf27a1ed1021693477606ea9dc7f` |
| `V7_MANUS_INTERFACE_DEFINITION_AUDIT.md` | `DERIVED_SUPPORTING_ARTIFACT` | 4823 | `4e963c256d3459ed5c2154b7b12c2f2fe161569eac57962e6e6ff28e346ba009` | `TRACKED_COMMITTED_REACHABLE` | `False` | `b44e1b2590dabf27a1ed1021693477606ea9dc7f` |
| `V7_MANUS_RC005_FORENSIC_AUDIT.json` | `VALID_HISTORICAL_ARTIFACT` | 10681 | `d2489c241b3a9493cca37151614f3ea7b7c50259d8f311e59fa8bceb1723beec` | `TRACKED_COMMITTED_REACHABLE` | `False` | `b44e1b2590dabf27a1ed1021693477606ea9dc7f` |
| `V7_MANUS_RC005_FORENSIC_AUDIT.md` | `DERIVED_SUPPORTING_ARTIFACT` | 3749 | `d477b8a89f0fcd2e3ae66143366353a58f0a686856ea15324b7a1070299891e5` | `TRACKED_COMMITTED_REACHABLE` | `False` | `b44e1b2590dabf27a1ed1021693477606ea9dc7f` |
| `artifacts/interface_semantics_closure/V7_INTERFACE_SEMANTICS_CLOSURE.json` | `VALID_HISTORICAL_ARTIFACT` | 24101 | `29f25f59cbe9ca1aa2fee6b86e39f99967861dd3aedbc739692f1b3052964022` | `TRACKED_COMMITTED_REACHABLE` | `False` | `b44e1b2590dabf27a1ed1021693477606ea9dc7f` |
| `artifacts/interface_semantics_closure/V7_INTERFACE_SEMANTICS_CLOSURE.md` | `DERIVED_SUPPORTING_ARTIFACT` | 6625 | `d75c824e2ac0ea9479dd00445ef88424b969799061b8949dac2e257b8993d4e2` | `TRACKED_COMMITTED_REACHABLE` | `False` | `b44e1b2590dabf27a1ed1021693477606ea9dc7f` |

## Duplicate analysis

No byte-identical duplicate was found against pre-existing tracked artifacts. The original historical package and reports were preserved under their original paths. No artifact was deleted.

## Baseline immutability

- V7-R3: unchanged.
- R4.1: unchanged.
- R4.2: absent.
- CAD/STEP/Geometry: unchanged.
- Design intent and engineering values: unchanged.
- RC-006: unchanged.

## Final Git state

No untracked files remain after the preservation commit/report scope. All seven discovered files are tracked, committed, and reachable from Git history.

**FULL_REPOSITORY_PERSISTENCE_VERIFIED**
