# V7 Manus Interface Definition Audit — Canonical RC-004 Re-Audit

**Mode:** Strict read-only repair/re-audit gate  
**Final status:** **VERIFIED_WITH_BLOCKERS**  
**Zero-bypass result:** ZIP↔manifest path/count/size/SHA equality verified with zero mismatches.

## Canonical package

| Field | Value |
|---|---|
| Filename | `V7-RC-004_INTERFACE_DEFINITION_PACKAGE.zip` |
| Package SHA-256 | `a2545e781afbdfd0a8241bb99a0d2eeb4fd9c3b4576476158e14d6cac4de0e75` |
| Manifest SHA-256 | `27509720eedc23229662036bcefce5b4499bacc1c83b3f7ff1a63cf3e2979f3b` |
| Member count | 4 |
| Manifest mismatch count | 0 |
| Canonical source | Authoritative original RC-004 artifacts with valid non-empty Markdown |

## Required claims

| Metric | Result |
|---|---:|
| DESIGN_INTENT_9_OF_9 | YES |
| SOURCE_SUPPORTED_COUNT | 2 |
| MODEL_DEFINED_COUNT | 9 |
| INFERRED_COUNT | 7 |
| PHYSICAL_VALIDATED_COUNT | 0 |
| UNSUPPORTED_CLAIMS | 0/FAIL for unsupported numerical physical claims |
| GEOMETRY_CHANGED | NO |
| CAE_READY | NO |
| MANUFACTURING_READY | NO |

## ZIP ↔ manifest verification

| Member | Size | SHA-256 | Result |
|---|---:|---|---|
| `V7_CURRENT_ENGINEERING_INTERFACE_DEFINITION.json` | 15477 | `4b41d70467482a2d7d35d5d33dd393bb346552def3e620e7a65914369a4ff572` | PASS |
| `V7_CURRENT_ENGINEERING_INTERFACE_DEFINITION.md` | 3571 | `b1e73ec95c87b5adfc5d5459b0782f393eb19e91df58b4def00730010df2c84d` | PASS |
| `V7_R4_DESIGN_INTENT_REVISION_MANIFEST.json` | 666 | `abe1ad151d921b6ec46a85cee7679346a2535c56ce5ffc11630ac1397acc1a9d` | PASS |
| `V7_R4_DESIGN_INTENT_SHA256_MANIFEST.json` | 607 | `27509720eedc23229662036bcefce5b4499bacc1c83b3f7ff1a63cf3e2979f3b` | PASS (not self-listed) |

The repaired canonical ZIP contains exactly the four semantic-definition members and no CAD/STEP payload. The source revision manifest states `geometry_changed=false`, `cad_changed=false`, `r4_1_changed=false`, and `r4_2_created=false`.

## Interface matrix

| Interface | Geometry | Model | Design intent | Source support | Physical validation | Final evidence level | Audit result |
|---|---|---|---|---|---|---|---|
| IF-R4-01 (120 ↔ 130L) | GEOMETRY_PROVEN | MODEL_DEFINED | DESIGN_INTENT_INFERRED | PARTIAL / NOT_PROVEN | NO | DESIGN_INTENT_INFERRED | PARTIAL |
| IF-R4-02 (120 ↔ 130R) | GEOMETRY_PROVEN | MODEL_DEFINED | DESIGN_INTENT_INFERRED | PARTIAL / NOT_PROVEN | NO | DESIGN_INTENT_INFERRED | PARTIAL |
| IF-R4-03 (120 ↔ 140) | GEOMETRY_PROVEN | MODEL_DEFINED | DESIGN_INTENT_INFERRED | PARTIAL / NOT_PROVEN | NO | DESIGN_INTENT_INFERRED | PARTIAL |
| IF-R4-04 (120 ↔ 200) | GEOMETRY_PROVEN | MODEL_DEFINED | DESIGN_INTENT_INFERRED | PARTIAL / NOT_PROVEN | NO | DESIGN_INTENT_INFERRED | PARTIAL |
| IF-R4-05 (130L ↔ 140) | GEOMETRY_PROVEN | MODEL_DEFINED | DESIGN_INTENT_INFERRED | PARTIAL / NOT_PROVEN | NO | DESIGN_INTENT_INFERRED | PARTIAL |
| IF-R4-06 (130R ↔ 140) | GEOMETRY_PROVEN | MODEL_DEFINED | DESIGN_INTENT_INFERRED | PARTIAL / NOT_PROVEN | NO | DESIGN_INTENT_INFERRED | PARTIAL |
| IF-R4-07 (140 ↔ 160) | GEOMETRY_PROVEN | MODEL_DEFINED | DESIGN_INTENT_INFERRED | PARTIAL / NOT_PROVEN | NO | DESIGN_INTENT_INFERRED | PARTIAL |
| IF-R4-08 (150L ↔ 160) | GEOMETRY_PROVEN | MODEL_DEFINED | MODEL_DEFINED | SOURCE_SUPPORTED | NO | MODEL_DEFINED | PASS |
| IF-R4-09 (160 ↔ J) | GEOMETRY_PROVEN | MODEL_DEFINED | MODEL_DEFINED | SOURCE_SUPPORTED | NO | MODEL_DEFINED | PASS |

## Evidence-level distinction

- `GEOMETRY_PROVEN`: all nine body pairs, intersection volumes, and face IDs match the authoritative R4 geometry evidence.
- `MODEL_DEFINED`: all nine are represented in the current kinematic/body model; IF-R4-08 and IF-R4-09 have direct joint mappings.
- `DESIGN_INTENT_INFERRED`: IF-R4-01 through IF-R4-07 have stated engineering roles consistent with the model, but no released mechanical source proves their detailed semantics.
- `SOURCE_SUPPORTED`: only IF-R4-08 and IF-R4-09 are counted, and only at model-artifact level for their explicit J_LINK_L/J_SEATBACK relations.
- `PHYSICAL_VALIDATED`: none. No physical validation, CAE, FE, crash, or manufacturing release was performed.

## IF-R4-08 / IF-R4-09

IF-R4-08 is consistent with `J_LINK_L`: parent `150L`, child `160`, q2-coupled link, constrained translations and q2-coupled rotation. IF-R4-09 is consistent with `J_SEATBACK`: parent `J`, child `160`, revolute axis `[0,1,0]`, and q2 limits `-10..25` degrees. These remain non-physical and do not establish pin/hinge clearance or hardware validation.

## Final status

**VERIFIED_WITH_BLOCKERS.** The canonical package is internally integrity-valid with zero mismatches. The engineering content remains design-intent/model definition only; physical validation and manufacturing readiness remain blocked.

Baselines preserved: V7-R3 immutable, R4.1 frozen, R4.2 not created.
