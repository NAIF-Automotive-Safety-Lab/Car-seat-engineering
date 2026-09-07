# P0 Parameter Count Reconciliation

| Counter | Count |
|---|---:|
| Raw numerical values | 846 |
| True engineering values | 836 |
| Non-engineering numbers | 3 |
| Unresolved raw values | 7 |
| Physical measurements | 0 |
| Test results | 0 |

## Primary classification partition

| Class | Count |
|---|---:|
| Design reference | 441 |
| Design target | 2 |
| Engineering assumption | 196 |
| Source-study value | 179 |
| Unresolved engineering value | 18 |
| **Total engineering values** | **836** |

The prior discrepancy `836 - (441 + 2 + 196 + 179) = 18` is resolved by an explicit `UNRESOLVED_ENGINEERING_VALUE` bucket. The raw partition is `836 + 3 + 7 = 846`. No record is counted in more than one primary class.

`PHYSICAL_MEASUREMENTS = 0` and `TEST_RESULTS = 0` remain unchanged.
