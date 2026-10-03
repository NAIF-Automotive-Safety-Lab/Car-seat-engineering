# T-OCS V7.0 REVISION-2 — Fail-Closed Execution Report

COMMAND:
aegis gate-status

ENVIRONMENT:
Python 3.13.5; CMake 3.31.6; Ninja 1.12.1; g++ 14.2.0; OCP available; pychrono unavailable; gmsh unavailable; ccx unavailable; ParaView/pvbatch unavailable; no Desktop Commander device; github.com DNS unavailable from the execution container.

ARTIFACTS:
Exact R4.1 STEP; deterministic Open CASCADE readback reports; Chrono capability/build/runtime status; G3 input closure/OpenRadioss status; G4 self-check/ledger; SHA-256 evidence; final gate result.

TESTS:
SHA-256 verification; STEP readback; 62-solid count; repeatability of baseline readback; native Chrono acquisition attempt; G3 fail-closed check; JSON validation; G4 hashing self-check; deterministic gate evaluation.

RESULTS:
R4.1 SHA-256 = fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68.
R4.1 size = 1,220,000 bytes.
STEP readback = IFSelect_RetDone; transferred roots = 1; solid count = 62.
Two baseline-readback reports are byte-identical by SHA-256: fb14843e8c4b35a0f1b8548854c4d3922029d2266fc7ea65c8425111b6edcc73.
Native Project Chrono execution = BLOCKED because github.com DNS failed and pychrono is absent. No T-OCS kinematic scenario executed.
G3 physical-input closure = BLOCKED; OpenRadioss run refused.
G4 self-check = PASS; G4 gate = NOT_PROVEN.

REGRESSION:
BLOCKED — source-readback repeatability PASS; engineering G2/G3 regression is not executable while G2/G3 are blocked.

EVIDENCE:
G0-R41-001
G2-CAP-001
G3-CLOSURE-001
G4-SELF-001

BLOCKERS:
G1 NOT_PROVEN — no separate G1 acceptance contract supplied.
G2 BLOCKED — native Chrono runtime unavailable; authoritative T-OCS joint/constraint executable closure incomplete.
G3 BLOCKED — mandatory physical inputs missing/unverified; OpenRadioss runtime unavailable.
G4 NOT_PROVEN — complete critical-artifact evidence review/hash closure incomplete.

G2:
BLOCKED

G3:
BLOCKED

G4:
NOT_PROVEN

OVERALL:
BLOCKED

NEXT EXECUTION UNIT:
Close the authoritative T-OCS joint/constraint model from the exact SHA-matched R4.1 artifact — frames, axes, DOF, limits, reaction definitions, and LOCK170 semantics — without inference; then rerun native Chrono G2.
