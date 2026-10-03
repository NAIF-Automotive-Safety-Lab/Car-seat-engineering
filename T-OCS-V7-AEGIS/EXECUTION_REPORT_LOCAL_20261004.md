# T-OCS V7.0 REVISION-2 — Superseding Execution Report

COMMAND:
Direct host qualification and strict baseline verification.

HOST:
Linux 6.18.44 x86_64; AMD EPYC 9V74 80-Core Processor; 3 logical CPUs; 6,091,040 KiB RAM; 31,263,448 KiB root free; NVIDIA runtime not detected.

DEPENDENCIES:
gcc/g++/gfortran 14.2.0 PASS; CMake 3.31.6 PASS; Ninja 1.12.1 PASS; Make 4.4.1 PASS; Git 2.47.3 PASS; Python 3.13.5 PASS; pip 25.1.1 PASS; pytest 9.0.2 PASS; numpy 2.3.5 PASS; scipy 1.17.0 PASS; pandas 2.2.3 PASS; jsonschema 4.26.0 PASS; pydantic 2.13.4 PASS; matplotlib 3.10.8 PASS; OCP 7.9.3.1 PASS; cadquery 2.8.0 PASS.

MISSING DEPENDENCIES:
Git LFS; native Project Chrono runtime/source; OpenRadioss Starter/Engine; Gmsh; ParaView/pvbatch; CalculiX; Code_Aster; Unity 6; Omniverse.

BASELINE:
BLOCKED.

BASELINE EVIDENCE:
Exact file size = 1,220,000 bytes.
Computed SHA-256 = fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68 (64 hexadecimal characters).
STEP read = PASS, IFSelect_RetDone.
B-Rep transfer roots = 1.
Solid count = 62.

CRITICAL STOP CONDITION:
The SHA string supplied in the work order has 65 hexadecimal characters, so it cannot be a SHA-256 digest. It also differs from the computed valid digest. No correction or substitution was made.

G2:
BLOCKED — native Chrono acquisition/runtime was not attempted beyond the allowed environment qualification boundary because GitHub DNS/HTTPS is unavailable and no cached native Chrono runtime exists. No T-OCS engineering kinematic simulation was executed.

G3:
BLOCKED — mandatory physical inputs and required OpenRadioss/Gmsh/ParaView runtimes are not closed/available. No dynamic engineering simulation was executed.

G4:
NOT_PROVEN — evidence infrastructure exists, but the complete engineering evidence gate cannot close while G2/G3 are blocked.

ENVIRONMENT:
BLOCKED.

TESTS EXECUTED:
Host probe PASS; compiler/build probes PASS; Python venv/import/schema/determinism PASS; exact R4.1 SHA computation PASS; exact STEP read PASS; B-Rep transfer PASS; 62-solid count PASS; GitHub DNS/HTTPS/Git probes BLOCKED; Git LFS probe BLOCKED.

ARTIFACTS:
ENVIRONMENT/environment_manifest.json
ENVIRONMENT/capability_matrix.json
ENVIRONMENT/runtime_probes/*
ENVIRONMENT/versions/toolchain_versions.json
ENVIRONMENT/logs/network.log
ENVIRONMENT/reports/ENVIRONMENT_FINAL_REPORT.md
ENVIRONMENT/hashes/SHA256SUMS.json
VALIDATION/FINAL_GATE_RESULT.json

HASHES:
Environment evidence manifest SHA-256 = f244cc7d9bfc97bcef8e420b0caa1f5f5040ccd8c60cb5d94d5c5af26c817ef4

REGRESSION:
BLOCKED.

BLOCKERS:
1. Malformed governing R4.1 SHA-256 reference (65 hex characters).
2. github.com DNS/HTTPS/Git transport unavailable.
3. Git LFS unavailable.
4. Native Project Chrono unavailable.
5. Required G3 runtimes unavailable.

NO-GO ITEMS:
No G2 engineering simulation.
No G3 dynamic physics.
No physical validation claim.
No mutation/regeneration of R4.1.
No silent correction of the governing SHA.

NEXT SINGLE EXECUTION UNIT:
Resolve the authoritative R4.1 SHA-256 reference. Once explicitly corrected/confirmed, rerun the baseline gate before any native Chrono/OpenRadioss acquisition or engineering simulation.
