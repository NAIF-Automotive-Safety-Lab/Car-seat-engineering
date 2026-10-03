# T-OCS V7.0 REVISION-2 — Strict Environment Qualification Execution

COMMAND:
Direct host qualification: host/toolchain probes, native C++/Fortran compile and runtime tests, CMake/Ninja build, Python venv/import/schema/determinism tests, R4.1 SHA/read/B-Rep verification, DNS/HTTPS/Git/Git-LFS probes.

HOST:
Linux 6.18.44 x86_64; AMD EPYC 9V74 80-Core Processor; 3 logical CPUs; 6,091,040 KiB RAM; 31,263,448 KiB root free; NVIDIA GPU runtime not detected.

DEPENDENCIES:
gcc 14.2.0 PASS
g++ 14.2.0 PASS
gfortran 14.2.0 PASS
CMake 3.31.6 PASS
Ninja 1.12.1 PASS
Make 4.4.1 PASS
Git 2.47.3 PASS
Python 3.13.5 PASS
pip 25.1.1 PASS
pytest 9.0.2 PASS
numpy 2.3.5 PASS
scipy 1.17.0 PASS
pandas 2.2.3 PASS
jsonschema 4.26.0 PASS
pydantic 2.13.4 PASS
matplotlib 3.10.8 PASS
OCP 7.9.3.1 PASS
Project Chrono native BLOCKED
OpenRadioss Starter BLOCKED
OpenRadioss Engine BLOCKED
Gmsh BLOCKED
ParaView/pvbatch BLOCKED
CalculiX BLOCKED
Code_Aster BLOCKED
Unity 6 BLOCKED
Omniverse BLOCKED
Git LFS BLOCKED

MISSING DEPENDENCIES:
Native Project Chrono runtime/source
OpenRadioss Starter/Engine
Gmsh
ParaView/pvbatch
CalculiX
Code_Aster
Unity 6
Omniverse
Git LFS

BASELINE:
BLOCKED
Reason: the SHA string supplied by the current execution order is 65 hexadecimal characters. A SHA-256 digest must be 64 hexadecimal characters. Computed actual digest is fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68 (64 chars). SHA match therefore FAILS under the exact supplied reference. STEP read PASS; B-Rep transfer roots 1; solid count 62.

G2:
BLOCKED

G3:
BLOCKED

G4:
NOT_PROVEN

ENVIRONMENT:
BLOCKED

TESTS EXECUTED:
Native C++ compile/runtime PASS.
Native Fortran compile/runtime PASS.
CMake/Ninja configure/build PASS.
Python venv dependency/import/schema/determinism PASS.
R4.1 byte SHA calculation PASS.
R4.1 STEP read PASS.
R4.1 B-Rep transfer PASS.
R4.1 solid count = 62 PASS.
github.com DNS BLOCKED.
GitHub HTTPS BLOCKED.
Git transport BLOCKED.
Git LFS BLOCKED.
Native Chrono runtime qualification BLOCKED.
No T-OCS engineering simulation executed.

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
BLOCKED

BLOCKERS:
Malformed required R4.1 SHA reference.
github.com DNS/HTTPS/Git transport unavailable.
Git LFS not installed.
Native Chrono unavailable.
Required G3 runtimes unavailable.

NO-GO ITEMS:
Do not start G2 engineering simulation.
Do not start G3 dynamic physics.
Do not promote historical joint mappings.
Do not alter or regenerate R4.1.
Do not substitute the computed SHA for the malformed governing reference without explicit authority.

NEXT SINGLE EXECUTION UNIT:
Resolve the authoritative R4.1 SHA-256 reference. Confirm whether the governing reference is intended to be the computed 64-character digest or whether another exact byte artifact is authoritative. Then rerun the baseline gate before engine acquisition.
