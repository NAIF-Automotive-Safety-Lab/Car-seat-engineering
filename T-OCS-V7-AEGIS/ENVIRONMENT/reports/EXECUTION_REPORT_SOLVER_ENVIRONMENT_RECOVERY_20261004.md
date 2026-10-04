# T-OCS V7.0 REVISION-2 — Solver Environment Recovery Report

## Status

**FINAL STATUS = BLOCKED**

### Step 1 — R4.1 baseline

PASS.

The actual immutable artifact remains:

`/mnt/data/files/T_OCS_V0_R4_1_SHOULDER_REPAIR.step`

Size: **1,220,000 bytes**  
SHA-256: **fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68**  
STEP read: **PASS — IFSelect_RetDone**  
B-Rep: **PASS**  
Solids: **62**  
Mutation: **FALSE**

The gate evaluator now reads `R4_1_BASELINE_LOCK.json → expected_sha256` dynamically.

`GATE-BASELINE-CONFIG-001` = **PASS**.  
The regression changes only the governing lock value in a temporary fixture and confirms the gate follows that field. A source scan also confirms no embedded 64-character R4.1 SHA literal remains in the evaluator.

### Step 2 — Source discovery

No local Git repositories containing Chrono/OpenRadioss/Gmsh/CalculiX source were found. No relevant solver archive or native solver binary was found. Git LFS is unavailable.

Connected repository discovery found:

* **OpenRadioss:** `NAIF-Automotive-Safety-Lab/OpenRadioss`, fork of `gabrielrussell/OpenRadioss`, main @ `bad6688602d770e019b75773d2d5a48b4f9bd918`. This is a real OpenRadioss source fork, but G3 qualification was not started because G2 did not pass.
* **Gmsh:** `NAIF-Automotive-Safety-Lab/pygmsh`, fork of `nschloe/pygmsh`, main @ `5415fb97cb07df9846c7bb178be2fc494f1f4ca3`. This is a Python wrapper/interface, not native Gmsh source.
* **CalculiX:** `NAIF-Automotive-Safety-Lab/CalculiX-Examples`, fork of `calculix/CalculiX-Examples`, master @ `316273e9105e44ce7e3ee05059dac1bc3f256a69`. This is an examples repository, not the solver source.

### Step 3/4 — Project Chrono

Official upstream metadata was verified:

`projectchrono/chrono`  
main @ `9e95a025adbd454f04f263986865b88e38cadd43`

However, the native source tree is not present in the execution host, no local/cache archive exists, and no native Chrono executable/library exists.

Shell GitHub DNS/HTTPS/Git transport remains unavailable.

Therefore:

**CHRONO_SOURCE = BLOCKED**  
**CHRONO_BUILD = BLOCKED**  
**CHRONO_RUNTIME = BLOCKED**

No native Chrono probe was executed because the required source/runtime could not first be provisioned. Probe #1 and Probe #2 are therefore **BLOCKED/NOT_EXECUTED**, not PASS.

### Hard stop

Because native Chrono qualification did not pass:

* G2 authoritative kinematic execution was not started.
* G2 repeatability was not started.
* G3 OpenRadioss qualification was not started.
* G3 Gmsh qualification was not started.
* G3 CalculiX qualification was not started.
* No mesh was generated.
* No physics result was generated.
* G4 engineering closure was not started.

No solver was substituted.

No engineering input was invented.

No R4.1 artifact was modified.
