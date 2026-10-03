# ENV-CHRONO-RECOVERY-001

## Final status
BLOCKED

The host was re-qualified before any solver execution. The R4.1 immutable artifact remains unchanged and the baseline gate is PASS.

## Environment
Linux 6.18.44, x86_64, AMD EPYC 9V74 80-Core Processor, 3 logical CPUs, 6,091,040 KiB RAM, approximately 31,235,008 KiB free on the root filesystem.

Required build tools are installed and executable: GCC/G++/GFortran 14.2.0, CMake 3.31.6, Ninja 1.12.1, Make 4.4.1, Git 2.47.3.

## Network
Shell DNS resolution for github.com fails. HTTPS access to GitHub fails with host-resolution error. Git transport to the official Project Chrono repository fails for the same reason.

Classification: NETWORK-BLOCKED.

No network bypass was attempted.

## Chrono discovery
No authoritative Chrono source tree, archive, native executable, native library, CMake package configuration, or PyChrono module was found under the scanned local roots.

The approved connected source is projectchrono/chrono, branch main, observed head 9e95a025adbd454f04f263986865b88e38cadd43. It is inspectable through the repository integration but cannot be provisioned into the native execution runtime from the current network-constrained host.

## Qualification decision
SOURCE = BLOCKED
BUILD = BLOCKED
RUNTIME = BLOCKED
ENV-CHRONO-001-A = BLOCKED
ENV-CHRONO-001-B = BLOCKED
CHRONO REPEATABILITY = BLOCKED

No native Chrono probe was executed and no substitute was counted.

## G2/G3 lock
Because native Chrono qualification did not pass:
G2 = BLOCKED
G3 = LOCKED

No T-OCS kinematic execution, mesh generation, OpenRadioss execution, CalculiX execution, or physics result was produced.

## Baseline protection
R4.1: 1,220,000 bytes; SHA-256 fbe6b17cdbf7282a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68; STEP read PASS; B-Rep PASS; 62 solids; mutation FALSE.

GATE-BASELINE-CONFIG-001 = PASS. The evaluator reads R4_1_BASELINE_LOCK.json -> expected_sha256 dynamically and contains no embedded R4.1 digest.

## Required next action
Provide an approved authoritative Project Chrono source tree/archive or precompiled native runtime to the executable environment. Then execute the two required native qualification runs before any G2 model execution.