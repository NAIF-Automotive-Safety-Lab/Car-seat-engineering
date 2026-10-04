# LOCAL RUNTIME SALVAGE — CASADI

LOCAL RUNTIME SALVAGE
PASS

AVAILABLE ENGINE
CasADi 3.7.2

ENGINE CLASS
SYMBOLIC_NUMERICAL_OPTIMIZATION_FRAMEWORK

QUALIFICATION
PASS — native CasADi Python binding and bundled shared library load; deterministic software probe repeated with identical trajectory (max absolute difference 0.0, tolerance 1e-12).

G2-ALT
BLOCKED

G2-KINEMATIC-COMPUTATIONAL
BLOCKED

REPEATABILITY
PASS for CasADi software runtime qualification; G2/model repeatability NOT_PROVEN because no T-OCS engineering model executed.

CHRONO
BLOCKED

R4.1
UNCHANGED

G3
LOCKED

ROOT CAUSE
CasADi is the first newly available local computational environment, but the authoritative T-OCS G2 model inputs are incomplete. The existing G2 status identifies joint/constraint closure as BLOCKED. The V7 verification record identifies kinematics as DEFINED NOT EXECUTABLE and mass/CG/inertia as null. Creating a G2 case would require inventing engineering inputs, which is prohibited.

NO-GO CLAIMS
Do not claim Native Chrono PASS. Do not claim G2 PASS. Do not claim physical/crash/manufacturing validation. Do not treat the synthetic CasADi software probe as T-OCS engineering evidence.

NEXT SINGLE EXECUTION UNIT
Close the authoritative T-OCS G2 joint/constraint and initial-condition input set from controlled project evidence; no solver execution until that closure is released.

FINAL STATUS
BLOCKED