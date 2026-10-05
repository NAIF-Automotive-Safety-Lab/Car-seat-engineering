# G2 Automatic Python Physics Fallback Report

## Result

```text
CHRONO
BLOCKED

SELECTED PYTHON ENGINE
MuJoCo

ENGINE QUALIFICATION
PASS

G2-ALT EXECUTION
BLOCKED

G2-ALT REPEATABILITY
PASS

R4.1
UNCHANGED

G3
LOCKED

FINAL STATUS
BLOCKED
```

## Actual runtime evidence

The required fallback order was followed. PyChrono and PyBullet failed local import. MuJoCo 3.14.0 was then provisioned in an isolated Python 3.12 environment through the ordinary package index, compiled a real rigid-body/contact model, advanced 400 deterministic `mj_step` increments, extracted state trajectories and contact-force capability, and executed the identical test twice. The two normalized trajectories, velocities, and contact states match exactly.

This is **G2-ALT-MUJOCO runtime qualification PASS**. It is not a Native Chrono PASS and it is not a T-OCS G2 engineering execution.

## Hard boundary

The selected engine does not repair missing T-OCS authority. The controlled model still lacks an authoritative current body decomposition, joint crosswalk, frames/axes, constraints, limits, lock/rebound law, and initial state. Therefore a T-OCS model cannot be created without inventing data. `G2-ALT EXECUTION` remains BLOCKED and G3 stays LOCKED.

A made-up model that produces a green solver plot would be trash evidence: it would be a MuJoCo demo, not a validated car-seat engineering result.
