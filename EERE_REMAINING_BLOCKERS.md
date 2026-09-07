# EERE Remaining Blockers

## BLOCKED

1. The approved `pychrono.core` binding is unavailable, so the SR11 dynamic execution path cannot be verified.
2. Gmsh and CalculiX are unavailable, so mesh and FE solver execution cannot be claimed.
3. The physical-input closure remains incomplete. The repository itself records open questions for occupant definition, restraint geometry, material cards, mass/CG/inertia, friction, absorber characterization, vehicle hardpoints, acceptance criteria, coordinate system, and filtering.
4. PMI/GD&T and feature-recognition layers are not yet implemented as deterministic evidence producers.

## NOT_AVAILABLE

STEPcode, step-p21, Analysis Situs, FreeCAD, CadQuery, and BrepMFR are not installed or integrated. They are not represented as verified merely because their upstream repositories are known.

## VERIFIED BOUNDARY

The current verified boundary is byte-preserving artifact acquisition, SHA-256 identity, immutable replay, deterministic STEP entity scanning, and OCP B-Rep inventory. The exact R4.1 payload remains available and verified with SHA-256 `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68` and 62 solids.

## NEXT

Restore an approved, real Chrono binding or provide a controlled build specification. Then run the Chrono smoke test and the existing SR11 adapter without modifying unknown parameters or fabricating dynamics evidence.
