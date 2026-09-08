# EERE Remaining Blockers

## BLOCKED

1. The existing SR11 adapter requests `pychrono.cascade`; only the verified `pychrono.core` binding was built, so the adapter's Cascade STEP-loading path remains BLOCKED.
2. CalculiX was not available as an executable in the inspected environment, so FE solver execution remains NOT_AVAILABLE.
3. STEPcode and step-p21 are not available as executable independent parsers; the current forensic result is explicitly a single deterministic entity scanner and not a dual-parser comparison.
4. The physical-input closure remains incomplete. Occupant definition, restraint geometry, material cards, mass/CG/inertia, friction, absorber characterization, vehicle hardpoints, acceptance criteria, coordinate system, and filtering remain open in the repository's engineering registers.
5. The authoritative V7 native package remains absent. Therefore native V7 recovery, released BOM, PMI/GD&T, materials, fasteners, interfaces, and manufacturing authority remain blocked.

## VERIFIED BOUNDARY

The current verified boundary is byte-preserving acquisition, SHA-256 identity, immutable replay, deterministic STEP scanning, OCP B-Rep inventory, deterministic planar/cylindrical surface counts, explicit PMI absence reporting, Project Chrono core runtime smoke, and Gmsh executable detection.

The exact R4.1 payload is verified with SHA-256 `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`, 62 solids, 402 faces, 1,550 edges, and valid B-Rep status.

## NOT_AVAILABLE

Semantic hole/slot/fillet/chamfer recognition, PMI/GD&T extraction, independent STEPcode/step-p21 comparison, full CAD-to-test traceability graph, FreeCAD/CadQuery/Analysis Situs integration, and CalculiX solver execution are not claimed.

## NEXT

Build and verify `pychrono.cascade` only if the required OpenCascade/pythonocc integration can be completed without replacing the verified OCP path. Separately acquire and verify the authoritative V7 native source package and physical test records before attempting any release-gate closure.
