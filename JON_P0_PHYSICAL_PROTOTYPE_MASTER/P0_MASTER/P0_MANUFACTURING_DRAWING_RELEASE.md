# P0 Manufacturing Definition - controlled release content

P0 is a REAL HARDWARE prototype/test article. This package defines the build intent, but it is **not** a substitute for the authoritative native R4.1 CAD/drawing release.

## Manufacturing release rules
1. Use stable V7 reference numerals 100-210 where applicable.
2. Use R4.1/native CAD for all hardpoint coordinates, hole patterns, mating datums, critical edge distances, rail parallelism, hinge axes, and tolerance stacks.
3. Use source concept nominal dimensions only as envelope/setup references, never as pixel-derived hole locations.
4. Every load-bearing fabricated part receives material heat/lot traceability and inspection record.
5. Every replaceable functional item (absorber, lock, rebound unit, bushing/pin set) receives a serial number.
6. Every machined/welded part receives a P0 part revision and inspection record.
7. No production-intent claim is permitted; label the assembly TEST ARTICLE / PROTOTYPE.

## Fabrication package to issue to shop
- Part drawings for all P0-001..P0-016
- Assembly drawing with exploded view and section views
- Datum scheme and GD&T for rails/hinges/clevises
- Weld map and weld procedure reference
- Fastener schedule with certified grades and installation method
- Material certificates
- Inspection plan
- Test article serial plate drawing
- Sensor mounting details and cable routing
- Test fixture attachment drawing

## Release blocker
The uploaded V7 PDF does not expose enough machine-readable native CAD coordinates to honestly manufacture every critical interface from the PDF alone. That is a documentation limitation, not a request to invent geometry.
