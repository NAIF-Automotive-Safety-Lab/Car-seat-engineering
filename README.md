# T-OCS Final GitHub Actions Workflow

Place `.github/workflows/tocs-mbd.yml` in the repository.

The workflow is fail-closed and does not synthesize MBD results. It verifies the frozen R4.1 STEP before expensive solver setup, builds Chrono 10.0.0 with Python + Cascade support because the current prebuilt Conda PyChrono distribution does not provide Cascade, runs a real smoke test, then requires a real SR11 execution manifest and non-empty raw solver output before `V1_STATUS=PASS`.

It does not modify R4.1 and does not create R4.2.
