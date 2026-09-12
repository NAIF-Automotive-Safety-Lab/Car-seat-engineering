# Integration helper: DO NOT COMMIT SECRETS

This repository addition is for automation: CI should execute the following in a controlled runner (with no secrets in env):

- cd integrations/desktop-commander
- ./fetch_upstream.sh
- mkdir working && cp -r upstream-src/src working/
- cp upstream-src/package.json working/
- cd working
- npm ci
- npm run build
- npm test

All outputs (logs, artifacts) MUST be collected under artifacts/copilot/TASK-20260912-0001/
