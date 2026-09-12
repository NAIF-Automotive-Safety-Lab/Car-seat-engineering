INTEGRATION_PLAN

1) Audit (performed): license, package manager, dependencies, build/start, entrypoint, env
2) Fetch upstream at exact commit using fetch_upstream.sh
3) Selectively copy components listed in manifest.json into working/
4) Run CI: npm ci && npm run build && npm test
5) Run runtime smoke tests and MCP startup tests in isolated runner
6) Run security negative tests
7) Collect logs and artifact hashes under artifacts/copilot/TASK-20260912-0001/
8) Open PR with test results and artifacts
