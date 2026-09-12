Integrate DesktopCommanderMCP — selective transplant (metadata only)

What this commit contains:
- UPSTREAM_COMMIT.txt (source commit reference)
- NOTICE_LICENSE.md (attribution & license summary)
- README.md (integration instructions)
- manifest.json (declares what to fetch)
- fetch_upstream.sh (helper to fetch upstream source at exact commit)

What is NOT included in this commit:
- No upstream source files (.git ignored)
- No R4.1 or other immutable V7 artifacts modified
- No secrets or credentials

Intended workflow (manual / CI):
1. Clone this repository and checkout branch: aegis/desktop-commander-integration
2. Run: integrations/desktop-commander/fetch_upstream.sh
   - This will clone the upstream DesktopCommanderMCP at the recorded source commit into integrations/desktop-commander/upstream-src/
3. From the upstream-src, selectively copy the required components into integrations/desktop-commander/working/ as described in manifest.json
4. Run npm install (node >=18) and npm run build in the working directory
5. Run runtime and security tests as per INTEGRATION_TESTS.md

Rationale for this approach:
- Minimizes repository bloat
- Ensures precise provenance (we record the exact upstream commit)
- Allows CI to run installs/builds in an isolated environment without committing compiled artifacts into V7

Warnings:
- Do NOT commit build artifacts into main
- Do NOT modify or regenerate R4.1, STEP assets, or evidence manifest files
