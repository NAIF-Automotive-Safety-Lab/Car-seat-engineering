#!/usr/bin/env bash
set -euo pipefail

REPO="https://github.com/wonderwhy-er/DesktopCommanderMCP.git"
COMMIT="a781f5a4b8cfebac6638bc6fcbd38fca6326be53"
DEST_DIR="$(dirname "$0")/upstream-src"

if [ -d "$DEST_DIR" ]; then
  echo "upstream-src already exists at $DEST_DIR"
  echo "To refresh, remove the directory first."
  exit 0
fi

echo "Cloning upstream repository $REPO"
# shallow clone the single commit
git clone --no-checkout "$REPO" "$DEST_DIR"
cd "$DEST_DIR"
# fetch and checkout the specific commit
git fetch --depth 1 origin $COMMIT || git fetch --unshallow || true
git checkout $COMMIT

echo "Upstream source checked out to $DEST_DIR at commit $COMMIT"

echo "Next steps: copy selected files from $DEST_DIR into integrations/desktop-commander/working/ and follow README.md"
