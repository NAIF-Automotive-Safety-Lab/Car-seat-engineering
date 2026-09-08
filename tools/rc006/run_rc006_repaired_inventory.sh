#!/usr/bin/env bash
set -u
UPLOAD=/home/ubuntu/upload
OUT=/home/ubuntu/rc006-repaired-audit-work
rm -rf "$OUT"; mkdir -p "$OUT"
ZIP="$UPLOAD/RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_REPAIRED.zip"
sha256sum "$ZIP" > "$OUT/package.sha256"
unzip -tqq "$ZIP" > "$OUT/zip.integrity.txt" 2>&1
unzip -oq "$ZIP" -d "$OUT/package"
unzip -Z1 "$ZIP" | sort > "$OUT/zip.members.txt"
find "$OUT/package" -type f -exec sha256sum {} \; | sort > "$OUT/member.sha256.txt"
printf '%s\n' '--- package sha ---'; cat "$OUT/package.sha256"
printf '%s\n' '--- members ---'; cat "$OUT/zip.members.txt"
printf '%s\n' '--- integrity ---'; cat "$OUT/zip.integrity.txt"
printf '%s\n' '--- member hashes ---'; cat "$OUT/member.sha256.txt"
