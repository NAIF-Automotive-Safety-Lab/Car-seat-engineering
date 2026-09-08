#!/usr/bin/env bash
set -u
UPLOAD=/home/ubuntu/upload
OUT=/home/ubuntu/rc006-audit-work
rm -rf "$OUT"; mkdir -p "$OUT"
ZIP="$UPLOAD/RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE.zip"
sha256sum "$ZIP" > "$OUT/package.sha256"
unzip -tqq "$ZIP" > "$OUT/zip.integrity.txt" 2>&1
unzip -oq "$ZIP" -d "$OUT/package"
unzip -Z1 "$ZIP" | sort > "$OUT/zip.members.txt"
find "$OUT/package" -type f -exec sha256sum {} \; | sort > "$OUT/member.sha256.txt"
rg -n -i '331|evidence|record|material|density|mass|cg|inertia|fastener|preload|torque|clamp|friction|joint|absorber|18.?22|180|lock|rebound|vehicle|hardpoint|datum|h.?point|anchor|clearance|pulse|initial|ready|derivable|pending|blocking|physical|measured|validated|manufactur|geometry_modified|r4_1_modified|r4_2_created' "$OUT/package" > "$OUT/claim-matches.txt" || true
printf '%s\n' '--- package sha ---'; cat "$OUT/package.sha256"
printf '%s\n' '--- members ---'; cat "$OUT/zip.members.txt"
printf '%s\n' '--- integrity ---'; cat "$OUT/zip.integrity.txt"
printf '%s\n' '--- member hashes ---'; cat "$OUT/member.sha256.txt"
printf '%s\n' '--- claim matches ---'; sed -n '1,800p' "$OUT/claim-matches.txt"
