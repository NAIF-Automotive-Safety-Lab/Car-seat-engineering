#!/usr/bin/env bash
set -u
REPO=/home/ubuntu/audit-repo
OUT=/home/ubuntu/rc006-persistence-work
rm -rf "$OUT"; mkdir -p "$OUT"
cd "$REPO"
{
 echo '=== REPO ==='; pwd; git branch --show-current; git rev-parse HEAD; git status --short --branch
 echo '=== STATUS PORCELAIN V1 ==='; git status --porcelain=v1
 echo '=== STAGED ==='; git diff --cached --name-status
 echo '=== UNSTAGED ==='; git diff --name-status
 echo '=== UNTRACKED ==='; git ls-files --others --exclude-standard
 echo '=== IGNORED RELATED ==='; git status --ignored --short | grep -Ei 'RC-006|RC006|MANUS|audit|interface_semantics' || true
 echo '=== RECENT COMMITS ==='; git log -15 --format='%H%x09%ad%x09%s' --date=iso
 echo '=== BASELINE REFERENCES ==='; git tag --list; git log --all --oneline --decorate -20
} > "$OUT/repo-inventory.txt"
# Authoritative artifacts expected from the final RC-006 provenance gate.
for f in \
  V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.json \
  V7_MANUS_RC006_PROVENANCE_FINAL_FORENSIC_AUDIT.md \
  RC-006_EVIDENCE_PHYSICAL_INPUT_CLOSURE_PACKAGE_PROVENANCE_REPAIRED.zip \
  RC006_REPAIRED_CANONICAL_MANIFEST_FINAL.json \
  RC006_REPAIRED_PROVENANCE_SELF_VERIFICATION_FINAL.json \
  RC006_EMBEDDED_CONTENT_MANIFEST.json \
  RC006_PROVENANCE_REPAIR_FINAL_REPORT.json; do
  printf '%s\n' "$f" >> "$OUT/critical-files.txt"
done
printf '%s\n' '=== CRITICAL LOCAL / REPO PRESENCE ===' > "$OUT/critical-presence.txt"
while IFS= read -r f; do
  printf '%s\t' "$f"; if test -s "$REPO/$f"; then echo 'REPO_PRESENT'; elif test -s "/home/ubuntu/upload/$f"; then echo 'UPLOAD_PRESENT_ONLY'; else echo 'MISSING'; fi
done < "$OUT/critical-files.txt" >> "$OUT/critical-presence.txt"
cat "$OUT/repo-inventory.txt"
cat "$OUT/critical-presence.txt"
