#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

command -v docker >/dev/null 2>&1 || { echo 'BLOCKED: Docker is not installed.' >&2; exit 20; }
docker compose version >/dev/null 2>&1 || { echo 'BLOCKED: Docker Compose plugin is unavailable.' >&2; exit 21; }

if [[ -z "${AEGIS_DOCKER_IMAGE_DIGEST:-}" ]]; then
  echo 'BLOCKED: AEGIS_DOCKER_IMAGE_DIGEST is not attested.' >&2
  echo 'No solver execution is authorized until a real OCI image sha256 digest is recorded.' >&2
  exit 22
fi

if [[ ! "$AEGIS_DOCKER_IMAGE_DIGEST" =~ ^sha256:[a-f0-9]{64}$ ]]; then
  echo 'BLOCKED: malformed OCI image digest.' >&2
  exit 23
fi

mkdir -p data/runs data/ledger
touch data/runs/.keep data/ledger/.keep

echo 'AEGIS self-host bootstrap prerequisites satisfied.'
echo "Canonical Git SHA: ${AEGIS_CANONICAL_GIT_SHA:-UNSET}"
echo "R4.1 SHA-256: ${AEGIS_R41_SHA256:-UNSET}"
echo "Runtime image digest: ${AEGIS_DOCKER_IMAGE_DIGEST}"

docker compose -f deploy/docker-compose.yml up -d --build

echo 'Services started. Runtime remains evidence-gated.'
