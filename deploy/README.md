# V7-AEGIS Self-Hosted Migration — DigitalOcean + DeepSeek

This directory is the infrastructure bootstrap for the self-hosted control plane migration.

## Evidence boundary
- GitHub `main` remains the canonical engineering source of truth.
- R4.1 remains immutable.
- This migration branch contains infrastructure/bootstrap additions only; it does not alter V7 geometry, design intent, or engineering evidence.
- DeepSeek is an advisory AI provider only. It cannot create evidence, modify baselines, change truth state, or open release gates.
- Real CAE execution must use attested container images and independent hashes.

## Target stack
- Ubuntu 24.04 LTS on DigitalOcean
- Docker Engine + Compose
- PostgreSQL (self-hosted or managed later)
- DeepSeek API (OpenAI-compatible endpoint)
- OpenRadioss + Gmsh + CalculiX execution runtime as a pinned container

## Bootstrap
```bash
cp deploy/.env.example .env
# Fill secrets locally; never commit .env
bash deploy/install.sh
```

The installer is intentionally fail-closed when the runtime image digest is missing.
