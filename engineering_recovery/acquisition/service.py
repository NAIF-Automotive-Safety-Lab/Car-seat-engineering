from __future__ import annotations

import hashlib
import json
import mimetypes
import shutil
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


TOOL_VERSION = "EERE-0.1.0"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class ArtifactManifest:
    artifact_id: str
    artifact_name: str
    source: str
    source_type: str
    path: str
    file_size_bytes: int
    sha256: str
    media_type: str
    extension: str
    acquisition_timestamp: str
    tool_version: str
    repository_revision: str
    status: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class ArtifactAcquisitionService:
    def __init__(self, repository_root: Path, repository_revision: str = "UNKNOWN") -> None:
        self.repository_root = repository_root.resolve()
        self.repository_revision = repository_revision

    def calculate_sha256(self, path: Path) -> str:
        return sha256_file(path)

    def calculate_file_size(self, path: Path) -> int:
        return path.stat().st_size

    def detect_mime(self, path: Path) -> str:
        return mimetypes.guess_type(path.name)[0] or "application/octet-stream"

    def copy_without_mutation(self, source: Path, destination: Path) -> None:
        source = source.resolve()
        if not source.is_file():
            raise FileNotFoundError(str(source))
        destination.parent.mkdir(parents=True, exist_ok=True)
        before = sha256_file(source)
        if destination.exists():
            if sha256_file(destination) != before or destination.stat().st_size != source.stat().st_size:
                raise RuntimeError("IMMUTABLE_DESTINATION_IDENTITY_CONFLICT")
            return
        shutil.copyfile(source, destination)
        after = sha256_file(destination)
        if before != after:
            raise RuntimeError("COPY_MUTATION_DETECTED")

    def immutable_copy(self, source: Path, destination: Path) -> None:
        self.copy_without_mutation(source, destination)
        destination.chmod(0o444)

    def create_manifest(self, source: Path, artifact_path: Path | None = None) -> ArtifactManifest:
        source = source.resolve()
        if not source.is_file():
            raise FileNotFoundError(str(source))
        path = (artifact_path or source).resolve()
        digest = sha256_file(path)
        return ArtifactManifest(
            artifact_id=digest,
            artifact_name=source.name,
            source=str(source),
            source_type="local_file",
            path=str(path),
            file_size_bytes=path.stat().st_size,
            sha256=digest,
            media_type=self.detect_mime(path),
            extension=path.suffix.lower(),
            acquisition_timestamp=datetime.now(timezone.utc).isoformat(),
            tool_version=TOOL_VERSION,
            repository_revision=self.repository_revision,
            status="VERIFIED",
        )

    def acquire_file(self, source: Path, evidence_root: Path) -> tuple[Path, ArtifactManifest]:
        source = source.resolve()
        evidence_root = evidence_root.resolve()
        destination = evidence_root / source.name
        self.immutable_copy(source, destination)
        manifest = self.create_manifest(source, destination)
        (evidence_root / f"{source.name}.manifest.json").write_text(
            json.dumps(manifest.to_dict(), indent=2) + "\n", encoding="utf-8"
        )
        return destination, manifest

    def verify_identity(self, path: Path, manifest: ArtifactManifest | dict[str, Any]) -> bool:
        expected = manifest.sha256 if isinstance(manifest, ArtifactManifest) else manifest["sha256"]
        return path.is_file() and path.stat().st_size == int(manifest.file_size_bytes if isinstance(manifest, ArtifactManifest) else manifest["file_size_bytes"]) and sha256_file(path) == expected

    def replay_artifact(self, path: Path, manifest: ArtifactManifest | dict[str, Any]) -> dict[str, Any]:
        return {"status": "PASS" if self.verify_identity(path, manifest) else "BLOCKED", "path": str(path), "sha256": sha256_file(path) if path.is_file() else None}

    def record_source(self, source: Path) -> dict[str, str]:
        return {"source": str(source.resolve()), "source_type": "local_file"}

    def record_timestamp(self) -> str:
        return datetime.now(timezone.utc).isoformat()

    def acquire_directory(self, source: Path, evidence_root: Path) -> list[ArtifactManifest]:
        return [self.acquire_file(path, evidence_root)[1] for path in sorted(source.rglob("*")) if path.is_file()]

    def manifest_document(self, manifests: list[ArtifactManifest]) -> dict[str, Any]:
        return {"tool_version": TOOL_VERSION, "artifacts": [item.to_dict() for item in manifests]}


__all__ = ["ArtifactAcquisitionService", "ArtifactManifest", "sha256_file"]


if __name__ == "__main__":
    raise SystemExit("Use ArtifactAcquisitionService from Python; no implicit acquisition is performed.")
