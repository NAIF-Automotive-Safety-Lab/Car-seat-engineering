from __future__ import annotations

import hashlib
import json
import mimetypes
import os
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

TOOL_VERSION = "eere-acquisition/0.1.0"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


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
    absolute_or_project_relative_path: str
    file_size_bytes: int
    sha256: str
    media_type: str
    extension: str
    acquisition_timestamp: str
    tool_version: str
    repository_revision: str
    status: str

    def as_dict(self) -> dict[str, Any]:
        return self.__dict__.copy()


class ArtifactAcquisitionService:
    """Byte-preserving acquisition service. No artifact is authoritative without a manifest."""

    def __init__(self, project_root: Path, evidence_root: Path | None = None) -> None:
        self.project_root = project_root.resolve()
        self.evidence_root = (evidence_root or self.project_root / "artifacts/engineering-evidence").resolve()

    def calculate_sha256(self, path: Path) -> str:
        return sha256_file(path)

    def calculate_file_size(self, path: Path) -> int:
        return path.stat().st_size

    def detect_mime_type(self, path: Path) -> str:
        return mimetypes.guess_type(path.name)[0] or "application/octet-stream"

    def _revision(self) -> str:
        head = self.project_root / ".git/HEAD"
        return head.read_text(encoding="utf-8").strip() if head.exists() else "UNKNOWN"

    def copy_without_mutation(self, source: Path, destination: Path) -> None:
        if not source.is_file():
            raise FileNotFoundError(str(source))
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if sha256_file(source) != sha256_file(destination) or source.stat().st_size != destination.stat().st_size:
                raise ValueError("IMMUTABLE_COPY_EXISTING_MISMATCH")
            return
        with source.open("rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 20)
        os.chmod(destination, 0o444)
        if sha256_file(source) != sha256_file(destination) or source.stat().st_size != destination.stat().st_size:
            raise ValueError("IMMUTABLE_COPY_MISMATCH")

    def immutable_copy(self, source: Path, artifact_id: str | None = None) -> tuple[Path, ArtifactManifest]:
        source = source.resolve()
        if not source.is_file():
            raise FileNotFoundError(str(source))
        digest = sha256_file(source)
        artifact_id = artifact_id or digest
        destination = self.evidence_root / artifact_id / source.name
        self.copy_without_mutation(source, destination)
        manifest = ArtifactManifest(
            artifact_id=artifact_id,
            artifact_name=source.name,
            source=str(source),
            source_type="local_file",
            absolute_or_project_relative_path=str(destination.relative_to(self.project_root))
            if destination.is_relative_to(self.project_root) else str(destination),
            file_size_bytes=destination.stat().st_size,
            sha256=sha256_file(destination),
            media_type=self.detect_mime_type(source),
            extension=source.suffix.lower(),
            acquisition_timestamp=utc_now(),
            tool_version=TOOL_VERSION,
            repository_revision=self._revision(),
            status="ACQUIRED_IMMUTABLE",
        )
        (destination.parent / "manifest.json").write_text(json.dumps(manifest.as_dict(), indent=2) + "\n", encoding="utf-8")
        return destination, manifest

    def create_manifest(self, path: Path, source: str | None = None, status: str = "VERIFIED") -> ArtifactManifest:
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(str(path))
        digest = sha256_file(path)
        return ArtifactManifest(
            artifact_id=digest,
            artifact_name=path.name,
            source=source or str(path),
            source_type="local_file",
            absolute_or_project_relative_path=str(path.relative_to(self.project_root)) if path.is_relative_to(self.project_root) else str(path),
            file_size_bytes=path.stat().st_size,
            sha256=digest,
            media_type=self.detect_mime_type(path),
            extension=path.suffix.lower(),
            acquisition_timestamp=utc_now(),
            tool_version=TOOL_VERSION,
            repository_revision=self._revision(),
            status=status,
        )

    def verify_identity(self, path: Path, expected_sha256: str, expected_size: int | None = None) -> ArtifactManifest:
        manifest = self.create_manifest(path, status="VERIFIED")
        if manifest.sha256 != expected_sha256:
            raise ValueError(f"SHA256_MISMATCH:{manifest.sha256}")
        if expected_size is not None and manifest.file_size_bytes != expected_size:
            raise ValueError(f"FILE_SIZE_MISMATCH:{manifest.file_size_bytes}")
        return manifest

    def replay_artifact(self, manifest_path: Path) -> ArtifactManifest:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        path = self.project_root / manifest["absolute_or_project_relative_path"]
        return self.verify_identity(path, manifest["sha256"], manifest["file_size_bytes"])

    def acquire_file(self, source: Path, artifact_id: str | None = None) -> dict[str, Any]:
        destination, manifest = self.immutable_copy(source, artifact_id)
        return {"destination": str(destination), "manifest": manifest.as_dict()}

    def acquire_directory(self, source: Path) -> list[dict[str, Any]]:
        return [self.acquire_file(path) for path in sorted(source.rglob("*")) if path.is_file()]

    def record_source(self, source: Path) -> dict[str, Any]:
        return self.create_manifest(source).as_dict()

    def record_timestamp(self) -> str:
        return utc_now()


def verify_manifest_round_trip(manifest: dict[str, Any], root: Path) -> bool:
    path = root / manifest["absolute_or_project_relative_path"]
    return path.is_file() and path.stat().st_size == manifest["file_size_bytes"] and sha256_file(path) == manifest["sha256"]
