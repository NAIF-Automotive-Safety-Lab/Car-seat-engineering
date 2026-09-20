from __future__ import annotations
from dataclasses import dataclass
from pathlib import PurePosixPath
from .core import sha256

ALLOWED = {".step", ".stp", ".stl", ".obj", ".csv", ".json", ".pdf"}
MAGIC = {".pdf": b"%PDF", ".stl": None, ".obj": None, ".step": None, ".stp": None, ".csv": None, ".json": None}
@dataclass(frozen=True)
class IngestedFile: filename: str; size: int; sha256: str; media_type: str; parser_status: str

def validate_upload(filename: str, data: bytes, media_type: str, limit: int = 25 * 1024 * 1024) -> IngestedFile:
    name = PurePosixPath(filename).name
    if name != filename or not name or "\x00" in name: raise ValueError("path traversal rejected")
    suffix = PurePosixPath(name).suffix.lower()
    if suffix not in ALLOWED: raise ValueError("unsupported file format")
    if not data or len(data) > limit: raise ValueError("invalid file size")
    expected = MAGIC[suffix]
    if expected and not data.startswith(expected): raise ValueError("content signature mismatch")
    # Native semantic execution is deliberately not fabricated: available parsers only validate metadata.
    parser = "METADATA_VALIDATED" if suffix in {".json", ".csv", ".pdf"} else "EXTERNAL_RUNTIME_REQUIRED"
    return IngestedFile(name, len(data), sha256(data), media_type, parser)
