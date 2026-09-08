"""EERE: Engineering Evidence & Recovery Engine."""

__version__ = "0.1.0"

from .acquisition import ArtifactAcquisitionService

__all__ = ["ArtifactAcquisitionService"]


def main() -> int:
    from .cli import main as _main
    return _main()

if __name__ == "__main__":
    raise SystemExit(main())


# subpackages are intentionally imported lazily to keep acquisition usable without CAD runtimes.

