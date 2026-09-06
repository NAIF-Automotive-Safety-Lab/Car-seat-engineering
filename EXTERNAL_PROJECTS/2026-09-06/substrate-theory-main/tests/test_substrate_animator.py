"""Tests for substrate_animator — animated substrate visualization.

Each animation function is asserted to produce a non-empty file on disk.
ffmpeg is not required: when unavailable, the animator falls back to GIF
via PillowWriter and the requested ``.mp4`` path is rewritten to ``.gif``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest  # noqa: E402

from stiff_medium.substrate_animator import (  # noqa: E402
    animate_all,
    animate_cone_bouncing,
    animate_cosmology_desaturation,
    animate_kink_scattering,
    animate_lattice_2d,
)


# Tighter, faster animations for the test suite — full 5 s renders are
# triggered separately by the bulk driver test.
_FAST_FRAMES = 24
_FAST_FPS = 12


def _assert_file_nonempty(path: str, min_bytes: int = 1024) -> None:
    p = Path(path)
    assert p.exists(), f"Animation file not created: {path}"
    assert p.is_file(), f"Path is not a regular file: {path}"
    size = p.stat().st_size
    assert size >= min_bytes, (
        f"Animation file too small ({size} bytes < {min_bytes}): {path}"
    )


# ---------------------------------------------------------------------------
# Individual animations
# ---------------------------------------------------------------------------

def test_animate_kink_scattering(tmp_path):
    out = tmp_path / "kink_scattering.mp4"
    final = animate_kink_scattering(str(out), fps=_FAST_FPS,
                                    n_frames=_FAST_FRAMES)
    _assert_file_nonempty(final)


def test_animate_cosmology_desaturation(tmp_path):
    out = tmp_path / "cosmology_desaturation.mp4"
    final = animate_cosmology_desaturation(str(out), fps=_FAST_FPS,
                                           n_frames=_FAST_FRAMES)
    _assert_file_nonempty(final)


def test_animate_cone_bouncing(tmp_path):
    out = tmp_path / "cone_bouncing.mp4"
    final = animate_cone_bouncing(str(out), fps=_FAST_FPS,
                                  n_frames=_FAST_FRAMES)
    _assert_file_nonempty(final)


def test_animate_lattice_2d(tmp_path):
    out = tmp_path / "lattice_2d_evolution.mp4"
    final = animate_lattice_2d(str(out), fps=_FAST_FPS,
                               n_frames=_FAST_FRAMES, grid_n=32)
    _assert_file_nonempty(final)


# ---------------------------------------------------------------------------
# Bulk driver
# ---------------------------------------------------------------------------

def test_animate_all(tmp_path):
    results = animate_all(str(tmp_path), fps=_FAST_FPS,
                          n_frames=_FAST_FRAMES)
    expected = {"kink_scattering", "cosmology_desaturation",
                "cone_bouncing", "lattice_2d_evolution"}
    assert set(results) == expected
    for _, path in results.items():
        _assert_file_nonempty(path)


# ---------------------------------------------------------------------------
# Fallback behavior: when ffmpeg is unavailable, .mp4 → .gif
# ---------------------------------------------------------------------------

def test_mp4_request_falls_back_to_gif_when_no_ffmpeg(tmp_path):
    """Asking for .mp4 returns a .gif path when ffmpeg is missing."""
    from matplotlib.animation import FFMpegWriter

    out = tmp_path / "anim.mp4"
    final = animate_kink_scattering(str(out), fps=_FAST_FPS,
                                    n_frames=_FAST_FRAMES)
    final_p = Path(final)
    if FFMpegWriter.isAvailable():
        # If ffmpeg IS present, mp4 is honored
        assert final_p.suffix == ".mp4", final
    else:
        # No ffmpeg → must be promoted to .gif
        assert final_p.suffix == ".gif", final
    _assert_file_nonempty(final)


# ---------------------------------------------------------------------------
# Smoke: the produced files declare animation extension we expect
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("func,name", [
    (animate_kink_scattering, "kink_scattering"),
    (animate_cosmology_desaturation, "cosmology_desaturation"),
    (animate_cone_bouncing, "cone_bouncing"),
    (animate_lattice_2d, "lattice_2d_evolution"),
])
def test_extension_is_mp4_or_gif(tmp_path, func, name):
    out = tmp_path / f"{name}.mp4"
    final = func(str(out), fps=_FAST_FPS, n_frames=_FAST_FRAMES)
    assert Path(final).suffix.lower() in {".mp4", ".gif"}
    _assert_file_nonempty(final)
