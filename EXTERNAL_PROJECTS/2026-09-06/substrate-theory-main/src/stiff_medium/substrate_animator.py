"""Substrate Animator — animated visualization sequences of substrate dynamics.

Produces short MP4 (preferred when ffmpeg is available) or GIF fallback
animations for four canonical B3-substrate scenarios:

  1. animate_kink_scattering  — sine-Gordon kink-antikink collision
  2. animate_cosmology_desaturation — σ(t) sweep from 1/2 (saturated) → 0
  3. animate_cone_bouncing — substrate strain bouncing inside σ = ±1/2 cone
  4. animate_lattice_2d — 2D Gaussian bump radiating outward

A convenience driver `animate_all(output_dir)` writes every animation in
sequence using the requested filenames (with extension auto-promoted to
.gif when ffmpeg is absent).

Design notes
------------
- Uses matplotlib `FuncAnimation`.
- Writer selection: prefer FFMpegWriter (mp4) when available; otherwise
  fall back to PillowWriter (gif). When the caller asked for `.mp4` and
  only Pillow is available, the path is rewritten to `.gif` so the file
  is still produced and can be opened.
- Backend is forced to `Agg` so the module is safe in headless test
  environments.
- All four animations target ~5 seconds at 30 fps (150 frames).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter  # noqa: E402


# ---------------------------------------------------------------------------
# Writer selection helper
# ---------------------------------------------------------------------------

def _resolve_writer_and_path(output_path: str, fps: int = 30
                             ) -> Tuple[object, str]:
    """Pick (writer, final_path) given the requested output path.

    Strategy:
      - If ffmpeg is available *and* the user asked for .mp4 → use FFMpegWriter.
      - Otherwise → use PillowWriter and force the extension to .gif.

    Returns
    -------
    (writer, final_path)
        ``writer`` is an instantiated matplotlib writer; ``final_path`` is
        the path that will actually be written (extension may differ from
        the input).
    """
    p = Path(output_path)
    ext = p.suffix.lower()
    want_mp4 = ext == ".mp4"

    if want_mp4 and FFMpegWriter.isAvailable():
        return FFMpegWriter(fps=fps, bitrate=1800), str(p)

    # Fall back to GIF — promote extension if needed.
    if ext != ".gif":
        p = p.with_suffix(".gif")
    return PillowWriter(fps=fps), str(p)


def _ensure_parent(path: str) -> None:
    """Make sure the parent directory of ``path`` exists."""
    parent = Path(path).parent
    if parent and not parent.exists():
        parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Kink-antikink scattering (sine-Gordon)
# ---------------------------------------------------------------------------

def _sine_gordon_kink_pair(x: np.ndarray, t: float, v: float = 0.4,
                           xi: float = 1.0, x0: float = 12.0) -> np.ndarray:
    """Analytic kink (left) + antikink (right) superposition.

    Each soliton is a Lorentz-boosted sine-Gordon kink:
        φ_K(x,t)  = 4 arctan(exp(+γ (x + x0 − v t) / ξ))   (moving right)
        φ_A(x,t) = 4 arctan(exp(−γ (x − x0 + v t) / ξ))   (moving left)
    Their sum approximates a kink-antikink pair colliding head-on near t≈x0/v.
    """
    gamma = 1.0 / np.sqrt(max(1.0 - v * v, 1e-12))
    arg_k = np.clip(+gamma * (x + x0 - v * t) / xi, -50.0, 50.0)
    arg_a = np.clip(-gamma * (x - x0 + v * t) / xi, -50.0, 50.0)
    return 4.0 * np.arctan(np.exp(arg_k)) - 4.0 * np.arctan(np.exp(arg_a))


def animate_kink_scattering(output_path: str, *, fps: int = 30,
                            n_frames: int = 150) -> str:
    """Animate a sine-Gordon kink-antikink collision.

    The two solitons approach head-on, scatter elastically (sine-Gordon is
    integrable), and emerge with their original profile up to a phase shift.

    Parameters
    ----------
    output_path : str
        Where to write the animation. ``.mp4`` is preferred; falls back to
        ``.gif`` when ffmpeg is not installed.
    fps : int
        Frames per second (default 30).
    n_frames : int
        Total number of frames (default 150 → 5 s @ 30 fps).

    Returns
    -------
    str
        The actual path written (extension may differ from ``output_path``
        if the writer fell back to GIF).
    """
    _ensure_parent(output_path)
    writer, final_path = _resolve_writer_and_path(output_path, fps=fps)

    x = np.linspace(-25.0, 25.0, 600)
    v = 0.45
    x0 = 12.0
    t_min, t_max = 0.0, 2.0 * x0 / v + 6.0  # before, during, after collision
    times = np.linspace(t_min, t_max, n_frames)

    fig, ax = plt.subplots(figsize=(7.5, 4.0), dpi=100)
    line_phi, = ax.plot([], [], color="#1f77b4", lw=2.0, label=r"$\phi(x,t)$")
    line_dphi, = ax.plot([], [], color="#d62728", lw=1.0, alpha=0.6,
                          label=r"$\partial_x\phi$")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(-2.0 * np.pi - 1.0, 2.0 * np.pi + 1.0)
    ax.axhline(0.0, color="gray", lw=0.5, alpha=0.5)
    ax.set_xlabel("x / ξ")
    ax.set_ylabel("φ")
    ax.set_title("Sine-Gordon: kink-antikink scattering")
    ax.legend(loc="upper right", fontsize=9)
    time_txt = ax.text(0.02, 0.95, "", transform=ax.transAxes,
                       fontsize=10, va="top")

    def init():
        line_phi.set_data([], [])
        line_dphi.set_data([], [])
        time_txt.set_text("")
        return line_phi, line_dphi, time_txt

    def update(frame: int):
        t = times[frame]
        phi = _sine_gordon_kink_pair(x, t, v=v, x0=x0)
        dphi = np.gradient(phi, x)
        line_phi.set_data(x, phi)
        line_dphi.set_data(x, dphi)
        time_txt.set_text(f"t = {t:5.2f}   v = {v:.2f}c")
        return line_phi, line_dphi, time_txt

    anim = FuncAnimation(fig, update, init_func=init, frames=n_frames,
                         interval=1000.0 / fps, blit=True)
    anim.save(final_path, writer=writer)
    plt.close(fig)
    return final_path


# ---------------------------------------------------------------------------
# 2. Cosmological de-saturation: σ(t) from 1/2 → 0
# ---------------------------------------------------------------------------

def animate_cosmology_desaturation(output_path: str, *, fps: int = 30,
                                   n_frames: int = 150) -> str:
    """Animate substrate de-saturation σ(t): 1/2 → 0.

    The left panel shows the global saturation level σ(t) on a stiff cap
    σ_max = 1/2; the right panel shows the corresponding (toy) Hubble rate
    H(t) ∝ −σ̇(t)/σ(t). The desaturation profile is a smooth tanh sweep
    chosen for visual clarity, not for a quantitative cosmology fit.
    """
    _ensure_parent(output_path)
    writer, final_path = _resolve_writer_and_path(output_path, fps=fps)

    n = 400
    t = np.linspace(0.0, 1.0, n)
    sigma_full = 0.25 * (1.0 - np.tanh(8.0 * (t - 0.5))) + 1e-3
    # Hubble proxy: H = − σ̇ / σ
    sigma_dot = np.gradient(sigma_full, t)
    H_full = -sigma_dot / np.clip(sigma_full, 1e-3, None)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), dpi=100)
    ax_s, ax_h = axes

    ax_s.set_xlim(0.0, 1.0)
    ax_s.set_ylim(-0.05, 0.55)
    ax_s.axhline(0.5, color="gray", ls="--", lw=0.7, label=r"$\sigma_{\max}=1/2$")
    ax_s.axhline(0.0, color="gray", ls="-", lw=0.5)
    ax_s.set_xlabel("t / t_today")
    ax_s.set_ylabel(r"$\sigma(t)$")
    ax_s.set_title("Substrate saturation")
    line_s, = ax_s.plot([], [], color="#2ca02c", lw=2.2)
    dot_s, = ax_s.plot([], [], "o", color="#2ca02c", ms=6)
    ax_s.legend(loc="upper right", fontsize=9)

    ax_h.set_xlim(0.0, 1.0)
    ax_h.set_ylim(-0.5, max(H_full.max(), 1.0) * 1.1)
    ax_h.axhline(0.0, color="gray", lw=0.5)
    ax_h.set_xlabel("t / t_today")
    ax_h.set_ylabel(r"$H(t) \propto -\dot\sigma/\sigma$")
    ax_h.set_title("Effective Hubble rate")
    line_h, = ax_h.plot([], [], color="#9467bd", lw=2.2)
    dot_h, = ax_h.plot([], [], "o", color="#9467bd", ms=6)

    info_txt = fig.text(0.5, 0.96, "", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    def init():
        line_s.set_data([], [])
        dot_s.set_data([], [])
        line_h.set_data([], [])
        dot_h.set_data([], [])
        info_txt.set_text("")
        return line_s, dot_s, line_h, dot_h, info_txt

    def update(frame: int):
        # Map frame to a fractional progress through t.
        idx = int(np.clip(frame * (n - 1) / max(n_frames - 1, 1), 0, n - 1))
        line_s.set_data(t[: idx + 1], sigma_full[: idx + 1])
        dot_s.set_data([t[idx]], [sigma_full[idx]])
        line_h.set_data(t[: idx + 1], H_full[: idx + 1])
        dot_h.set_data([t[idx]], [H_full[idx]])
        info_txt.set_text(
            f"σ(t) = {sigma_full[idx]:.3f}    H(t) = {H_full[idx]:+.3f}"
        )
        return line_s, dot_s, line_h, dot_h, info_txt

    anim = FuncAnimation(fig, update, init_func=init, frames=n_frames,
                         interval=1000.0 / fps, blit=False)
    anim.save(final_path, writer=writer)
    plt.close(fig)
    return final_path


# ---------------------------------------------------------------------------
# 3. Cone bouncing: substrate strain inside |σ| ≤ 1/2
# ---------------------------------------------------------------------------

def animate_cone_bouncing(output_path: str, *, fps: int = 30,
                          n_frames: int = 150) -> str:
    """Animate a substrate strain bouncing inside the σ = ±1/2 cone.

    The cone walls lie at σ = ±1/2; the strain σ(t) executes a damped
    elastic oscillation, reflecting at each wall. This is a simple
    1-D illustration of the cone-bouncing protocol used for the lepton
    mass tower.
    """
    _ensure_parent(output_path)
    writer, final_path = _resolve_writer_and_path(output_path, fps=fps)

    sigma_max = 0.5
    omega = 2.0 * np.pi * 1.5  # 1.5 oscillations per t-unit
    gamma = 0.05
    duration = 5.0
    times = np.linspace(0.0, duration, n_frames)
    # Damped sinusoid clipped at the cone walls (reflects via abs+sign trick).
    raw = sigma_max * np.exp(-gamma * times) * np.sin(omega * times)

    # Build cone-wall geometry for the side panel: a 2D cone in (t, σ).
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=100)
    ax_t, ax_c = axes

    # Time-trace panel
    ax_t.set_xlim(0.0, duration)
    ax_t.set_ylim(-sigma_max - 0.1, sigma_max + 0.1)
    ax_t.axhline(+sigma_max, color="red", ls="--", lw=1.0, label=r"$+\sigma_{\max}$")
    ax_t.axhline(-sigma_max, color="red", ls="--", lw=1.0, label=r"$-\sigma_{\max}$")
    ax_t.axhline(0.0, color="gray", lw=0.5)
    ax_t.set_xlabel("t")
    ax_t.set_ylabel(r"$\sigma(t)$")
    ax_t.set_title("Strain trace (cone-bounded)")
    ax_t.legend(loc="upper right", fontsize=9)
    line_t, = ax_t.plot([], [], color="#1f77b4", lw=2.0)
    dot_t, = ax_t.plot([], [], "o", color="#1f77b4", ms=8)

    # Cone-geometry panel: cone opens with half-angle 45° at origin.
    # Walls: σ = ±r where r is the radial coord from the apex.
    r_ax = np.linspace(0.0, 1.2, 100)
    ax_c.plot(r_ax, +sigma_max * np.ones_like(r_ax), color="red", lw=1.5)
    ax_c.plot(r_ax, -sigma_max * np.ones_like(r_ax), color="red", lw=1.5)
    ax_c.fill_between(r_ax, -sigma_max, +sigma_max, color="#fff5e0", alpha=0.7)
    ax_c.set_xlim(0.0, 1.2)
    ax_c.set_ylim(-0.7, 0.7)
    ax_c.set_xlabel("r (radial)")
    ax_c.set_ylabel(r"$\sigma$")
    ax_c.set_title("Cone interior (|σ| ≤ 1/2)")
    pt_c, = ax_c.plot([], [], "o", color="#1f77b4", ms=10, zorder=5)
    trail_c, = ax_c.plot([], [], "-", color="#1f77b4", lw=1.0, alpha=0.4)

    info_txt = fig.text(0.5, 0.96, "", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    # Recent trail in the cone view: derive a synthetic radial coordinate
    # that grows monotonically with time so the bounce is visible.
    r_traj = 0.05 + 1.0 * times / duration

    def init():
        line_t.set_data([], [])
        dot_t.set_data([], [])
        pt_c.set_data([], [])
        trail_c.set_data([], [])
        info_txt.set_text("")
        return line_t, dot_t, pt_c, trail_c, info_txt

    def update(frame: int):
        line_t.set_data(times[: frame + 1], raw[: frame + 1])
        dot_t.set_data([times[frame]], [raw[frame]])
        pt_c.set_data([r_traj[frame]], [raw[frame]])
        # Trail of recent positions in cone view
        lo = max(0, frame - 30)
        trail_c.set_data(r_traj[lo: frame + 1], raw[lo: frame + 1])
        info_txt.set_text(
            f"t = {times[frame]:.2f}    σ = {raw[frame]:+.3f}   "
            f"|σ|/σ_max = {abs(raw[frame]) / sigma_max:.2f}"
        )
        return line_t, dot_t, pt_c, trail_c, info_txt

    anim = FuncAnimation(fig, update, init_func=init, frames=n_frames,
                         interval=1000.0 / fps, blit=False)
    anim.save(final_path, writer=writer)
    plt.close(fig)
    return final_path


# ---------------------------------------------------------------------------
# 4. 2D lattice: Gaussian bump radiating outward
# ---------------------------------------------------------------------------

def animate_lattice_2d(output_path: str, *, fps: int = 30,
                       n_frames: int = 150, grid_n: int = 64) -> str:
    """Animate a 2D Gaussian bump radiating outward as a circular wave.

    Initial condition: u(x,y,0) = exp(−r²/2σ²) at origin. Time evolution
    is illustrated by an outward-propagating radial wave + amplitude
    decay. This is a *visual* propagation; for the full Lagrangian solver
    see :mod:`stiff_medium.lattice_substrate_2d`.
    """
    _ensure_parent(output_path)
    writer, final_path = _resolve_writer_and_path(output_path, fps=fps)

    L = 6.0
    x = np.linspace(-L, L, grid_n)
    y = np.linspace(-L, L, grid_n)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    rr = np.sqrt(xx ** 2 + yy ** 2)

    duration = 5.0
    times = np.linspace(0.0, duration, n_frames)
    sigma0 = 0.6
    c_wave = 1.0  # propagation speed (units of grid)

    fig, ax = plt.subplots(figsize=(5.5, 5.0), dpi=100)
    # Pre-render a 2D field for frame 0 to set color limits.
    field0 = np.exp(-rr ** 2 / (2.0 * sigma0 ** 2))
    img = ax.imshow(field0, extent=(-L, L, -L, L), origin="lower",
                    cmap="RdBu_r", vmin=-1.0, vmax=1.0, animated=True)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("2D substrate: Gaussian bump radiating")
    cbar = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("u(x,y,t)")
    time_txt = ax.text(0.02, 0.96, "", transform=ax.transAxes,
                       fontsize=10, va="top",
                       bbox=dict(boxstyle="round", facecolor="white",
                                 edgecolor="gray", alpha=0.8))

    def init():
        img.set_data(field0)
        time_txt.set_text("t = 0.00")
        return img, time_txt

    def update(frame: int):
        t = times[frame]
        # Outgoing circular wave packet: amplitude decays as 1/√(1+r),
        # radial profile is a Gaussian shell at r = c·t with width sigma0.
        shell = np.exp(-((rr - c_wave * t) ** 2) / (2.0 * sigma0 ** 2))
        # Static residual bump near origin damps with time.
        bump = np.exp(-rr ** 2 / (2.0 * sigma0 ** 2)) * np.exp(-2.0 * t)
        amp = 1.0 / np.sqrt(1.0 + rr)
        field = bump + shell * amp * np.cos(2.0 * np.pi * (rr - c_wave * t))
        img.set_data(field)
        time_txt.set_text(f"t = {t:4.2f}")
        return img, time_txt

    anim = FuncAnimation(fig, update, init_func=init, frames=n_frames,
                         interval=1000.0 / fps, blit=False)
    anim.save(final_path, writer=writer)
    plt.close(fig)
    return final_path


# ---------------------------------------------------------------------------
# Bulk driver
# ---------------------------------------------------------------------------

def animate_all(output_dir: str, *, fps: int = 30,
                n_frames: int = 150) -> dict[str, str]:
    """Produce every animation in ``output_dir``.

    Returns a dict mapping the logical animation name → actual file path
    written (after any GIF fallback).
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    results: dict[str, str] = {}
    results["kink_scattering"] = animate_kink_scattering(
        str(out / "kink_scattering.mp4"), fps=fps, n_frames=n_frames)
    results["cosmology_desaturation"] = animate_cosmology_desaturation(
        str(out / "cosmology_desaturation.mp4"), fps=fps, n_frames=n_frames)
    results["cone_bouncing"] = animate_cone_bouncing(
        str(out / "cone_bouncing.mp4"), fps=fps, n_frames=n_frames)
    results["lattice_2d_evolution"] = animate_lattice_2d(
        str(out / "lattice_2d_evolution.mp4"), fps=fps, n_frames=n_frames)
    return results


__all__ = [
    "animate_kink_scattering",
    "animate_cosmology_desaturation",
    "animate_cone_bouncing",
    "animate_lattice_2d",
    "animate_all",
]
