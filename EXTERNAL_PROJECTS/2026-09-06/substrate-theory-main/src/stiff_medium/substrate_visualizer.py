"""Substrate Visualizer.

Unified matplotlib-based visualizer for B3 framework outputs.

Provides a single ``SubstrateVisualizer`` class that produces standard
plots used across the B3 stack: 2D field snapshots, kink trajectories,
phonon dispersion relations, mass-torque comparison bar charts, PMNS
matrix heatmaps, CMB TT power spectra, falsifier timelines, the
sine-Gordon-style substrate potential, and a bulk PNG dump helper.

All methods operate non-interactively (matplotlib ``Agg`` backend) and
return ``matplotlib.figure.Figure`` objects, so callers can either save
them, embed them in reports, or close them.
"""

from __future__ import annotations

import os
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


class SubstrateVisualizer:
    """Matplotlib-based plotter for B3 substrate observables.

    Each ``plot_*`` method returns a ``matplotlib.figure.Figure`` and
    also stores it on ``self.figures`` keyed by a short name. Use
    :meth:`save_all` to write every cached figure to a directory as
    PNG files.
    """

    DEFAULT_DPI = 120
    DEFAULT_FIGSIZE = (6.0, 4.0)

    def __init__(self, style: str | None = None, dpi: int | None = None):
        self.dpi = dpi or self.DEFAULT_DPI
        self.figures: dict[str, plt.Figure] = {}
        if style:
            try:
                plt.style.use(style)
            except (OSError, ValueError):
                pass

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------
    def _new_fig(self, name: str, figsize: tuple[float, float] | None = None):
        fig, ax = plt.subplots(
            figsize=figsize or self.DEFAULT_FIGSIZE, dpi=self.dpi
        )
        self.figures[name] = fig
        return fig, ax

    # ------------------------------------------------------------------
    # 1. 2D field snapshot
    # ------------------------------------------------------------------
    def plot_field_2d(self, field: np.ndarray, t: float) -> plt.Figure:
        """Plot a 2D snapshot ``u(x, y, t)`` of the substrate field."""
        field = np.asarray(field)
        if field.ndim != 2:
            raise ValueError(
                f"field must be 2D, got shape {field.shape}"
            )
        fig, ax = self._new_fig("field_2d", figsize=(5.5, 4.5))
        im = ax.imshow(
            field,
            origin="lower",
            cmap="RdBu_r",
            aspect="equal",
            interpolation="nearest",
        )
        ax.set_title(f"Substrate field u(x, y) at t = {t:.3g}")
        ax.set_xlabel("x [grid]")
        ax.set_ylabel("y [grid]")
        fig.colorbar(im, ax=ax, label="u")
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # 2. 1D kink evolution
    # ------------------------------------------------------------------
    def plot_kink_evolution(
        self, history: Sequence[Mapping[str, float]] | np.ndarray
    ) -> plt.Figure:
        """Plot kink position vs time from a simulation history."""
        if isinstance(history, np.ndarray):
            arr = np.asarray(history)
            if arr.ndim != 2 or arr.shape[1] < 2:
                raise ValueError(
                    "ndarray history must have shape (N, >=2): [t, x, ...]"
                )
            t = arr[:, 0]
            x = arr[:, 1]
        else:
            t = np.array([h["t"] for h in history], dtype=float)
            x = np.array([h["x"] for h in history], dtype=float)
        if t.size == 0:
            raise ValueError("history is empty")

        fig, ax = self._new_fig("kink_evolution")
        ax.plot(t, x, color="#1f77b4", lw=1.6, label="kink position")
        ax.set_xlabel("time")
        ax.set_ylabel("kink x")
        ax.set_title("Kink position vs time")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # 3. Phonon dispersion
    # ------------------------------------------------------------------
    def plot_dispersion(
        self,
        omega_k: Mapping[str, np.ndarray] | tuple[np.ndarray, np.ndarray],
    ) -> plt.Figure:
        """Plot phonon dispersion ``omega(k)``.

        Accepts either a ``{"k": ..., "omega": ...}`` mapping or a
        ``(k, omega)`` tuple.
        """
        if isinstance(omega_k, Mapping):
            k = np.asarray(omega_k["k"], dtype=float)
            omega = np.asarray(omega_k["omega"], dtype=float)
        else:
            k, omega = omega_k
            k = np.asarray(k, dtype=float)
            omega = np.asarray(omega, dtype=float)
        if k.shape != omega.shape:
            raise ValueError("k and omega must have the same shape")

        fig, ax = self._new_fig("dispersion")
        ax.plot(k, omega, color="#d62728", lw=1.6, label=r"$\omega(k)$")
        ax.set_xlabel("k")
        ax.set_ylabel(r"$\omega$")
        ax.set_title("Phonon dispersion relation")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # 4. Mass-torque table
    # ------------------------------------------------------------------
    def plot_mass_torque_table(self, engine: Any) -> plt.Figure:
        """Bar chart of predicted vs observed masses.

        ``engine`` may be either:
        - an object exposing ``.predictions()`` returning a list of
          ``{name, predicted, observed}`` dicts, or
        - a list of such dicts directly, or
        - a mapping ``{name: (predicted, observed)}``.
        """
        rows = self._normalize_mass_rows(engine)
        if not rows:
            raise ValueError("engine produced no predictions")

        names = [r["name"] for r in rows]
        pred = np.array([r["predicted"] for r in rows], dtype=float)
        obs = np.array([r["observed"] for r in rows], dtype=float)

        x = np.arange(len(names))
        width = 0.4

        fig, ax = self._new_fig(
            "mass_torque", figsize=(max(6.0, 0.6 * len(names) + 2), 4.0)
        )
        ax.bar(x - width / 2, pred, width, label="predicted", color="#2ca02c")
        ax.bar(x + width / 2, obs, width, label="observed", color="#7f7f7f")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45, ha="right")
        ax.set_ylabel("mass [MeV]")
        ax.set_title("Mass-torque predictions vs observed")
        ax.legend()
        ax.grid(alpha=0.3, axis="y")
        fig.tight_layout()
        return fig

    @staticmethod
    def _normalize_mass_rows(engine: Any) -> list[dict[str, Any]]:
        if hasattr(engine, "predictions") and callable(engine.predictions):
            data = engine.predictions()
        else:
            data = engine
        rows: list[dict[str, Any]] = []
        if isinstance(data, Mapping):
            for name, vals in data.items():
                p, o = vals
                rows.append(
                    {"name": str(name), "predicted": float(p), "observed": float(o)}
                )
        else:
            for r in data:
                rows.append(
                    {
                        "name": str(r["name"]),
                        "predicted": float(r["predicted"]),
                        "observed": float(r["observed"]),
                    }
                )
        return rows

    # ------------------------------------------------------------------
    # 5. PMNS matrix heatmap
    # ------------------------------------------------------------------
    def plot_pmns_matrix(self, matrix: np.ndarray) -> plt.Figure:
        """Heatmap of ``|U_alpha,i|**2`` for a 3x3 PMNS matrix."""
        m = np.asarray(matrix)
        if m.shape != (3, 3):
            raise ValueError(f"PMNS matrix must be 3x3, got {m.shape}")
        # If complex, square the magnitude; if real & looks already-squared
        # (all in [0,1]), keep as-is.
        if np.iscomplexobj(m):
            probs = np.abs(m) ** 2
        elif np.all((m >= 0) & (m <= 1.0 + 1e-9)):
            probs = m
        else:
            probs = m ** 2

        fig, ax = self._new_fig("pmns", figsize=(4.8, 4.2))
        im = ax.imshow(probs, cmap="viridis", vmin=0.0, vmax=1.0)
        flavors = [r"$\nu_e$", r"$\nu_\mu$", r"$\nu_\tau$"]
        masses = [r"$\nu_1$", r"$\nu_2$", r"$\nu_3$"]
        ax.set_xticks(range(3))
        ax.set_yticks(range(3))
        ax.set_xticklabels(masses)
        ax.set_yticklabels(flavors)
        for i in range(3):
            for j in range(3):
                ax.text(
                    j, i, f"{probs[i, j]:.3f}",
                    ha="center", va="center",
                    color="white" if probs[i, j] < 0.5 else "black",
                    fontsize=9,
                )
        ax.set_title(r"PMNS $|U_{\alpha i}|^2$")
        fig.colorbar(im, ax=ax, label="probability")
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # 6. CMB TT power spectrum
    # ------------------------------------------------------------------
    def plot_cmb_acoustic_peaks(
        self,
        C_l: np.ndarray | Mapping[str, np.ndarray],
    ) -> plt.Figure:
        """Plot CMB TT power spectrum ``D_l = l(l+1)C_l/2pi``."""
        if isinstance(C_l, Mapping):
            ell = np.asarray(C_l.get("l", C_l.get("ell")), dtype=float)
            cl = np.asarray(C_l["C_l"], dtype=float)
        else:
            cl = np.asarray(C_l, dtype=float)
            ell = np.arange(2, 2 + cl.size, dtype=float)
        if ell.shape != cl.shape:
            raise ValueError("ell and C_l must have matching shape")
        d_l = ell * (ell + 1.0) * cl / (2.0 * np.pi)

        fig, ax = self._new_fig("cmb")
        ax.plot(ell, d_l, color="#9467bd", lw=1.4)
        ax.set_xscale("log")
        ax.set_xlabel(r"multipole $\ell$")
        ax.set_ylabel(r"$D_\ell = \ell(\ell+1)C_\ell/2\pi$ [$\mu K^2$]")
        ax.set_title("CMB TT acoustic peaks")
        ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # 7. Falsifier timeline
    # ------------------------------------------------------------------
    DEFAULT_FALSIFIERS: list[tuple[int, str]] = [
        (2025, "DESI DR2 Sigma m_nu"),
        (2026, "KATRIN final m_beta"),
        (2027, "LIGO O5 SGWB"),
        (2028, "Hyper-K nu mass ordering"),
        (2029, "CMB-S4 r constraint"),
        (2030, "DUNE delta_CP"),
        (2032, "Einstein Telescope dispersion"),
        (2035, "LiteBIRD inflation B-mode"),
    ]

    def plot_falsifier_timeline(
        self,
        events: Sequence[tuple[int, str]] | None = None,
    ) -> plt.Figure:
        """Plot a horizontal experimental falsifier timeline 2025-2035."""
        evs = list(events or self.DEFAULT_FALSIFIERS)
        if not evs:
            raise ValueError("no events to plot")
        years = np.array([e[0] for e in evs], dtype=float)
        labels = [e[1] for e in evs]

        fig, ax = self._new_fig(
            "falsifiers", figsize=(8.0, max(3.5, 0.45 * len(evs) + 1.5))
        )
        ax.hlines(0, 2025, 2035, colors="#888", lw=1.5)
        for i, (y, lbl) in enumerate(zip(years, labels)):
            offset = 0.6 if i % 2 == 0 else -0.6
            ax.plot([y], [0], "o", color="#e377c2", ms=8)
            ax.annotate(
                f"{int(y)}: {lbl}",
                xy=(y, 0),
                xytext=(y, offset),
                ha="center",
                va="bottom" if offset > 0 else "top",
                fontsize=8,
                arrowprops=dict(arrowstyle="-", color="#aaa", lw=0.6),
            )
        ax.set_xlim(2024.5, 2035.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_yticks([])
        ax.set_xlabel("year")
        ax.set_title("B3 falsifier timeline 2025-2035")
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # 8. Substrate potential V(u)
    # ------------------------------------------------------------------
    def plot_substrate_potential(
        self,
        u_range: tuple[float, float] | np.ndarray = (-2.0 * np.pi, 2.0 * np.pi),
        K: float = 1.0,
        xi: float = 1.0,
        n: int = 400,
    ) -> plt.Figure:
        """Plot the substrate potential ``V(u) = (K/xi^2)(1 - cos(u))``."""
        if isinstance(u_range, tuple):
            u = np.linspace(u_range[0], u_range[1], n)
        else:
            u = np.asarray(u_range, dtype=float)
            if u.ndim != 1 or u.size < 2:
                raise ValueError("u_range must be 1D with >=2 samples")
        amp = K / (xi ** 2)
        V = amp * (1.0 - np.cos(u))

        fig, ax = self._new_fig("potential")
        ax.plot(u, V, color="#ff7f0e", lw=1.8)
        ax.set_xlabel("u")
        ax.set_ylabel("V(u)")
        ax.set_title(r"Substrate potential $V(u)=(K/\xi^2)(1-\cos u)$")
        ax.grid(alpha=0.3)
        ax.axhline(0, color="k", lw=0.5)
        fig.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # 9. Save all
    # ------------------------------------------------------------------
    def save_all(self, directory: str) -> list[str]:
        """Save every cached figure as ``<name>.png`` in ``directory``."""
        os.makedirs(directory, exist_ok=True)
        paths: list[str] = []
        for name, fig in self.figures.items():
            path = os.path.join(directory, f"{name}.png")
            fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
            paths.append(path)
        return paths

    def close_all(self) -> None:
        """Close all cached figures and clear the cache."""
        for fig in self.figures.values():
            plt.close(fig)
        self.figures.clear()


__all__ = ["SubstrateVisualizer"]
