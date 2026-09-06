"""Generation tower visualizer — paired GEOMETRY + SIM + VIZ for the
3-generation lepton/quark tower.

WHY EXACTLY 3 GENERATIONS?
--------------------------
The substrate has D = 3 spatial dimensions, and a lepton's mass-bearing mode
is one orientation of the K_4 braided-anchor cell along a substrate axis.
With three orthogonal substrate axes (x, y, z), there are EXACTLY three
distinct stack-orientations of the K_4 cell, and therefore exactly three
charged-lepton mass-modes:

    generation 1 (electron) -> K_4 stacked along x
    generation 2 (muon)     -> K_4 stacked along y
    generation 3 (tau)      -> K_4 stacked along z

This is *not* a fit; it is the geometric statement that 3 generations are
forced from D = 3. A fourth generation would need a fourth orthogonal
substrate axis, which D = 3 cannot provide.

MASS LADDER FROM SUBSTRATE TOPOLOGY
-----------------------------------
The mass jump between generations is set by the substrate integers
(n_M, K_pair, n_R, R) of the B3 rigidity grid:

    log(m_mu / m_e) = n_M / (K_pair**4 * pi)                  # = 268/(16 pi)
    log(m_tau/m_mu) = n_M / (K_pair**4 * pi)
                       - (n_R - R) / (K_pair * (K_pair + 1))   # = 268/(16 pi) - 15/6

At canonical integers (268, 2, 18, 3) these reproduce m_mu/m_e to <0.01%
and m_tau/m_mu to ~0.93% (open NLO; see ``tau_mass_unified``).

The same generation map applies to up-type and down-type quarks with
sector-dependent prefactors; we display the full 6-lepton + 6-quark ladder
on a log-scale axis for visual comparison.

VISUALS PRODUCED
----------------
22_generation_tower.png  -- 3D plot of three orthogonal K_4 cells (one per
                            generation), each labelled with its substrate
                            axis and predicted mass.
23_lepton_quark_ladder.png -- log-scale mass ladder for all 6 leptons +
                            6 quarks, with B3 generation jumps annotated.

REFERENCES
----------
src/stiff_medium/generation_map.py        -- full lepton/quark/neutrino map
src/stiff_medium/tau_mass_unified.py      -- canonical m_tau/m_mu formula
src/stiff_medium/k4_face_pair_geometry.py -- K_4 cell geometry
src/stiff_medium/b3_constants.py          -- canonical integers
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from .b3_constants import (
    K_pair as _K_PAIR,
    K_rank as _K_RANK,
    R as _R,
    n_M as _N_M,
    n_R as _N_R,
)
from .generation_map import (
    DOWN_QUARK_MASSES_MEV,
    GenerationMap,
    LEPTON_MASSES_MEV,
    UP_QUARK_MASSES_MEV,
)
from .k4_face_pair_geometry import K4Geometry
from .tau_mass_unified import M_MUON_MEV, M_TAU_MEV, log_m_tau_over_m_mu


# ---------------------------------------------------------------------------
# Spatial axes of the substrate
# ---------------------------------------------------------------------------

D_SPATIAL: int = 3
"""Number of spatial dimensions of the substrate. Forces the number of
charged-lepton generations: one stack-orientation per orthogonal axis."""

AXIS_NAMES: Tuple[str, str, str] = ("x", "y", "z")
"""Names of the three orthogonal substrate axes."""

GEN_LABELS: Tuple[str, str, str] = ("gen 1", "gen 2", "gen 3")
"""Labels of the three generations."""

LEPTON_LABELS: Tuple[str, str, str] = ("electron", "muon", "tau")
"""Charged-lepton names by generation."""

UP_QUARK_LABELS: Tuple[str, str, str] = ("up", "charm", "top")
DOWN_QUARK_LABELS: Tuple[str, str, str] = ("down", "strange", "bottom")
NEUTRINO_LABELS: Tuple[str, str, str] = ("nu_e", "nu_mu", "nu_tau")


# ---------------------------------------------------------------------------
# 1.  Generation tower geometry
# ---------------------------------------------------------------------------

@dataclass
class GenerationTowerGeometry:
    """Three orthogonal K_4 cells stacked along the substrate axes.

    The construction is forced by D = 3: there is one K_4 cell per
    orthogonal axis, and each cell carries a single mass-bearing
    orientation. The three cells together exhaust the available
    orientations -- there is no fourth orthogonal axis to host a
    generation 4.

    Attributes
    ----------
    spatial_dim : int
        Number of spatial dimensions (always 3 in this construction).
    n_generations : int
        Number of generations forced by ``spatial_dim`` (= 3).
    cells : list[K4Geometry]
        One K_4 cell per generation. Each cell is a regular tetrahedron
        translated along its assigned axis (x, y, z) by ``axis_spacing``
        and labelled with the corresponding generation index.
    axis_spacing : float
        Centre-to-centre separation between adjacent K_4 cells (units of
        d_face). Pure visualisation knob; not a physical parameter.
    map : GenerationMap
        Substrate generation map used for mass ratios.
    n_M, K_pair, n_R, R : int
        Canonical substrate integers (defaults from b3_constants).
    """

    spatial_dim: int = D_SPATIAL
    axis_spacing: float = 3.0
    n_M: int = _N_M
    K_pair: int = _K_PAIR
    n_R: int = _N_R
    R: int = _R
    map: GenerationMap = field(default_factory=GenerationMap)
    cells: List[K4Geometry] = field(init=False)
    axis_directions: NDArray[np.float64] = field(init=False)

    def __post_init__(self) -> None:
        if self.spatial_dim != D_SPATIAL:
            raise ValueError(
                f"GenerationTowerGeometry assumes D=3 spatial dimensions; "
                f"got spatial_dim={self.spatial_dim}"
            )
        # Three orthogonal axes (rows of identity).
        self.axis_directions = np.eye(self.spatial_dim, dtype=float)

        # Build one K_4 cell per generation, translated along its axis.
        self.cells = []
        for gen in range(self.spatial_dim):
            cell = K4Geometry()
            translation = self.axis_spacing * self.axis_directions[gen]
            cell.vertices = cell.vertices + translation
            cell.face_centroids = cell.face_centroids + translation
            cell.centre = cell.centre + translation
            self.cells.append(cell)

    # ------------------------------------------------------------------ count

    @property
    def n_generations(self) -> int:
        """Number of generations forced by the spatial dimension."""
        return self.spatial_dim

    def axis_for_generation(self, gen: int) -> str:
        """Return the substrate axis name (x, y, z) for generation ``gen`` (1-3)."""
        if gen not in (1, 2, 3):
            raise ValueError(f"generation must be 1, 2, or 3; got {gen}")
        return AXIS_NAMES[gen - 1]

    # ------------------------------------------------------------------ logs

    def log_ratio_12(self) -> float:
        """log(m_mu/m_e) = n_M / (K_pair**4 * pi)."""
        return self.n_M / (self.K_pair ** 4 * math.pi)

    def log_ratio_23(self) -> float:
        """log(m_tau/m_mu) = canonical n_R-active formula.

        Delegates to ``tau_mass_unified.log_m_tau_over_m_mu``.
        """
        return log_m_tau_over_m_mu(
            n_M=self.n_M, K_pair=self.K_pair, n_R=self.n_R, R=self.R,
        )

    def ratio_12(self) -> float:
        """Predicted m_mu/m_e from the 1->2 substrate jump."""
        return math.exp(self.log_ratio_12())

    def ratio_23(self) -> float:
        """Predicted m_tau/m_mu from the 2->3 substrate jump."""
        return math.exp(self.log_ratio_23())

    def ratio_13(self) -> float:
        """Predicted m_tau/m_e = ratio_12 * ratio_23."""
        return self.ratio_12() * self.ratio_23()

    # --------------------------------------------------------------- masses

    def predicted_lepton_masses_mev(self) -> Dict[str, float]:
        """Return the three charged-lepton masses (MeV) anchored at m_e.

        m_e is taken from PDG (``LEPTON_MASSES_MEV['e']``); m_mu and m_tau
        follow from the substrate ratios.
        """
        m_e = LEPTON_MASSES_MEV["e"]
        m_mu = m_e * self.ratio_12()
        m_tau = m_mu * self.ratio_23()
        return {"e": m_e, "mu": m_mu, "tau": m_tau}

    # --------------------------------------------------------------- forced

    def is_3gen_forced(self) -> bool:
        """True iff the spatial dimension forces n_generations == 3.

        Geometric statement: D = 3 admits exactly three orthogonal axes,
        so exactly three independent K_4 stack-orientations exist.
        """
        return self.spatial_dim == 3 and self.n_generations == 3

    # -------------------------------------------------------------- summary

    def summary(self) -> Dict[str, object]:
        """Compact summary used by tests / the renderer."""
        masses = self.predicted_lepton_masses_mev()
        return {
            "spatial_dim": self.spatial_dim,
            "n_generations": self.n_generations,
            "axes": list(AXIS_NAMES),
            "lepton_labels": list(LEPTON_LABELS),
            "log_ratio_12": self.log_ratio_12(),
            "log_ratio_23": self.log_ratio_23(),
            "ratio_12_pred": self.ratio_12(),
            "ratio_23_pred": self.ratio_23(),
            "ratio_13_pred": self.ratio_13(),
            "ratio_12_obs": LEPTON_MASSES_MEV["mu"] / LEPTON_MASSES_MEV["e"],
            "ratio_23_obs": LEPTON_MASSES_MEV["tau"] / LEPTON_MASSES_MEV["mu"],
            "ratio_13_obs": LEPTON_MASSES_MEV["tau"] / LEPTON_MASSES_MEV["e"],
            "predicted_masses_mev": masses,
            "n_M": self.n_M,
            "K_pair": self.K_pair,
            "n_R": self.n_R,
            "R": self.R,
        }


# ---------------------------------------------------------------------------
# 2.  Generation tower simulator (visualization helper)
# ---------------------------------------------------------------------------

@dataclass
class GenerationTowerSimulator:
    """Visualize the 3-generation lepton/quark tower.

    Two figures are produced:

    1. ``draw_geometric_tower`` -- 3D rendering of the three orthogonal
       K_4 cells (one per generation), each labelled with its substrate
       axis and predicted lepton mass.

    2. ``draw_full_ladder`` -- log-scale mass ladder for all 6 leptons
       (3 charged + 3 neutrinos) and 6 quarks (3 up-type + 3 down-type),
       with the substrate generation jumps annotated.

    Attributes
    ----------
    geometry : GenerationTowerGeometry
        The underlying tower geometry.
    """

    geometry: GenerationTowerGeometry = field(default_factory=GenerationTowerGeometry)

    # ---- error percentages ---------------------------------------------

    def lepton_errors_pct(self) -> Dict[str, float]:
        """Percent error of predicted charged-lepton ratios vs PDG."""
        s = self.geometry.summary()
        return {
            "mu/e": 100.0 * (s["ratio_12_pred"] - s["ratio_12_obs"]) / s["ratio_12_obs"],
            "tau/mu": 100.0 * (s["ratio_23_pred"] - s["ratio_23_obs"]) / s["ratio_23_obs"],
            "tau/e": 100.0 * (s["ratio_13_pred"] - s["ratio_13_obs"]) / s["ratio_13_obs"],
        }

    # ---- 22: 3D geometric tower ----------------------------------------

    def draw_geometric_tower(self, ax=None):
        """Render the three orthogonal K_4 cells on a 3D axis.

        Each cell is drawn at its translated centre with edges in black,
        face-centroid arrows in blue, and a generation label that lists
        the assigned substrate axis and the predicted lepton mass.

        Parameters
        ----------
        ax : matplotlib 3D axis, optional
            If None, a new figure + 3D axis is created.

        Returns
        -------
        matplotlib axis
            The axis the geometry was drawn into.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3D)

        if ax is None:
            fig = plt.figure(figsize=(11, 9))
            ax = fig.add_subplot(111, projection="3d")

        masses = self.geometry.predicted_lepton_masses_mev()
        lepton_keys = ("e", "mu", "tau")
        gen_colors = ("#d62728", "#1f77b4", "#2ca02c")  # red, blue, green

        # Draw each cell.
        for gen, cell in enumerate(self.geometry.cells):
            color = gen_colors[gen]
            v = cell.vertices
            # Vertices.
            ax.scatter(v[:, 0], v[:, 1], v[:, 2],
                       s=120, c=color, edgecolors="black", linewidths=1.5)
            # Edges.
            for a, b in cell.edges:
                ax.plot(*zip(v[a], v[b]), color=color, alpha=0.7, linewidth=1.8)
            # Face-centroid normals.
            fc = cell.face_centroids
            fn = cell.face_normals
            for i in range(4):
                ax.quiver(fc[i, 0], fc[i, 1], fc[i, 2],
                          fn[i, 0] * 0.3, fn[i, 1] * 0.3, fn[i, 2] * 0.3,
                          color=color, alpha=0.5, arrow_length_ratio=0.3)

            # Label the cell.
            centre = cell.centre
            axis_name = AXIS_NAMES[gen]
            mass_mev = masses[lepton_keys[gen]]
            mass_str = (f"{mass_mev:.4f} MeV" if mass_mev < 10
                        else f"{mass_mev:.2f} MeV" if mass_mev < 1000
                        else f"{mass_mev/1000:.3f} GeV")
            ax.text(centre[0], centre[1], centre[2] + 0.9,
                    f"gen {gen+1} ({LEPTON_LABELS[gen]})\n"
                    f"axis: {axis_name}\n"
                    f"m = {mass_str}",
                    fontsize=10, ha="center", color=color, weight="bold")

        # Draw the substrate axes themselves through the origin.
        amax = self.geometry.axis_spacing * 1.6
        for gen in range(D_SPATIAL):
            d = self.geometry.axis_directions[gen]
            ax.plot([0, amax * d[0]], [0, amax * d[1]], [0, amax * d[2]],
                    color="gray", linestyle="--", alpha=0.4, linewidth=1.2)
            ax.text(amax * d[0] * 1.05, amax * d[1] * 1.05, amax * d[2] * 1.05,
                    f"axis {AXIS_NAMES[gen]}", fontsize=9, color="gray")

        # Annotate the integer ratios.
        s = self.geometry.summary()
        title = (
            f"3-generation tower forced from D = {D_SPATIAL} substrate dims\n"
            f"K_4 stack orientations: {self.geometry.n_generations} "
            f"(one per orthogonal axis)\n"
            f"m_mu/m_e = exp(n_M/(K_pair^4 pi)) = exp(268/(16 pi)) = "
            f"{s['ratio_12_pred']:.3f}    "
            f"m_tau/m_mu = {s['ratio_23_pred']:.3f}"
        )
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_box_aspect([1, 1, 1])
        return ax

    # ---- 23: log-scale lepton + quark ladder ---------------------------

    def draw_full_ladder(self, ax=None):
        """Render a log-scale mass ladder with all 6 leptons + 6 quarks.

        Charged leptons (e, mu, tau) and the up-type and down-type quark
        towers are plotted on a log-y axis, grouped by generation. The
        substrate generation jumps log_ratio_12, log_ratio_23 are
        annotated in the figure title.

        Parameters
        ----------
        ax : matplotlib axis, optional
            If None, a new figure + axis is created.

        Returns
        -------
        matplotlib axis
            The axis the ladder was drawn into.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(11, 7))

        # Generation indices on the x axis.
        x = np.array([1, 2, 3], dtype=float)

        # Charged leptons (3).
        m_e, m_mu, m_tau = (LEPTON_MASSES_MEV[k] for k in ("e", "mu", "tau"))
        leptons = np.array([m_e, m_mu, m_tau])

        # Neutrinos (3) -- B3 NH lightest at 2.26 meV; m2, m3 from oscillation.
        DM21_SQ_EV2 = 7.41e-5
        DM31_SQ_EV2 = 2.511e-3
        m1_eV = 2.26e-3
        m2_eV = math.sqrt(m1_eV ** 2 + DM21_SQ_EV2)
        m3_eV = math.sqrt(m1_eV ** 2 + DM31_SQ_EV2)
        neutrinos = np.array([m1_eV, m2_eV, m3_eV]) * 1e-6  # eV -> MeV

        # Up-type and down-type quarks.
        up_q = np.array([UP_QUARK_MASSES_MEV[k] for k in ("u", "c", "t")])
        down_q = np.array([DOWN_QUARK_MASSES_MEV[k] for k in ("d", "s", "b")])

        # Predictions anchored at the lightest of each tower.
        s = self.geometry.summary()
        r12, r23 = s["ratio_12_pred"], s["ratio_23_pred"]
        leptons_pred = np.array([m_e, m_e * r12, m_e * r12 * r23])
        up_q_pred = np.array([up_q[0], up_q[0] * r12, up_q[0] * r12 * r23])
        down_q_pred = np.array([down_q[0], down_q[0] * r12, down_q[0] * r12 * r23])

        # Plot observed values as filled markers, predicted as open markers.
        ax.semilogy(x, leptons, "o-", color="#d62728", markersize=11,
                    label="charged leptons (PDG)", linewidth=1.8)
        ax.semilogy(x, leptons_pred, "o--", color="#d62728", markersize=11,
                    markerfacecolor="white", markeredgewidth=1.8,
                    label="charged leptons (B3, anchored at m_e)", linewidth=1.5,
                    alpha=0.7)

        ax.semilogy(x, up_q, "s-", color="#1f77b4", markersize=11,
                    label="up-type quarks (PDG)", linewidth=1.8)
        ax.semilogy(x, up_q_pred, "s--", color="#1f77b4", markersize=11,
                    markerfacecolor="white", markeredgewidth=1.8,
                    label="up-type (B3, anchored at m_u)", linewidth=1.5,
                    alpha=0.7)

        ax.semilogy(x, down_q, "^-", color="#2ca02c", markersize=11,
                    label="down-type quarks (PDG)", linewidth=1.8)
        ax.semilogy(x, down_q_pred, "^--", color="#2ca02c", markersize=11,
                    markerfacecolor="white", markeredgewidth=1.8,
                    label="down-type (B3, anchored at m_d)", linewidth=1.5,
                    alpha=0.7)

        ax.semilogy(x, neutrinos, "v-", color="#9467bd", markersize=11,
                    label="neutrinos (NH, B3 m_1=2.26 meV)", linewidth=1.5)

        # Annotate each particle.
        all_labels = (
            (x, leptons, LEPTON_LABELS, "#d62728", 1.10),
            (x, up_q, UP_QUARK_LABELS, "#1f77b4", 1.40),
            (x, down_q, DOWN_QUARK_LABELS, "#2ca02c", 1.40),
            (x, neutrinos, NEUTRINO_LABELS, "#9467bd", 1.40),
        )
        for xs, ys, labs, color, dy in all_labels:
            for xi, yi, lab in zip(xs, ys, labs):
                ax.text(xi, yi * dy, lab,
                        fontsize=9, ha="center", color=color, weight="bold")

        # Decorations.
        ax.set_xticks(x)
        ax.set_xticklabels([f"gen {int(g)}\n(axis {a})"
                            for g, a in zip(x, AXIS_NAMES)],
                           fontsize=10)
        ax.set_xlim(0.5, 3.5)
        ax.set_ylabel("mass [MeV]  (log scale)")
        title = (
            "Lepton + quark generation ladder -- one mode per substrate axis "
            f"(D = {D_SPATIAL} -> 3 generations)\n"
            f"B3 jumps: log(m_2/m_1) = n_M/(K_pair^4 pi) = {s['log_ratio_12']:.4f}, "
            f"log(m_3/m_2) = ... - (n_R-R)/(K_pair(K_pair+1)) = {s['log_ratio_23']:.4f}"
        )
        ax.set_title(title, fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8, ncol=2)

        # Y range covers neutrinos at ~10^-9 MeV through top at ~10^5 MeV.
        ax.set_ylim(1e-10, 5e6)
        return ax


# ---------------------------------------------------------------------------
# Convenience entry point used by render_all_visuals.py
# ---------------------------------------------------------------------------

def make_tower_figures():
    """Return ``(fig_tower, fig_ladder)`` matplotlib Figures.

    The renderer in ``scripts/render_all_visuals.py`` uses these.
    """
    import matplotlib.pyplot as plt

    geom = GenerationTowerGeometry()
    sim = GenerationTowerSimulator(geometry=geom)

    fig_tower = plt.figure(figsize=(11, 9))
    ax_tower = fig_tower.add_subplot(111, projection="3d")
    sim.draw_geometric_tower(ax=ax_tower)

    fig_ladder, ax_ladder = plt.subplots(figsize=(11, 7))
    sim.draw_full_ladder(ax=ax_ladder)

    return fig_tower, fig_ladder


__all__ = [
    "AXIS_NAMES",
    "D_SPATIAL",
    "DOWN_QUARK_LABELS",
    "GEN_LABELS",
    "GenerationTowerGeometry",
    "GenerationTowerSimulator",
    "LEPTON_LABELS",
    "NEUTRINO_LABELS",
    "UP_QUARK_LABELS",
    "make_tower_figures",
]
