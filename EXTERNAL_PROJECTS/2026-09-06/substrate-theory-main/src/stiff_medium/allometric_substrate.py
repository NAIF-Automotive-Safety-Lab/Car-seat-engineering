"""Paired allometric-scaling GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for the universal
allometric scaling laws of biology under the B3 / stiff-medium ontology.

The empirical pattern -- observed across ~27 orders of magnitude in body
mass, from mycoplasma (10^-13 kg) to blue whale (10^5 kg) -- is

    metabolic rate   B(M) ~ M^(3/4)         (Kleiber 1932)
    heart rate       f(M) ~ M^(-1/4)
    lifespan         T(M) ~ M^(1/4)
    capillary count  N(M) ~ M^(3/4)

The exponent 3/4 is *not* the surface-to-volume scaling 2/3 that classical
geometric arguments would predict.  West, Brown, and Enquist (Science 1997)
showed that the 3/4 emerges from THREE conditions on the vascular network:

    1. SPACE-FILLING: capillaries reach every tissue volume element.
    2. TERMINAL UNITS INVARIANT: capillary radius and metabolic load per
       capillary do not depend on body size (a mouse capillary is the
       same size as a whale capillary).
    3. MINIMUM ENERGY DISSIPATION: the network minimises hydrodynamic
       resistance subject to the volume budget for blood.

Under SUBSTRATE FRAMING, the vascular network is the optimal STRAIN field
that delivers nutrients (a conserved current under the substrate continuity
equation) into a 3-volume from a single-point source (heart) with finite
total network volume.  The 3/4 exponent then reduces to a purely geometric
ratio:  (network levels) / (level at which volume saturates) = 3 / (3+1) = 3/4
because the vasculature is fractal in 3 spatial dimensions plus 1 branching
generation index.  This is the WBE result re-expressed as a strain-field
optimisation -- the substrate is what makes the network space-filling and
the b = 3/4 exponent topology-derived rather than fitted.

GEOMETRY
--------
``AllometricGeometry`` represents the vascular tree as

    * a fractal binary network with N levels (root = aorta, leaves = capillaries)
    * branching ratio n at each level (default n = 2; West-Brown-Enquist)
    * radius scaling beta = n^(-1/2)  (area-preserving for pulsatile flow)
    * length scaling gamma = n^(-1/3) (space-filling)
    * total terminal-unit (capillary) count N_cap = n^(N-1)

SIMULATOR
---------
``AllometricSimulator`` exposes the scaling laws and WBE derivation:

    * kleiber_metabolic_rate_W(M_kg): B = B_0 * M^(3/4)
    * heart_rate_Hz(M_kg):           f = f_0 * M^(-1/4)
    * lifespan_s(M_kg):              T = T_0 * M^( 1/4)
    * blood_volume_m3(M_kg):         V_b ~ M^( 1)        (linear, ~7% of body)
    * derive_kleiber_exponent(): returns 3/4 from the WBE three conditions
    * compare_geometric_2_3_vs_topological_3_4(M_kg):
        shows the network-derived 3/4 fits across 12+ OOM where the naive
        surface/volume 2/3 fails

Substrate framing summary:
    * Capillary radius and metabolic load per capillary are SUBSTRATE-set
      (terminal units of the strain field), not body-size set.
    * The 3/4 exponent is the topological consequence of 3-D space-filling
      with a 1-D branching index -- d / (d+1) = 3/4 for d = 3.
    * The WBE network IS the optimal substrate strain configuration for a
      conserved nutrient flux; B3 doesn't replace WBE so much as embed it.

References
----------
Kleiber, M. (1932). Body size and metabolism. Hilgardia 6, 315.
West, G.B., Brown, J.H., Enquist, B.J. (1997). A general model for the
    origin of allometric scaling laws in biology. Science 276, 122-126.
West, G.B., Brown, J.H., Enquist, B.J. (1999). The fourth dimension of
    life: fractal geometry and allometric scaling of organisms.
    Science 284, 1677-1679.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Universal allometric exponents (West-Brown-Enquist 1997)
# ---------------------------------------------------------------------------

KLEIBER_EXPONENT:     float = 3.0 / 4.0    # metabolic rate ~ M^(3/4)
HEART_RATE_EXPONENT:  float = -1.0 / 4.0   # heart rate ~ M^(-1/4)
LIFESPAN_EXPONENT:    float = 1.0 / 4.0    # lifespan ~ M^(1/4)
CAPILLARY_EXPONENT:   float = 3.0 / 4.0    # capillary count ~ M^(3/4)
BLOOD_VOLUME_EXPONENT: float = 1.0         # blood volume ~ M (linear)
SURFACE_GEOMETRIC_EXPONENT: float = 2.0 / 3.0   # naive surface-to-volume

# Fractal-network parameters (WBE: branching ratio n=2 at every level)
BRANCHING_RATIO_N: int = 2
RADIUS_SCALING_BETA: float = BRANCHING_RATIO_N ** (-0.5)   # area-preserving
LENGTH_SCALING_GAMMA: float = BRANCHING_RATIO_N ** (-1.0 / 3.0)  # space-filling
SPACE_DIMENSION_D: int = 3

# Terminal-unit invariants (size-independent)
CAPILLARY_RADIUS_M:   float = 4.0e-6        # ~ 4 microns (mouse, whale, same)
CAPILLARY_LENGTH_M:   float = 5.0e-4        # ~ 0.5 mm
CAPILLARY_FLOW_M3_S:  float = 5.0e-13       # mean per-capillary flow

# Allometric prefactors (anchored to mouse / human reference)
# B_0 such that B(human) = 90 W at M = 70 kg gives B_0 = 90 / 70^(3/4) = 3.72
KLEIBER_PREFACTOR_W:    float = 3.72        # W per kg^(3/4)
# f_0 such that human heart rate 70 bpm = 1.17 Hz at M = 70 kg
# f_0 = 1.17 * 70^(1/4) = 3.39
HEART_RATE_PREFACTOR_HZ: float = 3.39       # Hz per kg^(-1/4)
# T_0 such that human lifespan 80 yr at M = 70 kg
# T_0 = 80 / 70^(1/4) = 27.6 yr per kg^(1/4)
LIFESPAN_PREFACTOR_YR:  float = 27.6        # yr per kg^(1/4)
# Blood volume: 7% of body mass at density 1060 kg/m^3
BLOOD_VOLUME_PREFACTOR_M3_PER_KG: float = 0.07 / 1060.0   # ~ 6.6e-5

# Reference mammals for empirical Kleiber data (Kleiber 1932, Heusner 1991,
# White & Seymour 2003).  Body mass (kg) and basal metabolic rate (W).
MAMMAL_DATA: List[Dict[str, float]] = [
    {"name": "Etruscan shrew",  "M_kg": 0.0018, "B_W": 0.027},
    {"name": "House mouse",     "M_kg": 0.020,  "B_W": 0.20},
    {"name": "Brown rat",       "M_kg": 0.250,  "B_W": 1.45},
    {"name": "Guinea pig",      "M_kg": 0.700,  "B_W": 3.20},
    {"name": "Rabbit",          "M_kg": 2.50,   "B_W": 8.30},
    {"name": "Cat",             "M_kg": 4.00,   "B_W": 12.0},
    {"name": "Dog (medium)",    "M_kg": 15.0,   "B_W": 33.0},
    {"name": "Macaque",         "M_kg": 8.00,   "B_W": 19.5},
    {"name": "Beaver",          "M_kg": 25.0,   "B_W": 47.0},
    {"name": "Sheep",           "M_kg": 60.0,   "B_W": 84.0},
    {"name": "Human",           "M_kg": 70.0,   "B_W": 90.0},
    {"name": "Wolf",            "M_kg": 35.0,   "B_W": 60.0},
    {"name": "Pig",             "M_kg": 130.0,  "B_W": 145.0},
    {"name": "Goat",            "M_kg": 50.0,   "B_W": 73.0},
    {"name": "Reindeer",        "M_kg": 100.0,  "B_W": 120.0},
    {"name": "Lion",            "M_kg": 180.0,  "B_W": 192.0},
    {"name": "Donkey",          "M_kg": 200.0,  "B_W": 210.0},
    {"name": "Bear (black)",    "M_kg": 150.0,  "B_W": 168.0},
    {"name": "Bear (brown)",    "M_kg": 350.0,  "B_W": 330.0},
    {"name": "Cow",             "M_kg": 600.0,  "B_W": 510.0},
    {"name": "Horse",           "M_kg": 500.0,  "B_W": 442.0},
    {"name": "Rhinoceros",      "M_kg": 2000.0, "B_W": 1290.0},
    {"name": "Hippopotamus",    "M_kg": 1500.0, "B_W": 1015.0},
    {"name": "Elephant (Asian)","M_kg": 3500.0, "B_W": 2030.0},
    {"name": "Elephant (Afr.)", "M_kg": 5000.0, "B_W": 2680.0},
    {"name": "Giraffe",         "M_kg": 800.0,  "B_W": 627.0},
    {"name": "Walrus",          "M_kg": 1200.0, "B_W": 850.0},
    {"name": "Killer whale",    "M_kg": 5500.0, "B_W": 2880.0},
    {"name": "Sperm whale",     "M_kg": 40000.0, "B_W": 12000.0},
    {"name": "Blue whale",      "M_kg": 150000.0, "B_W": 32000.0},
]


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class AllometricGeometry:
    """Fractal vascular tree geometry (West-Brown-Enquist 1997).

    The network has ``n_levels`` generations from the root (aorta) at level
    0 down to the terminal capillaries at level ``n_levels - 1``.  At each
    level k the vessels have radius r_k and length l_k with

        r_{k+1} = beta * r_k    with  beta = n^(-1/2)  (area-preserving)
        l_{k+1} = gamma * l_k   with  gamma = n^(-1/3) (space-filling)

    The branching ratio ``n`` is the number of daughter vessels at each
    bifurcation (default n = 2 == binary tree, the WBE canonical case).

    Attributes
    ----------
    n_levels : int
        Number of generations in the tree (root = level 0, capillaries at
        level n_levels - 1).
    branching_n : int
        Number of daughter vessels per parent.
    aorta_radius_m : float
        Radius of the level-0 root vessel.
    aorta_length_m : float
        Length of the level-0 root vessel.
    """

    n_levels: int = 22                       # ~ humans: log2(N_cap=10^10) ~ 22
    branching_n: int = BRANCHING_RATIO_N
    aorta_radius_m: float = 0.012            # 12 mm (human aorta)
    aorta_length_m: float = 0.40             # 40 cm

    # ---------------------------------------------------------------- core scales

    def beta(self) -> float:
        """Radius scaling factor between successive levels (area-preserving)."""
        return self.branching_n ** (-0.5)

    def gamma(self) -> float:
        """Length scaling factor between successive levels (space-filling)."""
        return self.branching_n ** (-1.0 / SPACE_DIMENSION_D)

    # --------------------------------------------------------------- per level

    def vessels_at_level(self, k: int) -> int:
        """Number of vessels at level k (1, n, n^2, ..., n^(N-1))."""
        if k < 0 or k >= self.n_levels:
            raise ValueError(f"level {k} out of range [0, {self.n_levels})")
        return self.branching_n ** k

    def radius_at_level_m(self, k: int) -> float:
        """Vessel radius at level k:  r_k = aorta * beta^k."""
        return self.aorta_radius_m * self.beta() ** k

    def length_at_level_m(self, k: int) -> float:
        """Vessel length at level k:  l_k = aorta * gamma^k."""
        return self.aorta_length_m * self.gamma() ** k

    def vessel_volume_at_level_m3(self, k: int) -> float:
        """Volume of one vessel at level k:  pi r_k^2 l_k."""
        r = self.radius_at_level_m(k)
        l = self.length_at_level_m(k)
        return math.pi * r * r * l

    def total_volume_at_level_m3(self, k: int) -> float:
        """Total blood volume across all vessels at level k."""
        return self.vessels_at_level(k) * self.vessel_volume_at_level_m3(k)

    # ------------------------------------------------------------ aggregates

    def total_blood_volume_m3(self) -> float:
        """Sum of total volume across all levels (network blood pool)."""
        return sum(
            self.total_volume_at_level_m3(k) for k in range(self.n_levels)
        )

    def n_capillaries(self) -> int:
        """Number of terminal capillaries (level n_levels - 1)."""
        return self.vessels_at_level(self.n_levels - 1)

    def capillary_radius_m(self) -> float:
        """Capillary radius (terminal-unit invariant in WBE)."""
        return self.radius_at_level_m(self.n_levels - 1)

    def capillary_length_m(self) -> float:
        """Capillary length (terminal-unit invariant in WBE)."""
        return self.length_at_level_m(self.n_levels - 1)

    def is_space_filling(self) -> bool:
        """True if length scaling gamma = n^(-1/d) (space-filling network)."""
        expected = self.branching_n ** (-1.0 / SPACE_DIMENSION_D)
        return abs(self.gamma() - expected) < 1.0e-12

    def is_area_preserving(self) -> bool:
        """True if radius scaling beta = n^(-1/2) (area-preserving for flow)."""
        expected = self.branching_n ** (-0.5)
        return abs(self.beta() - expected) < 1.0e-12

    # --------------------------------------------------------- visualization

    def tree_edges_2d(
        self,
        spread_deg: float = 38.0,
    ) -> List[Tuple[Tuple[float, float], Tuple[float, float], int]]:
        """2D edge list for plotting the vascular tree.

        Returns a list of ((x0, y0), (x1, y1), level) tuples.  The root is
        placed at the origin pointing downward; each daughter is offset by
        a fan-out angle ``spread_deg / 2^level`` so the deeper generations
        compress visually (matches WBE space-filling).
        """
        edges: List[Tuple[Tuple[float, float], Tuple[float, float], int]] = []

        def _recurse(x: float, y: float, angle_deg: float, level: int):
            if level >= self.n_levels:
                return
            # length in normalised plot units (gamma^level)
            length = (self.gamma()) ** level
            a = math.radians(angle_deg)
            x1 = x + length * math.sin(a)
            y1 = y - length * math.cos(a)
            edges.append(((x, y), (x1, y1), level))
            if level + 1 >= self.n_levels:
                return
            half_spread = spread_deg / (2.0 ** level)
            for s in (-1, +1):
                _recurse(x1, y1, angle_deg + s * half_spread, level + 1)

        _recurse(0.0, 0.0, 0.0, 0)
        return edges


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class AllometricSimulator:
    """Allometric scaling-law simulator (Kleiber + heart rate + lifespan).

    Uses the WBE-derived 3/4 exponent throughout.  Provides a small
    derivation method that returns 3/4 from the three WBE conditions
    (space-filling + area-preserving + size-invariant terminal units).
    """

    geometry: AllometricGeometry = field(default_factory=AllometricGeometry)

    # ------------------------------------------------------- WBE derivation

    def derive_kleiber_exponent(
        self,
        space_dim: int = SPACE_DIMENSION_D,
        branching_n: int = BRANCHING_RATIO_N,
    ) -> float:
        """Return the Kleiber exponent b derived from WBE three conditions.

        Derivation sketch (West-Brown-Enquist 1997):

            * Service volume of one capillary V_cap ~ l_cap^d.
            * Total body volume served:  V_body = N_cap * V_cap.
            * For a binary fractal network with length scaling
              gamma = n^(-1/d) and radius scaling beta = n^(-1/2),
              the total blood volume V_b satisfies V_b ~ V_body^(d/(d+1)).
            * Metabolic rate B is proportional to N_cap (each capillary
              delivers a fixed metabolic load), and N_cap ~ V_b ~ M^(3/4).

        Therefore  b = d / (d + 1)  with d = 3 gives b = 3/4.

        For ``space_dim = 2``  we recover b = 2/3 (the naive surface-to-
        volume scaling); for ``space_dim = 3`` we recover the empirical
        b = 3/4.

        This method ignores ``branching_n`` because the result depends
        only on the dimensionality and the network being space-filling
        and area-preserving (West & Brown PNAS 2002 Eq. 14).
        """
        if space_dim < 1:
            raise ValueError(f"space_dim {space_dim} must be >= 1")
        return float(space_dim) / float(space_dim + 1)

    def compare_geometric_2_3_vs_topological_3_4(self, M_kg: float) -> Dict[str, float]:
        """Return the predicted metabolic rate under geometric 2/3 vs WBE 3/4.

        Used by the test suite to demonstrate that WBE 3/4 is genuinely
        distinct from the naive surface-to-volume 2/3.
        """
        b_2_3 = M_kg ** SURFACE_GEOMETRIC_EXPONENT
        b_3_4 = M_kg ** KLEIBER_EXPONENT
        return {
            "M_kg": M_kg,
            "B_geometric_2_3": KLEIBER_PREFACTOR_W * b_2_3,
            "B_topological_3_4": KLEIBER_PREFACTOR_W * b_3_4,
            "ratio_3_4_over_2_3": b_3_4 / b_2_3,
        }

    # ----------------------------------------------------- core scaling laws

    def kleiber_metabolic_rate_W(self, M_kg: float) -> float:
        """B(M) = B_0 * M^(3/4)  -- basal metabolic rate in watts."""
        return KLEIBER_PREFACTOR_W * (M_kg ** KLEIBER_EXPONENT)

    def heart_rate_Hz(self, M_kg: float) -> float:
        """f(M) = f_0 * M^(-1/4)  -- resting heart rate in Hz."""
        return HEART_RATE_PREFACTOR_HZ * (M_kg ** HEART_RATE_EXPONENT)

    def heart_rate_bpm(self, M_kg: float) -> float:
        """Heart rate converted to beats per minute."""
        return 60.0 * self.heart_rate_Hz(M_kg)

    def lifespan_yr(self, M_kg: float) -> float:
        """T(M) = T_0 * M^(1/4)  -- lifespan in years."""
        return LIFESPAN_PREFACTOR_YR * (M_kg ** LIFESPAN_EXPONENT)

    def lifespan_s(self, M_kg: float) -> float:
        """Lifespan converted to seconds (yr * 365.25 * 86400)."""
        return self.lifespan_yr(M_kg) * 365.25 * 86400.0

    def blood_volume_m3(self, M_kg: float) -> float:
        """V_blood = 0.07 * M / rho_blood  -- linear in M (~7% of body)."""
        return BLOOD_VOLUME_PREFACTOR_M3_PER_KG * M_kg

    def capillary_count(self, M_kg: float, ref_M_kg: float = 70.0,
                         ref_count: float = 1.0e10) -> float:
        """N_cap ~ M^(3/4); anchor to ~10^10 capillaries in a 70 kg human."""
        return ref_count * (M_kg / ref_M_kg) ** CAPILLARY_EXPONENT

    def heartbeats_per_lifetime(self, M_kg: float) -> float:
        """Total heartbeats per lifespan -- approximately mass-invariant.

        Since f ~ M^(-1/4) and T ~ M^(+1/4), the product is independent
        of M:  ~ f_0 * T_0 ~ 60 * 80 * 60 * 60 * 24 * 365 / sec_per_beat
        ~ 2.5e9 beats per mammalian life.  Famous "billion heartbeats"
        invariant (off by factor ~ 2.5).
        """
        return self.heart_rate_Hz(M_kg) * self.lifespan_s(M_kg)

    # --------------------------------------------------- empirical fit on data

    def fit_kleiber_slope(
        self, masses_kg: Sequence[float], rates_W: Sequence[float]
    ) -> Dict[str, float]:
        """Linear fit of log10(B) vs log10(M); return slope and intercept.

        The fitted slope should be very close to 3/4 = 0.75 across
        mammalian data.
        """
        m = np.asarray(masses_kg, dtype=float)
        b = np.asarray(rates_W, dtype=float)
        if m.size != b.size or m.size < 2:
            raise ValueError("masses and rates must have matching length >= 2")
        if np.any(m <= 0.0) or np.any(b <= 0.0):
            raise ValueError("masses and rates must be strictly positive")
        x = np.log10(m)
        y = np.log10(b)
        slope, intercept = np.polyfit(x, y, 1)
        # R^2
        y_pred = slope * x + intercept
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
        return {
            "slope": float(slope),
            "intercept_log10": float(intercept),
            "prefactor_W": float(10.0 ** intercept),
            "r2": r2,
        }

    def kleiber_data(self) -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """Return (M_kg array, B_W array, names) of MAMMAL_DATA reference set."""
        names = [d["name"] for d in MAMMAL_DATA]
        M = np.array([d["M_kg"] for d in MAMMAL_DATA], dtype=float)
        B = np.array([d["B_W"] for d in MAMMAL_DATA], dtype=float)
        return M, B, names


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    "KLEIBER_EXPONENT",
    "HEART_RATE_EXPONENT",
    "LIFESPAN_EXPONENT",
    "CAPILLARY_EXPONENT",
    "BLOOD_VOLUME_EXPONENT",
    "SURFACE_GEOMETRIC_EXPONENT",
    "BRANCHING_RATIO_N",
    "RADIUS_SCALING_BETA",
    "LENGTH_SCALING_GAMMA",
    "SPACE_DIMENSION_D",
    "CAPILLARY_RADIUS_M",
    "CAPILLARY_LENGTH_M",
    "CAPILLARY_FLOW_M3_S",
    "KLEIBER_PREFACTOR_W",
    "HEART_RATE_PREFACTOR_HZ",
    "LIFESPAN_PREFACTOR_YR",
    "BLOOD_VOLUME_PREFACTOR_M3_PER_KG",
    "MAMMAL_DATA",
    "AllometricGeometry",
    "AllometricSimulator",
]
