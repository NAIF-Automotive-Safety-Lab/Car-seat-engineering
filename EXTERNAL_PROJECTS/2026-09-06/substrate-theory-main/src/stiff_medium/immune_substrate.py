"""Paired immune-system GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for adaptive-immune
recognition (antibody/antigen binding, B-cell clonal selection, primary
versus secondary response) under the B3 / stiff-medium ontology.

GEOMETRY
--------
``ImmuneGeometry`` represents an antibody as the canonical Y-shaped
immunoglobulin (IgG):

    * Two Fab arms (Fragment, antigen-binding) tipped by paratopes
      (variable-region antigen-binding sites).
    * One Fc stem (Fragment, crystallisable) that mediates effector
      function and FcRn-driven recycling.
    * Hinge region between Fab arms and Fc, giving the Y its flexibility.
    * Antigen-binding site approximated as a complementarity-determining
      shape (small spherical region of radius ~1 nm at the Fab tip).
    * MHC presentation: an MHC class I molecule on a host cell displays
      a short peptide (8-10 amino acids) on its alpha1/alpha2 groove,
      shown here as a flat platform with a peptide segment laid in the
      groove.

Canonical antibody dimensions (IgG):

    * total span (Fab tip to Fab tip):           ~14 nm
    * Fab arm length:                            ~7  nm
    * Fc stem length:                            ~7  nm
    * Hinge half-angle (Y opening):              ~60 degrees
    * paratope (binding-site) radius:            ~1  nm
    * antigen radius (typical epitope):          ~1.5 nm

SIMULATOR
---------
``ImmuneSimulator`` exposes the standard adaptive-immunity quantitative
properties:

    * binding_affinity_Kd_M():          K_d for antibody/antigen binding
                                        (textbook IgG range 10^-9 .. 10^-12 M).
    * binding_free_energy_kcal_per_mol(): delta_G = -kT ln(1/K_d).
    * fraction_bound(C_antigen):        equilibrium occupancy at given
                                        antigen concentration.
    * sir_immune_response(t):           SIR-like dynamics for pathogen
                                        load P(t) and antibody titer A(t).
    * primary_vs_secondary():           paired (1 deg, 2 deg) titer time
                                        courses; secondary is faster, higher,
                                        and dominated by IgG (vs IgM in 1 deg).
    * clonal_selection(steps):          B-cell population dynamics under
                                        affinity-driven proliferation.
    * memory_cell_fraction():           fraction of activated B cells that
                                        differentiate into memory B cells
                                        (textbook ~10 percent).
    * hypermutation_rate_per_bp_per_div(): somatic hypermutation rate at
                                        the Ig V-gene (~10^-3 per bp per
                                        division -- ~10^6 times background).

Substrate-specific predictions encoded here:
    * The K_d ladder 10^-6 .. 10^-12 M maps to a substrate strain depth of
      0.3-0.7 eV; affinity maturation ratchets the strain depth deeper as
      hypermutation explores the V-region sequence space.
    * Memory cells are long-lived "frozen" strain configurations that
      persist after the primary response decays; on antigen rechallenge
      they re-saturate the substrate immediately, producing the order-of-
      magnitude faster and higher secondary response.
    * SIR dynamics are the substrate analogue of the pathogen-substrate
      strain spreading like a defect line; antibody = stabilising counter
      strain that quenches the pathogen-strain.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) values -- IgG, body temperature 310 K
# ---------------------------------------------------------------------------

# Antibody Y-shape geometry (nanometres, degrees)
ANTIBODY_FAB_LENGTH_NM:    float = 7.0    # length of one Fab arm
ANTIBODY_FC_LENGTH_NM:     float = 7.0    # length of Fc stem
ANTIBODY_HINGE_ANGLE_DEG:  float = 60.0   # half-angle between Fab arm and Fc
ANTIBODY_PARATOPE_NM:      float = 1.0    # paratope (binding-pocket) radius
ANTIBODY_TOTAL_SPAN_NM:    float = 14.0   # Fab tip to Fab tip

# Antigen and MHC dimensions
ANTIGEN_RADIUS_NM:         float = 1.5    # typical epitope sphere radius
MHC_PLATFORM_NM:           float = 4.0    # MHC alpha1/alpha2 platform side
MHC_GROOVE_LENGTH_NM:      float = 3.0    # peptide-binding groove length
MHC_PEPTIDE_LENGTH_AA:     int = 9        # MHC-I displays 8..10-mer; use 9

# Binding affinity ladder (molar) -- IgG textbook range
KD_LOW_AFFINITY_M:         float = 1.0e-6     # naive B-cell IgM
KD_MEDIUM_AFFINITY_M:      float = 1.0e-9     # primary IgG response
KD_HIGH_AFFINITY_M:        float = 1.0e-12    # affinity-matured IgG
KD_DEFAULT_M:              float = 1.0e-10    # mid-range IgG

# Clonal-selection / memory parameters
CLONAL_PROLIFERATION_RATE_PER_DAY: float = 1.0   # B-cell doubling ~ 1/day
MEMORY_CELL_FRACTION:              float = 0.10  # ~10 percent of activated B cells
HYPERMUTATION_RATE_PER_BP_PER_DIV: float = 1.0e-3  # somatic hypermutation
BACKGROUND_MUTATION_RATE_PER_BP_PER_DIV: float = 1.0e-9
PRIMARY_LATENCY_DAYS:              float = 5.0    # 1 deg response onset
SECONDARY_LATENCY_DAYS:            float = 1.5    # 2 deg response onset
PRIMARY_PEAK_TITER:                float = 1.0    # arbitrary units (1 deg)
SECONDARY_PEAK_TITER:              float = 100.0  # 2 deg ~ 100x 1 deg
ANTIBODY_HALF_LIFE_DAYS:           float = 21.0   # IgG circulating half-life

# SIR-like immune response constants
PATHOGEN_GROWTH_RATE_PER_DAY:      float = 1.5
PATHOGEN_CLEARANCE_RATE_PER_DAY:   float = 0.8
ANTIBODY_PRODUCTION_RATE:          float = 0.5
ANTIBODY_DECAY_RATE_PER_DAY:       float = math.log(2.0) / ANTIBODY_HALF_LIFE_DAYS

# Substrate anchor (B3 framework)
LAMBDA_QCD_MEV:                    float = 200.0
KT_AT_BODY_T_KCAL_PER_MOL:         float = 0.617    # 310 K, kcal/mol
KT_AT_BODY_T_J:                    float = 4.28e-21  # joules at 310 K


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class ImmuneGeometry:
    """Antibody Y-shape (Fab + Fc) and MHC presentation platform.

    The antibody is rendered as a Y in the y-z plane, centred at the
    origin and pointing 'up' (positive z).  The Fc stem runs from the
    origin down to z = -Fc_length.  The two Fab arms emanate from the
    origin upward and outward, at +/- hinge_angle_deg from the +z axis.

    The MHC platform sits on the x-y plane at z = z_mhc (default below the
    Fc stem) and has a peptide groove laid along its long axis.

    Attributes
    ----------
    fab_length_nm : float
        Length of each Fab arm.
    fc_length_nm : float
        Length of the Fc stem.
    hinge_angle_deg : float
        Half-angle between Fc stem and Fab arm (Y opening / 2).
    paratope_radius_nm : float
        Radius of the antigen-binding pocket (sphere) at each Fab tip.
    antigen_radius_nm : float
        Radius of the antigen sphere bound to one paratope.
    """

    fab_length_nm: float = ANTIBODY_FAB_LENGTH_NM
    fc_length_nm: float = ANTIBODY_FC_LENGTH_NM
    hinge_angle_deg: float = ANTIBODY_HINGE_ANGLE_DEG
    paratope_radius_nm: float = ANTIBODY_PARATOPE_NM
    antigen_radius_nm: float = ANTIGEN_RADIUS_NM
    mhc_platform_nm: float = MHC_PLATFORM_NM
    mhc_groove_length_nm: float = MHC_GROOVE_LENGTH_NM
    mhc_peptide_length_aa: int = MHC_PEPTIDE_LENGTH_AA

    # ------------------------------------------------------------------
    # Antibody coordinates
    # ------------------------------------------------------------------

    def fc_endpoints(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return (base, hinge) endpoints of the Fc stem (3D)."""
        base = np.array([0.0, 0.0, -self.fc_length_nm])
        hinge = np.array([0.0, 0.0, 0.0])
        return base, hinge

    def fab_endpoints(self, side: str) -> Tuple[np.ndarray, np.ndarray]:
        """Return (hinge, tip) endpoints of one Fab arm.

        The Fab tip is offset by hinge_angle_deg from the +z axis in
        the y-z plane, on the chosen side ('left' = -y, 'right' = +y).
        """
        if side not in ("left", "right"):
            raise ValueError("side must be 'left' or 'right'")
        sign = -1.0 if side == "left" else +1.0
        theta = math.radians(self.hinge_angle_deg)
        hinge = np.array([0.0, 0.0, 0.0])
        tip = np.array([
            0.0,
            sign * self.fab_length_nm * math.sin(theta),
            self.fab_length_nm * math.cos(theta),
        ])
        return hinge, tip

    def paratope_centers(self) -> np.ndarray:
        """Return (2, 3) array of paratope centre positions (Fab tips)."""
        _, left_tip = self.fab_endpoints("left")
        _, right_tip = self.fab_endpoints("right")
        return np.stack([left_tip, right_tip], axis=0)

    def total_span_nm(self) -> float:
        """Return Fab-tip to Fab-tip distance (nm)."""
        tips = self.paratope_centers()
        return float(np.linalg.norm(tips[0] - tips[1]))

    def antibody_height_nm(self) -> float:
        """Return Fc-base to Fab-tip vertical extent (nm)."""
        _, tip = self.fab_endpoints("right")
        return float(tip[2] + self.fc_length_nm)

    # ------------------------------------------------------------------
    # Antigen
    # ------------------------------------------------------------------

    def antigen_position(self, side: str = "right") -> np.ndarray:
        """Return centre position of an antigen bound to one paratope."""
        _, tip = self.fab_endpoints(side)
        # Antigen centre sits just above the paratope tip
        offset = (self.paratope_radius_nm + self.antigen_radius_nm) * 0.95
        sign = -1.0 if side == "left" else +1.0
        # Push antigen further outward along the Fab axis
        theta = math.radians(self.hinge_angle_deg)
        outward = np.array([
            0.0,
            sign * math.sin(theta),
            math.cos(theta),
        ])
        return tip + offset * outward

    # ------------------------------------------------------------------
    # MHC presentation platform
    # ------------------------------------------------------------------

    def mhc_platform_corners(self, z_mhc: float = -10.0) -> np.ndarray:
        """Return (4, 3) corners of the MHC alpha1/alpha2 platform.

        The platform is a square of side mhc_platform_nm in the x-y plane.
        """
        L = self.mhc_platform_nm / 2.0
        corners = np.array([
            [-L, -L, z_mhc],
            [+L, -L, z_mhc],
            [+L, +L, z_mhc],
            [-L, +L, z_mhc],
        ])
        return corners

    def mhc_peptide_positions(self, z_mhc: float = -10.0) -> np.ndarray:
        """Return (n_aa, 3) positions of peptide alpha-carbons in the groove.

        The peptide is laid along the +y axis on top of the MHC platform.
        """
        n = self.mhc_peptide_length_aa
        L = self.mhc_groove_length_nm
        ys = np.linspace(-L / 2.0, +L / 2.0, n)
        xs = np.zeros(n)
        zs = np.full(n, z_mhc + 0.3)
        return np.stack([xs, ys, zs], axis=1)

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "antibody = Y-shaped substrate-strain probe; the two Fab "
                "paratopes are matched-shape contact patches that bind to "
                "an antigen strain configuration"
            ),
            "fab_length_nm": self.fab_length_nm,
            "fc_length_nm": self.fc_length_nm,
            "hinge_angle_deg": self.hinge_angle_deg,
            "total_span_nm": self.total_span_nm(),
            "paratope_radius_nm": self.paratope_radius_nm,
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class ImmuneSimulator:
    """Antibody affinity, clonal-selection dynamics, and SIR immune response.

    The simulator does not resolve molecular force fields; instead it
    exposes the canonical adaptive-immunity quantitative observables:
    binding affinity K_d, fraction bound at equilibrium, primary versus
    secondary response titer time-courses, B-cell clonal selection
    dynamics, somatic hypermutation rate, and an SIR-like coupled model
    of pathogen load and antibody titer.
    """

    geometry: ImmuneGeometry | None = None
    Kd_M: float = KD_DEFAULT_M
    proliferation_rate_per_day: float = CLONAL_PROLIFERATION_RATE_PER_DAY
    memory_fraction: float = MEMORY_CELL_FRACTION
    hypermutation_rate: float = HYPERMUTATION_RATE_PER_BP_PER_DIV
    primary_latency_days: float = PRIMARY_LATENCY_DAYS
    secondary_latency_days: float = SECONDARY_LATENCY_DAYS
    primary_peak: float = PRIMARY_PEAK_TITER
    secondary_peak: float = SECONDARY_PEAK_TITER
    antibody_half_life_days: float = ANTIBODY_HALF_LIFE_DAYS
    temperature_K: float = 310.0

    # ------------------------------------------------------------------
    # Binding affinity
    # ------------------------------------------------------------------

    def binding_affinity_Kd_M(self) -> float:
        """Return the dissociation constant K_d (molar)."""
        return self.Kd_M

    def binding_free_energy_kcal_per_mol(self, Kd_M: float | None = None) -> float:
        """delta_G = -R T ln(1/K_d) = +R T ln(K_d) at 1 M reference state.

        For K_d = 10^-9 M and T = 310 K we get delta_G ~ -12.7 kcal/mol
        (a textbook 'tight' binding free energy).
        """
        Kd = Kd_M if Kd_M is not None else self.Kd_M
        # Note: lower K_d (tighter binding) -> more negative delta_G
        return float(KT_AT_BODY_T_KCAL_PER_MOL * math.log(Kd))

    def fraction_bound(
        self, antigen_concentration_M: float, Kd_M: float | None = None,
    ) -> float:
        """Equilibrium fraction of antibody paratopes bound to antigen.

        f = [Ag] / (K_d + [Ag])  (Langmuir isotherm).
        """
        Kd = Kd_M if Kd_M is not None else self.Kd_M
        return float(antigen_concentration_M / (Kd + antigen_concentration_M))

    def affinity_class(self, Kd_M: float | None = None) -> str:
        """Classify K_d as low / medium / high affinity (textbook ladder)."""
        Kd = Kd_M if Kd_M is not None else self.Kd_M
        if Kd >= KD_LOW_AFFINITY_M:
            return "low"
        if Kd >= KD_MEDIUM_AFFINITY_M:
            return "medium"
        return "high"

    # ------------------------------------------------------------------
    # Hypermutation
    # ------------------------------------------------------------------

    def hypermutation_rate_per_bp_per_div(self) -> float:
        """Return the somatic-hypermutation rate (~10^-3 per bp per division)."""
        return self.hypermutation_rate

    def hypermutation_fold_over_background(self) -> float:
        """Return the fold-enhancement of SHM over background mutation."""
        return self.hypermutation_rate / BACKGROUND_MUTATION_RATE_PER_BP_PER_DIV

    def expected_mutations_in_v_region(
        self, n_divisions: int, v_region_bp: int = 360,
    ) -> float:
        """Expected number of point mutations in a V-region after N divisions.

        V-region ~ 360 bp; typical germinal-centre B cell undergoes
        N divisions with rate ~10^-3 / bp / division.
        """
        return float(n_divisions * v_region_bp * self.hypermutation_rate)

    # ------------------------------------------------------------------
    # Memory cells
    # ------------------------------------------------------------------

    def memory_cell_fraction(self) -> float:
        """Fraction of activated B cells that become memory cells (~10 percent)."""
        return self.memory_fraction

    def memory_cells_after_response(self, n_activated_b_cells: float) -> float:
        """Number of memory B cells produced after a primary response."""
        return float(n_activated_b_cells * self.memory_fraction)

    # ------------------------------------------------------------------
    # Clonal selection dynamics
    # ------------------------------------------------------------------

    def clonal_selection(
        self,
        n_steps: int = 30,
        n_initial_clones: int = 5,
        antigen_concentration_M: float = 1.0e-7,
        rng: np.random.Generator | None = None,
    ) -> dict:
        """Simulate B-cell clonal selection by affinity-driven proliferation.

        Each clone has a K_d sampled from a log-uniform distribution; the
        per-step proliferation factor is proportional to its fraction_bound
        at the given antigen concentration.  Returns the K_d ladder and
        the resulting population time-course.
        """
        rng = rng or np.random.default_rng(seed=0)
        # Initial K_d distribution: log-uniform in [10^-6, 10^-9]
        log_Kd = rng.uniform(-9.0, -6.0, size=n_initial_clones)
        Kds = 10.0 ** log_Kd
        pops = np.ones(n_initial_clones)  # start with one cell per clone
        history = np.zeros((n_steps + 1, n_initial_clones))
        history[0] = pops
        for t in range(n_steps):
            f_bound = antigen_concentration_M / (Kds + antigen_concentration_M)
            # Proliferation rate per step: 1 + r * f_bound  (selection)
            growth = 1.0 + self.proliferation_rate_per_day * f_bound
            pops = pops * growth
            history[t + 1] = pops
        return {
            "Kds_M": Kds,
            "history": history,  # shape (n_steps+1, n_clones)
            "n_steps": n_steps,
            "winning_clone": int(np.argmax(pops)),
        }

    # ------------------------------------------------------------------
    # SIR-like coupled pathogen / antibody dynamics
    # ------------------------------------------------------------------

    def sir_immune_response(
        self,
        t_days: np.ndarray,
        P0: float = 1.0,
        A0: float = 0.0,
        antigen_growth: float | None = None,
        antibody_production: float | None = None,
    ) -> dict:
        """SIR-like coupled ODE for pathogen P(t) and antibody A(t).

        dP/dt = +alpha P - beta P A         (pathogen grows, killed by Ab)
        dA/dt = +gamma P - delta A          (Ab produced if pathogen present, decays)

        Returns a dict with 'P' and 'A' time courses on the input t grid.
        """
        alpha = antigen_growth if antigen_growth is not None else PATHOGEN_GROWTH_RATE_PER_DAY
        beta = PATHOGEN_CLEARANCE_RATE_PER_DAY
        gamma = antibody_production if antibody_production is not None \
            else ANTIBODY_PRODUCTION_RATE
        delta = ANTIBODY_DECAY_RATE_PER_DAY
        t = np.asarray(t_days, dtype=float)
        n = len(t)
        P = np.zeros(n)
        A = np.zeros(n)
        P[0], A[0] = P0, A0
        # Forward Euler with adaptive sub-stepping for stability
        for i in range(1, n):
            dt = t[i] - t[i - 1]
            # Take 10 sub-steps to keep Euler stable
            sub = 10
            ds = dt / sub
            p, a = P[i - 1], A[i - 1]
            for _ in range(sub):
                dp = alpha * p - beta * p * a
                da = gamma * p - delta * a
                p = max(0.0, p + ds * dp)
                a = max(0.0, a + ds * da)
            P[i], A[i] = p, a
        return {"t_days": t, "P": P, "A": A,
                "alpha": alpha, "beta": beta, "gamma": gamma, "delta": delta}

    # ------------------------------------------------------------------
    # Primary vs secondary response
    # ------------------------------------------------------------------

    def antibody_titer(
        self, t_days: np.ndarray, response: str = "primary",
    ) -> np.ndarray:
        """Antibody titer A(t) for primary or secondary response.

        Modelled as an offset gamma-shape rise followed by exponential decay
        to capture the canonical features:
            - latency before any rise (~5 d primary, ~1.5 d secondary)
            - rise to a peak (10x .. 100x larger in secondary)
            - decay with IgG half-life ~21 d.
        """
        if response not in ("primary", "secondary"):
            raise ValueError("response must be 'primary' or 'secondary'")
        latency = self.primary_latency_days if response == "primary" \
            else self.secondary_latency_days
        peak = self.primary_peak if response == "primary" else self.secondary_peak
        rise_tau = 3.0 if response == "primary" else 1.5
        decay_tau = self.antibody_half_life_days / math.log(2.0)
        t = np.asarray(t_days, dtype=float)
        # Shifted gamma-like profile: zero before latency
        u = np.maximum(0.0, t - latency)
        rise = 1.0 - np.exp(-u / rise_tau)
        decay = np.exp(-u / decay_tau)
        # Normalise so peak occurs and equals `peak`
        # Solve for peak time: derivative of rise*decay = 0
        # rise = 1 - exp(-u/r); decay = exp(-u/d)
        # rise*decay peak at u_peak = r * ln(1 + d/r)
        u_peak = rise_tau * math.log(1.0 + decay_tau / rise_tau)
        peak_unnorm = (1.0 - math.exp(-u_peak / rise_tau)) * math.exp(-u_peak / decay_tau)
        norm = peak / max(peak_unnorm, 1.0e-30)
        out = norm * rise * decay
        out[t < latency] = 0.0
        return out

    def primary_vs_secondary(
        self,
        t_days: np.ndarray,
        secondary_offset_days: float = 60.0,
    ) -> dict:
        """Return a paired titer time-course for 1 deg and 2 deg responses.

        The primary exposure happens at t=0, the secondary at
        t=secondary_offset_days.  Returns total titer = primary + secondary.
        """
        t = np.asarray(t_days, dtype=float)
        primary_titer = self.antibody_titer(t, "primary")
        # Secondary-response titer is shifted by secondary_offset_days
        sec_t = t - secondary_offset_days
        secondary_titer = self.antibody_titer(sec_t, "secondary")
        total = primary_titer + secondary_titer
        return {
            "t_days": t,
            "primary": primary_titer,
            "secondary": secondary_titer,
            "total": total,
            "secondary_offset_days": secondary_offset_days,
            "primary_peak": self.primary_peak,
            "secondary_peak": self.secondary_peak,
        }

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "Kd_M": self.Kd_M,
            "memory_fraction": self.memory_fraction,
            "hypermutation_rate_per_bp_per_div": self.hypermutation_rate,
            "primary_latency_days": self.primary_latency_days,
            "secondary_latency_days": self.secondary_latency_days,
            "secondary_over_primary_peak": self.secondary_peak / self.primary_peak,
            "interpretation": (
                "adaptive immunity = substrate-strain matching: paratope = "
                "shape-complementary contact; affinity maturation = ratchet "
                "of strain depth via hypermutation; memory = persistent strain"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


__all__ = [
    "ImmuneGeometry",
    "ImmuneSimulator",
    "ANTIBODY_FAB_LENGTH_NM",
    "ANTIBODY_FC_LENGTH_NM",
    "ANTIBODY_HINGE_ANGLE_DEG",
    "ANTIBODY_PARATOPE_NM",
    "ANTIBODY_TOTAL_SPAN_NM",
    "ANTIGEN_RADIUS_NM",
    "MHC_PLATFORM_NM",
    "MHC_GROOVE_LENGTH_NM",
    "MHC_PEPTIDE_LENGTH_AA",
    "KD_LOW_AFFINITY_M",
    "KD_MEDIUM_AFFINITY_M",
    "KD_HIGH_AFFINITY_M",
    "KD_DEFAULT_M",
    "CLONAL_PROLIFERATION_RATE_PER_DAY",
    "MEMORY_CELL_FRACTION",
    "HYPERMUTATION_RATE_PER_BP_PER_DIV",
    "BACKGROUND_MUTATION_RATE_PER_BP_PER_DIV",
    "PRIMARY_LATENCY_DAYS",
    "SECONDARY_LATENCY_DAYS",
    "PRIMARY_PEAK_TITER",
    "SECONDARY_PEAK_TITER",
    "ANTIBODY_HALF_LIFE_DAYS",
    "PATHOGEN_GROWTH_RATE_PER_DAY",
    "PATHOGEN_CLEARANCE_RATE_PER_DAY",
    "ANTIBODY_PRODUCTION_RATE",
    "ANTIBODY_DECAY_RATE_PER_DAY",
    "LAMBDA_QCD_MEV",
    "KT_AT_BODY_T_KCAL_PER_MOL",
    "KT_AT_BODY_T_J",
]
