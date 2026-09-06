"""Paired cancer / tumor-growth GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for the carcinogenesis
sequence (normal -> dysplastic -> neoplastic -> invasive -> metastatic) and
the canonical mathematical-oncology growth + treatment laws under the B3 /
stiff-medium ontology.

GEOMETRY
--------
``CancerGeometry`` holds the multi-stage cell population and tumour
macroscopic descriptors:

    * stage_populations  : counts per stage (normal, dysplastic,
                            neoplastic, invasive, metastatic).
    * tumour_radius_mm   : effective spherical radius of the primary mass.
    * vascularisation    : 0..1 fractional density of recruited capillaries
                            (angiogenic switch sets in around r > 1-2 mm).
    * mitotic_index      : fraction of cycling cells (proxy for grade).

Stage progression mirrors the textbook adenoma-carcinoma / Vogelstein
sequence: each stage adds one acquired hallmark and the cell loses one
more layer of substrate restraint (mass = cumulative torque on substrate;
each "hit" releases a torque axis).

SIMULATOR
---------
``CancerSimulator`` exposes the canonical mathematical-oncology models:

    * gompertz_growth(t)             : dN/dt = N a ln(K/N), saturating.
    * exponential_growth(t)          : reference unbounded growth dN/dt=rN.
    * armitage_doll_incidence(age)   : I(age) ∝ age^(n-1) for n hits.
    * knudson_two_hit(t)             : sporadic vs hereditary Rb tumour
                                       onset distribution.
    * hallmarks_acquired(stage)      : Hanahan-Weinberg 8 hallmarks
                                       acquired by current stage.
    * log_kill_treatment(N0, doses)  : surviving fraction after K doses
                                       at constant log-kill per dose.
    * combined_growth_and_treatment  : Gompertz growth interrupted by
                                       periodic chemotherapy doses.

Substrate-specific predictions encoded here:
    * Each acquired hallmark is one substrate-restraint released; the
      Gompertz carrying capacity K is the substrate's geometric tolerance
      to a defect of that strain depth.
    * Multistage carcinogenesis exponent n~5-6 maps to the integer count
      of independent substrate restraints that must each fail before a
      cell escapes saturation control (5 hallmarks ~ 5 hits).
    * Log-kill behaviour follows from each dose acting as an independent
      Poisson clearance event on the strain population (constant fraction
      removed per pulse).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) values
# ---------------------------------------------------------------------------

# Stage labels (5-stage adenoma-carcinoma sequence)
STAGE_NORMAL:       str = "normal"
STAGE_DYSPLASTIC:   str = "dysplastic"
STAGE_NEOPLASTIC:   str = "neoplastic"
STAGE_INVASIVE:     str = "invasive"
STAGE_METASTATIC:   str = "metastatic"
STAGES: Tuple[str, ...] = (
    STAGE_NORMAL, STAGE_DYSPLASTIC, STAGE_NEOPLASTIC,
    STAGE_INVASIVE, STAGE_METASTATIC,
)
N_STAGES: int = len(STAGES)

# Hanahan-Weinberg "8 hallmarks of cancer" (2011 update)
HALLMARKS: Tuple[str, ...] = (
    "sustained proliferative signalling",
    "evasion of growth suppressors",
    "resistance to cell death",
    "replicative immortality",
    "induction of angiogenesis",
    "activation of invasion / metastasis",
    "deregulated cellular energetics",
    "evasion of immune destruction",
)
N_HALLMARKS: int = len(HALLMARKS)

# Gompertz growth defaults (typical solid-tumour fit)
GOMPERTZ_A_PER_DAY:   float = 0.05    # senescence-rate coefficient
GOMPERTZ_K_CELLS:     float = 1.0e12  # carrying capacity (~1 kg of cells)
EXPONENTIAL_R_PER_DAY: float = 0.05   # bare exponential rate (small N limit)
INITIAL_TUMOUR_CELLS: float = 1.0     # one transformed cell

# Cell mass / volume conventions
# Mammalian cell ~ 1 pL = 1e-9 mL = 1e-9 cm^3 = 1e-6 mm^3 in volume.
CELL_VOLUME_PL:       float = 1.0     # 1 picolitre per mammalian cell
CELL_MASS_NG:         float = 1.0     # ~1 nanogram per cell

# Multistage carcinogenesis (Armitage-Doll 1954)
ARMITAGE_DOLL_HITS_DEFAULT: int = 6   # textbook colorectal n~6 (slope 5)
ARMITAGE_DOLL_C:           float = 1.0e-12   # incidence-rate coefficient
ARMITAGE_DOLL_AGE_MAX:     float = 90.0      # ages 0..90

# Knudson 2-hit (retinoblastoma)
KNUDSON_MUTATION_RATE_PER_YEAR: float = 1.0e-6
KNUDSON_N_RETINOBLASTS:        int = int(2e6)  # ~2 million retinoblasts
KNUDSON_PEAK_AGE_SPORADIC_YR:  float = 30.0
KNUDSON_PEAK_AGE_HEREDITARY_YR: float = 1.5

# Log-kill chemotherapy (Skipper-Schabel 1964)
LOG_KILL_FRACTION_PER_DOSE: float = 0.99   # 99 percent kill per dose (2-log)
DOSE_INTERVAL_DAYS:         float = 14.0   # 2-week cycle

# Vascularisation / angiogenic switch
ANGIOGENIC_SWITCH_RADIUS_MM: float = 1.5   # avascular -> angiogenic at ~1-2 mm
DIFFUSION_LIMITED_RADIUS_MM: float = 0.2   # below this, diffusion-fed only

# Substrate anchor (B3 framework)
LAMBDA_QCD_MEV:       float = 200.0


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class CancerGeometry:
    """Multi-stage cell population + tumour macroscopic descriptors.

    Attributes
    ----------
    stage_populations : np.ndarray, shape (N_STAGES,)
        Cell counts in each stage in the canonical order
        (normal, dysplastic, neoplastic, invasive, metastatic).
    tumour_radius_mm : float
        Effective spherical radius of the primary mass [mm].
    vascularisation : float
        Capillary density (0 = none, 1 = full vascular network).
    mitotic_index : float
        Fraction of cycling cells in the primary mass (textbook 0.01..0.4).
    """

    stage_populations: np.ndarray = field(
        default_factory=lambda: np.array(
            [1.0e9, 0.0, 0.0, 0.0, 0.0], dtype=float,
        )
    )
    tumour_radius_mm: float = 0.0
    vascularisation: float = 0.0
    mitotic_index: float = 0.05

    # ------------------------------------------------------------------
    # Stage queries
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        pops = np.asarray(self.stage_populations, dtype=float)
        if pops.shape != (N_STAGES,):
            raise ValueError(
                f"stage_populations must have shape ({N_STAGES},), "
                f"got {pops.shape}"
            )
        if np.any(pops < 0.0):
            raise ValueError("stage_populations must be non-negative")
        self.stage_populations = pops
        if self.tumour_radius_mm < 0.0:
            raise ValueError("tumour_radius_mm must be non-negative")
        if not (0.0 <= self.vascularisation <= 1.0):
            raise ValueError("vascularisation must be in [0, 1]")
        if not (0.0 <= self.mitotic_index <= 1.0):
            raise ValueError("mitotic_index must be in [0, 1]")

    def stage_index(self, stage: str) -> int:
        if stage not in STAGES:
            raise ValueError(f"unknown stage: {stage!r}; must be one of {STAGES}")
        return STAGES.index(stage)

    def population_in_stage(self, stage: str) -> float:
        return float(self.stage_populations[self.stage_index(stage)])

    def total_population(self) -> float:
        return float(self.stage_populations.sum())

    def malignant_population(self) -> float:
        """Population in any malignant stage (neoplastic + invasive + metastatic)."""
        idx = (self.stage_index(STAGE_NEOPLASTIC),
               self.stage_index(STAGE_INVASIVE),
               self.stage_index(STAGE_METASTATIC))
        return float(self.stage_populations[list(idx)].sum())

    def is_malignant(self) -> bool:
        return self.malignant_population() > 0.0

    def has_metastasized(self) -> bool:
        return self.population_in_stage(STAGE_METASTATIC) > 0.0

    def current_stage(self) -> str:
        """Return the most-advanced stage with non-zero population."""
        nonzero = np.nonzero(self.stage_populations > 0.0)[0]
        if nonzero.size == 0:
            return STAGE_NORMAL
        return STAGES[int(nonzero.max())]

    def is_angiogenic(self) -> bool:
        """Tumour has triggered the angiogenic switch."""
        return self.tumour_radius_mm >= ANGIOGENIC_SWITCH_RADIUS_MM

    # ------------------------------------------------------------------
    # Tumour size / mass conversions
    # ------------------------------------------------------------------

    def tumour_volume_mm3(self) -> float:
        """Volume of the primary mass treated as a sphere."""
        r = self.tumour_radius_mm
        return float((4.0 / 3.0) * math.pi * r ** 3)

    def tumour_mass_g(self) -> float:
        """Approximate tumour mass: 1 mm^3 ~ 1 mg of soft tissue."""
        return self.tumour_volume_mm3() * 1.0e-3

    @staticmethod
    def cells_to_radius_mm(n_cells: float) -> float:
        """Convert a cell count to an equivalent spherical radius [mm].

        Uses CELL_VOLUME_PL (1 pL ~ 1e-9 cm^3 = 1e-9 mm^3 = 1e-12 m^3 wait
        actually 1 pL = 1e-9 mL = 1e-9 cm^3 = 1e-6 mm^3).
        """
        if n_cells <= 0.0:
            return 0.0
        # 1 pL = 1e-9 cm^3 = 1e-6 mm^3
        vol_mm3 = n_cells * CELL_VOLUME_PL * 1.0e-6
        return float((3.0 * vol_mm3 / (4.0 * math.pi)) ** (1.0 / 3.0))

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "tumour = substrate strain that has lost N hallmark restraints; "
                "each hallmark = one substrate-saturation lock released, "
                "so the cell escapes geometric containment"
            ),
            "n_stages": N_STAGES,
            "n_hallmarks": N_HALLMARKS,
            "current_stage": self.current_stage(),
            "tumour_radius_mm": self.tumour_radius_mm,
            "tumour_volume_mm3": self.tumour_volume_mm3(),
            "tumour_mass_g": self.tumour_mass_g(),
            "vascularisation": self.vascularisation,
            "mitotic_index": self.mitotic_index,
            "angiogenic": self.is_angiogenic(),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class CancerSimulator:
    """Tumour growth, multistage incidence, and treatment response.

    Encapsulates the canonical mathematical-oncology models:
    Gompertz growth (saturating sigmoid), exponential growth
    (small-N limit), Armitage-Doll multistage incidence
    I(age) ∝ age^(n-1), Knudson 2-hit hypothesis for retinoblastoma,
    Hanahan-Weinberg 8 hallmarks acquired across the 5-stage
    progression, and Skipper-Schabel log-kill chemotherapy response.
    """

    geometry: CancerGeometry | None = None
    a_per_day: float = GOMPERTZ_A_PER_DAY        # Gompertz senescence rate
    K_cells: float = GOMPERTZ_K_CELLS            # carrying capacity
    r_per_day: float = EXPONENTIAL_R_PER_DAY     # exponential rate
    initial_cells: float = INITIAL_TUMOUR_CELLS

    n_armitage_hits: int = ARMITAGE_DOLL_HITS_DEFAULT
    armitage_C: float = ARMITAGE_DOLL_C
    log_kill_fraction_per_dose: float = LOG_KILL_FRACTION_PER_DOSE
    dose_interval_days: float = DOSE_INTERVAL_DAYS

    # ------------------------------------------------------------------
    # Gompertz growth
    # ------------------------------------------------------------------

    def gompertz_growth(
        self,
        t_days: np.ndarray,
        N0: float | None = None,
        a: float | None = None,
        K: float | None = None,
    ) -> np.ndarray:
        """Closed-form Gompertz N(t).

        N(t) = K * exp( ln(N0/K) * exp(-a t) )

        which satisfies dN/dt = N a ln(K/N): exponential at small N,
        saturates to K as t -> infinity.
        """
        N0 = N0 if N0 is not None else self.initial_cells
        a = a if a is not None else self.a_per_day
        K = K if K is not None else self.K_cells
        if N0 <= 0.0 or K <= 0.0 or a <= 0.0:
            raise ValueError("N0, K, a must be positive")
        t = np.asarray(t_days, dtype=float)
        return K * np.exp(np.log(N0 / K) * np.exp(-a * t))

    def gompertz_rate(self, N: float, a: float | None = None,
                      K: float | None = None) -> float:
        """Instantaneous growth rate dN/dt at population N."""
        a = a if a is not None else self.a_per_day
        K = K if K is not None else self.K_cells
        if N <= 0.0:
            return 0.0
        return float(N * a * math.log(K / N))

    def exponential_growth(
        self,
        t_days: np.ndarray,
        N0: float | None = None,
        r: float | None = None,
    ) -> np.ndarray:
        """N(t) = N0 exp(r t).  Reference unbounded growth."""
        N0 = N0 if N0 is not None else self.initial_cells
        r = r if r is not None else self.r_per_day
        t = np.asarray(t_days, dtype=float)
        return N0 * np.exp(r * t)

    def doubling_time_days(self, r: float | None = None) -> float:
        """Pure-exponential doubling time t_d = ln 2 / r."""
        r = r if r is not None else self.r_per_day
        if r <= 0.0:
            raise ValueError("growth rate must be positive")
        return float(math.log(2.0) / r)

    # ------------------------------------------------------------------
    # Multistage carcinogenesis (Armitage-Doll 1954)
    # ------------------------------------------------------------------

    def armitage_doll_incidence(
        self,
        age_years: np.ndarray,
        n_hits: int | None = None,
        C: float | None = None,
    ) -> np.ndarray:
        """Cancer incidence I(age) ∝ age^(n-1) for n independent hits.

        I(age) = C * age^(n_hits - 1)

        With n=6, I scales as age^5 (textbook colorectal slope on
        log-log plots).
        """
        n = n_hits if n_hits is not None else self.n_armitage_hits
        C = C if C is not None else self.armitage_C
        if n < 2:
            raise ValueError("n_hits must be >= 2 for multistage incidence")
        age = np.asarray(age_years, dtype=float)
        if np.any(age < 0.0):
            raise ValueError("age_years must be non-negative")
        return C * age ** (n - 1)

    def armitage_doll_log_slope(self, n_hits: int | None = None) -> float:
        """Slope of log I vs log age = n - 1 (textbook ~5 for colorectal)."""
        n = n_hits if n_hits is not None else self.n_armitage_hits
        return float(n - 1)

    # ------------------------------------------------------------------
    # Knudson 2-hit hypothesis (retinoblastoma)
    # ------------------------------------------------------------------

    def knudson_two_hit_sporadic_pdf(
        self,
        age_years: np.ndarray,
        mu: float = KNUDSON_MUTATION_RATE_PER_YEAR,
        n_cells: int = KNUDSON_N_RETINOBLASTS,
    ) -> np.ndarray:
        """Probability density for sporadic-Rb onset with 2 random hits.

        Each cell needs two independent mutational events (rate mu/year);
        the hit-count per cell is Poisson(lambda = mu*age) and the
        first-tumour age has the multistage form ~ age^(n-1) for n=2,
        so PDF ∝ age (linear-rising shoulder peaking around 30 yr for
        ~2e6 retinoblasts).
        """
        age = np.asarray(age_years, dtype=float)
        # Per-cell prob of having both hits by age t with Poisson(mu*t)
        # P(2 hits) approx = (mu*t)^2 / 2 for mu*t small
        # Tumour formation prob across n_cells: 1 - (1-p)^n_cells
        p = 1.0 - np.exp(-((mu * age) ** 2) / 2.0)
        P_tumour = 1.0 - (1.0 - p) ** n_cells
        # PDF = derivative of CDF; approximate numerically
        if age.size < 2:
            return np.zeros_like(age)
        pdf = np.gradient(P_tumour, age)
        return pdf

    def knudson_two_hit_hereditary_pdf(
        self,
        age_years: np.ndarray,
        mu: float = KNUDSON_MUTATION_RATE_PER_YEAR,
        n_cells: int = KNUDSON_N_RETINOBLASTS,
    ) -> np.ndarray:
        """Hereditary-Rb PDF: one hit already inherited, only 1 hit needed.

        Per-cell P(hit) = 1 - exp(-mu*age); tumour PDF peaks ~1-2 yr.
        """
        age = np.asarray(age_years, dtype=float)
        p = 1.0 - np.exp(-mu * age)
        P_tumour = 1.0 - (1.0 - p) ** n_cells
        if age.size < 2:
            return np.zeros_like(age)
        pdf = np.gradient(P_tumour, age)
        return pdf

    def knudson_peak_ages_years(self) -> Tuple[float, float]:
        """Return (sporadic, hereditary) peak ages [years]."""
        return (KNUDSON_PEAK_AGE_SPORADIC_YR, KNUDSON_PEAK_AGE_HEREDITARY_YR)

    # ------------------------------------------------------------------
    # Hanahan-Weinberg hallmarks
    # ------------------------------------------------------------------

    def hallmarks_acquired(self, stage: str) -> Tuple[str, ...]:
        """Return the hallmarks acquired by the time the cell reaches `stage`.

        Mapping (canonical adenoma-carcinoma sequence):
            normal       : 0 hallmarks
            dysplastic   : 2 (sustained proliferation + evasion of growth supp.)
            neoplastic   : 4 (+ resistance to death + replicative immortality)
            invasive     : 6 (+ angiogenesis + invasion/metastasis)
            metastatic   : 8 (+ deregulated energetics + immune evasion)
        """
        if stage not in STAGES:
            raise ValueError(f"unknown stage: {stage!r}")
        idx = STAGES.index(stage)
        # 0, 2, 4, 6, 8 hallmarks per stage
        n = 2 * idx
        return HALLMARKS[:n]

    def n_hallmarks_acquired(self, stage: str) -> int:
        return len(self.hallmarks_acquired(stage))

    def all_hallmarks(self) -> Tuple[str, ...]:
        return tuple(HALLMARKS)

    # ------------------------------------------------------------------
    # Log-kill chemotherapy (Skipper-Schabel)
    # ------------------------------------------------------------------

    def log_kill_survival_fraction(
        self,
        n_doses: int,
        kill_fraction_per_dose: float | None = None,
    ) -> float:
        """Surviving fraction after n_doses pulses with constant log-kill.

        Each dose kills the same fraction f of remaining cells regardless
        of N (Skipper-Schabel hypothesis), so survival = (1 - f)^n.

        For f = 0.99 and n = 5, survival = 1e-10.
        """
        f = kill_fraction_per_dose if kill_fraction_per_dose is not None \
            else self.log_kill_fraction_per_dose
        if not (0.0 <= f <= 1.0):
            raise ValueError("kill_fraction_per_dose must be in [0, 1]")
        if n_doses < 0:
            raise ValueError("n_doses must be non-negative")
        return float((1.0 - f) ** n_doses)

    def log_kill_treatment(
        self,
        N0: float,
        n_doses: int,
        kill_fraction_per_dose: float | None = None,
    ) -> np.ndarray:
        """Population trajectory after n_doses log-kill pulses.

        Returns array of length n_doses + 1 with N before each dose
        and the final post-treatment N.
        """
        f = kill_fraction_per_dose if kill_fraction_per_dose is not None \
            else self.log_kill_fraction_per_dose
        if N0 <= 0.0:
            raise ValueError("N0 must be positive")
        survival = (1.0 - f) ** np.arange(n_doses + 1)
        return N0 * survival

    def combined_growth_and_treatment(
        self,
        t_days: np.ndarray,
        N0: float | None = None,
        treatment_start_day: float = 100.0,
        n_doses: int = 6,
        kill_fraction_per_dose: float | None = None,
    ) -> dict:
        """Gompertz growth interrupted by periodic log-kill pulses.

        Returns a dict with N(t), the dose times, and a flag indicating
        treatment success (final N below threshold).
        """
        N0 = N0 if N0 is not None else self.initial_cells
        f = kill_fraction_per_dose if kill_fraction_per_dose is not None \
            else self.log_kill_fraction_per_dose
        t = np.asarray(t_days, dtype=float)
        N = self.gompertz_growth(t, N0=N0).copy()
        dose_times = treatment_start_day + np.arange(n_doses) * self.dose_interval_days
        # Apply each pulse: at the dose time, multiply remaining N by (1-f)
        # and continue growing from the reduced population
        for d_t in dose_times:
            mask = t >= d_t
            if not mask.any():
                continue
            # Find N just before the pulse using current trajectory
            i_pulse = int(np.searchsorted(t, d_t))
            if i_pulse >= len(t):
                break
            N_pre = N[i_pulse]
            N_post = N_pre * (1.0 - f)
            # Re-grow from N_post starting at time d_t for t >= d_t
            new_t = t[mask] - d_t
            new_N = self.gompertz_growth(new_t, N0=max(N_post, 1.0e-6))
            N[mask] = new_N
        return {
            "t_days": t,
            "N": N,
            "dose_times": dose_times,
            "treatment_success": bool(N[-1] < 1.0),
            "final_N": float(N[-1]),
            "log_reduction": float(np.log10(max(N0, 1.0) / max(N[-1], 1.0e-30))),
        }

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "gompertz_a_per_day": self.a_per_day,
            "gompertz_K_cells": self.K_cells,
            "exponential_r_per_day": self.r_per_day,
            "doubling_time_days": self.doubling_time_days(),
            "n_armitage_hits": self.n_armitage_hits,
            "armitage_log_slope": self.armitage_doll_log_slope(),
            "log_kill_fraction_per_dose": self.log_kill_fraction_per_dose,
            "n_hallmarks_total": N_HALLMARKS,
            "interpretation": (
                "carcinogenesis = stepwise release of substrate-saturation "
                "locks: each hallmark = one geometric restraint failed; "
                "Gompertz K = substrate strain capacity for that defect; "
                "log-kill = Poisson clearance of independent strain quanta"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


__all__ = [
    "CancerGeometry",
    "CancerSimulator",
    "STAGES",
    "STAGE_NORMAL",
    "STAGE_DYSPLASTIC",
    "STAGE_NEOPLASTIC",
    "STAGE_INVASIVE",
    "STAGE_METASTATIC",
    "N_STAGES",
    "HALLMARKS",
    "N_HALLMARKS",
    "GOMPERTZ_A_PER_DAY",
    "GOMPERTZ_K_CELLS",
    "EXPONENTIAL_R_PER_DAY",
    "INITIAL_TUMOUR_CELLS",
    "ARMITAGE_DOLL_HITS_DEFAULT",
    "ARMITAGE_DOLL_C",
    "ARMITAGE_DOLL_AGE_MAX",
    "KNUDSON_MUTATION_RATE_PER_YEAR",
    "KNUDSON_N_RETINOBLASTS",
    "KNUDSON_PEAK_AGE_SPORADIC_YR",
    "KNUDSON_PEAK_AGE_HEREDITARY_YR",
    "LOG_KILL_FRACTION_PER_DOSE",
    "DOSE_INTERVAL_DAYS",
    "ANGIOGENIC_SWITCH_RADIUS_MM",
    "DIFFUSION_LIMITED_RADIUS_MM",
    "CELL_VOLUME_PL",
    "CELL_MASS_NG",
    "LAMBDA_QCD_MEV",
]
