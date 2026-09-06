"""Paired epidemiology GEOMETRY + SIMULATOR for the substrate framework.

Models the classical compartmental epidemics (SIR, SEIR) under the B3 /
stiff-medium ontology, where contagion propagates as a phase-transition
on the contact substrate: each susceptible cell that touches an infected
cell undergoes a saturation flip, which then relaxes (recovery) or holds
(removal).  The continuous mean-field ODE limit is the standard
Kermack-McKendrick (1927) SIR system, and adding an exposed (latent)
compartment recovers the SEIR variant used since the 1980s for
HIV/influenza/COVID-19 modelling.

GEOMETRY
--------
``EpidemiologyGeometry`` defines

    * a population of ``N`` agents
    * a contact network — Erdős-Rényi G(N, p) by default, with a
      configurable mean-degree ``mean_contacts``; the small-world or
      scale-free generators sit on top of the same adjacency abstraction
    * an age structure: K age groups with population fractions
      ``age_fractions`` and a K x K mixing matrix (Mossong 2008
      polymod-style "who acquires infection from whom")
    * a per-pair effective transmissibility τ that, multiplied by the
      contact rate ``c̄``, yields the macroscopic transmission rate

        β = τ · c̄                                                    (1)

    * a recovery rate γ = 1 / D_infectious, with D_infectious in days

The basic reproduction number is then

    R_0 = β / γ = τ · c̄ · D_infectious                              (2)

For an age-structured model, R_0 is the spectral radius of the next-
generation matrix

    K_ij = (β_ij / γ_j) · (N_i / N_j)                                (3)

SIMULATOR
---------
``EpidemiologySimulator`` integrates two compartmental ODE systems:

    SIR (3 compartments):
        dS/dt = -β S I / N
        dI/dt = +β S I / N - γ I
        dR/dt = +γ I

    SEIR (4 compartments):
        dS/dt = -β S I / N
        dE/dt = +β S I / N - σ E
        dI/dt = +σ E - γ I
        dR/dt = +γ I

with σ = 1 / D_latent the rate at which exposed agents become
infectious.

Two well-known closed-form results are provided as static methods:

    herd-immunity threshold  H = 1 - 1/R_0                          (4)
    final-size relation      R_∞ / N = 1 - exp(-R_0 · R_∞ / N)      (5)

Vaccination strategies modeled:
    * uniform (random)            — moves a fraction v of S to R at t=0
    * targeted (highest-degree)   — same effect on β via reduced ⟨k⟩
    * ring (around index cases)   — phenomenological, multiplies β by
                                    (1 - α_ring) for first ``D_latent``
                                    days

Substrate interpretation:
    * S cells = unsaturated substrate volume (free to flip)
    * I cells = saturated cells in active flip-emission state
    * R cells = relaxed cells with antibody / topological lock
    * E cells (SEIR) = saturating cells inside the σ ≤ 1/2 build-up
      phase (pre-emission latency)
    * vaccine = pre-emptive strain-lock that pulls S cells into R
      without passing through I (no flip emission)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference values (clinical / epidemiological literature)
# ---------------------------------------------------------------------------

# Classic measles R_0 ~ 12 - 18 (Anderson & May 1991), seasonal flu ~ 1.3
# COVID-19 ancestral ~ 2.5 (Liu et al., J. Travel Med. 2020).
# We pick R_0 = 2.5 as the canonical worked-example.
R0_DEFAULT: float = 2.5

# Typical infectious period for an acute respiratory virus.
D_INFECTIOUS_DAYS_DEFAULT: float = 5.0

# Typical latent period (exposed to infectious).
D_LATENT_DAYS_DEFAULT: float = 3.0

# Substrate anchor (B3 framework)
LAMBDA_QCD_MEV: float = 200.0


# ---------------------------------------------------------------------------
# Closed-form epidemiological results
# ---------------------------------------------------------------------------

def herd_immunity_threshold(R0: float) -> float:
    """Herd-immunity threshold H = 1 - 1/R_0  (Anderson & May 1991)."""
    if R0 <= 1.0:
        return 0.0
    return 1.0 - 1.0 / R0


def final_size_fraction(R0: float, tol: float = 1.0e-12,
                        max_iter: int = 200) -> float:
    """Final epidemic size fraction R_∞/N for SIR with all-susceptible IC.

    Solves the implicit Kermack-McKendrick relation

        x = 1 - exp(-R_0 · x)

    by fixed-point / bisection.  Returns 0.0 for R_0 ≤ 1, otherwise the
    unique positive root in (0, 1).
    """
    if R0 <= 1.0:
        return 0.0
    # Bisection on f(x) = x - 1 + exp(-R0 x) which has f(0)=0 and
    # one positive root in (0, 1).
    lo, hi = 1.0e-12, 1.0 - 1.0e-12
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f = mid - 1.0 + math.exp(-R0 * mid)
        if abs(f) < tol:
            return mid
        # f is monotonic increasing on (0, 1) for R0 > 1
        if f < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def peak_prevalence_fraction(R0: float) -> float:
    """Maximum I/N during an SIR epidemic, given S(0)/N = 1, I(0)=0+.

    Standard result obtained by setting dI/dt = 0:

        I_max / N = 1 - 1/R_0 - ln(R_0)/R_0

    Equivalently:  i_max = 1 - (1 + ln R_0) / R_0.
    """
    if R0 <= 1.0:
        return 0.0
    return 1.0 - (1.0 + math.log(R0)) / R0


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class ContactNetwork:
    """Compact descriptor of the contact-graph distribution.

    For a true agent-based simulation the user supplies an adjacency
    matrix; for the mean-field ODE compartmental model only the average
    degree ⟨k⟩ matters because it sets the contact rate c̄.
    """

    N: int = 10_000               # population size
    mean_contacts: float = 13.0   # mean daily distinct contacts
    network_type: str = "erdos-renyi"
    p_edge: float = 0.0           # edge probability for E-R; auto-set below
    seed: Optional[int] = 42

    def __post_init__(self) -> None:
        if self.N <= 1:
            raise ValueError("N must be > 1")
        if self.p_edge == 0.0:
            # Match ⟨k⟩ = (N - 1) p
            self.p_edge = float(self.mean_contacts) / max(1.0, self.N - 1.0)
        if self.network_type not in ("erdos-renyi", "small-world",
                                      "scale-free"):
            raise ValueError(f"unknown network_type: {self.network_type!r}")

    def expected_degree(self) -> float:
        return float(self.mean_contacts)


@dataclass
class AgeStructure:
    """K-group age structure with mixing matrix.

    Default: 4 broad bands (0-19, 20-39, 40-59, 60+) using globally
    averaged 2020 demographics.  The mixing matrix is the symmetric
    POLYMOD-style "homogeneous within band, half across bands" stencil
    used as a teaching default.
    """

    n_groups: int = 4
    age_fractions: np.ndarray = field(default_factory=lambda: np.array(
        [0.25, 0.30, 0.25, 0.20]
    ))
    mixing_matrix: Optional[np.ndarray] = None   # K x K, set in __post_init__

    def __post_init__(self) -> None:
        af = np.asarray(self.age_fractions, dtype=float)
        if af.shape != (self.n_groups,):
            raise ValueError(
                f"age_fractions shape {af.shape} != ({self.n_groups},)"
            )
        if not math.isclose(float(af.sum()), 1.0, abs_tol=1.0e-9):
            raise ValueError("age_fractions must sum to 1")
        self.age_fractions = af
        if self.mixing_matrix is None:
            M = 0.5 * np.ones((self.n_groups, self.n_groups))
            np.fill_diagonal(M, 1.0)
            self.mixing_matrix = M
        M = np.asarray(self.mixing_matrix, dtype=float)
        if M.shape != (self.n_groups, self.n_groups):
            raise ValueError(
                f"mixing_matrix shape {M.shape} != "
                f"({self.n_groups}, {self.n_groups})"
            )
        self.mixing_matrix = M


@dataclass
class EpidemiologyGeometry:
    """Population + contact network + age structure for an epidemic.

    The geometry is the *static* substrate over which the simulator
    evolves a strain (infection) field.  In B3 ontology each agent is
    a substrate cell; an edge of the contact network is a coupling
    that allows strain (saturation) to hop from a flipped cell to an
    unflipped neighbour.
    """

    N:               int   = 10_000
    mean_contacts:   float = 13.0
    network_type:    str   = "erdos-renyi"
    age_structure:   Optional[AgeStructure]    = None
    contact_network: Optional[ContactNetwork]  = None

    # Disease parameters
    R0:                  float = R0_DEFAULT
    D_infectious_days:   float = D_INFECTIOUS_DAYS_DEFAULT
    D_latent_days:       float = D_LATENT_DAYS_DEFAULT

    def __post_init__(self) -> None:
        if self.contact_network is None:
            self.contact_network = ContactNetwork(
                N=self.N, mean_contacts=self.mean_contacts,
                network_type=self.network_type,
            )
        if self.age_structure is None:
            self.age_structure = AgeStructure()

    # ------------------------------------------------------------------
    # Derived disease parameters
    # ------------------------------------------------------------------

    def gamma_per_day(self) -> float:
        """Recovery rate γ = 1 / D_infectious  (per day)."""
        return 1.0 / self.D_infectious_days

    def sigma_per_day(self) -> float:
        """Latency-exit rate σ = 1 / D_latent  (per day)."""
        return 1.0 / self.D_latent_days

    def beta_per_day(self) -> float:
        """Effective transmission rate β = R_0 · γ  (per day)."""
        return self.R0 * self.gamma_per_day()

    def transmissibility_per_contact(self) -> float:
        """Per-contact transmission probability τ = β / ⟨k⟩."""
        c = self.contact_network.expected_degree()
        return self.beta_per_day() / max(c, 1.0e-12)

    def herd_immunity_threshold(self) -> float:
        return herd_immunity_threshold(self.R0)

    def next_generation_matrix(self) -> np.ndarray:
        """K x K next-generation matrix; spectral radius = R_0_eff."""
        ag = self.age_structure
        K = ag.n_groups
        gamma = self.gamma_per_day()
        # Normalise the mixing matrix so that the spectral radius equals R_0
        # for the rank-1 limit (uniform mixing).
        M = ag.mixing_matrix.copy()
        beta_pair = self.R0 * gamma  # macroscopic β
        # Heterogeneous-mixing next-generation matrix
        # K_ij = (beta_pair · M_ij · n_i) / (gamma · M_total)
        M_total = float((ag.age_fractions @ M @ ag.age_fractions))
        Kmat = np.zeros((K, K))
        for i in range(K):
            for j in range(K):
                Kmat[i, j] = (beta_pair * M[i, j] * ag.age_fractions[i]
                              / (gamma * max(M_total, 1.0e-12)))
        return Kmat

    def R0_from_next_generation(self) -> float:
        """Spectral radius of next-generation matrix (effective R_0)."""
        Kmat = self.next_generation_matrix()
        eigvals = np.linalg.eigvals(Kmat)
        return float(np.max(np.real(eigvals)))

    # ------------------------------------------------------------------
    # Substrate signature
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "agent = substrate cell; contact edge = saturation hop; "
                "infection = flip emission (σ → 1/2); recovery = relax/lock"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "N":                self.N,
            "mean_contacts":    self.contact_network.expected_degree(),
            "R0":               self.R0,
            "gamma_per_day":    self.gamma_per_day(),
            "beta_per_day":     self.beta_per_day(),
            "tau_per_contact":  self.transmissibility_per_contact(),
            "H_threshold":      self.herd_immunity_threshold(),
        }


# ---------------------------------------------------------------------------
# SIMULATOR: SIR + SEIR ODE integrator with vaccination
# ---------------------------------------------------------------------------

@dataclass
class EpidemiologySimulator:
    """SIR / SEIR integrator with optional pre-epidemic vaccination.

    Uses the classical 4th-order Runge-Kutta scheme on the mean-field
    compartmental ODEs.  Output arrays vs t (days):

        SIR  -> S, I, R
        SEIR -> S, E, I, R

    All compartment fractions sum to 1 (population-conserving).
    """

    geometry: EpidemiologyGeometry = field(default_factory=EpidemiologyGeometry)
    dt_days: float = 0.1
    duration_days: float = 180.0

    # ------------------------------------------------------------------
    # SIR right-hand side
    # ------------------------------------------------------------------

    def _rhs_sir(self, S: float, I: float, R: float,
                 beta: float, gamma: float) -> Tuple[float, float, float]:
        N = S + I + R
        dS = -beta * S * I / N
        dI = +beta * S * I / N - gamma * I
        dR = +gamma * I
        return dS, dI, dR

    def _rhs_seir(self, S: float, E: float, I: float, R: float,
                  beta: float, sigma: float, gamma: float
                  ) -> Tuple[float, float, float, float]:
        N = S + E + I + R
        dS = -beta * S * I / N
        dE = +beta * S * I / N - sigma * E
        dI = +sigma * E - gamma * I
        dR = +gamma * I
        return dS, dE, dI, dR

    # ------------------------------------------------------------------
    # SIR simulation (RK4)
    # ------------------------------------------------------------------

    def simulate_sir(self, I0_fraction: float = 1.0e-4,
                     vaccinated_fraction: float = 0.0) -> dict:
        """Integrate SIR over ``duration_days``.

        Parameters
        ----------
        I0_fraction : float
            Initial infected fraction I(0)/N (small seed).
        vaccinated_fraction : float
            Fraction of the population pre-immunised at t=0.  Moved
            directly from S to R, never passes through I.
        """
        g = self.geometry
        N = float(g.N)
        beta = g.beta_per_day()
        gamma = g.gamma_per_day()

        n_steps = int(round(self.duration_days / self.dt_days)) + 1
        t = np.linspace(0.0, self.duration_days, n_steps)
        S = np.zeros(n_steps)
        I = np.zeros(n_steps)
        R = np.zeros(n_steps)

        v = max(0.0, min(1.0, float(vaccinated_fraction)))
        I[0] = I0_fraction * N
        R[0] = v * N
        S[0] = N - I[0] - R[0]

        h = self.dt_days
        for k in range(1, n_steps):
            s0, i0, r0 = S[k - 1], I[k - 1], R[k - 1]
            k1 = self._rhs_sir(s0, i0, r0, beta, gamma)
            k2 = self._rhs_sir(s0 + 0.5 * h * k1[0], i0 + 0.5 * h * k1[1],
                               r0 + 0.5 * h * k1[2], beta, gamma)
            k3 = self._rhs_sir(s0 + 0.5 * h * k2[0], i0 + 0.5 * h * k2[1],
                               r0 + 0.5 * h * k2[2], beta, gamma)
            k4 = self._rhs_sir(s0 + h * k3[0], i0 + h * k3[1],
                               r0 + h * k3[2], beta, gamma)
            S[k] = max(0.0, s0 + (h / 6.0) * (k1[0] + 2.0 * k2[0]
                                              + 2.0 * k3[0] + k4[0]))
            I[k] = max(0.0, i0 + (h / 6.0) * (k1[1] + 2.0 * k2[1]
                                              + 2.0 * k3[1] + k4[1]))
            R[k] = max(0.0, r0 + (h / 6.0) * (k1[2] + 2.0 * k2[2]
                                              + 2.0 * k3[2] + k4[2]))

        return {
            "t_days": t,
            "S": S, "I": I, "R": R,
            "S_frac": S / N, "I_frac": I / N, "R_frac": R / N,
            "N": N, "beta": beta, "gamma": gamma,
            "R0": g.R0,
        }

    # ------------------------------------------------------------------
    # SEIR simulation (RK4)
    # ------------------------------------------------------------------

    def simulate_seir(self, E0_fraction: float = 1.0e-4,
                      I0_fraction: float = 0.0,
                      vaccinated_fraction: float = 0.0) -> dict:
        g = self.geometry
        N = float(g.N)
        beta = g.beta_per_day()
        gamma = g.gamma_per_day()
        sigma = g.sigma_per_day()

        n_steps = int(round(self.duration_days / self.dt_days)) + 1
        t = np.linspace(0.0, self.duration_days, n_steps)
        S = np.zeros(n_steps)
        E = np.zeros(n_steps)
        I = np.zeros(n_steps)
        R = np.zeros(n_steps)

        v = max(0.0, min(1.0, float(vaccinated_fraction)))
        E[0] = E0_fraction * N
        I[0] = I0_fraction * N
        R[0] = v * N
        S[0] = N - E[0] - I[0] - R[0]

        h = self.dt_days
        for k in range(1, n_steps):
            s0, e0, i0, r0 = S[k - 1], E[k - 1], I[k - 1], R[k - 1]
            k1 = self._rhs_seir(s0, e0, i0, r0, beta, sigma, gamma)
            k2 = self._rhs_seir(s0 + 0.5 * h * k1[0], e0 + 0.5 * h * k1[1],
                                i0 + 0.5 * h * k1[2], r0 + 0.5 * h * k1[3],
                                beta, sigma, gamma)
            k3 = self._rhs_seir(s0 + 0.5 * h * k2[0], e0 + 0.5 * h * k2[1],
                                i0 + 0.5 * h * k2[2], r0 + 0.5 * h * k2[3],
                                beta, sigma, gamma)
            k4 = self._rhs_seir(s0 + h * k3[0], e0 + h * k3[1],
                                i0 + h * k3[2], r0 + h * k3[3],
                                beta, sigma, gamma)
            S[k] = max(0.0, s0 + (h / 6.0) * (k1[0] + 2 * k2[0]
                                              + 2 * k3[0] + k4[0]))
            E[k] = max(0.0, e0 + (h / 6.0) * (k1[1] + 2 * k2[1]
                                              + 2 * k3[1] + k4[1]))
            I[k] = max(0.0, i0 + (h / 6.0) * (k1[2] + 2 * k2[2]
                                              + 2 * k3[2] + k4[2]))
            R[k] = max(0.0, r0 + (h / 6.0) * (k1[3] + 2 * k2[3]
                                              + 2 * k3[3] + k4[3]))

        return {
            "t_days": t,
            "S": S, "E": E, "I": I, "R": R,
            "S_frac": S / N, "E_frac": E / N,
            "I_frac": I / N, "R_frac": R / N,
            "N": N, "beta": beta, "sigma": sigma, "gamma": gamma,
            "R0": g.R0,
        }

    # ------------------------------------------------------------------
    # Vaccination strategies
    # ------------------------------------------------------------------

    def vaccinate_uniform(self, fraction: float) -> dict:
        """Pre-immunise a fraction uniformly; returns full SIR run."""
        return self.simulate_sir(vaccinated_fraction=fraction)

    def vaccinate_targeted(self, fraction: float,
                           degree_amplification: float = 2.0) -> dict:
        """Targeted vaccination of high-degree nodes.

        Mean-field proxy: removing the highest-degree fraction f of nodes
        reduces the average remaining contact rate ⟨k⟩ by a factor
        ``(1 - degree_amplification · f)`` (clipped to ≥ 0.05) more
        strongly than uniform.  This yields a smaller effective R_0 for
        the same coverage, which is the standard theoretical motivation
        for ring/contact-based strategies (Pastor-Satorras & Vespignani
        2002).
        """
        g = self.geometry
        # Effective β reduction
        f = max(0.0, min(1.0, fraction))
        scale = max(0.05, 1.0 - degree_amplification * f)
        beta_eff = g.beta_per_day() * scale
        # Build a temporary geometry with reduced R_0
        R0_eff = beta_eff / g.gamma_per_day()
        new_geom = EpidemiologyGeometry(
            N=g.N, mean_contacts=g.mean_contacts,
            network_type=g.network_type,
            R0=R0_eff,
            D_infectious_days=g.D_infectious_days,
            D_latent_days=g.D_latent_days,
        )
        new_sim = EpidemiologySimulator(
            geometry=new_geom, dt_days=self.dt_days,
            duration_days=self.duration_days,
        )
        out = new_sim.simulate_sir(vaccinated_fraction=f)
        out["R0_effective"] = R0_eff
        return out

    def vaccinate_ring(self, fraction: float,
                       suppression: float = 0.6) -> dict:
        """Ring vaccination around index cases.

        Mean-field proxy: multiplies β by (1 - suppression · fraction)
        for the entire epidemic.  Conceptually closer to NPI + targeted
        contact tracing than to bulk vaccination.
        """
        g = self.geometry
        f = max(0.0, min(1.0, fraction))
        beta_eff = g.beta_per_day() * (1.0 - suppression * f)
        R0_eff = beta_eff / g.gamma_per_day()
        new_geom = EpidemiologyGeometry(
            N=g.N, mean_contacts=g.mean_contacts,
            network_type=g.network_type,
            R0=R0_eff,
            D_infectious_days=g.D_infectious_days,
            D_latent_days=g.D_latent_days,
        )
        new_sim = EpidemiologySimulator(
            geometry=new_geom, dt_days=self.dt_days,
            duration_days=self.duration_days,
        )
        out = new_sim.simulate_sir(vaccinated_fraction=0.0)
        out["R0_effective"] = R0_eff
        return out

    # ------------------------------------------------------------------
    # Analysis helpers
    # ------------------------------------------------------------------

    @staticmethod
    def peak_infected_fraction(result: dict) -> float:
        """Empirical max(I(t))/N from a simulation result."""
        return float(np.max(result["I_frac"]))

    @staticmethod
    def peak_time_days(result: dict) -> float:
        idx = int(np.argmax(result["I_frac"]))
        return float(result["t_days"][idx])

    @staticmethod
    def final_recovered_fraction(result: dict) -> float:
        return float(result["R_frac"][-1])

    @staticmethod
    def attack_rate(result: dict) -> float:
        """Cumulative attack rate = (R(∞) − R(0)) / N."""
        return float(result["R_frac"][-1] - result["R_frac"][0])


# ---------------------------------------------------------------------------
# Convenience: parametric R_0 sweep used by the visualizer
# ---------------------------------------------------------------------------

def r0_sweep(R0_values: np.ndarray,
             D_infectious_days: float = D_INFECTIOUS_DAYS_DEFAULT,
             N: int = 10_000,
             duration_days: float = 180.0,
             dt_days: float = 0.1,
             I0_fraction: float = 1.0e-4) -> dict:
    """Return a dict {R0: SIR run-result} for each R_0 in the sweep."""
    runs = {}
    for R0 in R0_values:
        geom = EpidemiologyGeometry(
            N=N, R0=float(R0),
            D_infectious_days=D_infectious_days,
        )
        sim = EpidemiologySimulator(
            geometry=geom, dt_days=dt_days, duration_days=duration_days,
        )
        runs[float(R0)] = sim.simulate_sir(I0_fraction=I0_fraction)
    return runs
