"""Paired brain-network dynamics GEOMETRY + SIMULATOR for the substrate
framework.

Models large-scale brain network activity as a substrate phenomenon.  Brain
regions are substrate cells coupled by white-matter tracts; cortical
oscillations are coherent strain-rate modes propagating across the network.
Synchronisation is described by the Kuramoto (1975) model on a graph; the
graph topology is small-world (Watts & Strogatz 1998), as documented for
human and macaque connectomes (Sporns 2011; Bassett & Bullmore 2017).

GEOMETRY
--------
``BrainNetworkGeometry`` defines

    * ``nodes``      -- a list of brain regions (with 3D coordinates so that
                        plots have a usable schematic layout)
    * ``edges``      -- white-matter tracts; the adjacency matrix is
                        constructed by a Watts-Strogatz small-world model
                        (``k_neighbors`` ring + rewiring probability ``p``)
    * a *Default Mode Network* sub-graph: medial PFC, posterior cingulate
      cortex, angular gyrus (L/R), hippocampus (L/R) — the canonical
      task-negative network of Raichle 2001 / Buckner 2008
    * Watts-Strogatz small-world coefficient

        sigma = (C / C_rand) / (L / L_rand)                              (1)

      A network is small-world when sigma > 1 (large clustering, short
      path length).  An additional measure (Telesford 2011) is

        omega = L_rand / L  -  C / C_lattice                              (2)

      with -1 ≤ omega ≤ 1 and |omega| ≪ 1 indicating small-worldness.

SIMULATOR
---------
``BrainNetworkSimulator`` implements

    * Kuramoto coupled phase oscillators on the geometry's adjacency
      matrix:
                d θ_i           K
                ----- = ω_i +  --- Σ_j A_ij sin(θ_j - θ_i)
                  dt            N
      with the Kuramoto order parameter
                r(t) = | (1/N) Σ_j exp(i θ_j(t)) |                       (3)
      and the mean-field critical coupling
                K_c = 2 / (π g(ω̄))                                       (4)
      for natural-frequency density g(ω).
    * Brain rhythm bands (Buzsáki 2006):
        delta 1-4 Hz, theta 4-8 Hz, alpha 8-12 Hz, beta 12-30 Hz,
        gamma 30-100 Hz.
    * 1/f^β power spectrum baseline (Voytek & Knight 2015) with band peaks.
    * Critical-brain neural avalanche size distribution
                P(s) ∝ s^{-3/2}                                          (5)
      (Beggs & Plenz 2003) — the slope is *exactly* -3/2 for a critical
      branching process (mean branching ratio sigma ≈ 1).

Substrate interpretation
------------------------
* node = mesoscopic substrate patch (population of ~10^4 - 10^5 neurons)
* edge = white-matter tract = bulk strain-coupling channel between patches
* phase θ_i = local oscillator phase of the patch's strain field
* coupling K = effective stiffness of the inter-patch substrate medium
* small-world topology = mixed local (lattice) + global (long-range) strain
  coupling, the optimum of cost vs information transfer
* alpha rhythm (8-12 Hz) ≈ 100 ms ≈ default substrate relaxation time
  for the cortical scale; gamma (30-100 Hz) is the high-K limit
* avalanches = saturation cascades on the substrate (σ ≤ ½ critical line)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import math

import numpy as np


# ---------------------------------------------------------------------------
# Substrate / framework constants
# ---------------------------------------------------------------------------

LAMBDA_QCD_MEV: float = 200.0


# ---------------------------------------------------------------------------
# Reference values from neuroscience (Buzsáki 2006; Buckner 2008;
# Beggs & Plenz 2003)
# ---------------------------------------------------------------------------

# Frequency bands (Hz)  -- (low, high)
BRAIN_BANDS_HZ: Dict[str, Tuple[float, float]] = {
    "delta": (1.0, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 12.0),
    "beta":  (12.0, 30.0),
    "gamma": (30.0, 100.0),
}

# Neural-avalanche slope (Beggs & Plenz 2003 / branching critical exponent)
AVALANCHE_SLOPE_EXACT: float = -1.5

# Default Watts-Strogatz parameters that give small-world (Bassett 2017)
WS_K_NEIGHBORS_DEFAULT: int = 6      # ring lattice degree
WS_P_REWIRE_DEFAULT: float = 0.10    # 0 = lattice, 1 = random


# ---------------------------------------------------------------------------
# Closed-form network metrics
# ---------------------------------------------------------------------------

def kuramoto_critical_coupling_lorentzian(gamma_omega: float) -> float:
    """Mean-field Kuramoto K_c for a Lorentzian g(ω) of half-width gamma.

    Standard result (Kuramoto 1984; Strogatz 2000):

        K_c = 2 * gamma_omega
    """
    if gamma_omega <= 0.0:
        raise ValueError("gamma_omega must be > 0")
    return 2.0 * gamma_omega


def avalanche_pdf(sizes: np.ndarray, slope: float = AVALANCHE_SLOPE_EXACT,
                  s_min: float = 1.0, s_max: Optional[float] = None
                  ) -> np.ndarray:
    """Power-law PDF P(s) ∝ s^slope on [s_min, s_max] (or s_min..max(sizes))."""
    s = np.asarray(sizes, dtype=float)
    s = np.where(s < s_min, s_min, s)
    if s_max is None:
        s_max = float(s.max())
    if abs(slope + 1.0) < 1.0e-12:
        norm = math.log(s_max / s_min)
        return (1.0 / norm) / s
    a = slope
    norm = (s_max ** (a + 1.0) - s_min ** (a + 1.0)) / (a + 1.0)
    return (s ** a) / norm


# ---------------------------------------------------------------------------
# Default Mode Network -- canonical region list
# ---------------------------------------------------------------------------

# Roughly normalised MNI-style coordinates (xyz in arbitrary units), one
# entry per region.  Layout is for visualization, not anatomical accuracy.
DEFAULT_DMN_REGIONS: List[Dict] = [
    {"name": "mPFC",            "xyz": ( 0.00,  0.95,  0.20)},
    {"name": "PCC",             "xyz": ( 0.00, -0.85,  0.30)},
    {"name": "L Angular Gyrus", "xyz": (-0.85, -0.20,  0.20)},
    {"name": "R Angular Gyrus", "xyz": ( 0.85, -0.20,  0.20)},
    {"name": "L Hippocampus",   "xyz": (-0.40, -0.10, -0.50)},
    {"name": "R Hippocampus",   "xyz": ( 0.40, -0.10, -0.50)},
]


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class BrainNetworkGeometry:
    """Whole-brain network: nodes = regions, edges = white-matter tracts.

    The default geometry is a Watts-Strogatz small-world graph with N
    nodes, ring degree ``k_neighbors``, and rewiring probability
    ``p_rewire``.  A canonical Default Mode Network sub-graph is
    constructed in addition, accessible via :attr:`dmn_regions` and
    :meth:`dmn_subgraph`.
    """

    n_nodes: int = 90              # ~Brodmann/AAL atlas count
    k_neighbors: int = WS_K_NEIGHBORS_DEFAULT
    p_rewire: float = WS_P_REWIRE_DEFAULT
    seed: Optional[int] = 42

    # Optional explicit DMN region list
    dmn_regions: List[Dict] = field(
        default_factory=lambda: [dict(r) for r in DEFAULT_DMN_REGIONS]
    )

    # Filled by ``__post_init__``
    adjacency: np.ndarray = field(init=False, repr=False)
    coords: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.n_nodes < 4:
            raise ValueError("n_nodes must be >= 4 for a meaningful graph")
        if self.k_neighbors < 2 or self.k_neighbors >= self.n_nodes:
            raise ValueError("k_neighbors must be in [2, n_nodes-1]")
        if not (0.0 <= self.p_rewire <= 1.0):
            raise ValueError("p_rewire must be in [0, 1]")

        rng = np.random.default_rng(self.seed)
        self.adjacency = self._watts_strogatz(
            self.n_nodes, self.k_neighbors, self.p_rewire, rng,
        )
        self.coords = self._spherical_layout(self.n_nodes)

    # ------------------------------------------------------------------
    # Watts-Strogatz construction
    # ------------------------------------------------------------------

    @staticmethod
    def _watts_strogatz(n: int, k: int, p: float,
                        rng: np.random.Generator) -> np.ndarray:
        """Construct the symmetric adjacency matrix of a WS small-world."""
        if k % 2 != 0:
            k = k + 1   # Watts-Strogatz requires even k (k/2 each side)
        A = np.zeros((n, n), dtype=np.uint8)
        # Step 1: ring lattice -- each node connected to k/2 on each side
        for i in range(n):
            for j in range(1, k // 2 + 1):
                A[i, (i + j) % n] = 1
                A[(i + j) % n, i] = 1
        # Step 2: rewire each edge with probability p, avoiding self-loops
        # and duplicate edges.  Iterate over the ring (i, (i+j) mod n).
        for i in range(n):
            for j in range(1, k // 2 + 1):
                if rng.random() < p:
                    nbr = (i + j) % n
                    # Pick a new target uniformly among non-neighbours of i
                    candidates = np.where(
                        (A[i] == 0) & (np.arange(n) != i)
                    )[0]
                    if candidates.size > 0:
                        new = int(rng.choice(candidates))
                        A[i, nbr] = 0
                        A[nbr, i] = 0
                        A[i, new] = 1
                        A[new, i] = 1
        return A

    @staticmethod
    def _spherical_layout(n: int) -> np.ndarray:
        """Approximately even points on a unit sphere (Fibonacci lattice)."""
        coords = np.zeros((n, 3))
        phi = math.pi * (3.0 - math.sqrt(5.0))   # golden angle
        for i in range(n):
            y = 1.0 - (i / max(1.0, n - 1.0)) * 2.0    # in [-1, 1]
            r = math.sqrt(max(0.0, 1.0 - y * y))
            theta = phi * i
            coords[i, 0] = math.cos(theta) * r
            coords[i, 1] = y
            coords[i, 2] = math.sin(theta) * r
        return coords

    # ------------------------------------------------------------------
    # Network metrics
    # ------------------------------------------------------------------

    def degree(self) -> np.ndarray:
        return self.adjacency.sum(axis=1)

    def mean_degree(self) -> float:
        return float(self.degree().mean())

    def edges(self) -> List[Tuple[int, int]]:
        """Undirected edge list (i < j)."""
        out: List[Tuple[int, int]] = []
        n = self.n_nodes
        for i in range(n):
            for j in range(i + 1, n):
                if self.adjacency[i, j]:
                    out.append((i, j))
        return out

    def n_edges(self) -> int:
        return int(self.adjacency.sum() // 2)

    def clustering_coefficient(self) -> float:
        """Mean local clustering coefficient C (Watts-Strogatz definition)."""
        A = self.adjacency
        n = A.shape[0]
        coeffs = np.zeros(n)
        for i in range(n):
            nbrs = np.where(A[i] > 0)[0]
            k_i = nbrs.size
            if k_i < 2:
                coeffs[i] = 0.0
                continue
            # Number of triangles through i
            sub = A[np.ix_(nbrs, nbrs)]
            triangles = sub.sum() / 2.0
            coeffs[i] = (2.0 * triangles) / (k_i * (k_i - 1))
        return float(coeffs.mean())

    def average_shortest_path_length(self,
                                     max_pairs: Optional[int] = None
                                     ) -> float:
        """Mean shortest path length L (BFS).  Skips disconnected pairs.

        ``max_pairs`` optionally caps the number of source nodes used to
        keep cost bounded (BFS from each is O(N + E)).  Default: all nodes.
        """
        A = self.adjacency
        n = A.shape[0]
        sources = np.arange(n)
        if max_pairs is not None and max_pairs < n:
            rng = np.random.default_rng(self.seed)
            sources = rng.choice(n, size=max_pairs, replace=False)
        total = 0.0
        count = 0
        for s in sources:
            # BFS
            dist = -np.ones(n, dtype=int)
            dist[s] = 0
            queue: List[int] = [int(s)]
            while queue:
                u = queue.pop(0)
                for v in np.where(A[u] > 0)[0]:
                    if dist[v] < 0:
                        dist[v] = dist[u] + 1
                        queue.append(int(v))
            mask = (dist > 0)
            total += float(dist[mask].sum())
            count += int(mask.sum())
        if count == 0:
            return float("inf")
        return total / count

    def small_world_sigma(self, n_random: int = 3,
                          max_pairs: Optional[int] = 80) -> float:
        """Watts-Strogatz small-world coefficient sigma = (C/C_r)/(L/L_r).

        Compares to Erdős-Rényi random graphs of equal size and density.
        """
        C = self.clustering_coefficient()
        L = self.average_shortest_path_length(max_pairs=max_pairs)

        n = self.n_nodes
        # Edge probability matched to current density
        p_er = 2.0 * self.n_edges() / max(1.0, n * (n - 1))

        Cs: List[float] = []
        Ls: List[float] = []
        for k in range(n_random):
            rng = np.random.default_rng((self.seed or 0) + 1009 + k)
            A_rand = (rng.random((n, n)) < p_er).astype(np.uint8)
            np.fill_diagonal(A_rand, 0)
            A_rand = np.maximum(A_rand, A_rand.T)
            geom_r = BrainNetworkGeometry.__new__(BrainNetworkGeometry)
            geom_r.n_nodes = n
            geom_r.k_neighbors = self.k_neighbors
            geom_r.p_rewire = 1.0
            geom_r.seed = (self.seed or 0) + 1009 + k
            geom_r.dmn_regions = []
            geom_r.adjacency = A_rand
            geom_r.coords = self.coords
            Cs.append(geom_r.clustering_coefficient())
            Ls.append(geom_r.average_shortest_path_length(max_pairs=max_pairs))
        C_rand = max(float(np.mean(Cs)), 1.0e-6)
        L_rand = max(float(np.mean(Ls)), 1.0e-6)
        return (C / C_rand) / (L / L_rand)

    # ------------------------------------------------------------------
    # Default Mode Network sub-graph
    # ------------------------------------------------------------------

    def dmn_subgraph(self) -> Tuple[List[str], np.ndarray, np.ndarray]:
        """Return (region_names, coords[K,3], adjacency[K,K]) for the DMN.

        Connectivity is the canonical DMN topology:  mPFC <-> PCC is the
        midline backbone; PCC <-> Angular Gyrus L/R; mPFC <-> Hippocampus
        L/R; PCC <-> Hippocampus L/R; AG-L <-> AG-R lateral coupling.
        """
        regs = self.dmn_regions
        names = [r["name"] for r in regs]
        coords = np.array([r["xyz"] for r in regs])
        K = len(regs)
        idx = {n: i for i, n in enumerate(names)}
        A = np.zeros((K, K), dtype=np.uint8)

        def link(a: str, b: str) -> None:
            if a in idx and b in idx:
                i, j = idx[a], idx[b]
                A[i, j] = 1
                A[j, i] = 1

        link("mPFC", "PCC")
        link("PCC", "L Angular Gyrus")
        link("PCC", "R Angular Gyrus")
        link("mPFC", "L Hippocampus")
        link("mPFC", "R Hippocampus")
        link("PCC", "L Hippocampus")
        link("PCC", "R Hippocampus")
        link("L Angular Gyrus", "R Angular Gyrus")
        link("L Angular Gyrus", "L Hippocampus")
        link("R Angular Gyrus", "R Hippocampus")
        return names, coords, A

    # ------------------------------------------------------------------
    # Default Mode Network — K_4-apex topology (substrate prediction)
    # ------------------------------------------------------------------

    def dmn_k4_apex_topology(self,
                             apex_name: str = "PCC",
                             ) -> Dict:
        """Identify the K_4-apex closure mode of the Default Mode Network.

        Substrate framework prediction.  The K_4 nucleon has 3 quark-vertex
        face + 1 apex closure vertex; the Waddington landscape has
        1 pluripotent + 3 germ-layer differentiations.  Both are "1+3"
        substrate structures.  The DMN is asserted to obey the same
        K_4-apex partition: PCC (posterior cingulate cortex) is the apex,
        the other 5 regions span a closed face-cluster, and PCC closes
        the network as the sum/closure mode.

        This method computes degree-based centrality for each DMN region
        and reports whether the apex carries strictly the highest score
        (the substrate-derived signature).  It also returns a ``closure``
        flag = True iff the apex's connectivity profile is the algebraic
        sum (uniform projection) of the face cluster — i.e. the apex
        connects to every other region with non-zero weight.

        Returns
        -------
        Dict with keys:
            apex_name            : str, the apex region (default "PCC")
            apex_degree          : int, degree of the apex node
            face_degrees         : Dict[str, int], degree per non-apex
            apex_centrality      : float, normalized degree centrality of apex
            face_centralities    : Dict[str, float], same for non-apex regions
            apex_is_max_central  : bool, True iff apex strictly ≥ all others
            apex_closure         : bool, True iff apex connects to every
                                   non-apex region (sum/closure pattern)
            connection_uniformity: float in [0,1], 1 = uniform connectivity
                                   to all face regions, 0 = single-target
        """
        names, _, A = self.dmn_subgraph()
        if apex_name not in names:
            raise ValueError(
                f"apex_name {apex_name!r} not in DMN regions {names}")
        K = len(names)
        idx = {n: i for i, n in enumerate(names)}
        a = idx[apex_name]
        face_idx = [i for i in range(K) if i != a]
        face_names = [names[i] for i in face_idx]

        # Degrees in the DMN sub-graph
        deg = A.sum(axis=1)
        apex_deg = int(deg[a])
        face_deg = {n: int(deg[idx[n]]) for n in face_names}
        # Degree centrality (normalized by K-1 max possible neighbours)
        norm = max(K - 1, 1)
        apex_cent = apex_deg / norm
        face_cent = {n: face_deg[n] / norm for n in face_names}

        apex_is_max = all(apex_deg >= face_deg[n] for n in face_names) \
            and any(apex_deg > face_deg[n] for n in face_names)

        # Closure: apex connects to every non-apex region
        connections_to_face = [int(A[a, j]) for j in face_idx]
        apex_closure = all(c == 1 for c in connections_to_face)

        # Uniformity: 1 - normalized variance of apex-to-face weights
        if len(connections_to_face) > 1 and sum(connections_to_face) > 0:
            mean_w = float(np.mean(connections_to_face))
            var_w = float(np.var(connections_to_face))
            # For binary {0,1} weights, max variance = 0.25 at 50/50 mix
            uniformity = 1.0 - min(var_w / 0.25, 1.0)
        else:
            uniformity = 0.0

        return {
            "apex_name": apex_name,
            "apex_degree": apex_deg,
            "face_degrees": face_deg,
            "apex_centrality": float(apex_cent),
            "face_centralities": face_cent,
            "apex_is_max_central": bool(apex_is_max),
            "apex_closure": bool(apex_closure),
            "connection_uniformity": float(uniformity),
            "k4_apex_pattern_detected": bool(apex_is_max and apex_closure),
        }

    def predict_pcc_role(self) -> Dict:
        """Substrate-framework prediction: PCC is the K_4-apex closure mode.

        Cross-domain alignment:
            * K_4 nucleon: 3 quark-vertex face + 1 apex closure vertex
            * Waddington landscape: 1 pluripotent state + 3 germ-layer
              differentiations (endoderm, mesoderm, ectoderm)
            * Qubit triality on the Bloch sphere: 3 mutually-unbiased
              bases (X, Y, Z) + 1 "no-orientation" closure
            * DMN: 1 PCC apex + 5 regional face-cluster (extended K_4 because
              the DMN is doubled by L/R lateralization on Angular Gyrus and
              Hippocampus, but the substrate-irreducible apex/face split is
              preserved)

        Returns a structured prediction with 4 testable consequences and
        their corresponding falsifiers.
        """
        topology = self.dmn_k4_apex_topology(apex_name="PCC")
        return {
            "claim": "PCC is the K_4-apex closure mode of the DMN",
            "substrate_basis": (
                "K_4 simplex topology with 1+N apex/face partition. "
                "Cross-domain reuse: K_4 nucleon, Waddington landscape, "
                "qubit triality. DMN inherits the same apex/face split."
            ),
            "topology_metrics": topology,
            "consequences": [
                {
                    "id": 1,
                    "statement": (
                        "PCC dynamics ≈ sum/closure of other 5 regions "
                        "(substrate-derived: apex closes the face)"
                    ),
                    "test": (
                        "Linear regression PCC(t) ~ Σ_i w_i · region_i(t) "
                        "should explain >70% variance with non-degenerate "
                        "(non-zero) weights on every face region"
                    ),
                    "falsifier": (
                        "PCC fits as independent oscillator with low "
                        "explained variance from face regions (<30%)"
                    ),
                },
                {
                    "id": 2,
                    "statement": (
                        "PCC connectivity uniformly distributed across "
                        "other 5 regions, not preferentially to subsets"
                    ),
                    "test": (
                        "Variance of PCC→region connection weights "
                        "(diffusion MRI tractography) is small relative "
                        "to mean: CV < 0.4"
                    ),
                    "falsifier": (
                        "PCC connectivity is concentrated on ≤2 of the "
                        "5 face regions (CV > 1.0)"
                    ),
                },
                {
                    "id": 3,
                    "statement": (
                        "PCC is LAST to lose activity in anesthesia "
                        "(loss of apex dissolves the bound network)"
                    ),
                    "test": (
                        "Track regional BOLD/EEG amplitudes during "
                        "propofol/sevoflurane induction; PCC should "
                        "remain coherent latest into deep anesthesia"
                    ),
                    "falsifier": (
                        "PCC loses coherence before any of the 5 face "
                        "regions during anesthesia induction"
                    ),
                },
                {
                    "id": 4,
                    "statement": (
                        "PCC lesions cause disproportionate disruption "
                        "of self-referential / autobiographical processing"
                    ),
                    "test": (
                        "Lesion-mapping studies (PCC vs equivalent-volume "
                        "mPFC or angular-gyrus lesions) show >2x deficit "
                        "on self-referential tasks for PCC lesions"
                    ),
                    "falsifier": (
                        "PCC lesions produce deficits comparable to or "
                        "smaller than other DMN-region lesions"
                    ),
                },
            ],
            "data_source": (
                "Human Connectome Project (HCP) S1200 release for "
                "tractography + resting-state BOLD; ds002785 "
                "(propofol-anesthesia BOLD) for consequence 3"
            ),
            "tier": "Tier 4 falsifiable substrate prediction",
        }

    # ------------------------------------------------------------------
    # Substrate signature
    # ------------------------------------------------------------------

    def substrate_signature(self) -> Dict:
        return {
            "interpretation": (
                "node = mesoscopic substrate patch; edge = white-matter "
                "tract = bulk strain-coupling channel; phase θ_i = local "
                "strain-field oscillator phase"
            ),
            "topology": "Watts-Strogatz small-world",
            "n_nodes": self.n_nodes,
            "n_edges": self.n_edges(),
            "mean_degree": self.mean_degree(),
            "lambda_qcd_mev": LAMBDA_QCD_MEV,
            "dmn_k4_apex": "PCC (posterior cingulate cortex)",
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class BrainNetworkSimulator:
    """Kuramoto + spectral + avalanche dynamics on a BrainNetworkGeometry."""

    geometry: BrainNetworkGeometry
    coupling_K: float = 1.0
    omega_center_hz: float = 10.0   # alpha-band center
    omega_spread_hz: float = 1.0    # half-width for Lorentzian g(ω)
    dt_s: float = 1.0e-3            # 1 ms integration step
    duration_s: float = 2.0
    seed: Optional[int] = 7

    # ------------------------------------------------------------------
    # Natural frequency draw and Kuramoto critical coupling
    # ------------------------------------------------------------------

    def natural_frequencies_rad_per_s(self) -> np.ndarray:
        """Draw N natural angular frequencies from a Lorentzian(ω0, gamma)."""
        rng = np.random.default_rng(self.seed)
        n = self.geometry.n_nodes
        u = rng.random(n)
        # Inverse-CDF sampling for Cauchy/Lorentzian
        f_hz = (self.omega_center_hz
                + self.omega_spread_hz * np.tan(math.pi * (u - 0.5)))
        return 2.0 * math.pi * f_hz

    def critical_coupling(self) -> float:
        """Mean-field K_c for the Lorentzian g(ω) (in rad/s units)."""
        gamma_rad = 2.0 * math.pi * self.omega_spread_hz
        return kuramoto_critical_coupling_lorentzian(gamma_rad)

    # ------------------------------------------------------------------
    # Kuramoto integrator
    # ------------------------------------------------------------------

    def simulate_kuramoto(self,
                          coupling_K: Optional[float] = None,
                          duration_s: Optional[float] = None,
                          dt_s: Optional[float] = None,
                          ) -> Dict:
        """Integrate Kuramoto dθ_i/dt = ω_i + (K/<k>) Σ A_ij sin(θ_j-θ_i).

        The K/<k> normalisation is the standard rescaling that lets the
        mean-field K_c keep its physical meaning on a sparse graph
        (Restrepo, Ott & Hunt 2005).  In the all-to-all limit A_ij = 1
        and <k> = N - 1, recovering the original Kuramoto K/N form.
        """
        K = float(self.coupling_K if coupling_K is None else coupling_K)
        T = float(self.duration_s if duration_s is None else duration_s)
        dt = float(self.dt_s if dt_s is None else dt_s)

        n = self.geometry.n_nodes
        A = self.geometry.adjacency.astype(float)
        omega = self.natural_frequencies_rad_per_s()

        rng = np.random.default_rng((self.seed or 0) + 31)
        theta = rng.uniform(-math.pi, math.pi, size=n)

        n_steps = int(round(T / dt)) + 1
        ts = np.linspace(0.0, T, n_steps)
        order_r = np.zeros(n_steps)
        theta_hist = np.zeros((n_steps, n))

        mean_k = max(self.geometry.mean_degree(), 1.0)
        kn = K / mean_k
        for s_idx in range(n_steps):
            # Order parameter
            z = np.exp(1j * theta).mean()
            order_r[s_idx] = float(np.abs(z))
            theta_hist[s_idx] = theta
            # RK1 / Euler step (sufficient at dt=1ms for these dynamics)
            sin_diff = np.sin(theta[None, :] - theta[:, None])
            interaction = (A * sin_diff).sum(axis=1)
            theta = theta + dt * (omega + kn * interaction)
            theta = (theta + math.pi) % (2.0 * math.pi) - math.pi

        return {
            "t_s": ts,
            "theta": theta_hist,
            "order_r": order_r,
            "omega_rad_s": omega,
            "coupling_K": K,
            "K_c": self.critical_coupling(),
        }

    def order_parameter_vs_K(self,
                             K_values: Sequence[float],
                             tail_fraction: float = 0.4,
                             duration_s: float = 1.0,
                             dt_s: float = 5.0e-3,
                             ) -> Tuple[np.ndarray, np.ndarray]:
        """Sweep K and return mean steady-state r(K)."""
        Ks = np.asarray(list(K_values), dtype=float)
        rs = np.zeros_like(Ks)
        for i, K in enumerate(Ks):
            res = self.simulate_kuramoto(coupling_K=float(K),
                                         duration_s=duration_s,
                                         dt_s=dt_s)
            n = res["order_r"].size
            tail = res["order_r"][int((1.0 - tail_fraction) * n):]
            rs[i] = float(tail.mean())
        return Ks, rs

    # ------------------------------------------------------------------
    # Brain rhythm bands and 1/f spectrum
    # ------------------------------------------------------------------

    def band_ranges_hz(self) -> Dict[str, Tuple[float, float]]:
        return dict(BRAIN_BANDS_HZ)

    def synthetic_lfp_spectrum(self,
                               n_freqs: int = 600,
                               beta: float = 1.5,
                               band_amps: Optional[Dict[str, float]] = None,
                               band_width_hz: float = 1.5,
                               ) -> Tuple[np.ndarray, np.ndarray]:
        """Synthetic broadband 1/f^β + 5 band peaks (Lorentzians).

        Returns (frequencies_hz, power_arb).  Useful as a stand-in for a
        single-channel local-field potential power spectrum without
        committing to a full microscopic neuron model.
        """
        f = np.linspace(0.5, 120.0, n_freqs)
        power = 1.0 / (f ** beta)
        if band_amps is None:
            band_amps = {b: 0.6 for b in BRAIN_BANDS_HZ}
        for band, (lo, hi) in BRAIN_BANDS_HZ.items():
            f0 = 0.5 * (lo + hi)
            amp = float(band_amps.get(band, 0.0))
            power = power + amp * (band_width_hz ** 2) / (
                (f - f0) ** 2 + band_width_hz ** 2
            )
        return f, power

    # ------------------------------------------------------------------
    # Critical avalanches (branching process)
    # ------------------------------------------------------------------

    def simulate_avalanches(self,
                            n_avalanches: int = 4000,
                            branching: float = 1.0,
                            max_size: int = 10000,
                            ) -> np.ndarray:
        """Galton-Watson critical branching (sigma = branching).

        Each ancestor produces a Poisson(branching) number of descendants;
        the avalanche terminates when the population dies out (or hits
        ``max_size`` as a safety cap).  At branching = 1 (critical) the
        avalanche-size distribution follows P(s) ∝ s^{-3/2}.
        """
        rng = np.random.default_rng((self.seed or 0) + 991)
        sizes = np.zeros(n_avalanches, dtype=int)
        for k in range(n_avalanches):
            pop = 1
            total = 1
            while pop > 0 and total < max_size:
                children = int(rng.poisson(branching * pop))
                pop = children
                total += pop
            sizes[k] = total
        return sizes

    @staticmethod
    def avalanche_log_log_slope(sizes: np.ndarray,
                                s_min: int = 2,
                                s_max: Optional[int] = None,
                                n_bins: int = 25,
                                ) -> Tuple[float, np.ndarray, np.ndarray]:
        """Fit a log-log slope to the binned avalanche-size PDF.

        Returns (slope, bin_centers, pdf).  Logarithmic binning is used.
        """
        s = np.asarray(sizes, dtype=int)
        s = s[s >= s_min]
        if s_max is None:
            s_max = int(s.max()) if s.size > 0 else s_min
        if s.size < 5 or s_max <= s_min:
            return 0.0, np.array([]), np.array([])
        edges = np.logspace(math.log10(s_min), math.log10(s_max), n_bins + 1)
        counts, edges = np.histogram(s, bins=edges)
        widths = np.diff(edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        pdf = counts / (s.size * np.maximum(widths, 1.0e-12))
        keep = (counts > 0) & (centers > 0)
        if int(keep.sum()) < 3:
            return 0.0, centers, pdf
        x = np.log10(centers[keep])
        y = np.log10(pdf[keep])
        # Linear regression slope
        slope, _ = np.polyfit(x, y, 1)
        return float(slope), centers, pdf

    # ------------------------------------------------------------------
    # Substrate signature
    # ------------------------------------------------------------------

    def substrate_signature(self) -> Dict:
        sig = self.geometry.substrate_signature()
        sig.update({
            "coupling_K": self.coupling_K,
            "omega_center_hz": self.omega_center_hz,
            "omega_spread_hz": self.omega_spread_hz,
            "K_c_mean_field": self.critical_coupling(),
            "bands_hz": dict(BRAIN_BANDS_HZ),
            "avalanche_slope_exact": AVALANCHE_SLOPE_EXACT,
        })
        return sig
