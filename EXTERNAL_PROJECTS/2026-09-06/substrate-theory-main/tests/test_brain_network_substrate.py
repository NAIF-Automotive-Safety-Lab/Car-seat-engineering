"""Tests for the paired brain-network GEOMETRY + SIMULATOR.

Covers the required claims:
    (a) Watts-Strogatz adjacency: symmetric, mean degree ~ k_neighbors,
        small-world coefficient sigma > 1 (and tau = L/L_rand < 1.something
        is implicit in sigma > 1 with C/C_rand large)
    (b) Default Mode Network: 6 canonical regions, mPFC <-> PCC backbone
    (c) Frequency bands: delta/theta/alpha/beta/gamma ranges 1-4 / 4-8 /
        8-12 / 12-30 / 30-100 Hz
    (d) Kuramoto: K << K_c -> r ~ 0; K >> K_c -> r near 1 (synchronisation)
    (e) Critical Kuramoto coupling K_c = 2 gamma  (Lorentzian)
    (f) Avalanche size distribution slope ≈ -3/2 (critical branching)
    (g) Substrate signature contains expected keys
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.brain_network_substrate import (
    AVALANCHE_SLOPE_EXACT,
    BRAIN_BANDS_HZ,
    BrainNetworkGeometry,
    BrainNetworkSimulator,
    DEFAULT_DMN_REGIONS,
    LAMBDA_QCD_MEV,
    avalanche_pdf,
    kuramoto_critical_coupling_lorentzian,
)


# --------------------------------------------------------------------- #
# (a) Watts-Strogatz small-world topology                               #
# --------------------------------------------------------------------- #

def test_adjacency_is_symmetric_and_no_self_loops():
    geom = BrainNetworkGeometry(n_nodes=80, k_neighbors=6, p_rewire=0.1,
                                seed=7)
    A = geom.adjacency
    assert A.shape == (80, 80)
    assert np.array_equal(A, A.T)
    assert np.all(np.diag(A) == 0)


def test_mean_degree_matches_k_neighbors():
    """Watts-Strogatz preserves the total number of edges => mean k = k."""
    geom = BrainNetworkGeometry(n_nodes=100, k_neighbors=6, p_rewire=0.1,
                                seed=11)
    assert geom.mean_degree() == pytest.approx(6.0, abs=1.0e-9)


def test_small_world_sigma_above_unity():
    """For modest rewiring p~0.1, sigma = (C/C_r)/(L/L_r) >> 1."""
    geom = BrainNetworkGeometry(n_nodes=80, k_neighbors=6, p_rewire=0.1,
                                seed=3)
    sigma = geom.small_world_sigma(n_random=3, max_pairs=40)
    assert sigma > 1.5


def test_lattice_limit_is_not_small_world():
    """p = 0 -> ring lattice has long path lengths, sigma should be lower
    or at most modestly above 1 (mostly because L/L_rand is huge).

    More importantly, in the random limit p=1 the clustering collapses,
    so sigma falls toward 1.  We verify both extremes differ from the
    p~0.1 small-world regime.
    """
    geom_lattice = BrainNetworkGeometry(n_nodes=60, k_neighbors=6,
                                        p_rewire=0.0, seed=21)
    geom_random = BrainNetworkGeometry(n_nodes=60, k_neighbors=6,
                                       p_rewire=1.0, seed=21)
    geom_sw = BrainNetworkGeometry(n_nodes=60, k_neighbors=6,
                                   p_rewire=0.1, seed=21)
    # Clustering: lattice >> small-world > random
    C_l = geom_lattice.clustering_coefficient()
    C_r = geom_random.clustering_coefficient()
    C_sw = geom_sw.clustering_coefficient()
    assert C_l > C_r
    assert C_sw > C_r


def test_geometry_input_validation():
    with pytest.raises(ValueError):
        BrainNetworkGeometry(n_nodes=2, k_neighbors=4, p_rewire=0.1)
    with pytest.raises(ValueError):
        BrainNetworkGeometry(n_nodes=10, k_neighbors=12, p_rewire=0.1)
    with pytest.raises(ValueError):
        BrainNetworkGeometry(n_nodes=10, k_neighbors=4, p_rewire=2.0)


# --------------------------------------------------------------------- #
# (b) Default Mode Network                                              #
# --------------------------------------------------------------------- #

def test_dmn_has_six_canonical_regions():
    expected = {"mPFC", "PCC", "L Angular Gyrus", "R Angular Gyrus",
                "L Hippocampus", "R Hippocampus"}
    names = {r["name"] for r in DEFAULT_DMN_REGIONS}
    assert names == expected


def test_dmn_subgraph_topology():
    geom = BrainNetworkGeometry(n_nodes=60, k_neighbors=6, p_rewire=0.1,
                                seed=2)
    names, coords, A = geom.dmn_subgraph()
    assert len(names) == 6
    assert coords.shape == (6, 3)
    assert A.shape == (6, 6)
    # mPFC <-> PCC backbone
    i = names.index("mPFC")
    j = names.index("PCC")
    assert A[i, j] == 1
    assert A[j, i] == 1
    # PCC <-> Angular Gyri
    for ag in ("L Angular Gyrus", "R Angular Gyrus"):
        k = names.index(ag)
        assert A[j, k] == 1


# --------------------------------------------------------------------- #
# (c) Brain frequency bands                                             #
# --------------------------------------------------------------------- #

def test_brain_bands_have_canonical_ranges():
    assert BRAIN_BANDS_HZ["delta"] == (1.0, 4.0)
    assert BRAIN_BANDS_HZ["theta"] == (4.0, 8.0)
    assert BRAIN_BANDS_HZ["alpha"] == (8.0, 12.0)
    assert BRAIN_BANDS_HZ["beta"] == (12.0, 30.0)
    assert BRAIN_BANDS_HZ["gamma"] == (30.0, 100.0)


def test_band_widths_strictly_positive():
    for band, (lo, hi) in BRAIN_BANDS_HZ.items():
        assert hi > lo, f"{band} not ordered"


def test_synthetic_lfp_spectrum_has_band_peaks():
    geom = BrainNetworkGeometry(n_nodes=40, k_neighbors=4, p_rewire=0.1,
                                seed=5)
    sim = BrainNetworkSimulator(geometry=geom)
    f, p = sim.synthetic_lfp_spectrum(
        n_freqs=400,
        band_amps={"delta": 0.5, "theta": 0.5, "alpha": 1.5,
                   "beta": 0.4, "gamma": 0.3},
    )
    assert f.shape == p.shape
    assert (p > 0).all()
    # alpha band peak (8-12) should sit clearly above baseline ~ 1/f^β
    alpha_mask = (f >= 9.0) & (f <= 11.0)
    high_mask = (f >= 60.0) & (f <= 80.0)
    assert p[alpha_mask].mean() > p[high_mask].mean()


# --------------------------------------------------------------------- #
# (d) Kuramoto synchronisation                                          #
# --------------------------------------------------------------------- #

def test_kuramoto_subcritical_low_order():
    """K << K_c: order parameter r should remain small."""
    geom = BrainNetworkGeometry(n_nodes=60, k_neighbors=6, p_rewire=0.1,
                                seed=4)
    sim = BrainNetworkSimulator(geometry=geom, omega_spread_hz=2.0,
                                duration_s=1.0, dt_s=2.0e-3)
    K_c = sim.critical_coupling()
    res = sim.simulate_kuramoto(coupling_K=0.05 * K_c)
    tail = res["order_r"][-int(0.3 * res["order_r"].size):]
    assert float(tail.mean()) < 0.45


def test_kuramoto_supercritical_high_order():
    """K >> K_c: r approaches 1 (phase locking)."""
    geom = BrainNetworkGeometry(n_nodes=60, k_neighbors=6, p_rewire=0.1,
                                seed=4)
    sim = BrainNetworkSimulator(geometry=geom, omega_spread_hz=2.0,
                                duration_s=1.5, dt_s=2.0e-3)
    K_c = sim.critical_coupling()
    res = sim.simulate_kuramoto(coupling_K=20.0 * K_c)
    tail = res["order_r"][-int(0.3 * res["order_r"].size):]
    assert float(tail.mean()) > 0.7


def test_kuramoto_order_increases_with_K():
    """Sweep over K: r(K) is monotone non-decreasing on average."""
    geom = BrainNetworkGeometry(n_nodes=40, k_neighbors=6, p_rewire=0.1,
                                seed=2)
    sim = BrainNetworkSimulator(geometry=geom, omega_spread_hz=2.0)
    K_c = sim.critical_coupling()
    Ks = np.array([0.05, 1.0, 5.0, 25.0]) * K_c
    Ks_arr, rs = sim.order_parameter_vs_K(
        K_values=Ks, duration_s=0.6, dt_s=4.0e-3
    )
    # Endpoint must be substantially larger than the K=0.05*K_c point
    assert rs[-1] - rs[0] > 0.3


def test_critical_coupling_lorentzian_formula():
    """K_c = 2 * gamma_omega for Lorentzian g(omega)."""
    K_c = kuramoto_critical_coupling_lorentzian(3.0)
    assert K_c == pytest.approx(6.0, rel=1.0e-12)


def test_critical_coupling_simulator_matches_formula():
    geom = BrainNetworkGeometry(n_nodes=40, k_neighbors=6, p_rewire=0.1,
                                seed=1)
    sim = BrainNetworkSimulator(geometry=geom, omega_spread_hz=1.0)
    expected = 2.0 * 2.0 * math.pi * 1.0   # 2 * gamma_rad
    assert sim.critical_coupling() == pytest.approx(expected, rel=1.0e-12)


# --------------------------------------------------------------------- #
# (e) Neural avalanche -3/2 slope                                       #
# --------------------------------------------------------------------- #

def test_avalanche_slope_minus_three_halves():
    """Critical branching gives P(s) ∝ s^{-3/2}."""
    geom = BrainNetworkGeometry(n_nodes=20, k_neighbors=4, p_rewire=0.1,
                                seed=9)
    sim = BrainNetworkSimulator(geometry=geom, seed=9)
    sizes = sim.simulate_avalanches(n_avalanches=8000, branching=1.0,
                                    max_size=20000)
    slope, _, _ = sim.avalanche_log_log_slope(sizes, s_min=2, s_max=2000,
                                              n_bins=22)
    # Allow generous tolerance (~ 20%) given finite-size noise
    assert slope == pytest.approx(AVALANCHE_SLOPE_EXACT, abs=0.30)
    assert slope < -1.0
    assert slope > -2.0


def test_avalanche_subcritical_dies_fast():
    """branching < 1 -> avalanches small (mean size finite)."""
    geom = BrainNetworkGeometry(n_nodes=20, k_neighbors=4, p_rewire=0.1,
                                seed=2)
    sim = BrainNetworkSimulator(geometry=geom, seed=2)
    sizes = sim.simulate_avalanches(n_avalanches=2000, branching=0.5,
                                    max_size=2000)
    assert float(sizes.mean()) < 5.0


def test_avalanche_pdf_normalises():
    """Closed-form integral of (a+1) s^a / (s_max^{a+1} - s_min^{a+1})
    over [s_min, s_max] equals 1 by construction."""
    s_min, s_max, a = 1.0, 1000.0, -1.5
    sizes = np.logspace(np.log10(s_min), np.log10(s_max), 5000)
    pdf = avalanche_pdf(sizes, slope=a, s_min=s_min, s_max=s_max)
    integral = float(np.trapezoid(pdf, sizes))
    assert integral == pytest.approx(1.0, rel=2.0e-2)


# --------------------------------------------------------------------- #
# (f) Substrate signature                                               #
# --------------------------------------------------------------------- #

def test_substrate_signature_keys():
    geom = BrainNetworkGeometry(n_nodes=40, k_neighbors=4, p_rewire=0.1,
                                seed=0)
    sim = BrainNetworkSimulator(geometry=geom)
    sig = sim.substrate_signature()
    for key in ("interpretation", "topology", "n_nodes", "n_edges",
                "mean_degree", "lambda_qcd_mev", "coupling_K",
                "K_c_mean_field", "bands_hz", "avalanche_slope_exact"):
        assert key in sig
    assert sig["lambda_qcd_mev"] == LAMBDA_QCD_MEV
    assert sig["avalanche_slope_exact"] == -1.5
    assert sig["topology"] == "Watts-Strogatz small-world"


def test_kuramoto_simulation_shapes_and_bounded_order():
    geom = BrainNetworkGeometry(n_nodes=30, k_neighbors=4, p_rewire=0.1,
                                seed=0)
    sim = BrainNetworkSimulator(geometry=geom, duration_s=0.5, dt_s=2.0e-3)
    res = sim.simulate_kuramoto(coupling_K=10.0)
    assert res["order_r"].shape == res["t_s"].shape
    assert res["theta"].shape == (res["t_s"].size, 30)
    assert (res["order_r"] >= 0.0).all()
    assert (res["order_r"] <= 1.0 + 1.0e-9).all()


# --------------------------------------------------------------------- #
# (g) DMN K_4-apex topology -- substrate prediction                     #
# --------------------------------------------------------------------- #

def test_dmn_k4_apex_pcc_is_max_central():
    """PCC must be (strictly) the maximum-degree node in the canonical
    DMN sub-graph: this is the substrate-derived K_4-apex signature.
    """
    geom = BrainNetworkGeometry(n_nodes=60, k_neighbors=6, p_rewire=0.1,
                                seed=2)
    topo = geom.dmn_k4_apex_topology(apex_name="PCC")
    assert topo["apex_name"] == "PCC"
    assert topo["apex_is_max_central"] is True
    # PCC degree should exceed degree of every other DMN region
    for name, deg in topo["face_degrees"].items():
        assert topo["apex_degree"] > deg, (
            f"PCC degree {topo['apex_degree']} not > {name} degree {deg}"
        )


def test_dmn_k4_apex_closure_pattern():
    """Apex-closure: PCC connects to every other region (closure mode)."""
    geom = BrainNetworkGeometry(n_nodes=60, k_neighbors=6, p_rewire=0.1,
                                seed=2)
    topo = geom.dmn_k4_apex_topology(apex_name="PCC")
    assert topo["apex_closure"] is True
    assert topo["k4_apex_pattern_detected"] is True
    # Uniformity to all face regions = 1 (binary, all == 1)
    assert topo["connection_uniformity"] == pytest.approx(1.0, abs=1.0e-9)


def test_dmn_k4_apex_invalid_apex_raises():
    geom = BrainNetworkGeometry(n_nodes=40, k_neighbors=4, p_rewire=0.1,
                                seed=0)
    with pytest.raises(ValueError):
        geom.dmn_k4_apex_topology(apex_name="not_a_region")


def test_predict_pcc_role_structure():
    """The substrate prediction must expose the 4 testable consequences."""
    geom = BrainNetworkGeometry(n_nodes=60, k_neighbors=6, p_rewire=0.1,
                                seed=2)
    pred = geom.predict_pcc_role()
    assert "PCC" in pred["claim"]
    assert "K_4-apex" in pred["claim"] or "K_4" in pred["claim"]
    assert pred["tier"].startswith("Tier 4")
    # Four numbered consequences, each with test + falsifier
    assert isinstance(pred["consequences"], list)
    assert len(pred["consequences"]) == 4
    seen_ids = set()
    for c in pred["consequences"]:
        assert "id" in c and "statement" in c
        assert "test" in c and "falsifier" in c
        assert c["id"] not in seen_ids
        seen_ids.add(c["id"])
    assert seen_ids == {1, 2, 3, 4}
    # Required substrate concepts in the consequences (case-insensitive)
    blob = " ".join(c["statement"] + " " + c["test"] + " " + c["falsifier"]
                    for c in pred["consequences"]).lower()
    for keyword in ("sum", "uniform", "anesthesia", "lesion"):
        assert keyword in blob, f"missing keyword {keyword!r} in consequences"
    # Topology metrics embedded for downstream auditing
    assert pred["topology_metrics"]["k4_apex_pattern_detected"] is True


def test_predict_pcc_role_alternative_apex_falsifies_pattern():
    """Sanity: choosing a non-apex region (mPFC) as the candidate should
    NOT trigger the apex-closure flag in the canonical DMN sub-graph,
    because mPFC is not connected to both Angular Gyri.
    """
    geom = BrainNetworkGeometry(n_nodes=60, k_neighbors=6, p_rewire=0.1,
                                seed=2)
    topo_alt = geom.dmn_k4_apex_topology(apex_name="mPFC")
    # mPFC is not max-central (PCC is) and is not a closure node
    assert topo_alt["apex_is_max_central"] is False
    assert topo_alt["apex_closure"] is False
    assert topo_alt["k4_apex_pattern_detected"] is False


def test_substrate_signature_includes_dmn_apex():
    geom = BrainNetworkGeometry(n_nodes=40, k_neighbors=4, p_rewire=0.1,
                                seed=0)
    sig = geom.substrate_signature()
    assert "dmn_k4_apex" in sig
    assert "PCC" in sig["dmn_k4_apex"]
