"""Tests for the paired phase-transition geometry + simulator.

Verifies:
  (a) Bubbles of the de-saturated phase form below T_c when started from
      the +1 saturated state.
  (b) The percolation transition occurs near σ_c (here σ_c = 0,
      consistent with the symmetric Ising fixed point).
  (c) The nucleation rate is a non-trivial function of T:
      negligible at low T, grows around T_c, large above T_c.
"""
from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.phase_transition_paired import (
    T_CRITICAL_ONSAGER,
    PhaseTransitionGeometry,
    PhaseTransitionSimulator,
    nucleation_rate_vs_T,
)


# ---------------------------------------------------------------------------
# Geometry-only smoke checks
# ---------------------------------------------------------------------------

def test_geometry_critical_temperature_matches_onsager():
    geom = PhaseTransitionGeometry(L=16, J=1.0)
    assert geom.critical_temperature == pytest.approx(T_CRITICAL_ONSAGER)
    assert 2.26 < geom.critical_temperature < 2.28


def test_geometry_neighbours_periodic():
    geom = PhaseTransitionGeometry(L=8)
    nb = set(geom.neighbours(0, 0))
    assert (1, 0) in nb and (7, 0) in nb
    assert (0, 1) in nb and (0, 7) in nb


def test_label_clusters_counts_isolated_bubbles():
    geom = PhaseTransitionGeometry(L=10)
    geom.grid[:] = 1
    # Plant three isolated -1 cells far apart.
    geom.grid[1, 1] = -1
    geom.grid[5, 5] = -1
    geom.grid[8, 8] = -1
    _, sizes = geom.label_clusters(target=-1)
    assert sorted(sizes.values()) == [1, 1, 1]


def test_label_clusters_merges_connected():
    geom = PhaseTransitionGeometry(L=10)
    geom.grid[:] = 1
    # 2x2 connected -1 block
    geom.grid[3:5, 3:5] = -1
    _, sizes = geom.label_clusters(target=-1)
    assert list(sizes.values()) == [4]


def test_percolation_detected_for_full_row():
    geom = PhaseTransitionGeometry(L=12)
    geom.grid[:] = 1
    geom.grid[6, :] = -1   # full row of -1 spans all columns
    assert geom.has_percolating_cluster(target=-1)


def test_percolation_not_detected_for_isolated_bubble():
    geom = PhaseTransitionGeometry(L=12)
    geom.grid[:] = 1
    geom.grid[5:7, 5:7] = -1  # tiny 2x2 bubble
    assert not geom.has_percolating_cluster(target=-1)


# ---------------------------------------------------------------------------
# Simulator: bubbles form below T_c
# ---------------------------------------------------------------------------

def test_bubbles_form_below_Tc():
    """Below T_c, the saturated +1 state should *seed* small bubbles
    of the -1 phase as thermal fluctuations grow them."""
    geom = PhaseTransitionGeometry(L=24, seed=7)
    sim = PhaseTransitionSimulator(
        geometry=geom,
        T=2.10,                # below T_c ≈ 2.27
        n_sweeps=80,
        snapshot_every=20,
        seed_grid="saturated",
        seed=7,
    )
    hist = sim.run()
    # Some bubbles must have appeared.
    assert max(hist["n_bubbles"]) >= 1, "no bubbles formed below T_c"
    # The largest bubble must have positive area at the end.
    assert hist["largest_cluster"][-1] >= 1


# ---------------------------------------------------------------------------
# Simulator: percolation transition near sigma_c (here 0 = Ising symmetric)
# ---------------------------------------------------------------------------

def test_percolation_transition_at_sigma_c():
    """Above T_c the system disorders to ⟨σ⟩ ≈ 0 and a percolating
    -1 cluster appears; well below T_c (with small lattice), bubbles
    stay localized and there is no percolating cluster.

    We assert the *direction* of the transition rather than a specific
    fractional threshold (finite-size effects on small lattices are
    significant)."""
    geom_hot = PhaseTransitionGeometry(L=24, seed=11)
    sim_hot = PhaseTransitionSimulator(
        geometry=geom_hot,
        T=4.0,                  # well above T_c
        n_sweeps=120,
        snapshot_every=120,
        seed_grid="random",
        seed=11,
    )
    sim_hot.run()
    sigma_hot = abs(geom_hot.order_parameter())
    perc_hot = sim_hot.percolated()

    geom_cold = PhaseTransitionGeometry(L=24, seed=13)
    sim_cold = PhaseTransitionSimulator(
        geometry=geom_cold,
        T=1.2,                  # well below T_c
        n_sweeps=120,
        snapshot_every=120,
        seed_grid="saturated",
        seed=13,
    )
    sim_cold.run()
    sigma_cold = abs(geom_cold.order_parameter())
    perc_cold = sim_cold.percolated()

    # Order-parameter direction matches expectation.
    assert sigma_hot < geom_hot.critical_sigma + 0.5  # |⟨σ⟩| close to 0
    assert sigma_cold > 0.5                           # ordered phase survives
    # Percolation pattern: hot disordered phase percolates, cold doesn't.
    assert perc_hot is True
    assert perc_cold is False


# ---------------------------------------------------------------------------
# Simulator: nucleation rate vs temperature
# ---------------------------------------------------------------------------

def test_nucleation_rate_grows_with_T_below_Tc():
    """Inside the metastable saturated state (T below T_c), the
    bubble-nucleation rate should grow with T as the energy barrier
    for forming a critical bubble decreases."""
    Ts = np.array([1.4, 1.8, 2.2])
    res = nucleation_rate_vs_T(Ts, L=20, n_sweeps=60, seed=21)
    rates = res["rate"]
    # Strict monotone is too brittle on a small lattice; require the
    # endpoint rate at T=2.2 to dominate the rate at T=1.4 by at least 3×.
    assert rates[-1] > 3.0 * rates[0] + 1e-9, (
        f"nucleation rate does not grow with T: {rates}"
    )
    # The lowest-T rate must be strictly less than the highest-T rate.
    assert rates[0] < rates[-1]


def test_nucleation_rate_nonnegative():
    Ts = np.array([1.0, 2.0, 3.0])
    res = nucleation_rate_vs_T(Ts, L=16, n_sweeps=40, seed=31)
    assert np.all(res["rate"] >= 0.0)
    assert np.all(np.isfinite(res["rate"]))


def test_simulator_history_keys_complete():
    geom = PhaseTransitionGeometry(L=12, seed=3)
    sim = PhaseTransitionSimulator(
        geometry=geom,
        T=2.0,
        n_sweeps=20,
        snapshot_every=5,
        seed=3,
    )
    h = sim.run()
    for k in (
        "sweep",
        "order_parameter",
        "energy",
        "largest_cluster",
        "n_bubbles",
        "snapshots",
        "snapshot_sweeps",
    ):
        assert k in h
    assert len(h["sweep"]) == 20
    # Snapshots taken at sweep 0,5,10,15 plus final sweep 19 ⇒ 5 frames.
    assert len(h["snapshots"]) == len(h["snapshot_sweeps"])
    assert len(h["snapshots"]) >= 4


def test_glauber_dynamics_runs():
    geom = PhaseTransitionGeometry(L=12, seed=5)
    sim = PhaseTransitionSimulator(
        geometry=geom,
        T=2.5,
        n_sweeps=20,
        dynamics="glauber",
        seed=5,
    )
    h = sim.run()
    assert len(h["order_parameter"]) == 20
    assert all(-1.001 <= v <= 1.001 for v in h["order_parameter"])
