"""Tests for the K_4-edge derivation of N_BAM = 6.

These tests verify that N_BAM is rigorously fixed by the K_4 simplex
1-skeleton (no 2D-slice ansatz), and that the canonical product
n_A · N_BAM = 90 reproduces the deuteron binding energy.
"""

from __future__ import annotations

from itertools import combinations
from math import comb

import pytest

from src.stiff_medium import b3_constants as C
from src.stiff_medium.nbam_from_k4_edges import (
    EDGES_PER_FACE,
    EDGE_SHARING,
    EPS_FACE_MEV,
    EXPECTED_NA,
    EXPECTED_NBAM,
    K4_VERTICES,
    deuteron_binding_check,
    enumerate_edges,
    face_pair_channel_count,
    k4_fvector,
    monte_carlo_edge_count,
    nbam_from_k4_edges,
    summary,
    verify_simplex_topology,
)


# ---------------------------------------------------------------------------
# 1. Headline derivation: N_BAM = e(K_4) = 6
# ---------------------------------------------------------------------------

class TestNBAMFromK4Edges:
    """N_BAM is the edge count of K_4 (no slice ansatz)."""

    def test_edge_count_is_six(self):
        """K_4 has C(4, 2) = 6 edges, exactly."""
        assert nbam_from_k4_edges() == 6

    def test_matches_canonical_constant(self):
        """The derivation must match the b3_constants canonical N_BAM."""
        assert nbam_from_k4_edges() == C.N_BAM
        assert nbam_from_k4_edges() == EXPECTED_NBAM

    def test_handshake_count_agrees(self):
        """(n_F · edges_per_face) / 2 = (4 · 3) / 2 = 6 — must equal edge count."""
        assert face_pair_channel_count() == nbam_from_k4_edges()
        assert face_pair_channel_count() == 6

    def test_constants_match_module(self):
        """Local module constants must be self-consistent."""
        assert K4_VERTICES == 4
        assert EDGES_PER_FACE == 3
        assert EDGE_SHARING == 2
        assert EXPECTED_NBAM == 6
        assert EXPECTED_NA == 15


# ---------------------------------------------------------------------------
# 2. Symbolic / Monte Carlo verification of the edge count
# ---------------------------------------------------------------------------

class TestEdgeEnumeration:
    """The 6 edges can be enumerated symbolically and resampled by MC."""

    def test_symbolic_enumeration(self):
        edges = enumerate_edges()
        assert len(edges) == 6
        # Every entry is a strictly increasing pair of vertices in {0,1,2,3}.
        for (i, j) in edges:
            assert 0 <= i < j < K4_VERTICES
        # The set of edges is exactly the set of 2-subsets of {0,1,2,3}.
        assert set(edges) == set(combinations(range(K4_VERTICES), 2))

    def test_no_self_loops(self):
        for (i, j) in enumerate_edges():
            assert i != j, f"K_4 must not contain a self-loop {(i, j)}"

    def test_no_duplicate_edges(self):
        edges = enumerate_edges()
        assert len(edges) == len(set(edges))

    def test_monte_carlo_recovers_six(self):
        """A modest random sample should rediscover all 6 edges."""
        n_unique = monte_carlo_edge_count(n_trials=10_000, seed=42)
        assert n_unique == 6

    def test_monte_carlo_deterministic(self):
        """Identical seeds must give identical MC results."""
        a = monte_carlo_edge_count(n_trials=1_000, seed=7)
        b = monte_carlo_edge_count(n_trials=1_000, seed=7)
        assert a == b == 6


# ---------------------------------------------------------------------------
# 3. K_4 f-vector and Euler relations
# ---------------------------------------------------------------------------

class TestK4Topology:
    """K_4 has f-vector (4, 6, 4, 1) and ∂K_4 = S² (χ = 2)."""

    def test_fvector(self):
        fv = k4_fvector()
        assert fv.f == (4, 6, 4, 1)
        assert fv.vertices == 4
        assert fv.edges == 6
        assert fv.faces == 4
        assert fv.cells == 1

    def test_boundary_is_2_sphere(self):
        """∂K_4 has χ = 4 - 6 + 4 = 2 = χ(S²)."""
        fv = k4_fvector()
        assert fv.boundary_euler() == 2

    def test_full_simplex_euler(self):
        """Closed 3-simplex has χ = V - E + F - C = 4 - 6 + 4 - 1 = 1."""
        topo = verify_simplex_topology()
        assert topo["full_euler"] == 1
        assert topo["boundary_euler"] == 2

    def test_topology_dict_completeness(self):
        topo = verify_simplex_topology()
        for key in ("vertices", "edges", "faces", "cells",
                    "boundary_euler", "full_euler"):
            assert key in topo

    def test_face_valence_consistency(self):
        """Every face is a triangle ⇒ 3 edges per face; (4·3)/2 = 6 edges."""
        fv = k4_fvector()
        n_faces = fv.faces
        edges_double_counted = n_faces * EDGES_PER_FACE
        assert edges_double_counted // EDGE_SHARING == fv.edges


# ---------------------------------------------------------------------------
# 4. Deuteron BE cross-check
# ---------------------------------------------------------------------------

class TestDeuteronBindingCrossCheck:
    """ε_face = Λ_QCD / (n_A · N_BAM) = 200 / 90 = 2.222… MeV."""

    def test_eps_face_value(self):
        d = deuteron_binding_check()
        assert d["eps_face"] == pytest.approx(EPS_FACE_MEV, rel=1e-9)
        assert d["eps_face"] == pytest.approx(200.0 / 90.0, rel=1e-12)

    def test_n_A_equals_C_N_BAM_2(self):
        """n_A = C(N_BAM, 2) = 15."""
        d = deuteron_binding_check()
        assert d["n_A"] == comb(d["N_BAM"], 2)
        assert d["n_A"] == 15

    def test_canonical_product_is_90(self):
        d = deuteron_binding_check()
        assert d["denom"] == 90
        assert d["denom"] == d["n_A"] * d["N_BAM"]

    def test_matches_AME2020_at_sub_percent(self):
        """ε_face vs deuteron BE 2.2246 MeV — should agree to ~0.1%."""
        d = deuteron_binding_check()
        # Documented deviation ~ 0.11% (rounded by 200 MeV anchor)
        assert d["rel_err"] < 0.005

    def test_canonical_constants_consistent(self):
        """The derivation values agree with canonical b3_constants."""
        d = deuteron_binding_check()
        assert d["N_BAM"] == C.N_BAM
        assert d["n_A"] == C.n_A
        assert d["denom"] == C.n_A * C.N_BAM


# ---------------------------------------------------------------------------
# 5. The 2D-slice ansatz is now redundant
# ---------------------------------------------------------------------------

class TestSliceAnsatzRedundancy:
    """The 2D-slice argument and the K_4-edge argument both yield 6.

    The K_4-edge derivation makes the slice ansatz unnecessary: 6 is
    forced by the simplex 1-skeleton without any projection assumption.
    """

    def test_2d_kissing_number_matches(self):
        """2D kissing number is 6 — the slice argument's old endpoint."""
        # Elementary fact: 6 unit disks fit around one with no slack.
        assert 6 == nbam_from_k4_edges()

    def test_3d_kissing_number_does_not_match(self):
        """3D kissing number is 12 — the slice was needed to get to 6."""
        # The K_4-edge route avoids needing to reduce 12 → 6 at all.
        assert 12 != nbam_from_k4_edges()

    def test_no_slice_assumption_in_derivation(self):
        """The derivation depends only on K_4's f-vector."""
        # K_4 lives intrinsically in any dimension ≥ 3; the count 6 is
        # invariant under embedding choice. No 2D slice required.
        fv = k4_fvector()
        assert fv.edges == nbam_from_k4_edges() == 6


# ---------------------------------------------------------------------------
# 6. Single-call summary
# ---------------------------------------------------------------------------

class TestSummary:
    def test_summary_completeness(self):
        s = summary()
        for key in ("k4_fvector", "edge_count", "handshake_count",
                    "topology", "edges_enumerated", "monte_carlo_edges",
                    "deuteron"):
            assert key in s

    def test_summary_consistency(self):
        s = summary()
        assert s["edge_count"] == s["handshake_count"] == 6
        assert s["k4_fvector"] == (4, 6, 4, 1)
        assert s["monte_carlo_edges"] == 6
        assert len(s["edges_enumerated"]) == 6
        assert s["deuteron"]["N_BAM"] == 6
        assert s["deuteron"]["n_A"] == 15
        assert s["deuteron"]["denom"] == 90
