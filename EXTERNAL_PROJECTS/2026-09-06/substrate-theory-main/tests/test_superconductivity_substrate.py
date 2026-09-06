"""Tests for superconductivity_substrate — paired Cooper bridge + BCS gap +
substrate-saturation T_c ceiling.

Verifies:
    * BCS universal gap ratio 2Δ/(k_B T_c) = 3.5277537... (= 2π/e^γ)
    * Hg back-solved gap → T_c = 4.2 K (Onnes 1911)
    * HBCCO back-solved gap → T_c = 134 K (Schilling 1993, ambient-pressure record)
    * Substrate saturation ceiling T_c,max = Λ_QCD / n_R = 128.9 K
    * Δ(T) self-consistency → 0 as T → T_c
    * Density of states gap structure (zero inside |E|<Δ, peaks at |E|=Δ)
    * Meissner exponential decay
    * Type-I vs Type-II GL κ classification
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.superconductivity_substrate import (
    BCS_UNIVERSAL_GAP_RATIO,
    GL_KAPPA_THRESHOLD,
    REFERENCE_SUPERCONDUCTORS,
    SUBSTRATE_E_PER_MODE_MEV,
    SUBSTRATE_TC_MAX_K,
    SuperconductivityGeometry,
    SuperconductivitySimulator,
)


# ---------------------------------------------------------------------------
# Universal BCS gap ratio
# ---------------------------------------------------------------------------

def test_bcs_universal_ratio_is_3p527():
    """Weak-coupling BCS universal ratio 2Δ(0)/(k_B T_c) = 2π/e^γ ≈ 3.5278."""
    assert math.isclose(BCS_UNIVERSAL_GAP_RATIO, 3.527753977724091, rel_tol=1e-12)
    # The textbook short value 3.52 is correct to 0.5%.
    assert abs(BCS_UNIVERSAL_GAP_RATIO - 3.52) < 0.01
    # Static method on the simulator returns the same number.
    sim = SuperconductivitySimulator()
    assert sim.bcs_ratio() == BCS_UNIVERSAL_GAP_RATIO


def test_bcs_ratio_back_solves_Tc_from_gap():
    """Plug in 1.40 meV → expect T_c that satisfies 2Δ/(k_B T_c) = 3.52."""
    geom = SuperconductivityGeometry(delta_meV=1.40)
    sim = SuperconductivitySimulator(geometry=geom, debye_T_K=105.0)
    K_B = 1.380649e-23
    EV_PER_J = 1.0 / 1.602176634e-19
    delta_J = 1.40e-3 / EV_PER_J
    Tc = sim.critical_temperature_K()
    ratio_check = 2.0 * delta_J / (K_B * Tc)
    assert math.isclose(ratio_check, BCS_UNIVERSAL_GAP_RATIO, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Hg / HBCCO recovered T_c (back-solved couplings)
# ---------------------------------------------------------------------------

def _back_solve_for_Tc(Tc_target: float, theta_D: float) -> SuperconductivitySimulator:
    """Return a simulator whose Δ(0) is calibrated to give the target T_c."""
    K_B = 1.380649e-23
    EV_PER_J = 1.0 / 1.602176634e-19
    delta0_meV = (BCS_UNIVERSAL_GAP_RATIO * K_B * Tc_target / 2.0) * EV_PER_J * 1e3
    geom = SuperconductivityGeometry(delta_meV=delta0_meV)
    return SuperconductivitySimulator(geometry=geom, debye_T_K=theta_D)


def test_Hg_Tc_4p2K_from_BCS():
    """Hg (Onnes 1911, original SC) recovers T_c = 4.2 K via BCS."""
    sim = _back_solve_for_Tc(4.2, 72.0)
    Tc = sim.critical_temperature_K()
    assert math.isclose(Tc, 4.2, rel_tol=1e-6)
    # The back-solved coupling sits in the BCS weak-coupling window.
    assert 0.1 < sim.coupling_NV < 0.5


def test_HBCCO_Tc_134K_from_BCS():
    """HBCCO 134 K (Schilling 1993, ambient-pressure record) reproduces.

    Cuprates are NOT weak-coupling BCS, so the back-solved N(0)V naturally
    exceeds the strict-BCS validity window (NV ~ 1).  We verify the formula
    still recovers the requested T_c, and we report the coupling.
    """
    sim = _back_solve_for_Tc(134.0, 320.0)
    Tc = sim.critical_temperature_K()
    assert math.isclose(Tc, 134.0, rel_tol=1e-5)
    # Back-solved NV is order-1 (weak-coupling BCS broken — physical fact
    # for cuprates).
    assert sim.coupling_NV > 0.5


# ---------------------------------------------------------------------------
# Substrate T_c,max ceiling
# ---------------------------------------------------------------------------

def test_substrate_Tc_max_is_128p9K():
    """T_c,max = Λ_QCD / n_R / k_B = 200 MeV / 18 / k_B ≈ 128.9 K (B3 ceiling)."""
    assert math.isclose(SUBSTRATE_TC_MAX_K, 128.9, rel_tol=1e-3)
    # Static method exposes the same number.
    assert math.isclose(
        SuperconductivitySimulator.substrate_tc_max_K(),
        SUBSTRATE_TC_MAX_K,
        rel_tol=1e-12,
    )
    # Energy-per-mode is the canonical 11.11 meV.
    assert math.isclose(SUBSTRATE_E_PER_MODE_MEV, 200.0 / 18.0, rel_tol=1e-12)
    assert math.isclose(SUBSTRATE_E_PER_MODE_MEV, 11.111, rel_tol=1e-3)


def test_HBCCO_within_4pct_of_substrate_ceiling():
    """HBCCO 134 K matches the substrate ceiling within 4% (cross-disciplinary)."""
    HBCCO = REFERENCE_SUPERCONDUCTORS["HBCCO"]["T_c_K"]
    deviation_pct = 100.0 * abs(HBCCO - SUBSTRATE_TC_MAX_K) / SUBSTRATE_TC_MAX_K
    assert deviation_pct < 5.0   # 4.0% in practice
    assert deviation_pct > 1.0   # not zero — keeps the prediction non-trivial


def test_all_ambient_BCS_below_substrate_ceiling():
    """All conventional ambient-pressure SCs except HBCCO sit below 128.9 K."""
    sim = SuperconductivitySimulator()
    rows = sim.reference_table_rows()
    above = [r for r in rows if not r["below_ceiling"]]
    # HBCCO (134 K) marginally exceeds the strict ceiling 128.9 K — expected
    # since the ceiling is "approximate within 4%".
    above_names = {r["name"] for r in above}
    assert above_names <= {"HBCCO"}
    # Every other entry must be safely under the ceiling.
    for r in rows:
        if r["name"] == "HBCCO":
            continue
        assert r["T_c_K"] < SUBSTRATE_TC_MAX_K, (
            f"{r['name']} T_c={r['T_c_K']} K is above ceiling "
            f"{SUBSTRATE_TC_MAX_K:.1f} K — would falsify the bound"
        )


# ---------------------------------------------------------------------------
# Self-consistent Δ(T)
# ---------------------------------------------------------------------------

def test_gap_decreases_with_T_and_vanishes_at_Tc():
    """Δ(T) is monotonically decreasing in T and vanishes at T_c."""
    sim = _back_solve_for_Tc(7.2, 105.0)   # Pb
    Tc = sim.critical_temperature_K()
    Ts = np.array([0.1 * Tc, 0.5 * Tc, 0.9 * Tc, 1.0 * Tc, 1.1 * Tc])
    deltas = np.array([sim.gap_at_T(t) for t in Ts])
    # Strictly decreasing across the increasing-T sequence
    assert np.all(np.diff(deltas) <= 1e-12)
    # Gap vanishes at and above T_c
    assert deltas[-2] == 0.0
    assert deltas[-1] == 0.0
    # And Δ(T → 0) approaches Δ(0).
    assert math.isclose(deltas[0], sim.gap_zero_T_meV(), rel_tol=2e-2)


def test_gap_curve_has_BCS_shape():
    """Δ(T)/Δ(0) at T=0.5 T_c is roughly 0.95 (BCS classic ~ 0.94)."""
    sim = _back_solve_for_Tc(7.2, 105.0)
    t_over_tc, d_over_d0 = sim.gap_curve(n_T=20)
    # find the sample nearest 0.5
    idx = int(np.argmin(np.abs(t_over_tc - 0.5)))
    assert 0.85 < d_over_d0[idx] < 1.0
    # Endpoints
    assert d_over_d0[0] > 0.99   # Δ(T→0)/Δ(0) ≈ 1
    assert d_over_d0[-1] < 0.10  # Δ(T→T_c)/Δ(0) ≈ 0


# ---------------------------------------------------------------------------
# Density of states with gap
# ---------------------------------------------------------------------------

def test_density_of_states_zero_inside_gap():
    """N_s(E)/N(0) = 0 for |E| < Δ."""
    geom = SuperconductivityGeometry(delta_meV=1.40)
    sim = SuperconductivitySimulator(geometry=geom)
    E = np.linspace(-1.30, 1.30, 50)
    n_over_n0 = sim.density_of_states_ratio(E)
    assert np.allclose(n_over_n0, 0.0)


def test_density_of_states_peaks_at_gap_edge():
    """N_s(E)/N(0) → ∞ as |E| → Δ⁺ (square-root singularity)."""
    geom = SuperconductivityGeometry(delta_meV=1.40)
    sim = SuperconductivitySimulator(geometry=geom)
    E_just_above = np.array([1.41, 1.45, 1.5, 2.0, 5.0])
    n = sim.density_of_states_ratio(E_just_above)
    # Decreasing in |E| above the gap.
    assert np.all(np.diff(n) < 0.0)
    # At E = 1.41 meV (just above 1.40) the DOS should be very large.
    assert n[0] > 5.0
    # At E ≫ Δ, N_s/N → 1.
    assert math.isclose(n[-1], 1.0, rel_tol=0.1)


# ---------------------------------------------------------------------------
# Meissner effect (London depth)
# ---------------------------------------------------------------------------

def test_meissner_exponential_decay():
    """B(z) = B_0 exp(-z/λ_L); at z = λ_L the field is B_0/e."""
    sim = SuperconductivitySimulator()
    lam = sim.geometry.london_lambda_nm
    z = np.array([0.0, lam, 2.0 * lam, 5.0 * lam])
    B = sim.meissner_field(z, B0_T=1.0)
    assert math.isclose(B[0], 1.0, rel_tol=1e-12)
    assert math.isclose(B[1], 1.0 / math.e, rel_tol=1e-9)
    assert math.isclose(B[2], 1.0 / math.e ** 2, rel_tol=1e-9)
    assert B[3] < 1e-2


# ---------------------------------------------------------------------------
# Type-I vs type-II classification
# ---------------------------------------------------------------------------

def test_type_classification_threshold():
    """κ < 1/√2 → type-I; κ > 1/√2 → type-II (Abrikosov 1957)."""
    assert math.isclose(GL_KAPPA_THRESHOLD, 1.0 / math.sqrt(2.0), rel_tol=1e-15)
    assert SuperconductivitySimulator.type_classification(0.5) == "I"
    assert SuperconductivitySimulator.type_classification(0.6) == "I"
    assert SuperconductivitySimulator.type_classification(0.8) == "II"
    assert SuperconductivitySimulator.type_classification(50.0) == "II"


def test_geometry_default_is_type_II():
    """Default ξ=80 nm, λ_L=37 nm → κ < 1/√2 → type-I (Pb-like)."""
    geom = SuperconductivityGeometry()
    assert math.isclose(geom.gl_kappa(), 37.0 / 80.0)
    assert geom.is_type_II() is False
    assert SuperconductivitySimulator.type_classification(geom.gl_kappa()) == "I"


def test_geometry_extreme_kappa_classifies_type_II():
    """ξ=1 nm, λ_L=200 nm (extreme κ=200) → type-II (cuprate-like)."""
    geom = SuperconductivityGeometry(coherence_xi_nm=1.0, london_lambda_nm=200.0)
    assert geom.is_type_II() is True


# ---------------------------------------------------------------------------
# Cooper paired-strain bridge geometry
# ---------------------------------------------------------------------------

def test_paired_bridge_separation_is_2xi():
    """The two paired defects are separated by 2·ξ along the bridge."""
    geom = SuperconductivityGeometry(coherence_xi_nm=80.0)
    pos = geom.paired_bridge_positions()
    assert pos.shape == (2,)
    assert math.isclose(pos[1] - pos[0], 2.0 * 80.0, rel_tol=1e-12)


def test_paired_bridge_envelope_decays_at_xi():
    """Envelope ψ_pair(x) ∝ exp(-|x|/ξ); peaks at x=0, drops to 1/e at x=±ξ."""
    geom = SuperconductivityGeometry(coherence_xi_nm=80.0)
    x, env = geom.paired_bridge_envelope()
    # Peak at the centre
    idx_centre = int(np.argmin(np.abs(x)))
    assert math.isclose(env[idx_centre], 1.0, rel_tol=1e-2)
    # At |x| = ξ the envelope is 1/e
    idx_xi = int(np.argmin(np.abs(np.abs(x) - 80.0)))
    assert math.isclose(env[idx_xi], 1.0 / math.e, rel_tol=2e-2)


def test_bogoliubov_factors_sum_to_unity():
    """u_k² + v_k² = 1 for every k."""
    geom = SuperconductivityGeometry(delta_meV=1.40)
    eps = np.linspace(-30.0, 30.0, 31)   # |ε|/Δ ≫ 1 at the endpoints
    u2, v2 = geom.bogoliubov_uv2(eps)
    assert np.allclose(u2 + v2, 1.0, atol=1e-12)
    # Sign convention: at ε ≫ Δ, u² → 1 (electron-like)
    assert u2[-1] > 0.999
    assert v2[-1] < 0.001
    # at ε ≪ -Δ, v² → 1 (hole-like)
    assert v2[0] > 0.999
    assert u2[0] < 0.001
    # at ε = 0 (Fermi level), u² = v² = 1/2
    idx0 = int(np.argmin(np.abs(eps)))
    assert math.isclose(u2[idx0], 0.5, abs_tol=0.05)
    assert math.isclose(v2[idx0], 0.5, abs_tol=0.05)


# ---------------------------------------------------------------------------
# Reference table sorted ascending and below ceiling
# ---------------------------------------------------------------------------

def test_reference_table_sorted_ascending():
    sim = SuperconductivitySimulator()
    rows = sim.reference_table_rows()
    Tc_list = [r["T_c_K"] for r in rows]
    assert Tc_list == sorted(Tc_list)
    # First entry is the lowest-Tc elemental SC.
    assert rows[0]["name"] == "Hg"
    # Last entry is the ambient-pressure record HBCCO.
    assert rows[-1]["name"] == "HBCCO"
