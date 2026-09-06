"""Tests for catalysis_substrate.

Covers:
  - CatalysisGeometry: Sabatier optimum near E_bind ~ 0; BEP barrier symmetry.
  - Sabatier volcano: Pt is the optimal metal for H2 dissociation
    (peak TOF among Au/Cu/Pt/Pd/Ni/Ru/Re/W).
  - Michaelis-Menten saturation: v -> k_cat [E]_0 as [S] -> infinity;
    v -> k_cat [E]_0 / 2 at [S] = K_M.
  - Lineweaver-Burk linearization: 1/v vs 1/[S] is straight.
  - Carbonic anhydrase k_cat / K_M ~ 1.25e8 within order of magnitude
    of reference 1.5e8.
  - Hill coefficient fit recovers Hb n_H = 2.8.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from stiff_medium.catalysis_substrate import (
    CatalysisGeometry,
    CatalysisSimulator,
    H_BINDING_eV,
    ENZYME_KCAT_KM,
    HILL_REFERENCE,
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def sim():
    s = CatalysisSimulator()
    s.compare_to_reference()
    return s


def _close(a, b, pct):
    return abs(a - b) / abs(b) * 100.0 < pct


# ---------------------------------------------------------------------------
# 1. CatalysisGeometry: Sabatier and BEP
# ---------------------------------------------------------------------------

def test_geometry_default_pt_is_optimal():
    g = CatalysisGeometry()           # Pt-like
    assert g.is_sabatier_optimal(tol=0.15)


def test_geometry_au_too_weak_not_optimal():
    g = CatalysisGeometry(binding_energy=H_BINDING_eV["Au"])
    assert not g.is_sabatier_optimal(tol=0.15)


def test_geometry_w_too_strong_not_optimal():
    g = CatalysisGeometry(binding_energy=H_BINDING_eV["W"])
    assert not g.is_sabatier_optimal(tol=0.15)


def test_bep_barrier_symmetry():
    """E_a is minimized at E_bind = 0; rises on both sides."""
    g = CatalysisGeometry(bep_alpha=0.5, bep_intercept=0.75)
    e_a_zero = g.activation_energy(0.0)
    e_a_pos = g.activation_energy(+0.5)
    e_a_neg = g.activation_energy(-0.5)
    assert e_a_pos > e_a_zero, (e_a_pos, e_a_zero)
    assert e_a_neg > e_a_zero, (e_a_neg, e_a_zero)


def test_transition_state_position_above_site():
    g = CatalysisGeometry()
    ts = g.transition_state_position()
    assert ts[2] > g.site_position[2]


def test_strain_potential_minimum_at_re():
    g = CatalysisGeometry(binding_energy=-0.10, site_width=1.5)
    # at r = r_e + site, V = -D_e (the well depth)
    r_well = g.site_position + np.array([0.0, 0.0, g.site_width])
    V_well = float(g.strain_potential(r_well))
    # at r far away, V -> 0
    r_far = g.site_position + np.array([0.0, 0.0, 20.0])
    V_far = float(g.strain_potential(r_far))
    assert V_well < V_far
    assert V_well == pytest.approx(-abs(g.binding_energy), abs=1e-9)


# ---------------------------------------------------------------------------
# 2. Sabatier volcano
# ---------------------------------------------------------------------------

def test_sabatier_pt_is_optimal_for_hydrogenation(sim):
    """Pt sits at the volcano peak for H2 dissociation."""
    table = sim.results["sabatier_table"]
    peak = sim.results["sabatier_peak_metal"]
    # Peak must be near the Sabatier optimum: |E_bind| < 0.15 eV.
    assert abs(table[peak]["E_bind_eV"]) < 0.15
    # Pt itself must be among the top three metals.
    sorted_metals = sorted(
        table.keys(), key=lambda m: table[m]["TOF_s_inv"], reverse=True
    )
    assert "Pt" in sorted_metals[:3], sorted_metals[:3]


def test_sabatier_au_and_w_both_slow():
    """Au binds too weakly; W binds too strongly. Both are slow."""
    s = CatalysisSimulator()
    v = s.sabatier_volcano(metals=["Au", "Pt", "W"])
    assert v["Pt"]["TOF_s_inv"] > v["Au"]["TOF_s_inv"]
    assert v["Pt"]["TOF_s_inv"] > v["W"]["TOF_s_inv"]


def test_sabatier_volcano_shape_two_sided():
    """The TOF as a function of E_bind is concave (volcano-shaped)."""
    s = CatalysisSimulator()
    e_binds = np.linspace(-0.8, +0.4, 25)
    rates = []
    for eb in e_binds:
        g = CatalysisGeometry(binding_energy=float(eb))
        rates.append(s.eyring_rate(g.activation_energy(float(eb))))
    rates = np.array(rates)
    peak_idx = int(np.argmax(rates))
    # Peak must lie strictly inside the scan, not at either edge
    assert 0 < peak_idx < len(rates) - 1
    # Rates must rise then fall around the peak
    assert rates[peak_idx] > rates[0]
    assert rates[peak_idx] > rates[-1]


# ---------------------------------------------------------------------------
# 3. Michaelis-Menten saturation
# ---------------------------------------------------------------------------

def test_mm_saturation_at_high_S():
    """As [S] >> K_M, v -> V_max = k_cat [E]_0."""
    k_cat, K_M, E0 = 100.0, 1e-3, 1e-6
    V_max = k_cat * E0
    v_sat = CatalysisSimulator.michaelis_menten(
        S=1.0, k_cat=k_cat, K_M=K_M, E_total=E0
    )
    assert _close(v_sat, V_max, 1.0)


def test_mm_half_max_at_KM():
    """At [S] = K_M, v = V_max / 2."""
    k_cat, K_M, E0 = 50.0, 2e-4, 1e-6
    V_max = k_cat * E0
    v = CatalysisSimulator.michaelis_menten(
        S=K_M, k_cat=k_cat, K_M=K_M, E_total=E0
    )
    assert _close(v, V_max / 2.0, 1e-6)


def test_mm_linear_at_low_S():
    """At [S] << K_M, v ~ (k_cat / K_M) * [E]_0 * [S] (first-order in S)."""
    k_cat, K_M, E0 = 80.0, 5e-3, 1e-6
    S = 1e-5     # well below K_M
    v_calc = CatalysisSimulator.michaelis_menten(
        S=S, k_cat=k_cat, K_M=K_M, E_total=E0
    )
    v_pred = k_cat * E0 * S / K_M
    assert _close(v_calc, v_pred, 1.0)


def test_lineweaver_burk_linear():
    """1/v vs 1/[S] is exactly linear with intercept 1/V_max and slope K_M/V_max."""
    k_cat, K_M, E0 = 100.0, 1e-3, 1e-6
    V_max = k_cat * E0
    S = np.logspace(-5, -1, 30)
    inv_S, inv_v = CatalysisSimulator.lineweaver_burk(S, k_cat, K_M, E0)
    slope, intercept = np.polyfit(inv_S, inv_v, 1)
    assert _close(intercept, 1.0 / V_max, 1.0)
    assert _close(slope, K_M / V_max, 1.0)


# ---------------------------------------------------------------------------
# 4. Enzyme efficiency
# ---------------------------------------------------------------------------

def test_carbonic_anhydrase_perfect_enzyme(sim):
    """k_cat / K_M for carbonic anhydrase is order-of-magnitude 10^8."""
    eff = sim.results["CA_k_cat_over_K_M"]
    assert 1e7 < eff < 1e9


def test_carbonic_anhydrase_k_cat_million(sim):
    """k_cat ~ 10^6 s^-1 for carbonic anhydrase."""
    k_cat = sim.results["CA_k_cat_s_inv"]
    assert 1e5 < k_cat < 1e7


def test_enzyme_efficiency_ratio():
    s = CatalysisSimulator()
    eff = s.enzyme_efficiency(k_cat=1e6, K_M=1e-2)
    assert eff == pytest.approx(1e8, rel=1e-9)


# ---------------------------------------------------------------------------
# 5. Hill cooperativity for hemoglobin
# ---------------------------------------------------------------------------

def test_hill_coefficient_hemoglobin_28(sim):
    """Hill plot fit recovers n_H = 2.8 for hemoglobin."""
    n_H_fit = sim.results["Hb_n_H_fit"]
    assert _close(n_H_fit, 2.8, 1.0)


def test_hill_K_half_hemoglobin(sim):
    """K_half (P_50) for Hb is ~ 26 mmHg."""
    k_half = sim.results["Hb_K_half_fit_mmHg"]
    assert _close(k_half, 26.0, 1.0)


def test_hill_n1_reduces_to_michaelis_menten():
    """n_H = 1 reduces to the Michaelis-Menten form."""
    S = np.array([0.5, 1.0, 2.0, 5.0])
    v_hill = CatalysisSimulator.hill_equation(S, V_max=1.0, K_half=1.0, n_H=1.0)
    v_mm = CatalysisSimulator.michaelis_menten(S, k_cat=1.0, K_M=1.0, E_total=1.0)
    np.testing.assert_allclose(v_hill, v_mm, rtol=1e-12)


def test_hill_cooperativity_steeper_than_mm():
    """At n_H > 1, the curve is steeper than MM at the inflection point."""
    S = np.linspace(0.1, 10.0, 100)
    v_hill = CatalysisSimulator.hill_equation(S, V_max=1.0, K_half=1.0, n_H=2.8)
    v_mm = CatalysisSimulator.hill_equation(S, V_max=1.0, K_half=1.0, n_H=1.0)
    # find slope at K_half (S=1)
    idx = int(np.argmin(np.abs(S - 1.0)))
    slope_hill = (v_hill[idx + 1] - v_hill[idx - 1]) / (S[idx + 1] - S[idx - 1])
    slope_mm = (v_mm[idx + 1] - v_mm[idx - 1]) / (S[idx + 1] - S[idx - 1])
    assert slope_hill > slope_mm


# ---------------------------------------------------------------------------
# 6. Eyring TST sanity check
# ---------------------------------------------------------------------------

def test_eyring_low_barrier_fast():
    s = CatalysisSimulator()
    rate_low = s.eyring_rate(0.1)        # easy barrier
    rate_high = s.eyring_rate(1.5)       # hard barrier
    assert rate_low > rate_high


def test_eyring_temperature_dependence():
    s_cold = CatalysisSimulator(temperature=200.0)
    s_hot = CatalysisSimulator(temperature=500.0)
    e_a = 0.5
    assert s_hot.eyring_rate(e_a) > s_cold.eyring_rate(e_a)


# ---------------------------------------------------------------------------
# 7. End-to-end report runs cleanly
# ---------------------------------------------------------------------------

def test_report_runs():
    s = CatalysisSimulator()
    txt = s.report()
    assert "Sabatier" in txt
    assert "Carbonic anhydrase" in txt
    assert "Hemoglobin" in txt


def test_compare_table_keys():
    s = CatalysisSimulator()
    table = s.compare_to_reference()
    for k in ("Sabatier_peak_metal", "CA_k_cat_over_K_M", "Hb_Hill_n"):
        assert k in table
