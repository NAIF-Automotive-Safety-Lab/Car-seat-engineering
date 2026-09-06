"""Unit tests for the substrate field solver.

Verifies the substrate Lagrangian (§18.11) actually produces the kink
geometry by:
  - confirming the analytic ansatz is a stationary point;
  - confirming gradient flow from a crude initial guess produces the kink;
  - confirming the kink rest energy converges to the analytic value 8;
  - confirming kink-antikink attraction with the predicted exponential.
"""

import numpy as np
import pytest

from stiff_medium.substrate_field_solver import (
    kink_profile_1d,
    antikink_profile_1d,
    energy_density_1d,
    total_energy_1d,
    static_residual_1d,
    validate_static_kink_1d,
    relax_kink_from_step_initial,
    cylindrical_kink_energy_per_length,
    axisymmetric_lump_energy,
    two_kink_interaction_unrelaxed_curve,
    two_kink_same_charge_unrelaxed,
    geometric_origin_of_8,
    kink_energy_with_exact_derivative,
)


# ---------------------------------------------------------------------------
# Static 1D kink — analytic ansatz
# ---------------------------------------------------------------------------

def test_analytic_kink_boundary_conditions():
    """φ(-∞) = 0 and φ(+∞) = 2π for the analytic kink."""
    x = np.array([-100.0, +100.0])
    phi = kink_profile_1d(x)
    assert phi[0] == pytest.approx(0.0, abs=1e-30)
    assert phi[1] == pytest.approx(2.0 * np.pi, abs=1e-30)


def test_analytic_kink_centre_value():
    """φ(0) = 4 arctan(1) = π."""
    phi0 = kink_profile_1d(np.array([0.0]))[0]
    assert phi0 == pytest.approx(np.pi, abs=1e-12)


def test_static_kink_energy_converges_to_8():
    """E_K should converge to 8 with refinement (via exact derivative)."""
    e = kink_energy_with_exact_derivative(L=40.0, N=801)
    assert e == pytest.approx(8.0, abs=1e-10)


def test_validate_static_kink_residual_small():
    """The analytic kink is a stationary point: residual is small."""
    res = validate_static_kink_1d(L=30.0, N=1601)
    assert res.energy_numeric == pytest.approx(8.0, rel=1e-3)
    assert res.max_residual < 1e-3


# ---------------------------------------------------------------------------
# Gradient flow — substrate produces the kink
# ---------------------------------------------------------------------------

def test_gradient_flow_from_step_recovers_kink_energy():
    """Starting from a step ansatz, gradient flow drives E → 8."""
    res = relax_kink_from_step_initial(L=30.0, N=401)
    assert res.energy_numeric == pytest.approx(8.0, rel=2e-3)
    # EOM residual should approach zero — proof of stationarity
    assert res.max_residual < 1e-6


# ---------------------------------------------------------------------------
# 3D cylindrical kink — kink line
# ---------------------------------------------------------------------------

def test_cylindrical_kink_matches_1d():
    """3D cylindrical kink line has same energy/length as 1D kink."""
    res = cylindrical_kink_energy_per_length(N_z=1601)
    assert res.asymptotic_to_1d == pytest.approx(1.0, abs=2e-3)


# ---------------------------------------------------------------------------
# Derrick's theorem
# ---------------------------------------------------------------------------

def test_derrick_3d_lump_collapses():
    """A 3D spherical lump should collapse to vacuum (E → 0)."""
    res = axisymmetric_lump_energy(L=6.0, N=121, n_iter=8000)
    assert res["energy"] < 0.5
    assert abs(res["phi_at_origin"]) < 0.5


# ---------------------------------------------------------------------------
# Two-kink interaction
# ---------------------------------------------------------------------------

def test_kink_antikink_is_attractive():
    """Kink-antikink potential is negative at large R (attraction)."""
    R = np.array([3.0, 4.0, 5.0, 6.0, 7.0])
    res = two_kink_interaction_unrelaxed_curve(R, L=30.0, N=1601)
    assert all(res["V_int_numeric"] < 0)


def test_kink_kink_is_repulsive():
    """Same-charge kink-kink potential is positive at large R."""
    R = np.array([3.0, 4.0, 5.0, 6.0, 7.0])
    res = two_kink_same_charge_unrelaxed(R, L=30.0, N=1601)
    assert all(res["V_int_numeric"] > 0)


def test_kink_antikink_exponential_falloff():
    """Asymptotic |V_int(R)| ∝ e^{-R/ξ}: log-fit slope ≈ -1."""
    R = np.linspace(4.0, 7.0, 7)
    res = two_kink_interaction_unrelaxed_curve(R, L=30.0, N=1601)
    V = res["V_int_numeric"]
    slope, _ = np.polyfit(R, np.log(-V), 1)
    assert slope == pytest.approx(-1.0, abs=0.06)


# ---------------------------------------------------------------------------
# Geometric origin of '8'
# ---------------------------------------------------------------------------

def test_eight_is_exact():
    """E_K = 8 exactly: numerical convergence at machine precision."""
    res = geometric_origin_of_8((201, 401, 801))
    # With exact derivative, Simpson reaches ~1e-13 by N=801
    assert res["rel_errors_with_exact_derivative"][-1] < 1e-12


def test_centred_difference_converges_quadratically():
    """Centred-difference E_K converges as O(N⁻²)."""
    res = geometric_origin_of_8((201, 401, 801, 1601))
    rel = res["rel_errors"]
    # Each refinement should reduce the error by roughly 4×
    for i in range(len(rel) - 1):
        ratio = rel[i] / rel[i + 1]
        assert ratio > 3.0, f"Convergence ratio {ratio:.2f} at N={res['N_values'][i]}"


# ---------------------------------------------------------------------------
# Energy density basics
# ---------------------------------------------------------------------------

def test_vacuum_energy_density_zero():
    """Energy density vanishes at φ = 0 (vacuum)."""
    phi = np.zeros(11)
    e = energy_density_1d(phi, dx=0.1)
    assert np.allclose(e, 0.0)


def test_vacuum_energy_density_zero_at_2pi():
    """Energy density also vanishes at φ = 2π (the other vacuum)."""
    phi = np.full(11, 2.0 * np.pi)
    e = energy_density_1d(phi, dx=0.1)
    assert np.allclose(e, 0.0, atol=1e-14)


def test_residual_zero_for_constant_field():
    """∇²(const) − sin(0) = 0, ∇²(const) − sin(2π) = 0."""
    phi = np.zeros(50)
    r = static_residual_1d(phi, dx=0.1)
    assert np.allclose(r, 0.0)
    phi = np.full(50, 2.0 * np.pi)
    r = static_residual_1d(phi, dx=0.1)
    assert np.allclose(r, 0.0, atol=1e-14)


def test_antikink_complements_kink():
    """kink(x) + antikink(x) = 2π (antikink is the topological complement)."""
    x = np.linspace(-5, 5, 21)
    s = kink_profile_1d(x) + antikink_profile_1d(x)
    assert np.allclose(s, 2.0 * np.pi)
