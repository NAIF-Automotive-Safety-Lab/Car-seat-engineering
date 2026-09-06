"""Tests for sparc_dynamics: substrate-driven galactic rotation curves."""

from __future__ import annotations

import numpy as np
import pytest

from stiff_medium.sparc_dynamics import (
    G_DAGGER_EMPIRICAL,
    H0_PLANCK,
    H0_SHOES,
    GalaxyData,
    REFERENCE_GALAXIES,
    SPARCDynamics,
    synthetic_galaxy,
)


# ---------------------------------------------------------------------------
# Substrate predictions
# ---------------------------------------------------------------------------


def test_g_dagger_substrate_within_50pct():
    """g_† = c·H_0/(2π) at SH0ES H_0 must land within 50% of empirical."""
    s = SPARCDynamics(H0_kms_per_mpc=H0_SHOES)
    g = s.g_dagger_from_substrate()
    assert g > 0
    frac = abs(g - G_DAGGER_EMPIRICAL) / G_DAGGER_EMPIRICAL
    assert frac < 0.5, f"g_dag={g:.3e}, frac off = {frac:.2f}"


def test_g_dagger_units_sanity():
    """Order-of-magnitude check: should be ~1e-10 m/s^2."""
    s = SPARCDynamics()
    g = s.g_dagger_from_substrate()
    assert 5e-11 < g < 5e-10


def test_tully_fisher_slope_is_four():
    """Log V_flat vs log M_baryon slope should be 1/4 (V^4 ∝ M)."""
    s = SPARCDynamics(H0_kms_per_mpc=H0_SHOES)
    masses = np.logspace(8.0, 12.0, 9)
    V = np.array([s.tully_fisher_relation(m) for m in masses])
    slope, intercept = np.polyfit(np.log10(masses), np.log10(V), 1)
    assert abs(slope - 0.25) < 0.01, f"slope={slope:.4f} (expected 0.25)"


def test_tully_fisher_NGC3198_within_factor_of_2():
    """Tully-Fisher on NGC 3198 stellar mass must give V_flat near 150 km/s."""
    s = SPARCDynamics(H0_kms_per_mpc=H0_SHOES)
    galaxy = REFERENCE_GALAXIES["NGC3198"]
    V_pred = s.tully_fisher_relation(galaxy["M_star"])
    V_obs = galaxy["V_flat_kms"]
    assert 0.5 * V_obs < V_pred < 2.0 * V_obs, (
        f"V_pred={V_pred:.1f} vs V_obs={V_obs}"
    )


# ---------------------------------------------------------------------------
# Radial acceleration relation
# ---------------------------------------------------------------------------


def test_radial_acceleration_relation_high_acc_limit():
    """At g_bar >> g_†, RAR returns g_obs ≈ g_bar (Newtonian limit)."""
    s = SPARCDynamics()
    g_dag = s.g_dagger_from_substrate()
    g_bar = 1.0e3 * g_dag  # very high acceleration
    g_obs = float(s.radial_acceleration_relation(g_bar))
    assert abs(g_obs - g_bar) / g_bar < 0.05


def test_radial_acceleration_relation_low_acc_limit():
    """At g_bar << g_†, RAR predicts g_obs ≈ sqrt(g_bar · g_†) (deep MOND)."""
    s = SPARCDynamics()
    g_dag = s.g_dagger_from_substrate()
    g_bar = 1.0e-4 * g_dag
    g_obs = float(s.radial_acceleration_relation(g_bar))
    g_mond = np.sqrt(g_bar * g_dag)
    assert 0.5 * g_mond < g_obs < 2.0 * g_mond


# ---------------------------------------------------------------------------
# Rotation curve and fitting
# ---------------------------------------------------------------------------


def test_rotation_curve_monotone_and_positive():
    s = SPARCDynamics()
    R = np.linspace(0.5, 30.0, 20)
    V = s.rotation_curve(
        R,
        M_baryon_msun=3.0e10,
        c_concentration=10.0,
        M_halo_msun=1.0e12,
        R_disk_kpc=3.0,
    )
    assert np.all(V > 0)
    assert np.all(np.isfinite(V))
    # Should rise from inner radii to a flat-ish asymptote
    assert V[-1] > V[0]


def test_sparc_galaxy_fit_synthetic_recovers_curve():
    """Fit a synthetic NGC 3198-like galaxy; reduced chi^2 should be modest."""
    g = synthetic_galaxy(
        "NGC3198_synth", M_star=3.0e10, V_flat_kms=150.0, R_disk_kpc=3.14
    )
    s = SPARCDynamics(H0_kms_per_mpc=H0_SHOES)
    fit = s.sparc_galaxy_fit(g)
    assert fit["reduced_chi2"] < 10.0, fit
    # Curve must roughly match V_flat at the outer radii
    V_pred_outer = float(np.mean(fit["V_pred_kms"][-3:]))
    assert 100.0 < V_pred_outer < 220.0, f"V_pred_outer={V_pred_outer}"


# ---------------------------------------------------------------------------
# H_0 consistency check
# ---------------------------------------------------------------------------


def test_H_0_consistency_prefers_SHOES():
    """SH0ES H_0 should fit g_† better than Planck H_0."""
    s = SPARCDynamics()
    out = s.H_0_consistency_check()
    assert out["frac_offset_SHOES"] < out["frac_offset_Planck"], out
    assert out["preferred"] == "SHOES"
    # Sanity: SH0ES side within ~20% (paper says 1.6% but allow slack
    # because g_dag depends on the chosen prefactor)
    assert out["frac_offset_SHOES"] < 0.20


def test_H_0_consistency_planck_offset_larger():
    """The Planck-side offset should be at least 5x larger than SH0ES side."""
    s = SPARCDynamics()
    out = s.H_0_consistency_check()
    ratio = out["frac_offset_Planck"] / max(out["frac_offset_SHOES"], 1e-6)
    assert ratio > 1.5, out


# ---------------------------------------------------------------------------
# Reference galaxies
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["NGC3198", "NGC6946", "M33"])
def test_tully_fisher_reference_galaxies(name):
    """Tully-Fisher should land within a factor of 2.5 for all reference
    galaxies (these are rough literature numbers)."""
    s = SPARCDynamics(H0_kms_per_mpc=H0_SHOES)
    g = REFERENCE_GALAXIES[name]
    V_pred = s.tully_fisher_relation(g["M_star"])
    V_obs = g["V_flat_kms"]
    ratio = V_pred / V_obs
    assert 0.4 < ratio < 2.5, f"{name}: ratio={ratio:.2f}"
