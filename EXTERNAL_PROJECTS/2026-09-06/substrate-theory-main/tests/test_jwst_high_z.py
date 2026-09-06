"""Tests for src.stiff_medium.jwst_high_z."""

import math

import pytest

from src.stiff_medium.jwst_high_z import JWSTHighZ


@pytest.fixture
def model():
    return JWSTHighZ()


def test_halo_mass_function_decreases_with_M(model):
    """At fixed z, dn/dlogM must decrease with halo mass M."""
    for z in (10.0, 12.0, 15.0):
        phi_prev = math.inf
        for logM in (9.0, 10.0, 11.0, 12.0, 13.0):
            phi = model.halo_mass_function(z, 10 ** logM)
            assert phi > 0
            assert phi < phi_prev, (
                f"HMF should decrease with M at z={z}; "
                f"got phi(logM={logM})={phi:.3e} >= prev {phi_prev:.3e}"
            )
            phi_prev = phi


def test_halo_mass_function_rejects_nonpositive(model):
    with pytest.raises(ValueError):
        model.halo_mass_function(z=12.0, M=0.0)


def test_high_z_galaxies_plausible(model):
    """At z>10, JADES-mass galaxies should have non-negligible predicted density."""
    rows = model.compare_to_data()
    high_z_rows = [r for r in rows if r["z"] > 10.0]
    assert len(high_z_rows) >= 3
    for r in high_z_rows:
        phi = r["phi_pred_Mpc-3_dex-1"]
        # Plausible high-z SMF densities lie in ~ 1e-7 .. 1e-2 Mpc^-3 dex^-1
        assert 1e-9 < phi < 1e-1, (
            f"Implausible phi for {r['name']} at z={r['z']}: {phi:.3e}"
        )


def test_stellar_mass_function_consistent(model):
    """SMF at z=10, M*=1e8.5 should be positive and below low-z normalisation."""
    phi_z10 = model.stellar_mass_function(10.0, 10 ** 8.5)
    phi_z6 = model.stellar_mass_function(6.0, 10 ** 8.5)
    assert phi_z10 > 0
    assert phi_z6 > 0
    assert phi_z6 > phi_z10  # density grows with cosmic time


def test_reionization_tau_in_planck_range(model):
    """tau should land in the broad observational window [0.05, 0.08]."""
    reion = model.reionization_history()
    tau = reion["tau"]
    assert 0.05 <= tau <= 0.08, f"tau out of range: {tau:.4f}"
    assert reion["z_50pct_reionised"] > 6.0
    assert reion["z_complete"] < reion["z_50pct_reionised"]


def test_jades_and_ceers_samples_have_high_z(model):
    jades = model.JADES_observations()
    ceers = model.CEERS_observations()
    assert len(jades) >= 4
    assert len(ceers) >= 3
    assert any(g["z_spec"] > 13 for g in jades)
    assert any(g["z_spec"] > 10 for g in ceers)


def test_cosmic_dawn_predictions_ordering(model):
    d = model.cosmic_dawn_predictions()
    assert d["z_first_stars"] > d["z_first_galaxies"] > d["z_50pct_reionised"]
    assert d["z_50pct_reionised"] > d["z_complete_reionisation"]


def test_tension_dict_nonempty(model):
    t = model.tension_with_LCDM()
    assert isinstance(t, dict)
    assert len(t) >= 4
    for v in t.values():
        assert isinstance(v, str) and len(v) > 0


def test_report_contains_key_sections(model):
    rep = model.report()
    assert "Halo mass function" in rep
    assert "Reionisation" in rep
    assert "Tension with LambdaCDM" in rep
    assert "JADES" in rep or "CEERS" in rep or "Maisie" in rep
