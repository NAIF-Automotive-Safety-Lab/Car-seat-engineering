"""Tests for stiff_medium.neutrino_oscillation."""

from __future__ import annotations

import numpy as np
import pytest

from stiff_medium.neutrino_oscillation import (
    NeutrinoOscillation,
    OscillationConfig,
    DELTA_M2_21,
    DELTA_M2_31_NH,
    M1_SUBSTRATE_EV,
)


@pytest.fixture(scope="module")
def osc() -> NeutrinoOscillation:
    return NeutrinoOscillation()


# ---------------------------------------------------------------------
# (a) PMNS unitarity
# ---------------------------------------------------------------------
def test_pmns_unitary(osc: NeutrinoOscillation) -> None:
    U = osc.pmns_matrix()
    assert U.shape == (3, 3)
    delta = U @ U.conj().T - np.eye(3)
    assert np.max(np.abs(delta)) < 1e-12


def test_pmns_columns_orthonormal(osc: NeutrinoOscillation) -> None:
    U = osc.pmns_matrix()
    delta = U.conj().T @ U - np.eye(3)
    assert np.max(np.abs(delta)) < 1e-12


def test_probabilities_sum_to_one(osc: NeutrinoOscillation) -> None:
    # Σ_β P(ν_α → ν_β) = 1 for each α
    for alpha in osc.flavours:
        s = sum(osc.oscillation_probability(alpha, beta, 0.6, 295.0)
                for beta in osc.flavours)
        assert s == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------
# Mass eigenvalues
# ---------------------------------------------------------------------
def test_mass_NH_consistency(osc: NeutrinoOscillation) -> None:
    m1, m2, m3 = osc.mass_eigenvalues_NH()
    assert m1 == pytest.approx(M1_SUBSTRATE_EV, rel=1e-12)
    assert m2 ** 2 - m1 ** 2 == pytest.approx(DELTA_M2_21, rel=1e-9)
    assert m3 ** 2 - m1 ** 2 == pytest.approx(DELTA_M2_31_NH, rel=1e-9)


def test_mass_IH_floor(osc: NeutrinoOscillation) -> None:
    # IH has hard floor Σm_ν ≳ 100 meV (incompatible with substrate 60.5 meV)
    s = osc.sigma_m_nu("IH")
    assert s > 0.095   # eV
    assert s < 0.110   # eV  (m3 = 0 baseline)
    m1, m2, m3 = osc.mass_eigenvalues_IH()
    assert m3 < m1 < m2
    assert m2 ** 2 - m1 ** 2 == pytest.approx(DELTA_M2_21, rel=1e-9)


# ---------------------------------------------------------------------
# (b) Atmospheric oscillation peak at L/E ~ 500 km/GeV
# ---------------------------------------------------------------------
def test_atmospheric_peak_at_500_km_per_GeV(osc: NeutrinoOscillation) -> None:
    res = osc.atmospheric_neutrino_oscillation()
    LoverE = res["first_min_L_over_E_km_per_GeV"]
    # First minimum: 1.267 · Δm²_31 · L/E = π/2  ⇒  L/E ≈ 489 km/GeV
    expected = (np.pi / 2.0) / (1.26693281 * DELTA_M2_31_NH)  # ~489 km/GeV
    assert LoverE == pytest.approx(expected, rel=0.05)
    # Disappearance is strong at the minimum
    assert res["P_mu_mu_at_min"] < 0.10


# ---------------------------------------------------------------------
# (c) Solar νe survival probability ~0.55 at low E (pp), ~0.30 at high E
# ---------------------------------------------------------------------
def test_solar_low_energy_survival(osc: NeutrinoOscillation) -> None:
    sol = osc.solar_neutrino_problem(E_MeV=0.3)
    # Vacuum-averaged for pp neutrinos: ~0.55 (Super-K / SNO consistent).
    assert 0.50 <= sol["P_ee_low_E"] <= 0.60


def test_solar_high_energy_survival(osc: NeutrinoOscillation) -> None:
    sol = osc.solar_neutrino_problem(E_MeV=10.0)
    # Matter-dominated ⁸B: ~0.30 (SNO).
    assert 0.25 <= sol["P_ee_high_E"] <= 0.35


# ---------------------------------------------------------------------
# T2K / DUNE smoke tests
# ---------------------------------------------------------------------
def test_T2K_appearance_nonzero(osc: NeutrinoOscillation) -> None:
    t = osc.T2K_results()
    assert 0.0 < t["P_mu_e"] < 0.20
    assert 0.0 < t["P_mu_mu"] < 1.0


def test_T2K_CP_asymmetry_sign_and_magnitude(osc: NeutrinoOscillation) -> None:
    t = osc.T2K_results()
    assert abs(t["CP_asymmetry"]) <= 1.0
    # δ_CP = 3π/4 ≠ 0, π → asymmetry should be measurably non-zero
    assert abs(t["CP_asymmetry"]) > 0.01


def test_dune_predictions_shape(osc: NeutrinoOscillation) -> None:
    d = osc.dune_predictions()
    n = len(d["E_GeV"])
    assert n == len(d["P_mu_e_NH"]) == len(d["P_mu_mu_NH"])
    assert n == len(d["CP_asymmetry"]) == len(d["P_mu_e_IH"])
    assert d["max_CP_asymmetry"] > 0.0


# ---------------------------------------------------------------------
# Hierarchy distinguishability
# ---------------------------------------------------------------------
def test_NH_IH_differ(osc: NeutrinoOscillation) -> None:
    P_NH = osc.oscillation_probability("mu", "e", 2.0, 1300.0, hierarchy="NH")
    P_IH = osc.oscillation_probability("mu", "e", 2.0, 1300.0, hierarchy="IH")
    assert P_NH != pytest.approx(P_IH, abs=1e-3)


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------
def test_report_returns_keys(osc: NeutrinoOscillation, capsys) -> None:
    out = osc.report()
    keys_required = {
        "hierarchy", "delta_CP_rad",
        "mass_eigenvalues_NH_eV", "mass_eigenvalues_IH_eV",
        "Sigma_m_nu_NH_eV", "Sigma_m_nu_IH_eV",
        "PMNS_unitarity_residual",
        "atm_first_min_L_over_E", "atm_P_mu_mu_at_min",
        "solar_P_ee_high_E", "solar_P_ee_low_E",
        "T2K_P_mu_e", "T2K_P_mu_mu", "T2K_CP_asymmetry",
    }
    assert keys_required.issubset(out.keys())
    captured = capsys.readouterr()
    assert "B3 Neutrino Oscillation Report" in captured.out
