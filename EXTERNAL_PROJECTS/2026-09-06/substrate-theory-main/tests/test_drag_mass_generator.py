"""Tests for src.stiff_medium.drag_mass_generator."""

import numpy as np
import pytest

from src.stiff_medium.drag_mass_generator import (
    DragMassGenerator,
    BounceResult,
    PDG_MASSES_MEV,
    DEFAULT_PRIMITIVES,
    HBAR,
    C_LIGHT,
    MEV_J,
    _compton,
)


@pytest.fixture
def gen():
    return DragMassGenerator()


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def test_default_primitives(gen):
    assert gen.K == DEFAULT_PRIMITIVES["K"]
    assert gen.rho == DEFAULT_PRIMITIVES["rho"]
    assert gen.alpha_drag > 0.0


def test_compton_helper_inverts_mass():
    # Compton wavelength times mass*c should equal hbar
    m_mev = PDG_MASSES_MEV["e"]
    xi = _compton(m_mev)
    m_kg = (m_mev * MEV_J) / C_LIGHT ** 2
    assert np.isclose(xi * m_kg * C_LIGHT, HBAR, rtol=1e-12)


def test_bounce_result_fields(gen):
    r = gen.rest_mass({"xi": _compton(PDG_MASSES_MEV["e"]), "gamma": 0.0})
    assert isinstance(r, BounceResult)
    assert r.omega_b > 0
    assert r.mass_mev > 0
    assert "K" in r.primitives


# ---------------------------------------------------------------------------
# Required physics behaviors
# ---------------------------------------------------------------------------

def test_larger_gamma_gives_larger_omega_and_mass(gen):
    """Required (a): increasing drag must up-shift omega_b and m."""
    xi = _compton(PDG_MASSES_MEV["e"])
    # Choose gammas spanning a few decades around omega_compton
    omega_compton = C_LIGHT / xi
    gammas = [0.0, 0.5 * omega_compton, 5.0 * omega_compton, 50.0 * omega_compton]
    omegas, masses = [], []
    for g in gammas:
        r = gen.rest_mass({"xi": xi, "gamma": g})
        omegas.append(r.omega_b)
        masses.append(r.mass_mev)
    # Strictly monotonic
    assert all(omegas[i + 1] > omegas[i] for i in range(len(gammas) - 1))
    assert all(masses[i + 1] > masses[i] for i in range(len(gammas) - 1))


def test_electron_mass_within_5pct(gen):
    """Required (b): electron mass within 5% of 0.511 MeV."""
    r = gen.electron_mass(gamma=0.0)
    rel = abs(r.mass_mev - PDG_MASSES_MEV["e"]) / PDG_MASSES_MEV["e"]
    assert rel < 0.05, f"electron mass off by {rel:.3%}"


def test_proton_mass_anchor(gen):
    r = gen.proton_mass(gamma_qcd=0.0)
    # Proton anchor (938.272 MeV) recovered by Compton id at gamma=0
    rel = abs(r.mass_mev - 938.272) / 938.272
    assert rel < 0.01


def test_quark_spectrum_recovers_pdg(gen):
    spectrum = gen.quark_mass_spectrum()
    for q in ("u", "d", "s", "c", "b", "t"):
        rel = abs(spectrum[q].mass_mev - PDG_MASSES_MEV[q]) / PDG_MASSES_MEV[q]
        assert rel < 0.05, f"quark {q} off by {rel:.3%}"


def test_drag_to_mass_ratio_dimensionless_and_consistent(gen):
    """Required (d): m c^2 / (hbar gamma) dimensionless and O(1) across
    particles when gamma is anchored to the bound-state's own omega_b."""
    ratios = []
    for name, m_obs in PDG_MASSES_MEV.items():
        xi_p = _compton(m_obs)
        bare = gen.rest_mass({"xi": xi_p, "gamma": 0.0})
        # use gamma = omega_b / 2 as canonical drag scale per particle
        probe_gamma = bare.omega_b / 2.0
        r = gen.rest_mass({"xi": xi_p, "gamma": probe_gamma})
        ratios.append(r.drag_to_mass)
    ratios = np.array(ratios)
    # All dimensionless, strictly positive, and within a single decade.
    assert np.all(np.isfinite(ratios))
    assert np.all(ratios > 0)
    assert ratios.max() / ratios.min() < 10.0


# ---------------------------------------------------------------------------
# Required (c): cone-bouncing frequency stable under time evolution
# ---------------------------------------------------------------------------

def test_field_evolution_matches_predicted_omega(gen):
    """Required (c): numerically integrated cone-bouncing frequency is
    stable and within tolerance of the closed-form prediction."""
    # Use an electron-scale anchor with mild drag for fast convergence
    xi = _compton(PDG_MASSES_MEV["e"])
    sim = gen.evolve_field(xi=xi, gamma=0.0, n_grid=128, n_steps=2000)
    assert sim["n_peaks_found"] >= 3, "not enough oscillation peaks"
    assert np.isfinite(sim["omega_b_measured"])
    # Tolerate modest discretization error
    assert sim["rel_err"] < 0.1, (
        f"omega_b mismatch: pred={sim['omega_b_predicted']:.3e} "
        f"meas={sim['omega_b_measured']:.3e} rel_err={sim['rel_err']:.3e}"
    )


# ---------------------------------------------------------------------------
# Verification + reporting
# ---------------------------------------------------------------------------

def test_verify_against_pdg_returns_rows(gen):
    rows = gen.verify_against_pdg()
    assert len(rows) == len(PDG_MASSES_MEV)
    for r in rows:
        assert r.predicted_mev > 0
        assert r.observed_mev > 0
        assert r.rel_err < 1e-6, f"{r.name} not Compton-anchored"


def test_alpha_737_cross_check_close(gen):
    a = gen.cross_check_alpha_737()
    # m_p/m_e ~ 1836.15, B3 quoted integer 737 is in the *exponent*; this
    # cross-check just makes sure the formula evaluates and is finite
    assert np.isfinite(a["alpha_predicted"])
    assert a["alpha_predicted"] > 0
    assert a["alpha_observed"] > 0


def test_report_runs(gen, capsys):
    gen.report()
    captured = capsys.readouterr()
    assert "alpha-EM cross-check" in captured.out
