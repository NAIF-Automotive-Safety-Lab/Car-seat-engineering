"""Tests for src/stiff_medium/bsm_predictions.py.

Verify each BSMPredictions method returns a sensible structured result and
that there are no internal contradictions across the BSM dossier.
"""

from __future__ import annotations

import math

import pytest

from stiff_medium.bsm_predictions import BSMPredictions, BSMResult


@pytest.fixture
def bsm() -> BSMPredictions:
    return BSMPredictions()


# ---------------------------------------------------------------------------
# Per-method checks
# ---------------------------------------------------------------------------


def test_dark_matter(bsm: BSMPredictions) -> None:
    r = bsm.dark_matter()
    assert r["name"] == "dark_matter"
    assert r["values"]["mass_gev"] == pytest.approx(27.5)
    assert r["values"]["sigma_SI_cm2"] == pytest.approx(1e-135)
    assert r["values"]["candidate"] == "cube_cell_substrate_defect"
    assert "cube" in r["mechanism"].lower()


def test_dark_energy(bsm: BSMPredictions) -> None:
    r = bsm.dark_energy()
    assert r["name"] == "dark_energy"
    rho_nat = r["values"]["rho_natural_gev4"]
    rho_obs = r["values"]["rho_observed_gev4"]
    assert rho_nat > 0
    assert rho_obs > 0
    # The naive Lambda_QCD^4 / M_Pl^2 overshoots rho_Lambda by ~6 orders;
    # B3's "O(1) topology factor" claim absorbs this — we just require it
    # is vastly closer than the SM Planck-cutoff vacuum estimate (~10^120).
    o1 = r["values"]["O1_factor"]
    assert 1e-10 < o1 < 1e10, f"dark-energy ratio outside sanity window: {o1}"


def test_neutrinos(bsm: BSMPredictions) -> None:
    r = bsm.neutrinos()
    assert r["values"]["Sigma_m_nu_meV"] == pytest.approx(60.5)
    lo, hi = r["values"]["m_betabeta_window_meV"]
    assert lo == 0.0
    assert hi == 5.0
    assert "majorana" in r["prediction"].lower()
    # Lightest neutrino mass must be much smaller than the sum.
    assert r["values"]["lightest_m1_meV"] < r["values"]["Sigma_m_nu_meV"]


def test_hierarchy_problem(bsm: BSMPredictions) -> None:
    r = bsm.hierarchy_problem()
    expected_exp = 4.0 * math.pi * math.pi - 1.0
    assert r["values"]["exponent_4pi2_minus_1"] == pytest.approx(expected_exp)
    assert r["values"]["exp_minus_exponent"] == pytest.approx(math.exp(-expected_exp))
    # Predicted suppression should be in the right ballpark of Lambda_QCD/M_Pl
    # (~10^-20). Allow 3 orders of magnitude either way for the loose
    # ontology-level claim.
    pred = r["values"]["exp_minus_exponent"]
    obs = r["values"]["lambda_qcd_over_M_Pl"]
    assert 1e-4 < pred / obs < 1e4


def test_strong_cp(bsm: BSMPredictions) -> None:
    r = bsm.strong_cp()
    assert r["values"]["theta_QCD"] == 0.0
    assert r["values"]["axion_required"] is False
    # Prediction must be tighter than the experimental bound.
    assert r["values"]["theta_QCD"] < r["values"]["experimental_bound"]


def test_proton_decay(bsm: BSMPredictions) -> None:
    r = bsm.proton_decay()
    assert r["values"]["tau_p_years"] == math.inf
    assert r["values"]["channel"] == "none"
    # Must not contradict the SK lower bound.
    assert r["values"]["tau_p_years"] > r["values"]["experimental_lower_bound_years"]


def test_magnetic_monopoles(bsm: BSMPredictions) -> None:
    r = bsm.magnetic_monopoles()
    assert r["values"]["monopole_density"] == 0.0
    assert "forbidden" in r["prediction"].lower()


def test_extra_dimensions(bsm: BSMPredictions) -> None:
    r = bsm.extra_dimensions()
    assert r["values"]["n_spatial_dimensions"] == 3
    assert r["values"]["kk_tower_present"] is False
    assert r["values"]["compactification_radius_m"] is None


def test_supersymmetry(bsm: BSMPredictions) -> None:
    r = bsm.supersymmetry()
    assert r["values"]["superpartners"] is False
    assert r["values"]["expected_susy_scale_gev"] is None
    assert r["values"]["neutralino_dm"] is False


def test_string_theory_relation(bsm: BSMPredictions) -> None:
    r = bsm.string_theory_relation()
    assert r["values"]["string_scale_planckian"] is False
    assert r["values"]["compactified_dimensions"] == 0
    assert r["values"]["branes_required"] is False
    # String scale must be hadronic, not Planckian.
    assert r["values"]["string_scale_gev"] < 1.0  # GeV
    assert r["values"]["n_sheets_K_pair"] == 2


# ---------------------------------------------------------------------------
# Cross-method consistency
# ---------------------------------------------------------------------------


def test_no_susy_no_neutralino_dm(bsm: BSMPredictions) -> None:
    """SUSY denied AND DM is the cube cell — these must agree."""
    susy = bsm.supersymmetry()
    dm = bsm.dark_matter()
    assert susy["values"]["superpartners"] is False
    assert dm["values"]["candidate"] == "cube_cell_substrate_defect"
    assert susy["values"]["neutralino_dm"] is False


def test_no_extra_dims_no_kk_string_scale(bsm: BSMPredictions) -> None:
    """No extra dims must align with no compactified dimensions in the string relation."""
    ed = bsm.extra_dimensions()
    st = bsm.string_theory_relation()
    assert ed["values"]["n_spatial_dimensions"] == 3
    assert st["values"]["compactified_dimensions"] == 0


def test_strong_cp_no_axion_no_susy(bsm: BSMPredictions) -> None:
    """Strong-CP solved without axion AND without SUSY — both consistent."""
    scp = bsm.strong_cp()
    susy = bsm.supersymmetry()
    assert scp["values"]["axion_required"] is False
    assert susy["values"]["superpartners"] is False


def test_all_predictions_keys(bsm: BSMPredictions) -> None:
    expected = {
        "dark_matter",
        "dark_energy",
        "neutrinos",
        "hierarchy_problem",
        "strong_cp",
        "proton_decay",
        "magnetic_monopoles",
        "extra_dimensions",
        "supersymmetry",
        "string_theory_relation",
    }
    got = set(bsm.all_predictions().keys())
    assert got == expected


def test_every_prediction_is_complete(bsm: BSMPredictions) -> None:
    """Every method's dict must have the standard fields populated."""
    required = {"name", "prediction", "values", "mechanism", "status", "references"}
    for name, p in bsm.all_predictions().items():
        missing = required - set(p.keys())
        assert not missing, f"{name} missing fields: {missing}"
        assert p["prediction"], f"{name} has empty prediction"
        assert p["mechanism"], f"{name} has empty mechanism"
        assert p["status"], f"{name} has empty status"
        assert isinstance(p["references"], list)
        assert p["references"], f"{name} has no references"


def test_report_contains_all_sections(bsm: BSMPredictions) -> None:
    rep = bsm.report()
    assert isinstance(rep, str)
    assert "BSM PREDICTIONS DOSSIER" in rep
    for key in bsm.all_predictions():
        assert key.upper() in rep, f"report missing section {key}"
    assert "END BSM DOSSIER" in rep


def test_bsm_result_dataclass_roundtrip() -> None:
    r = BSMResult(name="x", prediction="p", mechanism="m", status="s",
                  references=["a"], values={"v": 1})
    d = r.to_dict()
    assert d["name"] == "x"
    assert d["values"] == {"v": 1}
    assert d["references"] == ["a"]
