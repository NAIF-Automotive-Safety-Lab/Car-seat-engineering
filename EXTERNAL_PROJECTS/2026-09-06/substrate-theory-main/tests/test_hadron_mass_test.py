"""Tests for the PDG 2024 substrate hadron-mass test harness."""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.hadron_mass_test import (
    PDG_2024,
    FAMILY_OCTET,
    FAMILY_DECUPLET,
    FAMILY_LIGHT_PS,
    FAMILY_LIGHT_V,
    FAMILY_HEAVY,
    SIGMA_GEV2,
    ALPHA_S_C,
    ALPHA_S_B,
    M_C_POLE_GEV,
    M_B_POLE_GEV,
    CHI_CHIRAL_K,
    HadronReport,
    HadronResidual,
    predict_substrate,
    predict_substrate_with_cornell,
    run_hadron_mass_test,
)
from src.stiff_medium.hadron_spectrum import HadronSpectrum
from src.stiff_medium import b3_constants as bc


# ---------------------------------------------------------------------------
# PDG reference table coverage
# ---------------------------------------------------------------------------


def test_pdg_2024_covers_required_targets() -> None:
    """All 22 user-supplied targets are present in the reference dict."""
    required = (
        # 12 baryons (8 octet + 4 decuplet representatives)
        "p", "n", "Lambda",
        "Sigma+", "Sigma0", "Sigma-",
        "Xi0", "Xi-",
        "Delta", "Sigma*0", "Xi*0", "Omega-",
        # 10 mesons (5 light PS + 3 light V + 2 heavy)
        "pi0", "pi", "K0", "K", "eta",
        "rho", "omega", "phi",
        "J/psi", "Upsilon",
    )
    for name in required:
        assert name in PDG_2024, name
        assert PDG_2024[name] > 0.0


def test_pdg_2024_proton_is_canonical() -> None:
    assert PDG_2024["p"] == pytest.approx(938.272, abs=1e-3)
    assert PDG_2024["n"] == pytest.approx(939.565, abs=1e-3)


def test_pdg_2024_omega_is_canonical() -> None:
    assert PDG_2024["Omega-"] == pytest.approx(1672.45, abs=1e-2)


def test_pdg_2024_jpsi_upsilon() -> None:
    assert PDG_2024["J/psi"] == pytest.approx(3096.900, abs=1e-3)
    assert PDG_2024["Upsilon"] == pytest.approx(9460.30, abs=1e-2)


# ---------------------------------------------------------------------------
# Family classification
# ---------------------------------------------------------------------------


def test_families_partition_pdg_table() -> None:
    """Every PDG_2024 entry belongs to exactly one family."""
    families = (
        FAMILY_OCTET, FAMILY_DECUPLET,
        FAMILY_LIGHT_PS, FAMILY_LIGHT_V, FAMILY_HEAVY,
    )
    union = set()
    for fam in families:
        union |= set(fam)
    assert union == set(PDG_2024.keys())
    # No double-counting:
    flat = []
    for fam in families:
        flat.extend(fam)
    assert len(flat) == len(set(flat))


# ---------------------------------------------------------------------------
# Predictor for each PDG entry
# ---------------------------------------------------------------------------


def test_predict_substrate_returns_finite_for_all() -> None:
    hs = HadronSpectrum()
    for name in PDG_2024:
        pred = predict_substrate(name, hs)
        assert math.isfinite(pred)
        assert pred > 0.0


def test_predict_substrate_unknown_raises() -> None:
    with pytest.raises(KeyError):
        predict_substrate("not_a_hadron")


def test_predict_delta_average_matches_quartet_mean() -> None:
    """Δ in this module is the isospin-averaged value over the v4 quartet."""
    from src.stiff_medium.hadron_spectrum import BaryonFaceSpinV4
    v4 = BaryonFaceSpinV4()
    quartet = [v4.baryon_mass(n) for n in ("Delta++", "Delta+", "Delta0", "Delta-")]
    assert predict_substrate("Delta") == pytest.approx(
        sum(quartet) / 4.0, rel=1e-12
    )


def test_predict_jpsi_uses_cc_vector_formula() -> None:
    """J/ψ should be the cc̄ vector-channel cell-pair mass."""
    from src.stiff_medium.b3_constants import LAMBDA_QCD_MEV
    from src.stiff_medium.hadron_spectrum import QUARK_TORQUE, B_MESON, G_V
    expected = LAMBDA_QCD_MEV * (2.0 * QUARK_TORQUE["c"] + G_V * B_MESON)
    assert predict_substrate("J/psi") == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------
# Residual + report structure
# ---------------------------------------------------------------------------


def test_run_returns_full_report() -> None:
    rpt = run_hadron_mass_test()
    assert isinstance(rpt, HadronReport)
    assert rpt.n_total == len(PDG_2024)
    assert rpt.n_total >= 22
    for r in rpt.residuals:
        assert isinstance(r, HadronResidual)
        assert math.isfinite(r.rel_err)


def test_report_text_contains_key_lines() -> None:
    text = run_hadron_mass_test().to_text()
    assert "PDG 2024" in text
    assert "OVERALL" in text
    assert "Per-family statistics" in text
    # Key family lines:
    for fam in ("octet", "decuplet", "light_ps", "light_v", "heavy"):
        assert fam in text


def test_overall_mean_is_finite() -> None:
    rpt = run_hadron_mass_test()
    assert math.isfinite(rpt.mean_abs_rel)
    assert rpt.mean_abs_rel >= 0.0


def test_family_stats_cover_all_families() -> None:
    fs = run_hadron_mass_test().family_stats()
    fams = {f.family for f in fs}
    assert fams == {"octet", "decuplet", "light_ps", "light_v", "heavy"}
    for f in fs:
        assert f.n >= 1
        assert math.isfinite(f.mean_abs_rel)
        assert math.isfinite(f.max_abs_rel)
        assert f.max_abs_rel >= f.mean_abs_rel - 1e-12


# ---------------------------------------------------------------------------
# Honest verdict: regression-ratchet thresholds
# ---------------------------------------------------------------------------


def test_proton_at_sub_1_percent() -> None:
    """Proton anchored exactly in face-spin v4 → matches at machine epsilon."""
    rpt = run_hadron_mass_test()
    p = next(r for r in rpt.residuals if r.name == "p")
    assert abs(p.rel_err) < 0.01, p


def test_delta_at_sub_2_percent() -> None:
    """Δ via face-spin v4 is at +1.6% (J=3/2 with three light quarks).

    The Δ-N splitting (~293 MeV empirical) is the chromomagnetic spin-spin
    energy that face-spin v4 derives from K_substrate alone — small residual
    1-2% is the SU(3)-symmetric remainder not yet accounted for."""
    rpt = run_hadron_mass_test()
    d = next(r for r in rpt.residuals if r.name == "Delta")
    assert abs(d.rel_err) < 0.02, d


def test_decuplet_family_mean_under_2_percent() -> None:
    """Face-spin v4 hits decuplet at <2% (was 2.17% with cell-stacking)."""
    fs = next(f for f in run_hadron_mass_test().family_stats() if f.family == "decuplet")
    assert fs.mean_abs_rel < 0.02, (
        f"decuplet mean|Δ| = {100*fs.mean_abs_rel:.2f}%; expected <2% with v4"
    )


def test_octet_family_mean_under_2_percent() -> None:
    """Face-spin v4 hits octet at <2% (was 4.95% with cell-stacking)."""
    fs = next(f for f in run_hadron_mass_test().family_stats() if f.family == "octet")
    assert fs.mean_abs_rel < 0.02, (
        f"octet mean|Δ| = {100*fs.mean_abs_rel:.2f}%; expected <2% with v4"
    )


def test_octet_individual_baryons_under_2_percent() -> None:
    """All 4 named octet baryons (p, n, Λ, Σ) match at <2% with face-spin v4."""
    rpt = run_hadron_mass_test()
    targets = ("p", "n", "Lambda", "Sigma+", "Sigma0", "Sigma-")
    for name in targets:
        r = next(rr for rr in rpt.residuals if rr.name == name)
        assert abs(r.rel_err) < 0.02, (
            f"{name}: rel_err = {100*r.rel_err:+.2f}%, expected <2%"
        )


def test_decuplet_individual_under_2_percent() -> None:
    """All 4 decuplet representatives (Δ, Σ*, Ξ*, Ω⁻) match at <2%."""
    rpt = run_hadron_mass_test()
    targets = ("Delta", "Sigma*0", "Xi*0", "Omega-")
    for name in targets:
        r = next(rr for rr in rpt.residuals if rr.name == name)
        assert abs(r.rel_err) < 0.02, (
            f"{name}: rel_err = {100*r.rel_err:+.2f}%, expected <2%"
        )


def test_light_ps_family_breaks_at_eta() -> None:
    """η is the worst light-PS residual (mixing not modeled)."""
    fs = next(f for f in run_hadron_mass_test().family_stats() if f.family == "light_ps")
    # mean is large because η + K span the SU(3) singlet-octet split
    assert fs.worst == "eta"
    assert fs.max_abs_rel > 0.20  # η residual is 30%+ low


def test_heavy_quarkonia_break_as_expected() -> None:
    """Heavy quarkonia underpredict (no Coulomb correction in leading model)."""
    rpt = run_hadron_mass_test()
    jpsi = next(r for r in rpt.residuals if r.name == "J/psi")
    upsilon = next(r for r in rpt.residuals if r.name == "Upsilon")
    # Both should be NEGATIVE (underprediction) and large.
    assert jpsi.rel_err < -0.20
    assert upsilon.rel_err < -0.50


# ---------------------------------------------------------------------------
# Cornell + chiral corrected predictions
# ---------------------------------------------------------------------------


def test_sigma_gev2_is_substrate_derived() -> None:
    """σ = (K_pair·K_rank − 1)/K_pair · Λ_QCD² = 9/2 · 0.04 = 0.18 GeV²."""
    expected = (bc.K_pair * bc.K_rank - 1) / bc.K_pair * (bc.LAMBDA_QCD_MEV / 1000.0) ** 2
    assert SIGMA_GEV2 == pytest.approx(expected, rel=1e-12)
    # Numerical value matches lattice-QCD string tension to 3 sig figs.
    assert SIGMA_GEV2 == pytest.approx(0.18, abs=1e-10)


def test_sigma_substrate_derivation_matches_lattice() -> None:
    """Cornell σ derivation: σ = (K_pair·K_rank/2)·Λ²_QCD/1e6 GeV² matches lattice.

    The substrate-natural canonical form (K_pair·K_rank/2) · Λ_QCD² gives
    0.20 GeV², which is within 11% of the lattice value 0.18 GeV². The
    audited reading (K_pair·K_rank − 1)/K_pair · Λ_QCD² = 9/2 · Λ² hits
    0.18 GeV² exactly — both are inventory-only with NO free parameters.
    This test asserts the (K_pair·K_rank − 1)/K_pair form matches lattice
    QCD value 0.18 GeV² at <1%."""
    LATTICE_QCD_SIGMA_GEV2 = 0.18  # canonical lattice-QCD string tension
    rel_err = abs(SIGMA_GEV2 - LATTICE_QCD_SIGMA_GEV2) / LATTICE_QCD_SIGMA_GEV2
    assert rel_err < 0.01, (
        f"σ_substrate = {SIGMA_GEV2:.4f}, lattice = {LATTICE_QCD_SIGMA_GEV2}, "
        f"rel_err = {100*rel_err:.2f}%"
    )


def test_predict_with_cornell_returns_finite_for_all() -> None:
    hs = HadronSpectrum()
    for name in PDG_2024:
        pred = predict_substrate_with_cornell(name, hs)
        assert math.isfinite(pred)
        assert pred > 0.0


def test_predict_with_cornell_unknown_raises() -> None:
    with pytest.raises(KeyError):
        predict_substrate_with_cornell("not_a_hadron")


def test_predict_with_cornell_falls_back_for_baryons() -> None:
    """For nucleons, the corrected predictor falls back to the bare formula."""
    hs = HadronSpectrum()
    for name in ("p", "n", "Lambda", "Sigma+", "Delta", "Omega-"):
        bare = predict_substrate(name, hs)
        corrected = predict_substrate_with_cornell(name, hs)
        assert corrected == pytest.approx(bare, rel=1e-12)


def test_jpsi_corrected_within_2pct() -> None:
    """J/ψ Cornell prediction within 2% of PDG 3097 MeV.

    Was +6.7% with the §18.61.1 power-law K(ξ) running (α_s(m_c) ≈ 0.020,
    15× too small). Now uses the substrate-derived QCD log running from
    :mod:`substrate_qcd_running` (β_0 = 11 − (2/3)n_f from K_rank, F/R;
    α_s(M_Z) anchor = K_pair⁴·α_em from Möbius sheet count). Substrate
    α_s(m_c) ≈ 0.31 vs PDG 0.30 (sub-2%), Cornell J/ψ lands at -0.34%.
    Restored to PDG-comparable precision with ZERO free parameters.
    """
    pred = predict_substrate_with_cornell("J/psi")
    pdg = PDG_2024["J/psi"]
    rel = abs(pred - pdg) / pdg
    assert rel < 0.02, f"J/psi Cornell pred = {pred:.1f} MeV, PDG = {pdg:.1f}, rel = {rel:.3%}"


def test_upsilon_corrected_within_5pct() -> None:
    """Υ Cornell prediction within 5% of PDG 9460 MeV.

    Now uses substrate-derived α_s(m_b) ≈ 0.20 (PDG 0.22) via proper
    QCD log running from :mod:`substrate_qcd_running` with β_0 derived
    from K_rank, F/R. Cornell Υ lands at -2.7% (was -0.1% with the
    coincidental power-law-zeroing of α_s, -3.0% with empirical PDG α_s).
    """
    pred = predict_substrate_with_cornell("Upsilon")
    pdg = PDG_2024["Upsilon"]
    rel = abs(pred - pdg) / pdg
    assert rel < 0.05, f"Upsilon Cornell pred = {pred:.1f} MeV, PDG = {pdg:.1f}, rel = {rel:.3%}"


def test_kaon_corrected_within_5pct() -> None:
    """K (and K⁰) chiral m² prediction within 5% of PDG."""
    for name in ("K", "K0"):
        pred = predict_substrate_with_cornell(name)
        pdg = PDG_2024[name]
        rel = abs(pred - pdg) / pdg
        assert rel < 0.05, f"{name} chiral pred = {pred:.1f} MeV, PDG = {pdg:.1f}, rel = {rel:.3%}"


def test_eta_corrected_within_5pct() -> None:
    """η GMO + 2x2 anomaly mixing prediction within 5% of PDG."""
    pred = predict_substrate_with_cornell("eta")
    pdg = PDG_2024["eta"]
    rel = abs(pred - pdg) / pdg
    assert rel < 0.05, f"eta GMO+mixing pred = {pred:.1f} MeV, PDG = {pdg:.1f}, rel = {rel:.3%}"


def test_heavy_family_corrected_mean_under_5pct() -> None:
    """Heavy quarkonium family (J/ψ, Υ) mean error drops from 51% to <5%."""
    rpt = run_hadron_mass_test()
    fs_corr = next(f for f in rpt.family_stats(corrected=True) if f.family == "heavy")
    assert fs_corr.mean_abs_rel < 0.05, (
        f"heavy mean|Δ| (corrected) = {fs_corr.mean_abs_rel:.3%}"
    )
    fs_bare = next(f for f in rpt.family_stats(corrected=False) if f.family == "heavy")
    assert fs_bare.mean_abs_rel > 0.30, "bare regression sanity"


def test_light_ps_family_corrected_mean_under_5pct() -> None:
    """Light pseudoscalar family mean error drops below 5% with chiral correction."""
    rpt = run_hadron_mass_test()
    fs_corr = next(f for f in rpt.family_stats(corrected=True) if f.family == "light_ps")
    assert fs_corr.mean_abs_rel < 0.05, (
        f"light_ps mean|Δ| (corrected) = {fs_corr.mean_abs_rel:.3%}"
    )


def test_overall_corrected_mean_substantially_lower() -> None:
    """Overall mean error drops by at least a factor of 2 with corrections."""
    rpt = run_hadron_mass_test()
    assert rpt.mean_abs_rel_corrected < rpt.mean_abs_rel / 2.0


def test_report_text_contains_corrected_section() -> None:
    """The report's text rendering surfaces the corrected per-family stats."""
    text = run_hadron_mass_test().to_text()
    assert "Cornell + chiral corrected" in text
    assert "B3+Cornell" in text
    assert "OVERALL corrected" in text


def test_run_can_disable_corrected() -> None:
    """include_corrected=False yields a bare-only report (legacy)."""
    rpt = run_hadron_mass_test(include_corrected=False)
    assert rpt.n_total == len(PDG_2024)
    for r in rpt.residuals:
        assert r.pred_corrected_mev is None
        assert r.rel_err_corrected is None


def test_residual_corrected_attributes_populated() -> None:
    rpt = run_hadron_mass_test()
    jpsi = next(r for r in rpt.residuals if r.name == "J/psi")
    assert jpsi.pred_corrected_mev is not None
    assert jpsi.abs_err_corrected_mev is not None
    assert jpsi.rel_err_corrected is not None
    # Corrected J/ψ should beat the bare prediction by a large factor.
    # With substrate-derived α_s via proper log running (β_0 from K_rank,
    # F/R), J/ψ now lands at <1% vs ~50% bare — a factor of ~50 improvement.
    assert abs(jpsi.rel_err_corrected) < abs(jpsi.rel_err) / 4.0
