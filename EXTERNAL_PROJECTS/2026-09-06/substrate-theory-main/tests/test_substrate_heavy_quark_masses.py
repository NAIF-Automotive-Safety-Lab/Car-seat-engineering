"""Tests for substrate-derived heavy-quark pole masses (m_c, m_b).

The form

    m_pole = (K_pair² / K_rank) · T_q^(n_R / (K_pair · K_rank)) · Λ_QCD
           = (4/5) · T_q^(9/5) · 200 MeV

is substrate-composite (Cat-A*): zero free parameters, integer-rigid,
single form fits BOTH heavy quarks. The exponent 9/5 lacks an independent
derivation from substrate dynamics — open research target.

Tests check:
  * m_c, m_b match PDG MS-bar to <1%
  * Form is integer-rigid (perturbations break the prediction)
  * Alternative forms (additive, pure α_s, K_4 closure) all fail
  * Wilson conversion to pole scheme works with substrate α_s
  * Module-level constants stay self-consistent
"""

from __future__ import annotations

import math

import pytest

from src.stiff_medium import b3_constants as bc
from src.stiff_medium.hadron_spectrum import QUARK_TORQUE
from src.stiff_medium.substrate_heavy_quark_masses import (
    POLE_PREFACTOR_A,
    POLE_EXPONENT_P,
    heavy_quark_pole_mass_MeV,
    m_c_pole_substrate_GeV,
    m_b_pole_substrate_GeV,
    m_c_pole_wilson_substrate_GeV,
    m_b_pole_wilson_substrate_GeV,
    alternative_forms_tried,
    rigidity_under_integer_shifts,
    heavy_quark_predictions,
    compute_full_report,
    format_report,
)


# ---------------------------------------------------------------------------
# Form constants
# ---------------------------------------------------------------------------

def test_pole_prefactor_is_K_pair_sq_over_K_rank() -> None:
    """a = K_pair² / K_rank = 4/5 = 0.8, exact integer ratio."""
    assert POLE_PREFACTOR_A == bc.K_pair ** 2 / bc.K_rank
    assert math.isclose(POLE_PREFACTOR_A, 0.8, rel_tol=1e-12)


def test_pole_exponent_is_nR_over_KpairKrank() -> None:
    """p = n_R / (K_pair · K_rank) = 18/10 = 9/5 = 1.8, exact integer ratio."""
    assert POLE_EXPONENT_P == bc.n_R / (bc.K_pair * bc.K_rank)
    assert math.isclose(POLE_EXPONENT_P, 1.8, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Per-quark mass predictions match PDG MS-bar to <1%
# ---------------------------------------------------------------------------

def test_m_c_pole_matches_PDG_MSbar_under_1pct() -> None:
    """Substrate m_c MS-bar = 1.275 GeV exactly matches PDG."""
    m_c_GeV = m_c_pole_substrate_GeV()
    pdg_MSbar = 1.275
    assert abs(m_c_GeV - pdg_MSbar) / pdg_MSbar < 0.01


def test_m_b_pole_matches_PDG_MSbar_under_1pct() -> None:
    """Substrate m_b MS-bar ≈ 4.20 GeV matches PDG 4.18 GeV at <1%."""
    m_b_GeV = m_b_pole_substrate_GeV()
    pdg_MSbar = 4.18
    assert abs(m_b_GeV - pdg_MSbar) / pdg_MSbar < 0.01


def test_m_c_explicit_value() -> None:
    """Numerical regression: m_c = (4/5)·3.168^(9/5)·200 MeV ≈ 1275.0 MeV."""
    m_c_MeV = heavy_quark_pole_mass_MeV(QUARK_TORQUE["c"])
    assert math.isclose(m_c_MeV, 1275.0, abs_tol=2.0)  # within 2 MeV


def test_m_b_explicit_value() -> None:
    """Numerical regression: m_b = (4/5)·6.146^(9/5)·200 MeV ≈ 4203 MeV."""
    m_b_MeV = heavy_quark_pole_mass_MeV(QUARK_TORQUE["b"])
    assert math.isclose(m_b_MeV, 4203.0, abs_tol=5.0)  # within 5 MeV


# ---------------------------------------------------------------------------
# Single form fits BOTH heavy quarks at <1%
# ---------------------------------------------------------------------------

def test_single_form_fits_both_heavy_quarks() -> None:
    """No per-quark tuning: same (a, p) gives <1% on c AND b."""
    preds = heavy_quark_predictions()
    for name, p in preds.items():
        assert abs(p.rel_err) < 0.01, (
            f"{name}: rel_err = {p.rel_err*100:+.2f}%, must be < 1%"
        )


# ---------------------------------------------------------------------------
# Rigidity: integer shifts break the prediction
# ---------------------------------------------------------------------------

def test_canonical_perturbation_passes() -> None:
    """At canonical integers, prediction matches PDG."""
    rig = rigidity_under_integer_shifts()
    canonical = next(r for r in rig if "canonical" in r.perturbation)
    assert abs(canonical.rel_err_c) < 0.01
    assert abs(canonical.rel_err_b) < 0.01


def test_K_pair_perturbation_breaks_prediction() -> None:
    """Shift K_pair: 2 → 3 must break m_c, m_b by >5%."""
    rig = rigidity_under_integer_shifts()
    K_pair_shift = next(r for r in rig if "K_pair: 2 → 3" in r.perturbation)
    assert abs(K_pair_shift.rel_err_c) > 0.05 or abs(K_pair_shift.rel_err_b) > 0.05


def test_K_rank_perturbation_breaks_prediction() -> None:
    """Shift K_rank: 5 → 4 or 5 → 6 must break the prediction."""
    rig = rigidity_under_integer_shifts()
    for label in ("K_rank: 5 → 4", "K_rank: 5 → 6"):
        r = next(rg for rg in rig if label in rg.perturbation)
        assert abs(r.rel_err_c) > 0.05 or abs(r.rel_err_b) > 0.05, (
            f"{label}: m_c err {r.rel_err_c*100:+.2f}%, "
            f"m_b err {r.rel_err_b*100:+.2f}% — should be >5%"
        )


def test_n_R_perturbation_breaks_prediction() -> None:
    """Shift n_R: 18 → 12 must break the prediction badly."""
    rig = rigidity_under_integer_shifts()
    n_R_shift = next(r for r in rig if "n_R: 18 → 12" in r.perturbation)
    assert abs(n_R_shift.rel_err_c) > 0.10
    assert abs(n_R_shift.rel_err_b) > 0.10


# ---------------------------------------------------------------------------
# Alternative forms ALL fail (validates form selection)
# ---------------------------------------------------------------------------

def test_alternative_forms_mostly_fail() -> None:
    """Each rejected alternative form fails by >10%; the chosen form passes."""
    alts = alternative_forms_tried()
    chosen = next(a for a in alts if "[CHOSEN]" in a.label)
    rejected = [a for a in alts if "[CHOSEN]" not in a.label]

    # Chosen form is best
    assert chosen.max_abs_rel < 0.01

    # All rejected forms are off by >10% somewhere
    for a in rejected:
        assert a.max_abs_rel > 0.10, (
            f"Rejected form '{a.label}' was actually within 10% — "
            "form selection should have flagged this"
        )


def test_bare_torque_underpredicts_heavy_quarks() -> None:
    """T_q · Λ alone is off by >50% for c, b — motivates the lift."""
    alts = alternative_forms_tried()
    bare = next(a for a in alts if "Bare T_q" in a.label)
    assert bare.max_abs_rel > 0.5


# ---------------------------------------------------------------------------
# Wilson conversion to pole scheme
# ---------------------------------------------------------------------------

def test_wilson_pole_lifts_above_MSbar() -> None:
    """1-loop Wilson conversion gives m_pole > m_MSbar (positive shift)."""
    m_c_msbar = m_c_pole_substrate_GeV()
    m_c_wilson = m_c_pole_wilson_substrate_GeV(n_loops=1)
    m_b_msbar = m_b_pole_substrate_GeV()
    m_b_wilson = m_b_pole_wilson_substrate_GeV(n_loops=1)
    assert m_c_wilson > m_c_msbar
    assert m_b_wilson > m_b_msbar


def test_wilson_pole_charm_in_expected_range() -> None:
    """Substrate Wilson m_c pole ≈ 1.4-1.5 GeV (between MS-bar and full pole)."""
    m_c_pole = m_c_pole_wilson_substrate_GeV(n_loops=1)
    assert 1.40 < m_c_pole < 1.50


def test_wilson_pole_bottom_in_expected_range() -> None:
    """Substrate Wilson m_b pole ≈ 4.5-4.6 GeV."""
    m_b_pole = m_b_pole_wilson_substrate_GeV(n_loops=1)
    assert 4.50 < m_b_pole < 4.65


# ---------------------------------------------------------------------------
# Hadron-mass-test integration: substrate masses exposed at module level
# ---------------------------------------------------------------------------

def test_hadron_mass_test_exposes_substrate_masses() -> None:
    """hadron_mass_test re-exports M_C_MSBAR_SUBSTRATE_GEV and M_B_MSBAR_SUBSTRATE_GEV."""
    from src.stiff_medium.hadron_mass_test import (
        M_C_MSBAR_SUBSTRATE_GEV, M_B_MSBAR_SUBSTRATE_GEV,
    )
    assert math.isclose(M_C_MSBAR_SUBSTRATE_GEV, m_c_pole_substrate_GeV(), rel_tol=1e-12)
    assert math.isclose(M_B_MSBAR_SUBSTRATE_GEV, m_b_pole_substrate_GeV(), rel_tol=1e-12)


def test_hadron_mass_test_keeps_empirical_pole_for_cornell() -> None:
    """M_C_POLE_GEV, M_B_POLE_GEV remain empirical (Cat-C) for Cornell scheme.

    Cornell phenomenology uses kinetic/pole-strip mass which is HIGHER than
    MS-bar by an O(α_s) Wilson coefficient. The substrate's MS-bar prediction
    is correct but in a different scheme; using it directly degrades the
    Cornell J/ψ and Υ predictions. Honest decision: keep empirical Cornell
    pole-strip values (Cat-C) for the Cornell solver, expose the substrate
    MS-bar prediction (Cat-A*) separately.
    """
    from src.stiff_medium.hadron_mass_test import M_C_POLE_GEV, M_B_POLE_GEV
    # Empirical Cornell-scheme values
    assert math.isclose(M_C_POLE_GEV, 1.32, abs_tol=0.01)
    assert math.isclose(M_B_POLE_GEV, 4.50, abs_tol=0.01)


# ---------------------------------------------------------------------------
# Cornell J/ψ, Υ verification: predictions remain solid with empirical scheme
# ---------------------------------------------------------------------------

def test_cornell_jpsi_remains_under_one_percent() -> None:
    """J/ψ from Cornell with substrate α_s + empirical Cornell pole stays <1%."""
    from src.stiff_medium.hadron_mass_test import (
        _quarkonium_mass_MeV, M_C_POLE_GEV, ALPHA_S_C, SIGMA_GEV2,
    )
    jpsi_pred = _quarkonium_mass_MeV(M_C_POLE_GEV, ALPHA_S_C, SIGMA_GEV2)
    jpsi_pdg = 3096.9
    assert abs(jpsi_pred - jpsi_pdg) / jpsi_pdg < 0.01


def test_cornell_upsilon_remains_under_three_percent() -> None:
    """Υ from Cornell with substrate α_s + empirical Cornell pole stays <3%."""
    from src.stiff_medium.hadron_mass_test import (
        _quarkonium_mass_MeV, M_B_POLE_GEV, ALPHA_S_B, SIGMA_GEV2,
    )
    upsilon_pred = _quarkonium_mass_MeV(M_B_POLE_GEV, ALPHA_S_B, SIGMA_GEV2)
    upsilon_pdg = 9460.3
    assert abs(upsilon_pred - upsilon_pdg) / upsilon_pdg < 0.03


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def test_report_renders_without_error() -> None:
    """format_report produces a multi-line text report."""
    report = compute_full_report()
    text = format_report(report)
    assert "Substrate heavy-quark pole-mass" in text
    assert "K_pair² / K_rank = 4/5" in text
    assert "n_R / (K_pair · K_rank)" in text


def test_verdict_calls_out_Cat_A_star() -> None:
    """Verdict honestly calls the form Cat-A*, not Cat-A."""
    report = compute_full_report()
    assert "Cat-A*" in report.verdict
    assert "no independent physical derivation" in report.verdict
