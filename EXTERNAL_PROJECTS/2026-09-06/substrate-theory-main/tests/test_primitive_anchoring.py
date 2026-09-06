"""Tests for src/stiff_medium/primitive_anchoring.py."""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.primitive_anchoring import (
    ALPHA_B3,
    ALPHA_OBS,
    ANCHOR_CHOICES,
    HIERARCHY_B3,
    MP_OVER_ME,
    AnchorResult,
    PrimitiveAnchoring,
    best_anchor,
    compare_anchors,
    report_all_anchors,
)


# ---------------------------------------------------------------------------
# Electron-Compton anchor: the canonical B3 substrate
# ---------------------------------------------------------------------------
def test_electron_compton_xi_value() -> None:
    pa = PrimitiveAnchoring("electron_compton")
    assert pa.xi == pytest.approx(3.8616e-13, rel=1e-3)


def test_electron_compton_K_order() -> None:
    pa = PrimitiveAnchoring("electron_compton")
    assert 1e23 < pa.K < 1e25  # ≈ 1.42e24 Pa


def test_electron_compton_rho_order() -> None:
    pa = PrimitiveAnchoring("electron_compton")
    assert 1e6 < pa.rho < 1e8  # ≈ 1.58e7 kg/m^3


def test_electron_compton_gamma_positive() -> None:
    pa = PrimitiveAnchoring("electron_compton")
    assert pa.gamma > 0
    assert pa.gamma == pytest.approx(7.76e20, rel=1e-2)


def test_alpha_reproduces_CODATA() -> None:
    """ℏ c = K ξ⁴ identity must reconstruct α to machine precision."""
    pa = PrimitiveAnchoring("electron_compton")
    assert pa.derived_alpha() == pytest.approx(ALPHA_OBS, rel=1e-9)


def test_alpha_close_to_B3_formula() -> None:
    """B3 formula 11/(48π³)·exp(-3π/737) ≈ α to ~10⁻⁵."""
    res = PrimitiveAnchoring("electron_compton").compute_cross_checks()
    assert abs(res.residuals["alpha_vs_B3_formula"]) < 1e-4


def test_kink_equals_8me_at_electron_compton() -> None:
    pa = PrimitiveAnchoring("electron_compton")
    assert pa.derived_kink_over_electron() == pytest.approx(8.0, rel=1e-12)


def test_kink_mass_4_MeV() -> None:
    pa = PrimitiveAnchoring("electron_compton")
    # 8 m_e c² = 8 × 0.511 MeV ≈ 4.088 MeV
    assert pa.derived_kink_mass_MeV() == pytest.approx(4.088, rel=1e-3)


def test_mp_over_me_consistency() -> None:
    pa = PrimitiveAnchoring("electron_compton")
    assert pa.derived_mp_over_me_via_kink() == pytest.approx(MP_OVER_ME, rel=1e-9)


# ---------------------------------------------------------------------------
# Alternative anchors
# ---------------------------------------------------------------------------
def test_proton_compton_anchor() -> None:
    pa = PrimitiveAnchoring("proton_compton")
    assert pa.xi == pytest.approx(2.103e-16, rel=1e-3)
    assert pa.K > 1e30
    assert pa.derived_alpha() == pytest.approx(ALPHA_OBS, rel=1e-9)


def test_planck_length_anchor() -> None:
    pa = PrimitiveAnchoring("planck_length")
    assert pa.xi == pytest.approx(1.6163e-35, rel=1e-3)
    # hierarchy ξ/ℓ_P = 1 exactly for this anchor
    assert pa.derived_hierarchy() == pytest.approx(1.0, rel=1e-12)


def test_invalid_anchor_raises() -> None:
    with pytest.raises(ValueError):
        PrimitiveAnchoring("not_an_anchor")


# ---------------------------------------------------------------------------
# Multi-anchor scan / best selection
# ---------------------------------------------------------------------------
def test_compare_anchors_returns_all() -> None:
    results = compare_anchors()
    assert set(results.keys()) == set(ANCHOR_CHOICES)
    for r in results.values():
        assert isinstance(r, AnchorResult)


def test_best_anchor_is_electron_compton() -> None:
    """Electron Compton anchor minimizes kink+mp/me consistency residuals."""
    assert best_anchor() == "electron_compton"


def test_alpha_anchor_invariant() -> None:
    """α from ℏc identity is anchor-independent (universality check)."""
    vals = [PrimitiveAnchoring(a).derived_alpha() for a in ANCHOR_CHOICES]
    for v in vals:
        assert v == pytest.approx(ALPHA_OBS, rel=1e-9)


def test_report_cross_checks_string(capsys: pytest.CaptureFixture) -> None:
    pa = PrimitiveAnchoring("electron_compton")
    txt = pa.report_cross_checks()
    captured = capsys.readouterr()
    assert "Anchor: electron_compton" in txt
    assert "ξ" in txt or "xi" in captured.out.lower() or "K" in txt
    assert "Residuals" in txt


def test_report_all_anchors_runs(capsys: pytest.CaptureFixture) -> None:
    txt = report_all_anchors()
    assert "Best self-consistent anchor" in txt
    for a in ANCHOR_CHOICES:
        assert a in txt


# ---------------------------------------------------------------------------
# Drag-as-mass-generator anchoring (new)
# ---------------------------------------------------------------------------
def test_solve_with_drag_mass_anchor_electron() -> None:
    """γ from drag-mass identity must give ℏγ = m_e c² for electron."""
    import scipy.constants as sc

    pa = PrimitiveAnchoring("electron_compton")
    gamma_drag = pa.solve_with_drag_mass_anchor("electron")
    expected = sc.m_e * sc.c ** 2 / sc.hbar
    assert gamma_drag == pytest.approx(expected, rel=1e-12)


def test_drag_to_mass_ratio_unity_for_electron_compton() -> None:
    """For electron_compton anchor, γ = c/ξ = m_e c²/ℏ so ratio = 1 exactly."""
    pa = PrimitiveAnchoring("electron_compton")
    assert pa.drag_to_mass_ratio("electron") == pytest.approx(1.0, rel=1e-12)


def test_electron_drag_anchor_yields_electron_mass() -> None:
    """electron_drag anchor reconstructs m_e to machine precision."""
    import scipy.constants as sc

    pa = PrimitiveAnchoring("electron_drag")
    m_kg = pa.derived_drag_mass_kg()
    assert m_kg == pytest.approx(sc.m_e, rel=1e-12)
    # cross_check_drag_alpha returns the same residual as the standard check
    drag = pa.cross_check_drag_alpha()
    assert abs(drag["residual_vs_B3"]) < 1e-4


def test_B3_constants_sane() -> None:
    """Sanity: B3 α formula and hierarchy magnitude."""
    assert 7.29e-3 < ALPHA_B3 < 7.30e-3
    assert HIERARCHY_B3 > 1e16
