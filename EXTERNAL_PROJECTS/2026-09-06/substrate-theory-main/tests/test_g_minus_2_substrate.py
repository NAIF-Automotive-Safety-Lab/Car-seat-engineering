"""Tests for src/stiff_medium/g_minus_2_substrate.py."""

from __future__ import annotations

import math

import pytest

from stiff_medium.g_minus_2_substrate import (
    A_MU_BMW_LATTICE,
    A_MU_FERMILAB_2023,
    A_MU_FERMILAB_2023_ERR,
    A_MU_SM_TI_2020,
    GMinus2,
)


@pytest.fixture
def gm2() -> GMinus2:
    return GMinus2()


def test_qed_anchor_is_aoyama_5loop(gm2: GMinus2) -> None:
    """QED component matches Aoyama et al. (Atoms 7, 28) within 1 × 10⁻¹¹."""
    assert abs(gm2.a_mu_QED() - 116_584_718.931) < 1.0


def test_hvp_branches_known_values(gm2: GMinus2) -> None:
    """Both HVP sources return their published anchors."""
    assert gm2.a_mu_hadronic_VP("TI") == pytest.approx(6_845.0, abs=1.0)
    assert gm2.a_mu_hadronic_VP("BMW") == pytest.approx(7_075.0, abs=1.0)
    with pytest.raises(ValueError):
        gm2.a_mu_hadronic_VP("nonsense")


def test_lbl_and_ew_within_published_range(gm2: GMinus2) -> None:
    assert 70.0 < gm2.a_mu_hadronic_LbL() < 120.0
    assert 150.0 < gm2.a_mu_EW() < 160.0


def test_substrate_contribution_is_small_and_positive(gm2: GMinus2) -> None:
    """B3 substrate δa_μ should sit in the few × 10⁻¹¹ band — large enough
    to matter, far too small to overshoot the SM."""
    a_sub = gm2.a_mu_substrate()
    # Reasonable-bound check: 0 < a_sub < 1000 × 10⁻¹¹.
    assert a_sub > 0.0
    assert a_sub < 1000.0
    # Order-of-magnitude sanity: should not exceed the EW component.
    assert a_sub < gm2.a_mu_EW()


def test_total_within_5sigma_of_fermilab_with_BMW_HVP(gm2: GMinus2) -> None:
    """Total = QED + HVP(BMW) + LbL + EW + B3 within 5σ of Fermilab 2023."""
    cmp = gm2.compare_to_fermilab(hvp_source="BMW")
    assert abs(cmp["n_sigma"]) < 5.0, cmp


def test_total_within_5sigma_of_fermilab_with_TI_HVP(gm2: GMinus2) -> None:
    """With TI'20 HVP the substrate contribution should help close the
    Fermilab gap; still pass the 5σ envelope."""
    cmp = gm2.compare_to_fermilab(hvp_source="TI")
    assert abs(cmp["n_sigma"]) < 5.0, cmp


def test_compare_to_BMW_lattice_returns_well_formed_dict(gm2: GMinus2) -> None:
    out = gm2.compare_to_BMW_lattice()
    for key in ("a_mu_B3", "a_mu_BMW", "diff", "combined_sigma", "n_sigma"):
        assert key in out
    assert out["a_mu_BMW"] == A_MU_BMW_LATTICE
    assert out["combined_sigma"] > 0
    # B3+TI vs BMW must be at least *closer* than TI alone vs BMW (substrate
    # adds positive contribution moving TI toward BMW).
    diff_pure_ti = A_MU_SM_TI_2020 - A_MU_BMW_LATTICE
    assert abs(out["diff"]) < abs(diff_pure_ti)


def test_total_is_sum_of_components(gm2: GMinus2) -> None:
    total = gm2.a_mu_total(hvp_source="BMW")
    expected = (
        gm2.a_mu_QED()
        + gm2.a_mu_hadronic_VP("BMW")
        + gm2.a_mu_hadronic_LbL()
        + gm2.a_mu_EW()
        + gm2.a_mu_substrate()
    )
    assert math.isclose(total, expected, rel_tol=1e-12)


def test_report_renders_and_mentions_key_anchors(gm2: GMinus2) -> None:
    text = gm2.report()
    assert isinstance(text, str)
    for needle in ("Fermilab", "BMW", "QED", "substrate", "HVP"):
        assert needle in text


def test_substrate_uses_topology_fixed_xi(gm2: GMinus2) -> None:
    """ξ_μ = 9/125 is hard-coded from braid-junction inventory."""
    assert gm2.xi_mu == pytest.approx(9.0 / 125.0)
