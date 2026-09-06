"""Tests for HadronSpectrum (K_4 cell-stacking hadron mass calculator)."""

from __future__ import annotations

import math

import pytest

from src.stiff_medium.hadron_spectrum import (
    HadronSpectrum,
    QUARK_TORQUE,
    OCTET,
    DECUPLET,
    PS_NONET,
    V_NONET,
    PDG,
)


@pytest.fixture
def hs() -> HadronSpectrum:
    return HadronSpectrum()


# ---------------------------------------------------------------------------
# Required acceptance criteria from the task spec
# ---------------------------------------------------------------------------


def test_proton_within_1pct(hs: HadronSpectrum) -> None:
    pred = hs.baryon_mass("p")
    pdg = PDG["p"]
    assert abs(pred - pdg) / pdg < 0.01, f"p={pred:.2f} vs PDG {pdg:.2f}"


def test_delta_within_5pct(hs: HadronSpectrum) -> None:
    for name in ("Delta++", "Delta+", "Delta0", "Delta-"):
        pred = hs.baryon_mass(name)
        pdg = PDG[name]
        assert abs(pred - pdg) / pdg < 0.05, f"{name}={pred:.2f} vs {pdg:.2f}"


def test_sigma0_within_5pct(hs: HadronSpectrum) -> None:
    pred = hs.baryon_mass("Sigma0")
    pdg = PDG["Sigma0"]
    assert abs(pred - pdg) / pdg < 0.05, f"Sigma0={pred:.2f} vs {pdg:.2f}"


def test_gell_mann_okubo_holds(hs: HadronSpectrum) -> None:
    g = hs.gell_mann_okubo()
    # GMO should hold to within ~5% in this leading inventory model.
    assert abs(g["rel_residual"]) < 0.05, g


# ---------------------------------------------------------------------------
# Spectrum coverage / sanity
# ---------------------------------------------------------------------------


def test_octet_full_coverage(hs: HadronSpectrum) -> None:
    spec = hs.octet_spectrum()
    assert set(spec.keys()) == set(OCTET)
    assert len(spec) == 8
    for m in spec.values():
        assert 800.0 < m < 1700.0


def test_decuplet_full_coverage(hs: HadronSpectrum) -> None:
    spec = hs.decuplet_spectrum()
    assert set(spec.keys()) == set(DECUPLET)
    assert len(spec) == 10
    for m in spec.values():
        assert 1100.0 < m < 2000.0


def test_meson_nonet_runs(hs: HadronSpectrum) -> None:
    spec = hs.meson_nonet()
    # at least the four canonical light vectors and the two pseudoscalars present
    for n in ("pi", "K", "rho", "K*", "phi"):
        assert n in spec
        assert spec[n] > 0.0


# ---------------------------------------------------------------------------
# Mass orderings (model-independent qualitative checks)
# ---------------------------------------------------------------------------


def test_mass_ordering_octet(hs: HadronSpectrum) -> None:
    assert hs.baryon_mass("p") < hs.baryon_mass("Lambda")
    assert hs.baryon_mass("Lambda") < hs.baryon_mass("Sigma0")
    assert hs.baryon_mass("Sigma0") < hs.baryon_mass("Xi0")


def test_mass_ordering_decuplet(hs: HadronSpectrum) -> None:
    assert hs.baryon_mass("Delta+") < hs.baryon_mass("Sigma*0")
    assert hs.baryon_mass("Sigma*0") < hs.baryon_mass("Xi*0")
    assert hs.baryon_mass("Xi*0") < hs.baryon_mass("Omega-")


def test_decuplet_above_octet_for_same_quarks(hs: HadronSpectrum) -> None:
    assert hs.baryon_mass("Delta+") > hs.baryon_mass("p")
    assert hs.baryon_mass("Sigma*0") > hs.baryon_mass("Sigma0")
    assert hs.baryon_mass("Xi*0") > hs.baryon_mass("Xi0")


def test_pseudoscalars_lighter_than_vectors(hs: HadronSpectrum) -> None:
    # Same quark content: pseudoscalar < vector
    assert hs.meson_mass("pi") < hs.meson_mass("rho")
    assert hs.meson_mass("K") < hs.meson_mass("K*")


# ---------------------------------------------------------------------------
# Compare-to-PDG plumbing & report
# ---------------------------------------------------------------------------


def test_compare_to_pdg_returns_residuals(hs: HadronSpectrum) -> None:
    res = hs.compare_to_pdg()
    assert len(res) > 0
    for r in res:
        assert r.pdg > 0.0
        assert math.isfinite(r.pred)
        assert math.isfinite(r.rel_err)


def test_report_contains_key_lines(hs: HadronSpectrum) -> None:
    text = hs.report()
    assert "B3 hadron spectrum" in text
    assert "GMO check" in text
    assert "p" in text and "Delta+" in text and "Omega-" in text


# ---------------------------------------------------------------------------
# Constant invariants
# ---------------------------------------------------------------------------


def test_quark_torque_ordering() -> None:
    assert QUARK_TORQUE["u"] < QUARK_TORQUE["s"]
    assert QUARK_TORQUE["s"] < QUARK_TORQUE["c"]
    assert QUARK_TORQUE["c"] < QUARK_TORQUE["b"]


def test_isospin_split_small() -> None:
    rel = abs(QUARK_TORQUE["d"] - QUARK_TORQUE["u"]) / QUARK_TORQUE["u"]
    assert rel < 0.01  # small d-u split (< 1%)


def test_unknown_meson_raises(hs: HadronSpectrum) -> None:
    with pytest.raises(KeyError):
        hs.meson_mass("not_a_meson")


def test_unknown_baryon_raises(hs: HadronSpectrum) -> None:
    with pytest.raises(KeyError):
        hs.baryon_mass("not_a_baryon")
