"""Tests for src/stiff_medium/quark_mass_spectrum.py"""

import math

import pytest

from stiff_medium.quark_mass_spectrum import (
    GENERATION_MAP,
    PDG_MASSES_MEV,
    QuarkMassSpectrum,
)


@pytest.fixture
def spec() -> QuarkMassSpectrum:
    return QuarkMassSpectrum()


# ------------------------------------------------------------------
#  Mass ordering
# ------------------------------------------------------------------

def test_mass_ordering_strict(spec: QuarkMassSpectrum) -> None:
    """m_t > m_b > m_c > m_s > m_d > m_u (strict)."""
    masses = [spec.mass(q) for q in ("u", "d", "s", "c", "b", "t")]
    for lo, hi in zip(masses, masses[1:]):
        assert hi > lo, f"ordering violated: {hi} not > {lo}"


def test_mass_ordering_explicit(spec: QuarkMassSpectrum) -> None:
    assert spec.mass("t") > spec.mass("b")
    assert spec.mass("b") > spec.mass("c")
    assert spec.mass("c") > spec.mass("s")
    assert spec.mass("s") > spec.mass("d")
    assert spec.mass("d") > spec.mass("u")


# ------------------------------------------------------------------
#  PDG agreement (acknowledge: quark sector partial, allow 50%)
# ------------------------------------------------------------------

def test_masses_within_50pct_of_pdg(spec: QuarkMassSpectrum) -> None:
    for q, pdg in PDG_MASSES_MEV.items():
        m = spec.mass(q)
        rel = abs(m - pdg) / pdg
        assert rel < 0.50, f"{q}: substrate {m} vs PDG {pdg}, off by {rel:.1%}"


def test_ratios_within_50pct_of_pdg(spec: QuarkMassSpectrum) -> None:
    for label, d in spec.mass_ratios().items():
        rel = abs(d["err_pct"]) / 100.0
        assert rel < 0.50, f"ratio {label} off by {rel:.1%}"


# ------------------------------------------------------------------
#  API surface
# ------------------------------------------------------------------

def test_compare_to_pdg_returns_six(spec: QuarkMassSpectrum) -> None:
    comps = spec.compare_to_PDG()
    assert len(comps) == 6
    assert {c.quark for c in comps} == set("udscbt")


def test_running_mass_decreases_with_scale(spec: QuarkMassSpectrum) -> None:
    """For light quarks, MSbar mass DECREASES with rising mu (alpha_s shrinks)."""
    m_low = spec.running_mass("s", 2.0)
    m_high = spec.running_mass("s", 100.0)
    assert m_high < m_low


def test_running_mass_at_reference_matches(spec: QuarkMassSpectrum) -> None:
    # at the reference scale running_mass should equal mass()
    m_ref = spec.running_mass("s", 2.0)
    assert math.isclose(m_ref, spec.mass("s"), rel_tol=1e-6)


def test_pole_above_msbar_for_heavy(spec: QuarkMassSpectrum) -> None:
    pv = spec.pole_vs_msbar()
    for q in ("c", "b", "t"):
        assert pv[q]["pole_gev"] > pv[q]["msbar_gev"]
        assert pv[q]["shift_pct"] > 0


def test_constituent_above_current_for_light(spec: QuarkMassSpectrum) -> None:
    cv = spec.current_vs_constituent()
    for q in ("u", "d", "s"):
        assert cv[q]["constituent_mev"] > cv[q]["current_mev"]
    # u, d constituent close to canonical ~300 MeV
    assert 250.0 < cv["u"]["constituent_mev"] < 350.0
    assert 250.0 < cv["d"]["constituent_mev"] < 350.0
    # s constituent close to canonical ~500 MeV
    assert 450.0 < cv["s"]["constituent_mev"] < 550.0


def test_koide_returns_both_towers(spec: QuarkMassSpectrum) -> None:
    k = spec.koide_quark_relation()
    assert "up_type" in k and "down_type" in k
    for d in k.values():
        assert 0.0 < d["Q"] < 1.0  # Q always in (0, 1) for positive masses


def test_report_contains_key_sections(spec: QuarkMassSpectrum) -> None:
    r = spec.report()
    for token in ("QUARK MASS SPECTRUM", "ratios", "MSbar", "Constituent", "Koide"):
        assert token.lower() in r.lower()


def test_unknown_quark_raises(spec: QuarkMassSpectrum) -> None:
    with pytest.raises(KeyError):
        spec.mass("g")
