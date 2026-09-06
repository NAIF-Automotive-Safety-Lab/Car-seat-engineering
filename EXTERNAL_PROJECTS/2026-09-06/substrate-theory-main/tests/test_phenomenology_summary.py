"""Tests for PhenomenologySummary aggregator."""

import os
import tempfile

import pytest

from src.stiff_medium.phenomenology_summary import PhenomenologySummary


@pytest.fixture
def summary():
    return PhenomenologySummary()


# --------------------------------------------------------------------
# Each table generates and has the expected keys / shape
# --------------------------------------------------------------------
def test_table_1_basic_constants(summary):
    t = summary.paper_table_1_basic_constants()
    for k in ("alpha_inv", "m_e_MeV", "m_p_MeV", "c_m_per_s", "hbar_eVs"):
        assert k in t, f"missing {k}"
        row = t[k]
        assert "prediction" in row and "observed" in row and "source" in row


def test_table_2_lepton_masses(summary):
    t = summary.paper_table_2_lepton_masses()
    for name in ("electron", "muon", "tau"):
        assert name in t
        assert t[name]["prediction_MeV"] > 0
        assert t[name]["rel_err"] is not None
        assert t[name]["rel_err"] < 0.01  # all within 1%


def test_table_3_quark_masses(summary):
    t = summary.paper_table_3_quark_masses()
    assert "up_MeV" in t and "top_GeV" in t
    assert t["top_GeV"]["prediction"] > 100  # ~172 GeV


def test_table_4_mixing_matrices(summary):
    t = summary.paper_table_4_mixing_matrices()
    assert "CKM" in t and "PMNS" in t
    assert "Vud" in t["CKM"]
    assert "theta12_deg" in t["PMNS"]
    # Cabibbo 1/sqrt(20) prediction
    cab = t["CKM"]["Cabibbo_lambda_pred"]
    assert abs(cab["prediction"] - 0.2236) < 1e-3
    # within ~2% of Vus
    assert cab["rel_err"] < 0.02


def test_table_5_cosmology(summary):
    t = summary.paper_table_5_cosmology()
    for k in ("H_0", "Sigma_mnu", "Omega_m", "Omega_L"):
        assert k in t
    assert 60 < t["H_0"]["prediction"] < 80
    assert 0 < t["Omega_m"]["prediction"] < 1


def test_table_6_dark_sector(summary):
    t = summary.paper_table_6_dark_sector()
    assert "m_DM_GeV" in t
    assert "sigma_SI_cm2" in t
    assert t["m_DM_GeV"]["prediction"] > 0


def test_table_7_falsifiers(summary):
    t = summary.paper_table_7_falsifiers()
    assert isinstance(t, list)
    assert len(t) >= 5
    for row in t:
        assert {"test", "experiment", "status", "result"} <= set(row.keys())


# --------------------------------------------------------------------
# Report
# --------------------------------------------------------------------
def test_report_contains_all_sections(summary):
    r = summary.report()
    assert "Table 1" in r
    assert "Table 2" in r
    assert "Table 3" in r
    assert "Table 4a" in r
    assert "Table 4b" in r
    assert "Table 5" in r
    assert "Table 6" in r
    assert "Table 7" in r
    assert "B3 Substrate" in r


# --------------------------------------------------------------------
# LaTeX export
# --------------------------------------------------------------------
def test_latex_export_creates_file_and_valid_skeleton(summary, tmp_path):
    fp = tmp_path / "out.tex"
    body = summary.latex_export(str(fp))
    assert fp.exists()
    text = fp.read_text(encoding="utf-8")
    assert text == body
    # LaTeX skeleton
    assert r"\documentclass{article}" in text
    assert r"\begin{document}" in text
    assert r"\end{document}" in text
    # All seven tables represented
    for cap in ("Table 1", "Table 2", "Table 3",
                "Table 4a", "Table 4b", "Table 5",
                "Table 6", "Table 7"):
        assert cap in text, f"missing caption {cap}"
    # Booktabs lines
    assert r"\toprule" in text and r"\bottomrule" in text
    # Balanced begin/end tabular
    assert text.count(r"\begin{tabular}") == text.count(r"\end{tabular}")
    assert text.count(r"\begin{table}") == text.count(r"\end{table}")


def test_latex_escapes_special_chars():
    from src.stiff_medium.phenomenology_summary import _tex_escape
    assert _tex_escape("a_b") == r"a\_b"
    assert _tex_escape("50%") == r"50\%"
    assert _tex_escape("x&y") == r"x\&y"
