"""Tests for stiff_medium.ckm_complete."""

import math
import pytest

from stiff_medium.ckm_complete import CKMComplete, PDG, WolfensteinParams


@pytest.fixture(scope="module")
def ckm():
    return CKMComplete()


def test_wolfenstein_returns_dataclass(ckm):
    w = ckm.wolfenstein()
    assert isinstance(w, WolfensteinParams)
    for attr in ("lam", "A", "rho_bar", "eta_bar"):
        assert hasattr(w, attr)


def test_lambda_matches_pdg_within_1pct(ckm):
    """Test (a): lambda = 1/sqrt(20) = 0.2236, within 1% of PDG 0.225."""
    lam = ckm.wolfenstein().lam
    assert abs(lam - 0.2236) < 1e-3
    frac_err = abs(lam - PDG["lambda"]) / PDG["lambda"]
    assert frac_err < 0.01, f"|lambda - PDG|/PDG = {frac_err:.4f} > 1%"


def test_A_substrate_ratio(ckm):
    """A = K_pair/K_rank = 16/19 ~ 0.842."""
    A = ckm.wolfenstein().A
    assert abs(A - 16.0 / 19.0) < 1e-12
    # within ~3% of PDG
    assert abs(A - PDG["A"]) / PDG["A"] < 0.05


def test_unitarity_to_1e12(ckm):
    """Test (b): V V^dag = I to better than 1e-12."""
    res = ckm.unitarity_residual()
    assert res < 1e-12, f"unitarity residual {res:.3e} exceeds 1e-12"


def test_jarlskog_order_of_magnitude(ckm):
    """Test (c): J ~ 3e-5, within order of magnitude."""
    J = ckm.jarlskog_invariant()
    assert J > 0
    # Order-of-magnitude window: PDG / 10  <  J  <  PDG * 10
    assert PDG["J"] / 10 < J < PDG["J"] * 10
    # Stronger: within factor 2
    assert 0.5 * PDG["J"] < J < 2.0 * PDG["J"]


def test_matrix_explicit_shape_and_diagonal(ckm):
    M = ckm.matrix_explicit()
    assert len(M) == 3 and all(len(r) == 3 for r in M)
    # Diagonal should be near 1
    for i in range(3):
        assert M[i][i] > 0.97
    # |Vub| should be small (~ few x 1e-3)
    assert 1e-3 < M[0][2] < 1e-2


def test_unitarity_triangle_returns_three_angles(ckm):
    a, b, g = ckm.unitarity_triangle()
    for ang in (a, b, g):
        assert 0 < ang < math.pi
    # Angles should sum to ~pi (closure of the triangle)
    assert abs((a + b + g) - math.pi) < 1e-9


def test_B_meson_mixing_keys_and_positivity(ckm):
    bm = ckm.B_meson_mixing()
    for k in ("DeltaM_d", "DeltaM_s", "ratio_s_d"):
        assert k in bm and bm[k] > 0
    # Hierarchy: dM_s >> dM_d
    assert bm["DeltaM_s"] > 10 * bm["DeltaM_d"]


def test_compare_to_PDG_returns_tuples(ckm):
    cmp = ckm.compare_to_PDG()
    for key, val in cmp.items():
        assert len(val) == 3, f"{key} should be (pred, pdg, frac_err)"
        pred, pdg, err = val
        assert isinstance(pred, float)
        assert isinstance(pdg, float)
        assert isinstance(err, float)


def test_report_is_nonempty_string(ckm):
    r = ckm.report()
    assert isinstance(r, str)
    assert "Wolfenstein" in r
    assert "Jarlskog" in r
    assert "lambda" in r


def test_Vus_close_to_lambda(ckm):
    """|Vus| should equal lambda to good precision."""
    M = ckm.matrix_explicit()
    Vus = M[0][1]
    assert abs(Vus - ckm.lam) / ckm.lam < 1e-3
