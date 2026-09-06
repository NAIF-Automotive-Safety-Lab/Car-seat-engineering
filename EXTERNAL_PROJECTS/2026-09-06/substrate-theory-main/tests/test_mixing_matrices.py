"""Tests for CKM/PMNS mixing matrix calculator."""

import numpy as np
import pytest

from stiff_medium.mixing_matrices import MixingMatrices, PDG


@pytest.fixture
def mm():
    return MixingMatrices()


def test_ckm_unitary(mm):
    V = mm.ckm_matrix()
    assert mm.unitarity_check(V, tol=1e-12)
    assert mm.unitarity_residual(V) < 1e-12


def test_pmns_unitary(mm):
    U = mm.pmns_matrix()
    assert mm.unitarity_check(U, tol=1e-12)
    assert mm.unitarity_residual(U) < 1e-12


def test_cabibbo_within_2pct(mm):
    pred = mm.lam
    obs = PDG["sin_theta_C"]
    assert abs(pred - obs) / obs < 0.02


def test_pmns_angles_within_5pct(mm):
    for key, sub in [("sin2_theta12_pmns", mm.sin2_theta12_pmns),
                     ("sin2_theta13_pmns", mm.sin2_theta13_pmns),
                     ("sin2_theta23_pmns", mm.sin2_theta23_pmns)]:
        obs = PDG[key]
        assert abs(sub - obs) / obs < 0.05, f"{key}: pred={sub}, obs={obs}"


def test_jarlskog_ckm_order_of_magnitude(mm):
    J = mm.jarlskog_invariant_ckm()
    assert 1e-6 < abs(J) < 1e-4
    # within order of magnitude of measured 3.18e-5
    ratio = abs(J) / PDG["J_CKM"]
    assert 0.1 < ratio < 10.0


def test_jarlskog_pmns_order_of_magnitude(mm):
    J = mm.jarlskog_invariant_pmns()
    assert abs(J) > 1e-3
    ratio = abs(J) / PDG["J_PMNS"]
    assert 0.1 < ratio < 10.0


def test_pmns_jarlskog_much_larger_than_ckm(mm):
    Jc = abs(mm.jarlskog_invariant_ckm())
    Jp = abs(mm.jarlskog_invariant_pmns())
    assert Jp > 100.0 * Jc


def test_ckm_magnitudes_match_wolfenstein(mm):
    V = mm.ckm_matrix()
    # |V_us| ≈ λ
    assert abs(abs(V[0, 1]) - mm.lam) < 1e-3
    # |V_cb| ≈ A λ²
    assert abs(abs(V[1, 2]) - mm.A * mm.lam ** 2) < 1e-3
    # |V_ub| in correct ballpark
    assert 1e-3 < abs(V[0, 2]) < 1e-2


def test_compare_to_pdg_runs(mm, capsys):
    rep = mm.compare_to_pdg()
    out = capsys.readouterr().out
    assert "sin_theta_C" in out
    assert "J_PMNS" in out
    assert isinstance(rep, dict)


def test_params_dataclass(mm):
    p = mm.params()
    assert p.lam == mm.lam
    assert p.delta_CP_pmns == mm.delta_CP_pmns
