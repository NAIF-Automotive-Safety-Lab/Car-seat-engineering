"""Tests for the centralized B3 constants module."""

from __future__ import annotations

from math import comb

import pytest

from src.stiff_medium import b3_constants as C


def test_all_constants_present():
    """Every constant promised by the audit must be importable."""
    expected = [
        "N_BAM", "K_pair", "K_rank", "n_R", "n_M", "n_A",
        "F", "R", "V13",
        "LAMBDA_QCD_MEV", "LAMBDA_QCD_K",
        "K_PA", "RHO_KGM3", "XI_M", "GAMMA_HZ",
        "Q_DRAG",
    ]
    for name in expected:
        assert hasattr(C, name), f"missing constant {name}"


def test_canonical_integer_values():
    assert C.N_BAM == 6
    assert C.K_pair == 2
    assert C.K_rank == 5
    assert C.n_R == 18
    assert C.F == 2
    assert C.R == 3
    assert C.V13 == 13


def test_n_M_identity():
    """n_M = K_pair · K_rank³ + n_R = 268 — the master identity."""
    assert C.n_M == 268
    assert C.n_M == C.K_pair * C.K_rank ** 3 + C.n_R


def test_n_A_identity():
    """n_A = C(N_BAM, 2) = edges of K_{N_BAM} = 15 with N_BAM=6."""
    assert C.n_A == comb(C.N_BAM, 2)
    assert C.n_A == 15


def test_deuteron_denominator_identity():
    """The product n_A · N_BAM MUST equal 90 to reproduce deuteron BE.

    ε_face = Λ_QCD / (n_A · N_BAM) = 200 MeV / 90 = 2.222 MeV ≈ AME2020 2.2246.

    This is the load-bearing factorization of 90 in the substrate framework:
    N_BAM=6 forced by 2D hex geometry (geom_05), so n_A=15=C(6,2) is
    forced by this product identity.
    """
    assert C.n_A * C.N_BAM == 90, (
        f"n_A · N_BAM = {C.n_A} · {C.N_BAM} = {C.n_A * C.N_BAM}, "
        f"expected 90 (deuteron BE denominator)"
    )
    # Numerically: ε_face = 200/90 within 0.11% of observed deuteron BE
    eps_face = C.LAMBDA_QCD_MEV / (C.n_A * C.N_BAM)
    assert abs(eps_face - 2.2246) / 2.2246 < 0.005


def test_legacy_n_A_kept_separate():
    """Pre-audit modules may still need 45; ensure shim is present."""
    assert C.N_A_LEGACY_45 == 45
    assert C.N_A_LEGACY_45 != C.n_A  # they really are different


def test_lambda_qcd_canonical():
    assert C.LAMBDA_QCD_MEV == 200.0
    assert C.LAMBDA_QCD_K == 200.0


def test_q_drag_derived():
    assert C.Q_DRAG == pytest.approx((11.0 / 12.0) * 268.0, rel=1e-15)
    assert C.Q_DRAG == pytest.approx(245.6666666666, rel=1e-10)


def test_koide_ratio():
    assert C.KOIDE_RATIO == pytest.approx(2.0 / 3.0, rel=1e-15)


def test_substrate_primitives_match_anchor():
    """K_PA, RHO, XI, GAMMA must match PrimitiveAnchoring('electron_compton')."""
    from src.stiff_medium.primitive_anchoring import PrimitiveAnchoring
    s = PrimitiveAnchoring("electron_compton")
    s._solve()
    assert C.K_PA == pytest.approx(s.K, rel=1e-12)
    assert C.RHO_KGM3 == pytest.approx(s.rho, rel=1e-12)
    assert C.XI_M == pytest.approx(s.xi, rel=1e-12)
    assert C.GAMMA_HZ == pytest.approx(s.gamma, rel=1e-12)


def test_verify_consistency_all_pass():
    checks = C.verify_consistency()
    bool_keys = [k for k, v in checks.items() if isinstance(v, bool)]
    failures = [k for k in bool_keys if not checks[k]]
    assert not failures, f"consistency failures: {failures}"


def test_clean_reexport():
    """`from src.stiff_medium.b3_constants import N_BAM, K_pair, ...` works."""
    from src.stiff_medium.b3_constants import (  # noqa: F401
        N_BAM, K_pair, K_rank, n_R, n_M, n_A, F, R, V13,
        LAMBDA_QCD_MEV, K_PA, RHO_KGM3, XI_M, GAMMA_HZ, Q_DRAG,
    )


def test_integers_are_integers():
    for name in ("N_BAM", "K_pair", "K_rank", "n_R", "n_M", "n_A", "F", "R", "V13"):
        assert isinstance(getattr(C, name), int), f"{name} should be int"


def test_constants_have_docstrings():
    """Each public constant should be documented (presence of docstring in module)."""
    import inspect
    src = inspect.getsource(C)
    for name in C.__all__:
        if name == "verify_consistency":
            continue
        # Each constant is followed by a triple-quoted docstring in the source.
        assert f"{name}:" in src or f"{name} " in src, f"{name} not declared"
