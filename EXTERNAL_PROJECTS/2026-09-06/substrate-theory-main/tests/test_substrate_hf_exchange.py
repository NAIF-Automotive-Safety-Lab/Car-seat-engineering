"""Tests for substrate_hf_exchange — substrate-derived HF exchange kernel.

Verifies that:
  * The three exchange coefficients are pure-integer ratios of K_pair=2
    and K_rank=5 (no per-element fitting).
  * The configuration helpers (k_p, k_s, doubly-occupied, has_inner_core)
    match Aufbau ground-state expectations for H..Ar.
  * The half-shell, closed-s, and s/p factors apply ONLY to the elements
    they should (no spurious corrections).
  * Applying the substrate exchange to the Roothaan-HF Koopmans IE drops
    the H..Ar mean error from 6.4% to 2.78%.
  * The predicted half-shell anomaly N→O has the correct sign.
"""

from __future__ import annotations

import pytest

from src.stiff_medium.b3_constants import K_pair, K_rank
from src.stiff_medium.substrate_hf_exchange import (
    ANGULAR_BUDGET,
    CLOSED_S_PAIR_FACTOR,
    P_SHELL_ANGULAR,
    SHEETS_PAULI,
    ExchangeAudit,
    apply_substrate_exchange,
    audit_exchange,
    closed_s_pair_factor,
    doubly_occupied_p,
    half_shell_anomaly_eV,
    half_shell_factor,
    has_inner_core,
    k_p_least_bound,
    k_s_least_bound,
    predict_substrate_HF_exchange,
    s_p_separation_factor,
    same_spin_pairs_p,
    substrate_exchange_factor,
)


# --------------------------------------------------------------------------- #
# Substrate constants are exact integer ratios from K_pair / K_rank           #
# --------------------------------------------------------------------------- #

def test_constants_are_integer_ratios():
    """SHEETS_PAULI, ANGULAR_BUDGET, P_SHELL_ANGULAR are pure integers from
    K_pair=2 and K_rank=5; CLOSED_S_PAIR_FACTOR = 1/10 exactly."""
    assert SHEETS_PAULI == K_pair == 2
    assert ANGULAR_BUDGET == K_rank == 5
    assert P_SHELL_ANGULAR == K_pair * K_rank == 10
    assert CLOSED_S_PAIR_FACTOR == 1.0 / 10.0
    assert CLOSED_S_PAIR_FACTOR == 1.0 / (K_pair * K_rank)


# --------------------------------------------------------------------------- #
# Configuration helpers (Aufbau ground state for H..Ar)                       #
# --------------------------------------------------------------------------- #

def test_k_p_least_bound_for_p_block():
    """Boron (2p¹) through Ne (2p⁶), and Al (3p¹) through Ar (3p⁶)."""
    expected_p = {
        5: 1, 6: 2, 7: 3, 8: 4, 9: 5, 10: 6,
        13: 1, 14: 2, 15: 3, 16: 4, 17: 5, 18: 6,
    }
    for Z, k in expected_p.items():
        assert k_p_least_bound(Z) == k, f"Z={Z}: expected k_p={k}"


def test_k_p_zero_for_s_targets():
    """k_p = 0 when the least-bound subshell is an s-shell."""
    s_targets = [1, 2, 3, 4, 11, 12]
    for Z in s_targets:
        assert k_p_least_bound(Z) == 0


def test_k_s_least_bound_for_s_block():
    """k_s = 1 for H, Li, Na (one s electron); k_s = 2 for He, Be, Mg."""
    assert k_s_least_bound(1) == 1   # H
    assert k_s_least_bound(2) == 2   # He
    assert k_s_least_bound(3) == 1   # Li
    assert k_s_least_bound(4) == 2   # Be
    assert k_s_least_bound(11) == 1  # Na
    assert k_s_least_bound(12) == 2  # Mg


def test_k_s_zero_for_p_targets():
    p_targets = [5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18]
    for Z in p_targets:
        assert k_s_least_bound(Z) == 0


def test_doubly_occupied_p():
    """Hund's rule: d = max(0, k_p - 3)."""
    expected = {
        5: 0, 6: 0, 7: 0,        # B, C, N — all unpaired
        8: 1, 9: 2, 10: 3,       # O, F, Ne — 1, 2, 3 doubly-occupied
        13: 0, 14: 0, 15: 0,     # Al, Si, P — all unpaired
        16: 1, 17: 2, 18: 3,     # S, Cl, Ar — 1, 2, 3 doubly-occupied
    }
    for Z, d in expected.items():
        assert doubly_occupied_p(Z) == d, f"Z={Z}: expected d={d}"


def test_same_spin_pairs_p():
    """Same-spin pair count for the p-shell under Hund's rule.

    p^k → (n_α, n_β) with n_α=min(k,3), n_β=max(0,k-3).
    ssp = C(n_α,2) + C(n_β,2).
    """
    expected = {
        5: 0,  6: 1,  7: 3,   # p^1, p^2, p^3
        8: 3,  9: 4,  10: 6,  # p^4, p^5, p^6
    }
    for Z, ssp in expected.items():
        assert same_spin_pairs_p(Z) == ssp, f"Z={Z}: expected ssp={ssp}"


def test_has_inner_core():
    """H, He: no inner core (n_t=1).  Li onwards: yes."""
    assert not has_inner_core(1)
    assert not has_inner_core(2)
    for Z in range(3, 19):
        assert has_inner_core(Z), f"Z={Z}: expected inner core"


# --------------------------------------------------------------------------- #
# Half-shell factor — only applies to broken half-shells                      #
# --------------------------------------------------------------------------- #

def test_half_shell_factor_one_for_intact_or_no_p_shell():
    """k_p ≤ 3 (B, C, N, Al, Si, P) and non-p targets all give factor=1."""
    intact_or_nonp = [1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15]
    for Z in intact_or_nonp:
        assert half_shell_factor(Z) == 1.0, (
            f"Z={Z}: expected half-shell factor = 1 for intact/non-p"
        )


def test_half_shell_factor_correct_for_broken():
    """For p^4/p^5/p^6 the formula is f = 1 - K_pair/k_p² = 1 - 2/k_p²."""
    expected = {
        8:  1.0 - 2.0 / 16.0,   # p^4: 1 - 2/16 = 7/8 = 0.875
        9:  1.0 - 2.0 / 25.0,   # p^5: 1 - 2/25 = 23/25 = 0.92
        10: 1.0 - 2.0 / 36.0,   # p^6: 1 - 2/36 = 17/18 ≈ 0.944
        16: 1.0 - 2.0 / 16.0,   # S analogous to O
        17: 1.0 - 2.0 / 25.0,   # Cl analogous to F
        18: 1.0 - 2.0 / 36.0,   # Ar analogous to Ne
    }
    for Z, f_expected in expected.items():
        assert half_shell_factor(Z) == pytest.approx(f_expected, abs=1e-12)


def test_half_shell_factor_below_one_for_broken_shells():
    """Half-shell factor must be < 1 for k_p > 3 (HF over-predicts there)."""
    for Z in [8, 9, 10, 16, 17, 18]:
        assert half_shell_factor(Z) < 1.0


# --------------------------------------------------------------------------- #
# Closed-s-pair factor                                                        #
# --------------------------------------------------------------------------- #

def test_closed_s_pair_factor_only_for_be_mg():
    """f = 11/10 for Be, Mg; f = 1 for everything else (including He, Li, Na)."""
    expected = {
        1: 1.0,                # H: no s-pair
        2: 1.0,                # He: no inner core
        3: 1.0,                # Li: s¹
        4: 1.0 + 1.0 / 10.0,   # Be: closed s² + inner core (1s²)
        11: 1.0,               # Na: s¹
        12: 1.0 + 1.0 / 10.0,  # Mg: closed s² + inner core (Ne)
    }
    for Z, f_expected in expected.items():
        assert closed_s_pair_factor(Z) == pytest.approx(f_expected, abs=1e-12)


def test_closed_s_pair_factor_one_for_p_targets():
    """Any p-target has factor = 1 from this term."""
    p_targets = [5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18]
    for Z in p_targets:
        assert closed_s_pair_factor(Z) == 1.0


# --------------------------------------------------------------------------- #
# s/p separation (placeholder — currently 1.0; reserved for future research)  #
# --------------------------------------------------------------------------- #

def test_s_p_separation_factor_is_unity():
    """Currently the s/p separation is captured by the K_rank screening
    in atom_substrate, so the dedicated factor here is 1.0."""
    for Z in range(1, 19):
        assert s_p_separation_factor(Z) == 1.0


# --------------------------------------------------------------------------- #
# Total exchange factor                                                       #
# --------------------------------------------------------------------------- #

def test_total_exchange_factor_one_for_simple_cases():
    """H, He, Li, Na, B, C, N, Al, Si, P all have total factor = 1
    (no half-shell anomaly + no closed-s with inner core)."""
    no_correction = [1, 2, 3, 5, 6, 7, 11, 13, 14, 15]
    for Z in no_correction:
        assert substrate_exchange_factor(Z) == 1.0, (
            f"Z={Z}: expected total exchange factor = 1"
        )


def test_total_exchange_factor_below_one_for_broken_p():
    for Z in [8, 9, 10, 16, 17, 18]:
        assert substrate_exchange_factor(Z) < 1.0


def test_total_exchange_factor_above_one_for_be_mg():
    for Z in [4, 12]:
        assert substrate_exchange_factor(Z) > 1.0
        assert substrate_exchange_factor(Z) == pytest.approx(11.0 / 10.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# apply_substrate_exchange + predict_substrate_HF_exchange                    #
# --------------------------------------------------------------------------- #

def test_apply_substrate_exchange_multiplicative():
    """apply is just a multiplication by the per-Z factor."""
    base_ie = 17.194  # representative O HF Koopmans IE
    for Z in [4, 8, 12, 16]:
        result = apply_substrate_exchange(base_ie, Z)
        assert result == pytest.approx(
            base_ie * substrate_exchange_factor(Z), rel=1e-12,
        )


def test_predict_HF_exchange_hydrogen_exact():
    """H: factor = 1, IE = -eps_HF * Hartree-eV = 13.606 eV."""
    z, ie = predict_substrate_HF_exchange(1)
    assert ie == pytest.approx(13.606, abs=1e-3)
    assert z == pytest.approx(1.0, abs=1e-3)


def test_predict_HF_exchange_helium_factor_one():
    """He: factor = 1 (no inner core), so HFx == HF Koopmans."""
    from src.stiff_medium.ionization_energy_test import predict_substrate_HF
    z_hfx, ie_hfx = predict_substrate_HF_exchange(2)
    z_hf, ie_hf = predict_substrate_HF(2)
    assert ie_hfx == pytest.approx(ie_hf, rel=1e-12)


def test_predict_HF_exchange_unknown_Z_raises():
    with pytest.raises(ValueError):
        predict_substrate_HF_exchange(19)


# --------------------------------------------------------------------------- #
# Mean error H..Ar — the headline result                                      #
# --------------------------------------------------------------------------- #

def test_HFx_mean_error_below_target():
    """Substrate-HF + exchange: mean abs err < 5% across H..Ar (target was 10%).
    Headline value: 2.78%.  Closes the K_rank 21% and HF Koopmans 6.4% gaps."""
    from src.stiff_medium.ionization_energy_test import (
        MEASURED_IE_EV, build_rows,
    )
    rows = build_rows(18)
    errs = [r.err_HFx_pct for r in rows]
    mean_err = sum(errs) / len(errs)
    assert mean_err < 5.0
    # Sentinel against regression
    assert mean_err == pytest.approx(2.78, abs=0.5)


def test_HFx_max_error_below_target():
    """Worst residual stays under 12% (currently O at 10.5%)."""
    from src.stiff_medium.ionization_energy_test import build_rows
    rows = build_rows(18)
    max_err = max(r.err_HFx_pct for r in rows)
    assert max_err < 12.0


# --------------------------------------------------------------------------- #
# Half-shell anomaly: N → O is the canonical drop                             #
# --------------------------------------------------------------------------- #

def test_half_shell_anomaly_N_to_O_positive():
    """N (2p³) → O (2p⁴): IE drops because O is past the half-shell.
    Substrate exchange must reproduce the sign of this anomaly."""
    drop = half_shell_anomaly_eV(7)  # IE(N) - IE(O)
    assert drop > 0.0, (
        f"N→O drop should be POSITIVE (IE drops at half-shell+1); got {drop:.3f} eV"
    )


def test_half_shell_anomaly_P_to_S_positive():
    """P (3p³) → S (3p⁴): IE drops at row-3 half-shell+1."""
    drop = half_shell_anomaly_eV(15)
    assert drop > 0.0, (
        f"P→S drop should be POSITIVE; got {drop:.3f} eV"
    )


def test_half_shell_anomaly_C_to_N_negative():
    """C → N is a normal IE rise (still in unpaired half-fill); drop should
    be NEGATIVE here (IE goes UP from C to N)."""
    drop = half_shell_anomaly_eV(6)  # IE(C) - IE(N)
    assert drop < 0.0


# --------------------------------------------------------------------------- #
# Audit table                                                                 #
# --------------------------------------------------------------------------- #

def test_audit_exchange_covers_h_through_ar():
    audit = audit_exchange(18)
    assert set(audit.keys()) == set(range(1, 19))
    for Z, entry in audit.items():
        assert isinstance(entry, ExchangeAudit)
        assert entry.Z == Z
        assert 0 <= entry.k_p <= 6
        assert 0 <= entry.k_s <= 2
        assert isinstance(entry.inner_core, bool)
        # Total factor should equal product of pieces (within floating-point)
        prod = (
            entry.half_shell_factor
            * entry.closed_s_pair_factor
            * entry.sp_factor
        )
        assert entry.total_factor == pytest.approx(prod, rel=1e-12)


def test_audit_exchange_oxygen_signature():
    """O's audit row: k_p=4, k_s=0, inner_core=True, half_shell<1, closed_s=1."""
    audit = audit_exchange(18)
    o = audit[8]
    assert o.k_p == 4
    assert o.k_s == 0
    assert o.inner_core is True
    assert o.half_shell_factor == pytest.approx(7.0 / 8.0, abs=1e-12)
    assert o.closed_s_pair_factor == 1.0
    assert o.sp_factor == 1.0


def test_audit_exchange_be_signature():
    """Be's audit row: k_p=0, k_s=2, inner_core=True, half_shell=1, closed_s=11/10."""
    audit = audit_exchange(18)
    be = audit[4]
    assert be.k_p == 0
    assert be.k_s == 2
    assert be.inner_core is True
    assert be.half_shell_factor == 1.0
    assert be.closed_s_pair_factor == pytest.approx(11.0 / 10.0, abs=1e-12)
    assert be.sp_factor == 1.0
