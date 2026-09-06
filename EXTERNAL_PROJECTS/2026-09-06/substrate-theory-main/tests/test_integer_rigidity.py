"""
test_integer_rigidity.py
========================

Tests for the B3 integer rigidity framework.

Verifies:
  (a) n_M = 268 uniquely minimises the chi^2 over a +/- window.
  (b) N_BAM = 6 vs 7 gives wildly different eps_face / deuteron BE.
  (c) K_pair = 2 vs 3 catastrophically breaks m_mu / m_e.
  (d) Canonical baseline chi^2 is small (predictions match observations).
  (e) Each named perturb_* method runs and produces a sensible result.
"""

import math
import pytest

from src.stiff_medium.integer_rigidity import (
    IntegerRigidity,
    OBSERVED,
    CANONICAL_INTEGERS,
    predict_m_mu_over_m_e,
    predict_m_tau_over_m_e,
    predict_deuteron_BE,
    predict_eps_face,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def rig():
    return IntegerRigidity()


# ---------------------------------------------------------------------------
# Baseline sanity
# ---------------------------------------------------------------------------

def test_canonical_chi2_is_small(rig):
    """The canonical integers should match observation to <1% per channel."""
    baseline = rig.baseline_chi2()
    # 5 observables * (0.01)^2 = 5e-4 worst case if every channel were 1% off
    assert baseline < 5e-3, f"baseline chi^2 = {baseline:.3e} too large"


def test_canonical_predictions_match_pdg():
    """Sanity: each canonical predictor lands within tolerance of observed."""
    pred_mu = predict_m_mu_over_m_e(CANONICAL_INTEGERS)
    obs_mu = OBSERVED["m_mu_over_m_e"]
    assert abs(pred_mu - obs_mu) / obs_mu < 5e-3

    pred_d = predict_deuteron_BE(CANONICAL_INTEGERS)
    obs_d = OBSERVED["deuteron_BE"]
    assert abs(pred_d - obs_d) / obs_d < 1e-2

    assert predict_eps_face(CANONICAL_INTEGERS) == pytest.approx(1.0 / 6.0)


# ---------------------------------------------------------------------------
# (a) n_M = 268 is the unique minimum
# ---------------------------------------------------------------------------

def test_n_M_uniquely_minimises_chi2(rig):
    u = rig.verify_uniqueness("n_M", window=5)
    assert u["argmin_delta"] == 0, (
        f"n_M argmin shifted to delta={u['argmin_delta']} (chi^2={u['argmin_chi2']:.3e})"
    )
    assert u["is_unique_minimum"], f"non-unique minimum, ties at {u['ties']}"


def test_n_M_perturbation_degrades_predictions(rig):
    """Shifting n_M by +/-1 must materially worsen chi^2."""
    baseline = rig.baseline_chi2()
    for d in (-1, +1):
        r = rig.perturb_n_M(d)
        assert r.total_chi2 > baseline * 5, (
            f"n_M delta={d}: chi^2={r.total_chi2:.3e} only "
            f"{r.degrade_factor(baseline):.2f}x baseline"
        )


# ---------------------------------------------------------------------------
# (b) N_BAM = 6 vs 7 changes eps_face / deuteron meaningfully
# ---------------------------------------------------------------------------

def test_N_BAM_perturbation_breaks_eps_face(rig):
    r = rig.perturb_N_BAM(+1)         # N_BAM = 7
    eps_row = next(o for o in r.observables if o["observable"] == "eps_face")
    # canonical baseline 1/6 = 0.16667 vs 1/7 = 0.14286 -> 14.3% off
    assert eps_row["rel_err"] > 0.10, (
        f"eps_face at N_BAM=7 only {eps_row['rel_err']:.3%} off canonical"
    )


def test_N_BAM_perturbation_breaks_deuteron(rig):
    r = rig.perturb_N_BAM(+1)         # N_BAM = 7
    d_row = next(o for o in r.observables if o["observable"] == "deuteron_BE")
    # 200 / (45*7/3) = 200/105 = 1.905 MeV vs observed 2.225 -> ~14% off
    assert d_row["rel_err"] > 0.10, (
        f"deuteron_BE at N_BAM=7 only {d_row['rel_err']:.3%} off"
    )


# ---------------------------------------------------------------------------
# (c) K_pair = 2 vs 3 catastrophically breaks the muon ratio
# ---------------------------------------------------------------------------

def test_K_pair_perturbation_breaks_muon(rig):
    """K_pair = 3 makes m_mu/m_e drop from 206 to ~exp(268/(81 pi)) ~ exp(1.05) ~ 2.86."""
    r = rig.perturb_K_pair(+1)
    mu_row = next(o for o in r.observables if o["observable"] == "m_mu/m_e")
    # observed 206.77; perturbed ~ 2.86 -> rel_err ~ 0.986
    assert mu_row["rel_err"] > 0.5, (
        f"K_pair=3: m_mu/m_e only {mu_row['rel_err']:.3%} off"
    )
    # And K_pair=1 makes the exponent enormous (n_M/pi ~ 85 -> exp blows up)
    r2 = rig.perturb_K_pair(-1)
    mu_row2 = next(o for o in r2.observables if o["observable"] == "m_mu/m_e")
    assert mu_row2["rel_err"] > 1.0, (
        f"K_pair=1: m_mu/m_e only {mu_row2['rel_err']:.3%} off"
    )


def test_K_pair_uniqueness(rig):
    u = rig.verify_uniqueness("K_pair", window=2)
    # K_pair must minimise at canonical 2 (delta=0)
    assert u["argmin_delta"] == 0, (
        f"K_pair argmin shifted to delta={u['argmin_delta']}"
    )


# ---------------------------------------------------------------------------
# (d) K_rank perturbation breaks Cabibbo
# ---------------------------------------------------------------------------

def test_K_rank_perturbation_breaks_cabibbo(rig):
    r = rig.perturb_K_rank(+1)        # K_rank = 6 -> 1/sqrt(30) = 0.1826
    cab_row = next(o for o in r.observables if o["observable"] == "cabibbo")
    # observed 0.22534; perturbed 0.1826 -> ~19% off
    assert cab_row["rel_err"] > 0.10, (
        f"K_rank=6: cabibbo only {cab_row['rel_err']:.3%} off"
    )


# ---------------------------------------------------------------------------
# (d') n_R perturbation must change m_tau/m_e (not silently no-op)
# ---------------------------------------------------------------------------

def test_n_R_perturbation_changes_tau_prediction():
    """Audit guard: n_R must enter m_tau/m_e EXPLICITLY.

    Before the audit fix, predict_m_tau_over_m_e ignored its `n_R`
    argument and used K_rank/K_pair instead -- so perturbing n_R was a
    no-op. The reflection-counted reformulation
        log(tau/mu) = n_M/(K_pair^4 pi) - (n_R - R)/(K_pair*(K_pair+1))
    is numerically identical at the canonical integers (15/6 = 5/2)
    but DOES respond to a unit shift in n_R. This test guards against
    n_R silently becoming inert again.
    """
    canon = dict(CANONICAL_INTEGERS)
    canon["R"] = 3
    pred_canonical = predict_m_tau_over_m_e(canon)

    perturbed = dict(canon)
    perturbed["n_R"] = canon["n_R"] + 1   # 18 -> 19, leave n_M fixed
    pred_perturbed = predict_m_tau_over_m_e(perturbed)

    rel_change = abs(pred_perturbed - pred_canonical) / pred_canonical
    assert rel_change > 1e-3, (
        f"n_R perturbation produced rel change {rel_change:.3e} -- "
        f"n_R is not actually used in m_tau/m_e predictor (audit regression)."
    )


def test_n_R_perturbation_breaks_tau_in_engine(rig):
    """Engine-level guard: rig.perturb_n_R must shift the m_tau/m_e row."""
    r = rig.perturb_n_R(+1)
    tau_row = next(o for o in r.observables
                   if o["observable"] == "m_tau/m_e")
    canonical_pred = tau_row["canonical_pred"]
    perturbed_pred = tau_row["perturbed_pred"]
    rel_change = abs(perturbed_pred - canonical_pred) / canonical_pred
    assert rel_change > 1e-3, (
        f"n_R+1 perturbation only shifted m_tau/m_e by {rel_change:.3e} "
        f"({canonical_pred} -> {perturbed_pred}); n_R inactive."
    )


def test_n_R_canonical_matches_K_rank_form():
    """Verify the (n_R - R)/(K_pair*(K_pair+1)) rewrite preserves numerics.

    At canonical integers (n_R=18, R=3, K_pair=2, K_rank=5) we must have
    (n_R - R)/(K_pair*(K_pair+1)) = 15/6 = 2.5 = K_rank/K_pair. This is
    the load-bearing integer identity that lets us promote n_R to an
    active observable input WITHOUT changing the canonical prediction.
    """
    n_R, R_int, K_pair, K_rank = 18, 3, 2, 5
    lhs = (n_R - R_int) / (K_pair * (K_pair + 1))
    rhs = K_rank / K_pair
    assert lhs == rhs == 2.5


# ---------------------------------------------------------------------------
# (e) API smoke: every named perturb_* method works
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method,delta", [
    ("perturb_n_M",   +1),
    ("perturb_N_BAM", +1),
    ("perturb_K_pair", +1),
    ("perturb_K_rank", +1),
    ("perturb_n_R",    +1),
    ("perturb_n_A",    +1),
])
def test_perturb_methods_runnable(rig, method, delta):
    r = getattr(rig, method)(delta)
    assert r.delta == delta
    assert r.total_chi2 >= 0.0
    assert isinstance(r.observables, list)


def test_rigidity_report_runs(rig, capsys):
    text = rig.rigidity_report(deltas=(-1, +1))
    assert "B3 Integer Rigidity Report" in text
    assert "n_M" in text
