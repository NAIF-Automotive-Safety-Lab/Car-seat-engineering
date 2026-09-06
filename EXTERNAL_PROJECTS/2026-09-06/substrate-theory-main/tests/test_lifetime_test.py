"""Tests for lifetime_test — substrate vs PDG 2024 lifetime audit."""

from __future__ import annotations

import math

import pytest

from stiff_medium.lifetime_test import (
    PDG_TABLE,
    LifetimePrediction,
    average_residual,
    kaon_charged_lifetime,
    kaon_long_lifetime,
    kaon_short_lifetime,
    muon_lifetime,
    neutron_bottle_beam_verdict,
    neutron_lifetime,
    pion_charged_lifetime,
    pion_neutral_lifetime,
    run_full_test,
    tau_lifetime,
)


# ---------------------------------------------------------------------------
# (A) Each particle lifetime is within 5% of PDG
# ---------------------------------------------------------------------------

def test_muon_within_1pct() -> None:
    p = muon_lifetime()
    assert abs(p.fractional_error) < 0.01, (
        f"muon: τ_sub = {p.tau_substrate_s*1e6:.4f} μs, "
        f"err = {p.percent_error:.2f}%"
    )


def test_tau_within_5pct() -> None:
    p = tau_lifetime()
    assert abs(p.fractional_error) < 0.05, (
        f"tau: τ_sub = {p.tau_substrate_s*1e15:.2f} fs, "
        f"err = {p.percent_error:.2f}%"
    )


def test_pion_charged_within_3pct() -> None:
    p = pion_charged_lifetime()
    assert abs(p.fractional_error) < 0.03, (
        f"π±: τ_sub = {p.tau_substrate_s*1e9:.3f} ns, "
        f"err = {p.percent_error:.2f}%"
    )


def test_pion_neutral_within_3pct() -> None:
    p = pion_neutral_lifetime()
    assert abs(p.fractional_error) < 0.03, (
        f"π⁰: τ_sub = {p.tau_substrate_s*1e17:.3f} × 10⁻¹⁷ s, "
        f"err = {p.percent_error:.2f}%"
    )


def test_kaon_charged_within_5pct() -> None:
    p = kaon_charged_lifetime()
    assert abs(p.fractional_error) < 0.05, (
        f"K±: τ_sub = {p.tau_substrate_s*1e9:.3f} ns, "
        f"err = {p.percent_error:.2f}%"
    )


def test_kaon_short_within_5pct() -> None:
    p = kaon_short_lifetime()
    assert abs(p.fractional_error) < 0.05, (
        f"K_S: τ_sub = {p.tau_substrate_s*1e12:.3f} ps, "
        f"err = {p.percent_error:.2f}%"
    )


def test_kaon_long_within_5pct() -> None:
    p = kaon_long_lifetime()
    assert abs(p.fractional_error) < 0.05, (
        f"K_L: τ_sub = {p.tau_substrate_s*1e9:.3f} ns, "
        f"err = {p.percent_error:.2f}%"
    )


def test_neutron_within_2pct() -> None:
    p = neutron_lifetime()
    assert abs(p.fractional_error) < 0.02, (
        f"neutron: τ_sub = {p.tau_substrate_s:.2f} s, "
        f"err = {p.percent_error:.2f}%"
    )


# ---------------------------------------------------------------------------
# (B) Aggregate statistics
# ---------------------------------------------------------------------------

def test_aggregate_mean_residual_under_3pct() -> None:
    preds = run_full_test()
    stats = average_residual(preds)
    assert stats["mean_abs_pct_error"] < 3.0, (
        f"mean |% err| = {stats['mean_abs_pct_error']:.2f}%, expected < 3%"
    )


def test_all_8_within_5pct() -> None:
    preds = run_full_test()
    stats = average_residual(preds)
    assert stats["n_within_5pct"] == 8, (
        f"only {stats['n_within_5pct']} / 8 within 5% of PDG"
    )


def test_at_least_5_within_2pct() -> None:
    preds = run_full_test()
    stats = average_residual(preds)
    assert stats["n_within_2pct"] >= 5, (
        f"only {stats['n_within_2pct']} / 8 within 2% of PDG"
    )


def test_log_residual_bounded() -> None:
    """log10(τ_sub/τ_PDG) bounded — substrate gets the order-of-magnitude
    right across 22 decades of lifetime (8.4×10⁻¹⁷ s to 887 s)."""
    preds = run_full_test()
    stats = average_residual(preds)
    assert stats["mean_abs_log10_residual"] < 0.05, (
        f"mean |log10 residual| = {stats['mean_abs_log10_residual']:.4f}"
    )


# ---------------------------------------------------------------------------
# (C) Lifetime span — substrate covers 22 decades
# ---------------------------------------------------------------------------

def test_lifetime_dynamic_range_19_decades() -> None:
    """Substrate spans ~19 decades from π⁰ (8.5e-17 s) to neutron (892 s)."""
    preds = run_full_test()
    log_taus = [math.log10(p.tau_substrate_s) for p in preds]
    span = max(log_taus) - min(log_taus)
    assert span > 18.0, (
        f"substrate lifetime span = {span:.1f} decades, expected > 18"
    )


# ---------------------------------------------------------------------------
# (D) Bottle vs beam puzzle — substrate sits closer to BEAM
# ---------------------------------------------------------------------------

def test_neutron_substrate_favors_beam_method() -> None:
    bb = neutron_bottle_beam_verdict()
    assert bb["substrate_favors"] == "beam"
    assert bb["nsigma_to_beam"] < bb["nsigma_to_bottle"]


def test_neutron_bottle_beam_tension_at_least_3sigma() -> None:
    bb = neutron_bottle_beam_verdict()
    assert bb["bottle_beam_tension_sigma"] >= 3.0, (
        f"bottle/beam tension = {bb['bottle_beam_tension_sigma']:.1f}σ"
    )


# ---------------------------------------------------------------------------
# (E) PDG values are physical
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", list(PDG_TABLE.keys()))
def test_pdg_values_positive(name: str) -> None:
    pdg = PDG_TABLE[name]
    assert pdg.tau_central > 0
    assert pdg.tau_sigma > 0
    assert pdg.tau_sigma < pdg.tau_central  # σ < value (precision check)


# ---------------------------------------------------------------------------
# (F) Each prediction returns the right structure
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("fn", [
    muon_lifetime,
    tau_lifetime,
    pion_charged_lifetime,
    pion_neutral_lifetime,
    kaon_charged_lifetime,
    kaon_short_lifetime,
    kaon_long_lifetime,
    neutron_lifetime,
])
def test_prediction_structure(fn) -> None:
    p = fn()
    assert isinstance(p, LifetimePrediction)
    assert p.tau_substrate_s > 0
    assert p.tau_pdg_s > 0
    assert isinstance(p.method, str)
    assert isinstance(p.passes, bool)
    assert math.isfinite(p.fractional_error)
