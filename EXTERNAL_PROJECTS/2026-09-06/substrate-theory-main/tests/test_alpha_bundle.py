"""Tests for stiff_medium.alpha_bundle."""

from __future__ import annotations

import math

import pytest

from stiff_medium.alpha_bundle import (
    ALPHA_CODATA,
    E_HL_OBSERVED,
    PI,
    alpha_from_beta,
    audit,
    beta_from_alpha,
    candidate_beta_values,
)


def test_forward_inverse_consistency() -> None:
    """alpha_from_beta and beta_from_alpha must round-trip exactly."""
    alpha = ALPHA_CODATA
    beta = beta_from_alpha(alpha)
    alpha_back = alpha_from_beta(beta)
    assert math.isclose(alpha, alpha_back, rel_tol=1e-15), (
        f"Round-trip α → β → α failed: {alpha} → {beta} → {alpha_back}"
    )


def test_beta_observed_value() -> None:
    """β_observed = √(4π α_CODATA) should be ≈ 0.3028..."""
    expected = math.sqrt(4.0 * PI * ALPHA_CODATA)
    assert math.isclose(E_HL_OBSERVED, expected, rel_tol=1e-12)
    # Sanity: should be the standard Heaviside-Lorentz electric charge
    assert 0.302 < E_HL_OBSERVED < 0.304


def test_alpha_from_beta_matches_codata_at_observed_beta() -> None:
    """alpha_from_beta(β_observed) must reproduce α_CODATA."""
    alpha = alpha_from_beta(E_HL_OBSERVED)
    assert math.isclose(alpha, ALPHA_CODATA, rel_tol=1e-12)


def test_candidates_are_sorted_by_residual() -> None:
    """candidate_beta_values must return list sorted by ascending residual."""
    cands = candidate_beta_values()
    residuals = [c.residual_pct for c in cands]
    assert residuals == sorted(residuals), (
        f"Candidates not sorted: {residuals}"
    )


def test_candidate_target_is_zero_residual() -> None:
    """The 'β_observed' candidate must have zero residual by construction."""
    cands = candidate_beta_values()
    target = next(c for c in cands if c.name == "β_observed")
    assert target.residual_pct == 0.0


def test_audit_excludes_target_from_best_candidate() -> None:
    """audit().best_candidate must NOT be the tautological 'β_observed' entry."""
    a = audit()
    assert a.best_candidate.name != "β_observed", (
        "audit() should report best DERIVED candidate, not the tautology"
    )


def test_audit_best_candidate_within_known_bound() -> None:
    """The best derived candidate is currently 1/√11 or √(2/π)/√7, both <1% off.

    If a future commit improves the candidate list, this bound should tighten.
    """
    a = audit()
    assert a.best_residual_pct < 1.0, (
        f"Best derived candidate now has residual {a.best_residual_pct:.4f}% — "
        "either the candidate list regressed or a clean derivation is in reach. "
        "Update the test bound."
    )


def test_pi_inverse_is_about_10pct_off() -> None:
    """The simplest topological candidate, 1/π, should give ~10% residual.

    This is a sanity check that the basic Möbius-holonomy-only guess
    doesn't accidentally land at α — which would be too good a coincidence.
    """
    alpha_pi = alpha_from_beta(1.0 / PI)
    inv_alpha_pi = 1.0 / alpha_pi
    # 1/(4π³) ≈ 1/124.025
    assert math.isclose(inv_alpha_pi, 4.0 * PI**3, rel_tol=1e-12)
    residual_pct = 100.0 * abs(alpha_pi - ALPHA_CODATA) / ALPHA_CODATA
    assert 10.0 < residual_pct < 11.0


@pytest.mark.parametrize(
    "alpha_input",
    [1.0 / 137.036, 1.0 / 100.0, 1.0 / 200.0, 0.01],
)
def test_round_trip_at_various_alpha(alpha_input: float) -> None:
    """Round-trip must hold for any reasonable α value."""
    beta = beta_from_alpha(alpha_input)
    alpha_back = alpha_from_beta(beta)
    assert math.isclose(alpha_input, alpha_back, rel_tol=1e-14)
