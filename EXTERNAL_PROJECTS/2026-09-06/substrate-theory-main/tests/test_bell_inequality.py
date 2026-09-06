"""Tests for stiff_medium.bell_inequality."""

from __future__ import annotations

import math

import numpy as np
import pytest

from stiff_medium.bell_inequality import (
    CHSH_CLASSICAL_BOUND,
    TSIRELSON_BOUND,
    BellInequality,
)


@pytest.fixture
def bell() -> BellInequality:
    return BellInequality(seed=42)


def test_classical_bound_is_two(bell: BellInequality) -> None:
    assert bell.chsh_classical() == pytest.approx(2.0)
    assert bell.chsh_classical() <= 2.0


def test_quantum_bound_is_tsirelson(bell: BellInequality) -> None:
    assert bell.chsh_quantum() == pytest.approx(2.0 * math.sqrt(2.0))
    assert bell.chsh_quantum() == pytest.approx(2.8284271247, rel=1e-9)


def test_substrate_matches_quantum(bell: BellInequality) -> None:
    assert bell.chsh_substrate() == pytest.approx(bell.chsh_quantum())
    assert bell.chsh_substrate() == pytest.approx(TSIRELSON_BOUND)


def test_quantum_violates_classical() -> None:
    assert TSIRELSON_BOUND > CHSH_CLASSICAL_BOUND


@pytest.mark.parametrize(
    "theta_a,theta_b,expected",
    [
        (0.0, 0.0, -1.0),
        (0.0, math.pi, 1.0),
        (0.0, math.pi / 2.0, 0.0),
        (math.pi / 4.0, -math.pi / 4.0, -math.cos(math.pi / 2.0)),
    ],
)
def test_singlet_correlation_form(
    bell: BellInequality, theta_a: float, theta_b: float, expected: float
) -> None:
    assert bell.singlet_state_correlation(theta_a, theta_b) == pytest.approx(
        expected, abs=1e-12
    )


def test_singlet_correlation_matches_minus_cos(bell: BellInequality) -> None:
    rng = np.random.default_rng(0)
    a = rng.uniform(-np.pi, np.pi, size=50)
    b = rng.uniform(-np.pi, np.pi, size=50)
    for ta, tb in zip(a, b):
        assert bell.singlet_state_correlation(float(ta), float(tb)) == pytest.approx(
            -math.cos(ta - tb), abs=1e-12
        )


def test_optimal_angles_saturate_tsirelson(bell: BellInequality) -> None:
    a1, a2, b1, b2 = bell.optimal_chsh_angles()
    s = bell.chsh_value(a1, a2, b1, b2)
    assert abs(s) == pytest.approx(TSIRELSON_BOUND, abs=1e-12)


def test_bell_test_simulation_recovers_tsirelson(bell: BellInequality) -> None:
    out = bell.bell_test_simulation(n_pairs=200_000)
    # ~5 sigma window using analytic SE estimate.
    assert abs(out["S"] - TSIRELSON_BOUND) < 5.0 * out["S_error"]
    assert out["S"] > CHSH_CLASSICAL_BOUND  # violates classical


def test_loophole_free_experiments_violate_classical(bell: BellInequality) -> None:
    exps = bell.loophole_free_compliance()
    assert len(exps) >= 4
    for ex in exps:
        assert ex.s_value > CHSH_CLASSICAL_BOUND
        assert ex.sigma_above_classical >= 2.0


def test_compare_to_experiments_structure(bell: BellInequality) -> None:
    cmp = bell.compare_to_experiments()
    assert cmp["classical_bound"] == pytest.approx(2.0)
    assert cmp["quantum_bound"] == pytest.approx(TSIRELSON_BOUND)
    assert cmp["substrate_prediction"] == pytest.approx(TSIRELSON_BOUND)
    rows = cmp["experiments"]
    assert isinstance(rows, list) and len(rows) >= 4
    for row in rows:
        assert row["violates_classical"] is True


def test_report_runs_and_mentions_key_terms(bell: BellInequality) -> None:
    text = bell.report()
    assert isinstance(text, str)
    for needle in ("Tsirelson", "Substrate", "Classical", "Monte Carlo"):
        assert needle.lower() in text.lower()
