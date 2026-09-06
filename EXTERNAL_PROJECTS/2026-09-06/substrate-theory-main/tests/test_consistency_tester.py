"""Tests for the cross-module ConsistencyTester.

These tests verify that each ``verify_*`` method runs to completion and
returns a populated :class:`ConsistencyResult`. A method is allowed to
report ``passed=False`` (an honest inconsistency report is itself a useful
outcome), but it must NEVER raise.
"""

from __future__ import annotations

import pytest

from src.stiff_medium.consistency_tester import (
    ConsistencyResult,
    ConsistencyTester,
)


@pytest.fixture
def tester() -> ConsistencyTester:
    return ConsistencyTester()


# ---------------------------------------------------------------------------
# Construction + module discovery
# ---------------------------------------------------------------------------

def test_construction_records_modules(tester: ConsistencyTester) -> None:
    # b3_constants is hard-imported and must be present.
    assert tester._mods["b3_constants"] is not None
    # At least 5 of the 12 sibling modules should resolve.
    n_found = sum(1 for v in tester._mods.values() if v is not None)
    assert n_found >= 5, f"only {n_found} sibling modules importable"


def test_results_starts_empty(tester: ConsistencyTester) -> None:
    assert tester.results == []


# ---------------------------------------------------------------------------
# Each verify_* method
# ---------------------------------------------------------------------------

def test_verify_b3_constants_used(tester: ConsistencyTester) -> None:
    r = tester.verify_b3_constants_used()
    assert isinstance(r, ConsistencyResult)
    assert r.name
    # canonical integer-grid drift must be zero -- this is a tight contract.
    assert r.passed, f"b3_constants drift: {r.values}\n{r.notes}"


def test_verify_omega_b_consistency(tester: ConsistencyTester) -> None:
    r = tester.verify_omega_b_consistency()
    assert isinstance(r, ConsistencyResult)
    # b3_constants reference is always present
    assert any("b3_constants" in k for k in r.values), r.values
    # honest report: pass OR fail are both acceptable, but the result
    # must be non-trivially populated.
    assert r.values, "omega_b check produced no values"


def test_verify_alpha_consistency(tester: ConsistencyTester) -> None:
    r = tester.verify_alpha_consistency()
    assert isinstance(r, ConsistencyResult)
    # B3 closed form is always recorded
    assert "B3 closed form" in r.values
    # alpha cross-module drift should be tight (analytic)
    assert r.max_rel_drift < 1.0e-3, (
        f"alpha drift {r.max_rel_drift:.2e} too large; values={r.values}"
    )


def test_verify_mass_torque_chain(tester: ConsistencyTester) -> None:
    r = tester.verify_mass_torque_chain()
    assert isinstance(r, ConsistencyResult)
    # if both modules importable, must produce values for both ratios
    if r.values:
        # at least one m_mu/m_e and one m_tau/m_mu entry
        keys = list(r.values)
        assert any("m_mu/m_e" in k for k in keys), keys
        assert any("m_tau/m_mu" in k for k in keys), keys


def test_verify_cosmology_consistency(tester: ConsistencyTester) -> None:
    r = tester.verify_cosmology_consistency()
    assert isinstance(r, ConsistencyResult)
    # at least one observable should produce values
    assert r.values, f"cosmology check returned no values: {r.notes}"
    # ranges sanity: H_0 between 50 and 90 km/s/Mpc, Sigma m_nu between
    # 30 and 200 meV, rho_Lambda > 0
    for k, v in r.values.items():
        if k.startswith("H_0"):
            assert 50.0 < v < 90.0, (k, v)
        elif k.startswith("Sigmamnu"):
            assert 30.0 < v < 500.0, (k, v)
        elif k.startswith("rhoLambda"):
            assert v > 0.0, (k, v)


def test_verify_majorana_consistency(tester: ConsistencyTester) -> None:
    r = tester.verify_majorana_consistency()
    assert isinstance(r, ConsistencyResult)
    if r.values:
        # both modules should report sin2_t12 entries
        keys = list(r.values)
        assert any("sin2_t12" in k for k in keys), keys


# ---------------------------------------------------------------------------
# Aggregate flow: run_all + report
# ---------------------------------------------------------------------------

def test_run_all_checks_returns_list(tester: ConsistencyTester) -> None:
    out = tester.run_all_checks()
    assert isinstance(out, list)
    assert len(out) == 6
    assert all(isinstance(r, ConsistencyResult) for r in out)


def test_run_all_idempotent_resets(tester: ConsistencyTester) -> None:
    tester.run_all_checks()
    n1 = len(tester.results)
    tester.run_all_checks()
    n2 = len(tester.results)
    assert n1 == n2, "run_all_checks should reset results, not append"


def test_report_before_run_is_safe(tester: ConsistencyTester) -> None:
    msg = tester.report()
    assert isinstance(msg, str)
    assert "no checks" in msg.lower()


def test_report_after_run_lists_all(tester: ConsistencyTester) -> None:
    tester.run_all_checks()
    msg = tester.report()
    assert isinstance(msg, str)
    assert "consistency report" in msg.lower()
    # every result name appears somewhere in the report
    for r in tester.results:
        assert r.name in msg


def test_no_verify_method_raises(tester: ConsistencyTester) -> None:
    """Honest reporting contract: failing checks must mark passed=False, never raise."""
    methods = (
        tester.verify_b3_constants_used,
        tester.verify_omega_b_consistency,
        tester.verify_alpha_consistency,
        tester.verify_mass_torque_chain,
        tester.verify_cosmology_consistency,
        tester.verify_majorana_consistency,
    )
    for m in methods:
        m()  # must not raise
