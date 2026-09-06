"""Tests for src/stiff_medium/mobius_k4_numerical.py."""

import math

import numpy as np
import pytest

from src.stiff_medium.mobius_k4_numerical import (
    ANALYTIC,
    combinatorial_ratio,
    k4_edges,
    k4_faces,
    k4_vertices,
    mobius_section,
    monte_carlo_ratio,
    quadrature_ratio,
    run_all,
    tetrahedron_volume,
    verify_mobius_twist,
)


def test_k4_geometry():
    V = k4_vertices()
    assert V.shape == (4, 3)
    assert len(k4_edges()) == 6
    assert len(k4_faces()) == 4
    # Regular tetrahedron with edge length sqrt(2) -> volume = sqrt(2)^3/(6 sqrt(2)) ... use computed
    assert tetrahedron_volume(V) > 0.0


def test_mobius_twist_explicit():
    assert verify_mobius_twist(n=33)
    theta = np.linspace(-3.0, 3.0, 11)
    assert np.allclose(mobius_section(theta + 2 * math.pi), -mobius_section(theta), atol=1e-14)


def test_combinatorial_method():
    ratio, info = combinatorial_ratio()
    assert info["n_total"] == 12
    assert info["n_removed"] == 1
    assert ratio == pytest.approx(ANALYTIC, rel=1e-15)
    assert abs(ratio - ANALYTIC) / ANALYTIC < 1e-3


def test_quadrature_method():
    ratio, info = quadrature_ratio()
    assert ratio == pytest.approx(ANALYTIC, rel=1e-12)
    assert abs(ratio - ANALYTIC) / ANALYTIC < 1e-3
    assert info["trivial_total"] > info["mobius_total"]


def test_monte_carlo_method():
    ratio, info = monte_carlo_ratio(n_samples=2_000_000, seed=42)
    assert abs(ratio - ANALYTIC) / ANALYTIC < 1e-3
    assert info["std_err"] < 5e-3


def test_methods_agree():
    report = run_all(n_samples=4_000_000, seed=2026)
    methods = [report.combinatorial, report.monte_carlo, report.quadrature]
    spread = max(methods) - min(methods)
    # Combinatorial and quadrature agree exactly; MC is statistical.
    assert abs(report.combinatorial - report.quadrature) / ANALYTIC < 1e-12
    # All three within 0.1% of analytic.
    for m in methods:
        assert abs(m - ANALYTIC) / ANALYTIC < 1e-3
    # Spread is bounded by the MC standard error (mild tolerance: 5 sigma).
    assert spread / ANALYTIC < 5e-3


def test_report_pretty_contains_values():
    report = run_all(n_samples=200_000, seed=7)
    text = report.pretty()
    assert "11/12" in text
    assert "combinatorial" in text
    assert "Monte Carlo" in text
    assert "quadrature" in text
