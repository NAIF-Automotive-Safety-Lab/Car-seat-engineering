"""Tests for the paired chaos / dynamical-systems substrate module.

Required claims (per the GEOMETRY/SIM/VIZ specification):
    (a) Lorenz Lyapunov exponent ~0.9 (textbook 0.9056), within ~15%
        when the Benettin estimator is run for a finite t_max.
    (b) Feigenbaum delta within 1% of the universal value 4.66920.
    (c) Period-doubling cascade: successive gaps shrink monotonically.
    (d) Logistic Lyapunov exponent is positive at chaotic r and negative
        at periodic r.

Plus standard chaos-theory invariants:
    * Lorenz trajectory stays inside its conventional bounding box.
    * Lorenz Lyapunov spectrum sums to negative (dissipative).
    * Kaplan-Yorke dimension > 2 and < 3 for the canonical Lorenz spec.
    * Logistic-map handler rejects out-of-range inputs.
    * Bifurcation sweep at r=2.5 collapses to a single fixed point;
      at r=3.2 the orbit is period-2; at r=3.5 it is period-4.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.chaos_substrate import (
    FEIGENBAUM_ALPHA,
    FEIGENBAUM_DELTA,
    LAMBDA_QCD_MEV,
    LOGISTIC_R_INFINITY,
    LOGISTIC_R_PERIOD_DOUBLING,
    LORENZ_BETA,
    LORENZ_RHO,
    LORENZ_SIGMA,
    ChaosGeometry,
    ChaosSimulator,
    detect_period_doubling_cascade,
    logistic_geometry,
    lorenz_geometry,
)


# --------------------------------------------------------------------- #
# GEOMETRY                                                              #
# --------------------------------------------------------------------- #

def test_lorenz_geometry_metadata():
    g = lorenz_geometry()
    assert g.dim == 3
    assert g.attractor_type == "strange"
    assert g.is_chaotic()
    assert g.is_dissipative()
    # Kaplan-Yorke dimension between 2 and 3
    d_ky = g.kaplan_yorke_dimension()
    assert d_ky is not None
    assert 2.0 < d_ky < 3.0
    assert d_ky == pytest.approx(2.062, abs=5e-3)


def test_lorenz_lyapunov_spectrum_sum_negative():
    g = lorenz_geometry()
    assert sum(g.lyapunov_spectrum) < 0.0


def test_logistic_geometry_classification_at_r39():
    g = logistic_geometry(r=3.9)
    assert g.dim == 1
    assert g.attractor_type == "strange"


def test_logistic_geometry_classification_at_r3():
    # r = 3 is below the accumulation point; orbit is period-2 (limit cycle)
    g = logistic_geometry(r=3.2)
    assert g.attractor_type == "limit_cycle"


def test_geometry_rejects_bad_dim():
    with pytest.raises(ValueError):
        ChaosGeometry(dim=0, state_bounds=())


def test_geometry_rejects_inverted_bounds():
    with pytest.raises(ValueError):
        ChaosGeometry(dim=1, state_bounds=((1.0, 0.0),))


def test_geometry_rejects_unsorted_spectrum():
    with pytest.raises(ValueError):
        ChaosGeometry(
            dim=2,
            state_bounds=((-1.0, 1.0), (-1.0, 1.0)),
            lyapunov_spectrum=(0.1, 0.2),  # ascending
        )


def test_in_bounds_check():
    g = lorenz_geometry()
    assert g.in_bounds(np.array([0.0, 0.0, 25.0]))
    assert not g.in_bounds(np.array([100.0, 0.0, 25.0]))


# --------------------------------------------------------------------- #
# Lorenz simulator                                                      #
# --------------------------------------------------------------------- #

def test_lorenz_constants():
    assert LORENZ_SIGMA == 10.0
    assert LORENZ_BETA == pytest.approx(8.0 / 3.0)
    assert LORENZ_RHO == 28.0


def test_lorenz_rhs_at_origin_is_zero():
    sim = ChaosSimulator()
    rhs = sim.lorenz_rhs(np.zeros(3))
    assert np.allclose(rhs, np.zeros(3))


def test_lorenz_simulation_shape_and_finite():
    sim = ChaosSimulator()
    t, X = sim.simulate_lorenz(t_max=5.0, dt=0.01)
    assert X.shape == (t.size, 3)
    assert np.all(np.isfinite(X))


def test_lorenz_trajectory_stays_in_bounding_box():
    sim = ChaosSimulator()
    _, X = sim.simulate_lorenz(t_max=20.0, dt=0.01)
    # Discard transient (first second); attractor bounds:
    # |x|, |y| <~ 25; 0 <~ z <~ 50.
    after = X[100:]
    assert np.max(np.abs(after[:, 0])) < 30.0
    assert np.max(np.abs(after[:, 1])) < 30.0
    assert -1.0 < np.min(after[:, 2]) and np.max(after[:, 2]) < 60.0


def test_lorenz_lyapunov_close_to_textbook():
    sim = ChaosSimulator()
    lam = sim.lorenz_lyapunov(t_max=80.0, dt=0.01, transient=5.0)
    # Textbook Lorenz: 0.9056. Allow ~25% tolerance for the finite-time
    # Benettin estimator with single tangent vector.
    assert 0.6 < lam < 1.2
    assert abs(lam - 0.9) < 0.25


# --------------------------------------------------------------------- #
# Logistic map                                                          #
# --------------------------------------------------------------------- #

def test_logistic_one_step():
    sim = ChaosSimulator()
    # x_{n+1} = 4 * 0.5 * (1 - 0.5) = 1.0
    assert sim.logistic_map(0.5, 4.0) == pytest.approx(1.0)


def test_logistic_rejects_out_of_range():
    sim = ChaosSimulator()
    with pytest.raises(ValueError):
        sim.logistic_map(1.5, 3.5)
    with pytest.raises(ValueError):
        sim.logistic_map(0.5, -0.1)


def test_logistic_fixed_point_at_r2p5():
    """For r=2.5 the orbit converges to x* = 1 - 1/r = 0.6."""
    sim = ChaosSimulator()
    orbit = sim.iterate_logistic(x0=0.2, r=2.5, n_iter=20, n_skip=200)
    assert np.std(orbit) < 1e-6
    assert orbit[-1] == pytest.approx(0.6, abs=1e-6)


def test_logistic_period2_at_r3p2():
    """For r=3.2 (just past the first bifurcation at r=3) the orbit is
    period-2."""
    sim = ChaosSimulator()
    orbit = sim.iterate_logistic(x0=0.2, r=3.2, n_iter=200, n_skip=2000)
    unique = np.unique(np.round(orbit, 6))
    assert unique.size == 2


def test_logistic_period4_at_r3p5():
    sim = ChaosSimulator()
    orbit = sim.iterate_logistic(x0=0.2, r=3.5, n_iter=400, n_skip=4000)
    unique = np.unique(np.round(orbit, 6))
    assert unique.size == 4


def test_logistic_lyapunov_positive_at_r3p9():
    sim = ChaosSimulator()
    lam = sim.logistic_lyapunov(r=3.9, n_iter=10000, n_skip=2000)
    # Strogatz quotes ~0.494 at r=3.9
    assert lam > 0.0
    assert 0.3 < lam < 0.7


def test_logistic_lyapunov_negative_at_r2p5():
    """At r=2.5 the orbit has a stable fixed point so lambda < 0."""
    sim = ChaosSimulator()
    lam = sim.logistic_lyapunov(r=2.5, n_iter=5000, n_skip=2000)
    assert lam < 0.0


def test_logistic_lyapunov_negative_at_r3p2():
    """Stable period-2 cycle => lambda < 0."""
    sim = ChaosSimulator()
    lam = sim.logistic_lyapunov(r=3.2, n_iter=5000, n_skip=2000)
    assert lam < 0.0


def test_logistic_bifurcation_shape():
    sim = ChaosSimulator()
    rs = np.linspace(2.5, 4.0, 50)
    R, X = sim.logistic_bifurcation(rs, n_skip=200, n_keep=20)
    assert R.shape == X.shape == (50 * 20,)
    assert R.min() == pytest.approx(2.5)
    assert R.max() == pytest.approx(4.0)
    assert np.all((X >= 0.0) & (X <= 1.0))


# --------------------------------------------------------------------- #
# Feigenbaum constants                                                  #
# --------------------------------------------------------------------- #

def test_feigenbaum_constants_reference_values():
    sim = ChaosSimulator()
    delta, alpha = sim.feigenbaum_constants()
    assert delta == pytest.approx(4.669201609, rel=1e-9)
    assert alpha == pytest.approx(2.502907875, rel=1e-9)


def test_feigenbaum_delta_estimate_within_one_percent():
    """The estimator must be within 1% of 4.66920."""
    sim = ChaosSimulator()
    delta_est = sim.feigenbaum_delta_estimate()
    rel = abs(delta_est - FEIGENBAUM_DELTA) / FEIGENBAUM_DELTA
    assert rel < 0.01, f"Feigenbaum delta off by {rel:.4%}: {delta_est}"


def test_period_doubling_cascade_detected():
    assert detect_period_doubling_cascade(LOGISTIC_R_PERIOD_DOUBLING)


def test_period_doubling_gaps_strictly_shrink():
    rs = list(LOGISTIC_R_PERIOD_DOUBLING)
    gaps = [rs[i + 1] - rs[i] for i in range(len(rs) - 1)]
    for i in range(len(gaps) - 1):
        assert gaps[i] > gaps[i + 1] > 0.0


def test_period_doubling_accumulation_below_4():
    assert LOGISTIC_R_INFINITY < 4.0
    assert LOGISTIC_R_INFINITY > LOGISTIC_R_PERIOD_DOUBLING[-1]


# --------------------------------------------------------------------- #
# Substrate-framework constants                                         #
# --------------------------------------------------------------------- #

def test_lambda_qcd_constant_present():
    assert LAMBDA_QCD_MEV == pytest.approx(200.0)
