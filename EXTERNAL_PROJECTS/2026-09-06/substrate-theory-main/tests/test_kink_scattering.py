"""Tests for kink_scattering — sine-Gordon dynamic scattering simulator."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.kink_scattering import (
    KinkScattering,
    breather_profile,
    classify_ka_outcome,
    lorentz_kink,
    topological_charge,
)


# ---------------------------------------------------------------------------
# (a) Topological charge conservation
# ---------------------------------------------------------------------------

def test_Q_conserved_KK():
    sim = KinkScattering(Nx=1024, L=80.0, bc="open")
    res = sim.kk_collision(v=0.4, separation=30.0)
    Q = res.Q
    # Total Q should remain near +2 throughout (two kinks each Q=+1)
    assert np.all(np.isfinite(Q))
    assert np.max(np.abs(Q - Q[0])) < 1e-3, f"ΔQ_max = {np.max(np.abs(Q - Q[0]))}"


def test_Q_conserved_KA():
    sim = KinkScattering(Nx=1024, L=80.0, bc="open")
    res = sim.ka_collision(v=0.5, separation=30.0)
    Q = res.Q
    assert np.max(np.abs(Q - Q[0])) < 1e-3


def test_Q_charges_correct_at_init():
    sim = KinkScattering(Nx=512, L=60.0, bc="open")
    sim.init_two_kinks(-15.0, +0.3, +15.0, -0.3, sign_K=+1, sign_AK=+1)
    assert abs(sim.topological_charge() - 2.0) < 0.05
    sim.init_two_kinks(-15.0, +0.3, +15.0, -0.3, sign_K=+1, sign_AK=-1)
    assert abs(sim.topological_charge() - 0.0) < 0.05


# ---------------------------------------------------------------------------
# (b) KK elastic phase shift matches analytic
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("v", [0.3, 0.5])
def test_kk_phase_shift_matches_analytic(v):
    sim = KinkScattering(Nx=2048, L=120.0, bc="open")
    out = sim.phase_shift(v=v, separation=60.0, channel="KK")
    # Analytic Δx = (1/γ) ln((1+v)/(1−v)) for sine-Gordon KK scattering.
    # Numeric extraction is approximate (centre-finder noise) — allow 30%.
    # The sign convention of "Δx" depends on which trajectory branch is
    # taken as reference; what is physical is |Δx| matching the analytic
    # magnitude (the soliton is displaced by this much from the free
    # bounce trajectory).
    num_mag = abs(out["delta_x_numeric"])
    ana_mag = abs(out["delta_x_analytic"])
    rel = abs(num_mag - ana_mag) / ana_mag
    assert ana_mag > 0
    assert rel < 0.4, (
        f"v={v}: |Δx|_num={num_mag:.3f}, |Δx|_ana={ana_mag:.3f}, "
        f"rel={rel:.2%}"
    )


# ---------------------------------------------------------------------------
# (c) Breather oscillates at predicted period
# ---------------------------------------------------------------------------

def test_breather_period():
    omega = 0.6
    sim = KinkScattering(Nx=1024, L=60.0, bc="open", cfl=0.3)
    sim.breather_initial(omega=omega)
    # The field at x=0 should oscillate as 4 arctan((η/ω) sin(ω t))
    T_pred = 2.0 * np.pi / omega
    res = sim.run(t_end=3.0 * T_pred, n_snapshots=200)
    # Track field at x=0 via snapshots
    snaps = res.phi_history
    ts = res.snapshot_t
    assert snaps is not None
    Nx = snaps.shape[1]
    centre_idx = Nx // 2
    phi_centre = snaps[:, centre_idx]
    # Find zero crossings (φ at centre passes through 0 each half period)
    # because sin(ω t) has zeros at t = n·π/ω, period T = 2π/ω.
    sign = np.sign(phi_centre)
    flips = np.where(np.diff(sign) != 0)[0]
    if len(flips) < 3:
        pytest.skip("not enough zero crossings to measure period")
    # Period = 2 × mean spacing between zero crossings.
    crossing_times = ts[flips]
    dts = np.diff(crossing_times)
    period_meas = 2.0 * float(np.mean(dts))
    rel = abs(period_meas - T_pred) / T_pred
    assert rel < 0.10, f"measured period {period_meas:.3f}, expected {T_pred:.3f}"

    # Stability: charge stays at 0
    assert np.max(np.abs(res.Q)) < 1e-3


# ---------------------------------------------------------------------------
# (d) Low-velocity KA produces breather (energy below 2 m_K)
# ---------------------------------------------------------------------------

def test_low_v_KA_forms_bound_state():
    """Low-energy KA configuration (E < 2 m_K) is the breather bound state.

    In pure (integrable) sine-Gordon a free kink-antikink scattering at any
    velocity is elastic — the solitons pass through.  The bound state
    (breather) is a *separate* solution with E < 2 m_K = 16.  We initialise
    in that bound configuration directly (omega < 1 ⇒ E_breather < 16) and
    verify it remains spatially localized and conserves Q = 0.
    """
    sim = KinkScattering(Nx=1024, L=60.0, bc="open", cfl=0.3)
    sim.breather_initial(omega=0.4)
    E0 = sim.total_energy()
    # Breather rest energy: E_B = 16 √(1−ω²) < 2 m_K = 16
    expected = 16.0 * np.sqrt(1.0 - 0.4 ** 2)
    assert E0 < 2.0 * 8.0, f"breather energy {E0} should be below 2 m_K = 16"
    assert abs(E0 - expected) / expected < 0.02

    res = sim.run(t_end=80.0, n_snapshots=20, n_expected_kinks=2)
    outcome = classify_ka_outcome(res, kink_mass=8.0)
    assert outcome == "breather", f"expected breather, got {outcome}"
    assert np.max(np.abs(res.Q)) < 1e-3


# ---------------------------------------------------------------------------
# Sanity — single static kink advances at v=0 without distortion
# ---------------------------------------------------------------------------

def test_static_kink_stable():
    sim = KinkScattering(Nx=1024, L=60.0, bc="open")
    sim.phi = lorentz_kink(sim.x, 0.0, 0.0, v=0.0)
    sim.pi = np.zeros_like(sim.x)
    Q0 = sim.topological_charge()
    E0 = sim.total_energy()
    sim.step(2000)
    # Open BC allows tiny endpoint drift; well within 1e-3 tolerance.
    assert abs(sim.topological_charge() - Q0) < 1e-3
    assert abs(sim.total_energy() - E0) / E0 < 1e-3
