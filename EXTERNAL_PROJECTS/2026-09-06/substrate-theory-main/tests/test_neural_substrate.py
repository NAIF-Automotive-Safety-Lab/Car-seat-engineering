"""Tests for the paired neural action-potential GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) AP threshold ≈ -55 mV  (sub-threshold stimulus -> no firing)
    (b) Resting potential ≈ -65 mV
    (c) Peak voltage  ≈ +40 mV  (within reasonable HH-1952 tolerance)
    (d) Absolute refractory period ~ 1 ms; relative refractory ~ 4 ms
    (e) Myelinated conduction velocity ≈ 100 m/s,
        unmyelinated ≈ 1 m/s
    (f) Nernst potentials reproduce textbook E_Na ≈ +60 mV at body T

Plus geometric / electrophysiological sanity:
    * three channel species (Na+, K+, leak) present
    * gating variables stay in [0, 1]
    * action potential overshoots 0 mV when superthreshold
    * single AP, not bursting, for short pulse
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.neural_substrate import (
    ABSOLUTE_REFRACTORY_MS,
    RELATIVE_REFRACTORY_MS,
    F_FARADAY_C_PER_MOL,
    G_K_MAX_DEFAULT,
    G_L_DEFAULT,
    G_NA_MAX_DEFAULT,
    K_IN_MM,
    K_OUT_MM,
    LAMBDA_QCD_MEV,
    NA_IN_MM,
    NA_OUT_MM,
    R_GAS_CONST_J_PER_MOL_K,
    T_BODY_K,
    V_MYELINATED_M_PER_S,
    V_PEAK_TYP,
    V_REST,
    V_THRESHOLD,
    V_UNMYELINATED_M_PER_S,
    IonChannel,
    NeuralGeometry,
    NeuralSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> NeuralGeometry:
    return NeuralGeometry()


@pytest.fixture
def geom_myelinated() -> NeuralGeometry:
    return NeuralGeometry(myelinated=True)


@pytest.fixture
def sim() -> NeuralSimulator:
    return NeuralSimulator()


# --------------------------------------------------------------------- #
# (a) Threshold is around -55 mV                                        #
# --------------------------------------------------------------------- #

def test_threshold_constant_is_minus_55_mV():
    assert V_THRESHOLD == pytest.approx(-55.0, abs=1e-9)


def test_subthreshold_stimulus_does_not_fire(sim):
    """A tiny stimulus should depolarize but NOT fire an AP (V stays << 0)."""
    out = sim.simulate(I_inj_uA_per_cm2=1.0, stim_duration_ms=0.5)
    # Voltage must NOT cross 0 mV (no AP)
    assert sim.peak_voltage_mV(out) < 0.0
    assert not sim.fired_action_potential(out, V_peak_min=0.0)


def test_superthreshold_stimulus_fires_ap(sim):
    """A larger stimulus crosses threshold and overshoots 0 mV."""
    out = sim.simulate(I_inj_uA_per_cm2=15.0, stim_duration_ms=1.0)
    assert sim.crossed_threshold(out, V_th=V_THRESHOLD)
    assert sim.fired_action_potential(out, V_peak_min=0.0)


# --------------------------------------------------------------------- #
# (b) Resting potential ≈ -65 mV                                        #
# --------------------------------------------------------------------- #

def test_resting_potential_constant_is_minus_65_mV():
    assert V_REST == pytest.approx(-65.0, abs=1e-9)


def test_resting_state_is_stable(sim):
    """With zero stimulus, V stays near V_rest for the whole window."""
    out = sim.simulate(I_inj_uA_per_cm2=0.0, stim_duration_ms=0.0)
    # Membrane should not depolarize on its own
    assert np.max(out["V_mV"]) < V_REST + 1.0
    assert np.min(out["V_mV"]) > V_REST - 1.0


def test_resting_gating_steady_state(sim):
    """At V_rest, gating variables are all between 0 and 1, h is dominant."""
    V0, m0, h0, n0 = sim.resting_state()
    assert V0 == pytest.approx(V_REST, abs=1e-9)
    for x in (m0, h0, n0):
        assert 0.0 <= x <= 1.0
    # At V_rest, h ~ 0.6 (Na+ inactivation gate is mostly OPEN), m << h
    assert h0 > m0
    assert h0 > 0.4
    assert m0 < 0.1


# --------------------------------------------------------------------- #
# (c) Peak voltage ≈ +40 mV                                             #
# --------------------------------------------------------------------- #

def test_peak_voltage_close_to_plus_40_mV(sim):
    """The AP peak should land in the +30..+55 mV window."""
    out = sim.simulate(I_inj_uA_per_cm2=15.0, stim_duration_ms=1.0)
    peak = sim.peak_voltage_mV(out)
    # HH 1952 with these conductances yields ~ +40 mV; allow generous band
    assert peak == pytest.approx(V_PEAK_TYP, abs=15.0), (
        f"peak={peak:.1f} mV vs target {V_PEAK_TYP} mV"
    )
    # Sanity: peak is below E_Na (cannot exceed Na+ reversal potential)
    geom = NeuralGeometry()
    assert peak < geom.E_Na() + 1.0


def test_single_ap_for_short_pulse(sim):
    """A single short pulse should produce exactly one spike (no bursting)."""
    out = sim.simulate(I_inj_uA_per_cm2=15.0, stim_duration_ms=1.0)
    spikes = sim.spike_times_ms(out, V_threshold=0.0)
    assert len(spikes) == 1


# --------------------------------------------------------------------- #
# (d) Refractory period                                                 #
# --------------------------------------------------------------------- #

def test_absolute_refractory_constant_is_about_1_ms():
    assert 0.5 <= ABSOLUTE_REFRACTORY_MS <= 2.0


def test_relative_refractory_constant_is_about_4_ms():
    assert 2.0 <= RELATIVE_REFRACTORY_MS <= 6.0


def test_two_pulses_within_absolute_refractory_no_second_ap(sim):
    """A second pulse 1 ms after the first should fail to fire (absolute refractory)."""
    # Use longer simulation to capture both
    sim2 = NeuralSimulator(duration_ms=60.0)
    out = sim2.simulate(
        I_inj_uA_per_cm2=15.0,
        stim_start_ms=5.0,
        stim_duration_ms=1.0,
        second_stim_ms=6.0,         # only 1 ms after first stimulus start
        second_stim_duration_ms=1.0,
    )
    spikes = sim2.spike_times_ms(out, V_threshold=0.0)
    # During absolute refractory, only one AP fires
    assert len(spikes) == 1


def test_two_pulses_well_separated_fire_two_aps(sim):
    """Two pulses 20 ms apart should fire two APs (well past refractory)."""
    sim2 = NeuralSimulator(duration_ms=60.0)
    out = sim2.simulate(
        I_inj_uA_per_cm2=15.0,
        stim_start_ms=5.0,
        stim_duration_ms=1.0,
        second_stim_ms=25.0,
        second_stim_duration_ms=1.0,
    )
    spikes = sim2.spike_times_ms(out, V_threshold=0.0)
    assert len(spikes) == 2


def test_measure_refractory_returns_finite_min_gap():
    """Refractory probe should report some finite minimum-fireable gap."""
    sim2 = NeuralSimulator(duration_ms=60.0)
    out = sim2.measure_refractory_ms(
        I_inj_uA_per_cm2=15.0,
        stim_duration_ms=1.0,
        gap_ms=(1.0, 5.0, 10.0, 15.0, 20.0),
    )
    assert out["min_fireable_gap_ms"] < float("inf")
    # Smallest gap that fires must be larger than absolute refractory
    assert out["min_fireable_gap_ms"] >= 1.0


# --------------------------------------------------------------------- #
# (e) Conduction velocities                                             #
# --------------------------------------------------------------------- #

def test_unmyelinated_conduction_velocity(geom):
    assert geom.conduction_velocity_m_per_s() == pytest.approx(
        V_UNMYELINATED_M_PER_S, abs=1e-9
    )


def test_myelinated_conduction_velocity(geom_myelinated):
    assert geom_myelinated.conduction_velocity_m_per_s() == pytest.approx(
        V_MYELINATED_M_PER_S, abs=1e-9
    )


def test_myelin_speedup_ratio(sim):
    speedup = sim.conduction_speedup_myelin()
    assert speedup == pytest.approx(100.0, abs=1e-9)


# --------------------------------------------------------------------- #
# (f) Nernst potentials                                                 #
# --------------------------------------------------------------------- #

def test_nernst_potential_Na_around_plus_60_mV(geom):
    """E_Na = (RT/zF) ln(145/15) at body T should be ~ +60 mV."""
    E = geom.E_Na_nernst()
    # Textbook value ~ +60..+70 mV at 37 deg C
    assert 55.0 < E < 75.0


def test_nernst_potential_K_around_minus_95_mV(geom):
    """E_K = (RT/zF) ln(4/140) ~ -95 mV at body T."""
    E = geom.E_K_nernst()
    assert -100.0 < E < -85.0


def test_hh_reversal_potentials_are_calibrated_values(geom):
    """Simulator uses HH 1952 calibrated reversals (+50, -77, -54.4) mV."""
    assert geom.E_Na() == pytest.approx(50.0, abs=1e-9)
    assert geom.E_K()  == pytest.approx(-77.0, abs=1e-9)
    assert geom.E_L()  == pytest.approx(-54.4, abs=1e-9)


def test_nernst_potential_chloride_negative(geom):
    """E_Cl with mammalian concentrations is negative."""
    assert geom.E_Cl() < -50.0


def test_nernst_zero_when_concentrations_equal():
    """Nernst potential vanishes when [in] = [out]."""
    E = NeuralGeometry.nernst_potential_mV(100.0, 100.0, +1)
    assert E == pytest.approx(0.0, abs=1e-9)


def test_nernst_inverts_when_valence_flips():
    """E_anion = -E_cation for the same concentration ratio."""
    E_pos = NeuralGeometry.nernst_potential_mV(10.0, 100.0, +1)
    E_neg = NeuralGeometry.nernst_potential_mV(10.0, 100.0, -1)
    assert E_pos == pytest.approx(-E_neg, abs=1e-9)


# --------------------------------------------------------------------- #
# Channel inventory                                                     #
# --------------------------------------------------------------------- #

def test_three_channel_species(geom):
    """Geometry exposes Na+, K+, leak."""
    chs = geom.channels()
    assert len(chs) == 3
    names = {c.name for c in chs}
    assert names == {"Na+", "K+", "leak"}


def test_two_voltage_gated_one_leak(geom):
    chs = geom.channels()
    n_vgated = sum(1 for c in chs if c.voltage_gated)
    n_leak   = sum(1 for c in chs if not c.voltage_gated)
    assert n_vgated == 2
    assert n_leak == 1


def test_default_max_conductances(geom):
    assert geom.g_Na_max == pytest.approx(G_NA_MAX_DEFAULT)
    assert geom.g_K_max  == pytest.approx(G_K_MAX_DEFAULT)
    assert geom.g_L      == pytest.approx(G_L_DEFAULT)


# --------------------------------------------------------------------- #
# Gating-variable invariants                                            #
# --------------------------------------------------------------------- #

def test_gating_variables_in_unit_interval(sim):
    out = sim.simulate(I_inj_uA_per_cm2=15.0, stim_duration_ms=1.0)
    for key in ("m", "h", "n"):
        assert np.all(out[key] >= 0.0)
        assert np.all(out[key] <= 1.0)


def test_h_gate_inactivates_during_spike(sim):
    """h drops sharply during the AP (Na+ inactivation)."""
    out = sim.simulate(I_inj_uA_per_cm2=15.0, stim_duration_ms=1.0)
    h0 = out["h"][0]
    h_min = float(np.min(out["h"]))
    assert h_min < h0 - 0.2   # h drops by at least 0.2


def test_n_gate_activates_during_spike(sim):
    """n rises sharply during the AP (delayed-rectifier K+ activation)."""
    out = sim.simulate(I_inj_uA_per_cm2=15.0, stim_duration_ms=1.0)
    n0 = out["n"][0]
    n_max = float(np.max(out["n"]))
    assert n_max > n0 + 0.2


def test_m_gate_activates_first(sim):
    """m peaks before n (fast Na+ activation precedes delayed K+ activation)."""
    out = sim.simulate(I_inj_uA_per_cm2=15.0, stim_duration_ms=1.0)
    t_m_peak = out["t_ms"][int(np.argmax(out["m"]))]
    t_n_peak = out["t_ms"][int(np.argmax(out["n"]))]
    assert t_m_peak < t_n_peak


# --------------------------------------------------------------------- #
# Substrate signature                                                   #
# --------------------------------------------------------------------- #

def test_geometry_substrate_signature(geom):
    sig = geom.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["channel_species"] == 3
    assert sig["voltage_gated"] == 2
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV


def test_simulator_substrate_signature(sim):
    sig = sim.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["myelin_speedup"] == pytest.approx(100.0, abs=1e-9)
    assert "AP" in sig["interpretation"]
    assert "refractory" in sig["interpretation"]
