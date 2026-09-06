"""Tests for photosynthesis_substrate.

Covers:
  - PhotosynthesisGeometry: chlorophyll a/b absorption peaks within 5%
    of the textbook values 430+662 / 453+642 nm.
  - Reaction-center potentials P680 and P700 are positive (oxidizing).
  - Excited-state RCs are strong reductants (negative E_h).
  - Z-scheme: PS-II to PS-I electron flow goes downhill in E_h between
    the two photo-driven jumps.
  - Quantum efficiency of light harvesting ~95%.
  - 8-photon stoichiometry per O2 (4 PS-II + 4 PS-I).
  - Calvin cycle stoichiometry: 3 ATP + 2 NADPH per CO2;
    18 ATP + 12 NADPH + 6 CO2 per hexose.
  - Photon energy hc/lambda is correctly computed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from stiff_medium.photosynthesis_substrate import (
    CHLOROPHYLL_PEAKS_NM,
    PhotosynthesisGeometry,
    PhotosynthesisSimulator,
    Z_SCHEME_REDOX_eV,
    PHOTONS_PER_O2,
    ATP_PER_HEXOSE,
    NADPH_PER_HEXOSE,
    CO2_PER_HEXOSE,
    ATP_PER_CO2,
    NADPH_PER_CO2,
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def sim():
    s = PhotosynthesisSimulator()
    s.compare_to_reference()
    return s


def _close_pct(a, b, pct):
    return abs(a - b) / abs(b) * 100.0 < pct


# ---------------------------------------------------------------------------
# 1. PhotosynthesisGeometry: chlorophyll peaks within 5%
# ---------------------------------------------------------------------------

def test_chl_a_soret_peak_within_5pct():
    g = PhotosynthesisGeometry(chlorophyll_type="chl_a")
    peaks = g.absorption_peaks_nm()
    assert _close_pct(peaks["soret"], 430.0, 5.0)


def test_chl_a_qy_peak_within_5pct():
    g = PhotosynthesisGeometry(chlorophyll_type="chl_a")
    peaks = g.absorption_peaks_nm()
    assert _close_pct(peaks["qy"], 662.0, 5.0)


def test_chl_b_soret_peak_within_5pct():
    g = PhotosynthesisGeometry(chlorophyll_type="chl_b")
    peaks = g.absorption_peaks_nm()
    assert _close_pct(peaks["soret"], 453.0, 5.0)


def test_chl_b_qy_peak_within_5pct():
    g = PhotosynthesisGeometry(chlorophyll_type="chl_b")
    peaks = g.absorption_peaks_nm()
    assert _close_pct(peaks["qy"], 642.0, 5.0)


def test_chl_a_qy_red_of_soret():
    """Qy band is in the red, Soret in the blue — Qy > Soret in nm."""
    g = PhotosynthesisGeometry(chlorophyll_type="chl_a")
    p = g.absorption_peaks_nm()
    assert p["qy"] > p["soret"]


def test_unknown_chlorophyll_raises():
    g = PhotosynthesisGeometry(chlorophyll_type="chl_X")
    with pytest.raises(ValueError):
        g.absorption_peaks_nm()


def test_qy_band_energy_is_red_photon_eV():
    """Qy at 662 nm = ~1.87 eV (red photon)."""
    g = PhotosynthesisGeometry(chlorophyll_type="chl_a")
    e_qy = g.absorption_peak_eV("qy")
    assert 1.7 < e_qy < 2.0


def test_soret_band_energy_is_blue_photon_eV():
    """Soret at 430 nm = ~2.88 eV (blue photon)."""
    g = PhotosynthesisGeometry(chlorophyll_type="chl_a")
    e_soret = g.absorption_peak_eV("soret")
    assert 2.7 < e_soret < 3.1


def test_absorption_band_has_two_peaks():
    """Absorption spectrum has clear maxima near both Soret and Qy bands."""
    g = PhotosynthesisGeometry(chlorophyll_type="chl_a")
    wl = np.linspace(350.0, 750.0, 401)
    spectrum = g.absorption_band(wl)
    # peak detection: Soret at ~430, Qy at ~662
    soret_idx = np.argmax(spectrum[wl < 500])
    qy_window = (wl > 600) & (wl < 720)
    qy_idx_local = np.argmax(spectrum[qy_window])
    qy_wl = wl[qy_window][qy_idx_local]
    soret_wl = wl[wl < 500][soret_idx]
    assert _close_pct(soret_wl, 430.0, 5.0)
    assert _close_pct(qy_wl, 662.0, 5.0)


# ---------------------------------------------------------------------------
# 2. Reaction-center potentials
# ---------------------------------------------------------------------------

def test_p680_is_strong_oxidant():
    """P680 ground state is a strong oxidant (E_h > +1 V) to split water."""
    g = PhotosynthesisGeometry(reaction_center="P680")
    e = g.reaction_center_potential()
    assert e > +1.0


def test_p700_is_oxidant():
    """P700 ground state is a moderate oxidant (E_h > 0)."""
    g = PhotosynthesisGeometry(reaction_center="P700")
    e = g.reaction_center_potential()
    assert e > 0.0


def test_p680_excited_is_reductant():
    """P680* (after photon) is a reductant (E_h < 0)."""
    g = PhotosynthesisGeometry(reaction_center="P680")
    e = g.reaction_center_excited_potential()
    assert e < 0.0


def test_p700_excited_is_strongest_reductant():
    """P700* is the strongest reductant in the chain (most negative E_h)."""
    g = PhotosynthesisGeometry(reaction_center="P700")
    e_p700_ex = g.reaction_center_excited_potential()
    g2 = PhotosynthesisGeometry(reaction_center="P680")
    e_p680_ex = g2.reaction_center_excited_potential()
    assert e_p700_ex < e_p680_ex


def test_photon_excitation_makes_rc_more_reducing():
    """Photon absorption: E_h(P680*) < E_h(P680)."""
    g = PhotosynthesisGeometry(reaction_center="P680")
    assert g.reaction_center_excited_potential() < g.reaction_center_potential()


# ---------------------------------------------------------------------------
# 3. Z-scheme structure
# ---------------------------------------------------------------------------

def test_z_scheme_chain_starts_at_water():
    chain = PhotosynthesisSimulator.z_scheme()
    assert chain[0][0] == "H2O/O2"


def test_z_scheme_chain_ends_at_NADPH():
    chain = PhotosynthesisSimulator.z_scheme()
    assert chain[-1][0] == "NADPH"


def test_z_scheme_has_two_photon_jumps():
    """Two upward jumps where E_h goes more negative (PS-II and PS-I)."""
    chain = PhotosynthesisSimulator.z_scheme()
    energies = [e for _, e in chain]
    # find positions where E_h drops abruptly (more negative): photon pumps
    photon_pumps = []
    for i in range(1, len(energies)):
        if energies[i] - energies[i - 1] < -1.0:
            photon_pumps.append(i)
    assert len(photon_pumps) == 2, f"expected 2 photon jumps, got {photon_pumps}"


def test_z_scheme_overall_drop_water_to_nadph_negative():
    """The full chain lifts e- from H2O (+0.82) to NADPH (-0.32) - photons drive this."""
    chain = PhotosynthesisSimulator.z_scheme()
    e_first = chain[0][1]
    e_last = chain[-1][1]
    assert e_last < e_first
    delta = e_first - e_last
    assert delta > 1.0   # > 1 V driven uphill by 2 photons


# ---------------------------------------------------------------------------
# 4. Quantum efficiency ~95%
# ---------------------------------------------------------------------------

def test_quantum_yield_is_about_95pct(sim):
    qy = sim.light_harvest_efficiency()
    assert _close_pct(qy, 0.95, 1.0)


def test_quantum_yield_in_valid_range(sim):
    qy = sim.light_harvest_efficiency()
    assert 0.0 <= qy <= 1.0


def test_n_excitations_to_rc_is_high():
    """100 photons * 95% = 95 excitations to RC."""
    s = PhotosynthesisSimulator()
    n_excited = s.n_excitations_to_reach_rc(100)
    assert _close_pct(n_excited, 95.0, 1.0)


# ---------------------------------------------------------------------------
# 5. 8-photon stoichiometry per O2
# ---------------------------------------------------------------------------

def test_eight_photons_per_o2(sim):
    """Z-scheme: 4 e- per O2 * 2 photosystems = 8 photons / O2."""
    n = sim.photons_for_o2()
    assert n == 8


def test_one_o2_per_eight_photons(sim):
    o2 = sim.o2_per_n_photons(8)
    assert _close_pct(o2, 1.0, 1e-6)


def test_eight_photons_minimum_quantum_requirement():
    """The well-known 'minimum quantum requirement' of photosynthesis = 8."""
    assert PHOTONS_PER_O2 == 8


def test_stoichiometry_table_per_o2(sim):
    s = sim.stoichiometry_table()
    assert s["photons"] == 8.0
    assert s["NADPH"] == 2.0
    assert s["O2"] == 1.0
    assert s["H2O_split"] == 2.0


# ---------------------------------------------------------------------------
# 6. Calvin cycle stoichiometry
# ---------------------------------------------------------------------------

def test_calvin_per_co2_three_atp_two_nadph(sim):
    d = sim.calvin_cycle_demand_per_co2()
    assert d["ATP"] == pytest.approx(3.0, abs=1e-9)
    assert d["NADPH"] == pytest.approx(2.0, abs=1e-9)


def test_calvin_per_hexose_18_atp_12_nadph(sim):
    d = sim.calvin_cycle_demand_per_hexose()
    assert d["ATP"] == pytest.approx(18.0, abs=1e-9)
    assert d["NADPH"] == pytest.approx(12.0, abs=1e-9)
    assert d["CO2"] == pytest.approx(6.0, abs=1e-9)


def test_calvin_demand_scales_linearly():
    """Hexose demand = 6 * per-CO2 demand."""
    s = PhotosynthesisSimulator()
    per_co2 = s.calvin_cycle_demand_per_co2()
    per_hex = s.calvin_cycle_demand_per_hexose()
    assert per_hex["ATP"] == pytest.approx(6 * per_co2["ATP"], abs=1e-9)
    assert per_hex["NADPH"] == pytest.approx(6 * per_co2["NADPH"], abs=1e-9)


def test_photons_per_hexose_is_48(sim):
    """6 CO2/hexose * 8 photons/O2 = 48 photons/hexose (since 1 O2 per CO2)."""
    n = sim.photons_per_hexose()
    assert n == 48


def test_calvin_module_constants():
    """Module-level constants encode the textbook stoichiometry."""
    assert ATP_PER_CO2 == 3.0
    assert NADPH_PER_CO2 == 2.0
    assert CO2_PER_HEXOSE == 6
    assert ATP_PER_HEXOSE == 18
    assert NADPH_PER_HEXOSE == 12


# ---------------------------------------------------------------------------
# 7. Photon energy hc/lambda
# ---------------------------------------------------------------------------

def test_photon_energy_red_680nm(sim):
    """680 nm photon = 1.823 eV."""
    e = sim.ps2_photon_energy_eV()
    assert _close_pct(e, 1.823, 1.0)


def test_photon_energy_red_700nm(sim):
    """700 nm photon = 1.771 eV (slightly less than 680nm)."""
    e = sim.ps1_photon_energy_eV()
    assert _close_pct(e, 1.771, 1.0)
    assert e < sim.ps2_photon_energy_eV()


def test_photon_energy_inverse_wavelength():
    """E ~ 1/lambda: E(400) > E(700)."""
    s = PhotosynthesisSimulator()
    assert s.photon_energy_eV(400.0) > s.photon_energy_eV(700.0)


# ---------------------------------------------------------------------------
# 8. ATP synthesis from proton gradient
# ---------------------------------------------------------------------------

def test_atp_per_proton_about_quarter():
    """Chloroplast c14 ATP synthase: ~0.24 ATP/H+ (14/3 ~ 4.67 H+/ATP)."""
    s = PhotosynthesisSimulator()
    atp = s.atp_from_protons(12.0)
    # 12 H+ / ~4.14 H+/ATP ~ 2.9 ATP
    assert 2.5 < atp < 3.3


def test_pmf_from_pH_gradient():
    """PMF ~ 0.18 V from delta_pH = 3 (light-driven gradient)."""
    s = PhotosynthesisSimulator()
    pmf = s.proton_motive_force_eV(delta_pH=3.0)
    assert _close_pct(pmf, 0.177, 5.0)


def test_pmf_includes_membrane_potential():
    """Adding 100 mV delta_psi increases PMF."""
    s = PhotosynthesisSimulator()
    pmf_no_psi = s.proton_motive_force_eV(delta_pH=3.0)
    pmf_with_psi = s.proton_motive_force_eV(delta_pH=3.0, delta_psi_mV=100.0)
    assert pmf_with_psi > pmf_no_psi


# ---------------------------------------------------------------------------
# 9. Comparison & report
# ---------------------------------------------------------------------------

def test_compare_table_keys(sim):
    table = sim.results["_compare"]
    for k in ("chl_a Soret nm", "chl_a Qy nm",
              "chl_b Soret nm", "chl_b Qy nm",
              "photons / O2", "NADPH / O2",
              "light_harvest_qy"):
        assert k in table


def test_compare_chl_peaks_within_5pct(sim):
    """All four chlorophyll absorption peaks must be within 5% of reference."""
    table = sim.results["_compare"]
    for k in ("chl_a Soret nm", "chl_a Qy nm",
              "chl_b Soret nm", "chl_b Qy nm"):
        assert abs(table[k]["rel_err_pct"]) < 5.0, (k, table[k])


def test_report_runs():
    s = PhotosynthesisSimulator()
    txt = s.report()
    assert "Chlorophyll" in txt
    assert "Z-scheme" in txt
    assert "Calvin" in txt
    assert "NADPH" in txt
