"""Tests for spectroscopy_substrate.

All predicted peaks must lie within 5% of the literature/handbook
reference values.  In practice the closed-form models reproduce the
references to better than 1%.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from stiff_medium.spectroscopy_substrate import (
    BOND_FORCE_CONSTANT,
    IR_REF,
    NMR_REF,
    RAMAN_REF,
    UV_VIS_REF,
    Molecule,
    SpectroscopyGeometry,
    SpectroscopySimulator,
    make_simulator,
    reduced_mass,
)


# --------------------------------------------------------------------------- #
# fixtures                                                                    #
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def sim() -> SpectroscopySimulator:
    return make_simulator()


@pytest.fixture(scope="module")
def geom() -> SpectroscopyGeometry:
    return SpectroscopyGeometry()


# --------------------------------------------------------------------------- #
# helpers                                                                     #
# --------------------------------------------------------------------------- #

def _within_pct(calc: float, ref: float, pct: float) -> bool:
    return abs(calc - ref) / abs(ref) * 100.0 < pct


# --------------------------------------------------------------------------- #
# 1. GEOMETRY                                                                 #
# --------------------------------------------------------------------------- #

def test_geometry_has_required_molecules(geom: SpectroscopyGeometry) -> None:
    for name in ("benzene", "beta_carotene", "water", "methanol",
                 "ethanol", "acetone"):
        assert name in geom.molecules
        assert isinstance(geom.molecule(name), Molecule)


def test_geometry_vibrational_modes_unique(geom: SpectroscopyGeometry) -> None:
    modes = geom.vibrational_modes("ethanol")
    assert len(modes) == len(set(modes)), modes
    assert "C-H_sp3" in modes


def test_geometry_box_length_positive(geom: SpectroscopyGeometry) -> None:
    assert geom.electronic_box_length_m("benzene") > 0


def test_geometry_nuclear_spin_environments(geom: SpectroscopyGeometry) -> None:
    envs = geom.nuclear_spin_environments("ethanol")
    labels = {e[0] for e in envs}
    assert labels == {"CH3", "CH2", "OH"}
    integration_total = sum(e[2] for e in envs)
    assert integration_total == pytest.approx(6.0)  # 3+2+1


# --------------------------------------------------------------------------- #
# 2. UV-Vis HOMO->LUMO                                                        #
# --------------------------------------------------------------------------- #

def test_uvvis_benzene_within_5pct(sim: SpectroscopySimulator) -> None:
    lam = sim.uv_vis_lambda_nm("benzene")
    assert _within_pct(lam, UV_VIS_REF["benzene"], 5.0), lam


def test_uvvis_beta_carotene_within_5pct(sim: SpectroscopySimulator) -> None:
    lam = sim.uv_vis_lambda_nm("beta_carotene")
    assert _within_pct(lam, UV_VIS_REF["beta_carotene"], 5.0), lam


def test_uvvis_curve_peaks_at_lambda_max(sim: SpectroscopySimulator) -> None:
    waves = np.linspace(150.0, 600.0, 5000)
    curve = sim.uv_vis_curve("benzene", waves)
    lam_max_pred = waves[int(np.argmax(curve))]
    assert _within_pct(lam_max_pred, UV_VIS_REF["benzene"], 5.0)


# --------------------------------------------------------------------------- #
# 3. IR vibrational                                                           #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("bond,ref_key", [
    ("C=O",     "C=O_stretch_ketone"),
    ("O-H",     "O-H_stretch_alcohol_centre"),
    ("C-H_sp3", "C-H_stretch_sp3"),
    ("C-H_sp2", "C-H_stretch_aromatic"),
    ("C=C",     "C=C_stretch_aromatic"),
    ("C-O",     "C-O_stretch_alcohol"),
])
def test_ir_bond_within_5pct(
    sim: SpectroscopySimulator, bond: str, ref_key: str
) -> None:
    nu = sim.ir_wavenumber_cm(bond)
    assert _within_pct(nu, IR_REF[ref_key], 5.0), (bond, nu)


def test_ir_co_in_carbonyl_window(sim: SpectroscopySimulator) -> None:
    nu = sim.ir_wavenumber_cm("C=O")
    assert 1680.0 <= nu <= 1750.0, nu  # canonical C=O window


def test_ir_oh_in_alcohol_envelope(sim: SpectroscopySimulator) -> None:
    nu = sim.ir_wavenumber_cm("O-H")
    assert 3200.0 <= nu <= 3600.0, nu


def test_ir_ch_centre_in_2900_window(sim: SpectroscopySimulator) -> None:
    nu = sim.ir_wavenumber_cm("C-H_sp3")
    assert 2850.0 <= nu <= 3000.0, nu


def test_ir_curve_has_oh_and_ch_features(sim: SpectroscopySimulator) -> None:
    nus = np.linspace(500.0, 4000.0, 4000)
    curve = sim.ir_curve("ethanol", nus)
    # Find local maxima as features (any sample greater than both neighbours)
    # and confirm O-H around 3400 and C-H around 2960 are present.
    oh_window = (nus > 3200.0) & (nus < 3600.0)
    ch_window = (nus > 2800.0) & (nus < 3050.0)
    assert curve[oh_window].max() > 0.3
    assert curve[ch_window].max() > 0.3


# --------------------------------------------------------------------------- #
# 4. NMR ^1H                                                                  #
# --------------------------------------------------------------------------- #

def test_nmr_ethanol_three_peaks(sim: SpectroscopySimulator) -> None:
    peaks = sim.nmr_peaks("ethanol")
    assert len(peaks) == 3
    labels = {p[0] for p in peaks}
    assert labels == {"CH3", "CH2", "OH"}


@pytest.mark.parametrize("env,key", [
    ("CH3", "CH3_delta_ppm"),
    ("CH2", "CH2_delta_ppm"),
    ("OH",  "OH_delta_ppm"),
])
def test_nmr_ethanol_shifts_within_5pct(
    sim: SpectroscopySimulator, env: str, key: str
) -> None:
    peaks = {p[0]: p for p in sim.nmr_peaks("ethanol")}
    delta = peaks[env][1]
    assert _within_pct(delta, NMR_REF["ethanol"][key], 5.0), (env, delta)


def test_nmr_ethanol_integration_ratios(sim: SpectroscopySimulator) -> None:
    peaks = {p[0]: p[2] for p in sim.nmr_peaks("ethanol")}
    assert peaks["CH3"] == pytest.approx(3.0)
    assert peaks["CH2"] == pytest.approx(2.0)
    assert peaks["OH"]  == pytest.approx(1.0)


def test_nmr_split_pattern_ethanol(sim: SpectroscopySimulator) -> None:
    pat = sim.nmr_split_pattern("ethanol")
    assert pat["CH3"].startswith("triplet")
    assert pat["CH2"].startswith("quartet")
    assert pat["OH"].startswith("singlet")


def test_nmr_curve_integration_monotone(sim: SpectroscopySimulator) -> None:
    ppm = np.linspace(0.0, 10.0, 5000)
    spec, cum = sim.nmr_curve("ethanol", ppm)
    # cumulative integration should be monotone increasing along increasing ppm
    order = np.argsort(ppm)
    cum_sorted = cum[order]
    assert np.all(np.diff(cum_sorted) >= -1e-10)
    assert spec.max() > 0


# --------------------------------------------------------------------------- #
# 5. Raman                                                                    #
# --------------------------------------------------------------------------- #

def test_raman_benzene_ring_breathing_close(sim: SpectroscopySimulator) -> None:
    modes = dict(sim.raman_modes("benzene"))
    assert "ring_breathing" in modes
    assert _within_pct(modes["ring_breathing"],
                       RAMAN_REF["benzene_ring_breathing"], 5.0)


def test_raman_acetone_co_close_to_ir(sim: SpectroscopySimulator) -> None:
    """Raman C=O of acetone should be close to its IR C=O (~1715 cm^-1)."""
    modes = dict(sim.raman_modes("acetone"))
    assert "C=O_stretch" in modes
    assert _within_pct(modes["C=O_stretch"], 1715.0, 5.0)


# --------------------------------------------------------------------------- #
# 6. Comparison & report                                                      #
# --------------------------------------------------------------------------- #

def test_compare_returns_full_table(sim: SpectroscopySimulator) -> None:
    tbl = sim.compare()
    for v in tbl.values():
        assert {"calc", "ref", "abs_err", "rel_err_pct"} <= set(v)


def test_compare_all_within_5pct(sim: SpectroscopySimulator) -> None:
    """Every predicted peak in compare() must be within 5% of reference."""
    tbl = sim.compare()
    bad = {k: v for k, v in tbl.items() if v["rel_err_pct"] > 5.0}
    assert not bad, bad


def test_report_string_has_summary(sim: SpectroscopySimulator) -> None:
    s = sim.report()
    assert "Substrate Molecular Spectroscopy" in s
    assert "worst" in s


# --------------------------------------------------------------------------- #
# 7. Sanity / error paths                                                     #
# --------------------------------------------------------------------------- #

def test_unknown_molecule_raises(sim: SpectroscopySimulator) -> None:
    with pytest.raises(KeyError):
        sim.geometry.molecule("phlogiston")


def test_no_chromophore_raises(sim: SpectroscopySimulator) -> None:
    # acetone has no n_pi_electrons set for HOMO->LUMO box model
    with pytest.raises(ValueError):
        sim.uv_vis_lambda_nm("acetone")


def test_reduced_mass_symmetric() -> None:
    mu_ab = reduced_mass("C", "O")
    mu_ba = reduced_mass("O", "C")
    assert mu_ab == pytest.approx(mu_ba)


def test_force_constant_table_complete() -> None:
    for label in ("C=O", "O-H", "C-H_sp3", "C=C", "C-O"):
        assert label in BOND_FORCE_CONSTANT
        assert BOND_FORCE_CONSTANT[label] > 0
