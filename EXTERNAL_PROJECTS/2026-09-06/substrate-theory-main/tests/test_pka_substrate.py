"""Tests for the paired pKa GEOMETRY + SIMULATOR module.

Covers the four requested claims:
    (1) PKaGeometry: dissociation = bond breaking + shell release
    (2) PKaSimulator: predict pKa for the named acids within 1 unit
    (3) Henderson-Hasselbalch buffer pH
    (4) Titration curves (strong and weak acid)

The "8 common pKa values within 1 unit" requirement is enforced
via parametrize over the canonical eight reference acids.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.pka_substrate import (
    PKaGeometry,
    PKaSimulator,
    REFERENCE_ACIDS,
    R_GAS,
    T_STD,
    LN10,
    PKW,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def sim():
    return PKaSimulator()


@pytest.fixture
def acetic_geom():
    e = REFERENCE_ACIDS["Acetic acid"]
    return PKaGeometry(
        name="Acetic acid",
        D_XH=e["D_XH"],
        E_solv=e["E_solv"],
        S_conf=e.get("S_conf", 0.0),
    )


# --------------------------------------------------------------------- #
# (1) PKaGeometry: dissociation = bond breaking + shell release         #
# --------------------------------------------------------------------- #

def test_geometry_bridge_strain_positive(acetic_geom):
    """Bridge strain energy is the bond enthalpy."""
    assert acetic_geom.bridge_strain_energy() == pytest.approx(469.0)
    assert acetic_geom.bridge_strain_energy() > 0.0


def test_geometry_shell_release_positive(acetic_geom):
    """Shell release energy is positive (favors dissociation)."""
    assert acetic_geom.shell_release_energy() > 0.0


def test_geometry_delta_g_dissociation_components(acetic_geom):
    """DeltaG = D_XH - E_solv - S_conf."""
    dg = acetic_geom.delta_g_dissociation()
    expected = (acetic_geom.D_XH
                - acetic_geom.E_solv
                - acetic_geom.S_conf)
    assert dg == pytest.approx(expected)


def test_geometry_classifies_strong_acid():
    """HCl: delta_G < 0 -> saturated regime (strong acid)."""
    e = REFERENCE_ACIDS["HCl"]
    geom = PKaGeometry(name="HCl", D_XH=e["D_XH"],
                       E_solv=e["E_solv"], S_conf=e["S_conf"])
    assert geom.delta_g_dissociation() < 0.0
    assert "saturated" in geom.saturation_regime().lower()


def test_geometry_classifies_water_as_weak():
    """Water: delta_G ~ 90 kJ/mol -> barely-acid regime."""
    e = REFERENCE_ACIDS["Water"]
    geom = PKaGeometry(name="Water", D_XH=e["D_XH"],
                       E_solv=e["E_solv"], S_conf=e["S_conf"])
    dg = geom.delta_g_dissociation()
    assert dg > 70.0
    assert "barely" in geom.saturation_regime().lower()


def test_geometry_signature_has_required_keys(acetic_geom):
    sig = acetic_geom.substrate_signature()
    for k in ("interpretation",
              "bridge_strain_kJ_per_mol",
              "shell_release_kJ_per_mol",
              "delta_g_diss_kJ_per_mol",
              "regime"):
        assert k in sig


# --------------------------------------------------------------------- #
# (2) PKaSimulator: predict pKa within 1 unit for 8 common acids        #
# --------------------------------------------------------------------- #

CANONICAL_EIGHT = [
    "HCl",
    "HF",
    "HNO3",
    "Acetic acid",
    "Formic acid",
    "Phosphoric",
    "Ammonium",
    "Water",
]


@pytest.mark.parametrize("name", CANONICAL_EIGHT)
def test_predicted_pka_within_1_unit(sim, name):
    pred = sim.predict_pka(name)
    ref = sim.reference_pka(name)
    assert abs(pred - ref) < 1.0, (
        f"{name}: pred={pred:.2f}, ref={ref:.2f}, |err|={abs(pred-ref):.2f}")


def test_full_table_predictions_within_1_unit(sim):
    """All 15 tabulated acids must be within 1.0 pKa unit."""
    for name, pred, ref in sim.all_predictions():
        assert abs(pred - ref) < 1.0, (
            f"{name}: pred={pred:.2f}, ref={ref:.2f}")


def test_strong_acids_predicted_negative(sim):
    """HCl, HBr, HI, H2SO4, HNO3 should all be predicted pKa < 0."""
    for name in ("HCl", "HBr", "HI", "H2SO4", "HNO3"):
        assert sim.predict_pka(name) < 0.0


def test_weak_acids_predicted_positive(sim):
    """Acetic, water, ammonium, methanol -> pKa > 0."""
    for name in ("Acetic acid", "Water", "Ammonium", "Methanol"):
        assert sim.predict_pka(name) > 0.0


def test_pka_formula_consistent_with_RT_ln10(sim):
    """pKa = DeltaG / (R T ln 10)  identity check."""
    geom = sim._geom_for("Acetic acid")
    dg_j = geom.delta_g_dissociation() * 1.0e3
    expected = dg_j / (R_GAS * T_STD * LN10)
    pred = sim.predict_pka("Acetic acid")
    assert pred == pytest.approx(expected, rel=1e-9)


def test_unknown_acid_raises(sim):
    with pytest.raises(KeyError):
        sim.predict_pka("XYZ-not-an-acid")


# --------------------------------------------------------------------- #
# (3) Buffer / Henderson-Hasselbalch                                    #
# --------------------------------------------------------------------- #

def test_henderson_hasselbalch_at_equimolar_returns_pka(sim):
    """[A-]/[HA] = 1 -> pH = pKa exactly."""
    pka = sim.reference_pka("Acetic acid")
    assert sim.henderson_hasselbalch("Acetic acid", 1.0) == pytest.approx(pka)


def test_henderson_hasselbalch_factor_10_shifts_one_unit(sim):
    """[A-]/[HA] = 10 -> pH = pKa + 1; [A-]/[HA] = 0.1 -> pH = pKa - 1."""
    pka = sim.reference_pka("Ammonium")
    assert sim.henderson_hasselbalch("Ammonium", 10.0) == pytest.approx(pka + 1.0)
    assert sim.henderson_hasselbalch("Ammonium", 0.1) == pytest.approx(pka - 1.0)


def test_buffer_pH_equimolar(sim):
    pka = sim.reference_pka("Acetic acid")
    assert sim.buffer_pH("Acetic acid", 0.10, 0.10) == pytest.approx(pka)


def test_henderson_hasselbalch_negative_ratio_raises(sim):
    with pytest.raises(ValueError):
        sim.henderson_hasselbalch("Acetic acid", -1.0)


def test_fraction_dissociated_at_pka_is_half(sim):
    """alpha(pH = pKa) = 0.5"""
    pka = sim.reference_pka("Acetic acid")
    assert sim.fraction_dissociated("Acetic acid", pka) == pytest.approx(0.5)


def test_fraction_dissociated_far_above_pka_is_one(sim):
    pka = sim.reference_pka("Acetic acid")
    assert sim.fraction_dissociated("Acetic acid", pka + 5.0) > 0.999


def test_fraction_dissociated_far_below_pka_is_zero(sim):
    pka = sim.reference_pka("Acetic acid")
    assert sim.fraction_dissociated("Acetic acid", pka - 5.0) < 0.001


def test_water_pH_acetic_acid_0p1_M(sim):
    """0.1 M acetic acid: textbook pH ~ 2.87"""
    pH = sim.water_pH("Acetic acid", 0.10)
    assert pH == pytest.approx(2.87, abs=0.05)


# --------------------------------------------------------------------- #
# (4) Titration curves                                                  #
# --------------------------------------------------------------------- #

def test_titration_strong_acid_starts_low(sim):
    """0.10 M HCl titrated -> initial pH ~ 1."""
    v_mL, pH = sim.titration_curve("HCl", c_acid=0.10,
                                   v_acid_L=0.025, c_base=0.10)
    assert pH[0] < 1.5


def test_titration_strong_acid_ends_high(sim):
    """0.10 M HCl with NaOH -> final pH > 11 at 2x equivalence."""
    v_mL, pH = sim.titration_curve("HCl", c_acid=0.10,
                                   v_acid_L=0.025, c_base=0.10)
    assert pH[-1] > 11.0


def test_titration_weak_acid_buffer_at_half_equivalence(sim):
    """At half-equivalence the pH equals pKa for a weak acid."""
    v_mL, pH = sim.titration_curve("Acetic acid", c_acid=0.10,
                                   v_acid_L=0.025, c_base=0.10,
                                   n_points=400)
    v_eq = sim.equivalence_volume_mL(0.10, 0.025, 0.10)
    idx = int(np.argmin(np.abs(v_mL - 0.5 * v_eq)))
    pka = sim.reference_pka("Acetic acid")
    assert pH[idx] == pytest.approx(pka, abs=0.1)


def test_titration_weak_acid_equivalence_pH_above_7(sim):
    """Equivalence point of weak acid + strong base has pH > 7."""
    v_mL, pH = sim.titration_curve("Acetic acid", c_acid=0.10,
                                   v_acid_L=0.025, c_base=0.10,
                                   n_points=401)
    v_eq = sim.equivalence_volume_mL(0.10, 0.025, 0.10)
    idx = int(np.argmin(np.abs(v_mL - v_eq)))
    # Textbook value for 0.05 M sodium acetate is pH ~ 8.72;
    # require well above 7 (above 7.5 is comfortably above neutral).
    assert pH[idx] > 7.5


def test_titration_strong_acid_equivalence_pH_near_7(sim):
    v_mL, pH = sim.titration_curve("HCl", c_acid=0.10,
                                   v_acid_L=0.025, c_base=0.10,
                                   n_points=401)
    v_eq = sim.equivalence_volume_mL(0.10, 0.025, 0.10)
    idx = int(np.argmin(np.abs(v_mL - v_eq)))
    assert 6.5 < pH[idx] < 7.5


def test_titration_curve_monotone_increasing(sim):
    v_mL, pH = sim.titration_curve("Acetic acid", c_acid=0.10,
                                   v_acid_L=0.025, c_base=0.10,
                                   n_points=400)
    diffs = np.diff(pH)
    assert np.all(diffs >= -1.0e-6)


def test_titration_curve_shape(sim):
    v_mL, pH = sim.titration_curve("Acetic acid", n_points=123)
    assert v_mL.shape == (123,)
    assert pH.shape == (123,)


def test_equivalence_volume_arithmetic(sim):
    v_eq = sim.equivalence_volume_mL(0.10, 0.025, 0.10)
    assert v_eq == pytest.approx(25.0)


# --------------------------------------------------------------------- #
# Substrate metadata                                                    #
# --------------------------------------------------------------------- #

def test_simulator_signature_keys(sim):
    sig = sim.substrate_signature()
    for k in ("ontology", "n_acids_tabulated",
              "T_K", "kT_kJ_per_mol", "RT_ln10_kJ_per_mol"):
        assert k in sig


def test_simulator_default_table_has_15_acids(sim):
    assert len(sim.table) == 15
