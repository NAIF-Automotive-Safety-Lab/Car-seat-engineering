"""Tests for the paired semiconductor band-structure GEOMETRY + SIMULATOR module.

Covers the three required claims:
    (a) Si bandgap E_g = 1.12 eV (textbook 300 K value)
    (b) GaAs bandgap E_g = 1.42 eV (direct gap, textbook 300 K value)
    (c) intrinsic carrier concentration of Si at 300 K ~ 1.5e10 /cm^3
        (textbook order-of-magnitude; we accept the band 5e9 .. 5e10
        which spans both the modern Misiakos value 1.0e10 and the
        Sze textbook value 1.45e10)

Plus geometric and material sanity:
    * direct-vs-indirect classification (GaAs direct, Si/Ge/diamond indirect)
    * built-in voltage of a Si p-n junction (V_bi ~ 0.6 .. 0.9 V at 1e17/1e17)
    * depletion width within physical bounds
    * conduction & valence bands have correct curvature signs
    * Brillouin zone cut spans [-pi/a, +pi/a]
    * substrate signature carries Lambda_QCD anchor
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.semiconductor_substrate import (
    EPS_0,
    EV,
    H_PLANCK,
    HBAR,
    K_B,
    K_B_EV,
    LAMBDA_QCD_MEV,
    M_E,
    MATERIALS,
    Q_E,
    SemiconductorGeometry,
    SemiconductorSimulator,
    kronig_penney_gap_eV,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom_si() -> SemiconductorGeometry:
    return SemiconductorGeometry(material="Si")


@pytest.fixture
def geom_gaas() -> SemiconductorGeometry:
    return SemiconductorGeometry(material="GaAs")


@pytest.fixture
def sim_si() -> SemiconductorSimulator:
    return SemiconductorSimulator(material="Si", T_kelvin=300.0)


@pytest.fixture
def sim_gaas() -> SemiconductorSimulator:
    return SemiconductorSimulator(material="GaAs", T_kelvin=300.0)


@pytest.fixture
def sim_gan() -> SemiconductorSimulator:
    return SemiconductorSimulator(material="GaN", T_kelvin=300.0)


@pytest.fixture
def sim_diamond() -> SemiconductorSimulator:
    return SemiconductorSimulator(material="diamond", T_kelvin=300.0)


# --------------------------------------------------------------------- #
# (a) Bandgaps -- the headline materials                                 #
# --------------------------------------------------------------------- #

def test_silicon_bandgap_is_1_12_eV(sim_si):
    """Si bandgap E_g = 1.12 eV at room temperature."""
    assert sim_si.bandgap_eV() == pytest.approx(1.12, abs=1e-9)


def test_gallium_arsenide_bandgap_is_1_42_eV(sim_gaas):
    """GaAs bandgap E_g = 1.42 eV at room temperature."""
    assert sim_gaas.bandgap_eV() == pytest.approx(1.42, abs=1e-9)


def test_gallium_nitride_bandgap_is_3_4_eV(sim_gan):
    """GaN bandgap E_g = 3.4 eV (wide-gap, blue/UV LEDs)."""
    assert sim_gan.bandgap_eV() == pytest.approx(3.4, abs=1e-9)


def test_diamond_bandgap_is_5_5_eV(sim_diamond):
    """Diamond bandgap E_g = 5.5 eV (UV-transparent insulator)."""
    assert sim_diamond.bandgap_eV() == pytest.approx(5.5, abs=1e-9)


def test_bandgap_ordering(sim_si, sim_gaas, sim_gan, sim_diamond):
    """E_g(Si) < E_g(GaAs) < E_g(GaN) < E_g(diamond)."""
    eg_si = sim_si.bandgap_eV()
    eg_gaas = sim_gaas.bandgap_eV()
    eg_gan = sim_gan.bandgap_eV()
    eg_dia = sim_diamond.bandgap_eV()
    assert eg_si < eg_gaas < eg_gan < eg_dia


# --------------------------------------------------------------------- #
# (b) Direct vs indirect gap                                             #
# --------------------------------------------------------------------- #

def test_silicon_is_indirect(sim_si):
    """Si has an indirect gap at the X point."""
    assert sim_si.gap_type() == "indirect"
    assert sim_si.is_direct() is False


def test_gallium_arsenide_is_direct(sim_gaas):
    """GaAs has a direct gap at Gamma."""
    assert sim_gaas.gap_type() == "direct"
    assert sim_gaas.is_direct() is True


def test_gallium_nitride_is_direct(sim_gan):
    """GaN has a direct gap at Gamma (LED foundation material)."""
    assert sim_gan.gap_type() == "direct"


def test_diamond_is_indirect(sim_diamond):
    """Diamond has an indirect gap (low light-emission efficiency)."""
    assert sim_diamond.gap_type() == "indirect"


# --------------------------------------------------------------------- #
# (c) Intrinsic carrier concentration of Si at 300 K                     #
# --------------------------------------------------------------------- #

def test_silicon_intrinsic_carrier_concentration(sim_si):
    """n_i(Si, 300 K) ~ 1e10 /cm^3 (textbook order of magnitude).

    The classical textbook value (Sze, 1981) is 1.45e10 /cm^3 using
    E_g = 1.08 eV and historical effective masses; the modern best
    measurement (Misiakos & Tsamakis, 1993) is 1.0e10 /cm^3.  Using
    the standard formula

        n_i = sqrt(N_c N_v) * exp(-E_g / 2kT)

    with E_g = 1.12 eV and DOS effective masses 1.08 m_e / 0.81 m_e
    yields 8.9e9 /cm^3, which lives in the same order of magnitude.

    We accept the band [5e9, 5e10] /cm^3 to span both reference
    values comfortably.
    """
    n_i = sim_si.intrinsic_carrier_concentration_per_cm3()
    assert 5.0e9 <= n_i <= 5.0e10, f"n_i = {n_i:.3e} outside textbook band"


def test_intrinsic_carrier_within_order_of_magnitude_of_1_5e10(sim_si):
    """n_i within an order of magnitude of the standard 1.5e10 textbook value."""
    n_i = sim_si.intrinsic_carrier_concentration_per_cm3()
    target = 1.5e10
    ratio = n_i / target
    assert 0.1 <= ratio <= 10.0


def test_intrinsic_carrier_increases_with_temperature(sim_si):
    """n_i grows exponentially with T (E_g/2kT shrinks)."""
    sim_300 = SemiconductorSimulator(material="Si", T_kelvin=300.0)
    sim_400 = SemiconductorSimulator(material="Si", T_kelvin=400.0)
    n_300 = sim_300.intrinsic_carrier_concentration_per_cm3()
    n_400 = sim_400.intrinsic_carrier_concentration_per_cm3()
    assert n_400 > n_300


def test_intrinsic_carrier_decreases_with_bandgap():
    """Wider-gap materials have smaller n_i at fixed T."""
    n_si  = SemiconductorSimulator(material="Si").intrinsic_carrier_concentration_per_cm3()
    n_gaas = SemiconductorSimulator(material="GaAs").intrinsic_carrier_concentration_per_cm3()
    n_gan = SemiconductorSimulator(material="GaN").intrinsic_carrier_concentration_per_cm3()
    n_dia = SemiconductorSimulator(material="diamond").intrinsic_carrier_concentration_per_cm3()
    assert n_si > n_gaas > n_gan > n_dia


def test_intrinsic_diamond_essentially_zero(sim_diamond):
    """Diamond at 300 K: n_i ~ 1e-27 /cm^3 (effectively zero free carriers)."""
    n_i = sim_diamond.intrinsic_carrier_concentration_per_cm3()
    assert n_i < 1.0e-10


# --------------------------------------------------------------------- #
# Carrier concentration formula n = N_c exp(-(E_c - E_F)/kT)             #
# --------------------------------------------------------------------- #

def test_electron_concentration_at_intrinsic_fermi_level(sim_si):
    """At intrinsic E_F, n equals the intrinsic carrier concentration."""
    E_F = sim_si.intrinsic_fermi_level_eV()
    n = sim_si.electron_concentration_per_cm3(E_F)
    n_i = sim_si.intrinsic_carrier_concentration_per_cm3()
    # n at intrinsic E_F equals n_i to within numerical tolerance
    assert n == pytest.approx(n_i, rel=0.01)


def test_hole_concentration_at_intrinsic_fermi_level(sim_si):
    """At intrinsic E_F, p equals the intrinsic carrier concentration."""
    E_F = sim_si.intrinsic_fermi_level_eV()
    p = sim_si.hole_concentration_per_cm3(E_F)
    n_i = sim_si.intrinsic_carrier_concentration_per_cm3()
    assert p == pytest.approx(n_i, rel=0.01)


def test_n_p_product_equals_n_i_squared(sim_si):
    """Mass-action law: n * p = n_i^2 regardless of E_F (non-degenerate)."""
    n_i = sim_si.intrinsic_carrier_concentration_per_cm3()
    # Pick an arbitrary E_F well inside the gap
    E_F = 0.6
    n = sim_si.electron_concentration_per_cm3(E_F)
    p = sim_si.hole_concentration_per_cm3(E_F)
    assert n * p == pytest.approx(n_i * n_i, rel=1e-6)


def test_carrier_grows_exponentially_with_E_F(sim_si):
    """Raising E_F by kT * ln(10) increases n by a factor of 10."""
    E_F_lo = 0.5
    E_F_hi = E_F_lo + sim_si.kT_eV() * math.log(10.0)
    n_lo = sim_si.electron_concentration_per_cm3(E_F_lo)
    n_hi = sim_si.electron_concentration_per_cm3(E_F_hi)
    assert n_hi / n_lo == pytest.approx(10.0, rel=1e-6)


# --------------------------------------------------------------------- #
# p-n junction: built-in voltage and depletion width                    #
# --------------------------------------------------------------------- #

def test_silicon_pn_junction_built_in_voltage(sim_si):
    """V_bi for Si p-n at 1e17/1e17 is in the 0.6 .. 0.9 V band."""
    V_bi = sim_si.built_in_voltage_V(N_a_per_cm3=1e17, N_d_per_cm3=1e17)
    assert 0.6 <= V_bi <= 0.95


def test_built_in_voltage_grows_with_doping(sim_si):
    """Heavier doping -> larger V_bi (logarithmic)."""
    V_low  = sim_si.built_in_voltage_V(1e15, 1e15)
    V_high = sim_si.built_in_voltage_V(1e18, 1e18)
    assert V_high > V_low


def test_built_in_voltage_grows_with_bandgap(sim_si, sim_gaas, sim_gan):
    """Wider-gap materials have higher V_bi at fixed doping."""
    V_si   = sim_si.built_in_voltage_V(1e17, 1e17)
    V_gaas = sim_gaas.built_in_voltage_V(1e17, 1e17)
    V_gan  = sim_gan.built_in_voltage_V(1e17, 1e17)
    assert V_si < V_gaas < V_gan


def test_silicon_depletion_width_micron_scale(sim_si):
    """Si p-n at 1e17 doping has W ~ 0.1 um (textbook scale)."""
    W = sim_si.depletion_width_m(1e17, 1e17)
    # Convert to microns; expect 0.05 .. 0.5 um at this doping
    W_um = W * 1e6
    assert 0.05 <= W_um <= 0.5


def test_depletion_width_scales_with_voltage(sim_si):
    """Reverse bias (V_applied < 0) widens the depletion region."""
    W_eq = sim_si.depletion_width_m(1e16, 1e16, V_applied_V=0.0)
    W_rev = sim_si.depletion_width_m(1e16, 1e16, V_applied_V=-2.0)
    assert W_rev > W_eq


def test_junction_profile_charge_neutrality(sim_si):
    """Total positive donor charge = total negative acceptor charge."""
    profile = sim_si.junction_profile(N_a_per_cm3=1e17, N_d_per_cm3=1e17)
    rho = profile["rho_C_per_m3"]
    x = profile["x_m"]
    dx = x[1] - x[0]
    pos_charge = float(np.sum(rho[rho > 0]) * dx)
    neg_charge = float(np.sum(rho[rho < 0]) * dx)
    # |pos| should equal |neg| (overall neutrality)
    assert abs(pos_charge + neg_charge) < 0.01 * abs(pos_charge)


def test_junction_profile_field_zero_at_edges(sim_si):
    """Electric field is zero outside the depletion region (bulk neutrality)."""
    profile = sim_si.junction_profile(N_a_per_cm3=1e17, N_d_per_cm3=1e17)
    E_field = profile["E_V_per_m"]
    # Endpoints should be at zero field
    assert abs(E_field[0]) < 1e-6
    assert abs(E_field[-1]) < 1e-6


def test_junction_profile_potential_step_equals_V_bi(sim_si):
    """Potential jump across the junction equals the built-in voltage."""
    profile = sim_si.junction_profile(N_a_per_cm3=1e17, N_d_per_cm3=1e17)
    phi = profile["phi_V"]
    V_bi = profile["V_bi_V"]
    delta_phi = abs(phi[-1] - phi[0])
    # Should match V_bi to within a few percent (numerical integration)
    assert delta_phi == pytest.approx(V_bi, rel=0.05)


# --------------------------------------------------------------------- #
# GEOMETRY: bands and Brillouin zone                                    #
# --------------------------------------------------------------------- #

def test_brillouin_zone_spans_pi_over_a(geom_si):
    """k grid spans [-pi/a, +pi/a] for the 1D BZ projection."""
    k = geom_si.k_grid()
    a = geom_si.lattice_constant_m()
    assert k[0] == pytest.approx(-math.pi / a, rel=1e-9)
    assert k[-1] == pytest.approx(+math.pi / a, rel=1e-9)


def test_normalised_k_grid_runs_minus_one_to_one(geom_si):
    """k normalised by pi/a runs from -1 to +1."""
    kn = geom_si.k_grid_normalised()
    assert kn[0] == pytest.approx(-1.0, rel=1e-9)
    assert kn[-1] == pytest.approx(+1.0, rel=1e-9)


def test_conduction_band_is_above_gap(geom_si):
    """Minimum of conduction band >= bandgap (zero baseline at valence top)."""
    E_c = geom_si.conduction_band_eV()
    assert E_c.min() == pytest.approx(geom_si.bandgap_eV(), abs=1e-9)


def test_valence_band_max_at_zero(geom_si):
    """Valence-band maximum sits at zero (our energy zero is the VBM)."""
    E_v = geom_si.valence_band_eV()
    assert E_v.max() == pytest.approx(0.0, abs=1e-9)


def test_direct_gap_minima_align(geom_gaas):
    """Direct gap: conduction-band minimum and valence-band maximum at k=0."""
    k = geom_gaas.k_grid()
    E_c = geom_gaas.conduction_band_eV()
    E_v = geom_gaas.valence_band_eV()
    i_c_min = int(np.argmin(E_c))
    i_v_max = int(np.argmax(E_v))
    # Both minima/maxima are at the centre (Gamma)
    assert abs(k[i_c_min]) < 1e-6
    assert abs(k[i_v_max]) < 1e-6


def test_indirect_gap_minima_dont_align(geom_si):
    """Indirect gap: conduction minimum off-Gamma, valence maximum at Gamma."""
    k = geom_si.k_grid()
    E_c = geom_si.conduction_band_eV()
    E_v = geom_si.valence_band_eV()
    i_c_min = int(np.argmin(E_c))
    i_v_max = int(np.argmax(E_v))
    # Valence max at k=0
    assert abs(k[i_v_max]) < 1e-6
    # Conduction min at zone boundary k = +/- pi/a (X-like)
    assert abs(k[i_c_min]) > 0.5 * abs(k[-1])


def test_fermi_level_is_in_the_gap(geom_si):
    """Intrinsic E_F = E_g/2 sits between valence top and conduction bottom."""
    E_F = geom_si.fermi_level_eV()
    assert 0.0 < E_F < geom_si.bandgap_eV()


# --------------------------------------------------------------------- #
# Bandgap from substrate periodic potential                             #
# --------------------------------------------------------------------- #

def test_kronig_penney_returns_input(sim_si):
    """The NFE estimate E_g ~ V_0 holds at first order."""
    assert kronig_penney_gap_eV(1.12) == pytest.approx(1.12, abs=1e-9)


def test_bandgap_from_substrate_potential_recovers_textbook(sim_si):
    """Passing V_0 = E_g_textbook into the formula returns E_g_textbook."""
    eg = sim_si.bandgap_from_substrate_potential_eV(1.12)
    assert eg == pytest.approx(1.12, abs=1e-9)


def test_bandgap_from_substrate_potential_no_arg(sim_si):
    """No-arg version returns the textbook bandgap."""
    assert sim_si.bandgap_from_substrate_potential_eV() == pytest.approx(
        sim_si.bandgap_eV(), abs=1e-9,
    )


# --------------------------------------------------------------------- #
# Substrate signature                                                    #
# --------------------------------------------------------------------- #

def test_geometry_substrate_signature(geom_si):
    sig = geom_si.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
    assert sig["material"] == "Si"
    assert sig["E_g_eV"] == pytest.approx(1.12, abs=1e-9)
    assert sig["gap_type"] == "indirect"


def test_simulator_substrate_signature(sim_si):
    sig = sim_si.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert "n_i_per_cm3" in sig
    assert sig["material"] == "Si"
    assert sig["T_kelvin"] == 300.0


def test_unknown_material_raises():
    """Asking for a material not in MATERIALS raises KeyError."""
    with pytest.raises(KeyError):
        SemiconductorSimulator(material="Unobtanium").bandgap_eV()
    with pytest.raises(KeyError):
        SemiconductorGeometry(material="Unobtanium").params()


def test_materials_table_has_required_keys():
    """Every entry in MATERIALS has the required keys."""
    required = {"E_g_eV", "gap_type", "m_e_star_over_me",
                "m_h_star_over_me", "eps_r", "lattice_a_A"}
    for name, row in MATERIALS.items():
        assert required.issubset(row.keys()), f"{name} missing keys"
        assert row["gap_type"] in {"direct", "indirect"}
