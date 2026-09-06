"""Tests for the paired biological-membrane GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) bilayer thickness = 4.0 nm (within 5 percent of textbook)
    (b) bending modulus = 20 kT (POPC, body temperature)
    (c) lateral diffusion coefficient ~ 1 micron^2/s
    (d) self-assembly delta_G < 0 (favourable bilayer formation)
    (e) Helfrich curvature energy reproduces 8 pi kappa for a closed sphere
    (f) vesicle stability above the critical radius
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.membrane_substrate import (
    BENDING_MODULUS_KT,
    BILAYER_THICKNESS_NM,
    HEAD_THICKNESS_PER_LEAFLET_NM,
    LATERAL_DIFFUSION_UM2_PER_S,
    LIPID_AREA_NM2,
    TAIL_THICKNESS_PER_LEAFLET_NM,
    MembraneGeometry,
    MembraneSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> MembraneGeometry:
    return MembraneGeometry(n_lipids_per_leaflet=64)


@pytest.fixture
def sim() -> MembraneSimulator:
    return MembraneSimulator()


# --------------------------------------------------------------------- #
# (a) Bilayer thickness = 4.0 nm                                        #
# --------------------------------------------------------------------- #

def test_bilayer_thickness_4nm(geom):
    """Total head-to-head bilayer thickness must equal 4.0 nm."""
    d = geom.total_thickness_nm()
    assert abs(d - 4.0) / 4.0 < 0.05, (
        f"bilayer thickness {d:.2f} nm vs target 4.0 nm"
    )


def test_thickness_is_sum_of_layers(geom):
    """4.0 nm should equal 2 head + 2 tail thickness."""
    expected = 2.0 * (HEAD_THICKNESS_PER_LEAFLET_NM + TAIL_THICKNESS_PER_LEAFLET_NM)
    assert abs(BILAYER_THICKNESS_NM - expected) < 1e-9


def test_hydrophobic_thickness_two_leaflet(geom):
    """Hydrophobic core = 2x tail thickness per leaflet."""
    h = geom.hydrophobic_thickness_nm()
    assert abs(h - 2.0 * TAIL_THICKNESS_PER_LEAFLET_NM) < 1e-9


# --------------------------------------------------------------------- #
# Lipid placement (heads + tails, two leaflets)                         #
# --------------------------------------------------------------------- #

def test_two_leaflets_with_correct_z(geom):
    upper = geom.lipid_head_positions("upper")
    lower = geom.lipid_head_positions("lower")
    assert upper.shape[0] == geom.n_lipids_per_leaflet
    assert lower.shape[0] == geom.n_lipids_per_leaflet
    assert np.allclose(upper[:, 2], +geom.bilayer_thickness_nm / 2.0)
    assert np.allclose(lower[:, 2], -geom.bilayer_thickness_nm / 2.0)


def test_tails_point_to_midplane(geom):
    upper_tails = geom.tail_endpoint_positions("upper")
    lower_tails = geom.tail_endpoint_positions("lower")
    # Upper-leaflet tails sit BELOW their heads, in [0, +d/2]
    assert np.all(upper_tails[:, 2] < +geom.bilayer_thickness_nm / 2.0)
    assert np.all(upper_tails[:, 2] >= 0.0 - 1e-9)
    # Lower-leaflet tails sit ABOVE their heads, in [-d/2, 0]
    assert np.all(lower_tails[:, 2] > -geom.bilayer_thickness_nm / 2.0)
    assert np.all(lower_tails[:, 2] <= 0.0 + 1e-9)


def test_invalid_leaflet_raises(geom):
    with pytest.raises(ValueError):
        geom.lipid_head_positions("middle")


def test_patch_area_equals_lipid_count_times_area(geom):
    A = geom.patch_area_nm2()
    expected = geom.n_lipids_per_leaflet * geom.area_per_lipid_nm2
    assert abs(A - expected) < 1e-9


def test_lipid_area_textbook(geom):
    """POPC fluid-phase area per lipid ~ 0.65 nm^2."""
    assert abs(LIPID_AREA_NM2 - 0.65) < 0.05


def test_integral_proteins_lie_in_patch(geom):
    """Fluid-mosaic: integral proteins centred on bilayer midplane."""
    pos = geom.integral_protein_positions(n_proteins=8)
    L = geom.patch_side_length_nm()
    assert pos.shape == (8, 3)
    assert np.all(np.abs(pos[:, 0]) <= L / 2.0 + 1e-9)
    assert np.all(np.abs(pos[:, 1]) <= L / 2.0 + 1e-9)
    assert np.allclose(pos[:, 2], 0.0)


# --------------------------------------------------------------------- #
# (b) Bending modulus = 20 kT                                           #
# --------------------------------------------------------------------- #

def test_bending_modulus_20_kT(sim):
    """POPC bilayer at body temperature has kappa ~ 20 kT."""
    kappa = sim.bending_modulus_kT_value()
    assert abs(kappa - 20.0) < 1e-9


def test_bending_modulus_in_joules_positive(sim):
    """Bending modulus must be a positive energy (J)."""
    kappa_J = sim.bending_modulus_J()
    assert kappa_J > 0
    # 20 * kT_310 ~ 20 * 4.28e-21 = 8.6e-20 J
    assert 1e-20 < kappa_J < 1e-19


def test_helfrich_zero_on_flat_membrane(sim):
    """Flat membrane (c1 = c2 = 0) has zero curvature energy."""
    e = sim.helfrich_curvature_energy_density(c1=0.0, c2=0.0)
    assert abs(e) < 1e-12


def test_helfrich_sphere_8pi_kappa(sim):
    """Closed spherical vesicle has bending energy 8 pi kappa."""
    expected = 8.0 * math.pi * sim.bending_modulus_kT
    e = sim.vesicle_bending_energy_kT()
    assert abs(e - expected) / expected < 1e-9


def test_helfrich_quadratic_in_curvature(sim):
    """E should scale as c^2 for a uniformly curved patch."""
    e1 = sim.helfrich_curvature_energy_density(c1=0.1, c2=0.1)
    e2 = sim.helfrich_curvature_energy_density(c1=0.2, c2=0.2)
    # Doubling c -> 4x energy (only mean-curvature term, since K = c^2 also 4x;
    # ratio of total energy is 4 because both terms scale as c^2)
    assert abs(e2 / e1 - 4.0) < 1e-6


# --------------------------------------------------------------------- #
# (c) Lateral diffusion ~ 1 um^2/s                                      #
# --------------------------------------------------------------------- #

def test_lateral_diffusion_1_um2_per_s(sim):
    """Fluid-phase POPC bilayer has D ~ 1 um^2 / s."""
    D = sim.lateral_diffusion_coefficient_um2_per_s()
    assert abs(D - LATERAL_DIFFUSION_UM2_PER_S) < 1e-9
    assert 0.5 < D < 5.0  # textbook range for fluid bilayer lipids


def test_msd_linear_in_time(sim):
    """Mean-square displacement is linear in t for a 2D random walk."""
    r1 = sim.msd_um2(1.0)
    r2 = sim.msd_um2(2.0)
    assert abs(r2 / r1 - 2.0) < 1e-9


def test_diffusion_distance_one_micron_in_one_second(sim):
    """D = 1 um^2/s => RMS displacement ~ 2 um in 1 s (4 D t = 4 um^2)."""
    d = sim.diffusion_distance_nm(t_s=1.0)
    assert 1.5e3 < d < 2.5e3, f"diffusion distance {d:.0f} nm in 1 s"


# --------------------------------------------------------------------- #
# (d) Self-assembly free energy                                         #
# --------------------------------------------------------------------- #

def test_self_assembly_is_favourable(sim):
    """Bilayer formation must release free energy (delta G < 0)."""
    dG = sim.self_assembly_free_energy_per_lipid_kcal_per_mol()
    assert dG < 0.0, f"delta G assembly = {dG:.2f} kcal/mol (should be negative)"


def test_self_assembly_dominates_kT(sim):
    """Hydrophobic free energy per lipid >> kT -> robust assembly."""
    from src.stiff_medium.membrane_substrate import KT_AT_BODY_T_KCAL_PER_MOL
    dG = abs(sim.self_assembly_free_energy_per_lipid_kcal_per_mol())
    assert dG > 10.0 * KT_AT_BODY_T_KCAL_PER_MOL


def test_critical_micelle_concentration_mM_range(sim):
    """CMC for typical detergents is in the mM range; check it's positive."""
    cmc = sim.critical_micelle_concentration_M()
    assert cmc > 0.0
    # Order-of-magnitude only: CMC for n_ch2=12 detergents ~ 1 mM -- 1 M
    assert 1e-7 < cmc < 1e2


# --------------------------------------------------------------------- #
# (e) Vesicle stability                                                 #
# --------------------------------------------------------------------- #

def test_vesicle_critical_radius_nm_range(sim):
    """R_crit = 4 kappa / gamma_edge -- expect tens of nm for POPC."""
    R = sim.vesicle_critical_radius_nm()
    assert 5.0 < R < 100.0, f"R_crit = {R:.1f} nm outside expected 5-100 nm"


def test_vesicle_above_critical_is_stable(sim):
    """Liposome of 100 nm radius must be stable."""
    assert sim.vesicle_is_stable(radius_nm=100.0) is True


def test_vesicle_below_critical_unstable(sim):
    """A 1 nm vesicle is far below R_crit -- unstable."""
    assert sim.vesicle_is_stable(radius_nm=1.0) is False


# --------------------------------------------------------------------- #
# Capillary-wave undulation spectrum                                    #
# --------------------------------------------------------------------- #

def test_undulation_spectrum_diverges_at_long_wavelength(sim):
    """<|h_q|^2> ~ 1/q^4 -> diverges as q -> 0."""
    q = np.array([0.01, 0.1, 1.0])
    s = sim.undulation_spectrum(q)
    # Smaller q -> larger amplitude
    assert s[0] > s[1] > s[2]
    # Slope: log s ~ -4 log q
    slope = (math.log(s[2]) - math.log(s[0])) / (math.log(q[2]) - math.log(q[0]))
    assert abs(slope + 4.0) < 0.01, f"undulation slope {slope:.3f} (target -4)"


def test_undulation_amplitude_scales_with_patch(sim):
    """RMS undulation amplitude grows with patch size."""
    a_small = sim.rms_undulation_amplitude_nm(patch_side_nm=10.0)
    a_large = sim.rms_undulation_amplitude_nm(patch_side_nm=100.0)
    assert a_large > a_small
    assert a_small > 0.0


# --------------------------------------------------------------------- #
# Bending mode shape                                                    #
# --------------------------------------------------------------------- #

def test_bending_mode_shape_extrema(sim):
    """A (1,1) bending mode peaks at the patch centre."""
    X, Y, Z = sim.bending_mode_shape(side_nm=20.0, n_grid=41, nx=1, ny=1, amplitude_nm=1.0)
    # Centre cell is at index (20, 20)
    assert abs(Z[20, 20] - 1.0) < 1e-6
    # Edges are clamped to zero
    assert abs(Z[0, 0]) < 1e-6
    assert abs(Z[-1, -1]) < 1e-6


# --------------------------------------------------------------------- #
# Substrate metadata                                                    #
# --------------------------------------------------------------------- #

def test_geometry_signature_has_lambda_qcd(geom):
    sig = geom.substrate_signature()
    assert "lambda_qcd_MeV" in sig
    assert sig["lambda_qcd_MeV"] == 200.0
    assert sig["leaflets"] == 2


def test_simulator_signature_records_constants(sim):
    sig = sim.substrate_signature()
    assert sig["bending_modulus_kT"] == BENDING_MODULUS_KT
    assert sig["lateral_diffusion_um2_per_s"] == LATERAL_DIFFUSION_UM2_PER_S
