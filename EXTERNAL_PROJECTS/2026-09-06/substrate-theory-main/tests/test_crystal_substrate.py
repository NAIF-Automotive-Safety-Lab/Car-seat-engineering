"""Tests for crystal_substrate (paired GEOMETRY+SIM, B3 / Stiff-Medium)."""

import math

import numpy as np
import pytest

from stiff_medium.crystal_substrate import (
    REFERENCE,
    CrystalGeometry,
    CrystalSimulator,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _close(a: float, b: float, pct: float) -> bool:
    return abs(a - b) / abs(b) * 100.0 < pct


# ---------------------------------------------------------------------------
# 1. Geometry — Bravais lattice catalogue
# ---------------------------------------------------------------------------

def test_bravais_catalogue_has_14_lattices():
    geom = CrystalGeometry()
    names = geom.BRAVAIS_NAMES
    assert len(names) == 14, f"expected 14 Bravais lattices, got {len(names)}"


def test_bravais_cubic_cF_motif_is_fcc():
    L, motif = CrystalGeometry.bravais_basis("cF", a=1.0)
    assert motif.shape == (4, 3)
    # Origin always present
    assert np.allclose(motif[0], [0.0, 0.0, 0.0])
    # Other three at face centres
    expected = {(0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)}
    given = {tuple(m) for m in motif[1:]}
    assert given == expected


def test_bravais_cubic_cI_motif_is_bcc():
    L, motif = CrystalGeometry.bravais_basis("cI", a=1.0)
    assert motif.shape == (2, 3)
    assert np.allclose(motif[0], [0.0, 0.0, 0.0])
    assert np.allclose(motif[1], [0.5, 0.5, 0.5])


def test_bravais_hexagonal_vectors():
    L, _ = CrystalGeometry.bravais_basis("hP", a=1.0, c=1.633, gamma=120.0)
    # |a1| = |a2| = a; |a3| = c; angle(a1,a2) = 120°
    assert _close(np.linalg.norm(L[0]), 1.0, 0.01)
    assert _close(np.linalg.norm(L[1]), 1.0, 0.01)
    assert _close(np.linalg.norm(L[2]), 1.633, 0.01)
    cos_g = float(np.dot(L[0], L[1]) /
                  (np.linalg.norm(L[0]) * np.linalg.norm(L[1])))
    assert _close(cos_g, math.cos(math.radians(120.0)), 0.5)


# ---------------------------------------------------------------------------
# 2. Geometry — NaCl, diamond, graphite, CsCl, ZnS, perovskite
# ---------------------------------------------------------------------------

def test_nacl_geometry_4_na_4_cl_per_cell():
    nacl = CrystalGeometry.nacl()
    assert nacl["Na"].shape == (4, 3)
    assert nacl["Cl"].shape == (4, 3)
    assert nacl["z"] == 4


def test_nacl_geometry_a_matches_xrd():
    nacl = CrystalGeometry.nacl()
    # Using the default lattice constant
    assert _close(nacl["a"], REFERENCE["a_NaCl"], 1.0)


def test_diamond_geometry_8_atoms_per_cell():
    d = CrystalGeometry.diamond()
    assert d["C"].shape == (8, 3)
    assert d["z"] == 8


def test_diamond_nearest_neighbour_distance():
    d = CrystalGeometry.diamond()
    a = d["a"]
    pos = d["C"]
    # Compute the minimum pairwise distance
    diffs = pos[:, None, :] - pos[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    np.fill_diagonal(dists, np.inf)
    nn = dists.min()
    # Expected NN = a*sqrt(3)/4 ≈ 1.544 Å for diamond
    expected = a * math.sqrt(3.0) / 4.0
    assert _close(nn, expected, 0.5), f"diamond NN got {nn} expected {expected}"


def test_graphite_4_carbons_per_cell():
    g = CrystalGeometry.graphite()
    assert g["C"].shape == (4, 3)
    assert g["z"] == 4


def test_graphite_in_plane_C_C_distance():
    g = CrystalGeometry.graphite()
    a = g["a"]
    # In-plane NN C-C distance is a/sqrt(3) ≈ 1.42 Å for a=2.46 Å
    expected = a / math.sqrt(3.0)
    # Find within-layer minimum distance (z=0 atoms)
    z0 = g["C"][np.abs(g["C"][:, 2]) < 1e-6]
    if len(z0) >= 2:
        d01 = float(np.linalg.norm(z0[0] - z0[1]))
        assert _close(d01, expected, 1.0), f"got {d01} expected {expected}"


def test_cscl_geometry_one_each():
    cs = CrystalGeometry.cscl()
    assert cs["Cs"].shape == (1, 3)
    assert cs["Cl"].shape == (1, 3)


def test_zns_geometry_4_each():
    z = CrystalGeometry.zns()
    assert z["Zn"].shape == (4, 3)
    assert z["S"].shape == (4, 3)


def test_perovskite_ABX3():
    p = CrystalGeometry.perovskite()
    assert p["Sr"].shape == (1, 3)
    assert p["Ti"].shape == (1, 3)
    assert p["O"].shape == (3, 3)


# ---------------------------------------------------------------------------
# 3. Lattice constants from substrate ξ
# ---------------------------------------------------------------------------

def test_lattice_constant_NaCl_from_xi_default():
    sim = CrystalSimulator()
    a_m = sim.lattice_constant_from_xi("NaCl")
    a_ang = a_m * 1e10
    assert _close(a_ang, REFERENCE["a_NaCl"], 0.1)


def test_lattice_constant_scales_with_xi():
    """Doubling ξ should double every lattice constant."""
    sim1 = CrystalSimulator()
    a1 = sim1.lattice_constant_from_xi("NaCl")
    sim2 = CrystalSimulator(xi_m=2.0 * sim1.xi_m)
    a2 = sim2.lattice_constant_from_xi("NaCl")
    assert _close(a2, 2.0 * a1, 0.001)


def test_lattice_constants_all_supported_crystals():
    sim = CrystalSimulator()
    for name in ("NaCl", "CsCl", "ZnS", "diamond", "graphite", "SrTiO3"):
        a = sim.lattice_constant_from_xi(name)
        assert a > 1e-11 and a < 1e-8, f"{name}: a={a}"


# ---------------------------------------------------------------------------
# 4. Madelung constants
# ---------------------------------------------------------------------------

def test_madelung_NaCl_within_1pct():
    M = CrystalSimulator.madelung_nacl(N_max=15)
    assert _close(M, REFERENCE["M_NaCl"], 1.0), f"got {M}, ref {REFERENCE['M_NaCl']}"


def test_madelung_NaCl_high_precision():
    """N_max=20 should give the full literature value 1.747565."""
    M = CrystalSimulator.madelung_nacl(N_max=20)
    assert abs(M - 1.747565) < 1e-3


def test_madelung_CsCl_within_1pct():
    M = CrystalSimulator.madelung_cscl(N_real=8, N_recip=8)
    assert _close(M, REFERENCE["M_CsCl"], 1.0), f"got {M}, ref {REFERENCE['M_CsCl']}"


# ---------------------------------------------------------------------------
# 5. Packing fractions
# ---------------------------------------------------------------------------

def test_fcc_packing_pi_over_3_sqrt2():
    phi = CrystalSimulator.packing_fraction_fcc()
    assert abs(phi - math.pi / (3.0 * math.sqrt(2.0))) < 1e-15


def test_fcc_packing_numerical_value():
    phi = CrystalSimulator.packing_fraction_fcc()
    assert _close(phi, 0.7405, 0.01)


def test_bcc_packing_pi_sqrt3_over_8():
    phi = CrystalSimulator.packing_fraction_bcc()
    assert abs(phi - math.pi * math.sqrt(3.0) / 8.0) < 1e-15
    assert _close(phi, 0.6802, 0.01)


def test_sc_packing_pi_over_6():
    phi = CrystalSimulator.packing_fraction_sc()
    assert abs(phi - math.pi / 6.0) < 1e-15
    assert _close(phi, 0.5236, 0.01)


def test_diamond_packing_pi_sqrt3_over_16():
    phi = CrystalSimulator.packing_fraction_diamond()
    assert abs(phi - math.pi * math.sqrt(3.0) / 16.0) < 1e-15
    assert _close(phi, 0.3401, 0.01)


def test_hcp_equals_fcc_close_packed():
    phi_fcc = CrystalSimulator.packing_fraction_fcc()
    phi_hcp = CrystalSimulator.packing_fraction_hcp()
    assert abs(phi_hcp - phi_fcc) < 1e-15


# ---------------------------------------------------------------------------
# 6. X-ray diffraction (Bragg peaks)
# ---------------------------------------------------------------------------

def test_braggs_law_simple_case():
    # d=2.82 Å, Cu Kα (1.5406 Å), n=1: 2θ should be ~31.7° (NaCl 200 peak)
    two_theta = CrystalSimulator.bragg_two_theta_deg(d_spacing_ang=2.82)
    assert _close(two_theta, 31.7, 1.0)


def test_braggs_law_returns_nan_when_forbidden():
    # d much smaller than λ/2 → sin θ > 1
    two_theta = CrystalSimulator.bragg_two_theta_deg(d_spacing_ang=0.5)
    assert math.isnan(two_theta)


def test_cubic_d_spacing_111():
    # d_{111} = a / sqrt(3)
    a = 5.6402
    d = CrystalSimulator.cubic_d_spacing(a, 1, 1, 1)
    assert _close(d, a / math.sqrt(3.0), 0.001)


def test_fcc_selection_rule():
    sim = CrystalSimulator
    # All-even (200) allowed
    assert sim.fcc_allowed(2, 0, 0)
    # All-odd (111) allowed
    assert sim.fcc_allowed(1, 1, 1)
    # Mixed (210) forbidden
    assert not sim.fcc_allowed(2, 1, 0)
    # Mixed (100) forbidden
    assert not sim.fcc_allowed(1, 0, 0)


def test_diamond_selection_rule():
    sim = CrystalSimulator
    # All-odd (111) allowed
    assert sim.diamond_allowed(1, 1, 1)
    # All-even h+k+l=4 (220) allowed
    assert sim.diamond_allowed(2, 2, 0)
    # All-even h+k+l=2 (200) FORBIDDEN in diamond
    assert not sim.diamond_allowed(2, 0, 0)
    # All-even h+k+l=8 (400) allowed
    assert sim.diamond_allowed(4, 0, 0)
    # All-even h+k+l=6 (222) FORBIDDEN
    assert not sim.diamond_allowed(2, 2, 2)


def test_NaCl_xrd_first_strong_peak_is_200():
    sim = CrystalSimulator()
    peaks = sim.xrd_peaks("NaCl")
    # The strongest peak with high intensity should be (200) at ~31.7°
    strong = [p for p in peaks if p["intensity"] > 50.0]
    assert len(strong) >= 1
    # First strong peak should be (200)
    p0 = strong[0]
    assert p0["hkl"] == (0, 0, 2) or p0["hkl"] == (2, 0, 0) or p0["hkl"] == (0, 2, 0)
    assert _close(p0["two_theta_deg"], 31.7, 1.0)


def test_diamond_xrd_111_at_44deg():
    sim = CrystalSimulator()
    peaks = sim.xrd_peaks("diamond")
    # First peak should be (111) at ~43.9°
    p0 = peaks[0]
    assert p0["hkl"] == (1, 1, 1)
    assert _close(p0["two_theta_deg"], 43.9, 1.0)


def test_diamond_xrd_220_at_75deg():
    sim = CrystalSimulator()
    peaks = sim.xrd_peaks("diamond")
    # Find the (220) peak (== (022) in the canonical ordering used)
    p220 = next(
        (p for p in peaks if sorted(p["hkl"]) == [0, 2, 2]),
        None,
    )
    assert p220 is not None
    assert _close(p220["two_theta_deg"], 75.3, 1.0)


def test_xrd_peaks_sorted_by_two_theta():
    sim = CrystalSimulator()
    peaks = sim.xrd_peaks("NaCl")
    angles = [p["two_theta_deg"] for p in peaks]
    assert angles == sorted(angles)


# ---------------------------------------------------------------------------
# 7. Mass density (consistency with XRD reference)
# ---------------------------------------------------------------------------

def test_NaCl_mass_density_matches():
    sim = CrystalSimulator()
    rho = sim.density_NaCl()
    assert _close(rho, REFERENCE["rho_NaCl"], 1.0)


def test_diamond_mass_density_matches():
    sim = CrystalSimulator()
    rho = sim.density_diamond()
    assert _close(rho, REFERENCE["rho_diamond"], 1.0)


# ---------------------------------------------------------------------------
# 8. Convenience report
# ---------------------------------------------------------------------------

def test_report_NaCl_full_set():
    sim = CrystalSimulator()
    r = sim.report("NaCl")
    assert "a_metres" in r
    assert "a_ang" in r
    assert "madelung" in r
    assert "packing" in r
    assert "density_gcm3" in r
    assert _close(r["madelung"], REFERENCE["M_NaCl"], 1.0)
    assert _close(r["packing"], REFERENCE["phi_FCC"], 0.01)


def test_report_diamond():
    sim = CrystalSimulator()
    r = sim.report("diamond")
    assert _close(r["packing"], REFERENCE["phi_diamond"], 0.01)
    assert _close(r["density_gcm3"], REFERENCE["rho_diamond"], 1.0)
