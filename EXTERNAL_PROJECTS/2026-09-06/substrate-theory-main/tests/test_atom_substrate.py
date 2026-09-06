"""Tests for atom_substrate — atoms as K_4 nuclei + substrate-strain orbitals."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.atom_substrate import (
    AUFBAU_ORDER,
    MEASURED_IE_EV,
    RYDBERG_EV,
    SUBSHELL_CAPACITY,
    AtomGeometry,
    AtomSimulator,
    aufbau_configuration,
    configuration_string,
    orbital_density_grid,
    predicted_vs_measured,
)


# ---------------------------------------------------------------------------
# Aufbau filling
# ---------------------------------------------------------------------------

def test_aufbau_hydrogen():
    cfg = aufbau_configuration(1)
    assert cfg == [(1, 0, 1)]
    assert configuration_string(cfg) == "1s1"


def test_aufbau_helium():
    cfg = aufbau_configuration(2)
    assert cfg == [(1, 0, 2)]
    assert configuration_string(cfg) == "1s2"


def test_aufbau_carbon():
    """Carbon Z=6 must fill 1s^2 2s^2 2p^2."""
    cfg = aufbau_configuration(6)
    assert cfg == [(1, 0, 2), (2, 0, 2), (2, 1, 2)]
    assert configuration_string(cfg) == "1s2 2s2 2p2"


def test_aufbau_neon():
    """Neon Z=10 closes the n=2 shell: 1s^2 2s^2 2p^6."""
    cfg = aufbau_configuration(10)
    assert cfg == [(1, 0, 2), (2, 0, 2), (2, 1, 6)]
    assert configuration_string(cfg) == "1s2 2s2 2p6"


def test_subshell_capacity():
    assert SUBSHELL_CAPACITY[0] == 2
    assert SUBSHELL_CAPACITY[1] == 6
    assert SUBSHELL_CAPACITY[2] == 10


# ---------------------------------------------------------------------------
# AtomGeometry: nucleus is a K_4 stack
# ---------------------------------------------------------------------------

def test_hydrogen_geometry():
    atom = AtomGeometry(Z=1)
    assert atom.Z == 1
    assert atom.N == 1                      # default: most-stable isotope
    assert atom.n_nucleons == 2
    assert len(atom.nucleon_cells) == 2
    # Each K_4 cell has 4 vertices
    for cell in atom.nucleon_cells:
        assert cell.vertices.shape == (4, 3)
    assert atom.electron_config == [(1, 0, 1)]


def test_carbon_geometry():
    atom = AtomGeometry(Z=6)
    assert atom.Z == 6
    assert atom.n_nucleons == atom.Z + atom.N
    # Carbon has 6 electrons in 1s^2 2s^2 2p^2
    assert atom.configuration_string() == "1s2 2s2 2p2"
    # Least bound subshell is 2p
    assert atom.least_bound_subshell() == (2, 1)


def test_oxygen_geometry():
    atom = AtomGeometry(Z=8)
    assert atom.configuration_string() == "1s2 2s2 2p4"
    assert atom.least_bound_subshell() == (2, 1)


def test_nucleus_radius_matches_R0_A_third():
    """R = 1.2 * A^(1/3) fm — empirical nuclear-charge-radius scaling."""
    atom = AtomGeometry(Z=6, N=6)            # carbon-12
    A = 12
    expected = 1.2 * A ** (1.0 / 3.0)
    assert atom.nucleus_radius_fm == pytest.approx(expected)


def test_z_invalid():
    with pytest.raises(ValueError):
        AtomGeometry(Z=0)


# ---------------------------------------------------------------------------
# AtomSimulator: hydrogen ground state
# ---------------------------------------------------------------------------

def test_hydrogen_ground_energy_exact():
    """E_1(H) = -Rydberg = -13.605693 eV  (closed form, 0.001% match)."""
    e1 = AtomSimulator.hydrogen_ground_energy_ev()
    assert e1 == pytest.approx(-13.605693, rel=1e-5)


def test_hydrogen_orbital_energy_n2():
    """E_2(H) = -Ry/4 = -3.4014 eV."""
    e2 = AtomSimulator.orbital_energy_ev(Z_eff=1.0, n=2)
    assert e2 == pytest.approx(-RYDBERG_EV / 4.0, rel=1e-12)


def test_hydrogen_orbital_energy_n3():
    """E_3(H) = -Ry/9 = -1.5117 eV."""
    e3 = AtomSimulator.orbital_energy_ev(Z_eff=1.0, n=3)
    assert e3 == pytest.approx(-RYDBERG_EV / 9.0, rel=1e-12)


# ---------------------------------------------------------------------------
# AtomSimulator: ionization energies for H..Ne
# ---------------------------------------------------------------------------

def test_ionization_hydrogen():
    """H IE = 13.606 eV  (0.05% off measured because we use bare Rydberg)."""
    ie = AtomSimulator.ionization_energy_ev(1)
    assert ie == pytest.approx(13.606, rel=1e-3)


def test_ionization_helium_within_1pct():
    """He IE = 24.587 eV  (calibrated Z_eff hits the published value exactly)."""
    ie = AtomSimulator.ionization_energy_ev(2)
    assert ie == pytest.approx(24.587, rel=1e-2)


def test_ionization_lithium_within_5pct():
    """Li IE = 5.392 eV."""
    ie = AtomSimulator.ionization_energy_ev(3)
    assert ie == pytest.approx(5.392, rel=5e-2)


def test_ionization_beryllium_within_5pct():
    ie = AtomSimulator.ionization_energy_ev(4)
    assert ie == pytest.approx(9.323, rel=5e-2)


def test_ionization_carbon_within_5pct():
    ie = AtomSimulator.ionization_energy_ev(6)
    assert ie == pytest.approx(11.260, rel=5e-2)


def test_ionization_neon_within_5pct():
    ie = AtomSimulator.ionization_energy_ev(10)
    assert ie == pytest.approx(21.565, rel=5e-2)


def test_all_ionization_within_5pct():
    """Every element H..Ne reproduces the measured first IE within 5%."""
    rows = predicted_vs_measured(Z_max=10)
    for Z, sym, pred, meas, err_pct in rows:
        assert err_pct <= 5.0, (
            f"{sym} (Z={Z}): predicted {pred:.3f}, measured {meas:.3f}, "
            f"error {err_pct:.2f}%"
        )


# ---------------------------------------------------------------------------
# Substrate-Schroedinger numerical eigensolve cross-check
# ---------------------------------------------------------------------------

def test_radial_schrodinger_hydrogen_1s():
    """Numerical radial eigensolve recovers H 1s = -13.606 eV (within 0.5%)."""
    e1 = AtomSimulator.radial_schrodinger_energy(n=1, ell=0, Z_eff=1.0)
    assert e1 == pytest.approx(-RYDBERG_EV, rel=5e-3)


def test_radial_schrodinger_hydrogen_2p():
    """Numerical radial eigensolve recovers H 2p = -3.401 eV (within 1%)."""
    e2p = AtomSimulator.radial_schrodinger_energy(n=2, ell=1, Z_eff=1.0)
    assert e2p == pytest.approx(-RYDBERG_EV / 4.0, rel=1e-2)


# ---------------------------------------------------------------------------
# Orbital density grids — substrate strain patterns
# ---------------------------------------------------------------------------

def test_orbital_1s_density_centred():
    """1s density peaks at the origin (s-orbital is spherically symmetric)."""
    X, Y, Z, rho = orbital_density_grid(1, 0, 0, extent_bohr=8.0, grid=24)
    centre_idx = rho.shape[0] // 2
    assert rho[centre_idx, centre_idx, centre_idx] == pytest.approx(rho.max(), rel=1e-6)


def test_orbital_2p_z_has_two_lobes():
    """2p_z density has two lobes along z, zero at z=0."""
    X, Y, Z, rho = orbital_density_grid(2, 1, 0, extent_bohr=10.0, grid=25)
    # midplane z=0 has zero density (within numerical noise from grid)
    z_mid = rho.shape[2] // 2
    assert rho[:, :, z_mid].max() < 1e-10
    # but density on z > 0 axis is positive
    assert rho[X.shape[0] // 2, X.shape[1] // 2, :].max() > 0.0


def test_orbital_density_normalised_to_max():
    """Density grids are normalised so max == 1.0 for plotting convenience."""
    for (n, ell, m) in [(1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 1), (3, 2, 0)]:
        _, _, _, rho = orbital_density_grid(n, ell, m, extent_bohr=12.0, grid=20)
        assert rho.max() == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# predicted_vs_measured table sanity
# ---------------------------------------------------------------------------

def test_predicted_vs_measured_returns_ten_rows():
    rows = predicted_vs_measured(Z_max=10)
    assert len(rows) == 10
    assert rows[0][1] == "H"
    assert rows[5][1] == "C"
    assert rows[9][1] == "Ne"


def test_measured_ie_table_has_h_through_ne():
    """The reference table covers H..Ne."""
    for Z in range(1, 11):
        assert Z in MEASURED_IE_EV
