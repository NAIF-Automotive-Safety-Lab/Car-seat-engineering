"""Tests for the topological-defect zoo (kinks, vortices, monopoles, textures)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.topological_defect_zoo import (
    DEFAULT_XI,
    DefectVisualizer,
    Kink1D,
    Monopole3D,
    Texture3D,
    Vortex2D,
    classify_defects,
    zoo_summary,
)


# =============================================================================
# Kink1D
# =============================================================================


def test_kink_charge_quantised_plus_one():
    k = Kink1D(L=40.0, n_points=2001)
    Q = k.topological_charge()
    assert Q == pytest.approx(+1.0, abs=5e-3)


def test_antikink_charge_minus_one():
    k = Kink1D(chirality=-1, L=40.0, n_points=2001)
    Q = k.topological_charge()
    assert Q == pytest.approx(-1.0, abs=5e-3)


def test_kink_energy_matches_8_over_xi():
    for xi in (0.5, 1.0, 2.0):
        k = Kink1D(xi=xi, L=40.0 * xi, n_points=4001)
        E_num = k.energy()
        E_an = Kink1D.analytic_energy(xi)
        assert E_num == pytest.approx(E_an, rel=5e-3)
        assert math.isfinite(E_num)


def test_kink_field_asymptotes():
    k = Kink1D(L=40.0)
    x = np.array([-30.0, 30.0])
    phi = k.field(x)
    assert phi[0] == pytest.approx(0.0, abs=1e-10)
    assert phi[1] == pytest.approx(2.0 * math.pi, abs=1e-10)


# =============================================================================
# Vortex2D
# =============================================================================


@pytest.mark.parametrize("n", [-2, -1, 1, 2, 3])
def test_vortex_winding_quantised(n):
    v = Vortex2D(n=n)
    w = v.winding_number()
    assert w == pytest.approx(float(n), abs=1e-6)


def test_vortex_energy_finite_and_positive():
    v = Vortex2D(n=1)
    E = v.energy_per_length()
    assert math.isfinite(E)
    assert E > 0.0


def test_vortex_energy_grows_with_winding_squared():
    e1 = Vortex2D(n=1).energy_per_length()
    e3 = Vortex2D(n=3).energy_per_length()
    # for log-divergent global vortex the n^2 scaling dominates the gradient term
    assert e3 > 4.0 * e1


def test_vortex_amplitude_vanishes_at_core():
    v = Vortex2D(n=1)
    amp = v.amplitude()
    centre = amp[v.n_points // 2, v.n_points // 2]
    assert centre == pytest.approx(0.0, abs=1e-2)


# =============================================================================
# Monopole3D
# =============================================================================


def test_monopole_analytic_charge_one():
    m = Monopole3D()
    assert m.topological_charge_analytic() == 1


def test_monopole_numeric_charge_one():
    m = Monopole3D()
    Q = m.topological_charge_numeric()
    assert Q == pytest.approx(1.0, rel=1e-3)


def test_monopole_energy_finite_and_positive():
    m = Monopole3D()
    E = m.energy()
    assert math.isfinite(E)
    assert E > 0.0


def test_monopole_orientability_forbidden_default():
    m = Monopole3D()
    assert m.is_orientability_forbidden() is True
    assert "FORBIDDEN" in m.substrate_role()


def test_monopole_orientability_optional_off():
    m = Monopole3D(enforce_orientability=False)
    assert m.is_orientability_forbidden() is False


def test_monopole_unit_vector_unit_length():
    m = Monopole3D(n_points=21, L=2.0)
    n_hat = m.unit_vector()
    norms = np.sqrt(np.sum(n_hat ** 2, axis=-1))
    # except possibly at the centre, every vector should have unit length
    centre = (m.n_points // 2,) * 3
    mask = np.ones_like(norms, dtype=bool)
    mask[centre] = False
    assert np.allclose(norms[mask], 1.0, atol=1e-6)


# =============================================================================
# Texture3D
# =============================================================================


def test_texture_baryon_number_one():
    t = Texture3D()
    B = t.baryon_number()
    assert B == pytest.approx(1.0, abs=1e-3)


def test_texture_energy_finite_and_positive():
    t = Texture3D()
    E = t.energy()
    assert math.isfinite(E)
    assert E > 0.0


def test_texture_F_profile_boundary_conditions():
    t = Texture3D()
    r = np.array([0.0, 1e3 * t.R])
    F = t.F_profile(r)
    assert F[0] == pytest.approx(math.pi, abs=1e-12)
    assert F[1] == pytest.approx(0.0, abs=1e-3)


def test_texture_field_normalised_to_unity():
    t = Texture3D(n_points=15)
    comps = t.field_components()
    s2 = (comps["n0"] ** 2 + comps["n1"] ** 2 +
          comps["n2"] ** 2 + comps["n3"] ** 2)
    assert np.allclose(s2, 1.0, atol=1e-6)


# =============================================================================
# Classification & summary
# =============================================================================


def test_classification_has_all_four_classes():
    cls = classify_defects()
    assert "pi_0(Z_2)" in cls
    assert "pi_1(U(1))" in cls
    assert "pi_2(S^2)" in cls
    assert "pi_3(S^3)" in cls


def test_zoo_summary_runs():
    s = zoo_summary()
    assert s["Kink1D"]["Q"] == pytest.approx(1.0, abs=1e-2)
    assert s["Vortex2D"]["winding_n"] == 1
    assert s["Monopole3D"]["Q_analytic"] == 1
    assert s["Texture3D"]["B"] == pytest.approx(1.0, abs=1e-2)
    assert "pi_0(Z_2)" in s["classification"]


# =============================================================================
# Visualizer (smoke tests — produce figures without erroring)
# =============================================================================


def test_visualizer_kink_smoke():
    vis = DefectVisualizer(dpi=60)
    fig = vis.plot_kink(Kink1D())
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_visualizer_vortex_smoke():
    vis = DefectVisualizer(dpi=60)
    fig = vis.plot_vortex(Vortex2D(n=2, n_points=41))
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_visualizer_monopole_smoke():
    vis = DefectVisualizer(dpi=60)
    fig = vis.plot_monopole(Monopole3D(n_points=15, L=2.0))
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_visualizer_texture_smoke():
    vis = DefectVisualizer(dpi=60)
    fig = vis.plot_texture(Texture3D(n_points=15, L=2.0))
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close(fig)


# =============================================================================
# Substrate connections — documentation-level assertions
# =============================================================================


def test_substrate_role_strings():
    assert "fermion" in Kink1D().substrate_role().lower()
    assert "string" in Vortex2D().substrate_role().lower()
    assert "forbidden" in Monopole3D().substrate_role().lower()
    assert "baryon" in Texture3D().substrate_role().lower()


def test_default_xi_constant():
    assert DEFAULT_XI == 1.0
