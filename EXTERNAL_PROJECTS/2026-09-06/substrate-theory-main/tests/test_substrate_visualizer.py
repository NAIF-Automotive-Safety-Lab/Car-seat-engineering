"""Tests for SubstrateVisualizer."""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from stiff_medium.substrate_visualizer import SubstrateVisualizer  # noqa: E402


@pytest.fixture
def viz():
    v = SubstrateVisualizer()
    yield v
    v.close_all()


def test_plot_field_2d(viz):
    field = np.random.RandomState(0).randn(32, 48)
    fig = viz.plot_field_2d(field, t=1.5)
    assert isinstance(fig, plt.Figure)
    ax = fig.axes[0]
    # imshow places an AxesImage; the data shape should match
    img = ax.images[0]
    assert img.get_array().shape == (32, 48)
    assert "field_2d" in viz.figures


def test_plot_field_2d_rejects_non_2d(viz):
    with pytest.raises(ValueError):
        viz.plot_field_2d(np.zeros(10), t=0.0)


def test_plot_kink_evolution_dict(viz):
    history = [{"t": float(i), "x": float(i) * 0.5} for i in range(20)]
    fig = viz.plot_kink_evolution(history)
    line = fig.axes[0].lines[0]
    xs, ys = line.get_data()
    assert len(xs) == 20 and len(ys) == 20
    assert ys[-1] == pytest.approx(9.5)


def test_plot_kink_evolution_array(viz):
    arr = np.column_stack([np.linspace(0, 1, 11), np.linspace(0, 2, 11)])
    fig = viz.plot_kink_evolution(arr)
    assert len(fig.axes[0].lines[0].get_xdata()) == 11


def test_plot_dispersion_mapping(viz):
    k = np.linspace(0, 2.0, 50)
    omega = np.sqrt(1.0 + k ** 2)
    fig = viz.plot_dispersion({"k": k, "omega": omega})
    line = fig.axes[0].lines[0]
    assert line.get_xdata().shape == (50,)


def test_plot_dispersion_tuple(viz):
    k = np.linspace(0, 1.0, 20)
    omega = k * 2.0
    fig = viz.plot_dispersion((k, omega))
    assert fig.axes[0].lines[0].get_ydata()[-1] == pytest.approx(2.0)


def test_plot_dispersion_shape_mismatch(viz):
    with pytest.raises(ValueError):
        viz.plot_dispersion((np.arange(5), np.arange(6)))


def test_plot_mass_torque_table_list(viz):
    rows = [
        {"name": "p", "predicted": 938.2, "observed": 938.27},
        {"name": "n", "predicted": 939.4, "observed": 939.57},
        {"name": "Lambda", "predicted": 1115.0, "observed": 1115.68},
    ]
    fig = viz.plot_mass_torque_table(rows)
    ax = fig.axes[0]
    # 3 names * 2 bars
    assert len(ax.patches) == 6
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert labels == ["p", "n", "Lambda"]


def test_plot_mass_torque_table_engine_object(viz):
    class Engine:
        def predictions(self):
            return [
                {"name": "pi", "predicted": 139.0, "observed": 139.57},
                {"name": "K", "predicted": 494.0, "observed": 493.68},
            ]

    fig = viz.plot_mass_torque_table(Engine())
    assert len(fig.axes[0].patches) == 4


def test_plot_pmns_matrix_real(viz):
    # PDG-like |U|^2 row sums to ~1
    m = np.array(
        [[0.681, 0.297, 0.022],
         [0.144, 0.354, 0.502],
         [0.175, 0.349, 0.476]]
    )
    fig = viz.plot_pmns_matrix(m)
    img = fig.axes[0].images[0]
    assert img.get_array().shape == (3, 3)


def test_plot_pmns_matrix_complex(viz):
    rng = np.random.RandomState(1)
    u = rng.randn(3, 3) + 1j * rng.randn(3, 3)
    fig = viz.plot_pmns_matrix(u)
    arr = fig.axes[0].images[0].get_array()
    assert arr.shape == (3, 3)
    assert np.all(arr >= 0)


def test_plot_pmns_matrix_wrong_shape(viz):
    with pytest.raises(ValueError):
        viz.plot_pmns_matrix(np.eye(4))


def test_plot_cmb_acoustic_peaks_array(viz):
    cl = np.linspace(1e-10, 1e-9, 200)
    fig = viz.plot_cmb_acoustic_peaks(cl)
    line = fig.axes[0].lines[0]
    assert line.get_xdata().shape == (200,)
    assert fig.axes[0].get_xscale() == "log"


def test_plot_cmb_acoustic_peaks_mapping(viz):
    ell = np.arange(2, 50, dtype=float)
    cl = 1e-10 / (ell ** 1.5)
    fig = viz.plot_cmb_acoustic_peaks({"l": ell, "C_l": cl})
    assert fig.axes[0].lines[0].get_xdata().shape == ell.shape


def test_plot_falsifier_timeline_default(viz):
    fig = viz.plot_falsifier_timeline()
    ax = fig.axes[0]
    # one marker per default event
    assert len(ax.lines) >= len(SubstrateVisualizer.DEFAULT_FALSIFIERS)


def test_plot_falsifier_timeline_custom(viz):
    fig = viz.plot_falsifier_timeline([(2026, "X"), (2030, "Y")])
    assert "falsifiers" in viz.figures


def test_plot_substrate_potential_default(viz):
    fig = viz.plot_substrate_potential()
    line = fig.axes[0].lines[0]
    xs, ys = line.get_data()
    assert xs.size == 400
    # V >= 0 for the sine-Gordon potential
    assert np.all(ys >= -1e-12)
    # Maximum near 2K/xi^2 = 2
    assert ys.max() == pytest.approx(2.0, rel=1e-3)


def test_plot_substrate_potential_custom_array(viz):
    u = np.linspace(-1.0, 1.0, 50)
    fig = viz.plot_substrate_potential(u, K=2.0, xi=2.0)
    line = fig.axes[0].lines[0]
    assert line.get_xdata().size == 50


def test_save_all(tmp_path, viz):
    viz.plot_substrate_potential()
    viz.plot_dispersion((np.linspace(0, 1, 10), np.linspace(0, 1, 10)))
    viz.plot_field_2d(np.zeros((4, 4)), t=0.0)

    out = tmp_path / "plots"
    paths = viz.save_all(str(out))
    assert len(paths) == 3
    for p in paths:
        assert os.path.exists(p)
        assert p.endswith(".png")
        assert os.path.getsize(p) > 0


def test_save_all_creates_directory(tmp_path, viz):
    viz.plot_substrate_potential()
    target = tmp_path / "nested" / "deeper"
    paths = viz.save_all(str(target))
    assert target.exists()
    assert len(paths) == 1


def test_close_all_clears_cache(viz):
    viz.plot_substrate_potential()
    assert len(viz.figures) == 1
    viz.close_all()
    assert viz.figures == {}


def test_full_pipeline_smoke(tmp_path):
    """End-to-end: every plot method runs and saves."""
    v = SubstrateVisualizer()
    try:
        v.plot_field_2d(np.random.RandomState(2).randn(16, 16), t=0.1)
        v.plot_kink_evolution(
            [{"t": float(i), "x": float(i)} for i in range(5)]
        )
        v.plot_dispersion((np.arange(10), np.arange(10) ** 0.5))
        v.plot_mass_torque_table(
            [{"name": "a", "predicted": 1.0, "observed": 1.1}]
        )
        v.plot_pmns_matrix(np.eye(3) * 0.9 + 0.05)
        v.plot_cmb_acoustic_peaks(np.linspace(1e-10, 1e-9, 30))
        v.plot_falsifier_timeline()
        v.plot_substrate_potential()
        paths = v.save_all(str(tmp_path / "out"))
        assert len(paths) == 8
        for p in paths:
            assert os.path.exists(p) and os.path.getsize(p) > 0
    finally:
        v.close_all()
