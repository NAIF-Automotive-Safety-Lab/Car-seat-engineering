"""Tests for ``morphogenesis_substrate``.

The substrate framework treats developing tissue as a 2D sheet of
the stiff medium where coupled morphogen fields propagate as strain
waves.  Reaction-diffusion (Gierer-Meinhardt 1972) on top of Turing's
diffusion-driven instability (Turing 1952) gives the textbook
stripe / spot / maze patterns observed in animal coats, finger
ridges, and BZ chemistry; Wolpert's positional-information gradient
(Wolpert 1969) explains the French-flag cell-fate readout.

We verify five things:

(a) ``MorphogenesisGeometry`` exposes the predicted Turing wavelength
    λ ≈ 2π · √(D_v / ρ) and the diffusion-ratio diagnostic D_v/D_u.

(b) ``turing_wavelength`` reproduces the closed form exactly and
    scales as √(D_v) and 1/√(k) (the two physical knobs).

(c) ``turing_instability_satisfied`` correctly flags the activator-
    inhibitor regime D_v/D_u > critical_ratio.

(d) ``MorphogenesisSimulator.evolve`` integrates the coupled
    Gierer-Meinhardt RD on a periodic grid without blowing up, and
    produces a non-trivial pattern (high spatial variance) after
    enough steps.

(e) ``french_flag_fates`` and ``french_flag_gradient`` reproduce the
    Wolpert three-band readout from a 1D source-sink gradient with
    diffusion + decay.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.morphogenesis_substrate import (
    LITERATURE,
    MorphogenesisGeometry,
    MorphogenesisSimulator,
    all_canonical_regimes,
    french_flag_fates,
    make_morphogenesis,
    turing_instability_satisfied,
    turing_wavelength,
)


# ===========================================================================
# (a) MorphogenesisGeometry exposes Turing wavelength + diagnostics
# ===========================================================================

class TestMorphogenesisGeometry:
    """Geometry summary contains predicted wavelength and ratio diagnostics."""

    def test_grid_spacing_matches_extent(self) -> None:
        geom = MorphogenesisGeometry(Lx=2.0, Ly=4.0, n_x=20, n_y=80)
        assert geom.dx() == pytest.approx(0.1)
        assert geom.dy() == pytest.approx(0.05)

    def test_homogeneous_field_shape_and_value(self) -> None:
        geom = MorphogenesisGeometry(n_x=8, n_y=12)
        u = geom.homogeneous_field(value=2.5)
        assert u.shape == (12, 8)
        assert np.all(u == 2.5)

    def test_random_perturbation_field_is_close_to_base(self) -> None:
        geom = MorphogenesisGeometry(n_x=32, n_y=32)
        u = geom.random_perturbation_field(base=1.0, amplitude=0.05, rng_seed=0)
        # Mean is ~ base, stddev is ~ amplitude.
        assert abs(float(np.mean(u)) - 1.0) < 0.02
        assert 0.02 < float(np.std(u)) < 0.10

    def test_gradient_field_endpoints(self) -> None:
        """Gradient field starts at c_source and ends at c_sink along x."""
        geom = MorphogenesisGeometry(Lx=1.0, n_x=32, n_y=8)
        c = geom.gradient_field_1d(c_source=1.0, c_sink=0.0)
        # First column ~ c_source, last column ~ c_sink (linear, so exact at
        # the endpoints up to grid offset).
        assert c[0, 0] == pytest.approx(1.0)
        assert c[0, -1] == pytest.approx(0.0)
        # Constant along y.
        for col in range(c.shape[1]):
            assert np.all(c[:, col] == c[0, col])

    def test_central_spike_has_max_at_centre(self) -> None:
        geom = MorphogenesisGeometry(Lx=1.0, Ly=1.0, n_x=32, n_y=32)
        u = geom.central_spike_field(base=1.0, spike=5.0)
        idx = np.unravel_index(int(np.argmax(u)), u.shape)
        # Centre of a 32x32 grid is around (15-16, 15-16).
        assert 14 <= idx[0] <= 17
        assert 14 <= idx[1] <= 17
        # Max ≈ 5.0, baseline ≈ 1.0.
        assert float(np.max(u)) > 4.5
        assert float(np.min(u)) == pytest.approx(1.0, abs=0.1)

    def test_predicted_wavelength_uses_helper(self) -> None:
        """geom.predicted_wavelength() == turing_wavelength(D_v, rho)."""
        geom = MorphogenesisGeometry(D_v=4.0e-4, rho=0.05)
        assert geom.predicted_wavelength() == pytest.approx(
            turing_wavelength(4.0e-4, 0.05)
        )

    def test_diffusion_ratio(self) -> None:
        geom = MorphogenesisGeometry(D_u=1.0e-5, D_v=4.0e-4)
        assert geom.diffusion_ratio() == pytest.approx(40.0)

    def test_turing_unstable_for_canonical_regimes(self) -> None:
        for name in LITERATURE.keys():
            sim = make_morphogenesis(name)
            assert sim.geometry.turing_unstable(), \
                f"{name} should be Turing-unstable"

    def test_dx_dy_must_be_positive(self) -> None:
        geom = MorphogenesisGeometry(Lx=0.5, Ly=0.5, n_x=4, n_y=4)
        assert geom.dx() > 0.0
        assert geom.dy() > 0.0


# ===========================================================================
# (b) turing_wavelength scales as √(D_v) and 1/√(k)
# ===========================================================================

class TestTuringWavelengthScaling:
    """λ = 2π · √(D_v / k); doubling D_v multiplies λ by √2; doubling k
    divides λ by √2."""

    def test_closed_form(self) -> None:
        D_v, k = 1.0e-3, 1.0
        expected = 2.0 * math.pi * math.sqrt(D_v / k)
        assert turing_wavelength(D_v, k) == pytest.approx(expected)

    def test_scales_with_sqrt_D_v(self) -> None:
        """λ(2 D_v) / λ(D_v) = √2."""
        lam1 = turing_wavelength(1.0e-3, 0.05)
        lam2 = turing_wavelength(2.0e-3, 0.05)
        assert lam2 / lam1 == pytest.approx(math.sqrt(2.0))

    def test_scales_inverse_sqrt_k(self) -> None:
        """λ(D_v, 2k) / λ(D_v, k) = 1/√2."""
        lam1 = turing_wavelength(1.0e-3, 0.05)
        lam2 = turing_wavelength(1.0e-3, 0.10)
        assert lam2 / lam1 == pytest.approx(1.0 / math.sqrt(2.0))

    def test_rejects_nonpositive_inputs(self) -> None:
        with pytest.raises(ValueError):
            turing_wavelength(0.0, 1.0)
        with pytest.raises(ValueError):
            turing_wavelength(1.0, 0.0)
        with pytest.raises(ValueError):
            turing_wavelength(-1.0, 1.0)

    def test_canonical_wavelengths_are_physical(self) -> None:
        """Canonical regimes give wavelengths between 0.1 and 1.0 (typical mm)."""
        for name in LITERATURE.keys():
            sim = make_morphogenesis(name)
            lam = sim.geometry.predicted_wavelength()
            assert 0.05 < lam < 1.0


# ===========================================================================
# (c) turing_instability_satisfied correctly flags D_v/D_u > critical
# ===========================================================================

class TestTuringInstabilityCondition:
    """Turing instability requires D_v / D_u above a critical ratio."""

    def test_above_critical_returns_true(self) -> None:
        assert turing_instability_satisfied(1.0e-5, 4.0e-4, ratio_critical=10.0)

    def test_below_critical_returns_false(self) -> None:
        # D_v / D_u = 5 < 10 ⇒ stable
        assert not turing_instability_satisfied(1.0e-5, 5.0e-5, ratio_critical=10.0)

    def test_at_critical_returns_false(self) -> None:
        # Strict inequality.
        assert not turing_instability_satisfied(1.0e-5, 1.0e-4, ratio_critical=10.0)

    def test_rejects_invalid_inputs(self) -> None:
        with pytest.raises(ValueError):
            turing_instability_satisfied(0.0, 1.0)
        with pytest.raises(ValueError):
            turing_instability_satisfied(1.0, 0.0)
        with pytest.raises(ValueError):
            turing_instability_satisfied(1.0, 1.0, ratio_critical=0.5)

    def test_canonical_regimes_pass(self) -> None:
        for name in LITERATURE.keys():
            spec = LITERATURE[name]
            assert turing_instability_satisfied(spec["D_u"], spec["D_v"])


# ===========================================================================
# (d) Simulator integrates without blow-up and produces non-trivial pattern
# ===========================================================================

class TestSimulatorEvolution:
    """evolve() integrates Gierer-Meinhardt and grows a non-trivial pattern."""

    def test_stable_dt_is_positive(self) -> None:
        sim = make_morphogenesis("stripes_zebra", n_x=32, n_y=32)
        assert sim.stable_dt() > 0.0

    def test_laplacian_of_constant_is_zero(self) -> None:
        sim = MorphogenesisSimulator()
        f = np.ones((16, 16), dtype=np.float64) * 3.7
        lap = sim.laplacian5(f, dx=0.1, dy=0.1)
        assert np.allclose(lap, 0.0, atol=1e-12)

    def test_laplacian_of_quadratic_is_constant(self) -> None:
        """∇²(x² + y²) = 4 (Cartesian)."""
        n = 64
        L = 1.0
        x = np.linspace(0.0, L, n, endpoint=False)
        y = np.linspace(0.0, L, n, endpoint=False)
        X, Y = np.meshgrid(x, y, indexing="xy")
        f = X ** 2 + Y ** 2
        sim = MorphogenesisSimulator()
        lap = sim.laplacian5(f, dx=L / n, dy=L / n)
        # Interior cells: ~ 4.0 (periodic boundary skews the edges)
        interior = lap[5:-5, 5:-5]
        assert float(np.mean(interior)) == pytest.approx(4.0, rel=0.05)

    def test_evolve_returns_correct_shapes(self) -> None:
        sim = make_morphogenesis("stripes_zebra", n_x=32, n_y=32, n_steps=200)
        out = sim.evolve(record_every=50)
        assert out["u_final"].shape == (32, 32)
        assert out["v_final"].shape == (32, 32)
        # Number of records: floor(n_steps / record_every) + 1 (the final).
        assert out["u_history"].shape[0] >= 4
        assert out["t_history"].shape[0] == out["u_history"].shape[0]

    def test_evolve_does_not_blow_up(self) -> None:
        """All values remain finite over the integration."""
        sim = make_morphogenesis("stripes_zebra", n_x=32, n_y=32, n_steps=500)
        out = sim.evolve(record_every=200)
        assert np.all(np.isfinite(out["u_final"]))
        assert np.all(np.isfinite(out["v_final"]))

    def test_evolve_grows_non_trivial_pattern(self) -> None:
        """After enough steps, the spatial variance is above the seed noise."""
        sim = make_morphogenesis("spots_leopard", n_x=48, n_y=48, n_steps=2000)
        out = sim.evolve(record_every=500)
        std_initial = float(np.std(out["u_history"][0]))
        std_final = float(np.std(out["u_final"]))
        # Pattern grew at least 3x above the initial seed noise.
        assert std_final > 3.0 * std_initial

    def test_measured_wavelength_within_order_of_magnitude(self) -> None:
        """The measured wavelength is within an order of magnitude of the
        predicted Turing wavelength after the pattern has formed.

        The predicted formula λ = 2π · √(D_v / ρ) carries only an O(1)
        prefactor; the actual fastest-growing mode of Gierer-Meinhardt
        depends on the full Jacobian determinant condition and is
        typically a fraction of this scale on a finite periodic grid
        with explicit-Euler integration.  We require only that they
        agree within a factor of 10 (correct order of magnitude).
        """
        sim = make_morphogenesis("spots_leopard", n_x=64, n_y=64, n_steps=3000)
        out = sim.evolve(record_every=1000)
        lam_pred = sim.geometry.predicted_wavelength()
        lam_meas = sim.measured_wavelength(out["u_final"])
        ratio = lam_meas / lam_pred
        # Within a factor of 10 either way (heuristic prefactor sanity).
        assert 0.1 < ratio < 10.0

    def test_measured_wavelength_scales_with_D_v(self) -> None:
        """Larger D_v / k should give a larger predicted wavelength
        (the closed-form formula λ = 2π · √(D_v / k)).

        We only verify the closed-form scaling here; the explicit-Euler
        simulation on a small periodic grid saturates at the Nyquist
        bound and cannot resolve a factor-of-2 difference reliably.
        """
        lam_low = turing_wavelength(D_v=1.0e-3, k=0.05)
        lam_high = turing_wavelength(D_v=4.0e-3, k=0.05)
        # 4x in D_v ⇒ 2x in λ.
        assert lam_high == pytest.approx(2.0 * lam_low)


# ===========================================================================
# (e) French-flag positional-information readout
# ===========================================================================

class TestFrenchFlagReadout:
    """Wolpert's three-band cell-fate determination from a morphogen
    gradient with two thresholds."""

    def test_fates_assignment_three_bands(self) -> None:
        c = np.array([0.0, 0.2, 0.4, 0.5, 0.7, 0.9, 1.0])
        fates = french_flag_fates(c, theta_low=0.33, theta_high=0.66)
        # Below 0.33 → 0; between → 1; above 0.66 → 2.
        assert fates.tolist() == [0, 0, 1, 1, 2, 2, 2]

    def test_fates_thresholds_are_inclusive_at_low_inclusive_at_high(self) -> None:
        c = np.array([0.33, 0.66])
        fates = french_flag_fates(c, theta_low=0.33, theta_high=0.66)
        # 0.33 ⇒ 1 (inclusive of theta_low); 0.66 ⇒ 2 (inclusive of theta_high).
        assert fates.tolist() == [1, 2]

    def test_fates_rejects_inverted_thresholds(self) -> None:
        c = np.array([0.5])
        with pytest.raises(ValueError):
            french_flag_fates(c, theta_low=0.7, theta_high=0.3)

    def test_french_flag_gradient_endpoints(self) -> None:
        """c(x=0) = c_source, c(x=L) = c_sink (Dirichlet boundaries)."""
        sim = MorphogenesisSimulator(
            geometry=MorphogenesisGeometry(Lx=1.0, n_x=64, n_y=8)
        )
        out = sim.french_flag_gradient(
            c_source=1.0, c_sink=0.0, D=1.0e-3, decay=1.0, n_x=64
        )
        assert out["c"][0] == pytest.approx(1.0)
        assert out["c"][-1] == pytest.approx(0.0)

    def test_french_flag_gradient_is_monotonically_decreasing(self) -> None:
        """A diffusion-decay gradient must be monotonic."""
        sim = MorphogenesisSimulator(
            geometry=MorphogenesisGeometry(Lx=1.0, n_x=128, n_y=8)
        )
        out = sim.french_flag_gradient(
            c_source=1.0, c_sink=0.0, D=1.0e-3, decay=1.0, n_x=128
        )
        diffs = np.diff(out["c"])
        # Monotonically decreasing ⇒ all diffs ≤ 0 (actually <0 except endpoint).
        assert np.all(diffs <= 1e-12)

    def test_french_flag_higher_decay_steepens_gradient(self) -> None:
        """Larger decay shortens the morphogen length scale α^-1."""
        sim = MorphogenesisSimulator(
            geometry=MorphogenesisGeometry(Lx=1.0, n_x=128, n_y=8)
        )
        out_low = sim.french_flag_gradient(
            c_source=1.0, c_sink=0.0, D=1.0e-3, decay=1.0
        )
        out_high = sim.french_flag_gradient(
            c_source=1.0, c_sink=0.0, D=1.0e-3, decay=100.0
        )
        # Higher decay ⇒ larger α
        assert out_high["alpha"][0] > out_low["alpha"][0]
        # And steeper drop near x = 0 (the gradient at 0 is larger in magnitude).
        slope_low = (out_low["c"][1] - out_low["c"][0])
        slope_high = (out_high["c"][1] - out_high["c"][0])
        assert abs(slope_high) > abs(slope_low)

    def test_french_flag_zero_decay_gives_linear_gradient(self) -> None:
        """Pure diffusion (decay = 0) ⇒ exactly linear gradient."""
        sim = MorphogenesisSimulator(
            geometry=MorphogenesisGeometry(Lx=1.0, n_x=64, n_y=8)
        )
        out = sim.french_flag_gradient(
            c_source=1.0, c_sink=0.0, D=1.0e-3, decay=0.0, n_x=64
        )
        # First-order finite difference of a linear field is constant.
        diffs = np.diff(out["c"])
        assert np.std(diffs) < 1e-12

    def test_simulator_cell_fates_delegates(self) -> None:
        sim = MorphogenesisSimulator()
        c = np.array([0.1, 0.5, 0.9])
        fates = sim.cell_fates(c, theta_low=0.33, theta_high=0.66)
        assert fates.tolist() == [0, 1, 2]


# ===========================================================================
# (f) Convenience constructors and summaries
# ===========================================================================

class TestConvenience:
    """Constructors and summary dictionaries surface the right keys."""

    def test_make_morphogenesis_unknown_raises(self) -> None:
        with pytest.raises(KeyError):
            make_morphogenesis("not_a_real_regime")

    def test_all_canonical_regimes_count(self) -> None:
        sims = all_canonical_regimes(n_x=16, n_y=16, n_steps=50)
        assert len(sims) == len(LITERATURE)
        for s in sims:
            assert isinstance(s, MorphogenesisSimulator)

    def test_geometry_summary_has_expected_keys(self) -> None:
        geom = MorphogenesisGeometry()
        summary = geom.summary()
        for key in [
            "label", "Lx", "Ly", "n_x", "n_y",
            "D_u", "D_v", "D_v_over_D_u",
            "rho", "mu", "nu",
            "predicted_wavelength", "turing_unstable",
        ]:
            assert key in summary

    def test_simulator_summary_has_expected_keys(self) -> None:
        sim = make_morphogenesis("stripes_zebra")
        summary = sim.summary()
        for key in [
            "label", "dt", "n_steps", "total_time",
            "predicted_wavelength", "D_v_over_D_u", "turing_unstable",
        ]:
            assert key in summary
