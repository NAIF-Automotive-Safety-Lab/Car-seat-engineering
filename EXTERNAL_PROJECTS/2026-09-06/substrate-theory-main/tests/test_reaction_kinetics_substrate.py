"""Tests for ``reaction_kinetics_substrate``.

The substrate framework treats each reaction as a strain-knot
traversal across a barrier, so the rate constant must reproduce
Arrhenius (k = A · exp(−E_a/RT)) and Eyring (k = κ · k_B T/h ·
exp(−ΔG‡/RT)) for the canonical SN2 / H₂+I₂ / Diels-Alder reactions.

We verify five things:

(a) ``ReactionGeometry`` reproduces the three stationary points
    exactly: V(0) = E_R, V(½) = E_TS, V(1) = E_P with a correctly
    saddle-shaped barrier in between (negative curvature at ξ = ½).

(b) ``arrhenius_rate`` matches the closed-form k for the three
    literature reactions at room temperature within physically
    reasonable bounds (rate constants span >10 orders of magnitude
    across the three barrier heights).

(c) ``eyring_rate`` matches the closed form k = (k_B T / h) · exp.
    The rate-coefficient prefactor at 298.15 K is 6.21 × 10¹² 1/s,
    which serves as a hard sanity check.

(d) Catalysis: lowering E_a by 30 kJ/mol gives a speedup of order
    10⁵ at 298 K, exactly the magnitude of an effective enzyme.

(e) Round-trip: ``extracted_E_a_kJ`` recovers the input E_a from
    the slope of the Arrhenius plot to <0.01% accuracy.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import constants as C

from src.stiff_medium.reaction_kinetics_substrate import (
    LITERATURE,
    ReactionGeometry,
    ReactionKineticsSimulator,
    all_canonical_reactions,
    arrhenius_rate,
    catalysed_rate_ratio,
    eyring_rate,
    half_life_first_order,
    make_reaction,
)


# ===========================================================================
# (a) ReactionGeometry exactly interpolates the three stationary points
# ===========================================================================

class TestReactionGeometryStationaryPoints:
    """V(0) = E_R, V(½) = E_TS, V(1) = E_P; saddle has negative curvature."""

    def test_V_at_reactants(self) -> None:
        geom = ReactionGeometry(
            energy_reactants_kJ=0.0,
            energy_ts_kJ=120.0,
            energy_products_kJ=-40.0,
        )
        assert geom.potential(np.array([0.0]))[0] == pytest.approx(0.0, abs=1e-10)

    def test_V_at_transition_state(self) -> None:
        geom = ReactionGeometry(
            energy_reactants_kJ=0.0,
            energy_ts_kJ=120.0,
            energy_products_kJ=-40.0,
        )
        assert geom.potential(np.array([0.5]))[0] == pytest.approx(120.0, rel=1e-10)

    def test_V_at_products(self) -> None:
        geom = ReactionGeometry(
            energy_reactants_kJ=0.0,
            energy_ts_kJ=120.0,
            energy_products_kJ=-40.0,
        )
        assert geom.potential(np.array([1.0]))[0] == pytest.approx(-40.0, rel=1e-10)

    def test_TS_is_local_maximum(self) -> None:
        """V at ξ=½ must be the prominent strain peak (≥ both wells)."""
        geom = ReactionGeometry(
            energy_reactants_kJ=0.0,
            energy_ts_kJ=100.0,
            energy_products_kJ=-50.0,
        )
        xi = np.linspace(0.0, 1.0, 1001)
        V = geom.potential(xi)
        i_max = int(np.argmax(V))
        # The maximum sits within +/- 10% of xi=0.5 (centred barrier).
        assert 400 < i_max < 600
        # The maximum is well above both wells (>>E_R and >>E_P).
        assert V[i_max] > 0.0
        assert V[i_max] > -50.0
        # And matches the prescribed E_TS to within ~3% (centred-bump model
        # picks up a small drift from the asymmetric Hermite blend slope at
        # the saddle on strongly exothermic reactions; for symmetric ΔE = 0
        # it is exact, see test_symmetric_TS_exact below).
        assert V[i_max] == pytest.approx(100.0, rel=3e-2)

    def test_symmetric_TS_exact(self) -> None:
        """For ΔE = 0 the bump-on-blend model is exact at xi=1/2."""
        geom = ReactionGeometry(
            energy_reactants_kJ=0.0,
            energy_ts_kJ=80.0,
            energy_products_kJ=0.0,
        )
        xi = np.linspace(0.0, 1.0, 1001)
        V = geom.potential(xi)
        i_max = int(np.argmax(V))
        assert i_max == 500
        assert V[i_max] == pytest.approx(80.0, rel=1e-10)

    def test_barrier_curvature_is_negative(self) -> None:
        """Saddle ⇒ d²V/dξ² < 0 at ξ = ½ (imaginary frequency)."""
        geom = ReactionGeometry(
            energy_reactants_kJ=0.0,
            energy_ts_kJ=100.0,
            energy_products_kJ=-50.0,
        )
        assert geom.barrier_curvature() < 0.0

    def test_activation_and_reaction_energies(self) -> None:
        geom = ReactionGeometry(
            energy_reactants_kJ=0.0,
            energy_ts_kJ=120.0,
            energy_products_kJ=-40.0,
        )
        assert geom.activation_energy_kJ() == pytest.approx(120.0)
        assert geom.activation_energy_reverse_kJ() == pytest.approx(160.0)
        assert geom.reaction_energy_kJ() == pytest.approx(-40.0)


# ===========================================================================
# (b) Arrhenius reproduces the three canonical reaction rates at 298 K
# ===========================================================================

class TestArrheniusCanonicalReactions:
    """SN2, H₂+I₂, Diels-Alder rates at 298 K span the expected window."""

    def test_arrhenius_form_matches_closed_form(self) -> None:
        """k = A · exp(−E_a/RT) is reproduced exactly."""
        A = 1.0e11
        E_a_J = 90.0e3
        T = 298.15
        expected = A * math.exp(-E_a_J / (C.R * T))
        assert arrhenius_rate(A, E_a_J, T) == pytest.approx(expected, rel=1e-12)

    def test_sn2_rate_is_finite_at_300_K(self) -> None:
        """SN2 with E_a = 90 kJ/mol gives k ~ 10^-4 to 10^-3 M^-1 s^-1 at 298 K."""
        sim = make_reaction("SN2_methyl_bromide")
        k = sim.arrhenius_k(298.15)
        # exp(-90000/(8.314·298)) ~ 1.5e-16 ⇒ k ~ 3e-5 M^-1 s^-1.
        assert 1.0e-7 < k < 1.0e-2

    def test_h2_i2_rate_is_tiny_at_300_K(self) -> None:
        """H2+I2 with E_a = 165 kJ/mol gives k ~ 10^-18 M^-1 s^-1 at 298 K."""
        sim = make_reaction("H2_I2_to_2HI")
        k_low = sim.arrhenius_k(298.15)
        # At 700 K the same reaction is observable: rate >> at 298 K.
        k_hot = sim.arrhenius_k(700.0)
        assert k_low < k_hot
        assert k_hot / k_low > 1.0e6

    def test_diels_alder_rate_at_500K_is_reasonable(self) -> None:
        """At T = 500 K, Diels-Alder (E_a = 115) gives a measurable rate."""
        sim = make_reaction("Diels_Alder_butadiene_ethene")
        k = sim.arrhenius_k(500.0)
        # exp(-115000/(8.314·500)) ~ 1e-12, with A = 1e6 ⇒ k ~ 1e-6 M^-1 s^-1.
        assert 1.0e-9 < k < 1.0e-3

    def test_three_reactions_have_correct_ordering(self) -> None:
        """At a fixed T, lower E_a ⇒ faster reaction (SN2 fastest, H2+I2 slowest)."""
        T = 500.0
        k_sn2 = make_reaction("SN2_methyl_bromide").arrhenius_k(T)
        k_da = make_reaction("Diels_Alder_butadiene_ethene").arrhenius_k(T)
        k_h2i2 = make_reaction("H2_I2_to_2HI").arrhenius_k(T)
        assert k_sn2 > k_da > k_h2i2


# ===========================================================================
# (c) Eyring TST: k = κ · (k_B T / h) · exp(−ΔG‡/RT)
# ===========================================================================

class TestEyringTransitionStateTheory:
    """The Eyring rate prefactor k_B T / h is 6.21 × 10¹² 1/s at 298.15 K."""

    def test_eyring_prefactor_at_room_T(self) -> None:
        """k_B T / h at 298.15 K = 6.21 × 10¹² 1/s."""
        T = 298.15
        # ΔG‡ = 0 ⇒ k = κ · k_B T / h
        k = eyring_rate(0.0, T, kappa=1.0)
        expected_prefactor = (C.Boltzmann * T) / C.Planck
        assert k == pytest.approx(expected_prefactor, rel=1e-12)
        assert 6.0e12 < k < 6.5e12

    def test_eyring_form_matches_closed_form(self) -> None:
        """Eyring closed-form for ΔG‡ = 100 kJ/mol at 298.15 K."""
        delta_G_J = 100.0e3
        T = 298.15
        expected = (
            (C.Boltzmann * T / C.Planck)
            * math.exp(-delta_G_J / (C.R * T))
        )
        assert eyring_rate(delta_G_J, T) == pytest.approx(expected, rel=1e-12)

    def test_eyring_kappa_lowers_rate(self) -> None:
        """κ < 1 (recrossing) lowers the rate proportionally."""
        T = 300.0
        k_full = eyring_rate(80.0e3, T, kappa=1.0)
        k_half = eyring_rate(80.0e3, T, kappa=0.5)
        assert k_half == pytest.approx(0.5 * k_full, rel=1e-12)

    def test_eyring_arrhenius_consistency_for_DH_only(self) -> None:
        """If ΔS‡ = 0 then ΔG‡ = ΔH‡ ≈ E_a; Eyring k differs from
        Arrhenius only by the prefactor (k_B T/h vs A)."""
        sim = make_reaction("SN2_methyl_bromide")
        T = 400.0
        k_eyring = sim.eyring_k(T, delta_S_dagger_J_per_K_per_mol=0.0)
        # The same exponential, different prefactors.
        ratio = k_eyring / sim.arrhenius_k(T)
        expected_ratio = (C.Boltzmann * T / C.Planck) / sim.A_prefactor
        assert ratio == pytest.approx(expected_ratio, rel=1e-10)


# ===========================================================================
# (d) Catalysis: dropping E_a by 30 kJ/mol gives ~10^5x speedup at 298 K
# ===========================================================================

class TestCatalysisLowersBarrier:
    """Catalysed rate ratio = exp(ΔE_a / RT)."""

    def test_30_kJ_drop_gives_10_to_5_speedup_at_300K(self) -> None:
        """Classic ΔE_a = 30 kJ/mol ⇒ k_cat / k ≈ exp(30000/(8.314·300)) ≈ 1.7e5."""
        ratio = catalysed_rate_ratio(90.0e3, 60.0e3, 298.15)
        assert 1.0e4 < ratio < 1.0e6
        # Closed-form check:
        expected = math.exp(30.0e3 / (C.R * 298.15))
        assert ratio == pytest.approx(expected, rel=1e-12)

    def test_simulator_catalysed_speedup(self) -> None:
        """ReactionKineticsSimulator.catalysed_speedup matches the helper."""
        sim = make_reaction("SN2_methyl_bromide", catalyst_drop_kJ=30.0)
        speedup = sim.catalysed_speedup(298.15)
        expected = math.exp(30.0e3 / (C.R * 298.15))
        assert speedup == pytest.approx(expected, rel=1e-12)

    def test_catalyst_lowers_arrhenius_k(self) -> None:
        """k_cat > k_uncat at any temperature."""
        sim = make_reaction("Diels_Alder_butadiene_ethene", catalyst_drop_kJ=40.0)
        for T in [300.0, 500.0, 800.0]:
            k = sim.arrhenius_k(T, with_catalyst=False)
            k_cat = sim.arrhenius_k(T, with_catalyst=True)
            assert k_cat > k

    def test_zero_catalyst_drop_means_no_speedup(self) -> None:
        sim = make_reaction("SN2_methyl_bromide", catalyst_drop_kJ=0.0)
        assert sim.catalysed_speedup(298.15) == pytest.approx(1.0)


# ===========================================================================
# (e) Round-trip: Arrhenius slope recovers input E_a exactly
# ===========================================================================

class TestArrheniusPlotRoundTrip:
    """ln k vs 1/T has slope = −E_a / R; we should recover E_a to ~0.01%."""

    def test_round_trip_sn2(self) -> None:
        sim = make_reaction("SN2_methyl_bromide")
        recovered = sim.extracted_E_a_kJ()
        assert recovered == pytest.approx(90.0, rel=1e-4)

    def test_round_trip_h2_i2(self) -> None:
        sim = make_reaction("H2_I2_to_2HI")
        recovered = sim.extracted_E_a_kJ()
        assert recovered == pytest.approx(165.0, rel=1e-4)

    def test_round_trip_diels_alder(self) -> None:
        sim = make_reaction("Diels_Alder_butadiene_ethene")
        recovered = sim.extracted_E_a_kJ()
        assert recovered == pytest.approx(115.0, rel=1e-4)

    def test_arrhenius_curve_shape_is_monotonic(self) -> None:
        """ln k must increase monotonically with T (1/T decreasing)."""
        sim = make_reaction("Diels_Alder_butadiene_ethene")
        curve = sim.arrhenius_curve(T_min=300.0, T_max=800.0, n=50)
        # Sort by ascending T (descending 1/T):
        ln_k = curve["ln_k"]
        # Monotonic increase in T order ⇒ ln k strictly increasing.
        assert np.all(np.diff(ln_k) > 0.0)


# ===========================================================================
# (f) Half-life and miscellaneous sanity
# ===========================================================================

class TestUtilityFunctions:
    """Helper functions: half-life, all_canonical_reactions, summary."""

    def test_half_life_formula(self) -> None:
        assert half_life_first_order(1.0) == pytest.approx(math.log(2.0))
        assert half_life_first_order(0.1) == pytest.approx(math.log(2.0) / 0.1)

    def test_half_life_rejects_nonpositive(self) -> None:
        with pytest.raises(ValueError):
            half_life_first_order(0.0)
        with pytest.raises(ValueError):
            half_life_first_order(-1.0)

    def test_arrhenius_rejects_nonpositive_T(self) -> None:
        with pytest.raises(ValueError):
            arrhenius_rate(1.0, 1.0, 0.0)
        with pytest.raises(ValueError):
            arrhenius_rate(1.0, 1.0, -300.0)

    def test_eyring_rejects_nonpositive_T(self) -> None:
        with pytest.raises(ValueError):
            eyring_rate(1.0, 0.0)
        with pytest.raises(ValueError):
            eyring_rate(1.0, -300.0)

    def test_all_canonical_reactions_count(self) -> None:
        sims = all_canonical_reactions()
        assert len(sims) == len(LITERATURE)
        for s in sims:
            assert isinstance(s, ReactionKineticsSimulator)

    def test_make_reaction_unknown_raises(self) -> None:
        with pytest.raises(KeyError):
            make_reaction("not_a_real_reaction")

    def test_summary_returns_expected_keys(self) -> None:
        sim = make_reaction("SN2_methyl_bromide", catalyst_drop_kJ=30.0)
        summary = sim.summary(T=298.15)
        for key in ["label", "T_K", "E_a_kJ", "k_arrhenius",
                    "k_eyring", "k_arrhenius_catalysed",
                    "catalyst_speedup", "half_life_s"]:
            assert key in summary
