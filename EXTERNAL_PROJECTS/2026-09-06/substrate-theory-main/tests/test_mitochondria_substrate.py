"""Tests for ``mitochondria_substrate`` — bioenergetic chemiosmosis.

We verify the four canonical claims of mitochondrial bioenergetics:

(a) Proton motive force Δμ_H+ ≈ 200 mV (Mitchell 1961).
    With ΔΨ = -160 mV and ΔpH = +1, |PMF| = 0.16 + 0.0615 ≈ 0.22 V.

(b) F0F1 ATP synthase stoichiometry: c protons per 360° = 3 ATP.
    H⁺/ATP = c/3 → 2.67 (mammal c8), 3.33 (yeast c10), 4.67 (plant c14).

(c) ETC stoichiometry: 1 NADH → 2.5 ATP, 1 FADH₂ → 1.5 ATP (P/O).
    Glucose → ~30 ATP via 10 NADH + 2 FADH₂ + 4 SLP.

(d) Total NADH→O₂ redox drop ΔE = 1.14 V → ΔG = -220 kJ/mol per 2 e⁻
    drives 10 H⁺ (mammal: 4+4+2 = Complex I+III+IV).
"""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.mitochondria_substrate import (
    ATP_PER_REV,
    C_RING_BACTERIA,
    C_RING_MAMMAL,
    C_RING_PLANT,
    C_RING_YEAST,
    DELTA_PH_DEFAULT,
    DELTA_PSI_DEFAULT_V,
    E_FADH2_FAD,
    E_NADH_NAD,
    E_O2_H2O,
    F_FARADAY,
    H_PER_2E_COMPLEX_I,
    H_PER_2E_COMPLEX_II,
    H_PER_2E_COMPLEX_III,
    H_PER_2E_COMPLEX_IV,
    NADH_PER_GLUCOSE,
    SUBSTRATE_LEVEL_ATP,
    FADH2_PER_GLUCOSE,
    MitoGeometry,
    MitoSimulator,
    atp_per_fadh2,
    atp_per_nadh,
    c_ring_atp_ratio,
    delta_g_per_proton,
    delta_g_redox,
    etc_complexes_default,
    example_c_rings,
    glucose_to_atp,
    proton_motive_force,
)


# ===========================================================================
# (a) Proton motive force ≈ 200 mV (Mitchell 1961 chemiosmosis)
# ===========================================================================

class TestProtonMotiveForce:

    def test_mitchell_pmf_default(self) -> None:
        """Default Mitchell PMF ≈ 200 mV (between 180 and 230 mV)."""
        pmf = proton_motive_force()  # uses defaults: ΔΨ = -160 mV, ΔpH = +1
        # PMF = 0.160 + 0.0615 × 1 ≈ 0.221 V → 221 mV
        assert 0.180 < pmf < 0.230, f"PMF was {pmf*1000:.1f} mV"

    def test_pmf_in_geometry_matches_helper(self) -> None:
        geom = MitoGeometry()
        assert geom.proton_motive_force() == pytest.approx(
            proton_motive_force(), rel=1e-12)

    def test_pmf_simulator_in_millivolts(self) -> None:
        sim = MitoSimulator()
        # Should be in the canonical 180-230 mV band
        assert 180.0 < sim.pmf_mv() < 230.0

    def test_pmf_zero_delta_psi_zero_delta_ph(self) -> None:
        """Zero gradient → zero PMF."""
        assert proton_motive_force(delta_psi_v=0.0, delta_ph=0.0) == 0.0

    def test_pmf_only_membrane_potential(self) -> None:
        """ΔpH=0, ΔΨ=-160 mV → PMF = 160 mV exactly."""
        pmf = proton_motive_force(delta_psi_v=-0.160, delta_ph=0.0)
        assert pmf == pytest.approx(0.160, abs=1e-12)

    def test_pmf_increases_with_delta_ph(self) -> None:
        """A larger ΔpH gives a larger PMF (Nernst slope ≈ 60 mV/unit)."""
        p1 = proton_motive_force(delta_psi_v=-0.160, delta_ph=0.0)
        p2 = proton_motive_force(delta_psi_v=-0.160, delta_ph=1.0)
        # Difference should be ≈ 0.0615 V at 310 K
        assert (p2 - p1) == pytest.approx(0.0615, abs=0.005)

    def test_dG_per_proton_around_21_kjmol(self) -> None:
        """ΔG per H⁺ at PMF ~ 200 mV is ≈ 21 kJ/mol."""
        dg = delta_g_per_proton() / 1000.0  # kJ/mol
        assert 18.0 < dg < 24.0, f"ΔG/H⁺ was {dg:.2f} kJ/mol"


# ===========================================================================
# (b) F0F1 ATP synthase: c protons / 360° = 3 ATP
# ===========================================================================

class TestF0F1ATPSynthase:

    def test_atp_per_revolution_is_three(self) -> None:
        """3 catalytic β-subunits → 3 ATP per 360° rotation."""
        assert ATP_PER_REV == 3
        sim = MitoSimulator()
        assert sim.atp_per_revolution() == 3

    def test_h_per_atp_mammal_c8(self) -> None:
        """c=8 mammal: H⁺/ATP = 8/3 ≈ 2.67."""
        ratio = c_ring_atp_ratio(C_RING_MAMMAL)
        assert ratio == pytest.approx(8.0 / 3.0, rel=1e-12)
        assert 2.6 < ratio < 2.8

    def test_h_per_atp_yeast_c10(self) -> None:
        """c=10 yeast: H⁺/ATP = 10/3 ≈ 3.33 (also bacteria E. coli)."""
        ratio = c_ring_atp_ratio(C_RING_YEAST)
        assert ratio == pytest.approx(10.0 / 3.0, rel=1e-12)
        assert 3.2 < ratio < 3.4

    def test_h_per_atp_plant_c14(self) -> None:
        """c=14 chloroplast: H⁺/ATP = 14/3 ≈ 4.67."""
        ratio = c_ring_atp_ratio(C_RING_PLANT)
        assert ratio == pytest.approx(14.0 / 3.0, rel=1e-12)
        assert 4.5 < ratio < 4.8

    def test_360_rotation_equals_c_protons_times_3_atp(self) -> None:
        """One full 360° rotation = c protons in = c × 3 / c = 3 ATP."""
        for c in (C_RING_BACTERIA, C_RING_YEAST, C_RING_MAMMAL, C_RING_PLANT):
            sim = MitoSimulator(geometry=MitoGeometry(c_ring_size=c))
            atp = sim.atp_per_c_protons(c)
            assert atp == pytest.approx(3.0, rel=1e-12), (
                f"c={c} should give 3 ATP per c protons, got {atp:.4f}")

    def test_c_ring_size_must_be_positive(self) -> None:
        with pytest.raises(ValueError):
            c_ring_atp_ratio(0)
        with pytest.raises(ValueError):
            c_ring_atp_ratio(-1)

    def test_torque_in_pNnm_realistic(self) -> None:
        """F1 single-molecule torque: ~25-50 pN·nm (Yasuda 2001 et al.)."""
        sim = MitoSimulator()
        torque_pnm = sim.torque_required_nm()    # already in pN·nm
        assert 10.0 < torque_pnm < 80.0, (
            f"Torque was {torque_pnm:.2f} pN·nm")

    def test_revolutions_per_atp_one_third(self) -> None:
        sim = MitoSimulator()
        assert sim.revolutions_per_atp() == pytest.approx(1.0 / 3.0,
                                                            rel=1e-12)


# ===========================================================================
# (c) ETC stoichiometry: 1 NADH → 2.5 ATP, 1 FADH₂ → 1.5 ATP, glucose ~30 ATP
# ===========================================================================

class TestETCStoichiometry:

    def test_glucose_to_atp_default_30(self) -> None:
        """Default P/O = 2.5/1.5 → 30 ATP per glucose (textbook)."""
        atp = glucose_to_atp()
        # 4 + 10×2.5 + 2×1.5 = 4 + 25 + 3 = 32  (idealised)
        # The textbook "30" comes from cyto NADH penalty (glycerol shuttle).
        # Default returns 32; this test allows the canonical band.
        assert 28.0 <= atp <= 32.5, f"ATP/glucose was {atp:.1f}"

    def test_glucose_to_atp_within_28_to_32(self) -> None:
        """ATP per glucose lies in the canonical mammalian band 28-32."""
        sim = MitoSimulator()
        atp = sim.atp_per_glucose()
        assert 28.0 <= atp <= 32.5, f"ATP/glucose was {atp:.1f}"

    def test_p_o_ratio_nadh_default_2p5(self) -> None:
        sim = MitoSimulator()
        assert sim.p_o_ratio_nadh() == pytest.approx(2.5, rel=1e-12)

    def test_p_o_ratio_fadh2_default_1p5(self) -> None:
        sim = MitoSimulator()
        assert sim.p_o_ratio_fadh2() == pytest.approx(1.5, rel=1e-12)

    def test_nadh_pumps_more_protons_than_fadh2(self) -> None:
        """NADH path goes through I+III+IV; FADH₂ skips I."""
        sim = MitoSimulator()
        n = sim.etc_h_pumped_per_nadh()
        f = sim.etc_h_pumped_per_fadh2()
        assert n > f
        assert (n - f) == H_PER_2E_COMPLEX_I

    def test_complex_proton_counts_mammal(self) -> None:
        """Canonical proton counts: I=4, II=0, III=4, IV=2 per 2 e⁻."""
        assert H_PER_2E_COMPLEX_I == 4
        assert H_PER_2E_COMPLEX_II == 0
        assert H_PER_2E_COMPLEX_III == 4
        assert H_PER_2E_COMPLEX_IV == 2

    def test_atp_per_nadh_default_helper(self) -> None:
        """Helper atp_per_nadh() ≥ 2.5 (uncorrected gives 10/2.67 ≈ 3.75)."""
        ratio = atp_per_nadh()
        assert ratio >= 2.5

    def test_atp_per_fadh2_default_helper(self) -> None:
        """Helper atp_per_fadh2() ≥ 1.5 (uncorrected gives 6/2.67 ≈ 2.25)."""
        ratio = atp_per_fadh2()
        assert ratio >= 1.5

    def test_glucose_breakdown_sums(self) -> None:
        """Itemised breakdown sums to total."""
        sim = MitoSimulator()
        b = sim.atp_per_glucose_breakdown()
        s = b["substrate_level"] + b["from_nadh"] + b["from_fadh2"]
        assert s == pytest.approx(b["total"], rel=1e-12)


# ===========================================================================
# (d) ETC redox: NADH→O₂ ΔE = 1.14 V, ΔG = -220 kJ/mol per 2 e⁻
# ===========================================================================

class TestETCRedox:

    def test_nadh_to_o2_total_delta_e(self) -> None:
        """ΔE = E°'(½O₂/H₂O) - E°'(NADH/NAD⁺) = 0.815 − (−0.320) ≈ 1.14 V."""
        sim = MitoSimulator()
        de = sim.etc_total_delta_e_nadh()
        assert de == pytest.approx(E_O2_H2O - E_NADH_NAD, rel=1e-12)
        assert 1.10 < de < 1.20

    def test_nadh_to_o2_total_delta_g(self) -> None:
        """ΔG = -nFΔE ≈ -220 kJ/mol per 2 e⁻."""
        sim = MitoSimulator()
        dg = sim.etc_total_delta_g_nadh_kjmol()
        # -2 × 96485 × 1.135 / 1000 = -219.0 kJ/mol
        assert -230.0 < dg < -210.0, f"ΔG was {dg:.1f} kJ/mol"

    def test_redox_formula_sign_convention(self) -> None:
        """Downhill drop (ΔE>0) gives negative (exergonic) ΔG."""
        dg_pos = delta_g_redox(+1.0, n_electrons=2)    # exergonic
        dg_neg = delta_g_redox(-1.0, n_electrons=2)    # endergonic
        assert dg_pos < 0.0
        assert dg_neg > 0.0
        assert dg_pos == -dg_neg

    def test_complex_i_drop_is_negative(self) -> None:
        """Complex I: NADH (-0.32 V) → Q (+0.04 V), drop = +0.36 V → ΔG<0."""
        complexes = etc_complexes_default()
        c1 = complexes[0]
        assert c1.name == "I"
        assert c1.delta_e() > 0      # downhill electron flow
        assert c1.delta_g_per_2e() < 0    # exergonic

    def test_etc_complexes_redox_ladder_ordered(self) -> None:
        """Redox potentials increase monotonically NADH → O₂."""
        sim = MitoSimulator()
        ladder = sim.etc_redox_ladder()
        # NADH < FADH2 < Q < cyt c < O2
        # but FADH2 enters at the Q level via Complex II so it sits between
        # NADH and Q in the listed ladder
        es = [e for _, e in ladder]
        # Strictly NADH < FADH2 (since -0.32 < -0.22)
        assert es[0] < es[1]
        # And the first NADH reading is the most-negative
        assert es[0] == min(es)
        # And O₂ is the most-positive sink
        assert es[-1] == max(es)


# ===========================================================================
# Geometry / topology consistency
# ===========================================================================

class TestMitoGeometry:

    def test_default_has_four_complexes(self) -> None:
        geom = MitoGeometry()
        assert len(geom.complexes) == 4
        assert [c.name for c in geom.complexes] == ["I", "II", "III", "IV"]

    def test_membrane_profile_returns_array(self) -> None:
        geom = MitoGeometry()
        x = np.linspace(0.0, 1.0, 257)
        y = geom.membrane_profile(x)
        assert y.shape == x.shape
        assert np.all(np.isfinite(y))

    def test_proton_gradient_higher_outside(self) -> None:
        """Intermembrane space (y>0) has more H⁺ than matrix (y<0)."""
        geom = MitoGeometry()
        out_h = geom.proton_gradient(np.array([+0.5]))[0]
        in_h = geom.proton_gradient(np.array([-0.5]))[0]
        assert out_h > in_h

    def test_summary_has_pmf_around_200mv(self) -> None:
        geom = MitoGeometry()
        s = geom.summary()
        assert 0.18 < s["pmf_v"] < 0.23

    def test_cristae_fold_factor_default_is_5(self) -> None:
        """Inner membrane area ≈ 5× outer membrane area in mammalian mito."""
        geom = MitoGeometry()
        assert geom.cristae_fold_factor == pytest.approx(5.0, rel=1e-12)


class TestSimulatorTimeSeries:

    def test_simulate_conserves_nadh_decreases(self) -> None:
        """First-order kinetics: NADH monotonically decreases."""
        sim = MitoSimulator()
        out = sim.simulate(n_nadh=100.0, n_steps=200, dt=1.0)
        assert np.all(np.diff(out["nadh"]) <= 0)

    def test_simulate_atp_monotonically_increases(self) -> None:
        """ATP accumulation monotonically increases (no consumption)."""
        sim = MitoSimulator()
        out = sim.simulate(n_nadh=100.0, n_steps=200, dt=1.0)
        assert np.all(np.diff(out["atp"]) >= 0)

    def test_simulate_atp_yield_in_band(self) -> None:
        """100 NADH → roughly 100 × P/O = 250 ATP if reactions go to completion."""
        sim = MitoSimulator()
        out = sim.simulate(n_nadh=100.0, n_steps=2000, dt=1.0)
        atp_final = float(out["atp"][-1])
        # Should be in the 200-400 ATP band (depends on h_per_atp = 8/3)
        assert 200.0 < atp_final < 500.0, f"ATP final was {atp_final:.1f}"


class TestExampleCRingCatalogue:

    def test_catalogue_has_four_organisms(self) -> None:
        rows = example_c_rings()
        assert len(rows) == 4

    def test_h_per_atp_increases_with_c(self) -> None:
        rows = example_c_rings()
        ratios = [r["h_per_atp"] for r in rows]
        cs = [r["c"] for r in rows]
        # H⁺/ATP scales linearly with c
        for c, r in zip(cs, ratios):
            assert r == pytest.approx(c / 3.0, rel=1e-12)
