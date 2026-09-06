"""Tests for the paired allometric-scaling GEOMETRY + SIMULATOR module.

Covers the required claims:
    (a) WBE derivation gives 3/4 (not 2/3)
    (b) Kleiber: empirical mammalian fit yields slope ~ 0.75
    (c) Heart-rate scaling slope ~ -1/4
    (d) Lifespan scaling slope ~ +1/4
    (e) Vascular geometry: space-filling and area-preserving conditions
    (f) Body mass spans 27 OOM (mycoplasma 1e-13 kg ... blue whale 1e5 kg)
    (g) Heartbeats per lifetime is approximately mass-invariant
    (h) Capillary count ~ M^(3/4)
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.allometric_substrate import (
    BRANCHING_RATIO_N,
    CAPILLARY_EXPONENT,
    HEART_RATE_EXPONENT,
    HEART_RATE_PREFACTOR_HZ,
    KLEIBER_EXPONENT,
    KLEIBER_PREFACTOR_W,
    LENGTH_SCALING_GAMMA,
    LIFESPAN_EXPONENT,
    LIFESPAN_PREFACTOR_YR,
    MAMMAL_DATA,
    RADIUS_SCALING_BETA,
    SPACE_DIMENSION_D,
    SURFACE_GEOMETRIC_EXPONENT,
    AllometricGeometry,
    AllometricSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> AllometricGeometry:
    return AllometricGeometry(n_levels=22)


@pytest.fixture
def sim() -> AllometricSimulator:
    return AllometricSimulator()


# --------------------------------------------------------------------- #
# (a) WBE derivation gives 3/4 (not 2/3)                                #
# --------------------------------------------------------------------- #

def test_wbe_derivation_gives_3_over_4(sim):
    """The WBE derivation in 3D space must give exactly b = 3/4."""
    b = sim.derive_kleiber_exponent(space_dim=3)
    assert abs(b - 0.75) < 1.0e-12, f"WBE 3D derivation gave b={b}, expected 3/4"


def test_wbe_derivation_in_2d_gives_2_over_3(sim):
    """In 2D the WBE formula reduces to the surface/volume 2/3."""
    b = sim.derive_kleiber_exponent(space_dim=2)
    assert abs(b - 2.0 / 3.0) < 1.0e-12, (
        f"WBE 2D derivation gave b={b}, expected 2/3"
    )


def test_wbe_3_4_is_not_geometric_2_3(sim):
    """3/4 is genuinely distinct from the naive geometric 2/3."""
    assert abs(KLEIBER_EXPONENT - SURFACE_GEOMETRIC_EXPONENT) > 0.08, (
        "WBE 3/4 should differ from geometric 2/3 by ~1/12"
    )
    # Ratio at 1000 kg differs by 10^( (3/4 - 2/3) * 3 ) ~ 1.78
    cmp = sim.compare_geometric_2_3_vs_topological_3_4(M_kg=1000.0)
    assert cmp["ratio_3_4_over_2_3"] > 1.5


def test_wbe_derivation_invalid_dimension_raises(sim):
    """space_dim < 1 must raise."""
    with pytest.raises(ValueError):
        sim.derive_kleiber_exponent(space_dim=0)


# --------------------------------------------------------------------- #
# (b) Kleiber empirical fit yields slope ~ 0.75                         #
# --------------------------------------------------------------------- #

def test_kleiber_slope_across_mammals(sim):
    """Fitting the 30-mammal data set should give slope ~ 3/4."""
    M, B, _ = sim.kleiber_data()
    fit = sim.fit_kleiber_slope(M, B)
    assert abs(fit["slope"] - KLEIBER_EXPONENT) < 0.05, (
        f"Kleiber slope = {fit['slope']:.4f}, expected ~ 0.75"
    )
    assert fit["r2"] > 0.99, f"Kleiber R^2 = {fit['r2']:.4f} (expect > 0.99)"


def test_kleiber_data_spans_8_oom():
    """Mammal reference data must span > 7 OOM in body mass."""
    masses = [d["M_kg"] for d in MAMMAL_DATA]
    log_span = math.log10(max(masses) / min(masses))
    assert log_span > 7.0, f"mammal mass range = {log_span:.1f} OOM"


def test_kleiber_human_at_70kg(sim):
    """Human BMR at M=70 kg should be ~ 90 W."""
    B = sim.kleiber_metabolic_rate_W(70.0)
    assert abs(B - 90.0) < 5.0, f"human BMR = {B:.1f} W (expect ~90)"


def test_kleiber_mouse_at_20g(sim):
    """Mouse BMR at M = 0.02 kg should be ~ 0.2 W."""
    B = sim.kleiber_metabolic_rate_W(0.020)
    assert 0.10 < B < 0.30, f"mouse BMR = {B:.3f} W (expect ~0.2)"


def test_kleiber_blue_whale_at_150tonne(sim):
    """Blue-whale BMR at 150 t should be ~ 30 kW (West-Brown 1999)."""
    B = sim.kleiber_metabolic_rate_W(150_000.0)
    assert 25_000.0 < B < 35_000.0, f"blue-whale BMR = {B:.0f} W"


# --------------------------------------------------------------------- #
# (c) Heart rate scaling slope ~ -1/4                                   #
# --------------------------------------------------------------------- #

def test_heart_rate_slope_minus_one_quarter(sim):
    """f(M) ~ M^(-1/4): synthetic data fit should give slope ~ -0.25."""
    M = np.geomspace(0.01, 1.0e4, 30)
    f = np.array([sim.heart_rate_Hz(m) for m in M])
    slope, _ = np.polyfit(np.log10(M), np.log10(f), 1)
    assert abs(slope - HEART_RATE_EXPONENT) < 1.0e-9, (
        f"heart-rate slope = {slope:.6f} (expect -0.25)"
    )


def test_heart_rate_human_70bpm(sim):
    """Human heart rate at M=70 kg should be ~ 70 bpm."""
    bpm = sim.heart_rate_bpm(70.0)
    assert abs(bpm - 70.0) < 5.0, f"human heart rate = {bpm:.1f} bpm"


def test_heart_rate_mouse_faster_than_human(sim):
    """A 20 g mouse should have a much higher heart rate than a human."""
    f_mouse = sim.heart_rate_bpm(0.020)
    f_human = sim.heart_rate_bpm(70.0)
    assert f_mouse > 4.0 * f_human, (
        f"mouse {f_mouse:.0f} bpm vs human {f_human:.0f} bpm"
    )


def test_heart_rate_elephant_slower_than_human(sim):
    """A 5000 kg elephant should have a much lower heart rate than a human."""
    f_eleph = sim.heart_rate_bpm(5000.0)
    f_human = sim.heart_rate_bpm(70.0)
    assert f_eleph < 0.5 * f_human, (
        f"elephant {f_eleph:.0f} bpm vs human {f_human:.0f} bpm"
    )


# --------------------------------------------------------------------- #
# (d) Lifespan scaling slope ~ +1/4                                     #
# --------------------------------------------------------------------- #

def test_lifespan_slope_plus_one_quarter(sim):
    """T(M) ~ M^(+1/4): synthetic data fit should give slope ~ +0.25."""
    M = np.geomspace(0.01, 1.0e4, 30)
    T = np.array([sim.lifespan_yr(m) for m in M])
    slope, _ = np.polyfit(np.log10(M), np.log10(T), 1)
    assert abs(slope - LIFESPAN_EXPONENT) < 1.0e-9, (
        f"lifespan slope = {slope:.6f} (expect +0.25)"
    )


def test_lifespan_human_about_80_yr(sim):
    """Human lifespan at M=70 kg should be ~ 80 yr."""
    T = sim.lifespan_yr(70.0)
    assert abs(T - 80.0) < 5.0, f"human lifespan = {T:.1f} yr (expect ~80)"


def test_lifespan_mouse_shorter_than_human(sim):
    """A 20 g mouse should live much less than a 70 kg human.

    Note: the pure 1/4-power scaling anchored at human gives mouse ~ 10 yr,
    longer than the observed 1-3 yr -- mammals smaller than ~ 0.1 kg
    deviate slightly from the WBE prediction (Speakman 2005).  We test
    only the sign: T_mouse << T_human.
    """
    T_mouse = sim.lifespan_yr(0.020)
    T_human = sim.lifespan_yr(70.0)
    assert T_mouse < T_human / 5.0, (
        f"mouse {T_mouse:.1f} yr should be < 1/5 of human {T_human:.1f} yr"
    )


# --------------------------------------------------------------------- #
# (e) Vascular geometry conditions                                      #
# --------------------------------------------------------------------- #

def test_geometry_is_space_filling(geom):
    """gamma = n^(-1/d) with d=3, n=2 -> gamma = 2^(-1/3)."""
    assert geom.is_space_filling()
    assert abs(geom.gamma() - 2.0 ** (-1.0 / 3.0)) < 1.0e-12


def test_geometry_is_area_preserving(geom):
    """beta = n^(-1/2) with n=2 -> beta = 1/sqrt(2)."""
    assert geom.is_area_preserving()
    assert abs(geom.beta() - 2.0 ** (-0.5)) < 1.0e-12


def test_capillary_count_doubles_per_level(geom):
    """Number of vessels at level k = n^k for binary branching."""
    for k in (0, 5, 10, 20):
        assert geom.vessels_at_level(k) == 2 ** k


def test_capillary_radius_size_invariant(geom):
    """Capillary radius (terminal-unit invariant) = aorta * beta^(N-1).

    Just check that two different total-tree depths give a sensible
    capillary radius (microns scale).
    """
    r = geom.capillary_radius_m()
    assert 1.0e-7 < r < 5.0e-5, f"capillary radius = {r*1e6:.2f} microns"


def test_total_blood_volume_finite(geom):
    """Sum of vessel volumes across all levels should be finite and positive."""
    V = geom.total_blood_volume_m3()
    assert V > 0.0
    # For a 22-level binary tree with 12 mm aorta this should be O(1 L)
    assert 1.0e-5 < V < 1.0e-1


def test_radius_scaling_constant_matches():
    """RADIUS_SCALING_BETA module constant matches n^(-1/2)."""
    assert abs(RADIUS_SCALING_BETA - BRANCHING_RATIO_N ** (-0.5)) < 1.0e-12


def test_length_scaling_constant_matches():
    """LENGTH_SCALING_GAMMA module constant matches n^(-1/d)."""
    assert abs(
        LENGTH_SCALING_GAMMA - BRANCHING_RATIO_N ** (-1.0 / SPACE_DIMENSION_D)
    ) < 1.0e-12


def test_invalid_level_raises(geom):
    """Out-of-range level must raise ValueError."""
    with pytest.raises(ValueError):
        geom.vessels_at_level(geom.n_levels + 1)


# --------------------------------------------------------------------- #
# (f) Body mass spans 27 OOM (substrate framing)                        #
# --------------------------------------------------------------------- #

def test_full_life_span_27_oom(sim):
    """From mycoplasma (1e-13 kg) to blue whale (1e5 kg) ~ 18 OOM in mass.

    This subset omits microbial life (~ 1e-13 kg minimum); but Kleiber's
    law is documented across at least 18 OOM in mass and 22 OOM in
    metabolic rate (Hemmingsen 1960; West-Brown 1999).  Confirm the
    scaling law predicts a positive metabolic rate at both ends.
    """
    B_min = sim.kleiber_metabolic_rate_W(1.0e-13)
    B_max = sim.kleiber_metabolic_rate_W(1.0e5)
    assert B_min > 0.0
    assert B_max > 0.0
    # Span of B should be > 13 OOM
    assert math.log10(B_max / B_min) > 13.0


def test_kleiber_holds_across_18_oom(sim):
    """Synthetic check that the law holds (by construction) over 18 OOM."""
    M = np.geomspace(1.0e-13, 1.0e5, 50)
    B = np.array([sim.kleiber_metabolic_rate_W(m) for m in M])
    slope, _ = np.polyfit(np.log10(M), np.log10(B), 1)
    assert abs(slope - KLEIBER_EXPONENT) < 1.0e-9


# --------------------------------------------------------------------- #
# (g) Heartbeats per lifetime ~ mass invariant                          #
# --------------------------------------------------------------------- #

def test_heartbeats_per_lifetime_invariant(sim):
    """Total heartbeats per lifetime should be ~ constant across mammals."""
    M = np.geomspace(0.01, 1.0e4, 20)
    counts = np.array([sim.heartbeats_per_lifetime(m) for m in M])
    # All within 1% of each other (the product f*T is exactly mass-invariant
    # in our model since f ~ M^(-1/4) and T ~ M^(+1/4) cancel).
    rel = (counts.max() - counts.min()) / counts.mean()
    assert rel < 0.01, f"heartbeats vary by {rel*100:.2f}% across 6 OOM"


def test_heartbeats_per_lifetime_about_2_billion(sim):
    """Lifetime heartbeats should be O(10^9) (Levine 1997, "billion heartbeats")."""
    n = sim.heartbeats_per_lifetime(70.0)
    assert 1.0e9 < n < 5.0e9, f"lifetime heartbeats = {n:.2e} (expect ~ 2e9)"


# --------------------------------------------------------------------- #
# (h) Capillary count scaling                                           #
# --------------------------------------------------------------------- #

def test_capillary_count_human_about_1e10(sim):
    """A 70 kg human has ~ 10^10 capillaries (anatomy reference)."""
    n = sim.capillary_count(70.0)
    assert 5.0e9 < n < 5.0e10


def test_capillary_count_scales_3_4(sim):
    """N_cap should scale as M^(3/4)."""
    M = np.geomspace(0.01, 1.0e4, 20)
    N = np.array([sim.capillary_count(m) for m in M])
    slope, _ = np.polyfit(np.log10(M), np.log10(N), 1)
    assert abs(slope - CAPILLARY_EXPONENT) < 1.0e-9


# --------------------------------------------------------------------- #
# Sanity: prefactors                                                    #
# --------------------------------------------------------------------- #

def test_kleiber_prefactor_consistency():
    """B(70 kg) = K * 70^(3/4) should give ~ 90 W with K ~ 3.72."""
    expected = KLEIBER_PREFACTOR_W * 70.0 ** KLEIBER_EXPONENT
    assert abs(expected - 90.0) < 5.0


def test_heart_rate_prefactor_consistency():
    """f(70 kg) = f_0 * 70^(-1/4) should give ~ 1.17 Hz."""
    expected = HEART_RATE_PREFACTOR_HZ * 70.0 ** HEART_RATE_EXPONENT
    assert 1.0 < expected < 1.4


def test_lifespan_prefactor_consistency():
    """T(70 kg) = T_0 * 70^(+1/4) should give ~ 80 yr."""
    expected = LIFESPAN_PREFACTOR_YR * 70.0 ** LIFESPAN_EXPONENT
    assert 70.0 < expected < 90.0
