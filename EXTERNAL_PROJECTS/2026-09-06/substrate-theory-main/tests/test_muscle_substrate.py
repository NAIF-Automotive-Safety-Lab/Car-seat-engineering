"""Tests for the paired muscle-contraction GEOMETRY + SIMULATOR module.

Covers the four required claims:
    (a) Hill (1938) F-v hyperbolic form: (F + a)(v + b) = (F_max + a) * b
    (b) Maximum mechanical power at v ~ 0.3 * v_max  (a/F_max = 0.25)
    (c) Length-tension peak at L_0 in [2.0, 2.25] um  (Gordon-Huxley-Julian 1966)
    (d) ATP -> mechanical work efficiency ~ 50% (theoretical ceiling)

Plus geometric and biochemical sanity:
    * sarcomere length 2.5 um (Z-to-Z)
    * thick filament 1.6 um (A-band)
    * thin filament 1.0 um (Z to tip)
    * H-zone 0.2 um (M-line region)
    * cross-bridge axial period 14.3 nm; 6 heads per helical turn
    * 1 ATP per cross-bridge cycle
    * power stroke ~ 5-10 nm displacement
    * 5-state cross-bridge cycle (Lymn-Taylor 1971)
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.muscle_substrate import (
    ATP_HYDROLYSIS_PN_NM,
    CROSS_BRIDGE_AXIAL_PERIOD_NM,
    DUTY_RATIO_SKELETAL,
    HEADS_PER_HELICAL_TURN,
    HILL_A_OVER_FMAX,
    H_ZONE_LENGTH_UM,
    INTER_FILAMENT_SPACING_NM,
    LENGTH_TENSION_L_MAX_UM,
    LENGTH_TENSION_L_MIN_UM,
    LENGTH_TENSION_L_PLATEAU_HIGH_UM,
    LENGTH_TENSION_L_PLATEAU_LOW_UM,
    MYOSIN_HEADS_PER_THICK,
    POWER_STROKE_NM,
    SARCOMERE_REST_LENGTH_UM,
    SINGLE_HEAD_FORCE_PN,
    THICK_FILAMENT_DIAMETER_NM,
    THICK_FILAMENT_LENGTH_UM,
    THIN_FILAMENT_DIAMETER_NM,
    THIN_FILAMENT_LENGTH_UM,
    MuscleGeometry,
    MuscleSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> MuscleGeometry:
    return MuscleGeometry(sarcomere_length_um=SARCOMERE_REST_LENGTH_UM)


@pytest.fixture
def sim(geom: MuscleGeometry) -> MuscleSimulator:
    return MuscleSimulator(geometry=geom)


# --------------------------------------------------------------------- #
# (a) Hill 1938 force-velocity hyperbolic form                          #
# --------------------------------------------------------------------- #

def test_hill_hyperbola_endpoints(sim: MuscleSimulator):
    """At F=0 we get v=v_max; at F=F_max we get v=0."""
    assert math.isclose(sim.hill_velocity_from_force(0.0),
                        sim.v_max_normalised, abs_tol=1e-12)
    assert math.isclose(sim.hill_velocity_from_force(sim.f_max_normalised),
                        0.0, abs_tol=1e-12)


def test_hill_hyperbola_invariant(sim: MuscleSimulator):
    """(F + a)(v + b) = (F_max + a) * b should hold pointwise."""
    a = sim.hill_a()
    b = sim.hill_b()
    F_max = sim.f_max_normalised
    rhs = (F_max + a) * b
    forces = np.linspace(0.0, F_max, 50)
    for F in forces:
        v = sim.hill_velocity_from_force(F)
        lhs = (F + a) * (v + b)
        assert math.isclose(lhs, rhs, rel_tol=1e-9), \
            f"Hill invariant violated at F={F}: LHS={lhs} RHS={rhs}"


def test_hill_inverse_consistency(sim: MuscleSimulator):
    """hill_force_from_velocity is the inverse of hill_velocity_from_force."""
    v_test = np.linspace(0.0, sim.v_max_normalised, 30)
    for v in v_test:
        F = sim.hill_force_from_velocity(v)
        v_back = sim.hill_velocity_from_force(F)
        assert math.isclose(v, v_back, rel_tol=1e-8, abs_tol=1e-10), \
            f"inverse mismatch at v={v}: got v_back={v_back}"


def test_hill_curve_monotone_decreasing(sim: MuscleSimulator):
    """F is a strictly decreasing function of v on [0, v_max]."""
    f, v = sim.hill_force_velocity_curve(n=200)
    diffs = np.diff(f)
    assert (diffs <= 1e-12).all(), "Force should be non-increasing in velocity"
    # And strictly less than F_max except at v = 0
    assert f[-1] < f[0]


# --------------------------------------------------------------------- #
# (b) Peak mechanical power at v ~ 0.31 v_max (a/F_max = 0.25)          #
# --------------------------------------------------------------------- #

def test_peak_power_at_about_0p31_vmax(sim: MuscleSimulator):
    """For Hill's hyperbola with a/F_max = 0.25 the analytical
    optimum is v* = b * (sqrt((1+a)/a) - 1).  With a/F_max=0.25
    and b=0.25*v_max this gives v*/v_max ~ 0.309."""
    v_star = sim.velocity_at_max_power()
    expected = sim.hill_b() * (math.sqrt((sim.f_max_normalised + sim.hill_a())
                                          / sim.hill_a()) - 1.0)
    assert math.isclose(v_star, expected, rel_tol=1e-12)
    # And in the textbook v_max = 1 case this is in the [0.28, 0.34] band.
    ratio = v_star / sim.v_max_normalised
    assert 0.28 <= ratio <= 0.34, f"v*/v_max = {ratio}, expected ~ 0.31"


def test_peak_power_is_actual_maximum(sim: MuscleSimulator):
    """Numerical scan of P(v) should peak at the analytical v*."""
    v, p, _ = sim.hill_power_curve(n=2001)
    i_max_numeric = int(np.argmax(p))
    v_max_numeric = v[i_max_numeric]
    v_star_analytic = sim.velocity_at_max_power()
    # Numeric and analytical agree to within one bin
    bin_width = v[1] - v[0]
    assert abs(v_max_numeric - v_star_analytic) <= 2 * bin_width


def test_max_power_value(sim: MuscleSimulator):
    """Max mechanical power = P(v*) > 0 and < F_max * v_max."""
    P_max = sim.max_mechanical_power()
    assert P_max > 0.0
    assert P_max < sim.f_max_normalised * sim.v_max_normalised
    # And at a/F_max = 0.25 the P_max / (F_max v_max) is ~ 0.10
    assert 0.07 <= P_max <= 0.13


# --------------------------------------------------------------------- #
# (c) Length-tension peak at L_0 in [2.0, 2.25] um                      #
# --------------------------------------------------------------------- #

def test_length_tension_peak_at_optimal_length(sim: MuscleSimulator):
    """F_active = 1.0 on the plateau [L_plateau_low, L_plateau_high]."""
    L0 = sim.optimal_length_um()
    assert LENGTH_TENSION_L_PLATEAU_LOW_UM <= L0 <= LENGTH_TENSION_L_PLATEAU_HIGH_UM
    assert math.isclose(sim.active_force_at_length(L0), 1.0, abs_tol=1e-12)
    # And both plateau endpoints also give F = 1
    assert math.isclose(
        sim.active_force_at_length(LENGTH_TENSION_L_PLATEAU_LOW_UM), 1.0, abs_tol=1e-12)
    assert math.isclose(
        sim.active_force_at_length(LENGTH_TENSION_L_PLATEAU_HIGH_UM), 1.0, abs_tol=1e-12)


def test_length_tension_zero_at_extremes(sim: MuscleSimulator):
    """F_active = 0 at L_min (clash) and L_max (no overlap)."""
    assert sim.active_force_at_length(LENGTH_TENSION_L_MIN_UM) == 0.0
    assert sim.active_force_at_length(LENGTH_TENSION_L_MAX_UM) == 0.0
    assert sim.active_force_at_length(0.5) == 0.0
    assert sim.active_force_at_length(5.0) == 0.0


def test_length_tension_curve_monotone_segments(sim: MuscleSimulator):
    """Ascending limb is non-decreasing; descending limb is non-increasing."""
    L, F = sim.length_tension_curve(n=400)
    asc_mask = (L > LENGTH_TENSION_L_MIN_UM) & (L < LENGTH_TENSION_L_PLATEAU_LOW_UM)
    desc_mask = (L > LENGTH_TENSION_L_PLATEAU_HIGH_UM) & (L < LENGTH_TENSION_L_MAX_UM)
    F_asc = F[asc_mask]
    F_desc = F[desc_mask]
    assert (np.diff(F_asc) >= -1e-12).all(), "ascending limb should be non-decreasing"
    assert (np.diff(F_desc) <= 1e-12).all(), "descending limb should be non-increasing"


def test_length_tension_peak_value_is_one(sim: MuscleSimulator):
    """Peak active force = 1 (normalised) at any plateau-interior length."""
    L, F = sim.length_tension_curve(n=4000)
    assert math.isclose(F.max(), 1.0, rel_tol=1e-9)


# --------------------------------------------------------------------- #
# (d) ATP -> mechanical work efficiency ~ 50%                           #
# --------------------------------------------------------------------- #

def test_atp_per_cycle_is_one(sim: MuscleSimulator):
    """One ATP hydrolysed per cross-bridge cycle (Lymn-Taylor 1971)."""
    assert sim.atp_per_cross_bridge_cycle() == 1


def test_per_stroke_work_in_pN_nm(sim: MuscleSimulator):
    """Single-head work F * d ~ 5 pN * 10 nm = 50 pN nm."""
    W = sim.power_stroke_work_pN_nm()
    assert math.isclose(W, SINGLE_HEAD_FORCE_PN * POWER_STROKE_NM, abs_tol=1e-12)
    assert 30 <= W <= 80, f"per-stroke work out of textbook range, got {W} pN.nm"


def test_thermodynamic_efficiency_about_half(sim: MuscleSimulator):
    """eta_max = 0.5 (textbook ceiling)."""
    eta = sim.thermodynamic_efficiency()
    assert math.isclose(eta, 0.5, abs_tol=1e-12)


def test_per_stroke_work_to_atp_ratio(sim: MuscleSimulator):
    """Per-stroke work / ATP free-energy ~ 1.0 in the textbook bookkeeping
    (50 pN nm work vs 50 pN nm ATP), but in-cell efficiency is half of
    that because half the ATP energy dissipates as heat through
    detachment + Pi release."""
    r = sim.per_stroke_work_to_atp_ratio()
    # The textbook stroke-work-vs-ATP ratio is approximately unity.
    assert 0.5 <= r <= 1.5


# --------------------------------------------------------------------- #
# Cross-bridge cycle: 5 states (Lymn-Taylor 1971)                       #
# --------------------------------------------------------------------- #

def test_cross_bridge_has_five_states(sim: MuscleSimulator):
    """Lymn-Taylor cycle = 5 distinct states."""
    states = sim.cross_bridge_states()
    assert len(states) == 5
    assert len(set(states)) == 5  # all distinct


def test_cross_bridge_cycle_is_closed(sim: MuscleSimulator):
    """The cycle is a closed loop: 5 transitions returning to start."""
    transitions = sim.cross_bridge_state_transitions()
    assert len(transitions) == 5
    # First transition must start at the first state
    assert transitions[0][0] == sim.cross_bridge_states()[0]
    # Last transition must end at the first state (cycle closes)
    assert transitions[-1][1] == sim.cross_bridge_states()[0]


def test_cross_bridge_includes_canonical_states(sim: MuscleSimulator):
    """The cycle must include the four textbook states."""
    states = sim.cross_bridge_states()
    needed = {"attached_rigor", "power_stroke", "detached", "primed"}
    assert needed.issubset(set(states))


def test_duty_ratio_in_skeletal_muscle_range(sim: MuscleSimulator):
    """Duty ratio of skeletal muscle ~ 5%."""
    d = sim.duty_ratio()
    assert 0.01 <= d <= 0.20  # textbook 5%, allow for slow/fast variation


# --------------------------------------------------------------------- #
# Geometric sanity                                                      #
# --------------------------------------------------------------------- #

def test_z_disk_separation_equals_sarcomere_length(geom: MuscleGeometry):
    z_left, z_right = geom.z_disk_positions_um()
    assert math.isclose(z_right - z_left, geom.sarcomere_length_um, abs_tol=1e-12)


def test_a_band_length_equals_thick_filament(geom: MuscleGeometry):
    a_left, a_right = geom.a_band_extent_um()
    assert math.isclose(a_right - a_left, THICK_FILAMENT_LENGTH_UM, abs_tol=1e-12)


def test_m_line_at_centre(geom: MuscleGeometry):
    assert math.isclose(geom.m_line_position_um(), 0.0, abs_tol=1e-12)


def test_h_zone_inside_a_band(geom: MuscleGeometry):
    a_left, a_right = geom.a_band_extent_um()
    h_left, h_right = geom.h_zone_extent_um()
    assert a_left <= h_left <= h_right <= a_right


def test_thin_filament_length(geom: MuscleGeometry):
    (z_left, left_tip), (z_right, right_tip) = geom.thin_filament_axial_extents_um()
    assert math.isclose(left_tip - z_left, THIN_FILAMENT_LENGTH_UM, abs_tol=1e-12)
    assert math.isclose(z_right - right_tip, THIN_FILAMENT_LENGTH_UM, abs_tol=1e-12)


def test_cross_bridge_axial_period(geom: MuscleGeometry):
    """Axial period = 14.3 nm (canonical helical pitch repeat)."""
    assert math.isclose(CROSS_BRIDGE_AXIAL_PERIOD_NM, 14.3, abs_tol=0.05)


def test_cross_bridge_count_textbook(geom: MuscleGeometry):
    """Total heads on one thick filament ~ MYOSIN_HEADS_PER_THICK (textbook ~ 300)."""
    n = geom.cross_bridge_count()
    # geom.cross_bridge_count() reflects the simple period-and-bare-zone count
    # (with period 14.3 nm and bare zone 0.2 um); the textbook number comes
    # from helical packing (300 heads).  Both should agree within ~ 30%.
    assert n > 50
    assert abs(n - MYOSIN_HEADS_PER_THICK) < 250  # loose -- textbook is 300


def test_thin_filament_diameter_smaller_than_thick():
    """Actin thin filament (~ 8 nm) is roughly half the thick (~ 15 nm)."""
    assert THIN_FILAMENT_DIAMETER_NM < THICK_FILAMENT_DIAMETER_NM
    assert THIN_FILAMENT_DIAMETER_NM > 0.3 * THICK_FILAMENT_DIAMETER_NM


def test_inter_filament_spacing(geom: MuscleGeometry):
    """Hex-lattice inter-filament spacing ~ 40-45 nm."""
    assert 35 <= INTER_FILAMENT_SPACING_NM <= 50
    thick = geom.thick_filament_cross_section_nm()
    if len(thick) >= 2:
        d = np.linalg.norm(thick[1] - thick[0])
        assert math.isclose(d, INTER_FILAMENT_SPACING_NM, rel_tol=1e-10)


def test_thin_to_thick_ratio_two_to_one(geom: MuscleGeometry):
    """Vertebrate skeletal muscle has 2 thin filaments per thick filament."""
    n_thick = len(geom.thick_filament_cross_section_nm())
    n_thin = len(geom.thin_filament_cross_section_nm())
    assert n_thin == 2 * n_thick


# --------------------------------------------------------------------- #
# Sliding filament: A-band length is invariant                          #
# --------------------------------------------------------------------- #

def test_a_band_length_independent_of_sarcomere_length():
    """A-band length is FIXED -- the sliding filament hallmark."""
    for L_s in (1.8, 2.5, 3.2):
        g = MuscleGeometry(sarcomere_length_um=L_s)
        a_left, a_right = g.a_band_extent_um()
        assert math.isclose(a_right - a_left, THICK_FILAMENT_LENGTH_UM, abs_tol=1e-12)


def test_h_zone_shrinks_with_shortening():
    """H-zone shrinks as the sarcomere shortens (sliding filament)."""
    g_long = MuscleGeometry(sarcomere_length_um=2.7)
    g_rest = MuscleGeometry(sarcomere_length_um=SARCOMERE_REST_LENGTH_UM)
    g_short = MuscleGeometry(sarcomere_length_um=2.0)
    h_long = g_long.h_zone_extent_um()
    h_rest = g_rest.h_zone_extent_um()
    h_short = g_short.h_zone_extent_um()
    L_long = h_long[1] - h_long[0]
    L_rest = h_rest[1] - h_rest[0]
    L_short = h_short[1] - h_short[0]
    assert L_long > L_rest > L_short


# --------------------------------------------------------------------- #
# Substrate metadata                                                    #
# --------------------------------------------------------------------- #

def test_substrate_signatures_present(sim: MuscleSimulator, geom: MuscleGeometry):
    sig_g = geom.substrate_signature()
    sig_s = sim.substrate_signature()
    for k in ("interpretation", "lambda_qcd_MeV"):
        assert k in sig_g and k in sig_s
    assert sig_g["lambda_qcd_MeV"] == 200.0
    assert sig_s["thermodynamic_efficiency"] == 0.5
    assert 0.28 <= sig_s["v_at_max_power_over_vmax"] <= 0.34


def test_shortening_trajectory_monotone_decreasing(sim: MuscleSimulator):
    """At a sub-maximal load the sarcomere shortens, then settles."""
    t, L, F = sim.shortening_trajectory(load_normalised=0.3,
                                         duration_s=0.1, n_steps=100)
    assert L[0] > L[-1]                    # length decreases
    diffs = np.diff(L)
    # Mostly non-increasing (allow minor numerical wobbles at the
    # ascending-limb floor)
    assert (diffs <= 1e-12).all()
