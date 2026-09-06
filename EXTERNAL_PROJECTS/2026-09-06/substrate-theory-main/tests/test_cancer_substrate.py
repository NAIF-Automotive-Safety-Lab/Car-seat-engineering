"""Tests for the paired cancer / tumor-growth substrate module.

Required claims (per the GEOMETRY/SIM/VIZ specification):
    (a) Gompertz growth saturates at the carrying capacity K and is
        approximately exponential at small N.
    (b) Armitage-Doll multistage incidence I(age) ∝ age^(n-1):
        with n_hits=6 the incidence on a log-log plot has slope 5.
    (c) Knudson 2-hit hypothesis: hereditary peak (~1-2 yr) is much
        earlier than sporadic peak (~30 yr) for retinoblastoma.
    (d) Skipper-Schabel log-kill: surviving fraction = (1-f)^n,
        kill_fraction_per_dose=0.99 gives 1e-10 after 5 doses.
    (e) 5 stages and 8 Hanahan-Weinberg hallmarks; hallmark count
        is monotone non-decreasing across the stage sequence.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.cancer_substrate import (
    ANGIOGENIC_SWITCH_RADIUS_MM,
    ARMITAGE_DOLL_HITS_DEFAULT,
    EXPONENTIAL_R_PER_DAY,
    GOMPERTZ_A_PER_DAY,
    GOMPERTZ_K_CELLS,
    HALLMARKS,
    KNUDSON_PEAK_AGE_HEREDITARY_YR,
    KNUDSON_PEAK_AGE_SPORADIC_YR,
    LAMBDA_QCD_MEV,
    LOG_KILL_FRACTION_PER_DOSE,
    N_HALLMARKS,
    N_STAGES,
    STAGE_DYSPLASTIC,
    STAGE_INVASIVE,
    STAGE_METASTATIC,
    STAGE_NEOPLASTIC,
    STAGE_NORMAL,
    STAGES,
    CancerGeometry,
    CancerSimulator,
)


# --------------------------------------------------------------------- #
# Fixtures                                                              #
# --------------------------------------------------------------------- #

@pytest.fixture
def geom() -> CancerGeometry:
    return CancerGeometry()


@pytest.fixture
def sim() -> CancerSimulator:
    return CancerSimulator()


# --------------------------------------------------------------------- #
# GEOMETRY                                                              #
# --------------------------------------------------------------------- #

def test_geometry_default_normal(geom):
    """A fresh CancerGeometry has all population in the normal stage."""
    assert geom.current_stage() == STAGE_NORMAL
    assert not geom.is_malignant()
    assert not geom.has_metastasized()


def test_geometry_stage_count():
    """Five canonical stages are defined."""
    assert N_STAGES == 5
    assert STAGES == (STAGE_NORMAL, STAGE_DYSPLASTIC, STAGE_NEOPLASTIC,
                      STAGE_INVASIVE, STAGE_METASTATIC)


def test_geometry_population_in_stage(geom):
    """population_in_stage returns the right entry from the array."""
    g = CancerGeometry(stage_populations=np.array(
        [1.0e9, 1.0e6, 1.0e4, 1.0e2, 1.0], dtype=float))
    assert g.population_in_stage(STAGE_NORMAL) == 1.0e9
    assert g.population_in_stage(STAGE_METASTATIC) == 1.0
    assert g.is_malignant()
    assert g.has_metastasized()
    assert g.current_stage() == STAGE_METASTATIC


def test_geometry_rejects_negative():
    """Negative populations are rejected at construction."""
    with pytest.raises(ValueError):
        CancerGeometry(stage_populations=np.array(
            [-1.0, 0.0, 0.0, 0.0, 0.0]))


def test_geometry_rejects_wrong_shape():
    """Wrong-shape stage_populations array is rejected."""
    with pytest.raises(ValueError):
        CancerGeometry(stage_populations=np.array([1.0, 2.0, 3.0]))


def test_geometry_volume_and_mass():
    """Sphere volume and ~1 mg/mm^3 mass conversion."""
    g = CancerGeometry(tumour_radius_mm=1.0)
    expected_vol = (4.0 / 3.0) * math.pi
    assert abs(g.tumour_volume_mm3() - expected_vol) < 1.0e-9
    # Mass ~ vol * 1e-3 g/mm^3
    assert abs(g.tumour_mass_g() - expected_vol * 1.0e-3) < 1.0e-12


def test_geometry_angiogenic_switch():
    """Angiogenic flag flips at the canonical 1.5 mm radius."""
    g_small = CancerGeometry(tumour_radius_mm=0.5)
    g_big   = CancerGeometry(tumour_radius_mm=2.0)
    assert not g_small.is_angiogenic()
    assert g_big.is_angiogenic()
    # At exactly the switch radius
    g_edge = CancerGeometry(tumour_radius_mm=ANGIOGENIC_SWITCH_RADIUS_MM)
    assert g_edge.is_angiogenic()


def test_geometry_cells_to_radius():
    """1e9 cells ~ 1 cm^3 ~ radius ~ 6.2 mm sphere."""
    r = CancerGeometry.cells_to_radius_mm(1.0e9)
    # 1e9 cells * 1 pL/cell = 1e9 pL = 1 mL = 1 cm^3 = 1000 mm^3
    # Sphere radius for V=1000 mm^3: r = (3*1000 / (4 pi))^(1/3)
    expected = (3.0 * 1000.0 / (4.0 * math.pi)) ** (1.0 / 3.0)
    assert abs(r - expected) / expected < 1.0e-6


def test_geometry_substrate_signature(geom):
    """substrate_signature returns the required keys."""
    sig = geom.substrate_signature()
    for key in ("interpretation", "n_stages", "n_hallmarks",
                "current_stage", "lambda_qcd_MeV"):
        assert key in sig
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
    assert sig["n_stages"] == N_STAGES
    assert sig["n_hallmarks"] == N_HALLMARKS


# --------------------------------------------------------------------- #
# (a) Gompertz growth saturates at K                                    #
# --------------------------------------------------------------------- #

def test_gompertz_saturates_at_K(sim):
    """Long-time limit of N(t) is the carrying capacity K."""
    t = np.array([0.0, 100.0, 1000.0, 10000.0])
    N = sim.gompertz_growth(t, N0=1.0)
    # Final N should be within 0.1% of K
    assert abs(N[-1] - sim.K_cells) / sim.K_cells < 1.0e-3


def test_gompertz_starts_at_N0(sim):
    """N(0) = N0 exactly."""
    N0 = 1.0
    N = sim.gompertz_growth(np.array([0.0]), N0=N0)
    assert abs(N[0] - N0) < 1.0e-9


def test_gompertz_monotone_increasing(sim):
    """For N0 < K, N(t) is strictly increasing."""
    t = np.linspace(0.0, 500.0, 200)
    N = sim.gompertz_growth(t, N0=1.0)
    diffs = np.diff(N)
    assert np.all(diffs > 0.0), "Gompertz must be monotone increasing for N0 < K"


def test_gompertz_exponential_at_small_N(sim):
    """At small N, Gompertz looks like exponential growth.

    dN/dt = N a ln(K/N).  At N << K, ln(K/N) is approximately a slowly
    changing constant ~ ln(K/N0), so dN/dt approx N * (a ln(K/N0)),
    giving an effective exponential rate r_eff = a*ln(K/N0).  We check
    that the early-time slope on a log-N plot matches r_eff to within
    a few percent over the first doubling.
    """
    N0 = 1.0
    # Small enough that ln(K/N) doesn't change much over one doubling
    a = sim.a_per_day
    K = sim.K_cells
    r_eff = a * math.log(K / N0)
    t = np.linspace(0.0, 0.05 / a, 20)  # very early
    N = sim.gompertz_growth(t, N0=N0)
    # Fit log N vs t
    slope = float(np.polyfit(t, np.log(N), 1)[0])
    assert abs(slope - r_eff) / r_eff < 0.05


def test_gompertz_rate_zero_at_K(sim):
    """At N = K, instantaneous growth rate is zero."""
    rate = sim.gompertz_rate(sim.K_cells)
    assert abs(rate) < 1.0e-9


def test_exponential_growth_doubles(sim):
    """N(t) doubles every t_d = ln(2)/r."""
    td = sim.doubling_time_days()
    assert abs(td - math.log(2.0) / sim.r_per_day) < 1.0e-12
    N0 = 100.0
    N = sim.exponential_growth(np.array([0.0, td, 2.0 * td]), N0=N0)
    assert abs(N[1] - 2.0 * N0) / N0 < 1.0e-9
    assert abs(N[2] - 4.0 * N0) / N0 < 1.0e-9


# --------------------------------------------------------------------- #
# (b) Armitage-Doll multistage incidence ~ age^(n-1)                    #
# --------------------------------------------------------------------- #

def test_armitage_doll_default_slope_5(sim):
    """log I vs log age has slope n-1 (= 5 for default n_hits=6)."""
    assert sim.armitage_doll_log_slope() == 5
    assert sim.n_armitage_hits == ARMITAGE_DOLL_HITS_DEFAULT == 6


def test_armitage_doll_incidence_age5(sim):
    """Incidence at 80 yr should be 32x incidence at 40 yr (2^5)."""
    I40 = sim.armitage_doll_incidence(np.array([40.0]))[0]
    I80 = sim.armitage_doll_incidence(np.array([80.0]))[0]
    ratio = I80 / I40
    assert abs(ratio - 32.0) < 0.5


def test_armitage_doll_loglog_slope_5(sim):
    """Linear regression of log I vs log age should yield slope ~5."""
    age = np.linspace(20.0, 90.0, 100)
    I = sim.armitage_doll_incidence(age)
    slope, _ = np.polyfit(np.log(age), np.log(I), 1)
    assert abs(slope - 5.0) < 1.0e-9


def test_armitage_doll_rejects_n_under_2(sim):
    """n_hits < 2 has no multistage interpretation."""
    with pytest.raises(ValueError):
        sim.armitage_doll_incidence(np.array([20.0, 40.0]), n_hits=1)


def test_armitage_doll_with_custom_n(sim):
    """Slope = n-1 for arbitrary n."""
    age = np.linspace(20.0, 90.0, 50)
    for n in (3, 5, 7):
        I = sim.armitage_doll_incidence(age, n_hits=n)
        slope, _ = np.polyfit(np.log(age), np.log(I), 1)
        assert abs(slope - (n - 1)) < 1.0e-9


# --------------------------------------------------------------------- #
# (c) Knudson 2-hit hypothesis                                          #
# --------------------------------------------------------------------- #

def test_knudson_hereditary_peak_earlier_than_sporadic(sim):
    """Hereditary peak (~1-2 yr) earlier than sporadic peak (~30 yr)."""
    sporadic, hereditary = sim.knudson_peak_ages_years()
    assert hereditary < sporadic
    assert hereditary == KNUDSON_PEAK_AGE_HEREDITARY_YR
    assert sporadic == KNUDSON_PEAK_AGE_SPORADIC_YR


def test_knudson_hereditary_pdf_peaks_early(sim):
    """Hereditary PDF peaks well below age 5 yr."""
    age = np.linspace(0.05, 10.0, 400)
    pdf = sim.knudson_two_hit_hereditary_pdf(age)
    peak_age = age[int(np.argmax(pdf))]
    assert peak_age < 5.0, f"hereditary peak at {peak_age:.2f} yr (expect <5)"


def test_knudson_sporadic_pdf_peaks_later(sim):
    """Sporadic PDF peaks above age 10 yr (much later than hereditary)."""
    age = np.linspace(0.5, 80.0, 400)
    pdf = sim.knudson_two_hit_sporadic_pdf(age)
    peak_age = age[int(np.argmax(pdf))]
    assert peak_age > 10.0, f"sporadic peak at {peak_age:.2f} yr (expect >10)"


# --------------------------------------------------------------------- #
# (d) Hanahan-Weinberg 8 hallmarks acquired across stages               #
# --------------------------------------------------------------------- #

def test_eight_hallmarks_total():
    """Hanahan-Weinberg defines 8 hallmarks (2011 update)."""
    assert N_HALLMARKS == 8
    assert len(HALLMARKS) == 8


def test_normal_stage_has_no_hallmarks(sim):
    """Normal cells have no acquired hallmarks."""
    h = sim.hallmarks_acquired(STAGE_NORMAL)
    assert len(h) == 0


def test_metastatic_stage_has_all_hallmarks(sim):
    """Metastatic cells have acquired all 8 hallmarks."""
    h = sim.hallmarks_acquired(STAGE_METASTATIC)
    assert len(h) == N_HALLMARKS
    assert tuple(h) == HALLMARKS


def test_hallmark_count_monotone(sim):
    """Number of acquired hallmarks is non-decreasing across stages."""
    counts = [sim.n_hallmarks_acquired(s) for s in STAGES]
    diffs = np.diff(counts)
    assert np.all(diffs >= 0)
    assert counts == [0, 2, 4, 6, 8]


def test_hallmark_unknown_stage(sim):
    """Asking for an unknown stage raises ValueError."""
    with pytest.raises(ValueError):
        sim.hallmarks_acquired("unknown_stage")


# --------------------------------------------------------------------- #
# (e) Skipper-Schabel log-kill survival                                 #
# --------------------------------------------------------------------- #

def test_log_kill_survival_one_dose(sim):
    """One dose at f=0.99 leaves 1% survival."""
    s = sim.log_kill_survival_fraction(1, kill_fraction_per_dose=0.99)
    assert abs(s - 0.01) < 1.0e-12


def test_log_kill_survival_five_doses(sim):
    """Five doses at f=0.99 leave 1e-10 survival."""
    s = sim.log_kill_survival_fraction(5, kill_fraction_per_dose=0.99)
    assert abs(s - 1.0e-10) < 1.0e-15


def test_log_kill_survival_constant_fraction(sim):
    """Survival decays as (1-f)^n geometrically."""
    f = LOG_KILL_FRACTION_PER_DOSE
    for n in (0, 1, 2, 3, 5, 10):
        s = sim.log_kill_survival_fraction(n)
        assert abs(s - (1.0 - f) ** n) < 1.0e-15


def test_log_kill_zero_doses(sim):
    """Zero doses gives full survival."""
    assert sim.log_kill_survival_fraction(0) == 1.0


def test_log_kill_treatment_trajectory(sim):
    """log_kill_treatment returns the geometric decay sequence."""
    N0 = 1.0e9
    n_doses = 4
    traj = sim.log_kill_treatment(N0, n_doses=n_doses,
                                  kill_fraction_per_dose=0.9)
    # Trajectory should have n_doses+1 entries: N0, 0.1 N0, 0.01 N0, ...
    assert len(traj) == n_doses + 1
    assert abs(traj[0] - N0) < 1.0e-9
    for i in range(n_doses + 1):
        assert abs(traj[i] - N0 * 0.1 ** i) / N0 < 1.0e-12


def test_log_kill_rejects_bad_fraction(sim):
    """Kill fraction outside [0, 1] is rejected."""
    with pytest.raises(ValueError):
        sim.log_kill_survival_fraction(3, kill_fraction_per_dose=1.5)
    with pytest.raises(ValueError):
        sim.log_kill_survival_fraction(3, kill_fraction_per_dose=-0.1)


def test_log_kill_rejects_negative_doses(sim):
    """Negative dose count is rejected."""
    with pytest.raises(ValueError):
        sim.log_kill_survival_fraction(-1)


# --------------------------------------------------------------------- #
# Combined growth + treatment                                           #
# --------------------------------------------------------------------- #

def test_combined_treatment_reduces_population(sim):
    """Combined growth + treatment leaves a visible drop just after each pulse.

    Compare the treated trajectory to the untreated trajectory immediately
    after the first pulse: the treated population must be strictly lower.
    """
    t = np.linspace(0.0, 250.0, 600)
    res = sim.combined_growth_and_treatment(
        t, N0=1.0, treatment_start_day=200.0,
        n_doses=6, kill_fraction_per_dose=0.99,
    )
    # Just after the first pulse the treated N must drop below untreated
    untreated = sim.gompertz_growth(t, N0=1.0)
    i_after_pulse = int(np.searchsorted(t, 200.5))
    assert res["N"][i_after_pulse] < untreated[i_after_pulse]
    # And the treated trajectory at the last pulse should be < untreated there
    last_pulse_t = res["dose_times"][-1]
    j = int(np.searchsorted(t, last_pulse_t + 0.5))
    if j < len(t):
        assert res["N"][j] < untreated[j]


def test_combined_treatment_dose_times(sim):
    """Dose times are spaced by dose_interval_days."""
    t = np.linspace(0.0, 200.0, 300)
    res = sim.combined_growth_and_treatment(
        t, N0=1.0, treatment_start_day=50.0,
        n_doses=4,
    )
    assert len(res["dose_times"]) == 4
    diffs = np.diff(res["dose_times"])
    for d in diffs:
        assert abs(d - sim.dose_interval_days) < 1.0e-9


# --------------------------------------------------------------------- #
# Substrate metadata                                                    #
# --------------------------------------------------------------------- #

def test_simulator_substrate_signature(sim):
    """Simulator substrate_signature includes lambda_qcd and key params."""
    sig = sim.substrate_signature()
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
    assert sig["n_hallmarks_total"] == N_HALLMARKS
    assert sig["gompertz_K_cells"] == GOMPERTZ_K_CELLS
    assert sig["gompertz_a_per_day"] == GOMPERTZ_A_PER_DAY
    assert sig["exponential_r_per_day"] == EXPONENTIAL_R_PER_DAY
    assert sig["log_kill_fraction_per_dose"] == LOG_KILL_FRACTION_PER_DOSE
    assert sig["armitage_log_slope"] == 5
    assert "interpretation" in sig
