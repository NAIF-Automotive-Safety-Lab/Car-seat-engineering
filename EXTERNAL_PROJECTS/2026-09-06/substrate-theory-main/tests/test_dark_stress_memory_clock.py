from stiff_medium.dark_stress_memory_clock import (
    assess_memory_clock,
    coherence_crossing_clock,
    free_fall_clock,
    relaxation_count,
    self_interaction_clock,
)


def test_relaxation_count_is_expected_scale():
    assert 48.8 < relaxation_count() < 49.0


def test_coherence_crossing_clock_is_close_to_required_memory():
    result = coherence_crossing_clock()

    assert abs(result.error_pct) < 5.0
    assert result.tau_clock_myr < 1.1


def test_self_interaction_clock_is_too_slow_at_cluster_density():
    result = self_interaction_clock(density_kg_m3=1.0e-22)

    assert result.tau_clock_myr > 1000.0
    assert result.verdict.startswith("too slow")


def test_free_fall_clock_required_density_is_cluster_core_scale():
    result = free_fall_clock()

    assert 1.0e7 <= result.required_density_msun_kpc3 <= 1.0e8
    assert "possible" in result.verdict


def test_memory_clock_assessment_selects_coherence_crossing():
    assessment = assess_memory_clock()

    assert "coherence crossing" in assessment.verdict
    assert "too slow" in assessment.self_interaction_clock.verdict
