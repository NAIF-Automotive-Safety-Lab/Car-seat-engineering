from stiff_medium.dark_stress_hybrid import (
    assess_hybrid_dark_stress,
    hybrid_cluster_check,
    hybrid_galaxy_check,
    hybrid_self_interaction,
    minimum_mobile_fraction_for_lensing,
    split_dark_stress,
)


def test_split_dark_stress_conserves_effective_ratio():
    split = split_dark_stress(5.36, 0.85)

    assert split.mobile_to_baryon + split.polarization_to_baryon == split.effective_dark_to_baryon
    assert split.total_to_baryon == 6.36


def test_minimum_mobile_fraction_is_high_for_cluster_lensing():
    minimum = minimum_mobile_fraction_for_lensing(5.36, 0.70)

    assert 0.82 < minimum < 0.84


def test_hybrid_cluster_check_passes_with_mostly_mobile_dark_stress():
    result = hybrid_cluster_check(mobile_fraction_of_dark=0.85)

    assert result.passes_mobile_lensing
    assert result.passes_memory_offset
    assert result.mobile_fraction_of_total_lensing > 0.70


def test_hybrid_cluster_check_fails_if_dark_stress_is_mostly_locked():
    result = hybrid_cluster_check(mobile_fraction_of_dark=0.20)

    assert not result.passes_mobile_lensing


def test_hybrid_galaxy_decomposition_retains_rotation_success():
    result = hybrid_galaxy_check()

    assert result.verdict.startswith("passes")
    assert result.outer_mobile_to_baryon > result.outer_polarization_to_baryon


def test_hybrid_self_interaction_lands_in_interesting_range():
    result = hybrid_self_interaction()

    assert 0.1 <= result.sigma_over_m_cm2_g <= 1.0


def test_hybrid_assessment_rejects_wimp_but_keeps_dark_stress():
    assessment = assess_hybrid_dark_stress()

    assert not assessment.fundamental_wimp_needed
    assert assessment.substrate_dark_stress_needed
    assert assessment.verdict.startswith("hybrid route passes")
