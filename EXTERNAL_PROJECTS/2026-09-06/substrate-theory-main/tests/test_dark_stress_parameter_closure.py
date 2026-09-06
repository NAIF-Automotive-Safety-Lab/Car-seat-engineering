from stiff_medium.dark_stress_parameter_closure import (
    assess_dark_stress_parameter_closure,
    dark_to_baryon_phase_space_closure,
    halo_radius_alpha_closure,
    mobile_fraction_phase_closure,
    polarization_memory_time_candidate,
)


def test_dark_to_baryon_phase_space_closure_is_subpercent():
    result = dark_to_baryon_phase_space_closure()

    assert abs(result.error_pct) < 0.3


def test_mobile_fraction_phase_closure_passes_cluster_threshold():
    result = mobile_fraction_phase_closure()

    assert result.mobile_fraction > result.minimum_required_mobile_fraction
    assert result.locked_fraction < 0.17


def test_memory_candidate_matches_required_scale_but_is_not_closed():
    result = polarization_memory_time_candidate()

    assert abs(result.error_pct) < 0.1
    assert "not derived" in result.verdict


def test_halo_radius_alpha_closure_stays_self_interacting():
    result = halo_radius_alpha_closure()

    assert 0.1 <= result.sigma_over_m_cm2_g <= 1.0


def test_dark_stress_parameter_closure_verdict_keeps_open_clock():
    assessment = assess_dark_stress_parameter_closure()

    assert "promising" in assessment.verdict
    assert "memory clock" in assessment.verdict
