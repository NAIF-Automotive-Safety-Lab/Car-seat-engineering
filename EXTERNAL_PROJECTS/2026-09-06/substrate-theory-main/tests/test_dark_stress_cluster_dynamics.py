from stiff_medium.dark_stress_cluster_dynamics import assess_cluster_dynamics


def test_cluster_dynamics_keeps_mobile_component_dominant():
    result = assess_cluster_dynamics()

    assert result.mobile_fraction_of_total_lensing > 0.70
    assert result.mobile_peak_to_gas_locked_ratio > 2.0
    assert result.passes_mobile_lensing


def test_derived_memory_clock_hits_cluster_offset():
    result = assess_cluster_dynamics()

    assert abs(result.offset_error_pct) < 1.0
    assert 145.0 < result.predicted_memory_offset_kpc < 155.0
    assert result.passes_offset_memory


def test_polarization_is_local_memory_not_cluster_transport():
    result = assess_cluster_dynamics()

    assert 45.0 < result.stress_horizon_kpc < 55.0
    assert result.stress_horizon_fraction_of_offset < 0.40
    assert not result.polarization_alone_spans_offset


def test_cluster_dynamics_verdict_requires_transport_equations():
    result = assess_cluster_dynamics()

    assert "mobile kink" in result.verdict
    assert "transport equations" in result.verdict
