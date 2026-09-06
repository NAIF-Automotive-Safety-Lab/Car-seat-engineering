from stiff_medium.dark_stress_scale_closure import (
    assess_dark_stress_scale_closure,
    dark_stress_speed_closure,
    hubble_length_kpc,
    polarization_coherence_length_closure,
)


def test_hubble_length_is_gpc_scale():
    assert 4.0e6 < hubble_length_kpc() < 5.0e6


def test_polarization_coherence_length_closure_is_kpc_scale():
    result = polarization_coherence_length_closure()

    assert 0.95 < result.ell_pol_kpc < 1.05


def test_dark_stress_speed_closure_is_cluster_scale():
    result = dark_stress_speed_closure()

    assert 900.0 < result.v_dark_km_s < 1100.0


def test_combined_scale_closure_hits_memory_target():
    result = assess_dark_stress_scale_closure()

    assert abs(result.tau_error_pct) < 1.0
    assert "strong" in result.verdict
