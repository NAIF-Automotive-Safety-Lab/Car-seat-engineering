from stiff_medium.dark_stress_cluster_dynamics import assess_cluster_dynamics
from stiff_medium.dark_stress_transport import (
    assess_transport_profiles,
    post_collision_profiles,
    profile_mass,
)


def test_transport_profiles_preserve_dark_stress_masses():
    profiles = post_collision_profiles()
    dynamics = assess_cluster_dynamics()

    assert abs(profile_mass(profiles.baryon, profiles.x_kpc) - 1.0) < 1.0e-6
    assert (
        abs(profile_mass(profiles.mobile_kink, profiles.x_kpc) - dynamics.mobile_to_baryon)
        < 1.0e-6
    )
    assert (
        abs(
            profile_mass(profiles.locked_polarization, profiles.x_kpc)
            - dynamics.locked_to_baryon
        )
        < 1.0e-6
    )


def test_finite_speed_transport_keeps_polarization_local():
    result = assess_transport_profiles()

    assert result.polarization_peak_kpc == 0.0
    assert result.polarization_density_at_mobile_peak == 0.0
    assert result.polarization_leakage_fraction < 1.0e-6


def test_total_lensing_peak_tracks_mobile_component():
    result = assess_transport_profiles()

    assert abs(result.total_peak_offset_error_kpc) < 2.0
    assert abs(result.total_peak_kpc - result.mobile_peak_kpc) < 1.0
    assert result.mobile_to_central_peak_ratio > 2.0


def test_transport_verdict_keeps_hybrid_boundary():
    result = assess_transport_profiles()

    assert "mobile neutral-stress peak" in result.verdict
    assert "polarization remains local" in result.verdict
