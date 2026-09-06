from stiff_medium.substrate_polarization_dm import (
    assess_pure_polarization,
    cluster_separation_check,
    cosmological_a0,
    rotation_curve_check,
    solar_system_check,
)


def test_cosmological_a0_is_mond_scale():
    a0 = cosmological_a0()

    assert 0.9e-10 < a0 < 1.2e-10


def test_pure_polarization_gives_flat_galaxy_curve():
    result = rotation_curve_check()

    assert result.flatness_fraction < 0.15
    assert result.mass_ratio_at_outer_radius > 5.0
    assert result.verdict.startswith("passes")


def test_high_acceleration_shutoff_is_solar_system_safe():
    result = solar_system_check()

    assert result.fractional_excess < 1.0e-8
    assert result.verdict.startswith("passes")


def test_instantaneous_pure_polarization_fails_cluster_offsets():
    result = cluster_separation_check()

    assert result.required_memory_myr > 10.0
    assert not result.instantaneous_local_passes
    assert result.memory_turns_into_dark_stress


def test_strict_no_dm_fails_but_hybrid_survives():
    assessment = assess_pure_polarization()

    assert not assessment.strict_no_dm_passes
    assert assessment.hybrid_required
