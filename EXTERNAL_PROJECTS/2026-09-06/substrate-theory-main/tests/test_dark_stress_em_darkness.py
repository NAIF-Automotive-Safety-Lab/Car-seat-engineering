from stiff_medium.dark_stress_em_darkness import assess_em_darkness


def test_heavy_neutral_component_is_em_dark_but_gravitational():
    result = assess_em_darkness()
    heavy = result.heavy_neutral_kink

    assert heavy.mass_gev > 10.0
    assert heavy.gate.em_dark
    assert heavy.gate.gravitationally_visible
    assert not heavy.gate.emits_photons
    assert not heavy.gate.absorbs_detector_photons
    assert not heavy.gate.reflects_em_fields


def test_heavy_component_has_non_detector_resonant_internal_scales():
    result = assess_em_darkness()
    heavy = result.heavy_neutral_kink

    assert heavy.compton_frequency_hz > 1.0e24
    assert heavy.shear_mode_frequency_hz > 1.0e19
    assert heavy.shear_mode_energy_ev > 1.0e4


def test_locked_polarization_is_ultra_low_frequency_coherent_stress():
    result = assess_em_darkness()
    coherent = result.coherent_polarization

    assert coherent.gate.em_dark
    assert coherent.gate.gravitationally_visible
    assert coherent.memory_frequency_hz < 1.0e-14
    assert coherent.memory_quantum_energy_ev < 1.0e-28
    assert coherent.light_wavelength_for_memory_frequency_kpc > 1.0e4


def test_hybrid_em_darkness_matches_heavy_plus_ultralight_constraint():
    result = assess_em_darkness()

    assert not result.ordinary_em_detection_expected
    assert result.gravitational_detection_expected
    assert "heavy neutral" in result.verdict
    assert "ultra-light coherent" in result.verdict
