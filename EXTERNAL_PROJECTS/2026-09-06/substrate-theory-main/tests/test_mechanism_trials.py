from stiff_medium.mechanism_trials import (
    ckm_overlap_trials,
    cosmology_biharmonic_opacity_trial,
    dark_matter_polarization_halo_trial,
    lepton_boundary_loop_trials,
    orientation_vortex_bias_trials,
    uv_phase_slip_trials,
)


def test_uv_phase_slip_has_action_level_hit():
    trials = uv_phase_slip_trials()

    assert min(abs(t.action_error_pct) for t in trials) < 0.1


def test_lepton_boundary_loop_correction_hits_subpercent_ratios():
    best = lepton_boundary_loop_trials()[1]

    assert abs(best.mu_error_pct) < 0.3
    assert abs(best.tau_error_pct) < 0.3


def test_ckm_two_axis_overlap_is_close_to_cabibbo():
    best = ckm_overlap_trials()[0]

    assert abs(best.sin_error_pct) < 0.3


def test_biharmonic_opacity_passes_visibility_requirement():
    trial = cosmology_biharmonic_opacity_trial()

    assert trial.f_galaxy < trial.required_f_vis
    assert trial.f_acoustic > 0.5


def test_dark_matter_polarization_halo_is_self_interacting_scale():
    trial = dark_matter_polarization_halo_trial()

    assert 0.1 <= trial.sigma_over_m <= 1.0


def test_orientation_vortex_determinant_reaches_antimatter_suppression():
    trials = orientation_vortex_bias_trials()

    assert trials[1].anti_fraction < 2.0e-18
