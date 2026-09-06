import math

from stiff_medium.missing_piece_hypotheses import (
    KOIDE_TARGET,
    dark_matter_geometric_cross_sections,
    evaluate_uv_candidate,
    find_empirical_foot_phase,
    foot_phase_point,
    orientation_bias_requirements,
    radius_for_self_interaction,
    required_visibility,
    scan_transfer_windows,
)


def test_uv_exp_16pi_prefactor_is_near_target_but_not_exact():
    candidate = evaluate_uv_candidate(n_inst=16, p_sigma_lat=2)

    assert abs(candidate.relative_error) < 0.1
    assert "needs a real UV action mechanism" in candidate.verdict


def test_foot_pi_over_six_is_signed_branch_not_positive_koide_branch():
    point = foot_phase_point("pi/6", 1.0 / 6.0)

    assert math.isclose(point.signed_koide_q, KOIDE_TARGET, rel_tol=1e-12)
    assert abs(point.positive_koide_q - KOIDE_TARGET) > 0.1
    assert -1 in point.root_signs


def test_empirical_foot_phase_has_expected_conjugate_branch():
    empirical = find_empirical_foot_phase()
    conjugate = (2.0 - empirical) % 2.0

    assert math.isclose(empirical, 0.5959313664, rel_tol=0.0, abs_tol=1e-6)
    assert math.isclose(conjugate, 1.4040686336, rel_tol=0.0, abs_tol=1e-6)


def test_transfer_windows_can_pass_math_but_need_physics():
    required = required_visibility()
    windows = scan_transfer_windows()

    assert required < 5.0e-4
    assert any(w.passes_visibility and w.keeps_acoustic_visible for w in windows)


def test_dark_matter_qcd_size_dimer_is_effectively_collisionless():
    rows = dark_matter_geometric_cross_sections(radii_fm=(1.0,))

    assert rows[0].sigma_over_m_cm2_per_g < 1.0e-2
    assert radius_for_self_interaction(1.0) > 10.0


def test_orientation_bias_requirement_is_large_for_tiny_antifraction():
    row = orientation_bias_requirements((1.0e-18,))[0]

    assert row.delta_e_over_t_eff > 40.0
    assert row.tau_over_epoch_max < 0.03
