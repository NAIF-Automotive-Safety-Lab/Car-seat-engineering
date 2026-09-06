import math

from stiff_medium.cone_lattice_microgeometry import (
    allowed_orientation_selectors,
    assess_cone_lattice_microgeometry,
    orientation_energy,
    parallel_perpendicular_squares,
    scan_orientation_minimum,
)


def test_parallel_perpendicular_split_matches_cone_geometry():
    parallel_sq, perpendicular_sq = parallel_perpendicular_squares([1.0, 0.0, 1.0])

    assert abs(parallel_sq - perpendicular_sq) < 1.0e-12


def test_axial_symmetry_alone_allows_quadratic_bias():
    terms = allowed_orientation_selectors(self_dual_exchange=False)

    assert any("m = p - q" in term for term in terms)


def test_self_dual_exchange_makes_quartic_first_selector():
    terms = allowed_orientation_selectors(self_dual_exchange=True)

    assert not any(term.startswith("m = p - q") for term in terms)
    assert any("m^2 = (p - q)^2" in term for term in terms)


def test_positive_quartic_selects_45_degrees():
    theta_min, energy_min = scan_orientation_minimum(linear_bias=0.0, beta=1.0)

    assert abs(math.degrees(theta_min) - 45.0) < 1.0e-3
    assert energy_min < 1.0e-12
    assert orientation_energy(0.0, beta=1.0) > energy_min


def test_allowed_linear_bias_spoils_exact_cone_angle():
    result = assess_cone_lattice_microgeometry()

    assert result.quadratic_bias_allowed_without_dual
    assert result.self_dual_exchange_required
    assert abs(result.biased_minimum_angle_deg - 45.0) > 1.0
    assert result.bias_shift_deg > 1.0
    assert not result.cone_forced_by_current_symmetry


def test_positive_beta_is_required_for_stable_45_degree_cone():
    result = assess_cone_lattice_microgeometry()

    assert result.beta_positive_required
    assert result.quartic_curvature_at_minimum > 0.0
    assert abs(result.negative_beta_minimum_angle_deg - 45.0) > 1.0
    assert "positive beta" in result.verdict
