from stiff_medium.dark_stress_factor_scan import assess_factor_scan


def test_combined_memory_constraint_has_factor_degeneracy():
    result = assess_factor_scan()

    assert result.subpercent_tau_candidates > 1
    assert result.best_tau_candidate.key != (3, 3, 1, 5)


def test_physical_windows_select_scale_closure_candidate():
    result = assess_factor_scan()

    assert result.physical_subpercent_tau_candidates == 1
    assert result.best_physical_candidate.key == (3, 3, 1, 5)
    assert abs(result.best_physical_candidate.tau_error_pct) < 1.0


def test_factor_scan_keeps_derivation_gap_visible():
    result = assess_factor_scan()

    assert "degenerate" in result.verdict
    assert "derive" in result.verdict
