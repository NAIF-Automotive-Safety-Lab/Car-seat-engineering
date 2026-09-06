import math

import pytest

from stiff_medium.alpha_audit import build_alpha_audit_summary
from stiff_medium.alpha_derivation import run_alpha_derivation_audit
from stiff_medium.physical_constants import (
    ALPHA_THOMSON_CODATA_2022,
    INV_ALPHA_THOMSON_CODATA_2022,
)


def test_alpha_derivation_audit_reports_no_self_consistent_derivation():
    audit = run_alpha_derivation_audit()

    assert audit.target_inv_alpha == INV_ALPHA_THOMSON_CODATA_2022
    assert audit.target_alpha == ALPHA_THOMSON_CODATA_2022
    assert audit.n_fully_consistent_scan_points == 0
    assert audit.higgs_w_has_positive_alpha is False
    assert audit.higgs_w_g_thirring < 0.0
    assert audit.beta_gap_pi == pytest.approx(0.637481590, abs=5.0e-10)
    assert "No first-principles alpha derivation" in audit.conclusion


def test_alpha_two_loop_audit_keeps_fixed_result_separate_from_fit():
    summary = build_alpha_audit_summary()
    audit = summary.two_loop

    assert audit.target_inv_alpha == INV_ALPHA_THOMSON_CODATA_2022
    assert audit.higgs_w_has_positive_alpha is False

    assert audit.fixed_beta.beta_squared / math.pi == pytest.approx(3.91, abs=1.0e-12)
    assert audit.fixed_beta.inv_alpha_2loop == pytest.approx(136.507203650, abs=5.0e-10)
    assert summary.fixed_delta_inv_alpha == pytest.approx(-0.528795527, abs=5.0e-10)
    assert summary.fixed_abs_pct_inv_alpha == pytest.approx(0.385880739, abs=5.0e-10)
    assert summary.fixed_delta_alpha == pytest.approx(2.826815940716512e-05, abs=5.0e-17)

    assert audit.fitted_beta is not None
    assert audit.fitted_beta.inv_alpha_2loop == pytest.approx(audit.target_inv_alpha, abs=1.0e-9)
    assert audit.fitted_beta.beta_squared / math.pi != pytest.approx(3.91, abs=1.0e-6)
    assert "exact agreement requires fitting beta" in audit.conclusion
    assert "No first-principles alpha derivation" in summary.conclusion
