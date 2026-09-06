"""Tests for the canonical cone-bouncing protocol module."""

from __future__ import annotations

import numpy as np
import pytest

from src.stiff_medium.cone_bouncing_protocol import (
    C, HBAR, M_E, XI_E,
    omega_b_canonical,
    mass_from_omega,
    dimensionless_drag_to_mass,
    omega_A, omega_B, omega_D,
    cross_compare,
    electron_anchor_primitives,
    predicted_electron_mass_mev,
)


# ---------------------------------------------------------------------------
# (a) Electron mass at the anchor
# ---------------------------------------------------------------------------

def test_electron_anchor_recovers_m_e():
    """At ξ = ξ_e and γ = 0 the canonical formula must give m_e exactly."""
    p = electron_anchor_primitives()
    omega = omega_b_canonical(**p)

    # ω_b should equal m_e c² / ℏ
    omega_target = M_E * C * C / HBAR
    assert omega == pytest.approx(omega_target, rel=1e-12)

    # And the inferred mass = m_e
    m_kg = mass_from_omega(omega)
    assert m_kg == pytest.approx(M_E, rel=1e-12)

    # In MeV: ≈ 0.510999
    assert predicted_electron_mass_mev() == pytest.approx(0.510999, rel=1e-4)


def test_drag_only_adds_mass():
    """γ > 0 must monotonically *increase* ω_b above the Compton baseline."""
    p = electron_anchor_primitives()
    base = omega_b_canonical(**p)
    p_drag = dict(p)
    p_drag["gamma"] = 1e-12  # arbitrary positive drag
    dressed = omega_b_canonical(**p_drag)
    assert dressed >= base


# ---------------------------------------------------------------------------
# (b) Dimensionless ratio O(1)
# ---------------------------------------------------------------------------

def test_dimensionless_ratio_O1_in_drag_dominated_limit():
    """In the γ-dominated regime, m c² / (ℏγ) → 1/2  (Γ = γ/2ρ).

    Choose primitives where the drag term dominates the Compton term.
    """
    rho = 1.0
    K   = rho * C * C
    # Make ξ huge so ω₀ = c/ξ is tiny, and γ moderate so Γ dominates.
    xi    = 1.0      # 1 m → ω₀ ≈ c
    gamma = 1e20     # huge drag → Γ ≫ ω₀

    omega = omega_b_canonical(K, rho, xi, gamma)
    m_kg  = mass_from_omega(omega)
    ratio = dimensionless_drag_to_mass(m_kg, gamma)

    # Should be close to 1/2 (since ω_b ≈ Γ = γ/2ρ ⇒ ℏω_b/c² · c²/(ℏγ) = 1/(2ρ))
    # With ρ = 1 (natural-baseline) the ratio is exactly 1/2.
    assert 0.1 < ratio < 10.0, f"ratio={ratio} not O(1)"


def test_canonical_reduces_to_compton_at_zero_drag():
    """ω_b(γ=0) = c/ξ exactly."""
    p = electron_anchor_primitives()
    omega = omega_b_canonical(**p)
    expected = np.sqrt(p["K"] / p["rho"]) / p["xi"]
    assert omega == pytest.approx(expected, rel=1e-15)


# ---------------------------------------------------------------------------
# (c) Cross-comparison report — show the 3 modules' formulas + canonical
# ---------------------------------------------------------------------------

def test_cross_comparison_report(capsys):
    """Print the comparison and assert the audit-found inconsistencies."""
    p = electron_anchor_primitives()
    # Use a non-trivial drag so the formulas can disagree.
    K, rho, xi = p["K"], p["rho"], p["xi"]
    gamma = 1e-12

    rows = cross_compare(K, rho, xi, gamma)

    print("\n=== ω_b cross-comparison @ electron anchor (γ = 1e-12) ===")
    for key, r in rows.items():
        print(f"  [{key:9s}] {r.label:42s}  ω_b={r.omega_b:.4e}  "
              f"Δ={r.rel_err_vs_canonical:.3e}")

    # Canonical and D must agree (D *is* the canonical form).
    assert rows["D"].rel_err_vs_canonical == pytest.approx(0.0, abs=1e-15)

    # A agrees with canonical at γ=0; at finite γ they may differ — but
    # both retain the c/ξ baseline so they agree to within O(γ̃) here.
    # We just require A is the same order of magnitude as canonical.
    assert 0.1 < rows["A"].omega_b / rows["canonical"].omega_b < 10.0

    # B *cannot* recover the Compton baseline; at small γ it underestimates
    # by many orders of magnitude.  This is the audit-flagged failure.
    assert rows["B"].omega_b < 1e-3 * rows["canonical"].omega_b

    captured = capsys.readouterr()
    assert "canonical" in captured.out
    assert "drag_mass_generator" in captured.out
    assert "mass_torque_engine" in captured.out
    assert "lattice_substrate_2d" in captured.out


def test_input_validation():
    with pytest.raises(ValueError):
        omega_b_canonical(K=0, rho=1, xi=1, gamma=0)
    with pytest.raises(ValueError):
        omega_b_canonical(K=1, rho=1, xi=1, gamma=-1.0)
