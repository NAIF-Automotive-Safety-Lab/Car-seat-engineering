"""Tests for CasimirEffect."""

from __future__ import annotations

import math

import pytest

from stiff_medium.casimir_effect import CasimirEffect, CasimirResult


@pytest.fixture
def casimir() -> CasimirEffect:
    return CasimirEffect()


def test_pressure_at_1_micron_matches_textbook(casimir: CasimirEffect) -> None:
    """At d = 1 um, P should be approximately -1.3 mPa (textbook value)."""
    p = casimir.pressure_parallel_plates(1.0e-6)
    # Negative -> attractive
    assert p < 0
    p_mPa = p * 1e3
    # Textbook value is about -1.30 mPa; allow 5% tolerance
    assert math.isclose(p_mPa, -1.30, rel_tol=0.05), f"got {p_mPa} mPa"


def test_force_per_area_alias(casimir: CasimirEffect) -> None:
    d = 1.0e-6
    assert casimir.force_per_area(d) == casimir.pressure_parallel_plates(d)


def test_d_minus_4_scaling(casimir: CasimirEffect) -> None:
    """Doubling d must reduce |P| by factor 2^4 = 16."""
    d1 = 1.0e-6
    d2 = 2.0e-6
    p1 = casimir.pressure_parallel_plates(d1)
    p2 = casimir.pressure_parallel_plates(d2)
    ratio = p1 / p2
    assert math.isclose(ratio, 16.0, rel_tol=1e-9)

    # Generalize to several scales
    p_half = casimir.pressure_parallel_plates(0.5e-6)
    assert math.isclose(p_half / p1, 16.0, rel_tol=1e-9)


def test_invalid_separation_raises(casimir: CasimirEffect) -> None:
    with pytest.raises(ValueError):
        casimir.pressure_parallel_plates(0.0)
    with pytest.raises(ValueError):
        casimir.pressure_parallel_plates(-1e-6)


def test_compute_returns_dataclass(casimir: CasimirEffect) -> None:
    r = casimir.compute(1e-6)
    assert isinstance(r, CasimirResult)
    assert r.distance_m == 1e-6
    assert r.pressure_Pa == r.force_per_area_Pa
    assert "pressure_Pa" in r.as_dict()


def test_casimir_polder_attractive_and_d5_scaling(casimir: CasimirEffect) -> None:
    alpha_H = 7.42e-41  # SI polarizability of hydrogen
    f1 = casimir.casimir_polder_atom_wall(1e-7, alpha_H)
    f2 = casimir.casimir_polder_atom_wall(2e-7, alpha_H)
    assert f1 < 0  # attractive
    # F ~ d^-5: doubling d -> factor 32 reduction
    assert math.isclose(f1 / f2, 32.0, rel_tol=1e-9)


def test_repulsive_casimir_munday_geometry(casimir: CasimirEffect) -> None:
    """Gold-bromobenzene-silica is the canonical repulsive case."""
    out = casimir.repulsive_casimir(("silica", "bromobenzene", "gold"))
    assert out["repulsive"] is True
    # Symmetric vacuum case is attractive (not repulsive)
    out2 = casimir.repulsive_casimir(("gold", "vacuum", "gold"))
    assert out2["repulsive"] is False


def test_repulsive_casimir_input_validation(casimir: CasimirEffect) -> None:
    with pytest.raises(ValueError):
        casimir.repulsive_casimir(("gold", "vacuum"))
    with pytest.raises(ValueError):
        casimir.repulsive_casimir(("gold", "unobtanium", "silica"))


def test_dynamical_casimir_estimate_positive(casimir: CasimirEffect) -> None:
    out = casimir.dynamical_casimir(
        omega_drive_hz=1e10, amplitude_m=1e-9, cavity_length_m=1e-2
    )
    assert out["photon_rate_Hz_estimate"] > 0
    assert math.isclose(out["emitted_photon_freq_Hz"], 5e9)


def test_compare_to_lamoreaux(casimir: CasimirEffect) -> None:
    res = casimir.compare_to_lamoreaux_experiment()
    # Reproduce -1.3 mPa at 1 um
    assert math.isclose(res["predicted_pressure_mPa"], -1.30, rel_tol=0.05)
    # Sphere-plate force should be attractive (negative) and microscopically small
    assert res["predicted_force_sphere_plate_N"] < 0
    # Within Lamoreaux's quoted 5% experimental band
    assert res["quoted_experimental_agreement_percent"] <= 5.0


def test_substrate_origin_keys(casimir: CasimirEffect) -> None:
    s = casimir.substrate_origin()
    for key in ("qft_picture", "substrate_picture", "scaling", "pressure_formula"):
        assert key in s
    assert "d^-4" in s["scaling"]


def test_report_structure(casimir: CasimirEffect) -> None:
    r = casimir.report()
    assert "scan" in r and len(r["scan"]) >= 3
    assert math.isclose(r["at_1_micron_mPa"], -1.30, rel_tol=0.05)
    assert r["lamoreaux_1997"]["quoted_experimental_agreement_percent"] == 5.0
    assert "substrate_origin" in r
