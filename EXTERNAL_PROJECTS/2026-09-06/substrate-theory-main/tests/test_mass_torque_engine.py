"""Tests for src.stiff_medium.mass_torque_engine."""

import math

import numpy as np
import pytest

from src.stiff_medium.mass_torque_engine import (
    MassTorque,
    TorqueResult,
    DEFAULT_INTEGERS,
    DEFAULT_PRIMITIVES,
    LAMBDA_QCD_MEV,
)


@pytest.fixture
def engine():
    return MassTorque()


# ---------------------------------------------------------------------------
# Basic plumbing
# ---------------------------------------------------------------------------

def test_default_configuration(engine):
    assert engine.lambda_qcd == LAMBDA_QCD_MEV
    for k, v in DEFAULT_PRIMITIVES.items():
        assert engine.primitives[k] == v
    for k, v in DEFAULT_INTEGERS.items():
        assert engine.integers[k] == v


def test_list_configs_contains_all_named(engine):
    configs = engine.list_configs()
    for name in ("deuteron", "alpha", "muon", "tau", "higgs",
                 "hierarchy", "fine_structure", "t_c_max"):
        assert name in configs


def test_compute_returns_torque_result(engine):
    r = engine.compute("deuteron")
    assert isinstance(r, TorqueResult)
    assert r.name == "deuteron"
    assert r.value_mev > 0
    assert "Lambda_QCD" in r.formula or "n_A" in r.formula


def test_unknown_config_raises(engine):
    with pytest.raises(KeyError):
        engine.compute("does_not_exist")


def test_callable_shorthand(engine):
    a = engine("muon")
    b = engine.compute("muon")
    assert a.value_mev == b.value_mev


# ---------------------------------------------------------------------------
# All named configurations match observed within tolerance
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "deuteron", "alpha", "muon", "tau", "higgs",
    "hierarchy", "fine_structure", "t_c_max",
])
def test_verify_named_configuration(engine, name):
    v = engine.verify(name)
    assert v["passed"], (
        f"{name}: predicted={v['predicted']:.6g}, "
        f"observed={v['observed']:.6g}, rel_err={v['rel_err']:.3e}, "
        f"tol={v['tolerance']:.1e}"
    )


def test_deuteron_canonical_value(engine):
    v = engine.verify("deuteron")
    # 200/90 = 2.2222 MeV
    assert math.isclose(v["predicted"], 200.0 / 90.0, rel_tol=1e-9)


def test_muon_lepton_ratio_form(engine):
    """m_mu / m_e = exp(n_M / (K_pair^4 * pi))."""
    n_M = engine.integers["n_M"]
    K_pair = engine.integers["K_pair"]
    expected_log = n_M / (K_pair ** 4 * np.pi)
    from src.stiff_medium.mass_torque_engine import M_ELECTRON_MEV
    r = engine.compute("muon")
    assert math.isclose(r.value_mev, M_ELECTRON_MEV * np.exp(expected_log),
                        rel_tol=1e-12)


def test_hierarchy_exponential_form(engine):
    r = engine.compute("hierarchy")
    assert math.isclose(r.torque, np.exp(4 * np.pi ** 2 - 1), rel_tol=1e-12)


def test_fine_structure_formula(engine):
    r = engine.compute("fine_structure")
    expected = (11.0 / (48.0 * np.pi ** 3)) * np.exp(-3.0 * np.pi / 737.0)
    assert math.isclose(r.torque, expected, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Custom configurations
# ---------------------------------------------------------------------------

def test_predict_with_explicit_torque(engine):
    r = engine.predict({"name": "manual", "torque": 0.5})
    assert r.value_mev == 0.5 * engine.lambda_qcd
    assert r.name == "manual"


def test_predict_with_callable_formula(engine):
    f = lambda prims, ints: 1.0 / (ints["n_A"] * ints["N_BAM"])
    r = engine.predict({"name": "deuteron_custom", "formula": f})
    assert math.isclose(r.value_mev, 200.0 / 90.0, rel_tol=1e-9)


def test_predict_overrides_integers(engine):
    f = lambda prims, ints: 1.0 / ints["custom_n"]
    r = engine.predict({
        "name": "custom",
        "formula": f,
        "integers": {"custom_n": 4},
    })
    assert math.isclose(r.value_mev, 50.0, rel_tol=1e-12)


def test_predict_requires_torque_or_formula(engine):
    with pytest.raises(ValueError):
        engine.predict({"name": "broken"})


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------

def test_engine_is_deterministic():
    eng1 = MassTorque()
    eng2 = MassTorque()
    for name in eng1.list_configs():
        r1 = eng1.compute(name)
        r2 = eng2.compute(name)
        assert r1.value_mev == r2.value_mev
        assert r1.torque == r2.torque


def test_repeated_calls_match(engine):
    r1 = engine.compute("tau")
    r2 = engine.compute("tau")
    assert r1.value_mev == r2.value_mev


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def test_report_runs_and_returns_string(engine, capsys):
    text = engine.report()
    assert isinstance(text, str)
    captured = capsys.readouterr()
    assert "deuteron" in captured.out
    assert "muon" in captured.out


# ---------------------------------------------------------------------------
# Drag-as-mass-generator API
# ---------------------------------------------------------------------------

def test_cone_bouncing_freq_baseline(engine):
    # Canonical:  ω_b² = (c/ξ)² + (γ/2ρ)²,  c = √(K/ρ).
    # With unit primitives K=ρ=ξ=γ=1: ω_b = √(1 + 1/4) = √(5/4).
    expected = math.sqrt(1.0 + 0.25)
    assert math.isclose(engine.cone_bouncing_freq(), expected, rel_tol=1e-12)


def test_cone_bouncing_freq_scales_with_drag(engine):
    # Canonical adds γ in quadrature, not linearly.  At baseline γ=1 we
    # have ω² = 1 + 1/4; doubling γ to 2 gives ω² = 1 + 1.  Both must
    # exceed the bare Compton ω₀ = 1, and ω(γ=2) > ω(γ=1) (monotone).
    base = engine.cone_bouncing_freq()
    doubled = engine.cone_bouncing_freq({"gamma": 2.0})
    assert doubled > base
    assert math.isclose(doubled, math.sqrt(1.0 + 1.0), rel_tol=1e-12)
    assert math.isclose(base, math.sqrt(1.0 + 0.25), rel_tol=1e-12)


def test_cone_bouncing_freq_sqrt_K(engine):
    # K=4 ⇒ c=√4=2, ω₀=2, Γ=1/2 ⇒ ω_b = √(4 + 1/4) = √4.25.
    val = engine.cone_bouncing_freq({"K": 4.0})
    assert math.isclose(val, math.sqrt(4.0 + 0.25), rel_tol=1e-12)


def test_drag_torque_equals_omega_b(engine):
    assert engine.drag_torque() == engine.cone_bouncing_freq()


def test_electron_registered(engine):
    assert "electron" in engine.list_configs()


def test_electron_anchor_pinned(engine):
    from src.stiff_medium.mass_torque_engine import M_ELECTRON_MEV
    r = engine.compute("electron")
    assert math.isclose(r.value_mev, M_ELECTRON_MEV, rel_tol=1e-12)


def test_electron_verifies(engine):
    v = engine.verify("electron")
    assert v["passed"]
    assert v["rel_err"] < 1e-9


def test_verify_includes_drag_keys(engine):
    v = engine.verify("muon")
    assert "drag_torque" in v
    assert "drag_share" in v
    assert v["drag_torque"] > 0


def test_report_includes_drag_column(engine, capsys):
    engine.report()
    captured = capsys.readouterr()
    assert "drag" in captured.out
    assert "electron" in captured.out
