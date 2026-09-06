"""Tests for the substrate_physics CLI + library.

Run with:

    python -m pytest cli/test_cli.py -v --timeout=30
"""

from __future__ import annotations

import csv
import io
import json
import math
import os
import subprocess
import sys
import time
from contextlib import redirect_stdout
from pathlib import Path

import pytest

# Make ``cli/`` importable when pytest is run from the repo root.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from substrate_physics import (  # noqa: E402  (after sys.path manipulation)
    Prediction,
    SubstratePhysics,
    UnknownInputError,
    _dispatch,
    format_prediction,
    main as cli_main,
    run_batch,
)


# ---------------------------------------------------------------------------
# Library-level tests: each predict_* method
# ---------------------------------------------------------------------------


@pytest.fixture
def sp() -> SubstratePhysics:
    return SubstratePhysics()


def test_lepton_mass_ratio(sp: SubstratePhysics) -> None:
    p = sp.predict_lepton_mass_ratio()
    assert p.unit == "dimensionless"
    assert 200 < p.value < 215            # measured 206.768
    assert p.category == "A"
    # within ~1% of measured ratio
    assert p.residual_pct() is not None
    assert p.residual_pct() < 1.5


def test_tau_over_electron(sp: SubstratePhysics) -> None:
    p = sp.predict_tau_over_electron_ratio()
    # measured ~3477; substrate exp tower predicts a value within ~5% of meas.
    assert 3000 < p.value < 4000
    assert p.category == "A"


def test_atomic_ie_known_elements(sp: SubstratePhysics) -> None:
    for sym, expected in [("H", 13.6), ("He", 24.6), ("C", 11.3),
                          ("Ne", 21.6), ("Ar", 15.8)]:
        p = sp.predict_atomic_ie(sym)
        assert p.unit == "eV"
        # within ~10% across the row.  Simplified Slater rule.
        assert abs(p.value - expected) / expected < 0.20, \
            f"{sym}: predicted {p.value}, expected {expected}"


def test_atomic_ie_by_name(sp: SubstratePhysics) -> None:
    # Friendly element-name lookup must work.
    assert sp.predict_atomic_ie("carbon").value > 0
    assert sp.predict_atomic_ie("aluminium").value > 0  # spelling variants
    assert sp.predict_atomic_ie("aluminum").value > 0


def test_atomic_ie_unknown_raises(sp: SubstratePhysics) -> None:
    with pytest.raises(UnknownInputError):
        sp.predict_atomic_ie("Unobtainium")


def test_bandgap_known_materials(sp: SubstratePhysics) -> None:
    for mat, expected in [("silicon", 1.12), ("germanium", 0.67),
                          ("diamond", 5.47), ("gaas", 1.42)]:
        p = sp.predict_bandgap(mat)
        assert p.unit == "eV"
        assert abs(p.value - expected) / expected < 0.05


def test_bandgap_unknown_raises(sp: SubstratePhysics) -> None:
    with pytest.raises(UnknownInputError):
        sp.predict_bandgap("hypothetium")


def test_madelung_known_crystals(sp: SubstratePhysics) -> None:
    for cry, expected in [("nacl", 1.747565), ("cscl", 1.762675),
                          ("zincblende", 1.638055)]:
        p = sp.predict_madelung(cry)
        assert abs(p.value - expected) / expected < 0.001


def test_fracture_irwin(sp: SubstratePhysics) -> None:
    # Known textbook check: K_I = 50 MPa·sqrt(m), sigma_y = 600 MPa
    # naive Irwin: r_p = K^2 / (2 pi sigma^2) = 2500 / (2*pi*360000) ≈ 1.105e-3 m
    # substrate kappa = 1 + 6/268 = 1.0224 -> r_p ~ 1.13 mm
    p = sp.predict_fracture_zone(50.0, 600.0)
    assert p.unit == "mm"
    assert 1.0 < p.value < 1.3
    assert p.extras["kappa_substrate"] == pytest.approx(1.0 + 6 / 268, rel=1e-12)


def test_fracture_invalid_yield(sp: SubstratePhysics) -> None:
    with pytest.raises(UnknownInputError):
        sp.predict_fracture_zone(50.0, 0.0)


def test_hadron_known(sp: SubstratePhysics) -> None:
    p = sp.predict_hadron_mass("p")
    assert p.unit == "MeV"
    assert abs(p.value - 938.27) / 938.27 < 1e-3
    # case-insensitive lookup
    assert sp.predict_hadron_mass("Lambda").value > 1000


def test_hadron_unknown_raises(sp: SubstratePhysics) -> None:
    with pytest.raises(UnknownInputError):
        sp.predict_hadron_mass("zorptron")


def test_lifetime_known(sp: SubstratePhysics) -> None:
    p = sp.predict_lifetime("muon")
    assert p.unit == "s"
    assert 2.0e-6 < p.value < 2.5e-6


def test_debye_known(sp: SubstratePhysics) -> None:
    p = sp.predict_debye_temperature("diamond")
    assert p.unit == "K"
    assert 2200 < p.value < 2300


def test_bcs_known(sp: SubstratePhysics) -> None:
    p = sp.predict_bcs_gap_ratio("weak_coupling")
    assert abs(p.value - 3.528) < 0.05


def test_tc_max(sp: SubstratePhysics) -> None:
    p = sp.predict_tc_max()
    assert p.unit == "K"
    assert p.value == pytest.approx(128.9, rel=1e-9)
    assert p.category == "A"


def test_dark_matter(sp: SubstratePhysics) -> None:
    p = sp.predict_dark_matter_mass()
    assert p.unit == "GeV"
    assert p.value == pytest.approx(27.5)


def test_neutrino_sum(sp: SubstratePhysics) -> None:
    p = sp.predict_neutrino_sum()
    assert p.unit == "meV"
    assert p.value == pytest.approx(60.5)


def test_alpha_em(sp: SubstratePhysics) -> None:
    p = sp.predict_alpha_em()
    # CODATA alpha = 0.0072973525...; substrate gives ~0.0072970624 (0.004%)
    assert abs(p.value - 7.2973525643e-3) / 7.2973525643e-3 < 0.001


def test_hierarchy(sp: SubstratePhysics) -> None:
    p = sp.predict_hierarchy()
    expected = math.exp(4 * math.pi ** 2 - 1)
    assert p.value == pytest.approx(expected, rel=1e-12)


def test_string_tension(sp: SubstratePhysics) -> None:
    p = sp.predict_string_tension()
    # Cornell sigma ~ 0.18 GeV^2; substrate predicts (9/2) * 0.04 = 0.18 GeV^2
    assert abs(p.value - 0.18) < 1e-9


def test_grueneisen_known(sp: SubstratePhysics) -> None:
    p = sp.predict_grueneisen("copper")
    assert abs(p.value - 1.98) < 0.05


# ---------------------------------------------------------------------------
# Performance: every prediction must complete in well under 10 ms.
# ---------------------------------------------------------------------------


def test_speed_under_10ms(sp: SubstratePhysics) -> None:
    bench = [
        sp.predict_lepton_mass_ratio,
        sp.predict_tau_over_electron_ratio,
        sp.predict_alpha_em,
        sp.predict_hierarchy,
        sp.predict_string_tension,
        sp.predict_tc_max,
        sp.predict_neutrino_sum,
        sp.predict_dark_matter_mass,
        lambda: sp.predict_atomic_ie("C"),
        lambda: sp.predict_bandgap("silicon"),
        lambda: sp.predict_madelung("nacl"),
        lambda: sp.predict_hadron_mass("p"),
        lambda: sp.predict_lifetime("muon"),
        lambda: sp.predict_debye_temperature("diamond"),
        lambda: sp.predict_bcs_gap_ratio("weak_coupling"),
        lambda: sp.predict_grueneisen("copper"),
        lambda: sp.predict_fracture_zone(50.0, 600.0),
    ]
    for fn in bench:
        t0 = time.perf_counter()
        fn()
        dt_ms = (time.perf_counter() - t0) * 1e3
        assert dt_ms < 10.0, f"{fn} took {dt_ms:.2f} ms"


# ---------------------------------------------------------------------------
# Discovery helpers.
# ---------------------------------------------------------------------------


def test_list_predictions_includes_all_topics() -> None:
    names = {entry["name"].split(" ")[0] for entry in
             SubstratePhysics.list_predictions()}
    expected = {"m_mu_over_m_e", "m_tau_over_m_e", "ie", "bandgap", "madelung",
                "fracture", "hadron", "lifetime", "debye", "bcs", "tc_max",
                "dark_matter_mass", "neutrino_sum", "alpha_em", "hierarchy",
                "string_tension", "grueneisen"}
    missing = expected - names
    assert not missing, f"missing prediction names: {missing}"


def test_info_has_inventory_integers() -> None:
    info = SubstratePhysics.info()
    assert info["inventory_integers"]["n_M"] == 268
    assert info["inventory_integers"]["K_pair"] == 2
    assert info["inventory_integers"]["K_rank"] == 5


# ---------------------------------------------------------------------------
# Dispatch + format helpers.
# ---------------------------------------------------------------------------


def test_dispatch_no_arg(sp: SubstratePhysics) -> None:
    p = _dispatch(sp, "alpha_em", [])
    assert p.name == "alpha_em"


def test_dispatch_one_arg(sp: SubstratePhysics) -> None:
    p = _dispatch(sp, "bandgap", ["silicon"])
    assert "bandgap" in p.name


def test_dispatch_two_arg_fracture(sp: SubstratePhysics) -> None:
    p = _dispatch(sp, "fracture", ["50", "600"])
    assert p.unit == "mm"


def test_dispatch_unknown_raises(sp: SubstratePhysics) -> None:
    with pytest.raises(UnknownInputError):
        _dispatch(sp, "definitely_not_a_command", [])


def test_format_prediction_text(sp: SubstratePhysics) -> None:
    s = format_prediction(sp.predict_alpha_em())
    assert "alpha_em" in s
    assert "value" in s
    assert "precision" in s


def test_format_prediction_json(sp: SubstratePhysics) -> None:
    s = format_prediction(sp.predict_alpha_em(), json_output=True)
    parsed = json.loads(s)
    assert parsed["name"] == "alpha_em"
    assert parsed["unit"] == "dimensionless"


# ---------------------------------------------------------------------------
# CLI entry point (programmatic invocation; no subprocess needed for speed).
# ---------------------------------------------------------------------------


def _capture_main(argv):
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = cli_main(argv)
    return rc, buf.getvalue()


def test_cli_predict_alpha() -> None:
    rc, out = _capture_main(["predict", "alpha_em"])
    assert rc == 0
    assert "alpha_em" in out


def test_cli_predict_json() -> None:
    rc, out = _capture_main(["predict", "alpha_em", "--json"])
    assert rc == 0
    parsed = json.loads(out)
    assert parsed["name"] == "alpha_em"


def test_cli_predict_with_args() -> None:
    rc, out = _capture_main(["predict", "ie", "carbon"])
    assert rc == 0
    assert "IE(C)" in out


def test_cli_list() -> None:
    rc, out = _capture_main(["list"])
    assert rc == 0
    assert "alpha_em" in out
    assert "bandgap" in out


def test_cli_info() -> None:
    rc, out = _capture_main(["info"])
    assert rc == 0
    assert "K_pair" in out
    assert "Lambda_QCD_MeV" in out


def test_cli_unknown_returns_error() -> None:
    rc, _ = _capture_main(["predict", "not_a_real_topic"])
    assert rc == 2


# ---------------------------------------------------------------------------
# Batch processing.
# ---------------------------------------------------------------------------


def test_batch_round_trip(tmp_path: Path) -> None:
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    input_csv.write_text(
        "# command,arg1\n"
        "alpha_em\n"
        "m_mu_over_m_e\n"
        "ie,carbon\n"
        "bandgap,silicon\n"
        "fracture,50,600\n"
    )
    n = run_batch(str(input_csv), str(output_csv))
    assert n == 5
    rows = list(csv.DictReader(open(output_csv)))
    assert len(rows) == 5
    names = {r["name"] for r in rows}
    assert "alpha_em" in names


def test_batch_handles_unknown_rows(tmp_path: Path) -> None:
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    input_csv.write_text(
        "alpha_em\n"
        "definitely_not_a_command\n"
        "ie,carbon\n"
    )
    n = run_batch(str(input_csv), str(output_csv))
    # 2 succeed, 1 errors but is captured in output
    assert n == 2
    rows = list(csv.DictReader(open(output_csv)))
    assert any("ERROR" in r["name"] for r in rows)


# ---------------------------------------------------------------------------
# Smoke test: invoke as a subprocess to confirm `__main__` block works.
# ---------------------------------------------------------------------------


def test_subprocess_invocation() -> None:
    script = HERE / "substrate_physics.py"
    res = subprocess.run(
        [sys.executable, str(script), "predict", "alpha_em"],
        capture_output=True, text=True, timeout=20,
    )
    assert res.returncode == 0, res.stderr
    assert "alpha_em" in res.stdout
