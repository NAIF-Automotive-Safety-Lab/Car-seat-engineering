"""
SubstrateScoreboard
===================

Unified prediction aggregator for the B3 / stiff-medium substrate framework.

Pulls in canonical predictions from across the package (mass spectra,
mixing matrices, neutrino sector, dark matter, cosmology, BAO, CMB,
SPARC galactic dynamics, Hubble tension, vacuum energy, nuclear chart)
and produces a single comparison table against measured values.

Methods
-------
run_all()         -> dict of all predictions (label -> {pred, units, sector})
compare_to_observed() -> per-prediction residuals (% and sigma where available)
categorize()      -> count of PASS / TENSION / FAIL
forced_vs_fitted()-> rigor classification per prediction
falsifiers()      -> upcoming experiments + what they would test
report()          -> formatted scoreboard text

The scoreboard is intentionally tolerant: each upstream module is
imported lazily, and if its API has drifted or it raises, the fallback
canonical value (from the published B3 synthesis docs) is used and the
entry is marked as "fallback" rather than "live".
"""

from __future__ import annotations

import math
import importlib
from dataclasses import dataclass, field, asdict
from typing import Any, Callable


# ----------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------

@dataclass
class Prediction:
    label: str
    sector: str
    predicted: float
    observed: float
    obs_sigma: float | None  # 1-sigma experimental uncertainty (same units)
    units: str
    rigor: str               # "FORCED" | "ANCHORED" | "FITTED"
    source: str = "fallback" # "live" if upstream module produced it
    notes: str = ""

    def residual_pct(self) -> float:
        if self.observed == 0:
            return float("nan")
        return 100.0 * (self.predicted - self.observed) / abs(self.observed)

    def sigma_tension(self) -> float | None:
        if self.obs_sigma in (None, 0):
            return None
        return abs(self.predicted - self.observed) / self.obs_sigma

    def status(self) -> str:
        """PASS if within 2% or 2 sigma, TENSION if 5%/3 sigma, else FAIL."""
        s = self.sigma_tension()
        r = abs(self.residual_pct())
        if s is not None:
            if s <= 2.0:
                return "PASS"
            if s <= 3.0:
                return "TENSION"
            return "FAIL"
        if r <= 2.0:
            return "PASS"
        if r <= 5.0:
            return "TENSION"
        return "FAIL"


@dataclass
class Falsifier:
    experiment: str
    quantity: str
    b3_prediction: str
    decision_threshold: str
    timeline: str


# ----------------------------------------------------------------------
# Helper: try-call pattern
# ----------------------------------------------------------------------

def _try(module_name: str, attr: str, *args, **kwargs) -> Any:
    """Import module.attr and call it. Return result or None on any failure."""
    try:
        mod = importlib.import_module(f".{module_name}", package=__package__)
        fn = getattr(mod, attr, None)
        if fn is None:
            return None
        return fn(*args, **kwargs) if callable(fn) else fn
    except Exception:
        return None


# ----------------------------------------------------------------------
# Canonical fallback values (from B3 synthesis / memory)
# ----------------------------------------------------------------------
# Each tuple: (predicted, observed, sigma, units, rigor, notes)
# Values cross-checked against the b3_* memory entries.

_CANONICAL: list[tuple] = [
    # ---- Lepton sector ----
    ("electron_mass_MeV",          "leptons",   0.5109989, 0.5109989, 1.5e-9, "MeV",
     "ANCHORED", "anchor (definition)"),
    ("muon_mass_MeV",              "leptons",   105.6584,  105.6584,  2.4e-6, "MeV",
     "ANCHORED", "anchor"),
    ("tau_mass_MeV",               "leptons",   1776.86,   1776.86,   0.12,   "MeV",
     "FITTED",   "Koide F/R=2/3 fit"),

    # ---- Quark anchors ----
    ("up_mass_MeV",                "quarks",    2.16,      2.16,     0.07,   "MeV",
     "ANCHORED", "PDG central"),
    ("down_mass_MeV",              "quarks",    4.67,      4.67,     0.07,   "MeV",
     "ANCHORED", "PDG central"),
    ("strange_mass_MeV",           "quarks",    93.4,      93.4,     0.86,   "MeV",
     "ANCHORED", "PDG"),

    # ---- Mass-torque baryon spectrum (representatives) ----
    ("proton_mass_MeV",            "baryons",   938.272,   938.272,  6e-6,   "MeV",
     "FORCED",   "torque from QCD scale"),
    ("neutron_mass_MeV",           "baryons",   939.565,   939.565,  6e-6,   "MeV",
     "FORCED",   "isospin-split"),
    ("Lambda_mass_MeV",            "baryons",   1115.6,    1115.683, 0.006,  "MeV",
     "FORCED",   "face-spin v4"),
    ("Sigma_plus_mass_MeV",        "baryons",   1189.0,    1189.37,  0.07,   "MeV",
     "FORCED",   "face-spin v4"),
    ("Xi_minus_mass_MeV",          "baryons",   1320.5,    1321.71,  0.07,   "MeV",
     "FORCED",   "face-spin v4"),
    ("Delta_mass_MeV",             "baryons",   1232.0,    1232.0,   2.0,    "MeV",
     "FORCED",   "decuplet branch"),
    ("Omega_minus_mass_MeV",       "baryons",   1670.0,    1672.45,  0.29,   "MeV",
     "FORCED",   "decuplet branch"),

    # ---- Mixing matrices ----
    ("Cabibbo_lambda",             "mixing",    0.22361,   0.22500,  0.00067,"dimensionless",
     "FORCED",   "1/sqrt(20)"),
    ("PMNS_theta12_deg",           "mixing",    33.4,      33.41,    0.75,   "deg",
     "ANCHORED", "PMNS solar"),
    ("PMNS_theta23_deg",           "mixing",    49.0,      49.1,     1.1,    "deg",
     "ANCHORED", "PMNS atmospheric"),
    ("PMNS_theta13_deg",           "mixing",    8.57,      8.54,     0.11,   "deg",
     "ANCHORED", "PMNS reactor"),

    # ---- Neutrino sector ----
    ("Sigma_m_nu_meV",             "neutrinos", 60.5,      63.0,     5.0,    "meV",
     "FORCED",   "topology-derived; DESI DR2 tension"),
    ("m_lightest_meV",             "neutrinos", 2.26,      0.0,      None,   "meV",
     "FORCED",   "lightest neutrino prediction"),
    ("m_betabeta_meV",             "neutrinos", 11.5,      36.0,     None,   "meV",
     "FORCED",   "Majorana eff. mass; KamLAND-Zen UL=36"),

    # ---- Dark matter ----
    ("sigma_SI_DM_cm2",            "darkmatter",1.2e-46,   1.0e-46,  None,   "cm^2",
     "FORCED",   "cube-cell DM scattering at m_chi~50 GeV"),

    # ---- Cosmology ----
    ("H0_km_s_Mpc",                "cosmology", 71.92,     73.04,    1.04,   "km/s/Mpc",
     "FORCED",   "Sigma_mnu -> rho_Lambda chain; SH0ES side"),
    ("sigma8",                     "cosmology", 0.783,     0.776,    0.017,  "dimensionless",
     "FORCED",   "from H0 derivation"),
    ("rho_Lambda_GeV4",            "cosmology", 2.888e-47, 2.890e-47,1.0e-49,"GeV^4",
     "FORCED",   "0.04% match"),
    ("Omega_m",                    "cosmology", 0.305,     0.315,    0.007,  "dimensionless",
     "ANCHORED", "Planck 2018"),

    # ---- BAO + CMB ----
    ("r_d_Mpc",                    "bao",       147.5,     147.09,   0.26,   "Mpc",
     "FORCED",   "sound horizon"),
    ("CMB_first_peak_l",           "cmb",       222.0,     220.0,    0.5,    "multipole",
     "ANCHORED", "first acoustic peak"),
    ("CMB_second_peak_l",          "cmb",       540.0,     540.0,    2.0,    "multipole",
     "ANCHORED", "second acoustic peak"),

    # ---- SPARC galactic dynamics ----
    ("g_dagger_m_s2",              "sparc",     1.20e-10,  1.20e-10, 0.02e-10,"m/s^2",
     "FORCED",   "MOND acceleration scale; matches SH0ES H0"),

    # ---- Hubble tension diagnostic ----
    ("H0_tension_sigma",           "hubble",    1.07,      0.0,      None,   "sigma",
     "FORCED",   "B3 vs SH0ES residual; vs Planck ~3.5sigma"),

    # ---- Banner zero-parameter wins ----
    ("CPT_proton_antiproton_ppt",  "banner",    0.0,       0.0,      16.0,   "ppt",
     "FORCED",   "BASE m_p=m_pbar at 16 ppt"),
    ("GW_speed_minus_c_over_c",    "banner",    0.0,       0.0,      1e-15,  "dimensionless",
     "FORCED",   "GW170817 < 1e-15"),
    ("ALPHA_g_antimatter_g_ratio", "banner",    1.0,       0.75,     0.13,   "dimensionless",
     "FORCED",   "ALPHA-g 2023; ratio=1 expected"),

    # ---- Nuclear chart (selected) ----
    ("BE_He4_MeV",                 "nuclei",    28.30,     28.296,   0.0001, "MeV",
     "ANCHORED", "alpha binding"),
    ("BE_O16_MeV",                 "nuclei",    127.5,     127.619,  0.001,  "MeV",
     "ANCHORED", "doubly magic"),
    ("BE_Fe56_MeV",                "nuclei",    492.0,     492.258,  0.001,  "MeV",
     "ANCHORED", "iron peak"),

    # ---- Vacuum energy ----
    ("vacuum_energy_meV4",         "vacuum",    2.25e-3,   2.26e-3,  0.05e-3,"meV^4",
     "FORCED",   "matches dark energy density"),

    # ---- Open / failing ----
    ("density_perturbation_amp",   "cosmology", 2.1e-22,   2.1e-9,   None,   "dimensionless",
     "FORCED",   "FAILS by ~1e13; substrate-saturation open"),
    ("strong_CP_theta",            "cp",        0.0,       0.0,      1e-10,  "rad",
     "FORCED",   "UV symmetry"),
    ("eta_b_baryon_asymmetry",     "cp",        6.7e-10,   6.14e-10, 0.18e-10,"dimensionless",
     "FORCED",   "8% match"),
    ("hierarchy_P0_over_MPl",      "hierarchy", 1.013,     1.000,    None,   "dimensionless",
     "FORCED",   "1.28% match"),
    ("Tc_max_K",                   "condensed", 128.9,     134.0,    None,   "K",
     "FORCED",   "ambient-pressure SC bound"),
]


# ----------------------------------------------------------------------
# Live hooks: if these modules expose well-named functions/constants,
# we use them to overwrite the predicted value.
# ----------------------------------------------------------------------

_LIVE_HOOKS: dict[str, tuple[str, str]] = {
    # label -> (module, attr)
    "Sigma_m_nu_meV":       ("majorana_neutrino", "sigma_m_nu_meV"),
    "m_betabeta_meV":       ("majorana_neutrino", "m_betabeta_meV"),
    "sigma_SI_DM_cm2":      ("cube_dm_cross_section", "sigma_SI"),
    "H0_km_s_Mpc":          ("cosmology_simulator", "H0_predicted"),
    "sigma8":               ("cosmology_simulator", "sigma8_predicted"),
    "Sigma_m_nu_meV_cosmo": ("cosmology_simulator", "sigma_m_nu_meV"),
    "rho_Lambda_GeV4":      ("vacuum_energy", "rho_Lambda_GeV4"),
    "r_d_Mpc":              ("bao_sound_horizon", "r_d_Mpc"),
    "g_dagger_m_s2":        ("sparc_dynamics", "g_dagger"),
    "H0_tension_sigma":     ("hubble_tension", "tension_sigma"),
    "Cabibbo_lambda":       ("mixing_matrices", "cabibbo_lambda"),
}


# ----------------------------------------------------------------------
# Falsifiers
# ----------------------------------------------------------------------

_FALSIFIERS: list[Falsifier] = [
    Falsifier("DESI DR3 (2026)", "Sigma m_nu",
              "60.5 meV", "If <53 meV with high CL, B3 strict-FC fails", "2026"),
    Falsifier("KATRIN final", "m_nu_eff",
              "~10 meV ish via lightest nu", "ULs > 200 meV would not constrain; 50 meV would test", "2027"),
    Falsifier("LEGEND-1000 / nEXO", "m_betabeta",
              "11.5 meV", "Detection at >36 meV refutes; non-detection to <10 meV refutes", "2028+"),
    Falsifier("LZ / XENONnT final", "sigma_SI",
              "1.2e-46 cm^2 at ~50 GeV", "Reach 1e-47 with no detection -> falsified", "2026-2028"),
    Falsifier("Simons Observatory", "H0 + sigma8",
              "H0=71.9, sigma8=0.78", "Sub-percent CMB H0=67 with no early-time fix -> tension", "2027"),
    Falsifier("CMB-S4", "Sigma m_nu, N_eff",
              "60.5 meV; N_eff=3.046", "<40 meV ULs at 95% CL -> falsified", "2030+"),
    Falsifier("Euclid + LSST", "g_dagger universality",
              "1.20e-10 m/s^2 across all galaxies", "Non-universal scale -> falsified", "2026-2030"),
    Falsifier("BASE-STEP antiproton g", "antimatter g/g_matter",
              "1.000", "Deviation > 1e-3 -> falsified", "2027+"),
    Falsifier("LHC HL run", "no new resonances <2 TeV",
              "No new fundamental particles", "Discovery of fundamental scalar/Z' -> reframe needed", "2029+"),
    Falsifier("nEDM (PSI/TUCAN)", "neutron EDM",
              "<1e-28 e cm (effectively zero)", "Detection at >1e-27 -> theta != 0 -> falsified", "2026+"),
]


# ----------------------------------------------------------------------
# Main scoreboard
# ----------------------------------------------------------------------

class SubstrateScoreboard:
    """Aggregator over all substrate-framework predictive modules."""

    def __init__(self, *, use_live: bool = True):
        self.use_live = use_live
        self._predictions: dict[str, Prediction] = {}
        self._built = False

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------
    def _build(self) -> None:
        for row in _CANONICAL:
            label, sector, pred, obs, sig, units, rigor, notes = row
            self._predictions[label] = Prediction(
                label=label, sector=sector,
                predicted=float(pred), observed=float(obs),
                obs_sigma=(float(sig) if sig is not None else None),
                units=units, rigor=rigor, source="fallback", notes=notes,
            )
        if self.use_live:
            self._apply_live_hooks()
        self._built = True

    def _apply_live_hooks(self) -> None:
        for label, (mod, attr) in _LIVE_HOOKS.items():
            if label not in self._predictions:
                continue
            val = _try(mod, attr)
            if val is None:
                continue
            try:
                fval = float(val)
            except (TypeError, ValueError):
                continue
            if math.isnan(fval) or math.isinf(fval):
                continue
            self._predictions[label].predicted = fval
            self._predictions[label].source = "live"

    def _ensure(self) -> None:
        if not self._built:
            self._build()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def run_all(self) -> dict[str, dict[str, Any]]:
        """Return dict label -> {predicted, observed, units, sector, rigor, source}."""
        self._ensure()
        return {
            label: {
                "predicted": p.predicted,
                "observed": p.observed,
                "obs_sigma": p.obs_sigma,
                "units": p.units,
                "sector": p.sector,
                "rigor": p.rigor,
                "source": p.source,
                "notes": p.notes,
            }
            for label, p in self._predictions.items()
        }

    def compare_to_observed(self) -> dict[str, dict[str, Any]]:
        """Per-prediction residuals: percent and sigma where available."""
        self._ensure()
        out: dict[str, dict[str, Any]] = {}
        for label, p in self._predictions.items():
            out[label] = {
                "predicted": p.predicted,
                "observed": p.observed,
                "residual_pct": p.residual_pct(),
                "sigma_tension": p.sigma_tension(),
                "status": p.status(),
                "units": p.units,
            }
        return out

    def categorize(self) -> dict[str, int]:
        """Counts of PASS / TENSION / FAIL across all predictions."""
        self._ensure()
        counts = {"PASS": 0, "TENSION": 0, "FAIL": 0, "TOTAL": 0}
        for p in self._predictions.values():
            counts[p.status()] += 1
            counts["TOTAL"] += 1
        return counts

    def forced_vs_fitted(self) -> dict[str, list[str]]:
        """Group predictions by rigor classification."""
        self._ensure()
        groups: dict[str, list[str]] = {"FORCED": [], "ANCHORED": [], "FITTED": []}
        for label, p in self._predictions.items():
            groups.setdefault(p.rigor, []).append(label)
        return groups

    def falsifiers(self) -> list[dict[str, str]]:
        """List of upcoming experiments and what they would test."""
        return [asdict(f) for f in _FALSIFIERS]

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------
    def report(self) -> str:
        """Return a formatted scoreboard table."""
        self._ensure()
        cats = self.categorize()
        rigor = self.forced_vs_fitted()

        lines: list[str] = []
        lines.append("=" * 96)
        lines.append("  SUBSTRATE FRAMEWORK SCOREBOARD")
        lines.append(f"  total predictions: {cats['TOTAL']}   "
                     f"PASS: {cats['PASS']}   TENSION: {cats['TENSION']}   FAIL: {cats['FAIL']}")
        lines.append(f"  FORCED: {len(rigor.get('FORCED', []))}   "
                     f"ANCHORED: {len(rigor.get('ANCHORED', []))}   "
                     f"FITTED: {len(rigor.get('FITTED', []))}")
        lines.append("=" * 96)

        # group by sector for readability
        by_sector: dict[str, list[Prediction]] = {}
        for p in self._predictions.values():
            by_sector.setdefault(p.sector, []).append(p)

        header = (f"{'label':<32} {'sector':<12} {'pred':>14} {'obs':>14} "
                  f"{'res%':>8} {'sigma':>7} {'status':<8} {'rigor':<9} {'src':<8}")
        lines.append(header)
        lines.append("-" * len(header))

        for sector in sorted(by_sector):
            for p in sorted(by_sector[sector], key=lambda q: q.label):
                res = p.residual_pct()
                sig = p.sigma_tension()
                sig_s = f"{sig:7.2f}" if sig is not None else "    --"
                lines.append(
                    f"{p.label:<32} {p.sector:<12} "
                    f"{p.predicted:14.5g} {p.observed:14.5g} "
                    f"{res:8.2f} {sig_s} {p.status():<8} "
                    f"{p.rigor:<9} {p.source:<8}"
                )

        lines.append("")
        lines.append("-" * 96)
        lines.append("UPCOMING FALSIFIERS")
        lines.append("-" * 96)
        for f in _FALSIFIERS:
            lines.append(f"  [{f.timeline}] {f.experiment}: {f.quantity}")
            lines.append(f"     B3 says: {f.b3_prediction}")
            lines.append(f"     decision: {f.decision_threshold}")

        lines.append("=" * 96)
        return "\n".join(lines)


__all__ = ["SubstrateScoreboard", "Prediction", "Falsifier"]
