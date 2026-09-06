"""Substrate sound-velocity prediction vs measured longitudinal sound speeds.

WHAT THE TEST IS ACTUALLY MEASURING
-----------------------------------
The substrate Lagrangian L = (ρ/2)(∂_t φ)² - (K/2)(∇φ)² gives a single
universal long-wavelength dispersion ω(k) = c_s · k with c_s = sqrt(K/ρ).
For an isotropic three-dimensional elastic medium this generalises to two
acoustic branches:

    longitudinal :  c_L = sqrt((B + 4G/3) / ρ)        (P-wave / compression)
    transverse   :  c_T = sqrt(G / ρ)                 (S-wave / shear)

This module asks the SHARP question: when we feed the substrate Lagrangian's
own (B, G) prediction (substrate_elasticity.bulk_modulus_substrate /
shear_modulus_substrate) into the longitudinal sound-speed formula, do we
recover the measured c_L across a benchmark set of 15 materials spanning
~12× in c_L (water 1480 m/s -> diamond ≈ 18000 m/s)?

This is a multiplicative test of the Lagrangian's elastic structure: BOTH
the moduli AND the wave-equation kinematic content must be consistent with
measured sound speeds. A failure here would falsify either the K_4 face-pair
spring construction or the substrate Lagrangian's kinematic content (or
both); a success means the chain

        K_4 face-pair coupling -> (B, G) -> c_L = sqrt((B + 4G/3)/ρ)

is consistent with measured P-wave speeds across metals, covalent solids,
glasses, and a liquid.

HONEST SCOPING
--------------
The 15-material benchmark has TWO regimes:

  (A) Materials in the substrate_elasticity database (8 materials):
      Cu, Al, Au, Fe, Ni, Pb, diamond, Si — for these, the substrate
      prediction uses the substrate-derived (B, G) from K_4 face-pair
      coupling with ε_coh as one per-element handbook anchor (Cat-B).

  (B) Materials NOT in the substrate_elasticity database (7 materials):
      Be, Ag, W, Ge, Steel (~Fe), glass (silica), quartz, water — for these
      the substrate elasticity has not been derived in this version; we
      report the prediction using MEASURED (B, G) so the Lagrangian's
      kinematic content is still tested in isolation.

Both regimes are reported with explicit flags ('uses_substrate_moduli'
True/False).

CAVEATS
-------
* Water is a LIQUID — it has no shear modulus (G = 0); only c_L is meaningful
  and the substrate central-force model is not designed for liquids. We
  include it as a kinematic-only check (c_L = sqrt(B/ρ)).
* Glass (fused silica) and quartz are AMORPHOUS / non-cubic; the substrate
  central-force model is averaged.
* Steel is NOT pure iron; its measured c_L (≈5930 m/s) is dominated by Fe
  but with carbide microstructure; we report it tagged as such.
* For diamond and Si, the central-force substrate prediction underpredicts
  G/B (3/5 = 0.6 vs measured 1.08 for diamond) — this is the documented
  Cauchy-violation regime for strong angular bonds. Hence we ALSO compare
  the substrate-prediction c_L to the measured-moduli c_L to separate the
  Lagrangian's wave-equation content from the moduli accuracy.

References
----------
* Ashcroft & Mermin, Solid State Physics, Ch. 22-23 (sound speeds, Debye).
* Auld, Acoustic Fields and Waves in Solids (1990), Vol. 1, Ch. 3.
* CRC Handbook of Chemistry and Physics, 90th ed., Section 14 (acoustic
  properties of materials), and Section 12 (elastic constants).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np


# --------------------------------------------------------------------------- #
# Physical constants (SI)                                                     #
# --------------------------------------------------------------------------- #

PASCAL_PER_GPA: float = 1.0e9


# --------------------------------------------------------------------------- #
# Material database                                                           #
# --------------------------------------------------------------------------- #
#
# Per-material entry:
#   M_molar      : molar mass [kg/mol]   (mostly informational; NOT used
#                                          for c_L because ρ already given)
#   rho_mass     : mass density [kg/m³]  (CRC, room temperature)
#   B_GPa_meas   : measured bulk modulus [GPa]  (CRC 90th ed., Section 12)
#   G_GPa_meas   : measured shear modulus [GPa] (CRC 90th ed., Section 12;
#                                                 0 for liquids)
#   c_L_meas     : measured longitudinal sound speed [m/s] (CRC Section 14;
#                                                            user spec)
#   in_elasticity_db : whether substrate_elasticity has a substrate
#                      prediction for (B, G); when True we use the
#                      substrate-derived (B, G); when False we fall back
#                      to measured (B, G) and flag the row.
#   notes        : short caveat / classification
#
# All measured longitudinal sound-speed values are exactly the user-supplied
# spec values; the moduli/density values are CRC handbook unless flagged.

@dataclass(frozen=True)
class SoundMaterial:
    """Per-material data for the sound-velocity comparison."""
    name:             str
    rho_mass:         float            # kg/m³
    B_GPa_meas:       float            # measured bulk modulus
    G_GPa_meas:       float            # measured shear modulus (0 = liquid)
    c_L_meas:         float            # measured longitudinal sound speed (m/s)
    in_elasticity_db: bool             # True if substrate_elasticity has it
    notes:            str = ""

    @property
    def c_L_measured_moduli(self) -> float:
        """Reference c_L computed from the measured (B, G) and ρ.

        This is the 'kinematic-content' check: with measured moduli, does
        the substrate Lagrangian's wave equation reproduce the directly
        measured P-wave speed? For materials with G > 0 the deviation
        between this and ``c_L_meas`` measures both rho/moduli table
        consistency and any anisotropy averaging error.
        """
        return longitudinal_speed(self.B_GPa_meas, self.G_GPa_meas, self.rho_mass)


MATERIALS: Dict[str, SoundMaterial] = {
    # ----- Materials covered by substrate_elasticity (substrate-predicted B,G) ----- #
    "Diamond":   SoundMaterial(
        "Diamond",   rho_mass=3515.0,  B_GPa_meas=442.0, G_GPa_meas=478.0,
        c_L_meas=18000.0, in_elasticity_db=True,
        notes="diamond cubic, very strong angular bonds (G/B=1.08 >> 0.6)",
    ),
    "Aluminum":  SoundMaterial(
        "Aluminum",  rho_mass=2700.0,  B_GPa_meas= 76.0, G_GPa_meas= 26.0,
        c_L_meas=6420.0,  in_elasticity_db=True,
        notes="FCC metal, central-force regime",
    ),
    "Iron":      SoundMaterial(
        "Iron",      rho_mass=7874.0,  B_GPa_meas=170.0, G_GPa_meas= 82.0,
        c_L_meas=5960.0,  in_elasticity_db=True,
        notes="BCC metal, near-central force",
    ),
    "Copper":    SoundMaterial(
        "Copper",    rho_mass=8960.0,  B_GPa_meas=140.0, G_GPa_meas= 48.0,
        c_L_meas=4760.0,  in_elasticity_db=True,
        notes="FCC metal, near-central force",
    ),
    "Gold":      SoundMaterial(
        "Gold",      rho_mass=19300.0, B_GPa_meas=180.0, G_GPa_meas= 27.0,
        c_L_meas=3240.0,  in_elasticity_db=True,
        notes="FCC metal, soft (G/B=0.15 << 0.6)",
    ),
    "Lead":      SoundMaterial(
        "Lead",      rho_mass=11340.0, B_GPa_meas= 46.0, G_GPa_meas=  5.6,
        c_L_meas=2160.0,  in_elasticity_db=True,
        notes="FCC metal, very soft (G/B=0.12, strong non-central)",
    ),
    "Silicon":   SoundMaterial(
        "Silicon",   rho_mass=2329.0,  B_GPa_meas= 98.0, G_GPa_meas= 66.0,
        c_L_meas=8433.0,  in_elasticity_db=True,
        notes="diamond cubic covalent, angular bonds (G/B=0.67)",
    ),

    # ----- Materials NOT in substrate_elasticity (measured B,G fallback) ----- #
    "Beryllium": SoundMaterial(
        "Beryllium", rho_mass=1850.0,  B_GPa_meas=130.0, G_GPa_meas=132.0,
        c_L_meas=12890.0, in_elasticity_db=False,
        notes="HCP metal, low-density / very stiff (high c_L)",
    ),
    "Silver":    SoundMaterial(
        "Silver",    rho_mass=10490.0, B_GPa_meas=100.0, G_GPa_meas= 30.0,
        c_L_meas=3650.0,  in_elasticity_db=False,
        notes="FCC metal, soft (G/B=0.30)",
    ),
    "Tungsten":  SoundMaterial(
        "Tungsten",  rho_mass=19250.0, B_GPa_meas=310.0, G_GPa_meas=161.0,
        c_L_meas=5180.0,  in_elasticity_db=False,
        notes="BCC metal, very dense + very stiff (cancels in c_L)",
    ),
    "Germanium": SoundMaterial(
        "Germanium", rho_mass=5323.0,  B_GPa_meas= 75.0, G_GPa_meas= 55.0,
        c_L_meas=5400.0,  in_elasticity_db=False,
        notes="diamond cubic covalent, angular bonds (G/B=0.73)",
    ),
    "Steel":     SoundMaterial(
        "Steel",     rho_mass=7850.0,  B_GPa_meas=160.0, G_GPa_meas= 79.3,
        c_L_meas=5930.0,  in_elasticity_db=False,
        notes="iron alloy (~mild steel); not pure Fe",
    ),
    "Quartz":    SoundMaterial(
        "Quartz",    rho_mass=2648.0,  B_GPa_meas= 37.5, G_GPa_meas= 44.3,
        c_L_meas=5750.0,  in_elasticity_db=False,
        notes="α-quartz crystal (anisotropic, Voigt-Reuss-Hill avg)",
    ),
    "Glass":     SoundMaterial(
        "Glass",     rho_mass=2200.0,  B_GPa_meas= 36.4, G_GPa_meas= 31.2,
        c_L_meas=5968.0,  in_elasticity_db=False,
        notes="fused silica (amorphous SiO_2)",
    ),
    "Water":     SoundMaterial(
        "Water",     rho_mass=1000.0,  B_GPa_meas=  2.2, G_GPa_meas=  0.0,
        c_L_meas=1480.0,  in_elasticity_db=False,
        notes="liquid (G=0): c_L = sqrt(B/ρ); central-force model n/a",
    ),
}


# --------------------------------------------------------------------------- #
# Wave-speed kinematics from the substrate Lagrangian                         #
# --------------------------------------------------------------------------- #


def longitudinal_speed(B_GPa: float, G_GPa: float, rho: float) -> float:
    """Longitudinal (P-wave) sound speed from the substrate wave equation.

        c_L = sqrt((B + 4G/3) / ρ)

    Inputs in GPa, kg/m³. Output in m/s. Reduces to c_L = sqrt(B/ρ) for
    liquids (G = 0).
    """
    if rho <= 0.0:
        raise ValueError("rho must be positive")
    if B_GPa < 0.0 or G_GPa < 0.0:
        raise ValueError("B and G must be non-negative")
    M_pa = (B_GPa + (4.0 / 3.0) * G_GPa) * PASCAL_PER_GPA
    return math.sqrt(M_pa / rho)


def transverse_speed(G_GPa: float, rho: float) -> float:
    """Transverse (S-wave) sound speed: c_T = sqrt(G/ρ).

    Returns 0 for G = 0 (liquids do not support shear waves).
    """
    if rho <= 0.0:
        raise ValueError("rho must be positive")
    if G_GPa < 0.0:
        raise ValueError("G must be non-negative")
    return math.sqrt(G_GPa * PASCAL_PER_GPA / rho)


# --------------------------------------------------------------------------- #
# Substrate (B, G) lookup                                                     #
# --------------------------------------------------------------------------- #


def substrate_moduli(name: str) -> Tuple[Optional[float], Optional[float], bool]:
    """Return (B_GPa, G_GPa, used_substrate_flag) for material ``name``.

    Looks up the material in substrate_elasticity.MATERIALS and returns the
    substrate-derived (B, G) from K_4 face-pair coupling.  For materials not
    in the elasticity database, returns (None, None, False).
    """
    try:
        from .substrate_elasticity import (
            MATERIALS as ELASTIC_MATS,
            bulk_modulus_substrate,
            shear_modulus_substrate,
        )
    except ImportError:
        return None, None, False

    if name not in ELASTIC_MATS:
        return None, None, False

    emat = ELASTIC_MATS[name]
    return bulk_modulus_substrate(emat), shear_modulus_substrate(emat), True


# --------------------------------------------------------------------------- #
# Per-material prediction                                                     #
# --------------------------------------------------------------------------- #


def predict_sound_speed(mat: SoundMaterial) -> Dict[str, float]:
    """Run the substrate prediction for one material.

    Returns a dict with substrate-predicted c_L (using substrate (B, G) when
    available, else falling back to measured (B, G)), the measured-moduli
    reference c_L, and the directly measured c_L from the user spec.
    """
    B_sub, G_sub, used_substrate = substrate_moduli(mat.name)
    if not used_substrate:
        B_use = mat.B_GPa_meas
        G_use = mat.G_GPa_meas
    else:
        assert B_sub is not None and G_sub is not None
        B_use = B_sub
        G_use = G_sub

    c_L_pred = longitudinal_speed(B_use, G_use, mat.rho_mass)
    c_T_pred = transverse_speed(G_use, mat.rho_mass)
    c_L_meas_moduli = mat.c_L_measured_moduli
    c_L_meas = mat.c_L_meas

    rel_err = (c_L_pred - c_L_meas) / c_L_meas
    rel_err_kinematic = (c_L_meas_moduli - c_L_meas) / c_L_meas

    return {
        "name":               mat.name,
        "rho_mass":           mat.rho_mass,
        "B_GPa_used":         B_use,
        "G_GPa_used":         G_use,
        "B_GPa_meas":         mat.B_GPa_meas,
        "G_GPa_meas":         mat.G_GPa_meas,
        "used_substrate_BG":  used_substrate,
        "c_L_pred":           c_L_pred,
        "c_T_pred":           c_T_pred,
        "c_L_measured_moduli": c_L_meas_moduli,
        "c_L_meas":           c_L_meas,
        "rel_err":            rel_err,
        "rel_err_kinematic":  rel_err_kinematic,
        "notes":              mat.notes,
    }


# --------------------------------------------------------------------------- #
# Run-test driver                                                             #
# --------------------------------------------------------------------------- #


def _log_pearson(pred: np.ndarray, meas: np.ndarray) -> float:
    """Pearson correlation in log-log space (sound speeds vary ~12×)."""
    lp = np.log(pred)
    lm = np.log(meas)
    lpm = lp - lp.mean()
    lmm = lm - lm.mean()
    denom = math.sqrt(float((lpm * lpm).sum()) * float((lmm * lmm).sum()))
    return float((lpm * lmm).sum()) / denom if denom > 0.0 else float("nan")


def run_test() -> Dict[str, object]:
    """Run substrate sound-velocity predictions on all 15 materials.

    Returns a dict with:
        rows    : per-material prediction + measurement record
        summary : aggregate statistics for the FULL 15-material set
        summary_substrate_BG : aggregate statistics restricted to the
                               materials with substrate-derived (B, G)
                               (a sharper test of the substrate stack)
    """
    rows: Dict[str, Dict[str, float]] = {}
    pred_arr: List[float] = []
    meas_arr: List[float] = []
    pred_arr_sub: List[float] = []
    meas_arr_sub: List[float] = []
    pred_arr_kin: List[float] = []      # measured-moduli c_L
    meas_arr_kin: List[float] = []      # measured c_L

    for name, mat in MATERIALS.items():
        row = predict_sound_speed(mat)
        rows[name] = row
        pred_arr.append(row["c_L_pred"])
        meas_arr.append(row["c_L_meas"])
        pred_arr_kin.append(row["c_L_measured_moduli"])
        meas_arr_kin.append(row["c_L_meas"])
        if row["used_substrate_BG"]:
            pred_arr_sub.append(row["c_L_pred"])
            meas_arr_sub.append(row["c_L_meas"])

    pred = np.asarray(pred_arr)
    meas = np.asarray(meas_arr)
    rel = (pred - meas) / meas

    pred_kin = np.asarray(pred_arr_kin)
    meas_kin = np.asarray(meas_arr_kin)
    rel_kin = (pred_kin - meas_kin) / meas_kin

    pred_s = np.asarray(pred_arr_sub)
    meas_s = np.asarray(meas_arr_sub)
    rel_s = (pred_s - meas_s) / meas_s

    summary = {
        "n_materials":              len(rows),
        "mean_abs_rel_err":         float(np.abs(rel).mean()),
        "median_abs_rel_err":       float(np.median(np.abs(rel))),
        "max_abs_rel_err":          float(np.abs(rel).max()),
        "within_10pct":             int(sum(1 for r in rel if abs(r) <= 0.10)),
        "within_25pct":             int(sum(1 for r in rel if abs(r) <= 0.25)),
        "within_50pct":             int(sum(1 for r in rel if abs(r) <= 0.50)),
        "loglog_pearson":           _log_pearson(pred, meas),
    }

    # Restricted to the materials with substrate-derived (B, G):
    summary_substrate_BG = {
        "n_materials_substrate_BG": int(len(pred_s)),
        "mean_abs_rel_err":         float(np.abs(rel_s).mean()) if len(pred_s) else float("nan"),
        "median_abs_rel_err":       float(np.median(np.abs(rel_s))) if len(pred_s) else float("nan"),
        "max_abs_rel_err":          float(np.abs(rel_s).max()) if len(pred_s) else float("nan"),
        "within_25pct":             int(sum(1 for r in rel_s if abs(r) <= 0.25)),
        "within_50pct":             int(sum(1 for r in rel_s if abs(r) <= 0.50)),
        "loglog_pearson":           _log_pearson(pred_s, meas_s) if len(pred_s) >= 3 else float("nan"),
    }

    # Kinematic-only test (measured moduli + substrate Lagrangian wave eq):
    summary_kinematic = {
        "n_materials":              len(rows),
        "mean_abs_rel_err":         float(np.abs(rel_kin).mean()),
        "median_abs_rel_err":       float(np.median(np.abs(rel_kin))),
        "max_abs_rel_err":          float(np.abs(rel_kin).max()),
        "within_10pct":             int(sum(1 for r in rel_kin if abs(r) <= 0.10)),
        "loglog_pearson":           _log_pearson(pred_kin, meas_kin),
    }

    return {
        "rows":                  rows,
        "summary":               summary,
        "summary_substrate_BG":  summary_substrate_BG,
        "summary_kinematic":     summary_kinematic,
    }


# --------------------------------------------------------------------------- #
# CLI entrypoint                                                              #
# --------------------------------------------------------------------------- #


def main() -> None:
    res = run_test()
    rows = res["rows"]
    print("Substrate sound-velocity predictions (c_L) vs measured P-wave speeds")
    print("=" * 105)
    print(
        f"{'Material':<11}{'sub?':>5}{'rho':>8}"
        f"{'B_use':>8}{'G_use':>8}"
        f"{'c_L_pred':>10}{'c_L_kin':>10}{'c_L_meas':>10}"
        f"{'rel_err':>10}{'rel_kin':>10}"
    )
    for name, row in rows.items():
        sub_str = "YES" if row["used_substrate_BG"] else "no"
        print(
            f"{name:<11}{sub_str:>5}{row['rho_mass']:>8.0f}"
            f"{row['B_GPa_used']:>8.1f}{row['G_GPa_used']:>8.1f}"
            f"{row['c_L_pred']:>10.0f}{row['c_L_measured_moduli']:>10.0f}"
            f"{row['c_L_meas']:>10.0f}"
            f"{row['rel_err']:>+10.2%}{row['rel_err_kinematic']:>+10.2%}"
        )

    print()
    print("Summary (all 15 materials, substrate prediction)")
    for k, v in res["summary"].items():
        if isinstance(v, float):
            print(f"  {k:<28} {v:+.5f}")
        else:
            print(f"  {k:<28} {v}")

    print()
    print("Summary (8 materials with substrate-derived B, G)")
    for k, v in res["summary_substrate_BG"].items():
        if isinstance(v, float):
            print(f"  {k:<28} {v:+.5f}")
        else:
            print(f"  {k:<28} {v}")

    print()
    print("Summary (kinematic check: measured B,G + substrate wave eq)")
    for k, v in res["summary_kinematic"].items():
        if isinstance(v, float):
            print(f"  {k:<28} {v:+.5f}")
        else:
            print(f"  {k:<28} {v}")


if __name__ == "__main__":
    main()
