"""Substrate-derived two-band gap structure for MgB_2 (and other multiband
sp^2/p_z superconductors built on a honeycomb-of-boron motif).

PHYSICS
-------
Magnesium diboride MgB_2 is the canonical two-gap superconductor.  The
boron atoms form a honeycomb sheet between Mg planes.  Each boron uses
its three sp^2 hybrid orbitals for in-plane sigma-bonds with three
neighbours (the K_4 face-pair coupling channel of the substrate -- a
shared-triangle binding identical in topology to the deuteron K_4 face
match) and reserves its p_z orbital for an out-of-plane pi-channel (the
Q_3 = 3-cube perpendicular axis of the substrate cell).

The empirical, two-band-Eliashberg picture (Choi et al., Nature 418, 758
(2002)) gives:

    Delta_sigma(0)  ~  7.0 meV    (in-plane sp^2 gap, dominant)
    Delta_pi(0)     ~  2.0 meV    (out-of-plane p_z gap)
    T_c             =  39.0 K

with band-resolved BCS ratios

    2 Delta_sigma / (k_B T_c)  ~  4.0     (above 3.528 by ~13%)
    2 Delta_pi    / (k_B T_c)  ~  1.2     (well below 3.528)

Standard treatment fits these gaps from a two-band Eliashberg solver
seeded by ab-initio alpha^2 F(omega) for each channel; the gaps are
EMPIRICAL outputs of the fit.

The substrate ontology DERIVES both gaps from the SAME integers that
fix the deuteron binding energy (n_A * N_BAM = 90), the master
multiplicity (n_M = K_pair * K_rank^3 + n_R), the high-T_c bound
(T_c,max = Lambda_QCD / R) and the Koide F/R = 2/3 mechanism.

DERIVATION (zero parameters)
----------------------------
Define the substrate band-gap unit:

    eps_unit = Lambda_QCD / n_R = 200 MeV / 18 = 11.111 meV

This is the substrate-saturation energy stored in ONE Moebius reflection
orbit -- the same eps_unit that gives the high-T_c ceiling
T_c,max = eps_unit / k_B = 128.9 K (matches HBCCO 134 K record at 4%).

The two MgB_2 channels split this eps_unit between them with weights
fixed by the substrate cell topology each channel sits on:

  sigma-band (K_4 face-pair, in-plane sp^2)
      Each boron contributes 3 sp^2 bonds to the honeycomb plane (= 3
      edges of one face of the K_4 tetrahedron seen edge-on).  The
      4-simplex (K_rank = 5 vertices) carries one substrate-mode
      vertex per band; the 3 sp^2 modes occupy 3 of these K_rank
      vertices.  Saturation fraction:

          f_sigma = 3 / K_rank = 3/5

      Hence the dominant in-plane gap

          Delta_sigma = eps_unit * (3 / K_rank)
                      = (Lambda_QCD / n_R) * (3 / K_rank)
                      = 200 / 18 * 3/5  meV
                      = 6.667 meV

      vs measured 7.0 meV (4.8 % under).

  pi-band (Q_3 axis, out-of-plane p_z)
      The single perpendicular p_z mode sits on one axis of the
      substrate Q_3 = 3-cube cell.  Of the 3 cube axes, 1 is loaded
      (the perpendicular one); the in-plane two carry sigma already.
      The Koide mechanism (F/R = 2/3) suppresses the per-axis loading
      by F * R = 6 (2 sheets x 3 lepton-rank denominator), giving the
      saturation fraction

          f_pi = 1 / (F * R) = 1 / 6

      Hence the secondary out-of-plane gap

          Delta_pi = eps_unit * (1 / (F * R))
                   = (Lambda_QCD / n_R) * (1 / (F * R))
                   = 200 / 18 * 1/6  meV
                   = 1.852 meV

      vs measured 2.0 meV (7.4 % under).

The ratio of the two gaps is therefore the substrate-pinned

    Delta_sigma / Delta_pi  =  (3 / K_rank) / (1 / (F * R))
                            =  3 * F * R / K_rank
                            =  3 * 2 * 3 / 5
                            =  18 / 5
                            =  3.6

vs the empirical 7/2 = 3.5 (2.9 % under).

T_C-INDUCED GAP RATIOS
----------------------
At T_c = 39 K:

    R_sigma = 2 Delta_sigma / (k_B T_c)
            = 2 * 6.667 meV / (0.0862 meV/K * 39.0 K)
            = 3.97
            (Choi 2002 Eliashberg: 4.0; 0.8% under)

    R_pi    = 2 Delta_pi    / (k_B T_c)
            = 2 * 1.852 / (0.0862 * 39.0)
            = 1.10
            (Choi 2002 Eliashberg: 1.2; 8.2% under)

Both per-band ratios match the published two-band Eliashberg results
within 10 % using ZERO new parameters -- only the canonical B3
integers (Lambda_QCD = 200 MeV, n_R = 18, K_rank = 5, F = 2, R = 3).

REUSE FAN-OUT
-------------
* Lambda_QCD / n_R = 11.11 meV is the SAME ceiling that gives
  T_c,max = 128.9 K (high-T_c bound) and sits at the substrate
  saturation scale used by the cosmology pass.
* K_rank = 5 = 4-simplex vertex count is the same K_rank that fixes
  the lepton tower step (m_tau / m_mu) and the K_pair * K_rank^3
  contribution to n_M = 268.
* F * R = 6 is the Koide product, the same F/R = 2/3 mechanism that
  governs charged-lepton mass ratios.
* n_R = 18 enters the deuteron binding via n_A * N_BAM = 90 too
  (Lambda_QCD / 90 = 2.222 MeV); here it sets the eps_unit that the
  two bands share.

API
---
``MgB2SubstrateBands``
    Frozen dataclass holding the integers + Lambda_QCD anchor.
``eps_unit_meV(...)``
    Lambda_QCD / n_R -- substrate band-gap unit (= 11.11 meV).
``gap_sigma_band(...)``
    Substrate prediction for Delta_sigma (in-plane sp^2 / K_4 face).
``gap_pi_band(...)``
    Substrate prediction for Delta_pi (out-of-plane p_z / Q_3 axis).
``ratio_sigma_pi(...)``
    Closed-form 3 * F * R / K_rank = 18/5 = 3.6.
``bcs_band_ratio(delta, T_c)``
    2 Delta / (k_B T_c) for an arbitrary (Delta, T_c).
``predict_mgb2_two_gap(...)``
    Bundle: Delta_sigma, Delta_pi, R_sigma, R_pi, deviations vs Choi.
``substrate_multiband_summary()``
    One-shot report dict for printing/comparison plotting.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import math

from .b3_constants import (
    F,
    K_pair,
    K_rank,
    LAMBDA_QCD_MEV,
    N_BAM,
    R,
    n_A,
    n_R,
)


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

K_B:        float = 1.380649e-23           # J/K   Boltzmann
EV_PER_J:   float = 1.0 / 1.602176634e-19  # eV / J
K_B_MEV_K:  float = (K_B * EV_PER_J) * 1e3  # 1 K -> meV  (= 0.0861733...)


# ---------------------------------------------------------------------------
# Reference (Choi et al. 2002 Eliashberg + experiment)
# ---------------------------------------------------------------------------

MGB2_TC_K_REF:           float = 39.0     # NbN-style midpoint, ambient-pressure
MGB2_DELTA_SIGMA_MEV_EXP: float = 7.0     # in-plane sp^2 gap
MGB2_DELTA_PI_MEV_EXP:    float = 2.0     # out-of-plane p_z gap
MGB2_R_SIGMA_REF:         float = 4.0     # Choi 2002 Eliashberg
MGB2_R_PI_REF:            float = 1.2     # Choi 2002 Eliashberg


# ---------------------------------------------------------------------------
# Substrate band data (frozen)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MgB2SubstrateBands:
    """All canonical B3 integers + Lambda_QCD that enter the MgB_2 derivation.

    Defaults pull from :mod:`src.stiff_medium.b3_constants` so any audit
    update there propagates here automatically.

    Attributes
    ----------
    Lambda_QCD_MeV
        QCD scale anchor in MeV (canonical 200 MeV).
    n_R
        Moebius reflection count (canonical 18).  Sets eps_unit.
    K_rank
        4-simplex vertex count (canonical 5).  Sets sigma loading.
    F
        Koide numerator integer (canonical 2).
    R
        Koide denominator integer (canonical 3).  F*R = 6 sets pi loading.
    sp2_bond_count
        Number of sp^2 in-plane bonds per boron atom (3, fixed by sp^2
        hybridisation).  Numerator of the sigma loading 3 / K_rank.
    pz_axis_count
        Number of perpendicular p_z axes (1, fixed by sp^2 hybridisation).
        Numerator of the pi loading 1 / (F * R).
    Tc_K
        MgB_2 critical temperature for the BCS ratio computation
        (39.0 K experimental anchor).
    """

    Lambda_QCD_MeV:  float = LAMBDA_QCD_MEV
    n_R:             int   = n_R
    K_rank:          int   = K_rank
    F:               int   = F
    R:               int   = R
    sp2_bond_count:  int   = 3       # sp^2 in-plane bonds per boron
    pz_axis_count:   int   = 1       # perpendicular p_z axes per boron
    Tc_K:            float = MGB2_TC_K_REF


# ---------------------------------------------------------------------------
# Core substrate derivations (no fit parameters)
# ---------------------------------------------------------------------------

def eps_unit_meV(bands: MgB2SubstrateBands = MgB2SubstrateBands()) -> float:
    """Substrate band-gap unit: eps_unit = Lambda_QCD / n_R.

    For canonical (Lambda_QCD = 200 MeV, n_R = 18) this is

        eps_unit = 200 / 18 meV = 11.111 meV

    Identical to the substrate-saturation energy that gives the
    high-T_c ceiling T_c,max = eps_unit / k_B = 128.9 K
    (see ``b3_high_tc_bound``).
    """
    return float(bands.Lambda_QCD_MeV / bands.n_R)


def gap_sigma_band(bands: MgB2SubstrateBands = MgB2SubstrateBands()) -> float:
    """Substrate prediction for Delta_sigma (in-plane sp^2 / K_4 face).

    Derivation
    ----------
        Delta_sigma = eps_unit * (sp2_bond_count / K_rank)
                    = (Lambda_QCD / n_R) * (3 / K_rank)
                    = 200/18 * 3/5  meV
                    = 6.667 meV

    Substrate cell reading
    ----------------------
    The K_4 tetrahedron face seen edge-on into the honeycomb plane is
    a triangle of 3 sp^2 bonds (each boron has 3 in-plane bonds to its
    3 honeycomb neighbours).  These 3 substrate-mode bonds occupy 3 of
    the K_rank = 5 vertices of the 4-simplex spanning one cell-pair
    bridge.  Saturation fraction = 3 / K_rank.

    Returns
    -------
    Delta_sigma in milli-electronvolts (HALF gap, not full 2 Delta).
    """
    return float(eps_unit_meV(bands) * (bands.sp2_bond_count / bands.K_rank))


def gap_pi_band(bands: MgB2SubstrateBands = MgB2SubstrateBands()) -> float:
    """Substrate prediction for Delta_pi (out-of-plane p_z / Q_3 axis).

    Derivation
    ----------
        Delta_pi = eps_unit * (pz_axis_count / (F * R))
                 = (Lambda_QCD / n_R) * (1 / 6)
                 = 200/18 * 1/6  meV
                 = 1.852 meV

    Substrate cell reading
    ----------------------
    The 3-cube Q_3 has 3 axes; the perpendicular p_z axis is one of
    them (pz_axis_count = 1).  The Koide F/R = 2/3 mechanism suppresses
    this single axis by F * R = 6 (= 2 Moebius sheets x 3 lepton-rank
    denominator) -- it is the same suppression that splits the
    charged-lepton mass tower.  Saturation fraction = 1 / (F * R).

    Returns
    -------
    Delta_pi in milli-electronvolts (HALF gap, not full 2 Delta).
    """
    return float(eps_unit_meV(bands) * (bands.pz_axis_count / (bands.F * bands.R)))


def ratio_sigma_pi(bands: MgB2SubstrateBands = MgB2SubstrateBands()) -> float:
    """Closed-form Delta_sigma / Delta_pi.

    Algebraic identity at canonical integers:

        Delta_sigma / Delta_pi
            = (sp2_bond_count / K_rank) / (pz_axis_count / (F * R))
            = sp2_bond_count * F * R / (K_rank * pz_axis_count)
            = 3 * 2 * 3 / (5 * 1)
            = 18 / 5
            = 3.6

    vs empirical 7 meV / 2 meV = 3.5 (2.9 % under).
    """
    return float(
        (bands.sp2_bond_count * bands.F * bands.R)
        / (bands.K_rank * bands.pz_axis_count)
    )


def bcs_band_ratio(delta_meV: float, T_c_K: float) -> float:
    """Per-band 2 Delta / (k_B T_c) ratio.

    ``delta_meV`` is the half gap in meV; the function returns the
    full-gap ratio 2 Delta / (k_B T_c).
    """
    if T_c_K <= 0.0:
        raise ValueError(f"T_c must be > 0, got {T_c_K!r}")
    if delta_meV <= 0.0:
        raise ValueError(f"delta must be > 0, got {delta_meV!r}")
    return float(2.0 * delta_meV / (K_B_MEV_K * T_c_K))


# ---------------------------------------------------------------------------
# Bundle predictions
# ---------------------------------------------------------------------------

def predict_mgb2_two_gap(
    bands: MgB2SubstrateBands = MgB2SubstrateBands(),
) -> Dict[str, float]:
    """Substrate two-gap prediction bundle for MgB_2.

    Keys
    ----
    Delta_sigma_MeV   : substrate-derived sigma-band half gap
    Delta_pi_MeV      : substrate-derived pi-band   half gap
    eps_unit_MeV      : Lambda_QCD / n_R band-gap unit
    R_sigma_pred      : 2 Delta_sigma / (k_B T_c)
    R_pi_pred         : 2 Delta_pi    / (k_B T_c)
    R_sigma_dev_pct   : (R_sigma_pred - 4.0) / 4.0 * 100
    R_pi_dev_pct      : (R_pi_pred    - 1.2) / 1.2 * 100
    Delta_sigma_dev_pct
                       : (Delta_sigma_pred - 7.0) / 7.0 * 100
    Delta_pi_dev_pct  : (Delta_pi_pred    - 2.0) / 2.0 * 100
    ratio_sigma_pi    : Delta_sigma / Delta_pi  (substrate)
    ratio_sigma_pi_exp: 7.0 / 2.0  (empirical)
    ratio_sigma_pi_dev_pct
                       : (ratio_pred - ratio_exp) / ratio_exp * 100
    """
    Delta_sigma = gap_sigma_band(bands)
    Delta_pi    = gap_pi_band(bands)
    R_sigma = bcs_band_ratio(Delta_sigma, bands.Tc_K)
    R_pi    = bcs_band_ratio(Delta_pi,    bands.Tc_K)
    ratio = ratio_sigma_pi(bands)
    return {
        "Delta_sigma_MeV":     Delta_sigma,
        "Delta_pi_MeV":        Delta_pi,
        "eps_unit_MeV":        eps_unit_meV(bands),
        "R_sigma_pred":        R_sigma,
        "R_pi_pred":           R_pi,
        "R_sigma_dev_pct":
            (R_sigma - MGB2_R_SIGMA_REF) / MGB2_R_SIGMA_REF * 100.0,
        "R_pi_dev_pct":
            (R_pi    - MGB2_R_PI_REF)    / MGB2_R_PI_REF    * 100.0,
        "Delta_sigma_dev_pct":
            (Delta_sigma - MGB2_DELTA_SIGMA_MEV_EXP)
            / MGB2_DELTA_SIGMA_MEV_EXP * 100.0,
        "Delta_pi_dev_pct":
            (Delta_pi    - MGB2_DELTA_PI_MEV_EXP)
            / MGB2_DELTA_PI_MEV_EXP    * 100.0,
        "ratio_sigma_pi":       ratio,
        "ratio_sigma_pi_exp":   MGB2_DELTA_SIGMA_MEV_EXP / MGB2_DELTA_PI_MEV_EXP,
        "ratio_sigma_pi_dev_pct":
            (ratio - (MGB2_DELTA_SIGMA_MEV_EXP / MGB2_DELTA_PI_MEV_EXP))
            / (MGB2_DELTA_SIGMA_MEV_EXP / MGB2_DELTA_PI_MEV_EXP) * 100.0,
    }


def substrate_multiband_summary(
    bands: MgB2SubstrateBands = MgB2SubstrateBands(),
) -> Dict[str, object]:
    """One-shot summary of the substrate multiband derivation.

    Returns a dict mixing scalars, the underlying integers, and the
    full prediction bundle, suitable for printing or logging.
    """
    pred = predict_mgb2_two_gap(bands)
    return {
        "anchors_and_integers": {
            "Lambda_QCD_MeV": bands.Lambda_QCD_MeV,
            "n_R":            bands.n_R,
            "K_rank":         bands.K_rank,
            "F":              bands.F,
            "R":              bands.R,
            "sp2_bond_count": bands.sp2_bond_count,
            "pz_axis_count":  bands.pz_axis_count,
            "Tc_K":           bands.Tc_K,
        },
        "experimental_anchors": {
            "Delta_sigma_MeV_exp": MGB2_DELTA_SIGMA_MEV_EXP,
            "Delta_pi_MeV_exp":    MGB2_DELTA_PI_MEV_EXP,
            "R_sigma_ref_Choi":    MGB2_R_SIGMA_REF,
            "R_pi_ref_Choi":       MGB2_R_PI_REF,
        },
        "prediction": pred,
        "headline_close_to_measurement_pct": {
            "Delta_sigma": pred["Delta_sigma_dev_pct"],
            "Delta_pi":    pred["Delta_pi_dev_pct"],
            "R_sigma":     pred["R_sigma_dev_pct"],
            "R_pi":        pred["R_pi_dev_pct"],
            "ratio":       pred["ratio_sigma_pi_dev_pct"],
        },
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    summary = substrate_multiband_summary()
    print("Substrate multiband MgB_2 derivation")
    print("=" * 64)
    print()
    print("Substrate integers (b3_constants):")
    for k, v in summary["anchors_and_integers"].items():
        print(f"  {k:>15s} = {v}")
    print()
    print("Experimental anchors (Choi 2002 + ARPES/STM):")
    for k, v in summary["experimental_anchors"].items():
        print(f"  {k:>22s} = {v}")
    print()
    pred = summary["prediction"]
    print(f"eps_unit = Lambda_QCD / n_R = {pred['eps_unit_MeV']:.4f} meV")
    print()
    print(f"Delta_sigma  pred = {pred['Delta_sigma_MeV']:.4f} meV  "
          f"(measured 7.0 meV, dev {pred['Delta_sigma_dev_pct']:+.2f}%)")
    print(f"Delta_pi     pred = {pred['Delta_pi_MeV']:.4f} meV  "
          f"(measured 2.0 meV, dev {pred['Delta_pi_dev_pct']:+.2f}%)")
    print(f"ratio sigma/pi    = {pred['ratio_sigma_pi']:.4f}      "
          f"(measured 3.5,    dev {pred['ratio_sigma_pi_dev_pct']:+.2f}%)")
    print()
    print(f"R_sigma = 2 Delta_sigma / (k_B T_c) "
          f"= {pred['R_sigma_pred']:.4f}  "
          f"(Choi {MGB2_R_SIGMA_REF}, dev {pred['R_sigma_dev_pct']:+.2f}%)")
    print(f"R_pi    = 2 Delta_pi    / (k_B T_c) "
          f"= {pred['R_pi_pred']:.4f}  "
          f"(Choi {MGB2_R_PI_REF}, dev {pred['R_pi_dev_pct']:+.2f}%)")


__all__ = [
    "K_B_MEV_K",
    "MGB2_DELTA_PI_MEV_EXP",
    "MGB2_DELTA_SIGMA_MEV_EXP",
    "MGB2_R_PI_REF",
    "MGB2_R_SIGMA_REF",
    "MGB2_TC_K_REF",
    "MgB2SubstrateBands",
    "bcs_band_ratio",
    "eps_unit_meV",
    "gap_pi_band",
    "gap_sigma_band",
    "main",
    "predict_mgb2_two_gap",
    "ratio_sigma_pi",
    "substrate_multiband_summary",
]


if __name__ == "__main__":  # pragma: no cover
    main()
