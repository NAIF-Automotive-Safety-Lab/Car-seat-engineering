"""BCS universal gap ratio: substrate prediction vs measured superconductors.

The substrate-paired Cooper-bridge ground state predicts a UNIVERSAL,
parameter-free weak-coupling ratio between the zero-temperature gap and
the critical temperature:

       2 Δ(0) / (k_B T_c)  =  2 π / e^γ   ≃   3.527754

where γ = 0.5772156649... is the Euler-Mascheroni constant.  The same
result is the BCS-1957 universal ratio, derived in the substrate
ontology from the paired-strain bridge unbinding condition (see
``superconductivity_substrate.SuperconductivitySimulator.bcs_ratio``).

This module compares that prediction against tabulated literature
values of (T_c, 2Δ(0)) for ten elemental and one multiband
superconductor.

Strong-coupling (Allen--Dynes) extension
----------------------------------------
For materials where lambda * (T_c / omega_log)^2 is non-negligible the
substrate / weak-coupling line is corrected by the Allen--Dynes
strong-coupling envelope (Allen & Dynes 1975; Carbotte 1990 review):

    2 Delta(0) / (k_B T_c)
        ~ (2 pi / e^gamma) * [ 1 + 12.5 * (T_c / omega_log)^2
                                       * ln(omega_log / (2 T_c)) ]

where omega_log is the logarithmic average of the phonon spectrum
alpha^2 F(omega) (in K).  The substrate ontology accommodates this
correction WITHOUT adding new parameters: lambda and omega_log are
material-specific phonon-spectrum inputs, the same inputs Eliashberg
theory takes (and the same inputs the substrate framework would have
to take to specify which medium-mode spectrum is being unbound).
The correction does NOT change the substrate-derived universal ratio
2 pi / e^gamma -- it adds the leading T_c/omega_log expansion term that
both Eliashberg and substrate strong-coupling treatments produce.

For multiband superconductors (MgB_2) the single-band Allen--Dynes form
fails: each band has its own gap and its own coupling.  The module
provides a separate ``mgb2_multiband_prediction`` that handles the
sigma- and pi-channel gaps independently, with an optional DOS-weighted
composite gap for visual comparison.

Materials
---------
The reference table here is independently sourced (Tinkham, Carbotte,
NIST, PDG element-card values) and intentionally distinct from the
back-solved Hg/HBCCO calibration table in
``superconductivity_substrate.REFERENCE_SUPERCONDUCTORS`` --- this is a
*test*, not a refit:

    name     T_c (K)     2Δ(0) (meV)     class
    -------  ----------  --------------  -----------------------------
    Hg        4.15        1.41            elemental, weak-coupling BCS
    Pb        7.20        2.74            elemental, strong-coupling BCS
    Sn        3.72        1.15            elemental, weak-coupling BCS
    Al        1.20        0.34            elemental, weak-coupling BCS
    Nb        9.25        3.05            elemental, weak/strong-borderline BCS
    V         5.40        1.55            elemental, weak-coupling BCS
    Ta        4.48        1.40            elemental, weak-coupling BCS
    In        3.40        1.05            elemental, weak-coupling BCS
    Tl        2.39        0.74            elemental, weak-coupling BCS
    MgB_2    39.0        14.0             multiband (σ + π), expected to deviate

The measured ratio R_meas = 2Δ(0) / (k_B T_c) is compared to the
substrate prediction R_pred = 2π / e^γ ≃ 3.528.  Per-material percent
deviation is reported.

API
---
``MATERIALS``                  : tuple of ``MaterialDatum`` (immutable).
``BCS_RATIO_PRED``             : float, 2π/e^γ.
``measured_ratio(T_c, gap)``   : compute R_meas from K and meV.
``deviation_percent(R)``       : (R - R_pred) / R_pred * 100.
``predict_with_allen_dynes(lam, omega_log_K, T_c_K)``
                               : strong-coupling-corrected gap ratio.
``MATERIAL_PHONON_PARAMS``     : dict of (lambda, omega_log_K) per material
                                 from McMillan / Allen-Dynes literature.
``mgb2_multiband_prediction()``: per-band MgB_2 sigma/pi ratios + DOS-
                                 weighted composite.
``run_test()``                 : returns a list of ``ResultRow`` with
                                 (name, T_c, gap_meV, R_meas, dev_pct,
                                 dev_pct_ad, within_5pct).
``summary_stats(rows)``        : mean / max / min |dev|, n_within_5pct
                                 (bare AND Allen-Dynes-corrected).
``render_bcs_gap_ratio_test()``: matplotlib bar-chart figure for the
                                 visuals/120 PNG (bare BCS / AD / measured).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Tuple

import math

import numpy as np


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

K_B:        float = 1.380649e-23           # J/K   Boltzmann
EV_PER_J:   float = 1.0 / 1.602176634e-19  # eV / J
K_B_MEV_K:  float = (K_B * EV_PER_J) * 1e3 # 1 K -> meV  (= 0.0861733...)


# ---------------------------------------------------------------------------
# Substrate / BCS prediction
# ---------------------------------------------------------------------------

EULER_GAMMA:    float = 0.5772156649015329
BCS_RATIO_PRED: float = 2.0 * math.pi / math.exp(EULER_GAMMA)
# = 3.527753977724091


# ---------------------------------------------------------------------------
# Reference data (literature)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MaterialDatum:
    """One measured superconductor entry.

    Attributes
    ----------
    name      : symbol / label.
    T_c_K     : critical temperature (K), measured at zero applied field
                and ambient pressure.
    gap_meV   : full gap 2 Δ(0) in milli-electronvolts at T -> 0.
    klass     : qualitative class -- ``"elemental"``, ``"multiband"``.
    note      : short physics note (coupling regime / known anomalies).
    """
    name:    str
    T_c_K:   float
    gap_meV: float
    klass:   str
    note:    str


MATERIALS: Tuple[MaterialDatum, ...] = (
    MaterialDatum("Hg",   4.15,  1.41, "elemental",
                  "Onnes 1911; mildly strong-coupling"),
    MaterialDatum("Pb",   7.20,  2.74, "elemental",
                  "Strong-coupling Eliashberg, ratio ~4.4 known"),
    MaterialDatum("Sn",   3.72,  1.15, "elemental",
                  "Weak-coupling, near-textbook BCS"),
    MaterialDatum("Al",   1.20,  0.34, "elemental",
                  "Weak-coupling reference, ratio ~3.3-3.5"),
    MaterialDatum("Nb",   9.25,  3.05, "elemental",
                  "Weak-/strong-coupling borderline; 2-band hints"),
    MaterialDatum("V",    5.40,  1.55, "elemental",
                  "Weak-coupling, near textbook"),
    MaterialDatum("Ta",   4.48,  1.40, "elemental",
                  "Weak-coupling, near textbook"),
    MaterialDatum("In",   3.40,  1.05, "elemental",
                  "Weak-coupling, near textbook"),
    MaterialDatum("Tl",   2.39,  0.74, "elemental",
                  "Weak-coupling, near textbook"),
    MaterialDatum("MgB2", 39.0, 14.00, "multiband",
                  "Two-band sigma+pi: dominant sigma gap, deviation expected"),
)


# ---------------------------------------------------------------------------
# Core arithmetic
# ---------------------------------------------------------------------------

def measured_ratio(T_c_K: float, gap_meV: float) -> float:
    """Return R_meas = 2Δ(0) / (k_B T_c), DIMENSIONLESS.

    ``gap_meV`` is the full gap 2Δ(0) in meV (the literature convention
    used throughout this module's table).  ``T_c_K`` is in Kelvin.

    R_meas = (gap_meV * 1e-3 / EV_PER_J) / (K_B * T_c_K)
           = gap_meV / (K_B_MEV_K * T_c_K).
    """
    if T_c_K <= 0.0:
        raise ValueError(f"T_c must be > 0, got {T_c_K!r}")
    if gap_meV <= 0.0:
        raise ValueError(f"gap must be > 0, got {gap_meV!r}")
    return float(gap_meV / (K_B_MEV_K * T_c_K))


def deviation_percent(R_meas: float, R_pred: float = BCS_RATIO_PRED) -> float:
    """Return (R_meas - R_pred) / R_pred * 100."""
    return float((R_meas - R_pred) / R_pred * 100.0)


# ---------------------------------------------------------------------------
# Allen-Dynes strong-coupling correction
# ---------------------------------------------------------------------------
#
# Allen & Dynes 1975, Phys. Rev. B 12, 905; Carbotte 1990, Rev. Mod. Phys.
# 62, 1027 (table I, equation 2.27).  The leading correction beyond the
# universal weak-coupling ratio is
#
#     2 Delta(0) / (k_B T_c) = (2 pi / e^gamma)
#         * [ 1 + 12.5 * (T_c / omega_log)^2 * ln(omega_log / (2 T_c)) ]
#
# valid for T_c << omega_log.  For T_c approaching omega_log the
# correction grows large and the truncated series is no longer a good
# approximation; callers are expected to use it only where the sign of
# the bracket is positive (i.e. omega_log > 2 T_c).

ALLEN_DYNES_PREFACTOR: float = 12.5


def predict_with_allen_dynes(
    lam:         float,
    omega_log_K: float,
    T_c_K:       float,
    R_pred_weak: float = BCS_RATIO_PRED,
) -> float:
    """Strong-coupling-corrected gap ratio per Allen--Dynes.

    Parameters
    ----------
    lam
        Electron-phonon coupling strength (mass enhancement parameter
        lambda).  Listed here for API completeness and future extension
        but does NOT enter the leading T_c/omega_log correction term;
        the lambda dependence is folded into T_c itself via McMillan /
        Allen-Dynes T_c formula.  Higher-order corrections in lambda
        (mu* renormalisation, strong-coupling shape factors) are NOT
        included in the leading-term form used here.
    omega_log_K
        Phonon log-average frequency in K (= hbar omega_log / k_B).
        Tabulated values come from McMillan 1968 and Allen--Dynes 1975
        re-fits of measured tunnelling alpha^2 F(omega) spectra.
    T_c_K
        Critical temperature in K.
    R_pred_weak
        Weak-coupling baseline (default = substrate / BCS = 2 pi/e^gamma).

    Returns
    -------
    float
        Predicted 2 Delta(0) / (k_B T_c).  Returns ``R_pred_weak``
        unchanged if omega_log <= 2 T_c (the leading-term correction
        breaks down outside that regime; signal that AD truncation is
        not applicable rather than report a negative correction).
    """
    if T_c_K <= 0.0:
        raise ValueError(f"T_c must be > 0, got {T_c_K!r}")
    if omega_log_K <= 0.0:
        raise ValueError(f"omega_log must be > 0, got {omega_log_K!r}")
    if lam < 0.0:
        raise ValueError(f"lambda must be >= 0, got {lam!r}")
    if omega_log_K <= 2.0 * T_c_K:
        # Leading-order Allen-Dynes is not valid here; do not extrapolate.
        return float(R_pred_weak)
    correction = (
        ALLEN_DYNES_PREFACTOR
        * (T_c_K / omega_log_K) ** 2
        * math.log(omega_log_K / (2.0 * T_c_K))
    )
    return float(R_pred_weak * (1.0 + correction))


# Material-specific phonon parameters from McMillan 1968 / Allen-Dynes 1975
# tabulations (Carbotte 1990 review, table I) and Choi et al. 2002 for MgB_2.
#
#     name : (lambda, omega_log [K])
#
# These are EMPIRICAL inputs to the strong-coupling correction; they are
# not derived from substrate axioms.  The substrate framework predicts
# the universal weak-coupling line 2 pi / e^gamma; the Allen-Dynes
# correction pulls in the same external phonon-spectrum data that
# Eliashberg theory takes.

MATERIAL_PHONON_PARAMS: dict = {
    # name: (lambda, omega_log_K)  -- Allen-Dynes / Carbotte literature
    "Hg":   (1.62,  72.0),
    "Pb":   (1.55,  56.0),
    "Sn":   (0.72, 152.0),
    "Al":   (0.43, 305.0),
    "Nb":   (0.82, 130.0),
    "V":    (0.63, 217.0),
    "Ta":   (0.69, 159.0),
    "In":   (0.81, 100.0),
    "Tl":   (0.71,  72.0),
    # MgB_2 dominant sigma channel (Choi et al. 2002, Nature 418, 758).
    # The single-band Allen-Dynes form is intentionally a poor fit here;
    # ``mgb2_multiband_prediction`` provides the per-band breakdown.
    "MgB2": (0.96, 778.0),
}


# ---------------------------------------------------------------------------
# Substrate-derived (lambda, omega_log) lookup
# ---------------------------------------------------------------------------
#
# The substrate Eliashberg derivation in
# :mod:`src.stiff_medium.substrate_eliashberg` builds alpha^2 F(omega) from
# the substrate phonon spectrum (Debye-quadratic acoustic background +
# Einstein peak at (2/3) omega_D from K_4-cell internal modes) plus a
# Hopfield-style electron-phonon coupling whose magnitude is anchored ONCE
# on Pb (lambda_Pb = 1.55).  The function below builds a substitute for
# MATERIAL_PHONON_PARAMS using ONLY substrate-derivable inputs (M_ion,
# elastic moduli, valence) for the SP-elemental subset; transition metals
# (Nb, V, Ta) are passed through with the empirical entries because the
# free-electron substrate N(0) does not capture their d-band Fermi-level
# DOS (a known substrate-Eliashberg honest-caveat).
#
# This lets ``run_test(phonon_params=substrate_phonon_params())`` rerun
# the full BCS gap-ratio test with substrate-derived strong-coupling
# corrections, providing a head-to-head comparison with the empirical
# baseline.

def substrate_phonon_params(
    fall_back_to_empirical: bool = True,
) -> dict:
    """Return a (lambda, omega_log_K) dict derived from the substrate
    Eliashberg derivation in :mod:`substrate_eliashberg`.

    Parameters
    ----------
    fall_back_to_empirical
        If True (default), names not in
        ``substrate_eliashberg.MATERIAL_DB`` are filled in from
        ``MATERIAL_PHONON_PARAMS`` so the returned dict has the same key
        set as ``MATERIAL_PHONON_PARAMS`` (typically the only fall-back is
        MgB_2, which the substrate single-band model does not handle).
        If False, the returned dict only contains substrate-supplied keys.

    Returns
    -------
    dict[str, tuple[float, float]]
        Same shape as MATERIAL_PHONON_PARAMS but with substrate-derived
        (lambda, omega_log_K) for SP-metal entries.
    """
    # Lazy import to avoid a circular import at module load time.
    from .substrate_eliashberg import MATERIAL_DB, predict_lambda_omega_log
    out: dict = {}
    for name in MATERIAL_DB:
        out[name] = predict_lambda_omega_log(name)
    if fall_back_to_empirical:
        for name, val in MATERIAL_PHONON_PARAMS.items():
            if name not in out:
                out[name] = val
    return out


# ---------------------------------------------------------------------------
# Multiband module: MgB_2 sigma + pi
# ---------------------------------------------------------------------------
#
# MgB_2 is the canonical two-band superconductor.  The sigma band hosts
# the dominant gap (Delta_sigma ~ 7 meV, full gap 2 Delta_sigma ~ 14 meV)
# and a moderate coupling (lambda_sigma ~ 1.0); the pi band hosts a
# secondary gap (Delta_pi ~ 2 meV, 2 Delta_pi ~ 4 meV) with weak
# coupling (lambda_pi ~ 0.4).  Choi et al. 2002 (Nature 418, 758)
# reports per-band gap ratios 2 Delta_sigma / k_B T_c ~ 4.0 and
# 2 Delta_pi / k_B T_c ~ 1.2 from two-band Eliashberg.  The single-band
# universal ratio simply does not apply to either band individually.

@dataclass(frozen=True)
class MgB2Bands:
    """Per-band parameters for the MgB_2 two-gap fit."""
    Tc_K:           float = 39.0
    Delta_sigma_meV: float = 7.0    # half-gap on sigma band (Choi 2002)
    Delta_pi_meV:    float = 2.0    # half-gap on pi band   (Choi 2002)
    w_sigma:         float = 0.55   # DOS-weighted sigma fraction
    lambda_sigma:    float = 0.96
    omega_log_sigma_K: float = 778.0
    lambda_pi:       float = 0.42
    omega_log_pi_K:  float = 696.0


def mgb2_multiband_prediction(bands: MgB2Bands = MgB2Bands()) -> dict:
    """Per-band gap-ratio breakdown for MgB_2.

    The MEASURED-side rows (R_sigma, R_pi, R_composite) use the Choi 2002
    half-gaps stored on the bands dataclass.  The SUBSTRATE-side rows
    (R_sigma_substrate_pred, R_pi_substrate_pred, etc.) are derived
    parameter-free from the K_4 face-pair (sigma) and Q_3 axis (pi)
    couplings encoded in :mod:`src.stiff_medium.substrate_multiband`.

    Returns a dict with:
        R_sigma            : 2 Delta_sigma_exp / (k_B T_c)         (measured)
        R_pi               : 2 Delta_pi_exp    / (k_B T_c)         (measured)
        R_composite        : (w_sigma * 2 Delta_sigma + w_pi * 2 Delta_pi)
                              / (k_B T_c)                          (measured)
        R_weighted_substrate
                          : DOS-weighted average of measured band ratios
                            (algebraically identical to R_composite for
                            a fixed T_c).
        R_pred_weak        : substrate / BCS weak-coupling line (3.528)
        R_pred_dominant_AD : Allen-Dynes correction applied to dominant
                              sigma band (with sigma lambda, omega_log)
        dev_composite_pct  : composite vs measured (using full sigma 2D=14)
        Choi_R_sigma_ref   : 4.0 (Choi et al. 2002 Eliashberg result)
        Choi_R_pi_ref      : 1.2 (Choi et al. 2002 Eliashberg result)

        --- substrate (K_4 / Q_3) parameter-free predictions ---
        Delta_sigma_substrate_pred_meV
                          : (Lambda_QCD / n_R) * (3 / K_rank)  =  6.667 meV
        Delta_pi_substrate_pred_meV
                          : (Lambda_QCD / n_R) * (1 / (F * R)) =  1.852 meV
        R_sigma_substrate_pred
                          : 2 Delta_sigma_substrate / (k_B T_c)
        R_pi_substrate_pred
                          : 2 Delta_pi_substrate    / (k_B T_c)
        dev_sigma_substrate_pct
                          : substrate sigma vs Choi 4.0
        dev_pi_substrate_pct
                          : substrate pi    vs Choi 1.2
        ratio_sigma_pi_substrate
                          : 3 * F * R / K_rank = 18/5 = 3.6
                            (vs 7/2 = 3.5 measured)
    """
    # Measured side -------------------------------------------------
    twoD_sigma = 2.0 * bands.Delta_sigma_meV
    twoD_pi    = 2.0 * bands.Delta_pi_meV
    R_sigma = measured_ratio(bands.Tc_K, twoD_sigma)
    R_pi    = measured_ratio(bands.Tc_K, twoD_pi)
    avg_gap = bands.w_sigma * twoD_sigma + (1.0 - bands.w_sigma) * twoD_pi
    R_composite = measured_ratio(bands.Tc_K, avg_gap)
    R_pred_dominant = predict_with_allen_dynes(
        bands.lambda_sigma, bands.omega_log_sigma_K, bands.Tc_K
    )

    # Substrate K_4 / Q_3 predictions (parameter-free) ----------------
    # Local import to keep the bcs_gap_ratio_test module importable even
    # when the substrate_multiband module is absent (graceful fallback).
    try:
        from .substrate_multiband import (
            MgB2SubstrateBands,
            gap_pi_band,
            gap_sigma_band,
            ratio_sigma_pi,
        )
        sub_bands = MgB2SubstrateBands(Tc_K=bands.Tc_K)
        Delta_sigma_sub = gap_sigma_band(sub_bands)
        Delta_pi_sub    = gap_pi_band(sub_bands)
        R_sigma_sub = measured_ratio(bands.Tc_K, 2.0 * Delta_sigma_sub)
        R_pi_sub    = measured_ratio(bands.Tc_K, 2.0 * Delta_pi_sub)
        ratio_sub   = ratio_sigma_pi(sub_bands)
        dev_sigma_substrate_pct = deviation_percent(
            R_sigma_sub, R_pred=4.0)
        dev_pi_substrate_pct = deviation_percent(
            R_pi_sub, R_pred=1.2)
    except ImportError:
        Delta_sigma_sub = float("nan")
        Delta_pi_sub    = float("nan")
        R_sigma_sub     = float("nan")
        R_pi_sub        = float("nan")
        ratio_sub       = float("nan")
        dev_sigma_substrate_pct = float("nan")
        dev_pi_substrate_pct    = float("nan")

    return {
        "R_sigma":            R_sigma,
        "R_pi":               R_pi,
        "R_composite":        R_composite,
        "R_weighted_substrate":
            bands.w_sigma * R_sigma + (1.0 - bands.w_sigma) * R_pi,
        "R_pred_weak":        BCS_RATIO_PRED,
        "R_pred_dominant_AD": R_pred_dominant,
        "dev_sigma_pct":      deviation_percent(R_sigma, R_pred_dominant),
        "dev_composite_pct":  deviation_percent(R_composite, BCS_RATIO_PRED),
        "Choi_R_sigma_ref":   4.0,
        "Choi_R_pi_ref":      1.2,
        # Substrate (K_4 / Q_3) zero-parameter predictions
        "Delta_sigma_substrate_pred_meV": Delta_sigma_sub,
        "Delta_pi_substrate_pred_meV":    Delta_pi_sub,
        "R_sigma_substrate_pred":         R_sigma_sub,
        "R_pi_substrate_pred":            R_pi_sub,
        "dev_sigma_substrate_pct":        dev_sigma_substrate_pct,
        "dev_pi_substrate_pct":           dev_pi_substrate_pct,
        "ratio_sigma_pi_substrate":       ratio_sub,
    }


# ---------------------------------------------------------------------------
# Test driver
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ResultRow:
    """Per-material BCS gap ratio test result.

    Attributes
    ----------
    name, T_c_K, gap_meV, klass, note
        Pass-through from ``MaterialDatum``.
    R_meas
        Measured 2 Delta(0) / (k_B T_c).
    R_pred
        Bare substrate / BCS weak-coupling prediction (3.528).
    R_pred_ad
        Allen--Dynes-corrected prediction using literature
        (lambda, omega_log) for this material.  Equals ``R_pred``
        when no phonon parameters are available or the AD correction
        is out of its validity range.
    R_pred_ad_substrate
        Allen--Dynes-corrected prediction using SUBSTRATE-DERIVED
        (lambda, omega_log) from :mod:`substrate_eliashberg`.  Equals
        ``R_pred`` when no substrate-derived phonon entry exists for
        the material.
    dev_pct
        Deviation of R_meas from bare R_pred, in percent.
    dev_pct_ad
        Deviation of R_meas from Allen-Dynes-corrected R_pred_ad
        (empirical phonon parameters).
    dev_pct_ad_substrate
        Deviation of R_meas from substrate-AD R_pred_ad_substrate.
    within_5pct
        ``abs(dev_pct) <= threshold_pct`` (BARE).
    within_5pct_ad
        ``abs(dev_pct_ad) <= threshold_pct`` (Allen-Dynes-corrected,
        EMPIRICAL phonon parameters).
    within_5pct_ad_substrate
        ``abs(dev_pct_ad_substrate) <= threshold_pct`` (Allen-Dynes-
        corrected, SUBSTRATE-derived phonon parameters).
    """
    name:           str
    T_c_K:          float
    gap_meV:        float
    klass:          str
    R_meas:         float
    R_pred:         float
    R_pred_ad:      float
    R_pred_ad_substrate: float
    dev_pct:        float
    dev_pct_ad:     float
    dev_pct_ad_substrate: float
    within_5pct:    bool
    within_5pct_ad: bool
    within_5pct_ad_substrate: bool
    note:           str


def run_test(
    materials: Iterable[MaterialDatum] = MATERIALS,
    R_pred:    float = BCS_RATIO_PRED,
    threshold_pct: float = 5.0,
    phonon_params: dict = MATERIAL_PHONON_PARAMS,
    phonon_params_substrate: dict | None = None,
) -> List[ResultRow]:
    """Compute R_meas, bare and Allen-Dynes-corrected deviations,
    and pass/fail at ``threshold_pct``.

    Returns one ``ResultRow`` per input material.  Order is preserved.

    For multiband materials (klass == "multiband") the Allen-Dynes
    correction is applied to the DOMINANT band (e.g. MgB_2 sigma) using
    the listed (lambda, omega_log).  The single-band AD line is known
    to underestimate the multiband gap ratio; see
    ``mgb2_multiband_prediction`` for the per-band breakdown.

    Parameters
    ----------
    phonon_params
        EMPIRICAL (lambda, omega_log) per material from Allen-Dynes /
        Carbotte literature.  Default ``MATERIAL_PHONON_PARAMS``.
    phonon_params_substrate
        SUBSTRATE-DERIVED (lambda, omega_log) per material from
        :func:`substrate_phonon_params`.  When ``None`` (default), the
        function builds the substrate-derived map automatically; pass
        an empty dict to disable substrate-AD computation entirely.
    """
    if phonon_params_substrate is None:
        try:
            phonon_params_substrate = substrate_phonon_params()
        except (ImportError, KeyError):
            # Graceful fallback: substrate module unavailable.
            phonon_params_substrate = {}

    rows: List[ResultRow] = []
    for m in materials:
        R = measured_ratio(m.T_c_K, m.gap_meV)
        d = deviation_percent(R, R_pred)
        if m.name in phonon_params:
            lam, wlog = phonon_params[m.name]
            R_pred_ad = predict_with_allen_dynes(lam, wlog, m.T_c_K, R_pred)
        else:
            R_pred_ad = R_pred
        if m.name in phonon_params_substrate:
            lam_s, wlog_s = phonon_params_substrate[m.name]
            R_pred_ad_s = predict_with_allen_dynes(
                lam_s, wlog_s, m.T_c_K, R_pred)
        else:
            R_pred_ad_s = R_pred
        d_ad   = deviation_percent(R, R_pred_ad)
        d_ad_s = deviation_percent(R, R_pred_ad_s)
        rows.append(ResultRow(
            name=m.name,
            T_c_K=m.T_c_K,
            gap_meV=m.gap_meV,
            klass=m.klass,
            R_meas=R,
            R_pred=R_pred,
            R_pred_ad=R_pred_ad,
            R_pred_ad_substrate=R_pred_ad_s,
            dev_pct=d,
            dev_pct_ad=d_ad,
            dev_pct_ad_substrate=d_ad_s,
            within_5pct=abs(d) <= threshold_pct,
            within_5pct_ad=abs(d_ad) <= threshold_pct,
            within_5pct_ad_substrate=abs(d_ad_s) <= threshold_pct,
            note=m.note,
        ))
    return rows


def summary_stats(rows: List[ResultRow]) -> dict:
    """Aggregate statistics across a result set.

    Returns a dict with:
        n                 : number of materials
        n_within_5pct     : count of |deviation| <= 5 % (BARE)
        n_within_5pct_ad  : count of |deviation| <= 5 % (Allen-Dynes-corrected)
        mean_abs_dev      : mean |deviation| in percent (BARE)
        mean_abs_dev_ad   : mean |deviation| (Allen-Dynes-corrected)
        max_abs_dev       : max  |deviation| (BARE)
        max_abs_dev_ad    : max  |deviation| (Allen-Dynes-corrected)
        min_abs_dev       : min  |deviation| (BARE)
        min_abs_dev_ad    : min  |deviation| (Allen-Dynes-corrected)
        elemental_only    : restricted to klass=="elemental" (BARE + AD)
    """
    if not rows:
        raise ValueError("empty rows")
    abs_devs       = np.array([abs(r.dev_pct)              for r in rows], dtype=float)
    abs_devs_ad    = np.array([abs(r.dev_pct_ad)           for r in rows], dtype=float)
    abs_devs_ad_s  = np.array([abs(r.dev_pct_ad_substrate) for r in rows], dtype=float)
    n_pass         = int(sum(1 for r in rows if r.within_5pct))
    n_pass_ad      = int(sum(1 for r in rows if r.within_5pct_ad))
    n_pass_ad_s    = int(sum(1 for r in rows if r.within_5pct_ad_substrate))
    elem           = [r for r in rows if r.klass == "elemental"]
    elem_abs       = np.array([abs(r.dev_pct)              for r in elem], dtype=float)
    elem_abs_ad    = np.array([abs(r.dev_pct_ad)           for r in elem], dtype=float)
    elem_abs_ad_s  = np.array([abs(r.dev_pct_ad_substrate) for r in elem], dtype=float)
    elem_pass      = int(sum(1 for r in elem if r.within_5pct))
    elem_pass_ad   = int(sum(1 for r in elem if r.within_5pct_ad))
    elem_pass_ad_s = int(sum(1 for r in elem if r.within_5pct_ad_substrate))
    return {
        "n":                          len(rows),
        "n_within_5pct":              n_pass,
        "n_within_5pct_ad":           n_pass_ad,
        "n_within_5pct_ad_substrate": n_pass_ad_s,
        "mean_abs_dev":               float(abs_devs.mean()),
        "mean_abs_dev_ad":            float(abs_devs_ad.mean()),
        "mean_abs_dev_ad_substrate":  float(abs_devs_ad_s.mean()),
        "max_abs_dev":                float(abs_devs.max()),
        "max_abs_dev_ad":             float(abs_devs_ad.max()),
        "max_abs_dev_ad_substrate":   float(abs_devs_ad_s.max()),
        "min_abs_dev":                float(abs_devs.min()),
        "min_abs_dev_ad":             float(abs_devs_ad.min()),
        "min_abs_dev_ad_substrate":   float(abs_devs_ad_s.min()),
        "elemental_only": {
            "n":                          len(elem),
            "n_within_5pct":              elem_pass,
            "n_within_5pct_ad":           elem_pass_ad,
            "n_within_5pct_ad_substrate": elem_pass_ad_s,
            "mean_abs_dev":               float(elem_abs.mean())     if len(elem) else 0.0,
            "mean_abs_dev_ad":            float(elem_abs_ad.mean())  if len(elem) else 0.0,
            "mean_abs_dev_ad_substrate":  float(elem_abs_ad_s.mean())if len(elem) else 0.0,
            "max_abs_dev":                float(elem_abs.max())      if len(elem) else 0.0,
            "max_abs_dev_ad":             float(elem_abs_ad.max())   if len(elem) else 0.0,
            "max_abs_dev_ad_substrate":   float(elem_abs_ad_s.max()) if len(elem) else 0.0,
            "min_abs_dev":                float(elem_abs.min())      if len(elem) else 0.0,
            "min_abs_dev_ad":             float(elem_abs_ad.min())   if len(elem) else 0.0,
            "min_abs_dev_ad_substrate":   float(elem_abs_ad_s.min()) if len(elem) else 0.0,
        },
    }


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def render_bcs_gap_ratio_test(out_path: str) -> str:
    """Render a 3-series bar chart per material: bare BCS prediction,
    Allen-Dynes-corrected prediction, and measured 2Delta/(k_B T_c).

    The bottom panel shows percent deviation (bare vs AD) so the user can
    see at a glance which materials the strong-coupling correction
    rescues.  Saves PNG to ``out_path`` and returns the path.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = run_test()
    names    = [r.name        for r in rows]
    R_meas   = [r.R_meas      for r in rows]
    R_pred   = [r.R_pred      for r in rows]   # bare 3.528
    R_pred_ad= [r.R_pred_ad   for r in rows]
    devs     = [r.dev_pct     for r in rows]
    devs_ad  = [r.dev_pct_ad  for r in rows]
    klasses  = [r.klass       for r in rows]

    color_meas = "#1f77b4"  # tab:blue
    color_bare = "#7f7f7f"  # gray
    color_ad   = "#d62728"  # tab:red
    color_mb   = "#9467bd"  # tab:purple, multiband annotation
    edge_for_class = {"elemental": "black", "multiband": color_mb}

    fig, (ax1, ax2) = plt.subplots(
        nrows=2, ncols=1, figsize=(13.0, 9.0),
        gridspec_kw=dict(height_ratios=[3, 2], hspace=0.40),
    )

    # ---- top: 3-series grouped bar chart per material ---------------
    n = len(names)
    x = np.arange(n)
    bw = 0.27

    bars_bare = ax1.bar(
        x - bw, R_pred, width=bw, color=color_bare, edgecolor="black",
        linewidth=0.6, label=r"bare BCS  $2\pi/e^\gamma=3.528$",
    )
    bars_ad = ax1.bar(
        x,      R_pred_ad, width=bw, color=color_ad,   edgecolor="black",
        linewidth=0.6, label="substrate + Allen-Dynes",
    )
    bars_meas = ax1.bar(
        x + bw, R_meas,    width=bw, color=color_meas, edgecolor="black",
        linewidth=0.6, label="measured",
    )
    # outline multiband-class bars in purple to flag MgB_2
    for k, b in zip(klasses, bars_meas):
        if k == "multiband":
            b.set_edgecolor(color_mb)
            b.set_linewidth(2.0)

    ax1.axhline(BCS_RATIO_PRED, color=color_bare, linestyle="--",
                linewidth=1.0, alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(names)
    ax1.set_ylabel(r"$2\Delta(0) / (k_B T_c)$  (dimensionless)",
                   fontsize=11)
    ax1.set_title("BCS gap ratio: bare substrate / BCS  vs  "
                  "Allen-Dynes-corrected  vs  measured",
                  fontsize=12)
    ax1.set_ylim(0.0, 5.2)
    ax1.grid(axis="y", alpha=0.3)
    # numeric labels above measured bars
    for b, r in zip(bars_meas, rows):
        ax1.text(b.get_x() + b.get_width() / 2.0,
                 b.get_height() + 0.05,
                 f"{r.R_meas:.2f}",
                 ha="center", va="bottom", fontsize=7.5)
    ax1.legend(loc="upper left", fontsize=9, framealpha=0.95, ncol=3)

    # ---- bottom: percent deviation (bare vs AD) --------------------
    bars2_bare = ax2.bar(
        x - bw / 2.0, devs, width=bw, color=color_bare,
        edgecolor="black", linewidth=0.6, label="bare BCS deviation",
    )
    bars2_ad = ax2.bar(
        x + bw / 2.0, devs_ad, width=bw, color=color_ad,
        edgecolor="black", linewidth=0.6, label="Allen-Dynes deviation",
    )
    ax2.axhline(0.0, color="black", linewidth=0.8)
    ax2.axhspan(-5.0, 5.0, color="crimson", alpha=0.08,
                label=r"$\pm 5\%$ match band")
    ax2.set_xticks(x)
    ax2.set_xticklabels(names)
    ax2.set_ylabel("deviation from\nprediction (%)", fontsize=11)
    ax2.grid(axis="y", alpha=0.3)
    ax2.legend(loc="lower left", fontsize=9, ncol=3)
    for b, r in zip(bars2_bare, rows):
        y = b.get_height()
        offset = 0.6 if y >= 0 else -1.4
        ax2.text(b.get_x() + b.get_width() / 2.0,
                 y + offset,
                 f"{r.dev_pct:+.0f}",
                 ha="center", va="bottom" if y >= 0 else "top",
                 fontsize=7)
    for b, r in zip(bars2_ad, rows):
        y = b.get_height()
        offset = 0.6 if y >= 0 else -1.4
        ax2.text(b.get_x() + b.get_width() / 2.0,
                 y + offset,
                 f"{r.dev_pct_ad:+.0f}",
                 ha="center", va="bottom" if y >= 0 else "top",
                 fontsize=7)

    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _print_table(rows: List[ResultRow]) -> None:
    print(f"\nBCS gap ratio test: substrate prediction (bare) = "
          f"{BCS_RATIO_PRED:.6f}")
    print(f" {'name':<6}  {'T_c[K]':>7}  {'2D[meV]':>8}  "
          f"{'R_meas':>7}  "
          f"{'R_AD_e':>7}  {'R_AD_s':>7}  "
          f"{'dev_b':>7}  {'dev_ADe':>7}  {'dev_ADs':>7}  "
          f"{'<=5%b':>5}  {'<=5%ADe':>7}  {'<=5%ADs':>7}  class")
    for r in rows:
        flag_b   = "Y" if r.within_5pct              else " "
        flag_ad  = "Y" if r.within_5pct_ad           else " "
        flag_ads = "Y" if r.within_5pct_ad_substrate else " "
        print(f" {r.name:<6}  {r.T_c_K:>7.2f}  {r.gap_meV:>8.3f}  "
              f"{r.R_meas:>7.3f}  "
              f"{r.R_pred_ad:>7.3f}  {r.R_pred_ad_substrate:>7.3f}  "
              f"{r.dev_pct:>+7.2f}  {r.dev_pct_ad:>+7.2f}  "
              f"{r.dev_pct_ad_substrate:>+7.2f}  "
              f"{flag_b:>5}  {flag_ad:>7}  {flag_ads:>7}  {r.klass}")


def main() -> None:
    rows = run_test()
    _print_table(rows)
    stats = summary_stats(rows)
    print("\nAggregate (all materials):")
    print(f"  n={stats['n']}, "
          f"bare n_within_5pct={stats['n_within_5pct']}, "
          f"AD-emp n_within_5pct={stats['n_within_5pct_ad']}, "
          f"AD-sub n_within_5pct={stats['n_within_5pct_ad_substrate']}")
    print(f"  bare    mean|dev|={stats['mean_abs_dev']:.2f}%  max|dev|="
          f"{stats['max_abs_dev']:.2f}%")
    print(f"  AD-emp  mean|dev|={stats['mean_abs_dev_ad']:.2f}%  max|dev|="
          f"{stats['max_abs_dev_ad']:.2f}%")
    print(f"  AD-sub  mean|dev|={stats['mean_abs_dev_ad_substrate']:.2f}%  max|dev|="
          f"{stats['max_abs_dev_ad_substrate']:.2f}%")
    elem = stats["elemental_only"]
    print("Aggregate (elemental only):")
    print(f"  n={elem['n']}, "
          f"bare n_within_5pct={elem['n_within_5pct']}, "
          f"AD-emp n_within_5pct={elem['n_within_5pct_ad']}, "
          f"AD-sub n_within_5pct={elem['n_within_5pct_ad_substrate']}")
    print(f"  bare    mean|dev|={elem['mean_abs_dev']:.2f}%  max|dev|="
          f"{elem['max_abs_dev']:.2f}%")
    print(f"  AD-emp  mean|dev|={elem['mean_abs_dev_ad']:.2f}%  max|dev|="
          f"{elem['max_abs_dev_ad']:.2f}%")
    print(f"  AD-sub  mean|dev|={elem['mean_abs_dev_ad_substrate']:.2f}%  max|dev|="
          f"{elem['max_abs_dev_ad_substrate']:.2f}%")
    # MgB_2 multiband breakdown
    mb = mgb2_multiband_prediction()
    print("\nMgB_2 multiband (per-band, MEASURED):")
    print(f"  R_sigma     = {mb['R_sigma']:.3f}  (Choi 2002 Eliashberg "
          f"~ {mb['Choi_R_sigma_ref']:.2f})")
    print(f"  R_pi        = {mb['R_pi']:.3f}  (Choi 2002 Eliashberg "
          f"~ {mb['Choi_R_pi_ref']:.2f})")
    print(f"  R_composite = {mb['R_composite']:.3f}  "
          f"(DOS-weighted full gap)")
    print("\nMgB_2 multiband (per-band, SUBSTRATE K_4/Q_3 prediction, "
          "zero parameters):")
    print(f"  Delta_sigma = {mb['Delta_sigma_substrate_pred_meV']:.3f} meV  "
          f"(measured 7.0, dev "
          f"{(mb['Delta_sigma_substrate_pred_meV'] - 7.0)/7.0*100:+.2f}%)")
    print(f"  Delta_pi    = {mb['Delta_pi_substrate_pred_meV']:.3f} meV  "
          f"(measured 2.0, dev "
          f"{(mb['Delta_pi_substrate_pred_meV'] - 2.0)/2.0*100:+.2f}%)")
    print(f"  R_sigma     = {mb['R_sigma_substrate_pred']:.3f}  (Choi "
          f"~ {mb['Choi_R_sigma_ref']:.2f}, dev "
          f"{mb['dev_sigma_substrate_pct']:+.2f}%)")
    print(f"  R_pi        = {mb['R_pi_substrate_pred']:.3f}  (Choi "
          f"~ {mb['Choi_R_pi_ref']:.2f}, dev "
          f"{mb['dev_pi_substrate_pct']:+.2f}%)")
    print(f"  ratio sigma/pi = {mb['ratio_sigma_pi_substrate']:.3f}  "
          f"(measured 3.5)")


__all__ = [
    "ALLEN_DYNES_PREFACTOR",
    "BCS_RATIO_PRED",
    "EULER_GAMMA",
    "K_B_MEV_K",
    "MATERIALS",
    "MATERIAL_PHONON_PARAMS",
    "MaterialDatum",
    "MgB2Bands",
    "ResultRow",
    "deviation_percent",
    "main",
    "measured_ratio",
    "mgb2_multiband_prediction",
    "predict_with_allen_dynes",
    "render_bcs_gap_ratio_test",
    "run_test",
    "substrate_phonon_params",
    "summary_stats",
]


if __name__ == "__main__":  # pragma: no cover
    main()
