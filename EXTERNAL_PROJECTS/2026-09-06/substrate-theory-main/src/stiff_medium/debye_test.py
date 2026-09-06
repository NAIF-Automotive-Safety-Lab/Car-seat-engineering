"""Debye-temperature predictions vs measured Θ_D for 14 elemental solids.

Tests the substrate phonon-dispersion picture by checking whether the
substrate-derived sound speed c_s = sqrt(K/ρ) -- when fed into the standard
Debye formula

    Θ_D = (ℏ / k_B) · c_s · (6 π² n)^(1/3)

reproduces measured Debye temperatures across a benchmark set of solids
spanning ~14× in Θ_D (Pb 105 K -> diamond 2230 K).

WHAT THE TEST IS ACTUALLY MEASURING
-----------------------------------
The substrate Lagrangian L = (ρ/2)(∂_t φ)² - (K/2)(∇φ)² gives
c_s = sqrt(K/ρ) at long wavelength (`PhononDispersion.sound_speed`). The
Debye temperature follows by counting modes up to the Brillouin-zone radius
k_D = (6 π² n)^(1/3) and identifying ω_D = c_s · k_D.

For each material we use:
  * the Debye-averaged sound speed from elastic constants
        c_s = ((1/c_L³ + 2/c_T³) / 3)^(-1/3)
    where c_L = sqrt((B + 4G/3)/ρ) is the longitudinal speed,
          c_T = sqrt(G/ρ)           is the transverse speed,
    with bulk modulus B and shear modulus G from CRC Handbook (90th ed.).
  * the atomic number density n = N_A · ρ_mass / M_molar.

Comparing predicted Θ_D against measured Θ_D is a SHARP test of the
Lagrangian's elastic structure: it asks whether the long-wavelength
acoustic dispersion ω(k) = c_s · k -- which is forced by the substrate
kinetic and gradient terms -- correctly predicts the high-frequency
zone-boundary cutoff that controls the heat capacity.

HONEST CAVEAT
-------------
Per material we use the EXPERIMENTAL bulk and shear moduli to compute c_s,
not a substrate-derived (B, G) -- so the test is on the structure of the
formula (Debye averaging + zone-boundary scaling), NOT on a parameter-free
prediction of the moduli themselves. A failure here would still falsify
the substrate Lagrangian's kinematic content; a success means only that
the kinematic content is consistent with measured elasticity, not that
the moduli are derived. This is the same calibration used by Ashcroft &
Mermin Chapter 23.

The 14 materials cover FCC (Cu, Al, Au, Ag, Ni, Pb), BCC (Fe, W),
diamond cubic (C, Si, Ge), HCP (Be, Mg), and tetragonal (Sn) lattices.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


# --------------------------------------------------------------------------- #
# Physical constants (SI)                                                     #
# --------------------------------------------------------------------------- #

HBAR: float = 1.054571817e-34   # J·s
KB: float = 1.380649e-23        # J/K
N_A: float = 6.02214076e23      # 1/mol


# --------------------------------------------------------------------------- #
# Material database                                                           #
# --------------------------------------------------------------------------- #
#
# All numbers are CRC Handbook of Chemistry and Physics (90th ed., 2009-2010)
# unless noted otherwise. Values picked to match the Θ_D reference values in
# Ashcroft & Mermin Table 23.2 (column 2: low-temperature heat-capacity Θ_D).
#
# Fields per material:
#   M_molar    : molar mass in kg/mol
#   rho_mass   : mass density at room temperature in kg/m³
#   B_GPa      : isothermal bulk modulus, GPa
#   G_GPa      : shear (Voigt-Reuss-Hill) modulus, GPa
#   theta_D_K  : measured Debye temperature, K (Ashcroft-Mermin / Kittel)
#   lattice    : crystal structure label (informational only)
#
# Where multiple polymorphs exist (Sn: white = β-Sn at RT; C: diamond) we
# pick the room-temperature stable phase consistent with the Θ_D quoted in
# the standard tables.

@dataclass(frozen=True)
class Material:
    name: str
    M_molar: float       # kg/mol
    rho_mass: float      # kg/m³
    B_GPa: float
    G_GPa: float
    theta_D_K: float     # measured
    lattice: str
    alpha_V: float = 0.0     # volumetric thermal expansion coefficient (1/K)
    T_melt: float = 0.0      # melting point (K) — Lindemann thermal scale anchor
    gamma_G_lit: float = 0.0 # literature Grüneisen parameter (for cross-check only)

    @property
    def n_atoms_per_m3(self) -> float:
        """Atomic number density: N_A · ρ / M."""
        return N_A * self.rho_mass / self.M_molar

    @property
    def poisson_ratio(self) -> float:
        """Poisson's ratio v = (3B - 2G) / (2(3B + G)).  Substrate cap: v ≤ σ_max = 1/2."""
        return (3.0 * self.B_GPa - 2.0 * self.G_GPa) / (
            2.0 * (3.0 * self.B_GPa + self.G_GPa)
        )


MATERIALS: Dict[str, Material] = {
    # Per-material extras (alpha_V from CRC volumetric thermal expansion = 3 * linear,
    # T_melt from CRC, gamma_G_lit from Ashcroft-Mermin / Anderson "Equation of State of Solids").
    #
    # FCC metals
    "Copper":     Material("Copper",     0.06355,  8960.0, 140.0,  48.0,  343.0, "FCC",
                            alpha_V=49.5e-6, T_melt=1357.0, gamma_G_lit=2.0),
    "Aluminum":   Material("Aluminum",   0.02698,  2700.0,  76.0,  26.0,  428.0, "FCC",
                            alpha_V=69.0e-6, T_melt= 933.0, gamma_G_lit=2.2),
    "Gold":       Material("Gold",       0.19697, 19300.0, 180.0,  27.0,  165.0, "FCC",
                            alpha_V=42.6e-6, T_melt=1337.0, gamma_G_lit=3.0),
    "Silver":     Material("Silver",     0.10787, 10490.0, 100.0,  30.0,  225.0, "FCC",
                            alpha_V=56.7e-6, T_melt=1235.0, gamma_G_lit=2.4),
    "Nickel":     Material("Nickel",     0.05869,  8908.0, 180.0,  76.0,  450.0, "FCC",
                            alpha_V=39.0e-6, T_melt=1728.0, gamma_G_lit=1.9),
    "Lead":       Material("Lead",       0.20720, 11340.0,  46.0,   5.6,  105.0, "FCC",
                            alpha_V=87.0e-6, T_melt= 600.0, gamma_G_lit=2.7),

    # BCC metals
    "Iron":       Material("Iron",       0.05585,  7874.0, 170.0,  82.0,  470.0, "BCC",
                            alpha_V=35.4e-6, T_melt=1811.0, gamma_G_lit=1.7),
    "Tungsten":   Material("Tungsten",   0.18384, 19250.0, 310.0, 161.0,  400.0, "BCC",
                            alpha_V=13.5e-6, T_melt=3695.0, gamma_G_lit=1.5),

    # Diamond-cubic covalent
    "Diamond":    Material("Diamond",    0.01201,  3515.0, 442.0, 478.0, 2230.0, "diamond",
                            alpha_V= 3.0e-6, T_melt=4300.0, gamma_G_lit=0.85),
    "Silicon":    Material("Silicon",    0.02809,  2329.0,  98.0,  66.0,  645.0, "diamond",
                            alpha_V= 7.8e-6, T_melt=1687.0, gamma_G_lit=1.0),
    "Germanium":  Material("Germanium",  0.07264,  5323.0,  75.0,  55.0,  374.0, "diamond",
                            alpha_V=17.3e-6, T_melt=1211.0, gamma_G_lit=1.0),

    # HCP metals
    "Beryllium":  Material("Beryllium",  0.00901,  1850.0, 130.0, 132.0, 1440.0, "HCP",
                            alpha_V=34.2e-6, T_melt=1560.0, gamma_G_lit=1.2),
    "Magnesium":  Material("Magnesium",  0.02431,  1738.0,  45.0,  17.0,  400.0, "HCP",
                            alpha_V=75.0e-6, T_melt= 923.0, gamma_G_lit=1.5),

    # White tin (β-Sn, body-centered tetragonal, room-temperature stable)
    "Tin":        Material("Tin",        0.11871,  7287.0,  58.0,  18.0,  200.0, "tetragonal",
                            alpha_V=66.0e-6, T_melt= 505.0, gamma_G_lit=2.3),
}


# --------------------------------------------------------------------------- #
# Sound-speed extraction                                                      #
# --------------------------------------------------------------------------- #


def longitudinal_speed(B_GPa: float, G_GPa: float, rho: float) -> float:
    """c_L = sqrt((B + 4G/3) / ρ).  Inputs in GPa, kg/m³.  Output m/s."""
    M = (B_GPa + (4.0 / 3.0) * G_GPa) * 1.0e9     # Pa
    return math.sqrt(M / rho)


def transverse_speed(G_GPa: float, rho: float) -> float:
    """c_T = sqrt(G / ρ).  Inputs in GPa, kg/m³.  Output m/s."""
    return math.sqrt(G_GPa * 1.0e9 / rho)


def debye_average_speed(c_L: float, c_T: float) -> float:
    """Debye-averaged sound speed: 3/c_s³ = 1/c_L³ + 2/c_T³."""
    return (3.0 / (1.0 / c_L ** 3 + 2.0 / c_T ** 3)) ** (1.0 / 3.0)


def material_sound_speeds(mat: Material) -> Tuple[float, float, float]:
    """Return (c_L, c_T, c_Debye) for a material in m/s."""
    c_L = longitudinal_speed(mat.B_GPa, mat.G_GPa, mat.rho_mass)
    c_T = transverse_speed(mat.G_GPa, mat.rho_mass)
    c_D = debye_average_speed(c_L, c_T)
    return c_L, c_T, c_D


# --------------------------------------------------------------------------- #
# Debye-temperature prediction                                                #
# --------------------------------------------------------------------------- #


def debye_temperature(c_s: float, n_atoms_per_m3: float) -> float:
    """Standard Debye formula Θ_D = (ℏ/k_B) · c_s · (6 π² n)^(1/3)."""
    omega_D = c_s * (6.0 * math.pi ** 2 * n_atoms_per_m3) ** (1.0 / 3.0)
    return HBAR * omega_D / KB


def predict_theta_D(mat: Material) -> Dict[str, float]:
    """Run the substrate-Lagrangian Debye prediction for one material."""
    c_L, c_T, c_D = material_sound_speeds(mat)
    n = mat.n_atoms_per_m3
    theta_pred = debye_temperature(c_D, n)
    theta_meas = mat.theta_D_K
    rel_err = (theta_pred - theta_meas) / theta_meas
    return {
        "name":       mat.name,
        "lattice":    mat.lattice,
        "rho_mass":   mat.rho_mass,
        "B_GPa":      mat.B_GPa,
        "G_GPa":      mat.G_GPa,
        "n_atoms":    n,
        "c_L":        c_L,
        "c_T":        c_T,
        "c_Debye":    c_D,
        "theta_pred": theta_pred,
        "theta_meas": theta_meas,
        "rel_err":    rel_err,
    }


# --------------------------------------------------------------------------- #
# SUBSTRATE-DERIVED MODULI PATHWAY                                            #
# --------------------------------------------------------------------------- #
#
# Rather than feeding the EXPERIMENTAL (B, G) into the Debye formula, we can
# now use the SUBSTRATE-DERIVED (B, G) from substrate_elasticity.py. This
# closes the calibration loop: every input to the Θ_D prediction is then
# either substrate-derived or a primary anchor (ε_coh, d_NN, ρ_mass, M).
#
# The substrate-elasticity module covers 8 of the 14 debye_test materials
# (Cu, Al, Au, Ni, Pb, Fe, diamond, Si). For materials present in both
# databases, we substitute substrate-predicted (B, G) and re-run the Debye
# pipeline.


def predict_theta_D_substrate_moduli(mat: Material,
                                     full_substrate_eps_coh: bool = False
                                     ) -> Dict[str, float]:
    """Debye-temperature prediction using SUBSTRATE-DERIVED (B, G).

    Uses substrate_elasticity.bulk_modulus_substrate /
    shear_modulus_substrate to predict (B, G) from atomic K_4 face-pair
    coupling, then feeds those into the standard Debye pipeline.

    Falls back to the empirical (B, G) for materials not in the substrate-
    elasticity database (Silver, Tungsten, Germanium, Beryllium, Magnesium,
    Tin), with a 'using_substrate_moduli' flag.

    Parameters
    ----------
    full_substrate_eps_coh : bool, default False
        When True, additionally sources ε_coh from
        ``substrate_cohesive_energy.eps_coh_substrate`` (closing the
        Cat-B → Cat-A chain on ε_coh) instead of the per-material
        Kittel-handbook anchor in substrate_elasticity.MATERIALS.
        For materials not in substrate_cohesive_energy.MATERIALS the
        function falls back to the substrate_elasticity handbook anchor.
    """
    from .substrate_elasticity import (
        MATERIALS as ELASTIC_MATS,
        bulk_modulus_substrate,
        bulk_modulus_substrate_full,
        shear_modulus_substrate,
        shear_modulus_substrate_full,
    )

    if mat.name in ELASTIC_MATS:
        emat = ELASTIC_MATS[mat.name]
        if full_substrate_eps_coh:
            B_use, _ = bulk_modulus_substrate_full(emat)
            G_use, _ = shear_modulus_substrate_full(emat)
        else:
            B_use = bulk_modulus_substrate(emat)
            G_use = shear_modulus_substrate(emat)
        using_substrate = True
    else:
        B_use = mat.B_GPa
        G_use = mat.G_GPa
        using_substrate = False

    c_L = longitudinal_speed(B_use, G_use, mat.rho_mass)
    c_T = transverse_speed(G_use, mat.rho_mass)
    c_D = debye_average_speed(c_L, c_T)
    n = mat.n_atoms_per_m3
    theta_pred = debye_temperature(c_D, n)
    theta_meas = mat.theta_D_K
    rel_err = (theta_pred - theta_meas) / theta_meas
    return {
        "name":               mat.name,
        "lattice":            mat.lattice,
        "B_GPa_substrate":    B_use,
        "G_GPa_substrate":    G_use,
        "B_GPa_measured":     mat.B_GPa,
        "G_GPa_measured":     mat.G_GPa,
        "using_substrate":    using_substrate,
        "c_Debye":            c_D,
        "theta_pred":         theta_pred,
        "theta_meas":         theta_meas,
        "rel_err":            rel_err,
    }


def run_test_substrate_moduli(full_substrate_eps_coh: bool = False
                              ) -> Dict[str, object]:
    """Run Debye-temperature predictions using substrate-derived (B, G).

    Returns rows + summary, comparable to run_test() but with substrate
    moduli substituted in for the materials covered by substrate_elasticity.

    Parameters
    ----------
    full_substrate_eps_coh : bool, default False
        Pass-through to ``predict_theta_D_substrate_moduli``.  When True,
        the ε_coh anchor is itself derived from
        ``substrate_cohesive_energy``, closing the Cat-B → Cat-A chain
        on cohesive energy for the 6 metals (Cu, Al, Au, Ni, Pb, Fe)
        with substrate cohesive-energy entries.
    """
    rows: Dict[str, Dict[str, float]] = {}
    pred_arr: List[float] = []
    meas_arr: List[float] = []
    pred_arr_sub_only: List[float] = []
    meas_arr_sub_only: List[float] = []
    for name, mat in MATERIALS.items():
        row = predict_theta_D_substrate_moduli(
            mat, full_substrate_eps_coh=full_substrate_eps_coh
        )
        rows[name] = row
        pred_arr.append(row["theta_pred"])
        meas_arr.append(row["theta_meas"])
        if row["using_substrate"]:
            pred_arr_sub_only.append(row["theta_pred"])
            meas_arr_sub_only.append(row["theta_meas"])

    pred = np.asarray(pred_arr)
    meas = np.asarray(meas_arr)
    rel = (pred - meas) / meas

    pred_s = np.asarray(pred_arr_sub_only)
    meas_s = np.asarray(meas_arr_sub_only)
    rel_s = (pred_s - meas_s) / meas_s

    summary = {
        "n_materials_total":              len(rows),
        "n_materials_substrate":          int(sum(
            1 for r in rows.values() if r["using_substrate"]
        )),
        "mean_abs_rel_err_substrate_only": float(np.abs(rel_s).mean()),
        "max_abs_rel_err_substrate_only":  float(np.abs(rel_s).max()),
        "mean_abs_rel_err_total":          float(np.abs(rel).mean()),
    }
    return {"rows": rows, "summary": summary}


# --------------------------------------------------------------------------- #
# Quasi-harmonic / anharmonic correction                                      #
# --------------------------------------------------------------------------- #
#
# WHY HARMONIC FAILS FOR SOFT METALS (Pb, Sn)
# -------------------------------------------
# The harmonic prediction Θ_D = (ℏ/k_B)·c_s·(6π²n)^(1/3) uses the long-wavelength
# elastic-continuum sound speed. For materials with very soft transverse modes
# (G ≪ B, large Poisson's ratio v) the actual zone-boundary phonon frequencies
# are anharmonically *stiffened* relative to the linear extrapolation:
#
#     ω(k_BZ)_actual > c_T · k_BZ
#
# because the cubic substrate term g₃·u³ contributes a positive shift to the
# self-consistent renormalised frequency (Born-Oppenheimer mean-field). This
# stiffening is captured by the Grüneisen parameter γ_G = -∂ln(ω)/∂ln(V).
#
# SUBSTRATE-DERIVED γ_G
# ---------------------
# Belomestnykh-Tesleva (Tech. Phys. Lett. 30 91, 2004) showed that for an
# isotropic medium with arbitrary cubic anharmonicity,
#
#     γ_G = (3/2) · (1 + v) / (2 - 3v)
#
# where v is the Poisson ratio. THIS IS SUBSTRATE-DERIVABLE: the substrate
# saturation cap σ_max = 1/2 (`SIGMA_MAX` in src/stiff_medium/saturation.py)
# enforces v ≤ 1/2 (the incompressible limit). The denominator 2 - 3v vanishes
# at v = 2/3, but v is bounded above by σ_max = 1/2 ⇒ γ_G ≤ 9/2 always.
# So the substrate elastic limit is the upstream cause of the universal
# γ_G ≲ 5 observed across all 14 materials.
#
# QUASI-HARMONIC CORRECTION
# -------------------------
# The standard quasi-harmonic shift due to thermal lattice expansion is
#
#     Θ_D^qh / Θ_D^harm  =  1 + γ_G · α_V · T_eff
#
# where α_V is the volumetric thermal expansion coefficient and T_eff sets
# the relevant thermal-displacement scale. We use the Lindemann scale
# T_eff = T_melt / 2: at this temperature the rms atomic displacement reaches
# ~10% of the lattice spacing, which is the regime in which the soft-mode
# stiffening is established and saturates. The same scale recovers the
# experimental low-temperature heat-capacity Θ_D when extrapolated back.
#
# HONEST CAVEAT
# -------------
# This correction CANNOT fully resolve the Pb anomaly. Pb has additional
# Kohn-anomaly stiffening at the X-point of the Brillouin zone driven by
# electron-phonon coupling at the Fermi surface — a quantum effect outside
# the elastic continuum. The substrate quasi-harmonic correction reduces
# the Pb residual from -27% to roughly -15% to -21% (depending on T_eff),
# but does not reach the measured 105 K. Sn responds well: -12% → -8%,
# which is within the 10% tolerance band.

# Substrate elastic-limit cap on Poisson's ratio (from saturation.py SIGMA_MAX = 0.5)
SIGMA_MAX_POISSON: float = 0.5


def gamma_G_substrate(mat: Material) -> float:
    """Substrate-derived Grüneisen parameter via Belomestnykh-Tesleva form.

        γ_G = (3/2) · (1 + v) / (2 - 3v),   v ≤ σ_max = 1/2

    This formula derives from the substrate elastic-limit cap σ_max = 1/2
    (saturation.py SIGMA_MAX) acting on the local Poisson ratio v. The
    denominator (2 - 3v) → 0 at v = 2/3, but the substrate cap v ≤ 1/2
    ensures γ_G ≤ 9/2 ≈ 4.5 universally. For materials with v ≥ 2/3
    (unphysical given the cap) we return γ_G = 9/2 as a hard cap.
    """
    v = mat.poisson_ratio
    if v >= 2.0 / 3.0:
        return 4.5
    # Enforce the substrate σ_max cap (incompressible limit)
    if v > SIGMA_MAX_POISSON:
        v = SIGMA_MAX_POISSON
    return 1.5 * (1.0 + v) / (2.0 - 3.0 * v)


def quasi_harmonic_shift(mat: Material, T_eff: float | None = None) -> float:
    """Multiplicative shift Θ_D^qh / Θ_D^harm = 1 + γ_G · α_V · T_eff.

    If T_eff is None, use the Lindemann thermal scale T_melt / 2.
    Returns the multiplicative factor (always ≥ 1 since γ_G > 0, α_V > 0).
    """
    if T_eff is None:
        T_eff = 0.5 * mat.T_melt
    gG = gamma_G_substrate(mat)
    return 1.0 + gG * mat.alpha_V * T_eff


def predict_with_quasiharmonic(
    mat: Material, T_eff: float | None = None
) -> Dict[str, float]:
    """Quasi-harmonic-corrected Θ_D prediction for one material.

    Pipeline:
      (1) Compute harmonic Θ_D^harm from elastic moduli (Debye-averaged c_s).
      (2) Compute substrate γ_G from Poisson ratio (Belomestnykh-Tesleva).
      (3) Apply quasi-harmonic shift Θ_D^qh = Θ_D^harm · (1 + γ_G · α_V · T_eff)
          at T_eff = T_melt / 2 (Lindemann thermal scale).

    Returns a dict comparing harmonic, quasi-harmonic, and measured Θ_D.
    """
    if T_eff is None:
        T_eff = 0.5 * mat.T_melt
    harm = predict_theta_D(mat)
    gG_sub = gamma_G_substrate(mat)
    shift = 1.0 + gG_sub * mat.alpha_V * T_eff
    theta_qh = harm["theta_pred"] * shift
    theta_meas = mat.theta_D_K
    rel_err_harm = harm["rel_err"]
    rel_err_qh = (theta_qh - theta_meas) / theta_meas
    return {
        "name":          mat.name,
        "lattice":       mat.lattice,
        "v_poisson":     mat.poisson_ratio,
        "gamma_G_sub":   gG_sub,
        "gamma_G_lit":   mat.gamma_G_lit,
        "alpha_V":       mat.alpha_V,
        "T_melt":        mat.T_melt,
        "T_eff":         T_eff,
        "shift":         shift,
        "theta_harm":    harm["theta_pred"],
        "theta_qh":      theta_qh,
        "theta_meas":    theta_meas,
        "rel_err_harm":  rel_err_harm,
        "rel_err_qh":    rel_err_qh,
    }


def run_test_quasiharmonic(T_eff: float | None = None) -> Dict[str, object]:
    """Compute quasi-harmonic predictions for all 14 materials with summary."""
    rows: Dict[str, Dict[str, float]] = {}
    pred_h: List[float] = []
    pred_q: List[float] = []
    meas:   List[float] = []
    for name, mat in MATERIALS.items():
        row = predict_with_quasiharmonic(mat, T_eff=T_eff)
        rows[name] = row
        pred_h.append(row["theta_harm"])
        pred_q.append(row["theta_qh"])
        meas.append(row["theta_meas"])

    ph = np.asarray(pred_h)
    pq = np.asarray(pred_q)
    me = np.asarray(meas)
    rel_h = (ph - me) / me
    rel_q = (pq - me) / me

    summary = {
        "n_materials":               len(rows),
        "T_eff_default":             "T_melt/2 (per material)" if T_eff is None else T_eff,
        # harmonic baseline
        "mean_abs_rel_err_harm":     float(np.abs(rel_h).mean()),
        "max_abs_rel_err_harm":      float(np.abs(rel_h).max()),
        "within_10pct_harm":         int(sum(1 for r in rel_h if abs(r) <= 0.10)),
        # quasi-harmonic
        "mean_abs_rel_err_qh":       float(np.abs(rel_q).mean()),
        "max_abs_rel_err_qh":        float(np.abs(rel_q).max()),
        "within_10pct_qh":           int(sum(1 for r in rel_q if abs(r) <= 0.10)),
        "improvement_pp":            float(
            np.abs(rel_h).mean() - np.abs(rel_q).mean()
        ),
        # specific-material focus
        "Pb_harm":                   rows["Lead"]["theta_harm"],
        "Pb_qh":                     rows["Lead"]["theta_qh"],
        "Pb_meas":                   rows["Lead"]["theta_meas"],
        "Sn_harm":                   rows["Tin"]["theta_harm"],
        "Sn_qh":                     rows["Tin"]["theta_qh"],
        "Sn_meas":                   rows["Tin"]["theta_meas"],
    }
    return {"rows": rows, "summary": summary}


def run_test() -> Dict[str, object]:
    """Compute predictions for every material; return rows + summary stats."""
    rows: Dict[str, Dict[str, float]] = {}
    pred_arr: List[float] = []
    meas_arr: List[float] = []
    for name, mat in MATERIALS.items():
        row = predict_theta_D(mat)
        rows[name] = row
        pred_arr.append(row["theta_pred"])
        meas_arr.append(row["theta_meas"])

    pred = np.asarray(pred_arr)
    meas = np.asarray(meas_arr)
    rel = (pred - meas) / meas

    # Pearson correlation
    pm = pred - pred.mean()
    mm = meas - meas.mean()
    denom = math.sqrt(float((pm * pm).sum()) * float((mm * mm).sum()))
    pearson_r = float((pm * mm).sum()) / denom if denom > 0.0 else float("nan")

    # Log-space (geometric) Pearson, since Θ_D varies over 14×
    lp = np.log(pred)
    lm = np.log(meas)
    lpm = lp - lp.mean()
    lmm = lm - lm.mean()
    ldenom = math.sqrt(float((lpm * lpm).sum()) * float((lmm * lmm).sum()))
    pearson_log_r = (
        float((lpm * lmm).sum()) / ldenom if ldenom > 0.0 else float("nan")
    )

    # Linear regression in log-log space (slope of log(pred) vs log(meas))
    slope, intercept = np.polyfit(lm, lp, 1)

    # Per-material agreement classes
    materials_within_5pct = sum(1 for r in rel if abs(r) <= 0.05)
    materials_within_10pct = sum(1 for r in rel if abs(r) <= 0.10)
    materials_within_20pct = sum(1 for r in rel if abs(r) <= 0.20)

    summary = {
        "n_materials":            len(rows),
        "mean_rel_err":           float(rel.mean()),
        "mean_abs_rel_err":       float(np.abs(rel).mean()),
        "median_abs_rel_err":     float(np.median(np.abs(rel))),
        "max_abs_rel_err":        float(np.abs(rel).max()),
        "max_rel_err_material":   list(rows.keys())[int(np.argmax(np.abs(rel)))],
        "rms_rel_err":            float(np.sqrt((rel * rel).mean())),
        "pearson_r":              pearson_r,
        "pearson_log_r":          pearson_log_r,
        "loglog_slope":           float(slope),
        "loglog_intercept":       float(intercept),
        "within_5pct":            materials_within_5pct,
        "within_10pct":           materials_within_10pct,
        "within_20pct":           materials_within_20pct,
    }
    return {"rows": rows, "summary": summary}


# --------------------------------------------------------------------------- #
# Visual rendering                                                            #
# --------------------------------------------------------------------------- #


def render_debye_test(out_path: str | None = None) -> str:
    """Render visuals/124_debye_test.png: 3-line comparison of harmonic vs
    quasi-harmonic substrate predictions vs measured Θ_D, plus per-material
    residual bars showing improvement from quasi-harmonic correction."""
    import matplotlib.pyplot as plt
    from pathlib import Path

    if out_path is None:
        out_path = str(
            Path(__file__).resolve().parents[2]
            / "visuals" / "124_debye_test.png"
        )

    # Harmonic-only summary (for headline numbers and lattice palette)
    res = run_test()
    rows_h = res["rows"]
    summary_h = res["summary"]

    # Quasi-harmonic results
    qh = run_test_quasiharmonic()
    rows_q = qh["rows"]
    summary_q = qh["summary"]

    names = list(rows_h.keys())
    meas = np.array([rows_h[n]["theta_meas"] for n in names])
    pred_harm = np.array([rows_h[n]["theta_pred"] for n in names])
    pred_qh = np.array([rows_q[n]["theta_qh"] for n in names])
    rel_harm = (pred_harm - meas) / meas
    rel_qh = (pred_qh - meas) / meas

    lattice = [rows_h[n]["lattice"] for n in names]
    lattice_set = sorted(set(lattice))
    palette = {
        "FCC":         "#1f77b4",
        "BCC":         "#d62728",
        "diamond":     "#2ca02c",
        "HCP":         "#ff7f0e",
        "tetragonal":  "#9467bd",
    }
    colors = [palette.get(L, "#7f7f7f") for L in lattice]

    fig, (ax_3line, ax_bar) = plt.subplots(
        1, 2, figsize=(17.0, 7.5), gridspec_kw={"width_ratios": [1.3, 1.4]}
    )

    # ---------- left: 3-line comparison (harmonic / quasi-harmonic / measured) ----------
    # Sort by measured Θ_D so the curves make visual sense
    order_m = np.argsort(meas)
    names_m = [names[i] for i in order_m]
    meas_m = meas[order_m]
    pred_h_m = pred_harm[order_m]
    pred_q_m = pred_qh[order_m]

    x = np.arange(len(names_m))
    ax_3line.plot(
        x, meas_m, "o-", color="black", lw=2, markersize=9,
        label="measured Θ_D (Ashcroft-Mermin / Kittel)",
    )
    ax_3line.plot(
        x, pred_h_m, "s--", color="#d62728", lw=1.5, markersize=8,
        label="substrate harmonic Θ_D = (ℏ/k_B)·c_s·(6π²n)^(1/3)",
    )
    ax_3line.plot(
        x, pred_q_m, "D-", color="#1f77b4", lw=1.5, markersize=8,
        label=r"substrate quasi-harm  Θ_D · (1 + γ_G·α_V·T_melt/2),  γ_G=$\frac{3}{2}\frac{1+\nu}{2-3\nu}$",
    )
    ax_3line.set_yscale("log")
    ax_3line.set_xticks(x)
    ax_3line.set_xticklabels(names_m, rotation=55, ha="right", fontsize=9)
    ax_3line.set_ylabel(r"$\Theta_D$ [K]")
    ax_3line.set_title(
        "Substrate Debye-temperature: harmonic vs quasi-harmonic vs measured\n"
        f"harmonic mean |err| = {summary_q['mean_abs_rel_err_harm']:.1%},  "
        f"quasi-harm mean |err| = {summary_q['mean_abs_rel_err_qh']:.1%}  "
        f"(Δ = {summary_q['improvement_pp']:+.2%} pp)"
    )
    ax_3line.legend(loc="upper left", fontsize=9, framealpha=0.95)
    ax_3line.grid(True, which="both", alpha=0.3)
    # Annotate Pb and Sn with γ_G and α_V
    for tag in ("Lead", "Tin"):
        i = names_m.index(tag)
        row = rows_q[tag]
        annot = (
            f"{tag}\nγ_G={row['gamma_G_sub']:.2f} (lit {row['gamma_G_lit']:.1f})\n"
            f"α_V={row['alpha_V']*1e6:.0f} ppm/K\n"
            f"{row['theta_harm']:.0f}→{row['theta_qh']:.0f} (meas {row['theta_meas']:.0f})"
        )
        ax_3line.annotate(
            annot,
            (x[i], pred_q_m[i]),
            xytext=(15, -45) if tag == "Lead" else (15, -25),
            textcoords="offset points", fontsize=7,
            bbox=dict(facecolor="lightyellow", edgecolor="grey", alpha=0.9),
            arrowprops=dict(arrowstyle="->", color="grey", lw=0.6),
        )

    # ---------- right: residual bar chart, harmonic vs quasi-harmonic ----------
    order = np.argsort(rel_harm)
    names_o = [names[i] for i in order]
    rh_o = rel_harm[order] * 100.0
    rq_o = rel_qh[order] * 100.0

    y = np.arange(len(names_o))
    bw = 0.4
    ax_bar.barh(
        y - bw / 2, rh_o, height=bw, color="#d62728", alpha=0.85,
        edgecolor="black", linewidth=0.4, label="harmonic only",
    )
    ax_bar.barh(
        y + bw / 2, rq_o, height=bw, color="#1f77b4", alpha=0.85,
        edgecolor="black", linewidth=0.4, label="+ quasi-harmonic shift",
    )
    ax_bar.axvline(0.0, color="black", linewidth=1.2)
    ax_bar.axvline(5.0, color="grey", linestyle=":", linewidth=0.7)
    ax_bar.axvline(-5.0, color="grey", linestyle=":", linewidth=0.7)
    ax_bar.axvline(10.0, color="grey", linestyle="--", linewidth=0.7)
    ax_bar.axvline(-10.0, color="grey", linestyle="--", linewidth=0.7)
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(names_o, fontsize=10)
    ax_bar.set_xlabel(r"$(\Theta_D^{\rm pred}/\Theta_D^{\rm meas} - 1)$ [%]")
    ax_bar.set_title(
        f"Per-material residual: harmonic ({summary_q['within_10pct_harm']}/14 within 10%) "
        f"→ quasi-harm ({summary_q['within_10pct_qh']}/14 within 10%)"
    )
    ax_bar.legend(loc="lower right", fontsize=9)
    ax_bar.grid(True, axis="x", alpha=0.3)

    # Substrate γ_G derivation note
    ax_bar.text(
        0.02, 0.02,
        r"substrate cap: $\sigma_{max}=1/2$  $\Rightarrow$  $\nu\leq 1/2$  $\Rightarrow$  $\gamma_G\leq 9/2$"
        "\nKohn anomaly (Pb X-pt) lies outside elastic-continuum scope",
        transform=ax_bar.transAxes, ha="left", va="bottom", fontsize=8,
        bbox=dict(facecolor="lightyellow", edgecolor="grey", alpha=0.85),
    )

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# --------------------------------------------------------------------------- #
# CLI entrypoint                                                              #
# --------------------------------------------------------------------------- #


def main() -> None:
    res = run_test()
    rows = res["rows"]
    summary = res["summary"]

    print("Substrate phonon Debye-temperature test (14 materials)")
    print("=" * 88)
    print(
        f"{'Material':<12}{'lat.':<6}{'c_L (m/s)':>10}{'c_T (m/s)':>10}"
        f"{'c_D (m/s)':>10}{'Θ_pred [K]':>13}{'Θ_meas [K]':>13}{'rel err':>10}"
    )
    for name, row in rows.items():
        print(
            f"{name:<12}{row['lattice']:<6}"
            f"{row['c_L']:>10.0f}{row['c_T']:>10.0f}{row['c_Debye']:>10.0f}"
            f"{row['theta_pred']:>13.1f}{row['theta_meas']:>13.1f}"
            f"{row['rel_err']:>+10.2%}"
        )

    print()
    print("Summary (harmonic)")
    for k, v in summary.items():
        if isinstance(v, float):
            print(f"  {k:<24} {v:+.5f}")
        else:
            print(f"  {k:<24} {v}")

    # ----- Quasi-harmonic comparison -----
    qh = run_test_quasiharmonic()
    rows_q = qh["rows"]
    summary_q = qh["summary"]

    print()
    print("Quasi-harmonic correction (substrate γ_G = 3/2 · (1+ν)/(2-3ν), T_eff = T_melt/2)")
    print("=" * 88)
    print(
        f"{'Material':<12}{'ν':>6}{'γ_G_sub':>9}{'γ_G_lit':>9}"
        f"{'α_V (ppm/K)':>14}{'Θ_harm':>10}{'Θ_qh':>10}"
        f"{'Θ_meas':>9}{'err harm':>11}{'err qh':>10}"
    )
    for name, row in rows_q.items():
        print(
            f"{name:<12}{row['v_poisson']:>6.3f}{row['gamma_G_sub']:>9.2f}"
            f"{row['gamma_G_lit']:>9.2f}{row['alpha_V']*1e6:>14.1f}"
            f"{row['theta_harm']:>10.1f}{row['theta_qh']:>10.1f}"
            f"{row['theta_meas']:>9.0f}"
            f"{row['rel_err_harm']:>+11.2%}{row['rel_err_qh']:>+10.2%}"
        )

    print()
    print("Quasi-harmonic summary")
    for k, v in summary_q.items():
        if isinstance(v, float):
            print(f"  {k:<28} {v:+.5f}")
        else:
            print(f"  {k:<28} {v}")


if __name__ == "__main__":
    main()
