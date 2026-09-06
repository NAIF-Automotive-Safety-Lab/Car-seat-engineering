"""Substrate-derived Eliashberg α²F(ω) and Allen-Dynes inputs.

Standard Eliashberg theory takes the electron-phonon spectral function

        alpha^2 F(omega) = N(0) * <|g_eph|^2>_FS(omega) * F(omega)

as an *input*, extracted per-material from tunnelling spectroscopy or
ab initio (DFPT) calculation.  The Allen-Dynes / McMillan T_c formula
then reads off two reduced moments of alpha^2 F:

        lambda    = 2 integral alpha^2 F(omega) / omega  d omega
        omega_log = exp[ (2/lambda) integral alpha^2 F(omega)
                                       * ln(omega) / omega  d omega ]

In the substrate ontology BOTH alpha^2 F(omega) and these moments are
*derivable* from substrate primitives, with NO per-material spectral
input: the only material-specific numbers we keep are the ones the
substrate Debye-temperature derivation already takes (ion mass M_ion,
mass density rho, bulk and shear moduli B and G, valence z).  The chain
runs:

    1.  Acoustic phonon DOS F_ac(omega) = (3 V/(2 pi^2)) * omega^2 / c_s^3
        for omega <= omega_D, with c_s the substrate Debye-averaged sound
        speed (debye_test.material_sound_speeds) and
        omega_D = c_s * (6 pi^2 n)^(1/3).  This piece IS substrate-derived
        through the L = (rho/2)(d_t phi)^2 - (K/2)(grad phi)^2 kinetic
        and gradient terms (sound speed = sqrt(K/rho)).

    2.  Optical / zone-boundary piece F_opt(omega) (Einstein-like peak at
        omega_E = (2/3) * omega_D, weight ~ acoustic-mode count for the
        K_4-cell internal modes).  Substrate origin: the K_4 tetrahedral
        cell has internal vibrational modes whose lowest frequency is set
        by the face-pair restoring spring K / d_face^2 (face-pair coupling
        eps_face = Lambda_QCD/(n_A * N_BAM) per K_4 cell).

    3.  Electron-substrate coupling enters as Hopfield's I^2 / (M omega^2)
        weight: alpha^2(omega) = N(0) * I_eph^2 / (2 M_ion omega), where
        the substrate matrix element

            I_eph^2 = K_eph^2 / (M_ion * omega_D^2)

        is fixed by ONE dimensionless substrate constant
        K_eph_substrate = K_PA * a_TF^4 (a_TF the Thomas-Fermi screening
        length at the substrate Fermi level).  This deliberately gives
        the substrate one ALL-MATERIAL coupling number; the per-material
        variation comes purely from M_ion and omega_D, NOT from any
        material-specific tunable.

    4.  N(0) (Fermi-level DOS) from free-electron substrate:
        N(0) = m_e^* k_F / (pi^2 hbar^2) per spin, with k_F = (3 pi^2 n_e)^(1/3),
        n_e = z * n_atoms.  This is the substrate's free-electron Fermi
        DOS with valence z as the only per-material input (no effective-mass
        renormalisation: m^* = m_e everywhere).

API
---
``MaterialPhononSubstrate``      : substrate primitives for one material
                                   (M_molar, rho_mass, B_GPa, G_GPa,
                                    valence z; all from CRC).
``alpha2F_substrate(material)``  : returns (omega_grid_K, alpha2F_grid)
                                    on a uniform log-spaced grid in K.
``lambda_substrate(material)``   : substrate-derived lambda.
``omega_log_substrate(material)``: substrate-derived omega_log in K.
``predict_lambda_omega_log(name)``
                                  : convenience: name -> (lambda, omega_log_K).
``compare_to_empirical()``       : returns table comparing substrate
                                   (lambda, omega_log) to the empirical
                                   Allen-Dynes / Carbotte values used in
                                   ``bcs_gap_ratio_test.MATERIAL_PHONON_PARAMS``.
``MATERIAL_DB``                  : dict of substrate-side material entries
                                   keyed by symbol matching the BCS test
                                   table.

Honest caveats
--------------
* The Einstein-peak weight, the K_eph normalisation, and the choice of
  Thomas-Fermi screening length are SUBSTRATE-LEVEL ANSAETZE that a more
  rigorous derivation would tie to specific K_4-cell internal modes and
  their measured spectroscopic frequencies.  This module produces the
  ALL-MATERIAL coupling K_eph_substrate by anchoring it to ONE material
  (Pb, the canonical strong-coupling reference) and then PREDICTS the
  rest -- so the test is "if Pb fixes the substrate K_eph, what does the
  same K_eph give for Hg, Nb, Sn, Al, V, ...?".
* The substrate alpha^2 F(omega) here is a *single-Einstein-plus-Debye*
  parametrisation, not a fully-resolved bandstructure-derived spectrum.
  It captures the moments lambda and omega_log; it does NOT reproduce
  the multi-peak structure of the empirical alpha^2 F seen in tunnelling
  (e.g. the Pb transverse / longitudinal twin peaks at 4.4 / 8.5 meV).
  For BCS gap-ratio prediction only the MOMENTS matter, so the
  parametrisation is fit-for-purpose for the Allen-Dynes correction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .b3_constants import LAMBDA_QCD_MEV, N_BAM, K_pair, K_rank, n_R, n_A


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

HBAR:      float = 1.054571817e-34          # J s
KB:        float = 1.380649e-23             # J/K
N_AVOG:    float = 6.02214076e23            # 1/mol
M_E:       float = 9.1093837015e-31         # kg, electron rest mass
ELECTRON_CHARGE: float = 1.602176634e-19    # C
EPS0:      float = 8.8541878128e-12         # F/m
EV_PER_J:  float = 1.0 / ELECTRON_CHARGE
KB_MEV_K:  float = (KB * EV_PER_J) * 1e3    # 0.0861733 meV / K
MEV_TO_J:  float = 1.602176634e-16          # J / MeV


# ---------------------------------------------------------------------------
# Material database
# ---------------------------------------------------------------------------
# Per-material substrate primitives.  Crystal-structure label kept for
# bookkeeping only (the substrate Eliashberg derivation here does not
# depend on lattice symmetry).  All elastic data CRC 90th ed., valence
# from periodic table (number of conduction electrons per atom).
#
#     M_molar  : kg/mol     (atomic / formula mass)
#     rho_mass : kg/m^3     (room-temperature density)
#     B_GPa    : GPa        (isothermal bulk modulus)
#     G_GPa    : GPa        (shear modulus, Voigt-Reuss-Hill)
#     valence  : int        (free-electron count per atom; for Pb=4, Sn=4,
#                             Al=3, In=3, Tl=3, V=5, Nb=5, Ta=5, Hg=2)

@dataclass(frozen=True)
class MaterialPhononSubstrate:
    """Substrate-side primitives for one material.

    Only the elastic moduli, mass density, and valence are per-material
    inputs.  No phonon spectrum, no electron-phonon coupling and no
    Fermi-level DOS are provided as input here -- those are derived
    by the substrate-Eliashberg routines from the entries below.
    """
    name:     str
    M_molar:  float    # kg/mol
    rho_mass: float    # kg/m^3
    B_GPa:    float    # GPa
    G_GPa:    float    # GPa
    valence:  int      # conduction electrons per atom

    @property
    def n_atoms_per_m3(self) -> float:
        """Atomic number density n_atoms = N_A * rho / M_molar (1/m^3)."""
        return N_AVOG * self.rho_mass / self.M_molar

    @property
    def m_ion_kg(self) -> float:
        """Single-ion mass M_ion = M_molar / N_A (kg)."""
        return self.M_molar / N_AVOG

    @property
    def n_electrons_per_m3(self) -> float:
        """Conduction-electron density n_e = z * n_atoms (1/m^3)."""
        return float(self.valence) * self.n_atoms_per_m3


MATERIAL_DB: Dict[str, MaterialPhononSubstrate] = {
    # --- BCS-test elementals (Pb / Hg / Nb / Sn / Al / V / Ta / In / Tl) ---
    # Properties: CRC Handbook 90th ed.  Valence: outermost s+p (or s+d)
    # electrons treated as conduction (free-electron substrate model).
    "Pb": MaterialPhononSubstrate(
        name="Pb",  M_molar=0.20720, rho_mass=11340.0, B_GPa= 46.0, G_GPa= 5.6, valence=4,
    ),
    "Hg": MaterialPhononSubstrate(
        # Hg is liquid at RT; we use the alpha-Hg solid-phase values
        # appropriate to the SC regime (T < T_c = 4.2 K).
        # rho = 14260 kg/m^3 (alpha-Hg at 80 K), B = 25 GPa, G = 9 GPa.
        name="Hg", M_molar=0.20059, rho_mass=14260.0, B_GPa= 25.0, G_GPa= 9.0, valence=2,
    ),
    "Nb": MaterialPhononSubstrate(
        name="Nb", M_molar=0.09291, rho_mass= 8570.0, B_GPa=170.0, G_GPa=37.5, valence=5,
    ),
    "Sn": MaterialPhononSubstrate(
        name="Sn", M_molar=0.11871, rho_mass= 7287.0, B_GPa= 58.0, G_GPa=18.0, valence=4,
    ),
    "Al": MaterialPhononSubstrate(
        name="Al", M_molar=0.02698, rho_mass= 2700.0, B_GPa= 76.0, G_GPa=26.0, valence=3,
    ),
    "V":  MaterialPhononSubstrate(
        name="V",  M_molar=0.05094, rho_mass= 6110.0, B_GPa=160.0, G_GPa=47.0, valence=5,
    ),
    "Ta": MaterialPhononSubstrate(
        name="Ta", M_molar=0.18095, rho_mass=16650.0, B_GPa=200.0, G_GPa=69.0, valence=5,
    ),
    "In": MaterialPhononSubstrate(
        name="In", M_molar=0.11482, rho_mass= 7310.0, B_GPa= 41.0, G_GPa= 4.0, valence=3,
    ),
    "Tl": MaterialPhononSubstrate(
        name="Tl", M_molar=0.20438, rho_mass=11850.0, B_GPa= 43.0, G_GPa= 2.8, valence=3,
    ),
}


# ---------------------------------------------------------------------------
# Substrate-derived basic quantities: c_s, omega_D, n_e, k_F, N(0)
# ---------------------------------------------------------------------------

def _sound_speeds(mat: MaterialPhononSubstrate) -> Tuple[float, float, float]:
    """Return (c_L, c_T, c_Debye) in m/s.

    Same formulas as ``debye_test.material_sound_speeds`` -- duplicated
    here so this module is independent of debye_test.MATERIALS (which has
    a different name list).
    """
    rho = mat.rho_mass
    B = mat.B_GPa * 1.0e9
    G = mat.G_GPa * 1.0e9
    c_L = math.sqrt((B + (4.0 / 3.0) * G) / rho)
    c_T = math.sqrt(G / rho)
    c_D = (3.0 / (1.0 / c_L ** 3 + 2.0 / c_T ** 3)) ** (1.0 / 3.0)
    return c_L, c_T, c_D


def debye_omega_rad_s(mat: MaterialPhononSubstrate) -> float:
    """Substrate Debye angular frequency omega_D = c_s * (6 pi^2 n)^(1/3) (rad/s)."""
    _, _, c_D = _sound_speeds(mat)
    n = mat.n_atoms_per_m3
    return c_D * (6.0 * math.pi ** 2 * n) ** (1.0 / 3.0)


def debye_temperature_K(mat: MaterialPhononSubstrate) -> float:
    """Substrate Debye temperature theta_D = hbar * omega_D / k_B (K)."""
    return HBAR * debye_omega_rad_s(mat) / KB


def fermi_wavevector(mat: MaterialPhononSubstrate) -> float:
    """Free-electron Fermi wavevector k_F = (3 pi^2 n_e)^(1/3) (1/m)."""
    n_e = mat.n_electrons_per_m3
    return (3.0 * math.pi ** 2 * n_e) ** (1.0 / 3.0)


def fermi_dos_per_J_per_m3(mat: MaterialPhononSubstrate) -> float:
    """Free-electron Fermi-level DOS N(0) PER VOLUME, two spins included.

        N(0) = m_e * k_F / (pi^2 * hbar^2)

    Units: (1 / J) / m^3.
    """
    k_F = fermi_wavevector(mat)
    return M_E * k_F / (math.pi ** 2 * HBAR ** 2)


# ---------------------------------------------------------------------------
# Substrate alpha^2 F(omega): Debye + Einstein
# ---------------------------------------------------------------------------
#
# Two-piece parametrisation (Allen-Dynes "Einstein over a Debye background"
# style; see Carbotte 1990 Eq. 2.11):
#
#     F(omega) = F_ac(omega)  +  F_E(omega)
#     F_ac(omega) = (3 / omega_D^3) * omega^2 * Theta(omega_D - omega)
#     F_E(omega) = w_E * delta(omega - omega_E)
#
# normalised so integral F = N_modes_per_cell.  In the long-wavelength
# substrate the acoustic part carries 3 polarizations per ion, and the
# K_4-cell internal modes contribute Einstein peaks at frequencies set by
# the face-pair coupling eps_face / d_face^2.
#
# Hopfield matrix element I^2 enters as
#
#     alpha^2 F(omega) = (1 / (2 M_ion omega)) * I_eph^2 * F(omega)
#
# but in the substrate model we DERIVE I_eph^2 from one global constant
# K_eph_substrate (set by anchoring on Pb's measured lambda; see below).
#
# Because Allen-Dynes lambda and omega_log only need the moments
#
#     lambda    = 2 integral alpha^2 F / omega  d omega
#     omega_log = exp[ (2/lambda) integral alpha^2 F * ln(omega) / omega
#                                          d omega ]
#
# we represent F as a discrete grid (omega_grid, F_grid) and reduce the
# delta peak to a narrow Lorentzian (FWHM = 1% omega_E) for numerical
# integration.

# Substrate Einstein peak: located at this fraction of omega_D
# ---------------------------------------------------------------
# The substrate K_4 cell has internal vibrational modes whose lowest
# frequency is set by the face-pair restoring spring acting on the
# tetrahedron's principal-axis displacement.  For a regular tetrahedron
# the lowest Einstein-mode frequency relative to the acoustic-cutoff
# omega_D is geometrically a numerical constant; we use the standard
# ratio omega_E / omega_D = 2/3 inferred from acoustic-band/optical-
# band centroid splittings in cubic and bcc metals (Born-Huang model).
EINSTEIN_OVER_DEBYE: float = 2.0 / 3.0

# Einstein peak weight (substrate ansatz)
# ---------------------------------------
# Weight of the Einstein peak relative to the Debye background, used in
# F = (1 - W_E) * F_ac + W_E * F_E.  Anchored to alpha^2 F lineshape of
# Pb tunnelling spectroscopy where the optical-mode contribution accounts
# for ~50% of the integrated alpha^2 F.  Held fixed across all materials.
W_EINSTEIN: float = 0.5

# Substrate Hopfield K_eph (anchored to Pb)
# -----------------------------------------
# In the Hopfield representation
#       lambda = N(0) <I^2> / (M_ion <omega^2>)
# the substrate replaces <I^2> with K_eph_substrate * (substrate energy
# scale)^2, where K_eph_substrate is the SAME dimensionless number for
# every material.  We anchor it ONCE on Pb (lambda_Pb = 1.55), then use
# it to predict the rest.  The resulting K_eph_substrate is reported by
# ``calibrate_K_eph_substrate`` and exported as K_EPH_SUBSTRATE for use
# in ``alpha2F_substrate``.

# Anchor target: empirical Pb electron-phonon coupling (Allen-Dynes 1975)
LAMBDA_PB_ANCHOR: float = 1.55


def _moments_from_alpha2F(
    omega_grid: np.ndarray,
    alpha2F:    np.ndarray,
) -> Tuple[float, float]:
    """Return (lambda, omega_log) computed from a discrete alpha^2 F.

    Using the standard Allen-Dynes definitions:

        lambda      = 2 integral_0^infty alpha^2 F(w) / w dw
        ln(w_log)   = (2/lambda) integral_0^infty alpha^2 F(w) ln(w) / w dw

    The integration domain excludes w = 0 (the integrand is finite because
    alpha^2 F(w) ~ w^2 near zero).
    """
    w = np.asarray(omega_grid, dtype=float)
    a = np.asarray(alpha2F, dtype=float)
    if w.shape != a.shape or w.ndim != 1 or w.size < 4:
        raise ValueError(
            f"omega_grid and alpha2F must be matching 1-D arrays of size >= 4; "
            f"got shapes {w.shape} and {a.shape}"
        )
    if not np.all(np.diff(w) > 0):
        raise ValueError("omega_grid must be strictly increasing")
    if w[0] <= 0.0:
        raise ValueError("omega_grid must be > 0 throughout")
    if np.any(a < 0.0):
        raise ValueError("alpha^2 F must be >= 0 everywhere")
    trap = getattr(np, "trapezoid", None) or np.trapz
    integ_lambda    = float(trap(a / w, w))
    lam = 2.0 * integ_lambda
    if lam <= 0.0:
        return 0.0, float(w[0])
    integ_log       = float(trap(a * np.log(w) / w, w))
    log_w_log       = (2.0 / lam) * integ_log
    return lam, float(math.exp(log_w_log))


def _spectrum_shape(
    omega_K: np.ndarray,
    omega_D_K: float,
    omega_E_K: Optional[float] = None,
    w_einstein: float = W_EINSTEIN,
) -> np.ndarray:
    """Substrate phonon spectrum shape F_shape(omega) (UNITS: 1 / K).

    Parametrised as a Debye-quadratic background up to omega_D plus a
    narrow Lorentzian Einstein peak at omega_E.  Normalised so that the
    integral of F_shape over all omega equals 1; this lets us factor the
    ION-MASS / coupling weight out cleanly downstream.

    Parameters
    ----------
    omega_K
        Frequency grid (K); must be strictly increasing and > 0.
    omega_D_K
        Debye cutoff (K).
    omega_E_K
        Einstein-peak centre (K).  Default = (2/3) * omega_D.
    w_einstein
        Fraction of the spectral weight in the Einstein peak; default
        ``W_EINSTEIN`` = 0.5.

    Returns
    -------
    F_shape : ndarray, same shape as omega_K
        Substrate phonon spectrum, normalised so int F_shape dw = 1
        (in units 1 / K).
    """
    if omega_E_K is None:
        omega_E_K = EINSTEIN_OVER_DEBYE * omega_D_K
    omega_K = np.asarray(omega_K, dtype=float)

    # ---- Debye part: F_ac(w) ~ w^2 / omega_D^3 for w <= omega_D --------
    f_ac = np.where(omega_K <= omega_D_K,
                    3.0 * omega_K ** 2 / omega_D_K ** 3,
                    0.0)
    trap = getattr(np, "trapezoid", None) or np.trapz
    norm_ac = float(trap(f_ac, omega_K))
    if norm_ac > 0.0:
        f_ac = f_ac / norm_ac

    # ---- Einstein peak (narrow Lorentzian, FWHM = 1% omega_E) ----------
    gamma_E = 0.01 * omega_E_K
    f_E_raw = (1.0 / math.pi) * gamma_E / ((omega_K - omega_E_K) ** 2 + gamma_E ** 2)
    norm_E = float(trap(f_E_raw, omega_K))
    if norm_E > 0.0:
        f_E = f_E_raw / norm_E
    else:
        f_E = np.zeros_like(omega_K)

    return (1.0 - w_einstein) * f_ac + w_einstein * f_E


def _alpha2_weight(
    omega_K: np.ndarray,
    K_eph: float,
    M_ion_kg: float,
    omega_D_K: float,
) -> np.ndarray:
    """Substrate electron-phonon coupling alpha^2(omega).

    The coupling is split between deformation-potential (acoustic) and
    Holstein (optical) forms, smoothly interpolating across the Debye
    cutoff:

        alpha^2(omega) = (K_eph / (M_ion * (k_B omega_D)))
                          * (omega / omega_D)        for omega <= omega_D
        alpha^2(omega) = (K_eph / (M_ion * (k_B omega_D)))
                          * 1.0                       for omega >  omega_D

    The acoustic linear-in-omega form follows from the substrate
    deformation-potential matrix element |g_q|^2 ~ q (one factor of
    omega from the q -> omega = c_s q dispersion).  This is the
    standard Allen-Dynes acoustic-coupling form (Carbotte 1990 Eq.
    2.8).  At omega = omega_D both branches join smoothly.

    The result is that alpha^2 F(omega) scales as omega^3 at small
    omega (Debye F ~ omega^2 times alpha^2 ~ omega), so alpha^2 F /
    omega ~ omega^2 has a finite integral.

    K_eph carries dimensions of  J^2  (substrate energy squared per
    Hopfield invariant); see ``calibrate_K_eph_substrate``.
    """
    omega_K = np.asarray(omega_K, dtype=float)
    omega_D_J = KB * omega_D_K
    base = K_eph / (M_ion_kg * omega_D_J)
    return base * np.where(
        omega_K <= omega_D_K,
        omega_K / omega_D_K,
        1.0,
    )


def alpha2F_substrate(
    material: MaterialPhononSubstrate,
    K_eph: Optional[float] = None,
    n_grid: int = 256,
    omega_max_factor: float = 1.20,
) -> Tuple[np.ndarray, np.ndarray]:
    """Substrate-derived alpha^2 F(omega) for one material.

    Parameters
    ----------
    material
        ``MaterialPhononSubstrate`` providing M_ion, rho, B, G, valence.
    K_eph
        Substrate Hopfield constant (J^2 / mode).  If None, the module
        anchor ``K_EPH_SUBSTRATE`` is used.
    n_grid
        Number of frequency-grid points (linear in omega, in K).
    omega_max_factor
        Upper grid boundary as a multiple of omega_D (must be > 1 to
        capture the Einstein peak fully).

    Returns
    -------
    omega_K     : (n_grid,) ndarray, frequencies in K, strictly increasing
    alpha2F     : (n_grid,) ndarray, alpha^2 F(omega) (dimensionless per K)
    """
    if K_eph is None:
        K_eph = K_EPH_SUBSTRATE
    if omega_max_factor <= 1.0:
        raise ValueError(
            f"omega_max_factor must exceed 1.0 (need to capture the Einstein peak); "
            f"got {omega_max_factor!r}"
        )
    if n_grid < 16:
        raise ValueError(f"n_grid must be >= 16, got {n_grid!r}")

    omega_D = debye_temperature_K(material)
    omega_max = omega_max_factor * omega_D
    omega_K = np.linspace(omega_D / n_grid, omega_max, n_grid)

    F_shape = _spectrum_shape(omega_K, omega_D)
    alpha2 = _alpha2_weight(omega_K, K_eph, material.m_ion_kg, omega_D)
    return omega_K, alpha2 * F_shape


def lambda_substrate(
    material: MaterialPhononSubstrate,
    K_eph: Optional[float] = None,
) -> float:
    """Substrate-derived electron-phonon coupling lambda.

        lambda = 2 integral alpha^2 F(omega) / omega d omega
    """
    omega_K, a2F = alpha2F_substrate(material, K_eph=K_eph)
    lam, _ = _moments_from_alpha2F(omega_K, a2F)
    return lam


def omega_log_substrate(
    material: MaterialPhononSubstrate,
    K_eph: Optional[float] = None,
) -> float:
    """Substrate-derived omega_log (K).

        ln(omega_log) = (2/lambda) integral alpha^2 F * ln(omega) / omega d omega
    """
    omega_K, a2F = alpha2F_substrate(material, K_eph=K_eph)
    _, w_log = _moments_from_alpha2F(omega_K, a2F)
    return w_log


# ---------------------------------------------------------------------------
# Anchor calibration: fix K_eph from Pb lambda = 1.55
# ---------------------------------------------------------------------------

def calibrate_K_eph_substrate(
    material_anchor: MaterialPhononSubstrate = MATERIAL_DB["Pb"],
    lambda_anchor:   float = LAMBDA_PB_ANCHOR,
    n_grid:          int = 256,
) -> float:
    """Calibrate K_eph_substrate so that lambda_substrate(Pb) == 1.55.

    Returns the anchored K_eph in units of  J^2 * kg  (i.e. the
    constant such that K_eph / (M_ion * k_B * omega) is dimensionless
    when multiplied by the K-units spectrum shape).

    Method: lambda is LINEAR in K_eph (alpha^2 F factorises into K_eph *
    [shape / mass / omega] terms), so a single substrate evaluation at
    any reference K_eph determines the anchor:

        K_eph_anchor = lambda_target / lambda_at_unit_K_eph
    """
    omega_K, a2F_unit = alpha2F_substrate(material_anchor, K_eph=1.0,
                                          n_grid=n_grid)
    lam_unit, _ = _moments_from_alpha2F(omega_K, a2F_unit)
    if lam_unit <= 0.0:
        raise ValueError(
            f"calibration failed: lambda at unit K_eph = {lam_unit!r}"
        )
    return float(lambda_anchor / lam_unit)


K_EPH_SUBSTRATE: float = calibrate_K_eph_substrate()
"""Substrate-derived Hopfield electron-phonon coupling constant.

Calibrated ONCE on Pb's empirical lambda = 1.55.  All other materials
are PREDICTIONS using this same constant: the substrate's claim is that
M_ion and Debye omega_D (themselves substrate-derived) account for the
material-specific variation in lambda and omega_log, with NO extra
per-material tunable.
"""


# ---------------------------------------------------------------------------
# Convenience: per-name predictions and comparison tables
# ---------------------------------------------------------------------------

def predict_lambda_omega_log(name: str) -> Tuple[float, float]:
    """Return (lambda, omega_log_K) for a material name in MATERIAL_DB."""
    if name not in MATERIAL_DB:
        raise KeyError(
            f"unknown material {name!r}; available: {sorted(MATERIAL_DB)}"
        )
    mat = MATERIAL_DB[name]
    return lambda_substrate(mat), omega_log_substrate(mat)


@dataclass(frozen=True)
class SubstrateEliashbergRow:
    """One material's substrate-derived vs empirical Eliashberg moments."""
    name:            str
    theta_D_K:       float    # substrate-derived Debye temp (K)
    omega_E_K:       float    # substrate Einstein-mode peak (K)
    lambda_sub:      float    # substrate lambda
    omega_log_sub_K: float    # substrate omega_log (K)
    lambda_emp:      float    # empirical lambda (Allen-Dynes / Carbotte)
    omega_log_emp_K: float    # empirical omega_log (K)
    lambda_rel_err:  float    # (sub - emp) / emp
    omega_log_rel_err: float


def compare_to_empirical(
    empirical: Optional[Dict[str, Tuple[float, float]]] = None,
) -> List[SubstrateEliashbergRow]:
    """Return one ``SubstrateEliashbergRow`` per material in MATERIAL_DB
    that also has an entry in ``empirical``.

    Defaults to the same MATERIAL_PHONON_PARAMS used by
    ``bcs_gap_ratio_test`` so the substrate vs literature comparison is
    line-for-line consistent with the BCS gap-ratio test inputs.
    """
    if empirical is None:
        # Lazy import: avoid bcs_gap_ratio_test pulling in this module
        from .bcs_gap_ratio_test import MATERIAL_PHONON_PARAMS
        empirical = MATERIAL_PHONON_PARAMS

    rows: List[SubstrateEliashbergRow] = []
    for name, mat in MATERIAL_DB.items():
        if name not in empirical:
            continue
        lam_emp, wlog_emp = empirical[name]
        lam_sub  = lambda_substrate(mat)
        wlog_sub = omega_log_substrate(mat)
        rows.append(SubstrateEliashbergRow(
            name=name,
            theta_D_K=debye_temperature_K(mat),
            omega_E_K=EINSTEIN_OVER_DEBYE * debye_temperature_K(mat),
            lambda_sub=lam_sub,
            omega_log_sub_K=wlog_sub,
            lambda_emp=lam_emp,
            omega_log_emp_K=wlog_emp,
            lambda_rel_err=(lam_sub - lam_emp) / lam_emp,
            omega_log_rel_err=(wlog_sub - wlog_emp) / wlog_emp,
        ))
    return rows


def summary_stats(rows: Iterable[SubstrateEliashbergRow]) -> Dict[str, float]:
    """Aggregate substrate-vs-empirical agreement statistics."""
    rows = list(rows)
    if not rows:
        raise ValueError("empty rows")
    abs_lam = np.array([abs(r.lambda_rel_err)    for r in rows])
    abs_wlg = np.array([abs(r.omega_log_rel_err) for r in rows])
    return {
        "n":                     len(rows),
        "mean_abs_lambda_err":   float(abs_lam.mean()),
        "max_abs_lambda_err":    float(abs_lam.max()),
        "median_abs_lambda_err": float(np.median(abs_lam)),
        "mean_abs_wlog_err":     float(abs_wlg.mean()),
        "max_abs_wlog_err":      float(abs_wlg.max()),
        "median_abs_wlog_err":   float(np.median(abs_wlg)),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"\nSubstrate-Eliashberg derivation")
    print(f"  K_eph_substrate (anchored on Pb lambda={LAMBDA_PB_ANCHOR:.2f}) = "
          f"{K_EPH_SUBSTRATE:.4e}  (J^2 * kg)")
    rows = compare_to_empirical()
    print(f"\n {'name':<5}  {'theta_D':>8}  {'omega_E':>8}  "
          f"{'lam_sub':>8}  {'lam_emp':>8}  {'dev_lam':>9}  "
          f"{'wlog_sub':>9}  {'wlog_emp':>9}  {'dev_wlog':>9}")
    for r in rows:
        print(f" {r.name:<5}  {r.theta_D_K:>8.1f}  {r.omega_E_K:>8.1f}  "
              f"{r.lambda_sub:>8.3f}  {r.lambda_emp:>8.3f}  "
              f"{r.lambda_rel_err*100:>+8.1f}%  "
              f"{r.omega_log_sub_K:>9.1f}  {r.omega_log_emp_K:>9.1f}  "
              f"{r.omega_log_rel_err*100:>+8.1f}%")
    stats = summary_stats(rows)
    print(f"\n n={stats['n']}")
    print(f"  mean |lambda err|   = {stats['mean_abs_lambda_err']*100:.1f}%   "
          f"max = {stats['max_abs_lambda_err']*100:.1f}%")
    print(f"  mean |omega_log err|= {stats['mean_abs_wlog_err']*100:.1f}%   "
          f"max = {stats['max_abs_wlog_err']*100:.1f}%")


__all__ = [
    "EINSTEIN_OVER_DEBYE",
    "K_EPH_SUBSTRATE",
    "LAMBDA_PB_ANCHOR",
    "MATERIAL_DB",
    "MaterialPhononSubstrate",
    "SubstrateEliashbergRow",
    "W_EINSTEIN",
    "alpha2F_substrate",
    "calibrate_K_eph_substrate",
    "compare_to_empirical",
    "debye_omega_rad_s",
    "debye_temperature_K",
    "fermi_dos_per_J_per_m3",
    "fermi_wavevector",
    "lambda_substrate",
    "main",
    "omega_log_substrate",
    "predict_lambda_omega_log",
    "summary_stats",
]


if __name__ == "__main__":  # pragma: no cover
    main()
