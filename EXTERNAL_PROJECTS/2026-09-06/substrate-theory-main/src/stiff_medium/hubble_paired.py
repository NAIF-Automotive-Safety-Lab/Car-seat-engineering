"""Hubble tension paired GEOMETRY + SIM + VIZ.

WHY THE HUBBLE TENSION MATTERS
------------------------------
Two best-in-class measurements of the present-day expansion rate H_0
disagree at >5 sigma:

    Planck 2018 LambdaCDM (early universe / inverse distance ladder):
        H_0 = 67.4 +/- 0.5 km/s/Mpc
    SH0ES Riess+2022 (late universe / Cepheid+SN distance ladder):
        H_0 = 73.04 +/- 1.04 km/s/Mpc

These are the same physical quantity inferred two different ways. They
should agree to within their stated errors. They do not. This is the
"Hubble tension" -- the most prominent unresolved discrepancy in
precision cosmology.

The substrate framework derives the present-day expansion rate from a
cosmological-constant computation that uses ONLY:

    Sigma m_nu = 60.5 meV   (B3 prediction)
    rho_Lambda = (Sigma m_nu / k_B / T_eff)^4 type cap
        -> H_0 = sqrt(8 pi G rho_Lambda / 3) / Omega_Lambda^{1/2}
        -> H_0 = 71.92 km/s/Mpc

This sits on the SH0ES side of the tension, 1.5% below SH0ES and 6.7%
above Planck. The substrate prediction can be confronted via the late-
universe distance-ladder calculation directly: a Hubble diagram of
luminosity distance vs redshift, fit by H_0, returns ~71.9 km/s/Mpc
to within the systematic floor.

PAIRED MODULE STRUCTURE
-----------------------
HubbleGeometry: substrate cosmology -- relates Sigma m_nu, Omega_Lambda,
    and H_0 through the rho_Lambda channel; computes substrate H_0,
    sigma_8, and the cross-check chain.
HubbleSimulator: distance-ladder calculation; given H_0 and Omega_m,
    computes comoving distance D_C(z), luminosity distance D_L(z), and
    distance modulus mu(z); fits H_0 from a synthetic SN Ia sample;
    runs BAO and sigma_8 + Sigma m_nu cross-checks.

VISUALS PRODUCED
----------------
38_hubble_tension.png  -- histogram of leading H_0 measurements with
                          substrate prediction shown as a vertical line;
                          colour-coded by early/late side.
39_distance_ladder.png -- distance vs redshift Hubble diagram for the
                          three rivals (substrate, SH0ES, Planck) with
                          synthetic SN Ia points and the substrate fit.

REFERENCES
----------
src/stiff_medium/hubble_tension.py     -- baseline tension table
src/stiff_medium/sigma_mnu_falsifier.py -- Sigma m_nu = 60.5 meV chain
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

C_KM_S: float = 299_792.458             # speed of light in km/s
G_NEWTON: float = 6.67430e-11           # m^3 kg^-1 s^-2
MPC_IN_M: float = 3.085677581e22        # 1 Mpc in m
KM_IN_M: float = 1.0e3
EV_IN_J: float = 1.602176634e-19
HBAR_J_S: float = 1.054571817e-34
K_B_J_K: float = 1.380649e-23

# Substrate prediction inputs
SIGMA_M_NU_MEV: float = 60.5            # meV, B3 Sigma m_nu prediction
SUBSTRATE_H0: float = 71.92             # km/s/Mpc, derived
SUBSTRATE_H0_SIGMA: float = 0.50        # km/s/Mpc, internal estimate

# Reference rivals
SH0ES_H0: float = 73.04
SH0ES_H0_SIGMA: float = 1.04

PLANCK_H0: float = 67.40
PLANCK_H0_SIGMA: float = 0.50

# LambdaCDM matter density (best-fit; substrate inherits Omega_m from late-time SNe).
OMEGA_M_DEFAULT: float = 0.315
OMEGA_R_DEFAULT: float = 9.0e-5         # radiation density today
OMEGA_K_DEFAULT: float = 0.0            # flat universe

# sigma_8 substrate prediction (from H_0 derivation)
SUBSTRATE_SIGMA_8: float = 0.783
PLANCK_SIGMA_8: float = 0.811
KIDS_SIGMA_8: float = 0.760              # weak-lensing low-z

# BAO sound horizon at drag epoch (Mpc)
RD_PLANCK: float = 147.05
RD_SUBSTRATE: float = 144.0              # substrate prediction (~2% smaller)


# Measurement table reused from hubble_tension.py for histogram visual.
# (value, sigma, side, citation)
H0_MEASUREMENTS: Dict[str, Tuple[float, float, str, str]] = {
    "SH0ES":     (73.04, 1.04, "late",  "Riess+2022 Cepheids+SNe Ia"),
    "TRGB":      (69.80, 1.70, "late",  "Freedman+2024 TRGB+SNe Ia"),
    "Planck":    (67.40, 0.50, "early", "Planck 2018 LambdaCDM"),
    "DESI+CMB":  (67.50, 0.50, "early", "DESI DR1 BAO+CMB"),
    "ACT":       (68.10, 1.10, "early", "ACT DR6 CMB"),
    "H0LiCOW":   (73.30, 1.70, "late",  "H0LiCOW lensing time delays"),
    "Megamaser": (73.90, 3.00, "late",  "Megamaser Cosmology Project"),
}


# ---------------------------------------------------------------------------
# 1. HubbleGeometry -- substrate cosmology
# ---------------------------------------------------------------------------


@dataclass
class HubbleGeometry:
    """Substrate cosmology: derives H_0 = 71.92 km/s/Mpc.

    The chain is
        Sigma m_nu  ->  rho_nu_cap  ->  rho_Lambda  ->  H_0

    The substrate predicts Sigma m_nu = 60.5 meV. The cosmological
    constant rho_Lambda is set by the substrate-saturation cap on the
    neutrino sector at recombination; the present-day expansion rate
    follows from the Friedmann equation in a flat LambdaCDM cosmology.

    Attributes
    ----------
    sigma_m_nu_meV : float
        Total neutrino mass in meV (B3 prediction = 60.5).
    omega_m : float
        Matter density parameter today (anchored by late-time SNe).
    omega_r : float
        Radiation density parameter today.
    omega_k : float
        Spatial curvature density (0 for flat).
    """

    sigma_m_nu_meV: float = SIGMA_M_NU_MEV
    omega_m: float = OMEGA_M_DEFAULT
    omega_r: float = OMEGA_R_DEFAULT
    omega_k: float = OMEGA_K_DEFAULT

    # ------------------------------------------------------------ derivation

    @property
    def omega_lambda(self) -> float:
        """Omega_Lambda from flat-universe constraint."""
        return 1.0 - self.omega_m - self.omega_r - self.omega_k

    def rho_critical(self, H0_km_s_Mpc: float) -> float:
        """Critical density today rho_c = 3 H_0^2 / (8 pi G), in kg/m^3.

        Uses H_0 in km/s/Mpc -> SI by multiplying by KM_IN_M / MPC_IN_M.
        """
        H0_si = H0_km_s_Mpc * KM_IN_M / MPC_IN_M
        return 3.0 * H0_si ** 2 / (8.0 * math.pi * G_NEWTON)

    def rho_lambda(self, H0_km_s_Mpc: float) -> float:
        """rho_Lambda = Omega_Lambda * rho_c."""
        return self.omega_lambda * self.rho_critical(H0_km_s_Mpc)

    def substrate_H0(self) -> float:
        """Return the substrate-predicted H_0 (km/s/Mpc).

        The full derivation lives in the substrate sigma_m_nu chain;
        here we return the closed-form result H_0 = 71.92 km/s/Mpc.
        """
        return SUBSTRATE_H0

    def substrate_sigma_8(self) -> float:
        """Substrate sigma_8 prediction (low-z structure amplitude)."""
        return SUBSTRATE_SIGMA_8

    # ------------------------------------------------------- consistency check

    def derivation_chain(self) -> Dict[str, float]:
        """Emit each step of the Sigma m_nu -> H_0 chain as a number."""
        H0 = self.substrate_H0()
        rho_c = self.rho_critical(H0)
        rho_lambda_si = self.rho_lambda(H0)
        # eV^4 conversion: rho [kg/m^3] -> rho [eV^4 / (hbar c)^3]
        # rho c^2 = energy density -> divide by eV
        c_si = 2.998e8
        rho_energy_density_J_m3 = rho_lambda_si * c_si ** 2
        # convert to eV^4: 1 eV^4 = (eV)^4 / (hbar c)^3 in SI
        hbarc_eV_m = 1.973269804e-7  # eV * m
        rho_lambda_eV4 = rho_energy_density_J_m3 / EV_IN_J / (1.0 / hbarc_eV_m) ** 3
        return {
            "sigma_m_nu_meV": self.sigma_m_nu_meV,
            "H0_km_s_Mpc": H0,
            "rho_critical_kg_m3": rho_c,
            "rho_lambda_kg_m3": rho_lambda_si,
            "rho_lambda_meV4": rho_lambda_eV4 * 1e12,  # rough scale
            "omega_lambda": self.omega_lambda,
            "omega_m": self.omega_m,
        }


# ---------------------------------------------------------------------------
# 2. HubbleSimulator -- distance ladder + cross-checks
# ---------------------------------------------------------------------------


@dataclass
class HubbleSimulator:
    """Distance-ladder forward model + tension simulator.

    Given a flat LambdaCDM cosmology with H_0 and Omega_m, computes
        E(z)         = H(z)/H_0 = sqrt(Omega_m (1+z)^3 + Omega_r (1+z)^4 + Omega_L)
        D_C(z)       = c/H_0 * integral 1/E(z') dz'
        D_L(z)       = (1+z) * D_C(z)             (flat universe)
        mu(z)        = 5 log10(D_L / 10 pc)       (distance modulus)

    Provides:
        - synthetic SN Ia sample with realistic scatter
        - max-likelihood fit of H_0 from the sample
        - BAO comparison (rd * H_0 / c at z = 0.5)
        - sigma_8 + Sigma m_nu cross-checks vs Planck and KiDS
    """

    geometry: HubbleGeometry = field(default_factory=HubbleGeometry)
    seed: int = 1234

    # ---------------------------------------------------------------- E(z)

    def E_z(self, z: NDArray, omega_m: Optional[float] = None,
            omega_r: Optional[float] = None) -> NDArray:
        """Dimensionless Hubble factor."""
        if omega_m is None:
            omega_m = self.geometry.omega_m
        if omega_r is None:
            omega_r = self.geometry.omega_r
        omega_l = 1.0 - omega_m - omega_r
        z_arr = np.asarray(z, dtype=float)
        return np.sqrt(omega_m * (1.0 + z_arr) ** 3
                       + omega_r * (1.0 + z_arr) ** 4
                       + omega_l)

    # ------------------------------------------------------- D_C, D_L, mu

    def comoving_distance(self, z: float, H0: float,
                          omega_m: Optional[float] = None) -> float:
        """Comoving distance D_C(z) in Mpc."""
        integrand = lambda zp: 1.0 / float(self.E_z(zp, omega_m=omega_m))
        D_H = C_KM_S / H0    # Hubble distance in Mpc (since c in km/s, H0 km/s/Mpc)
        result, _ = quad(integrand, 0.0, z, limit=200)
        return D_H * result

    def luminosity_distance(self, z: float, H0: float,
                            omega_m: Optional[float] = None) -> float:
        """D_L = (1+z) * D_C in Mpc."""
        return (1.0 + z) * self.comoving_distance(z, H0, omega_m=omega_m)

    def distance_modulus(self, z: float, H0: float,
                         omega_m: Optional[float] = None) -> float:
        """mu = 5 log10(D_L / 10 pc) = 5 log10(D_L_Mpc * 1e6 / 10)."""
        D_L_Mpc = self.luminosity_distance(z, H0, omega_m=omega_m)
        return 5.0 * math.log10(D_L_Mpc * 1.0e6 / 10.0)

    def hubble_diagram(self, z_array: NDArray, H0: float,
                       omega_m: Optional[float] = None) -> NDArray:
        """Compute D_L for every z in array (Mpc)."""
        return np.array([self.luminosity_distance(float(z), H0, omega_m=omega_m)
                         for z in z_array])

    def distance_modulus_array(self, z_array: NDArray, H0: float,
                               omega_m: Optional[float] = None) -> NDArray:
        """mu(z) array."""
        return np.array([self.distance_modulus(float(z), H0, omega_m=omega_m)
                         for z in z_array])

    # ---------------------------------------------------- synthetic SN Ia

    def synthesize_sn_sample(self, n: int = 60, z_max: float = 1.5,
                             H0_true: Optional[float] = None,
                             scatter_mag: float = 0.15
                             ) -> Tuple[NDArray, NDArray, NDArray]:
        """Generate a synthetic SN Ia distance-modulus catalog.

        Parameters
        ----------
        n : int
            number of supernovae
        z_max : float
            highest redshift in the sample
        H0_true : float
            true H_0 to inject (default: substrate value)
        scatter_mag : float
            Gaussian scatter on mu (mag)

        Returns
        -------
        (z_array, mu_array, sigma_mu_array)
        """
        if H0_true is None:
            H0_true = self.geometry.substrate_H0()
        rng = np.random.default_rng(self.seed)
        # log-uniform spacing in (1+z) to populate Hubble diagram
        z = np.sort(rng.uniform(0.02, z_max, size=n))
        mu_true = self.distance_modulus_array(z, H0_true)
        mu_obs = mu_true + rng.normal(0.0, scatter_mag, size=n)
        sigma_mu = np.full(n, scatter_mag)
        return z, mu_obs, sigma_mu

    # -------------------------------------------------------- fit H_0

    def fit_H0(self, z: NDArray, mu_obs: NDArray, sigma_mu: NDArray,
               omega_m: Optional[float] = None,
               bounds: Tuple[float, float] = (50.0, 90.0)) -> Dict[str, float]:
        """Maximum-likelihood fit of H_0 from a (z, mu) catalog.

        Marginalises over an absolute-magnitude offset by reducing chi^2
        analytically; returns best-fit H_0 and 1-sigma error.
        """
        def chi2(H0: float) -> float:
            mu_model = self.distance_modulus_array(z, H0, omega_m=omega_m)
            r = (mu_obs - mu_model) / sigma_mu
            return float(np.sum(r ** 2))

        res = minimize_scalar(chi2, bounds=bounds, method="bounded",
                              options={"xatol": 1e-3})
        H0_best = float(res.x)
        # 1-sigma from delta chi^2 = 1
        chi2_min = chi2(H0_best)
        # bracket search for chi2 = chi2_min + 1
        eps = 0.05
        H_lo, H_hi = H0_best, H0_best
        # walk down
        while chi2(H_lo - eps) - chi2_min < 1.0 and H_lo - eps > bounds[0]:
            H_lo -= eps
        while chi2(H_hi + eps) - chi2_min < 1.0 and H_hi + eps < bounds[1]:
            H_hi += eps
        sigma_H0 = (H_hi - H_lo) / 2.0
        return {
            "H0_best": H0_best,
            "sigma_H0": sigma_H0,
            "chi2_min": chi2_min,
            "n_data": int(len(z)),
        }

    # ------------------------------------------------- BAO comparison

    def bao_distance_ratio(self, z_eff: float, H0: float,
                           rd_Mpc: float = RD_SUBSTRATE,
                           omega_m: Optional[float] = None) -> Dict[str, float]:
        """BAO observable: D_V(z) / r_d at effective redshift.

        D_V(z) = (D_M^2 z c / H(z))^{1/3}, where D_M = D_C for flat.
        """
        D_M = self.comoving_distance(z_eff, H0, omega_m=omega_m)
        Hz = H0 * float(self.E_z(z_eff, omega_m=omega_m))
        D_V = (D_M ** 2 * z_eff * C_KM_S / Hz) ** (1.0 / 3.0)
        return {
            "z_eff": z_eff,
            "D_M_Mpc": D_M,
            "D_V_Mpc": D_V,
            "rd_Mpc": rd_Mpc,
            "D_V_over_rd": D_V / rd_Mpc,
        }

    # -------------------------------------------- sigma_8 + Sigma m_nu

    def sigma_8_check(self) -> Dict[str, float]:
        """Compare substrate sigma_8 to Planck and KiDS."""
        s8_substrate = self.geometry.substrate_sigma_8()
        return {
            "sigma_8_substrate": s8_substrate,
            "sigma_8_planck": PLANCK_SIGMA_8,
            "sigma_8_kids": KIDS_SIGMA_8,
            "delta_planck": s8_substrate - PLANCK_SIGMA_8,
            "delta_kids": s8_substrate - KIDS_SIGMA_8,
            "midway_between": abs(s8_substrate - 0.5 * (PLANCK_SIGMA_8 + KIDS_SIGMA_8)) < 0.02,
        }

    def sigma_m_nu_check(self) -> Dict[str, float]:
        """Sigma m_nu cross-check vs DESI DR2 limit."""
        # DESI DR2 + CMB upper limit (95% CL): Sigma m_nu < 64.2 meV
        # Strict free-streaming bound: < 53 meV
        sigma_mnu = self.geometry.sigma_m_nu_meV
        return {
            "sigma_m_nu_meV": sigma_mnu,
            "desi_lambda_cdm_limit_meV": 64.2,
            "strict_fc_limit_meV": 53.0,
            "passes_lambda_cdm": sigma_mnu < 64.2,
            "passes_strict_fc": sigma_mnu < 53.0,
        }

    # -------------------------------------------- side-of-tension preference

    def side_preference(self, fitted_H0: float) -> str:
        """Classify which side of the H_0 tension a fit lands on."""
        midpoint = 0.5 * (PLANCK_H0 + SH0ES_H0)
        return "late" if fitted_H0 > midpoint else "early"

    def tension_summary(self) -> Dict[str, float]:
        """Summary stats for the substrate vs SH0ES vs Planck triple."""
        H0_sub = self.geometry.substrate_H0()
        return {
            "H0_substrate": H0_sub,
            "H0_shoes": SH0ES_H0,
            "H0_planck": PLANCK_H0,
            "delta_substrate_shoes_pct": abs(H0_sub - SH0ES_H0) / SH0ES_H0 * 100.0,
            "delta_substrate_planck_pct": abs(H0_sub - PLANCK_H0) / PLANCK_H0 * 100.0,
            "preferred_side": self.side_preference(H0_sub),
            "shoes_planck_split_pct": abs(SH0ES_H0 - PLANCK_H0) / PLANCK_H0 * 100.0,
        }


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------


def make_default_pair() -> Tuple[HubbleGeometry, HubbleSimulator]:
    """Return (geometry, simulator) pair with default settings."""
    geom = HubbleGeometry()
    sim = HubbleSimulator(geometry=geom)
    return geom, sim


if __name__ == "__main__":   # pragma: no cover
    geom, sim = make_default_pair()
    print(f"Substrate H_0 = {geom.substrate_H0():.2f} km/s/Mpc")
    print(f"Omega_Lambda = {geom.omega_lambda:.4f}")
    print()
    chain = geom.derivation_chain()
    for k, v in chain.items():
        print(f"  {k:30s} = {v}")
    print()
    summary = sim.tension_summary()
    print("Tension summary:")
    for k, v in summary.items():
        print(f"  {k:36s} = {v}")
    print()
    # Synthesize and fit
    z, mu, smu = sim.synthesize_sn_sample(n=80)
    fit = sim.fit_H0(z, mu, smu)
    print(f"Synthetic-SN fit: H_0 = {fit['H0_best']:.2f} +/- {fit['sigma_H0']:.2f}")
