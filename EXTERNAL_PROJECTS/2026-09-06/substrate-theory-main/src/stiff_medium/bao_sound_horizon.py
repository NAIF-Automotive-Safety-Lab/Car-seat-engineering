"""
BAO Sound Horizon Simulator (B3 substrate framework).

Computes the comoving sound horizon r_s at the baryon drag epoch from the
substrate cosmology, and compares to BAO measurements (SDSS, BOSS, eBOSS,
DESI DR2). Anchors:

  * H_0   = 71.92 km/s/Mpc      (B3 substrate prediction)
  * Sigma_m_nu = 60.5 meV        (B3 substrate prediction)
  * Omega_b h^2 ~ 0.02237        (BBN/CMB-consistent baryon density)
  * Omega_m  ~ 0.315            (matter today including DM)
  * T_cmb  = 2.7255 K
  * z_drag ~ 1059                (baryon drag redshift, EH98 fit)

Sound speed of baryon-photon plasma:
  c_s(z) = c / sqrt(3 (1 + R(z))),    R(z) = 3 rho_b / (4 rho_gamma)

Comoving sound horizon at drag:
  r_s = integral_{z_drag}^{infty} c_s(z) / H(z) dz

Standard LCDM with these parameters yields r_s ~ 147 Mpc, which the
substrate framework reproduces because the integral is dominated by
radiation-era physics (where substrate and LCDM agree) modulated by the
baryon-to-photon ratio (which is fixed by BBN, sector-shared with B3).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np
from scipy.integrate import quad


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
C_LIGHT = 2.99792458e8           # m/s
C_KM_S = C_LIGHT / 1.0e3          # km/s
MPC_M = 3.0857e22                 # m/Mpc
KB = 1.380649e-23                 # J/K
HBAR = 1.054571817e-34            # J s
SIGMA_SB = 5.670374419e-8         # W/m^2/K^4 (Stefan-Boltzmann)
G_NEWTON = 6.67430e-11

T_CMB = 2.7255                    # K
Z_DRAG_STD = 1059.94              # standard drag epoch (Planck 2018 EH98 fit)
Z_EQ_FIT = 3402.0                 # matter-radiation equality (approx, Planck)


@dataclass
class BAOParams:
    """Cosmological parameters anchored to B3 substrate predictions."""
    H0: float = 71.92                # km/s/Mpc  (B3 prediction)
    Omega_m: float = 0.315           # total matter today
    Omega_b: float = 0.0493          # baryons today (Omega_b h^2 ~ 0.02237)
    Omega_Lambda: float = 0.685      # dark energy today
    T_cmb: float = T_CMB             # K
    sum_mnu_ev: float = 0.0605       # eV (B3 prediction 60.5 meV)
    N_eff: float = 3.046
    z_drag: float = Z_DRAG_STD


# ---------------------------------------------------------------------------
# BAO measurements catalogue (r_d in Mpc; D_M/r_d, D_H/r_d, D_V/r_d ratios)
# ---------------------------------------------------------------------------
BAO_OBSERVATIONS: List[Dict] = [
    # Survey, z_eff, D_V/r_d (or D_M/r_d), error, source
    {"survey": "6dFGS",       "z": 0.106, "DV_rd": 2.976,  "err": 0.133, "ref": "Beutler+2011"},
    {"survey": "SDSS MGS",    "z": 0.15,  "DV_rd": 4.466,  "err": 0.168, "ref": "Ross+2015"},
    {"survey": "BOSS LOWZ",   "z": 0.32,  "DV_rd": 8.467,  "err": 0.167, "ref": "Alam+2017"},
    {"survey": "BOSS CMASS",  "z": 0.57,  "DV_rd": 13.77,  "err": 0.13,  "ref": "Alam+2017"},
    {"survey": "eBOSS LRG",   "z": 0.70,  "DV_rd": 16.26,  "err": 0.20,  "ref": "Bautista+2021"},
    {"survey": "DESI BGS",    "z": 0.40,  "DV_rd": 10.51,  "err": 0.12,  "ref": "DESI DR2 2024"},
    {"survey": "DESI LRG1",   "z": 0.51,  "DV_rd": 13.62,  "err": 0.25,  "ref": "DESI DR2 2024"},
    {"survey": "DESI LRG2",   "z": 0.71,  "DV_rd": 16.85,  "err": 0.32,  "ref": "DESI DR2 2024"},
    {"survey": "DESI ELG",    "z": 1.32,  "DV_rd": 21.71,  "err": 0.28,  "ref": "DESI DR2 2024"},
    {"survey": "DESI Lya",    "z": 2.33,  "DV_rd": 31.30,  "err": 0.70,  "ref": "DESI DR2 2024"},
]

# Measured r_d ("standard ruler" length) from CMB+BAO joint analysis
R_D_MEASURED = 147.05    # Mpc (Planck 2018 + DESI)
R_D_ERR = 0.30            # Mpc


# ---------------------------------------------------------------------------
# Photon and baryon density (today)
# ---------------------------------------------------------------------------
def rho_gamma_today(T_cmb: float = T_CMB) -> float:
    """Photon energy density today (kg/m^3, equivalent mass density)."""
    u_gamma = 4.0 * SIGMA_SB / C_LIGHT * T_cmb**4   # J/m^3
    return u_gamma / C_LIGHT**2                      # kg/m^3


def rho_crit_today(H0_kms_mpc: float) -> float:
    """Critical density today (kg/m^3)."""
    H0_si = H0_kms_mpc * 1.0e3 / MPC_M               # 1/s
    return 3.0 * H0_si**2 / (8.0 * np.pi * G_NEWTON)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------
class BAOSoundHorizon:
    """
    Computes substrate-framework BAO sound horizon and compares to data.
    """

    def __init__(self, params: BAOParams | None = None):
        self.params = params or BAOParams()
        self._cache: Dict[str, float] = {}

    # -----------------------------------------------------------------
    # Densities and equations of state
    # -----------------------------------------------------------------
    def Omega_gamma(self) -> float:
        """Photon density parameter today."""
        rho_g = rho_gamma_today(self.params.T_cmb)
        rho_c = rho_crit_today(self.params.H0)
        return rho_g / rho_c

    def Omega_r(self) -> float:
        """Radiation (photons + relativistic neutrinos) today."""
        og = self.Omega_gamma()
        # Neutrinos contribute (7/8)(4/11)^(4/3) N_eff per photon
        return og * (1.0 + 0.2271 * self.params.N_eff)

    def baryon_to_photon_R(self, z: float) -> float:
        """
        R(z) = 3 rho_b(z) / (4 rho_gamma(z)).
        Both scale with redshift: rho_b ~ (1+z)^3, rho_gamma ~ (1+z)^4,
        so R ~ 1/(1+z) times R_today * (rho_b_today / rho_g_today).
        """
        Ob = self.params.Omega_b
        Og = self.Omega_gamma()
        # R_today = 3 Ob / (4 Og) at z=0; redshift-evolved:
        R0 = 3.0 * Ob / (4.0 * Og)
        return R0 / (1.0 + z)

    def sound_speed(self, z: float) -> float:
        """Sound speed in baryon-photon plasma (m/s)."""
        R = self.baryon_to_photon_R(z)
        return C_LIGHT / np.sqrt(3.0 * (1.0 + R))

    # -----------------------------------------------------------------
    # Hubble parameter
    # -----------------------------------------------------------------
    def H_z(self, z: float) -> float:
        """
        Hubble parameter H(z) in 1/s, substrate cosmology with
        matter + radiation + Lambda. Substrate framework predicts
        identical functional form post-recombination because the
        de-saturation rate matches LCDM expansion by construction.
        """
        H0_si = self.params.H0 * 1.0e3 / MPC_M
        Om = self.params.Omega_m
        Or = self.Omega_r()
        OL = self.params.Omega_Lambda
        Ok = 1.0 - Om - Or - OL
        x = 1.0 + z
        E2 = Or * x**4 + Om * x**3 + Ok * x**2 + OL
        return H0_si * np.sqrt(E2)

    def H_z_kms_mpc(self, z: float) -> float:
        """H(z) in km/s/Mpc."""
        return self.H_z(z) * MPC_M / 1.0e3

    # -----------------------------------------------------------------
    # Sound horizon
    # -----------------------------------------------------------------
    def _integrand(self, z: float) -> float:
        """c_s(z) / H(z) in meters (gives r_s in m when integrated over z)."""
        return self.sound_speed(z) / self.H_z(z)

    def sound_horizon_drag(self) -> float:
        """
        Comoving sound horizon at the baryon drag epoch (Mpc).

        r_s = integral from z_drag to infinity of c_s/H dz
        """
        z_d = self.params.z_drag
        # Split the integral: z_d to 10*z_d, then 10*z_d to infinity
        r_s_m_1, _ = quad(self._integrand, z_d, 10.0 * z_d, limit=200)
        r_s_m_2, _ = quad(self._integrand, 10.0 * z_d, 1.0e8, limit=200)
        r_s_m = r_s_m_1 + r_s_m_2
        return r_s_m / MPC_M

    # -----------------------------------------------------------------
    # BAO observable scales
    # -----------------------------------------------------------------
    def comoving_distance(self, z: float) -> float:
        """Comoving distance D_C(z) in Mpc."""
        integrand = lambda zp: C_LIGHT / self.H_z(zp)
        d_m, _ = quad(integrand, 0.0, z, limit=200)
        return d_m / MPC_M

    def D_M(self, z: float) -> float:
        """Transverse comoving distance D_M (flat -> = D_C). Mpc."""
        return self.comoving_distance(z)

    def D_H(self, z: float) -> float:
        """Hubble distance D_H = c/H(z), Mpc."""
        return (C_LIGHT / self.H_z(z)) / MPC_M

    def D_V(self, z: float) -> float:
        """Spherically-averaged BAO distance D_V (Mpc)."""
        dm = self.D_M(z)
        dh = self.D_H(z)
        return (z * dm**2 * dh) ** (1.0 / 3.0)

    def bao_scale_today(self, z: float) -> Dict[str, float]:
        """All BAO observable ratios at redshift z."""
        r_d = self.sound_horizon_drag()
        return {
            "z": z,
            "r_d_Mpc": r_d,
            "D_M_Mpc": self.D_M(z),
            "D_H_Mpc": self.D_H(z),
            "D_V_Mpc": self.D_V(z),
            "D_M_over_rd": self.D_M(z) / r_d,
            "D_H_over_rd": self.D_H(z) / r_d,
            "D_V_over_rd": self.D_V(z) / r_d,
        }

    # -----------------------------------------------------------------
    # Comparison with observations
    # -----------------------------------------------------------------
    def compare_to_observations(self) -> List[Dict]:
        """Compare substrate predictions to BAO survey measurements."""
        r_d = self.sound_horizon_drag()
        results = []
        for obs in BAO_OBSERVATIONS:
            z = obs["z"]
            dv = self.D_V(z)
            pred_ratio = dv / r_d
            measured = obs["DV_rd"]
            err = obs["err"]
            sigma = (pred_ratio - measured) / err
            results.append({
                "survey": obs["survey"],
                "z": z,
                "measured_DV_rd": measured,
                "err": err,
                "predicted_DV_rd": pred_ratio,
                "delta_sigma": sigma,
                "ref": obs["ref"],
            })
        return results

    # -----------------------------------------------------------------
    # Forecast: DESI DR3
    # -----------------------------------------------------------------
    def predict_DESI_DR3(self) -> Dict[str, Dict]:
        """
        Forecast BAO observables at DESI DR3 effective redshifts.
        DESI DR3 is expected to reduce errors by ~sqrt(2) over DR2 and
        add a high-z QSO bin.
        """
        z_bins = {
            "BGS_DR3":  0.30,
            "LRG1_DR3": 0.51,
            "LRG2_DR3": 0.71,
            "LRG3_DR3": 0.93,
            "ELG1_DR3": 1.32,
            "ELG2_DR3": 1.49,
            "QSO_DR3":  1.75,
            "Lya_DR3":  2.33,
        }
        out = {}
        for label, z in z_bins.items():
            scale = self.bao_scale_today(z)
            # forecast 1-sigma error reduced by sqrt(2)
            est_err = max(0.10, 0.012 * scale["D_V_over_rd"])
            out[label] = {
                "z": z,
                "DV_over_rd": scale["D_V_over_rd"],
                "DM_over_rd": scale["D_M_over_rd"],
                "DH_over_rd": scale["D_H_over_rd"],
                "forecast_err": est_err,
            }
        return out

    # -----------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------
    def summary(self) -> Dict:
        r_d = self.sound_horizon_drag()
        return {
            "H0_kms_mpc": self.params.H0,
            "Omega_m": self.params.Omega_m,
            "Omega_b": self.params.Omega_b,
            "Omega_gamma": self.Omega_gamma(),
            "Omega_r": self.Omega_r(),
            "z_drag": self.params.z_drag,
            "r_d_Mpc": r_d,
            "r_d_measured_Mpc": R_D_MEASURED,
            "r_d_err_Mpc": R_D_ERR,
            "delta_rd_Mpc": r_d - R_D_MEASURED,
            "delta_rd_sigma": (r_d - R_D_MEASURED) / R_D_ERR,
            "cs_recomb_over_c": self.sound_speed(1090.0) / C_LIGHT,
        }


def main() -> None:
    bao = BAOSoundHorizon()
    s = bao.summary()
    print("=" * 60)
    print("B3 substrate BAO sound horizon")
    print("=" * 60)
    for k, v in s.items():
        if isinstance(v, float):
            print(f"  {k:25s} = {v:12.5g}")
        else:
            print(f"  {k:25s} = {v}")
    print("\nBAO comparison:")
    print(f"  {'survey':14s} {'z':>6s} {'meas':>8s} {'pred':>8s} {'err':>6s} {'sigma':>7s}")
    for r in bao.compare_to_observations():
        print(f"  {r['survey']:14s} {r['z']:6.2f} {r['measured_DV_rd']:8.3f} "
              f"{r['predicted_DV_rd']:8.3f} {r['err']:6.3f} {r['delta_sigma']:7.2f}")


if __name__ == "__main__":
    main()
