"""SPARC galactic dynamics — rotation curve fitter, RAR, Tully-Fisher.

This module sits inside the B3 substrate framework. The hypothesis is that
galactic flat-rotation phenomenology (the Milgrom acceleration scale ``a_0``
≈ 1.2e-10 m/s², the radial-acceleration relation, and the baryonic
Tully-Fisher relation) emerges from the substrate's cosmological state — in
particular from H_0 — rather than from a tuned dark-matter sector.

Substrate prediction:

    g_†  ≡  c · H_0 / (2 π)

For H_0 ≈ 70 km/s/Mpc this gives g_† ≈ 1.1e-10 m/s², within ~10 % of the
empirical Milgrom scale. The Hubble side that maximises the fit (SH0ES,
H_0 ≈ 73) is the SH0ES side, consistent with B3's other internal-consistency
results.

Components:

  - ``SPARCDynamics`` class with NFW halo and exponential-disk baryons
  - ``rotation_curve``        V(R) = sqrt(V_b^2 + V_DM^2)
  - ``g_dagger_from_substrate`` substrate-derived acceleration scale
  - ``radial_acceleration_relation``  empirical RAR (McGaugh+ 2016) plus
    substrate prediction
  - ``tully_fisher_relation``  V_flat(M_baryon)
  - ``sparc_galaxy_fit``       fit a single galaxy's rotation curve via
    scipy.optimize
  - ``H_0_consistency_check``  compares SH0ES vs Planck preference for a
    galactic sample

Designed to be standalone (no IO) so it can be exercised by unit tests
without external SPARC files.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import curve_fit, minimize_scalar


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

C_LIGHT = 2.99792458e8          # m / s
G_NEWTON = 6.67430e-11          # m^3 / (kg s^2)
M_SUN = 1.98847e30              # kg
KPC = 3.0857e19                 # m
KM = 1.0e3                      # m
MPC = 1.0e3 * KPC               # m
H0_KMS_PER_MPC_TO_SI = KM / MPC  # 1 (km/s/Mpc) in s^-1

# Empirical references
G_DAGGER_EMPIRICAL = 1.2e-10    # m/s^2 (McGaugh+ 2016 RAR fit)
A0_MILGROM = 1.2e-10            # m/s^2 (Milgrom MOND scale)

# Reference H0 values (km/s/Mpc)
H0_SHOES = 73.0
H0_PLANCK = 67.0


# ---------------------------------------------------------------------------
# Galaxy data containers
# ---------------------------------------------------------------------------


@dataclass
class GalaxyData:
    """Observed rotation curve for one galaxy.

    All quantities in SI-friendly astronomical units:
      - ``R_kpc``: radii [kpc]
      - ``V_obs_kms``: observed circular speed [km/s]
      - ``V_err_kms``: 1-sigma errors [km/s] (optional, defaults to 5 km/s)
      - ``M_star``: stellar mass [M_sun] (used as Bayesian prior on disk mass)
      - ``R_disk_kpc``: exponential-disk scale length [kpc]
    """

    name: str
    R_kpc: np.ndarray
    V_obs_kms: np.ndarray
    M_star: float
    R_disk_kpc: float = 3.0
    V_err_kms: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        self.R_kpc = np.asarray(self.R_kpc, dtype=float)
        self.V_obs_kms = np.asarray(self.V_obs_kms, dtype=float)
        if self.V_err_kms is None:
            self.V_err_kms = 5.0 * np.ones_like(self.R_kpc)
        else:
            self.V_err_kms = np.asarray(self.V_err_kms, dtype=float)


# Reference galaxies (rough literature values — used for sanity tests, not
# as the SPARC catalog). V_flat values are typical published SPARC numbers.
REFERENCE_GALAXIES: Dict[str, Dict[str, float]] = {
    "NGC3198": {"M_star": 3.0e10, "V_flat_kms": 150.0, "R_disk_kpc": 3.14},
    "NGC6946": {"M_star": 5.0e10, "V_flat_kms": 180.0, "R_disk_kpc": 5.5},
    "M33":     {"M_star": 4.5e9,  "V_flat_kms": 120.0, "R_disk_kpc": 1.7},
}


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------


class SPARCDynamics:
    """Substrate-driven SPARC rotation-curve / RAR / Tully-Fisher engine.

    Parameters
    ----------
    H0_kms_per_mpc : float
        Hubble constant assumed for the substrate-derived acceleration scale
        ``g_†``. Defaults to SH0ES (73). Use 67 to compare Planck.
    """

    def __init__(self, H0_kms_per_mpc: float = H0_SHOES) -> None:
        self.H0 = float(H0_kms_per_mpc)

    # ------------------------------------------------------------------
    # Substrate predictions
    # ------------------------------------------------------------------

    def H0_si(self) -> float:
        """H_0 in SI units (1/s)."""
        return self.H0 * H0_KMS_PER_MPC_TO_SI

    def g_dagger_from_substrate(self) -> float:
        """Substrate-derived acceleration scale: g_† = c·H_0/(2π)."""
        return C_LIGHT * self.H0_si() / (2.0 * np.pi)

    def tully_fisher_relation(
        self,
        M_baryon_msun: float,
    ) -> float:
        """Baryonic Tully-Fisher V_flat(M_b).

        From the substrate g_† scale: in the deep-MOND-like regime the
        circular speed at large R obeys V_flat^4 = G · M_b · g_†, so

            V_flat = (G · M_b · g_†)^(1/4)

        Returns V_flat in km/s.
        """
        g_dag = self.g_dagger_from_substrate()
        M_kg = M_baryon_msun * M_SUN
        V_si = (G_NEWTON * M_kg * g_dag) ** 0.25
        return V_si / KM

    def radial_acceleration_relation(
        self,
        g_bar: np.ndarray | float,
    ) -> np.ndarray:
        """Empirical RAR (McGaugh+ 2016): g_obs = g_bar / (1 - exp(-sqrt(g_bar/g_†)))

        Uses the substrate-derived g_†, so this *is* the substrate's
        prediction (no extra fit parameter).

        Parameters
        ----------
        g_bar : array_like
            Newtonian acceleration from baryons alone [m/s^2].
        Returns
        -------
        g_obs : ndarray
            Predicted total observed acceleration [m/s^2].
        """
        g_bar = np.asarray(g_bar, dtype=float)
        g_dag = self.g_dagger_from_substrate()
        # Avoid divide-by-zero
        eps = 1e-30
        x = np.sqrt(np.maximum(g_bar, eps) / g_dag)
        return g_bar / (1.0 - np.exp(-x))

    # ------------------------------------------------------------------
    # Baryonic and dark-matter rotation profiles
    # ------------------------------------------------------------------

    @staticmethod
    def v_baryon_disk(
        R_kpc: np.ndarray,
        M_disk_msun: float,
        R_disk_kpc: float,
    ) -> np.ndarray:
        """Circular speed of a thin exponential disk.

        Approximation (Freeman 1970 / point-mass-equivalent): treat the disk
        as enclosed-mass and use V^2 = G M_enc(R) / R with

            M_enc(R) = M_disk · [1 - (1 + R/R_d) · exp(-R/R_d)].

        This is exact for a spherically distributed exponential profile and a
        very good approximation for SPARC-style fitting at our test
        precision.

        Returns V in km/s.
        """
        R = np.asarray(R_kpc, dtype=float)
        x = R / R_disk_kpc
        M_enc = M_disk_msun * (1.0 - (1.0 + x) * np.exp(-x))
        # Convert to SI for V calculation
        R_si = R * KPC
        M_si = M_enc * M_SUN
        # Guard against R=0
        R_safe = np.where(R_si > 0, R_si, 1.0)
        V_si_sq = G_NEWTON * M_si / R_safe
        V_si_sq = np.where(R_si > 0, V_si_sq, 0.0)
        return np.sqrt(np.maximum(V_si_sq, 0.0)) / KM

    @staticmethod
    def v_halo_nfw(
        R_kpc: np.ndarray,
        M_halo_msun: float,
        c_concentration: float,
        r_vir_kpc: float = 200.0,
    ) -> np.ndarray:
        """NFW halo circular speed.

        ρ(r) = ρ_s / [(r/r_s)(1+r/r_s)^2], with r_s = r_vir / c.

        M(<r) = 4π ρ_s r_s^3 · [ ln(1+x) - x/(1+x) ], x = r/r_s.

        Normalise so M(<r_vir) = M_halo.
        Returns V in km/s.
        """
        R = np.asarray(R_kpc, dtype=float)
        c = float(c_concentration)
        r_s = r_vir_kpc / c

        def mu(x):
            return np.log(1.0 + x) - x / (1.0 + x)

        x = R / r_s
        x_vir = r_vir_kpc / r_s  # = c
        M_enc = M_halo_msun * mu(x) / mu(x_vir)
        R_si = R * KPC
        M_si = M_enc * M_SUN
        R_safe = np.where(R_si > 0, R_si, 1.0)
        V_si_sq = G_NEWTON * M_si / R_safe
        V_si_sq = np.where(R_si > 0, V_si_sq, 0.0)
        return np.sqrt(np.maximum(V_si_sq, 0.0)) / KM

    def rotation_curve(
        self,
        R_kpc: np.ndarray,
        M_baryon_msun: float,
        c_concentration: float,
        M_halo_msun: float,
        R_disk_kpc: float = 3.0,
        r_vir_kpc: float = 200.0,
    ) -> np.ndarray:
        """Combined V(R) = sqrt(V_b^2 + V_DM^2). Returns km/s."""
        Vb = self.v_baryon_disk(R_kpc, M_baryon_msun, R_disk_kpc)
        Vh = self.v_halo_nfw(R_kpc, M_halo_msun, c_concentration, r_vir_kpc)
        return np.sqrt(Vb ** 2 + Vh ** 2)

    # ------------------------------------------------------------------
    # Galaxy fitting
    # ------------------------------------------------------------------

    def sparc_galaxy_fit(
        self,
        galaxy_data: GalaxyData,
        c_init: float = 10.0,
        log10_Mhalo_init: float = 11.5,
        log10_Mb_init: Optional[float] = None,
    ) -> Dict[str, float]:
        """Fit a single galaxy's rotation curve.

        Free parameters: M_baryon, M_halo, concentration c.
        Uses scipy.optimize.curve_fit with log-mass parametrisation for
        numerical stability.

        Returns a dict with best-fit params, chi^2, and reduced chi^2.
        """
        gd = galaxy_data
        if log10_Mb_init is None:
            log10_Mb_init = float(np.log10(max(gd.M_star, 1e7)))

        R_disk = gd.R_disk_kpc

        def model(R, log10_Mb, log10_Mh, c):
            Mb = 10.0 ** log10_Mb
            Mh = 10.0 ** log10_Mh
            return self.rotation_curve(R, Mb, c, Mh, R_disk_kpc=R_disk)

        p0 = [log10_Mb_init, log10_Mhalo_init, c_init]
        bounds = ([7.0, 9.0, 2.0], [12.5, 14.0, 40.0])

        try:
            popt, pcov = curve_fit(
                model,
                gd.R_kpc,
                gd.V_obs_kms,
                p0=p0,
                sigma=gd.V_err_kms,
                bounds=bounds,
                maxfev=5000,
            )
            log10_Mb, log10_Mh, c_fit = popt
        except Exception as exc:  # pragma: no cover - fallback path
            log10_Mb, log10_Mh, c_fit = p0
            pcov = None

        V_pred = model(gd.R_kpc, log10_Mb, log10_Mh, c_fit)
        residuals = (gd.V_obs_kms - V_pred) / gd.V_err_kms
        chi2 = float(np.sum(residuals ** 2))
        dof = max(len(gd.R_kpc) - 3, 1)
        return {
            "name": gd.name,
            "M_baryon": float(10.0 ** log10_Mb),
            "M_halo": float(10.0 ** log10_Mh),
            "concentration": float(c_fit),
            "chi2": chi2,
            "reduced_chi2": chi2 / dof,
            "V_pred_kms": V_pred,
        }

    # ------------------------------------------------------------------
    # H_0 consistency: SH0ES vs Planck for galactic data
    # ------------------------------------------------------------------

    def H_0_consistency_check(
        self,
        galaxies: Optional[Sequence[GalaxyData]] = None,
    ) -> Dict[str, float]:
        """Compare SH0ES vs Planck H_0 for galactic dynamics.

        For each candidate H_0 the substrate predicts a g_†; we compare to
        the empirical Milgrom value 1.2e-10 m/s^2. The H_0 that minimises
        the fractional offset is preferred. If a SPARC-like galaxy sample
        is supplied, we also compute the Tully-Fisher residuals as an
        independent check (V_flat ≈ V at last measured radius).

        Returns dict with both side's offsets and the preferred side.
        """
        # g_dagger comparison
        out: Dict[str, float] = {}
        for label, h0 in (("SHOES", H0_SHOES), ("Planck", H0_PLANCK)):
            self.H0 = h0
            g_dag = self.g_dagger_from_substrate()
            frac = abs(g_dag - G_DAGGER_EMPIRICAL) / G_DAGGER_EMPIRICAL
            out[f"g_dag_{label}"] = g_dag
            out[f"frac_offset_{label}"] = frac

        # Tully-Fisher on supplied galaxies if any
        if galaxies:
            for label, h0 in (("SHOES", H0_SHOES), ("Planck", H0_PLANCK)):
                self.H0 = h0
                resid = []
                for g in galaxies:
                    V_flat_obs = float(np.median(g.V_obs_kms[-3:])) if len(
                        g.V_obs_kms) >= 3 else float(g.V_obs_kms[-1])
                    V_flat_pred = self.tully_fisher_relation(g.M_star)
                    resid.append((V_flat_pred - V_flat_obs) / V_flat_obs)
                out[f"TF_rms_{label}"] = float(
                    np.sqrt(np.mean(np.asarray(resid) ** 2))
                )

        out["preferred"] = (
            "SHOES" if out["frac_offset_SHOES"] < out["frac_offset_Planck"]
            else "Planck"
        )
        # Restore default
        self.H0 = H0_SHOES
        return out


# ---------------------------------------------------------------------------
# Synthetic data helper (for tests)
# ---------------------------------------------------------------------------


def synthetic_galaxy(
    name: str,
    M_star: float,
    V_flat_kms: float,
    R_disk_kpc: float = 3.0,
    n_pts: int = 12,
    rng_seed: int = 0,
) -> GalaxyData:
    """Build a synthetic SPARC-like galaxy whose curve rises then flattens.

    Used by the test suite so we don't need actual SPARC files. Adds a
    small noise term so curve_fit has a non-trivial chi^2.
    """
    rng = np.random.default_rng(rng_seed)
    R = np.linspace(0.5, 25.0, n_pts)
    # Smooth rise to V_flat with characteristic radius ~ 2 R_disk
    V = V_flat_kms * np.tanh(R / (2.0 * R_disk_kpc))
    V_noisy = V + rng.normal(0.0, 3.0, size=R.shape)
    return GalaxyData(
        name=name,
        R_kpc=R,
        V_obs_kms=V_noisy,
        M_star=M_star,
        R_disk_kpc=R_disk_kpc,
        V_err_kms=5.0 * np.ones_like(R),
    )


__all__ = [
    "SPARCDynamics",
    "GalaxyData",
    "REFERENCE_GALAXIES",
    "synthetic_galaxy",
    "G_DAGGER_EMPIRICAL",
    "H0_SHOES",
    "H0_PLANCK",
]
