"""Gravitational lensing simulator.

Substrate-framework picture: gravity is a substrate-strain-induced effective
refractive index for transverse waves. A point mass M at the origin warps the
substrate; light propagating across the strain field acquires a path-length
deflection equivalent to Einstein's prediction

    alpha = 4 G M / (b c^2)

where b is the impact parameter. This module collects the standard lensing
observables (Einstein radius, two-image strong-lens equation, microlensing
magnification, weak-lensing tangential shear, Refsdal time delay, and a
crude cluster giant-arc geometry) so they can be cross-checked against
canonical systems (1919 Sun deflection, Q0957+561, G2237+0305, LMC microlensing).

Conventions
-----------
- All inputs in SI unless stated otherwise.
- Distances are angular-diameter distances in metres.
- Angles returned in radians; helpers `arcsec`/`mas` convert.
- Mass M in kilograms (helper `M_solar` for solar masses).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy import constants as const

# ----- unit helpers ---------------------------------------------------------
G = const.G                      # 6.674e-11 m^3 kg^-1 s^-2
C = const.c                      # 2.998e8 m/s
M_SUN = 1.98892e30               # kg
PC = const.parsec                # 3.0857e16 m
KPC = 1e3 * PC
MPC = 1e6 * PC
GPC = 1e9 * PC
ARCSEC = np.pi / (180.0 * 3600.0)
MAS = 1e-3 * ARCSEC


def to_arcsec(theta_rad: float) -> float:
    return theta_rad / ARCSEC


def to_mas(theta_rad: float) -> float:
    return theta_rad / MAS


# ----- main class -----------------------------------------------------------
@dataclass
class GravitationalLensing:
    """Point-mass gravitational lensing in the substrate-strain picture."""

    G: float = G
    c: float = C

    # ------------------------------------------------------------------
    # 1. Deflection
    # ------------------------------------------------------------------
    def deflection_angle(self, M: float, b: float) -> float:
        """Einstein deflection angle (radians) for mass M, impact parameter b.

        alpha = 4 G M / (b c^2)
        """
        if b <= 0:
            raise ValueError("Impact parameter b must be positive.")
        return 4.0 * self.G * M / (b * self.c ** 2)

    # ------------------------------------------------------------------
    # 2. Einstein radius
    # ------------------------------------------------------------------
    def einstein_radius(
        self, M_lens: float, D_lens: float, D_source: float, D_LS: float
    ) -> float:
        """Angular Einstein radius (radians).

        theta_E = sqrt( 4 G M / c^2 * D_LS / (D_L * D_S) )
        """
        if D_lens <= 0 or D_source <= 0 or D_LS <= 0:
            raise ValueError("Distances must be positive.")
        return np.sqrt(
            4.0 * self.G * M_lens / self.c ** 2 * D_LS / (D_lens * D_source)
        )

    # ------------------------------------------------------------------
    # 3. Strong lensing — two image positions
    # ------------------------------------------------------------------
    def strong_lensing_images(
        self, beta: float, theta_E: float
    ) -> Tuple[float, float]:
        """Return (theta_plus, theta_minus) image angular positions.

        Lens equation for a point mass:  beta = theta - theta_E^2 / theta
            ->  theta^2 - beta*theta - theta_E^2 = 0
            ->  theta_pm = 0.5*(beta +/- sqrt(beta^2 + 4 theta_E^2))
        """
        disc = np.sqrt(beta ** 2 + 4.0 * theta_E ** 2)
        return 0.5 * (beta + disc), 0.5 * (beta - disc)

    def image_magnifications(
        self, beta: float, theta_E: float
    ) -> Tuple[float, float]:
        """Signed magnifications of the two images."""
        u = beta / theta_E
        # mu_pm = 1/4 * ( u/sqrt(u^2+4) + sqrt(u^2+4)/u  +/- 2 )
        s = np.sqrt(u ** 2 + 4.0)
        if u == 0:
            return np.inf, np.inf
        base = 0.25 * (u / s + s / u)
        return base + 0.5, base - 0.5

    # ------------------------------------------------------------------
    # 4. Microlensing magnification
    # ------------------------------------------------------------------
    def microlensing_magnification(self, u: float) -> float:
        """Paczynski point-source point-lens magnification.

        A(u) = (u^2 + 2) / (u * sqrt(u^2 + 4))
        """
        u = np.asarray(u, dtype=float)
        if np.any(u <= 0):
            raise ValueError("u must be positive.")
        return (u ** 2 + 2.0) / (u * np.sqrt(u ** 2 + 4.0))

    # ------------------------------------------------------------------
    # 5. Weak lensing — tangential shear of a point mass
    # ------------------------------------------------------------------
    def weak_lensing_shear(self, M: float, R: float, Sigma_crit: float | None = None,
                           D_lens: float | None = None,
                           D_source: float | None = None,
                           D_LS: float | None = None) -> float:
        """Tangential shear gamma_t at projected radius R from a point mass.

        For a point mass the convergence kappa = 0 outside the lens and
        gamma_t(R) = M / (pi * R^2 * Sigma_crit)
        equivalently in dimensionless form
        gamma_t(theta) = theta_E^2 / theta^2.

        If a critical surface density Sigma_crit is supplied directly it is
        used; otherwise the user can supply distances and we compute it as
        Sigma_crit = c^2 / (4 pi G) * D_S / (D_L D_LS).
        """
        if R <= 0:
            raise ValueError("R must be positive.")
        if Sigma_crit is None:
            if None in (D_lens, D_source, D_LS):
                # Use the simple "deflection-angle-over-radius" weak-field form,
                # gamma_t ~ 2 G M / (R c^2 * R) = (1/2) * alpha / R, but the
                # truly canonical answer for a point mass needs Sigma_crit.
                # Fall back to the dimensionless ratio assuming Sigma_crit
                # absorbs the geometry: returns (Schwarzschild radius)/(2R).
                rs = 2.0 * self.G * M / self.c ** 2
                return rs / (2.0 * R)
            Sigma_crit = (
                self.c ** 2 / (4.0 * np.pi * self.G)
                * D_source / (D_lens * D_LS)
            )
        return M / (np.pi * R ** 2 * Sigma_crit)

    # ------------------------------------------------------------------
    # 6. Refsdal time delay between two images
    # ------------------------------------------------------------------
    def time_delay(
        self,
        theta1: float,
        theta2: float,
        beta: float,
        theta_E: float,
        z_lens: float,
        H_0: float,
    ) -> float:
        """Refsdal time delay (seconds) between two images.

        Uses the point-mass result
            c * dt = (1+z_L) * D_dt * [ phi(theta1) - phi(theta2) ]
        with the Fermat potential
            phi(theta) = 0.5*(theta - beta)^2 - theta_E^2 ln|theta|.
        D_dt is the time-delay distance, approximated here as
            D_dt ~ c / H_0  (a low-z stand-in; for cosmological lenses the
        user should pass H_0 in SI units, s^-1).

        H_0 must be in 1/s (e.g. 70 km/s/Mpc -> 70e3 / MPC).
        """
        def phi(theta):
            return 0.5 * (theta - beta) ** 2 - theta_E ** 2 * np.log(abs(theta))

        D_dt = self.c / H_0
        dt = (1.0 + z_lens) * D_dt / self.c * (phi(theta2) - phi(theta1))
        return dt

    # ------------------------------------------------------------------
    # 7. Cluster giant-arc geometry
    # ------------------------------------------------------------------
    def cluster_lensing_arc(
        self, M_total: float, R_core: float, beta: float
    ) -> dict:
        """Crude giant-arc geometry for a cored isothermal-like cluster.

        We treat the cluster as a softened point mass: the deflection
        saturates inside R_core and is Einstein-like outside. The arc
        radius and (linear) tangential length are returned.
        """
        # Effective Einstein radius assuming D_LS/(D_L D_S) ~ 1/R_core
        theta_E = np.sqrt(4.0 * self.G * M_total / self.c ** 2 / R_core)
        # Two image solutions of the softened lens equation
        disc = np.sqrt(beta ** 2 + 4.0 * theta_E ** 2)
        theta_plus = 0.5 * (beta + disc)
        # Tangential magnification factor at theta_plus (point-mass form)
        denom = 1.0 - (theta_E / theta_plus) ** 2
        mu_t = np.inf if denom == 0 else 1.0 / denom
        arc_radius = theta_plus * R_core
        arc_length = abs(mu_t) * 2.0 * np.sqrt(theta_E ** 2 - beta ** 2 / 4.0
                                               if theta_E ** 2 > beta ** 2 / 4
                                               else 0.0) * R_core
        return {
            "theta_E_rad": theta_E,
            "theta_E_arcsec": to_arcsec(theta_E),
            "image_angle_rad": theta_plus,
            "arc_radius_m": arc_radius,
            "tangential_magnification": mu_t,
            "arc_length_m": arc_length,
        }

    # ------------------------------------------------------------------
    # convenience: Schwarzschild radius
    # ------------------------------------------------------------------
    def schwarzschild_radius(self, M: float) -> float:
        return 2.0 * self.G * M / self.c ** 2


# ----- canonical reference systems -----------------------------------------
def sun_grazing_deflection() -> float:
    """1919 Eddington test: light grazing the Sun's limb (arcsec)."""
    M_sun = 1.98892e30
    R_sun = 6.957e8
    lens = GravitationalLensing()
    return to_arcsec(lens.deflection_angle(M_sun, R_sun))


def lmc_microlensing_einstein_radius(M_lens_solar: float = 0.5,
                                     D_lens_kpc: float = 10.0,
                                     D_source_kpc: float = 50.0) -> float:
    """Einstein radius (mas) of a stellar-mass lens midway to the LMC."""
    lens = GravitationalLensing()
    D_L = D_lens_kpc * KPC
    D_S = D_source_kpc * KPC
    D_LS = D_S - D_L
    theta_E = lens.einstein_radius(M_lens_solar * M_SUN, D_L, D_S, D_LS)
    return to_mas(theta_E)
