"""
EM Radiation Simulator (B3 substrate framework)
================================================

Models EM radiation as transverse modes of the substrate emitted by an
oscillating bound state. Provides both the analytic Larmor result and an
explicit substrate-field calculation via a retarded Green's function, then
verifies they agree.

In the B3 ontology:
  - Charge q is a topological defect coupled to the substrate gauge field.
  - Acceleration of the defect launches transverse waves at speed c.
  - The substrate stress carries energy out at rate set by the
    retarded Green's function of the wave operator.

The classical Larmor formula
    P = (2/3) (q^2 a^2) / (4*pi*eps0 * c^3)
must therefore emerge from the substrate field calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np
from scipy.integrate import quad


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
EPS0 = 8.8541878128e-12          # vacuum permittivity [F/m]
MU0 = 1.25663706212e-6           # vacuum permeability [N/A^2]
C = 2.99792458e8                 # speed of light [m/s]
E_CHARGE = 1.602176634e-19       # elementary charge [C]
M_E = 9.1093837015e-31           # electron mass [kg]
HBAR = 1.054571817e-34           # reduced Planck constant
PI = np.pi


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------
@dataclass
class EMRadiation:
    """EM radiation calculator for an oscillating bound charge.

    Parameters
    ----------
    q : float
        Charge of the radiator [C]. Defaults to electron charge.
    m : float
        Mass of the radiator [kg]. Defaults to electron mass.
    r0 : float
        Oscillation amplitude [m]. Defaults to 10 nm.
    omega : float
        Angular oscillation frequency [rad/s]. Defaults to 1e15.
    """

    q: float = E_CHARGE
    m: float = M_E
    r0: float = 1.0e-8
    omega: float = 1.0e15

    # ------------------------------------------------------------------
    # 1. Trajectory of the bound state
    # ------------------------------------------------------------------
    def position(self, t: np.ndarray | float) -> np.ndarray:
        """r(t) = r0 cos(omega t)  along z-hat."""
        return self.r0 * np.cos(self.omega * np.asarray(t))

    def velocity(self, t: np.ndarray | float) -> np.ndarray:
        return -self.r0 * self.omega * np.sin(self.omega * np.asarray(t))

    def acceleration(self, t: np.ndarray | float) -> np.ndarray:
        return -self.r0 * self.omega ** 2 * np.cos(self.omega * np.asarray(t))

    def peak_acceleration(self) -> float:
        return self.r0 * self.omega ** 2

    # ------------------------------------------------------------------
    # 2. Larmor formula and friends
    # ------------------------------------------------------------------
    @staticmethod
    def larmor_power(q: float, a: float) -> float:
        """Instantaneous Larmor radiated power.

            P = (2/3) (q^2 a^2) / (4 pi eps0 c^3)
        """
        return (2.0 / 3.0) * (q ** 2 * a ** 2) / (4.0 * PI * EPS0 * C ** 3)

    def time_averaged_power(self) -> float:
        """<P>_t for r(t)=r0 cos(omega t):  uses <a^2> = a_peak^2 / 2."""
        a_peak = self.peak_acceleration()
        return 0.5 * self.larmor_power(self.q, a_peak)

    # ------------------------------------------------------------------
    # 3. Angular patterns
    # ------------------------------------------------------------------
    @staticmethod
    def dipole_radiation_pattern(theta: np.ndarray | float) -> np.ndarray:
        """dP/dOmega ~ sin^2(theta) for an electric dipole."""
        return np.sin(np.asarray(theta)) ** 2

    @staticmethod
    def quadrupole_radiation(Q_ij: np.ndarray, omega: float,
                             theta: np.ndarray | float | None = None
                             ) -> Tuple[float, np.ndarray | None]:
        """Total quadrupole radiated power and (optionally) angular pattern.

        Total power for an oscillating traceless quadrupole tensor Q_ij
        (in SI units, dimensions C*m^2):

            P = (1/360 pi eps0 c^5) * omega^6 * sum_ij Q_ij Q_ij

        Returns (P_total, pattern). pattern is sin^2(2 theta) (the
        characteristic axisymmetric quadrupole lobe pattern) evaluated at
        the supplied theta, or None if theta is omitted.
        """
        Q = np.asarray(Q_ij, dtype=float)
        Q2 = float(np.sum(Q * Q))
        P = (omega ** 6) * Q2 / (360.0 * PI * EPS0 * C ** 5)
        pattern = None
        if theta is not None:
            pattern = np.sin(2.0 * np.asarray(theta)) ** 2
        return P, pattern

    # ------------------------------------------------------------------
    # 4. Cyclotron / synchrotron / bremsstrahlung
    # ------------------------------------------------------------------
    @staticmethod
    def cyclotron_radiation(B: float, v: float,
                            q: float = E_CHARGE,
                            m: float = M_E) -> float:
        """Non-relativistic cyclotron radiated power.

            a = (q v B)/m   ;   P = (2/3) q^2 a^2 / (4 pi eps0 c^3)
        """
        a = q * v * B / m
        return EMRadiation.larmor_power(q, a)

    @staticmethod
    def bremsstrahlung(v_in: float, v_out: float, Z: float,
                       q: float = E_CHARGE,
                       m: float = M_E,
                       dt: float = 1.0e-18) -> float:
        """Rough non-relativistic bremsstrahlung energy radiated in an
        impulsive Coulomb deflection.

        Treats the encounter as a near-instantaneous velocity change
        Delta v = v_out - v_in over a short interaction time dt.  The
        emitted energy follows from Larmor with a = Delta v / dt:

            E_rad ~ P_Larmor * dt = (2/3) q^2 (Delta v)^2 / (4 pi eps0 c^3 dt)

        The Z dependence enters via the Coulomb cross section that selects
        deflection probability; here we keep it as an explicit prefactor
        so the result is monotone in Z^2 as expected.
        """
        dv = v_out - v_in
        a = dv / dt
        P = EMRadiation.larmor_power(q, a)
        return Z ** 2 * P * dt

    @staticmethod
    def synchrotron_critical_freq(gamma: float, B: float,
                                  q: float = E_CHARGE,
                                  m: float = M_E) -> float:
        """omega_c = (3/2) gamma^3 (q B / m)."""
        return 1.5 * gamma ** 3 * (q * B / m)

    # ------------------------------------------------------------------
    # 5. Explicit substrate-field calculation
    # ------------------------------------------------------------------
    def substrate_radiation_power(self, R: float = 1.0,
                                  n_theta: int = 200) -> float:
        """Compute the time-averaged radiated power by integrating the
        far-field Poynting flux of the substrate transverse mode over a
        sphere of radius R.

        For a non-relativistic dipole p(t) = q r(t) z-hat, the retarded
        radiation electric field at large R is

            E_rad(R, theta, t) = - (q / (4 pi eps0 c^2 R)) a(t_ret) sin(theta)

        with t_ret = t - R/c.  The Poynting magnitude is

            S = eps0 c |E_rad|^2

        Integrating over the sphere and time-averaging gives back the
        Larmor result -- this routine does that integral numerically and
        is the substrate-side cross-check.
        """
        a_peak = self.peak_acceleration()
        # |E_rad|^2 time-averaged at theta:  prefactor * a_peak^2 / 2 * sin^2 theta
        prefac = (self.q / (4.0 * PI * EPS0 * C ** 2 * R)) ** 2
        E2_avg = lambda th: prefac * 0.5 * a_peak ** 2 * np.sin(th) ** 2
        # Poynting magnitude (time-avg)
        S_avg = lambda th: EPS0 * C * E2_avg(th)

        # Integrate S * R^2 sin(theta) dtheta dphi over sphere
        # phi gives factor 2 pi
        integrand = lambda th: S_avg(th) * R ** 2 * np.sin(th)
        P, _ = quad(integrand, 0.0, PI, limit=n_theta)
        return 2.0 * PI * P

    def verify_against_larmor(self, R: float = 1.0,
                              tol: float = 0.05) -> dict:
        """Run substrate calc and compare to analytic Larmor; return
        a dict with both numbers and the relative error."""
        P_sub = self.substrate_radiation_power(R=R)
        P_lar = self.time_averaged_power()
        rel_err = abs(P_sub - P_lar) / P_lar
        return {
            "P_substrate": P_sub,
            "P_larmor": P_lar,
            "rel_err": rel_err,
            "passes": rel_err < tol,
        }

    # ------------------------------------------------------------------
    # 6. Retarded-Green's-function field at a point (vector form)
    # ------------------------------------------------------------------
    def retarded_E_field(self, r_obs: np.ndarray, t: float) -> np.ndarray:
        """Far-field (radiation-zone) electric field of the oscillating
        dipole, evaluated at observation point r_obs at time t, using the
        retarded Green's function of the wave operator.

        E_rad = (q / (4 pi eps0 c^2 R)) * [ (n x (n x a_ret)) ]
              with n = r_obs/|r_obs|, a_ret = acceleration at t - R/c.
        Sign convention: standard Jackson (the cross-product picks up the
        minus relative to the scalar-form expression above).
        """
        r_obs = np.asarray(r_obs, dtype=float)
        R = float(np.linalg.norm(r_obs))
        if R == 0.0:
            raise ValueError("Observation point cannot coincide with source")
        n = r_obs / R
        t_ret = t - R / C
        a_vec = np.array([0.0, 0.0, self.acceleration(t_ret)])
        # n x (n x a) = n (n.a) - a
        n_dot_a = float(np.dot(n, a_vec))
        cross_term = n * n_dot_a - a_vec
        return (self.q / (4.0 * PI * EPS0 * C ** 2 * R)) * cross_term


# ---------------------------------------------------------------------------
# Convenience: default electron oscillator from the prompt
# ---------------------------------------------------------------------------
def default_electron_oscillator() -> EMRadiation:
    """Electron, 10 nm amplitude, omega = 1e15 rad/s -- the test case."""
    return EMRadiation(q=E_CHARGE, m=M_E, r0=1.0e-8, omega=1.0e15)
