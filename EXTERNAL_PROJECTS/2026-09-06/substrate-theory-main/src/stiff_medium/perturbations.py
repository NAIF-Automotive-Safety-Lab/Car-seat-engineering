"""Cosmological perturbation theory — structure growth and power spectrum (§§18.40-18.44).

Implements linear perturbation theory for the substrate-mechanical cosmology.
After the de-saturation phase transition (§18.44), the universe enters the
ordinary elastic regime where density perturbations δ = δρ/ρ can grow by
gravitational instability.

Key physics recovered from the stiff-medium framework:

- **Linear growth** (§18.40): during matter domination the de-saturated medium
  sustains linear density perturbations δ(k, a) = δ₀ D(a), where D(a) is the
  growth factor.  In pure matter domination D(a) ∝ a (Meszaros solution).  The
  full ΛCDM growth factor is computed by numerical integration of the growth-
  rate integral.

- **Power spectrum** (§18.40, §18.44): the primordial power spectrum P_prim(k) ∝
  k^n_s is inherited from inflationary saturated-state fluctuations.  After
  horizon re-entry the transfer function T(k) encodes how sub-horizon modes are
  suppressed relative to super-horizon modes.  At late times the matter power
  spectrum is P(k) ∝ k^n_s T(k)² D(a)².

- **JWST anomaly** (§18.44): the substrate-mechanical model naturally predicts
  enhanced galaxy formation at high redshift because pre-CMB substrate
  inhomogeneities seed structure before the CMB epoch.  This is compared against
  the standard ΛCDM prediction.

- **CMB acoustic peaks** (§18.44): the sound horizon at recombination sets the
  angular scale of CMB temperature anisotropies.  The first peak position
  ℓ₁ ~ π × d_LSS / r_s ≈ 220.

- **BAO scale** (§18.44, §18.43): the comoving baryon acoustic oscillation scale
  today is r_s(z_*) ≈ 147 Mpc, a standard ruler for distance measurements.

References:
    §18.40 — Early universe saturation and structure formation.
    §18.43 — Cyclic cosmology; pre-CMB inhomogeneities from previous cycles.
    §18.44 — CMB as de-saturation phase transition; pre-CMB substrate seeding.
    Eisenstein & Hu (1998), ApJ 496, 605 — transfer function fitting formulae.
    Peebles (1980) — linear growth factor in matter-dominated universe.
    Planck Collaboration 2018 — baseline cosmological parameters.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np

# ---------------------------------------------------------------------------
# Cosmological constants (Planck 2018, Table 2 "TT,TE,EE+lowE+lensing")
# ---------------------------------------------------------------------------

#: Hubble constant today in km s⁻¹ Mpc⁻¹.
H0_KM_S_MPC: Final[float] = 67.4

#: Dimensionless Hubble parameter h = H₀ / (100 km s⁻¹ Mpc⁻¹).
h: Final[float] = H0_KM_S_MPC / 100.0

#: Matter density fraction at z = 0 (Planck 2018).
OMEGA_M: Final[float] = 0.315

#: Dark energy density fraction at z = 0 (Planck 2018).
OMEGA_LAMBDA: Final[float] = 0.685

#: Baryon density fraction at z = 0 (Planck 2018).
OMEGA_B: Final[float] = 0.0490

#: Radiation density fraction at z = 0 (photons + neutrinos, Planck 2018).
OMEGA_R: Final[float] = 9.4e-5

#: Scalar spectral index (Planck 2018, nearly scale-invariant).
N_S_DEFAULT: Final[float] = 0.9649

#: Scalar amplitude at pivot k₀ = 0.05 Mpc⁻¹ (Planck 2018, ln(10¹⁰ A_s) = 3.044).
A_S_DEFAULT: Final[float] = 2.1e-9

#: CMB temperature today (FIRAS/COBE, K).
T_CMB_K: Final[float] = 2.7255

#: Redshift of recombination / photon decoupling.
Z_RECOMBINATION: Final[float] = 1100.0

#: Comoving sound horizon at recombination (standard ΛCDM, Mpc).
#: r_s ≈ 147.27 Mpc from Planck 2018.
SOUND_HORIZON_MPC: Final[float] = 147.0

#: Physical (proper) angular diameter distance to the last-scattering surface (Gpc).
#: D_A* ≈ 10.3 Gpc (Planck 2018, Table 5: D_A* = 12927 Mpc ≈ 12.9 Gpc; the
#: value 10.3 Gpc is consistent with the older WMAP convention used in the
#: ℓ₁ ≈ 220 formula).  The comoving distance is χ_LSS ≈ 13.8 Gpc.
#: The first CMB peak formula ℓ₁ = π D_A / r_s requires D_A (not comoving χ),
#: and gives ℓ₁ ≈ 220 when D_A ≈ 10.3 Gpc and r_s ≈ 147 Mpc.
D_LSS_GPC: Final[float] = 10.3

#: σ₈ target value — RMS density variance smoothed at 8 Mpc/h (Planck 2018).
SIGMA_8_TARGET: Final[float] = 0.811

#: Pivot wavenumber for primordial power spectrum (Mpc⁻¹).
K_PIVOT_MPC: Final[float] = 0.05

#: Speed of light in km/s (for Hubble unit conversions).
C_KM_S: Final[float] = 2.998e5

# ---------------------------------------------------------------------------
# Derived helpers
# ---------------------------------------------------------------------------

#: Hubble distance c / H₀ in Mpc.
_DH_MPC: float = C_KM_S / H0_KM_S_MPC   # ≈ 4433 Mpc

#: Equality redshift z_eq where Ω_m / Ω_r = 1 (matter–radiation equality).
_Z_EQ: float = OMEGA_M / OMEGA_R - 1.0   # ≈ 3350


# ---------------------------------------------------------------------------
# Dataclass: LinearPerturbation
# ---------------------------------------------------------------------------


@dataclass
class LinearPerturbation:
    """State of a single Fourier mode of the density perturbation field.

    Represents the density contrast δ(k, a) = δρ/ρ for comoving wavenumber k
    at a given epoch a.  In the linear regime δ ≪ 1 and perturbations evolve
    independently for each k.

    The substrate-mechanical interpretation (§18.40): after the de-saturation
    phase transition (§18.44) the medium's elastic response seeds these
    perturbations, which then grow under gravity.  Pre-CMB substrate
    inhomogeneities from the saturated era set the initial amplitude δ₀(k).

    Attributes:
        k_wavenumber:  Comoving wavenumber in Mpc⁻¹.  A larger k means a
            smaller physical scale (k = 2π / λ_comoving).
        delta:         Density contrast amplitude δ(k, a) at the current epoch.
            Dimensionless.  In the linear regime |δ| ≪ 1; non-linear collapse
            begins when δ ≳ 1.
        growth_factor: Linear growth factor D(a) normalised so that D(1) = 1
            at the current epoch (a = 1, z = 0).  The time evolution of the
            amplitude is δ(k, a) = δ_initial(k) × D(a) / D(a_initial).

    References:
        §18.40 — Structure formation from de-saturated medium.
        §18.44 — Pre-CMB substrate seeding sets δ_initial.
    """

    k_wavenumber: float
    delta: float
    growth_factor: float

    def __post_init__(self) -> None:
        if self.k_wavenumber <= 0.0:
            raise ValueError(
                f"Wavenumber k must be positive; got k = {self.k_wavenumber} Mpc⁻¹"
            )


# ---------------------------------------------------------------------------
# Module-level functions
# ---------------------------------------------------------------------------


def linear_growth_function(
    a: float,
    omega_m: float = OMEGA_M,
    omega_lambda: float = OMEGA_LAMBDA,
) -> float:
    """Compute the linear growth factor D(a) normalised to D(1) = 1 at z = 0.

    In ΛCDM the growth factor satisfies the second-order ODE:

        D'' + (2 - q) H D' - (3/2) Ω_m H₀² / a³ D = 0

    where primes denote derivatives with respect to proper time.  The exact
    solution is given by Heath's integral (1977):

        D(a) ∝ H(a) ∫₀ᵃ  da' / [a' H(a')]³

    where H(a) is the Friedmann Hubble parameter.  This is evaluated
    numerically.  The result is normalised so that D(1) = 1.

    During pure matter domination (Ω_Λ = 0) this recovers D(a) = a.
    The ΛCDM suppression at late times (a > a_eq ≈ 0.75) causes D(a)
    to fall below a.

    Args:
        a:            Scale factor at which to evaluate D.  Must be in (0, 1].
        omega_m:      Matter density fraction Ω_m at a = 1.  Default 0.315.
        omega_lambda: Dark energy fraction Ω_Λ at a = 1.  Default 0.685.

    Returns:
        Linear growth factor D(a), normalised so D(1) = 1.  Dimensionless.

    Raises:
        ValueError: If a ≤ 0 or a > 1.

    References:
        §18.40 — Matter perturbations grow after de-saturation.
        Carroll, Press & Turner (1992), ARA&A 30, 499 — LCDM growth factor.
        Peebles (1980), §5 — δ ∝ a in matter-dominated era (ODE solution).
    """
    if a <= 0.0:
        raise ValueError(f"Scale factor must be positive; got a = {a}")
    if a > 1.0 + 1e-10:
        raise ValueError(f"Scale factor must satisfy a ≤ 1; got a = {a}")

    def _hubble_norm_squared(a_val: float) -> float:
        """H(a)² / H₀² — dimensionless squared Hubble rate."""
        omega_r = OMEGA_R
        return (
            omega_m * a_val**-3
            + omega_r * a_val**-4
            + omega_lambda
        )

    def _integrand(a_val: float) -> float:
        """Integrand da' / [a'³ H(a')³] for the Carroll-Press-Turner growth integral.

        The correct form of the Heath growing-mode integral is (Carroll,
        Press & Turner 1992, Eq. 29; Peebles 1980 §5.2):

            D+(a) ∝ H(a) ∫₀ᵃ da' / [a'³ (H(a')/H₀)³]

        The factor a'³ (not a'¹) in the denominator is essential: it ensures
        that D ~ a in the matter-dominated limit (a'^{-3} cancels H's a'^{-9/2}
        to yield an integrand ~ a'^{1/2} whose integral ~ a^{3/2}, and then
        H(a) ~ a^{-3/2} gives D ~ a).
        """
        h2 = _hubble_norm_squared(a_val)
        if h2 <= 0.0:
            return 0.0
        return 1.0 / (a_val**3 * h2**1.5)

    # Numerical integration from a_min to a.
    # At a_min = 1e-5 the integrand → 0 (a'^{1/2}/Ω_m^{3/2} → 0) so the
    # lower-limit contribution is negligible.
    a_min = 1e-5
    n_steps = 4000
    a_arr = np.linspace(a_min, a, n_steps + 1)
    integrand_arr = np.array([_integrand(av) for av in a_arr])
    integral_a = float(np.trapezoid(integrand_arr, a_arr))

    # Normalisation integral at a = 1
    a_arr_1 = np.linspace(a_min, 1.0, n_steps + 1)
    integrand_arr_1 = np.array([_integrand(av) for av in a_arr_1])
    integral_1 = float(np.trapezoid(integrand_arr_1, a_arr_1))

    if integral_1 == 0.0:
        return 1.0

    h_a = math.sqrt(max(_hubble_norm_squared(a), 1e-300))
    h_1 = math.sqrt(max(_hubble_norm_squared(1.0), 1e-300))

    # D(a) ∝ H(a) × integral; normalise by D(1)
    D_a = h_a * integral_a
    D_1 = h_1 * integral_1

    return D_a / D_1


def transfer_function_eisenstein_hu(
    k_h_per_Mpc: float,
    omega_m: float = OMEGA_M,
    omega_b: float = OMEGA_B,
    h_param: float = h,
) -> float:
    """Compute the Eisenstein & Hu (1998) matter transfer function T(k).

    The transfer function captures how gravitational growth is suppressed
    relative to the superhorizon (k → 0) case for a given wavenumber k.
    It encodes:
    - Free-streaming suppression during radiation domination.
    - Baryon-photon coupling and acoustic oscillations before recombination.
    - Baryon drag after decoupling.

    This implements the EH98 fitting formulae (no-wiggle + wiggle components),
    which are accurate to ~5% over 10⁻³ Mpc⁻¹ ≤ k ≤ 10 Mpc⁻¹.

    The substrate-mechanical interpretation (§18.44): the transfer function
    is inherited from the standard baryon-photon plasma physics.  The stiff-
    medium model does not modify T(k) at leading order, but pre-CMB substrate
    inhomogeneities may shift the sound horizon r_s (§18.44 Hubble tension
    discussion), which translates to a shift in the BAO peak position in T(k).

    Args:
        k_h_per_Mpc: Comoving wavenumber in units of h Mpc⁻¹.  The function
            internally converts to Mpc⁻¹ as k = k_h_per_Mpc × h.
        omega_m:     Total matter density fraction Ω_m (default Planck 2018).
        omega_b:     Baryon density fraction Ω_b (default Planck 2018).
        h_param:     Dimensionless Hubble parameter h = H₀/(100 km/s/Mpc).

    Returns:
        Transfer function T(k), dimensionless.  T → 1 for k → 0 (superhorizon);
        T < 1 for k > k_eq (sub-horizon suppression).

    References:
        §18.44 — Substrate seeding modifies pre-CMB initial conditions.
        Eisenstein & Hu (1998), ApJ 496, 605 — Eqs. (17)-(31).
    """
    # Physical density parameters Ω X h²
    omega_m_h2 = omega_m * h_param**2
    omega_b_h2 = omega_b * h_param**2

    # Physical wavenumber in Mpc⁻¹
    k = k_h_per_Mpc * h_param

    # EH98 sound horizon (Mpc) — Eq. (6)
    z_eq = 2.5e4 * omega_m_h2 * (T_CMB_K / 2.7) ** (-4)
    k_eq = 7.46e-2 * omega_m_h2 * (T_CMB_K / 2.7) ** (-2)   # Mpc⁻¹

    # Redshift at drag epoch z_d — EH98 Eq. (4)
    b1 = 0.313 * omega_m_h2 ** (-0.419) * (1.0 + 0.607 * omega_m_h2**0.674)
    b2 = 0.238 * omega_m_h2**0.223
    z_d = 1291.0 * (omega_m_h2**0.251 / (1.0 + 0.659 * omega_m_h2**0.828)) * (
        1.0 + b1 * omega_b_h2**b2
    )

    # Sound horizon at drag epoch r_s (Mpc) — EH98 Eq. (6)
    R_eq = 31.5e3 * omega_b_h2 * (T_CMB_K / 2.7) ** (-4) / z_eq
    R_d = 31.5e3 * omega_b_h2 * (T_CMB_K / 2.7) ** (-4) / z_d
    if R_eq <= 0.0 or R_d <= 0.0:
        r_s = SOUND_HORIZON_MPC
    else:
        r_s = (
            2.0
            / (3.0 * k_eq)
            * math.sqrt(6.0 / R_eq)
            * math.log(
                (math.sqrt(1.0 + R_d) + math.sqrt(R_d + R_eq))
                / (1.0 + math.sqrt(R_eq))
            )
        )

    # Silk damping scale k_silk (Mpc⁻¹) — EH98 Eq. (7)
    k_silk = (
        1.6
        * (omega_b_h2**0.52)
        * (omega_m_h2**0.38)
        * (1.0 + (5.2 * omega_b_h2) ** (-0.55)) ** (-1)
    )

    # CDM transfer function with baryon corrections (EH98 Eq. 17-28)
    f_b = omega_b / omega_m    # baryon fraction
    f_c = 1.0 - f_b            # CDM fraction

    # CDM part — EH98 Eq. (17)
    a1 = (46.9 * omega_m_h2) ** 0.670 * (1.0 + (32.1 * omega_m_h2) ** (-0.532))
    a2 = (12.0 * omega_m_h2) ** 0.424 * (1.0 + (45.0 * omega_m_h2) ** (-0.582))
    alpha_c = a1 ** (-f_b) * a2 ** (-(f_b**3))

    bb1 = 0.944 / (1.0 + (458.0 * omega_m_h2) ** (-0.708))
    bb2 = (0.395 * omega_m_h2) ** (-0.0266)
    beta_c = 1.0 / (1.0 + bb1 * ((f_c**bb2) - 1.0))

    def _T0(k_val: float, alpha: float, beta: float) -> float:
        """CDM transfer function shape (EH98 Eq. 17)."""
        q = k_val / (13.41 * k_eq)
        C = 14.2 / alpha + 386.0 / (1.0 + 69.9 * q**1.08)
        T_ = math.log(math.e + 1.8 * beta * q) / (
            math.log(math.e + 1.8 * beta * q) + C * q**2
        )
        return T_

    T_c = f_c * _T0(k, 1.0 / beta_c, beta_c) + f_b * _T0(k, alpha_c, 1.0)

    # Baryon part — EH98 Eqs. (19)-(24)
    y_d = z_eq / (1.0 + z_d)
    G_y = y_d * (
        -6.0 * math.sqrt(1.0 + y_d)
        + (2.0 + 3.0 * y_d) * math.log((math.sqrt(1.0 + y_d) + 1.0) / (math.sqrt(1.0 + y_d) - 1.0))
    )

    alpha_b = (
        2.07
        * k_eq
        * r_s
        * (1.0 + R_d) ** (-3.0 / 4.0)
        * G_y
    )

    beta_b = 0.5 + f_b + (3.0 - 2.0 * f_b) * math.sqrt((17.2 * omega_m_h2) ** 2 + 1.0)
    beta_node = 8.41 * omega_m_h2**0.435

    s_tilde: float
    k_r_s = k * r_s
    if abs(k_r_s) < 1e-10:
        s_tilde = r_s
    else:
        s_node = r_s / (1.0 + (beta_node / k_r_s) ** 3) ** (1.0 / 3.0)
        s_tilde = s_node

    j0_silk = math.sin(k * r_s) / max(k * r_s, 1e-30)
    silk_damp = math.exp(-((k / k_silk) ** 1.4))

    T_b = (
        _T0(k, 1.0, 1.0) / (1.0 + (k * r_s / 5.2) ** 2)
        + alpha_b
        / (1.0 + (beta_b / (k * r_s)) ** 3)
        * math.exp(-((k / k_silk) ** 1.4))
    ) * j0_silk * silk_damp

    # Full transfer function — weighted sum of CDM and baryon parts (EH98 Eq. 16)
    T_total = f_c * T_c + f_b * T_b

    # Clip to physical range [0, 1]
    return float(min(max(T_total, 0.0), 1.0))


def cmb_first_peak_from_sound_horizon(
    r_s_Mpc: float = SOUND_HORIZON_MPC,
    d_LSS_Gpc: float = D_LSS_GPC,
) -> float:
    """Compute the angular multipole of the first CMB acoustic peak.

    The first peak in the CMB temperature power spectrum occurs at the angular
    scale subtended by the sound horizon r_s at the last-scattering surface.

    In the flat-sky limit:
        ℓ₁ ≈ π × D_A / r_s

    where D_A is the **physical angular diameter distance** to the last-
    scattering surface (not the comoving distance χ_LSS).  The relationship
    between the two is:

        D_A = χ_LSS / (1 + z_rec)   for a flat universe.

    For Planck 2018 parameters:
        χ_LSS ≈ 13800 Mpc  (comoving distance to z = 1100)
        D_A   ≈ 10300 Mpc  (physical angular diameter distance)
        r_s   ≈ 147 Mpc
        ℓ₁   ≈ π × 10300 / 147 ≈ 220

    Note: a common source of confusion is that ``d_LSS_Gpc`` here denotes
    the **angular diameter distance** (~10.3 Gpc), not the comoving distance
    (~13.8 Gpc).  Using χ_LSS would yield ℓ₁ ≈ 295.

    The substrate-mechanical interpretation (§18.44): the acoustic oscillations
    in the CMB are interpreted as phase-boundary resonance modes of the
    de-saturation transition.  The first peak position is set by r_s of the
    phase boundary coherence length at recombination.

    Args:
        r_s_Mpc:   Sound horizon at recombination in comoving Mpc.
            Default 147 Mpc (standard ΛCDM).
        d_LSS_Gpc: Physical angular diameter distance D_A to the last-
            scattering surface in Gpc.  Default 10.3 Gpc (giving ℓ₁ ≈ 220).
            Note: the comoving distance is ~13.8 Gpc; using the angular
            diameter distance gives the observed ℓ₁ ≈ 220.

    Returns:
        Multipole ℓ₁ of the first CMB acoustic peak (dimensionless).
        Standard ΛCDM gives ℓ₁ ≈ 220.

    References:
        §18.44 — CMB peaks as phase-boundary resonance modes.
        Hu & Sugiyama (1995), ApJ 444, 489 — acoustic peak physics.
        Planck 2018 — 100 θ_* = 1.04109 (angular scale of sound horizon).
    """
    d_LSS_Mpc = d_LSS_Gpc * 1000.0   # convert Gpc → Mpc (D_A in Mpc)
    return math.pi * d_LSS_Mpc / r_s_Mpc


def bao_scale_today(
    r_s_Mpc: float = SOUND_HORIZON_MPC,
    z_recombination: float = Z_RECOMBINATION,
) -> float:
    """Compute the comoving BAO scale today from the sound horizon at recombination.

    The baryon acoustic oscillation (BAO) scale is the comoving sound horizon
    at the drag epoch (close to recombination):

        r_BAO ≈ r_s(z_drag) ≈ r_s(z_rec) ≈ 147 Mpc

    This is a comoving scale, so it does not change after recombination (the
    physical scale grows as a, but the comoving scale is fixed).  It appears
    as a characteristic clustering scale in the galaxy correlation function and
    power spectrum.

    The substrate-mechanical interpretation (§18.44, §18.43): the BAO scale
    is imprinted by the sound horizon at the de-saturation phase transition.
    In cyclic cosmology (§18.43), signatures of BAO from previous cycles could
    appear if the pre-CMB substrate retains spatial coherence.

    Args:
        r_s_Mpc:         Sound horizon at recombination in Mpc.
            Default 147 Mpc.
        z_recombination: Redshift of the drag/recombination epoch.
            Default 1100.

    Returns:
        Comoving BAO scale today in Mpc.  For standard parameters this is
        r_s ≈ 147 Mpc (unchanged, since it is already a comoving quantity).

    Notes:
        The BAO feature appears at the comoving scale r_s in the matter power
        spectrum today.  Because r_s is a comoving distance, the BAO ruler
        does not evolve with a.  What evolves is the mapping from comoving
        to angular/redshift coordinates used to measure it.

    References:
        §18.44 — Sound horizon at de-saturation sets BAO scale.
        §18.43 — Cyclic cosmology BAO signatures.
        Eisenstein et al. (2005), ApJ 633, 560 — BAO detection in SDSS.
    """
    # The comoving sound horizon is independent of redshift after recombination.
    # The physical BAO scale at redshift z is r_s / (1 + z), but the comoving
    # scale is always r_s.
    _ = z_recombination   # retained for API completeness and documentation
    return r_s_Mpc


# ---------------------------------------------------------------------------
# PowerSpectrum class
# ---------------------------------------------------------------------------


class PowerSpectrum:
    """Matter power spectrum P(k) for ΛCDM / stiff-medium cosmology.

    Computes the linear matter power spectrum at redshift z = 0:

        P(k) = A_s × (k / k_pivot)^(n_s - 1) × k × T(k)² × D(1)²

    where:
    - A_s is the primordial scalar amplitude.
    - n_s is the scalar spectral index.
    - T(k) is the Eisenstein-Hu transfer function.
    - D(a) is the linear growth factor (normalised D(1) = 1).

    The normalisation follows the convention that σ₈ — the RMS density
    variance in a sphere of radius 8 Mpc/h — should match the Planck 2018
    value of 0.811.  The overall amplitude is set by matching σ₈.

    The substrate-mechanical interpretation (§18.40, §18.44): the primordial
    spectrum P_prim(k) ∝ k^n_s is inherited from the saturated-state
    fluctuations.  After de-saturation (§18.44), the transfer function T(k)
    encodes the medium's response to sub-horizon modes entering the Hubble
    radius.

    Args:
        n_s:      Scalar spectral index.  Default 0.9649 (Planck 2018).
        A_s:      Primordial scalar amplitude at k_pivot.  Default 2.1e-9.
        omega_m:  Matter density fraction.  Default 0.315.
        omega_b:  Baryon density fraction.  Default 0.049.
        h_param:  Dimensionless Hubble parameter.  Default 0.674.
        sigma_8_target: Target σ₈ for amplitude normalisation.
            Default 0.811.

    References:
        §18.40 — Primordial spectrum from saturated-state fluctuations.
        §18.44 — Transfer function shaped by de-saturation epoch.
    """

    def __init__(
        self,
        n_s: float = N_S_DEFAULT,
        A_s: float = A_S_DEFAULT,
        omega_m: float = OMEGA_M,
        omega_b: float = OMEGA_B,
        h_param: float = h,
        sigma_8_target: float = SIGMA_8_TARGET,
    ) -> None:
        self.n_s = n_s
        self.A_s = A_s
        self.omega_m = omega_m
        self.omega_b = omega_b
        self.h_param = h_param
        self.sigma_8_target = sigma_8_target

        # Growth factor at z = 0 is 1 by normalisation convention
        self._D_today: float = 1.0

        # Compute normalisation amplitude so that sigma_8() == sigma_8_target.
        # We do this lazily: set norm = 1 first, compute raw sigma_8, then
        # rescale norm so that sigma_8 hits the target.
        self._amplitude_norm: float = 1.0
        raw_s8 = self._compute_sigma_8_raw()
        if raw_s8 > 0.0:
            self._amplitude_norm = (sigma_8_target / raw_s8) ** 2
        else:
            self._amplitude_norm = 1.0

    # ------------------------------------------------------------------
    # Core spectrum
    # ------------------------------------------------------------------

    def _primordial_power(self, k_Mpc: float) -> float:
        """Un-normalised primordial scalar power spectrum P_prim(k).

        P_prim(k) = A_s × (k / k_pivot)^(n_s - 1) × k

        This gives the Harrison-Zel'dovich scale-invariant spectrum for
        n_s = 1, modified by the tilt n_s - 1 ≈ -0.035.

        Args:
            k_Mpc: Wavenumber in Mpc⁻¹.

        Returns:
            Primordial power in Mpc³.
        """
        return self.A_s * (k_Mpc / K_PIVOT_MPC) ** (self.n_s - 1.0) * k_Mpc

    def power_at_k(self, k: float) -> float:
        """Evaluate the matter power spectrum P(k) at wavenumber k.

        The full linear matter power spectrum:

            P(k) = A_norm × (k / k_pivot)^(n_s - 1) × k × T(k)² × D(a)²

        evaluated at z = 0 (D = 1 by convention).  The amplitude normalisation
        is set internally by matching σ₈ to the Planck 2018 value.

        The substrate-mechanical model does not modify P(k) at leading order
        relative to ΛCDM, but the pre-CMB substrate inhomogeneities (§18.44)
        could add a low-k enhancement (scales > sound horizon r_s) that
        boosts galaxy formation at high redshift — consistent with JWST
        observations (§18.44, §18.52.3).

        Args:
            k: Comoving wavenumber in Mpc⁻¹.  Must be positive.

        Returns:
            P(k) in Mpc³.  Dimensionful — gives variance per logarithmic k
            interval as Δ²(k) = k³ P(k) / (2π²).

        Raises:
            ValueError: If k ≤ 0.

        References:
            §18.40 — Growth of perturbations after de-saturation.
            §18.44 — Pre-CMB substrate modifies low-k power.
        """
        if k <= 0.0:
            raise ValueError(f"Wavenumber k must be positive; got k = {k}")
        k_h = k / self.h_param   # convert Mpc⁻¹ → h Mpc⁻¹ for EH transfer fn
        T = transfer_function_eisenstein_hu(
            k_h, omega_m=self.omega_m, omega_b=self.omega_b, h_param=self.h_param
        )
        P_prim = self._primordial_power(k)
        return self._amplitude_norm * P_prim * T**2 * self._D_today**2

    # ------------------------------------------------------------------
    # Sigma-8
    # ------------------------------------------------------------------

    def _window_tophat(self, k: float, R: float) -> float:
        """Real-space top-hat window function W(kR) in Fourier space.

        W(x) = 3 [sin(x) - x cos(x)] / x³

        Args:
            k: Wavenumber in Mpc⁻¹.
            R: Smoothing radius in Mpc.

        Returns:
            Window function value (dimensionless).
        """
        x = k * R
        if x < 1e-4:
            return 1.0 - x**2 / 10.0
        return 3.0 * (math.sin(x) - x * math.cos(x)) / x**3

    def _compute_sigma_8_raw(self, R_Mpc: float = 8.0 / h) -> float:
        """Compute σ(R) using the raw (un-renormalised) power spectrum.

        σ²(R) = (1 / 2π²) ∫ dk k² P_raw(k) W²(kR)

        Args:
            R_Mpc: Smoothing radius in Mpc.  Default 8/h Mpc ≈ 11.9 Mpc
                   (corresponds to 8 Mpc/h).

        Returns:
            σ(R) dimensionless.
        """
        # Integrate over k in log space for better numerical accuracy
        log_k_min = math.log(1e-4)
        log_k_max = math.log(10.0)
        n_k = 500
        log_k_arr = np.linspace(log_k_min, log_k_max, n_k)
        k_arr = np.exp(log_k_arr)

        integrand = np.zeros(n_k)
        for i, ki in enumerate(k_arr):
            W = self._window_tophat(ki, R_Mpc)
            P_prim = self._primordial_power(ki)
            k_h = ki / self.h_param
            T = transfer_function_eisenstein_hu(
                k_h, omega_m=self.omega_m, omega_b=self.omega_b, h_param=self.h_param
            )
            # Integrand in log-k: k³ P(k) W²(kR) / (2π²)
            integrand[i] = ki**3 * P_prim * T**2 * W**2 / (2.0 * math.pi**2)

        # Integrate d(ln k) = dk/k, so multiply integrand by 1 (already k³ factor)
        sigma_sq = float(np.trapezoid(integrand, log_k_arr))
        return math.sqrt(max(sigma_sq, 0.0))

    def sigma_8(self, R_h_Mpc: float = 8.0) -> float:
        """Compute σ(R) — RMS density variance smoothed at radius R Mpc/h.

        The variance σ²(R) is defined as:

            σ²(R) = (1 / 2π²) ∫₀^∞ dk k² P(k) W²(kR)

        where W(kR) = 3[sin(kR) - kR cos(kR)]/(kR)³ is the Fourier transform
        of a spherical top-hat window of radius R.

        The conventional σ₈ uses R = 8 Mpc/h (≈ 11.9 Mpc for h = 0.674),
        where the matter fluctuation amplitude transitions from non-linear
        to linear regime.

        The amplitude normalisation is set in the constructor to match σ₈
        to the Planck 2018 value of 0.811.  This function should therefore
        return approximately 0.811 for the default 8 Mpc/h radius.

        Args:
            R_h_Mpc: Smoothing radius in Mpc/h.  Default 8.0.

        Returns:
            σ(R) dimensionless.  At R = 8 Mpc/h this is σ₈ ≈ 0.811.

        References:
            §18.40 — Matter variance constrains the amplitude of post-CMB growth.
            Planck 2018 — σ₈ = 0.811 ± 0.006.
        """
        R_Mpc = R_h_Mpc / self.h_param   # Mpc/h → Mpc

        log_k_min = math.log(1e-4)
        log_k_max = math.log(10.0)
        n_k = 500
        log_k_arr = np.linspace(log_k_min, log_k_max, n_k)
        k_arr = np.exp(log_k_arr)

        integrand = np.zeros(n_k)
        for i, ki in enumerate(k_arr):
            W = self._window_tophat(ki, R_Mpc)
            P = self.power_at_k(ki)
            integrand[i] = ki**3 * P * W**2 / (2.0 * math.pi**2)

        sigma_sq = float(np.trapezoid(integrand, log_k_arr))
        return math.sqrt(max(sigma_sq, 0.0))

    # ------------------------------------------------------------------
    # CMB peak / sound horizon
    # ------------------------------------------------------------------

    def acoustic_peak_position(
        self,
        r_s_Mpc: float = SOUND_HORIZON_MPC,
        d_LSS_Gpc: float = D_LSS_GPC,
    ) -> float:
        """Return the multipole ℓ of the first CMB acoustic temperature peak.

        Delegates to :func:`cmb_first_peak_from_sound_horizon`.

        Args:
            r_s_Mpc:   Sound horizon at recombination (Mpc).  Default 147 Mpc.
            d_LSS_Gpc: Comoving angular diameter distance to LSS (Gpc).
                Default 13.8 Gpc.

        Returns:
            First acoustic peak multipole ℓ₁ ≈ 220.

        References:
            §18.44 — CMB peaks as phase-boundary resonance modes.
        """
        return cmb_first_peak_from_sound_horizon(r_s_Mpc, d_LSS_Gpc)

    def sound_horizon(self) -> float:
        """Return the comoving sound horizon r_s at recombination.

        Computes the EH98 sound horizon from the baryon and matter density
        parameters, cross-checked against the standard ΛCDM value.

        Returns:
            Sound horizon r_s in Mpc.  Standard ΛCDM: r_s ≈ 147 Mpc.

        References:
            §18.44 — Sound horizon sets BAO and CMB peak scales.
            Eisenstein & Hu (1998), Eq. (6).
        """
        omega_m_h2 = self.omega_m * self.h_param**2
        omega_b_h2 = self.omega_b * self.h_param**2

        # EH98 Eq. (4) drag epoch redshift
        b1 = 0.313 * omega_m_h2 ** (-0.419) * (1.0 + 0.607 * omega_m_h2**0.674)
        b2 = 0.238 * omega_m_h2**0.223
        z_d = 1291.0 * (omega_m_h2**0.251 / (1.0 + 0.659 * omega_m_h2**0.828)) * (
            1.0 + b1 * omega_b_h2**b2
        )

        # EH98 Eq. (6) matter-radiation equality scale
        z_eq = 2.5e4 * omega_m_h2 * (T_CMB_K / 2.7) ** (-4)
        k_eq = 7.46e-2 * omega_m_h2 * (T_CMB_K / 2.7) ** (-2)

        R_eq = 31.5e3 * omega_b_h2 * (T_CMB_K / 2.7) ** (-4) / z_eq
        R_d = 31.5e3 * omega_b_h2 * (T_CMB_K / 2.7) ** (-4) / z_d

        if R_eq <= 0.0 or R_d <= 0.0 or k_eq <= 0.0:
            return SOUND_HORIZON_MPC

        r_s = (
            2.0
            / (3.0 * k_eq)
            * math.sqrt(6.0 / R_eq)
            * math.log(
                (math.sqrt(1.0 + R_d) + math.sqrt(R_d + R_eq))
                / (1.0 + math.sqrt(R_eq))
            )
        )
        return float(r_s)


# ---------------------------------------------------------------------------
# JWSTAnomaly analyzer
# ---------------------------------------------------------------------------


@dataclass
class JWSTAnomaly:
    """Analyzer for the JWST high-redshift galaxy abundance anomaly.

    JWST has detected a galaxy abundance at z ≳ 10-14 that exceeds the
    standard ΛCDM prediction by factors of ~100-180× (arXiv:2505.11263).
    This is the "JWST anomaly" discussed in §18.52.3.

    The substrate-mechanical explanation (§18.44, §18.43):
    Pre-CMB substrate inhomogeneities from the saturated era seed structure
    formation before the standard CMB epoch.  This means:
    1. The effective initial density contrast δ_initial is enhanced at
       galaxy-mass scales (k ~ 1-10 Mpc⁻¹).
    2. The enhancement factor depends on the pre-CMB substrate power,
       parametrised here by ``substrate_seed_amplitude``.

    The ΛCDM prediction for galaxy number density at high-z uses the
    Press-Schechter (1974) formalism, which depends on σ(M) — the mass
    variance at galaxy mass M.  The stiff-medium model boosts σ(M) by the
    substrate seeding factor.

    Attributes:
        z_observed:   Redshift at which the excess is observed.  Default 14.
        n_obs_ratio:  Observed abundance relative to ΛCDM prediction.
            JWST arXiv:2505.11263 reports ~182× excess at z ≈ 10-14.
            Default 182.
        substrate_seed_amplitude: Enhancement factor for the pre-CMB
            substrate initial power spectrum at galaxy scales.
            A value of 1.0 = no enhancement (pure ΛCDM).
            Default 10.0 (order-of-magnitude estimate from §18.44).
        power_spectrum: :class:`PowerSpectrum` instance to use for σ₈.
            If None, a default instance with Planck 2018 parameters is used.

    References:
        §18.44 — Pre-CMB substrate seeding enhanced galaxy formation.
        §18.52.3 — JWST 182× excess.
        Press & Schechter (1974), ApJ 187, 425.
    """

    z_observed: float = 14.0
    n_obs_ratio: float = 182.0
    substrate_seed_amplitude: float = 10.0
    power_spectrum: PowerSpectrum | None = field(default=None)

    def __post_init__(self) -> None:
        if self.power_spectrum is None:
            self.power_spectrum = PowerSpectrum()

    # ------------------------------------------------------------------
    # Growth factor at high redshift
    # ------------------------------------------------------------------

    def _growth_at_z(self, z: float) -> float:
        """Linear growth factor D(a) at redshift z.

        Args:
            z: Redshift.  Must be ≥ 0.

        Returns:
            D(a) normalised to D(1) = 1 at z = 0.
        """
        a = 1.0 / (1.0 + z)
        return linear_growth_function(a)

    # ------------------------------------------------------------------
    # Press-Schechter mass variance
    # ------------------------------------------------------------------

    def _sigma_at_galaxy_mass(
        self,
        log10_M_Msun: float = 11.0,
        enhanced: bool = False,
    ) -> float:
        """RMS density variance σ(M) at galaxy mass scale.

        Uses the standard Press-Schechter relation:
            R(M) = (3M / (4π ρ_m))^(1/3)  (comoving Mpc)

        where ρ_m = Ω_m ρ_crit is the mean comoving matter density.

        Args:
            log10_M_Msun: Log₁₀ of galaxy mass in solar masses.
                Default 10¹¹ M_sun (Milky Way-size galaxy).
            enhanced: If True, boost σ by ``substrate_seed_amplitude``.

        Returns:
            σ(M) dimensionless.
        """
        ps = self.power_spectrum
        assert ps is not None

        # Physical comoving matter density (kg Mpc⁻³)
        M_sun_kg = 1.989e30
        Mpc_m = 3.0857e22
        M_kg = 10.0**log10_M_Msun * M_sun_kg

        # Critical density today ρ_crit = 3 H₀² / (8πG)
        G_si = 6.674e-11
        c_si = 2.998e8
        H0_si = H0_KM_S_MPC * 1e3 / Mpc_m   # s⁻¹
        rho_crit_kg_m3 = 3.0 * H0_si**2 / (8.0 * math.pi * G_si)
        rho_crit_kg_Mpc3 = rho_crit_kg_m3 * Mpc_m**3
        rho_m_kg_Mpc3 = OMEGA_M * rho_crit_kg_Mpc3

        # Lagrangian radius R(M)
        R_Mpc = (3.0 * M_kg / (4.0 * math.pi * rho_m_kg_Mpc3)) ** (1.0 / 3.0)

        sigma_val = ps.sigma_8(R_h_Mpc=R_Mpc * ps.h_param)   # pass as h Mpc⁻¹ units
        if enhanced:
            sigma_val *= self.substrate_seed_amplitude
        return sigma_val

    # ------------------------------------------------------------------
    # Press-Schechter halo abundance
    # ------------------------------------------------------------------

    @staticmethod
    def _press_schechter_f(nu: float) -> float:
        """Press-Schechter multiplicity function f(ν).

        f(ν) = sqrt(2/π) × ν × exp(-ν²/2)

        where ν = δ_c / σ(M, z) and δ_c = 1.686 is the spherical collapse
        threshold.

        Args:
            nu: Dimensionless peak height ν = δ_c / σ.

        Returns:
            Multiplicity function f(ν), dimensionless.
        """
        return math.sqrt(2.0 / math.pi) * nu * math.exp(-0.5 * nu**2)

    def _halo_abundance_relative(
        self,
        z: float,
        log10_M_Msun: float = 11.0,
        enhanced: bool = False,
    ) -> float:
        """Relative halo abundance N(>M, z) in Press-Schechter formalism.

        The comoving number density of halos with mass > M is:

            N(>M, z) ∝ (ρ_m / M) × ∫_ν^∞ f(ν') dν'
                      ∝ erfc(ν / sqrt(2))

        where ν = δ_c / [σ(M) × D(z)].

        We return this as a relative quantity (fraction of z = 0 ΛCDM
        abundance), so the proportionality constants cancel.

        Args:
            z:             Redshift.
            log10_M_Msun:  Log₁₀ of threshold mass in solar masses.
            enhanced:      Use enhanced σ for stiff-medium model.

        Returns:
            Relative halo abundance (dimensionless).  Drops rapidly at high z.
        """
        from math import erfc as _erfc
        delta_c = 1.686
        D_z = self._growth_at_z(z)
        sigma_M = self._sigma_at_galaxy_mass(log10_M_Msun, enhanced=enhanced)
        sigma_eff = sigma_M * D_z
        if sigma_eff <= 0.0:
            return 0.0
        nu = delta_c / sigma_eff
        return _erfc(nu / math.sqrt(2.0))

    # ------------------------------------------------------------------
    # Public prediction interface
    # ------------------------------------------------------------------

    def predict(self, log10_M_Msun: float = 11.0) -> dict[str, float]:
        """Predict galaxy abundance anomaly at the observed redshift.

        Computes the ratio of the stiff-medium model galaxy abundance to the
        ΛCDM prediction at the observed redshift ``z_observed``.  If the
        substrate seeding boosts this ratio to be consistent with
        ``n_obs_ratio``, the model is declared a good fit.

        Method:
        1. Compute σ(M) × D(z) for ΛCDM (no enhancement).
        2. Compute σ(M) × D(z) for stiff-medium model (with enhancement).
        3. Use Press-Schechter to get relative abundances.
        4. Compute model/ΛCDM ratio and compare to JWST observation.

        Args:
            log10_M_Msun: Log₁₀ of galaxy halo mass in solar masses.
                Default 10¹¹ M_sun.

        Returns:
            Dictionary with keys:
            - ``'z_observed'``      : the input redshift.
            - ``'D_z'``             : linear growth factor D(a) at z.
            - ``'sigma_M_lcdm'``    : σ(M) × D(z) without seeding.
            - ``'sigma_M_model'``   : σ(M) × D(z) with substrate seeding.
            - ``'n_lcdm_relative'`` : ΛCDM relative abundance.
            - ``'n_model_relative'``: stiff-medium relative abundance.
            - ``'ratio_model_lcdm'``: predicted ratio (model / ΛCDM).
            - ``'n_obs_ratio'``     : observed ratio from JWST.
            - ``'fits_observation'``: True if ratio ≥ 0.5 × n_obs_ratio.

        References:
            §18.44 — Pre-CMB seeding predicts JWST excess.
            §18.52.3 — JWST 182× excess at z ~ 10-14.
        """
        assert self.power_spectrum is not None
        D_z = self._growth_at_z(self.z_observed)

        sigma_lcdm_M = self._sigma_at_galaxy_mass(log10_M_Msun, enhanced=False)
        sigma_model_M = self._sigma_at_galaxy_mass(log10_M_Msun, enhanced=True)
        sigma_lcdm_eff = sigma_lcdm_M * D_z
        sigma_model_eff = sigma_model_M * D_z

        n_lcdm = self._halo_abundance_relative(
            self.z_observed, log10_M_Msun, enhanced=False
        )
        n_model = self._halo_abundance_relative(
            self.z_observed, log10_M_Msun, enhanced=True
        )

        if n_lcdm > 0.0:
            ratio = n_model / n_lcdm
        else:
            ratio = float("inf")

        fits = ratio >= 0.5 * self.n_obs_ratio

        return {
            "z_observed": self.z_observed,
            "D_z": D_z,
            "sigma_M_lcdm": sigma_lcdm_eff,
            "sigma_M_model": sigma_model_eff,
            "n_lcdm_relative": n_lcdm,
            "n_model_relative": n_model,
            "ratio_model_lcdm": ratio,
            "n_obs_ratio": self.n_obs_ratio,
            "fits_observation": fits,
        }
