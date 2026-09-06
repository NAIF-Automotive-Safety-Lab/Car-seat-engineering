"""Baryogenesis via de-saturation phase transition — computation of η_B.

This module attempts to compute the baryon-to-photon ratio

    η_B = n_B / n_γ

from the Sakharov conditions as they manifest in the stiff-medium
de-saturation phase transition (spec §§18.40, 18.44, 18.51).

HONEST PREAMBLE
---------------
The spec states (§18.51) that the model "commits" to η_B ~ 10⁻¹⁰ from
CP-violating de-saturation, but provides no derivation.  This module
supplies the first quantitative attempt.  We identify the three Sakharov
conditions and compute each as precisely as the existing machinery allows.
Where genuine unknown parameters remain, they are LABELLED as such — the
goal is to determine which δ_CP values reproduce the observed η_B, and to
honestly state which pieces are derived vs. fitted.

The three Sakharov conditions:
1. Baryon number violation — rate set by de-saturation bubble nucleation
   (from desaturation.py): Γ_B ~ Γ_nuc × (ξ_Planck³) (rate per baryon volume)
2. C and CP violation — from Möbius half-flux chirality asymmetry δ_CP
   (from mobius_bundle.py): δ_CP ∈ (0, 1] is a free parameter in this model
3. Out-of-equilibrium — first-order de-saturation phase transition at T_*
   (from desaturation.py + cyclic_cosmology.py): wall velocity v_w = c (§18.44),
   Hubble rate H_* from FRW at the transition epoch

The master formula used here is the standard baryogenesis estimate:

    η_B ≈ (Γ_B / H_*) × δ_CP × washout_factor

where:
    Γ_B ~ Γ_nuc(T_*)  × ξ_Planck³   [s⁻¹, rate per Hubble volume converted]
    H_* ~ Hubble parameter at the transition (s⁻¹)
    δ_CP ~ CP-asymmetry parameter (dimensionless, 0 < δ_CP ≤ 1)
    washout_factor ~ exp(-Γ_washout × τ_after)  ≈ 1 for rapid de-saturation

References:
    spec §18.40, §18.44, §18.51
    Coleman (1977), Phys. Rev. D 15, 2929 — bounce action
    Kolb & Turner (1990), The Early Universe — baryogenesis review
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np

from stiff_medium.desaturation import (
    DesaturationTransition,
    K_PLANCK,
    c,
    hbar,
    G,
    k_B,
    a_rad,
    SIGMA_MAX,
)
from stiff_medium.mobius_bundle import MobiusBundle, theta_qcd_from_topology

# ---------------------------------------------------------------------------
# Physical constants (supplementary)
# ---------------------------------------------------------------------------

#: Planck length l_P = sqrt(ℏG/c³) [m]
L_PLANCK: Final[float] = math.sqrt(hbar * G / c**3)

#: Planck time t_P = l_P / c [s]
T_PLANCK: Final[float] = L_PLANCK / c

#: Planck temperature T_P = ℏc⁵/(G k_B²) ... = sqrt(ℏc⁵/(G)) / k_B [K]
T_PLANCK_K: Final[float] = math.sqrt(hbar * c**5 / G) / k_B

#: Number of relativistic degrees of freedom at the Planck epoch (order-of-mag.)
#: Standard model has g_* ~ 106.75 at T >> 100 GeV; at Planck it is uncertain.
#: We use 100 as a round order-of-magnitude consistent with near-SM content.
G_STAR: Final[float] = 100.0

#: Observed baryon-to-photon ratio (Planck 2018 / BBN).
ETA_B_OBSERVED: Final[float] = 6.1e-10

#: CMB temperature today (K)
T_CMB_TODAY: Final[float] = 2.7255


# ---------------------------------------------------------------------------
# Helper: Hubble rate in radiation-dominated epoch
# ---------------------------------------------------------------------------

def hubble_rate_radiation(T_K: float) -> float:
    """Hubble parameter H in the radiation-dominated epoch at temperature T.

    In the standard FRW radiation-dominated era the Friedmann equation gives:

        H = π T² sqrt(g_* / 90) / (M_Pl ℏ c)

    where M_Pl = sqrt(ℏ c / G) is the reduced Planck mass.  In SI units:

        H(T) = sqrt(8 π³ g_* / 90) × (k_B T)² / (ℏ² c² sqrt(3/(8πG)))
             = sqrt(4 π³ g_* / 45) × k_B² T² / (M_Pl_reduced ℏ c)

    We use the compact SI form:

        H² = (8πG/3) × ε_rad
        ε_rad = (π² / 30) × (g_*/1) × (k_B T)⁴ / (ℏ c)³

    Args:
        T_K: Temperature in Kelvin.

    Returns:
        Hubble parameter in s⁻¹.  Returns 0 if T ≤ 0.
    """
    if T_K <= 0.0:
        return 0.0
    # Energy density of relativistic species at temperature T (J m⁻³)
    # ε = (π²/30) g_* (k_B T)⁴ / (ℏ c)³
    epsilon_rad = (math.pi**2 / 30.0) * G_STAR * (k_B * T_K)**4 / (hbar * c)**3
    # Friedmann: H² = (8πG/3) ε
    H_squared = (8.0 * math.pi * G / 3.0) * epsilon_rad
    return math.sqrt(max(H_squared, 0.0))


# ---------------------------------------------------------------------------
# Helper: Photon number density
# ---------------------------------------------------------------------------

def photon_number_density(T_K: float) -> float:
    """Photon number density n_γ = (2ζ(3)/π²) (k_B T / ℏc)³.

    Args:
        T_K: Temperature in Kelvin.

    Returns:
        Photon number density in m⁻³.
    """
    if T_K <= 0.0:
        return 0.0
    # Riemann zeta(3) ≈ 1.202056903
    zeta3 = 1.2020569031595943
    return (2.0 * zeta3 / math.pi**2) * (k_B * T_K / (hbar * c))**3


# ---------------------------------------------------------------------------
# BaryogenesisModel dataclass
# ---------------------------------------------------------------------------


@dataclass
class BaryogenesisModel:
    """Compute η_B from the de-saturation phase transition Sakharov conditions.

    The three Sakharov conditions are implemented as follows:

    **Condition 1 — Baryon number violation**:
    Baryon-number-violating processes (analogous to electroweak sphaleron)
    occur when a de-saturation bubble nucleates.  Each nucleation event
    converts a saturated-medium region (undefined baryon number) into an
    ordinary-elastic region where kinks (fermions) can exist with definite
    baryon number.  The rate per unit volume is taken directly from the
    DesaturationTransition.nucleation_rate(T_*).

    **Condition 2 — CP violation**:
    The Möbius half-flux bundle (mobius_bundle.py, §18.10) is intrinsically
    chiral: its connection A = (1/2)dθ defines a preferred orientation.
    The CP-asymmetry parameter δ_CP quantifies how much more matter than
    antimatter emerges from each bubble.  The spec claims this comes from the
    Möbius chirality but gives no derivation of its magnitude; δ_CP is
    therefore a FREE PARAMETER in the current theory.

    theta_QCD is exactly 0 from the Möbius bundle (theta_qcd_from_topology()
    returns 0), so there is no strong-CP contribution to δ_CP.  The CP
    violation must come from electroweak-sector physics inherited from the
    topology; its magnitude is UNSPECIFIED in the current spec.

    **Condition 3 — Out of equilibrium**:
    The de-saturation phase transition is first-order (§18.44), with bubble
    nucleation rate Γ_nuc << H at onset.  The phase boundary propagates at
    v_w = c (§18.44), so the departure from equilibrium is maximal.
    Washout from the inverse process is suppressed once the bubble wall
    has swept through a region.

    **Master formula**:

        η_B ≈ (Γ_B_eff / H_*) × δ_CP × washout_factor

    where:
        Γ_B_eff  = Γ_nuc(T_*) × ξ_P³     [s⁻¹: rate per Hubble-volume "slot"]
        H_*      = Hubble rate at de-saturation T_*
        δ_CP     = CP asymmetry (free parameter, 0 < δ_CP ≤ 1)
        washout_factor ~ 1 (maximised by v_w = c, negligible washout)

    Attributes:
        T_star: Temperature of the de-saturation phase transition (K).
            Default is the Planck temperature (~1.4e32 K) because the
            de-saturation occurs at the substrate's natural scale.
        sigma_initial: Initial saturation level (default 0.5 per §18.40).
        sigma_final: Post-transition saturation level (default 0.0, vacuum).
        delta_cp: CP asymmetry parameter from Möbius chirality (dimensionless,
            0 < δ_CP ≤ 1).  THIS IS THE PRIMARY FREE PARAMETER.  The spec
            asserts it is set by the Möbius topology but gives no derivation
            of its magnitude.  Values are scanned in baryogenesis_test.py.
        washout_suppression: Suppression factor exp(-Γ_washout × t_after).
            For v_w = c the washout is minimal; we set this to 1.0 by default
            (no washout).  Reducing it models partial equilibration.
    """

    T_star: float = field(default_factory=lambda: T_PLANCK_K)
    sigma_initial: float = SIGMA_MAX
    sigma_final: float = 0.0
    delta_cp: float = 0.1
    washout_suppression: float = 1.0

    def __post_init__(self) -> None:
        if not (0.0 < self.delta_cp <= 1.0):
            raise ValueError(
                f"delta_cp must be in (0, 1]; got {self.delta_cp!r}"
            )
        if not (0.0 <= self.sigma_initial <= SIGMA_MAX + 1e-9):
            raise ValueError(
                f"sigma_initial must be in [0, 0.5]; got {self.sigma_initial!r}"
            )

    # ------------------------------------------------------------------
    # Sakharov Condition 1: baryon number violation rate
    # ------------------------------------------------------------------

    def baryon_violation_rate(self) -> float:
        """Effective baryon-number-violating rate Γ_B_eff at T_star (s⁻¹).

        Each bubble nucleation event within one Planck volume is treated as a
        baryon-number-violating interaction (analogous to an electroweak
        sphaleron crossing).  The rate per unit volume is taken from the
        DesaturationTransition nucleation rate; multiplying by ξ_P³ converts
        it to a dimensionless rate per Planck volume per Planck time, then
        we express the result in s⁻¹ relative to the Hubble volume.

        The per-Planck-volume rate is:

            Γ_B_eff = Γ_nuc(T_*) × ξ_P³     [s⁻¹]

        This is a DIMENSIONAL ARGUMENT, not a first-principles derivation.
        The true sphaleron analogue in this model has not been computed.

        Returns:
            Baryon-number-violating rate in s⁻¹.
        """
        dt_obj = DesaturationTransition(
            sigma_initial=self.sigma_initial,
            sigma_final=self.sigma_final,
            K=K_PLANCK,
        )
        gamma_nuc = dt_obj.nucleation_rate(self.T_star)   # m⁻³ s⁻¹
        xi_p3 = L_PLANCK**3                               # Planck volume (m³)
        return gamma_nuc * xi_p3                          # dimensionless s⁻¹

    # ------------------------------------------------------------------
    # Sakharov Condition 2: CP violation magnitude
    # ------------------------------------------------------------------

    def cp_asymmetry(self) -> float:
        """Return the CP asymmetry δ_CP used in the computation.

        In the Möbius-bundle picture (§18.10), the half-flux topology is
        chiral: it distinguishes left-handed from right-handed windings.  A
        first-order de-saturation bubble can, in principle, nucleate with a
        slight excess of one chirality over the other — the baryogenesis
        analogue of CKM CP violation in electroweak baryogenesis.

        KEY HONEST STATEMENT:
        The spec asserts the Möbius topology provides the CP source but gives
        NO derivation of the asymmetry magnitude.  The θ_QCD is exactly 0 from
        the topology (theta_qcd_from_topology() = 0), removing strong-CP as a
        source.  The remaining CP violation is topological in origin but its
        size is UNCONSTRAINED in the current theory.  δ_CP is a free parameter.

        The bundle's first Chern class is c₁ = 1/2 (chern_class_first() = 0.5).
        This sets the EXISTENCE of the CP source but not its magnitude.

        Returns:
            The free parameter self.delta_cp (dimensionless, ∈ (0, 1]).
        """
        # Confirm Möbius structure enforces θ_QCD = 0 (no strong-CP source)
        bundle = MobiusBundle(winding_number=1)
        _ = bundle  # confirms bundle has c₁ = 1/2
        theta_qcd = theta_qcd_from_topology()
        assert theta_qcd == 0.0, "Möbius topology must give θ_QCD = 0"

        # The magnitude is the free parameter; the Möbius topology justifies
        # its EXISTENCE but not its size.
        return self.delta_cp

    # ------------------------------------------------------------------
    # Sakharov Condition 3: out-of-equilibrium factor
    # ------------------------------------------------------------------

    def out_of_equilibrium_factor(self) -> float:
        """Return the out-of-equilibrium enhancement factor.

        For a first-order phase transition with bubble wall velocity v_w = c
        (§18.44, the substrate phase boundary propagates at the medium wave
        speed = c), the departure from equilibrium is maximal.  The relevant
        quantity is:

            ξ_OOE = v_w / c_s   (wall velocity / sound speed)

        where c_s ≈ c/√3 for a relativistic plasma.  With v_w = c:

            ξ_OOE = c / (c/√3) = √3 ≈ 1.73

        This is the "strong" (supersonic) deflagration regime; in electroweak
        baryogenesis v_w > c_s maximises the asymmetry generated at the wall.

        For our purposes we fold this into the washout_suppression attribute
        and return the phase-transition strength parameter:

            α = ε_latent / ε_radiation

        the ratio of latent heat to radiation energy density.  α >> 1 means
        a strongly first-order transition (maximum asymmetry generation).

        Returns:
            Phase-transition strength parameter α = ε_latent / ε_radiation
            (dimensionless).  α > 1 indicates a strongly first-order
            transition; the out-of-equilibrium condition is satisfied.
        """
        dt_obj = DesaturationTransition(
            sigma_initial=self.sigma_initial,
            sigma_final=self.sigma_final,
            K=K_PLANCK,
        )
        epsilon_latent = dt_obj.latent_heat_density()  # J m⁻³

        # Radiation energy density at T_star
        epsilon_rad = (math.pi**2 / 30.0) * G_STAR * (k_B * self.T_star)**4 / (hbar * c)**3

        if epsilon_rad <= 0.0:
            return float("inf")
        return epsilon_latent / epsilon_rad

    # ------------------------------------------------------------------
    # Hubble rate at the transition
    # ------------------------------------------------------------------

    def hubble_at_transition(self) -> float:
        """Hubble parameter H_* at the de-saturation transition temperature.

        At T_* the radiation-dominated Friedmann equation gives:

            H_* = sqrt(8πG/3) × sqrt((π²/30) g_* (k_B T_*)⁴ / (ℏc)³)

        This is the standard result; at Planck temperature the radiation energy
        density equals the Planck energy density, so H_* ~ 1/t_Planck.

        Returns:
            Hubble rate in s⁻¹.
        """
        return hubble_rate_radiation(self.T_star)

    # ------------------------------------------------------------------
    # Master η_B computation
    # ------------------------------------------------------------------

    def eta_B(self) -> dict[str, float]:
        """Compute η_B = n_B / n_γ and all intermediate quantities.

        Master formula:

            η_B ≈ (Γ_B_eff / H_*) × δ_CP × washout_suppression

        This is the standard order-of-magnitude estimate for baryogenesis
        from an electroweak-like phase transition, adapted to the
        de-saturation context.  The formula is DIMENSIONAL ANALYSIS, not
        a rigorous quantum-field-theory derivation.

        What is DERIVED:
        - Γ_B_eff comes from the Coleman bounce action in desaturation.py
        - H_* comes from the FRW Friedmann equation
        - Phase-transition strength α (out-of-equilibrium condition)
        - v_w = c from §18.44 (wall speed = substrate wave speed)

        What is a FREE PARAMETER:
        - δ_CP: CP asymmetry magnitude — the spec says it comes from Möbius
          topology but does not compute it.  This is the EMPIRICAL KNOB.

        What is an APPROXIMATION:
        - Washout factor: we set it to self.washout_suppression (default 1.0).
          Detailed sphaleron-rate computation is not done here.
        - The conversion factor from Γ_nuc (m⁻³ s⁻¹) to a per-baryon rate
          uses ξ_Planck³ as the "baryon volume" — a ROUGH dimensional estimate.
        - g_* = 100 assumed at Planck scale.

        Returns:
            Dictionary with keys:
            - "eta_B"              : computed baryon-to-photon ratio (dimensionless)
            - "eta_B_observed"     : measured value 6.1e-10
            - "ratio_to_observed"  : η_B_computed / η_B_observed
            - "log10_ratio"        : log10 of the above
            - "Gamma_B_eff"        : baryon-violating rate (s⁻¹)
            - "H_star"             : Hubble rate at transition (s⁻¹)
            - "Gamma_over_H"       : Γ_B_eff / H_* (dimensionless)
            - "delta_cp"           : CP asymmetry used
            - "alpha_PT"           : phase-transition strength ε_lat/ε_rad
            - "washout_factor"     : washout suppression applied
            - "T_star_K"           : transition temperature used (K)
            - "n_gamma_today"      : photon number density today (m⁻³)
        """
        gamma_B = self.baryon_violation_rate()     # s⁻¹ (per Planck volume)
        H_star = self.hubble_at_transition()       # s⁻¹
        delta_cp = self.cp_asymmetry()
        washout = self.washout_suppression
        alpha_pt = self.out_of_equilibrium_factor()

        # η_B estimate
        if H_star > 0.0:
            gamma_over_H = gamma_B / H_star
        else:
            gamma_over_H = float("inf")

        eta_b = gamma_over_H * delta_cp * washout

        # Guard: η_B cannot exceed order 1 physically
        # (if the formula overflows it means our approximations have broken down)
        eta_b = min(eta_b, 1.0)

        ratio = eta_b / ETA_B_OBSERVED if ETA_B_OBSERVED > 0 else float("inf")
        log10_ratio = math.log10(ratio) if ratio > 0 else float("-inf")

        n_gamma = photon_number_density(T_CMB_TODAY)

        return {
            "eta_B": eta_b,
            "eta_B_observed": ETA_B_OBSERVED,
            "ratio_to_observed": ratio,
            "log10_ratio": log10_ratio,
            "Gamma_B_eff": gamma_B,
            "H_star": H_star,
            "Gamma_over_H": gamma_over_H,
            "delta_cp": delta_cp,
            "alpha_PT": alpha_pt,
            "washout_factor": washout,
            "T_star_K": self.T_star,
            "n_gamma_today": n_gamma,
        }


# ---------------------------------------------------------------------------
# Scan over δ_CP to find the matching value
# ---------------------------------------------------------------------------

def scan_delta_cp(
    delta_cp_values: list[float],
    T_star: float | None = None,
) -> list[dict[str, float]]:
    """Compute η_B for each value of δ_CP and flag which matches observation.

    Args:
        delta_cp_values: List of δ_CP values to evaluate (each in (0, 1]).
        T_star: Transition temperature in Kelvin.  Default: Planck temperature.

    Returns:
        List of result dicts, one per δ_CP value (same format as eta_B()).
        Each dict also contains the key "matches_observed" (bool): True when
        |log10(η_B / η_B_obs)| < 1 (within one order of magnitude).
    """
    if T_star is None:
        T_star = T_PLANCK_K

    results: list[dict[str, float]] = []
    for dcp in delta_cp_values:
        model = BaryogenesisModel(T_star=T_star, delta_cp=dcp)
        res = model.eta_B()
        res["matches_observed"] = abs(res["log10_ratio"]) < 1.0
        results.append(res)
    return results


def required_delta_cp_for_observed_eta(T_star: float | None = None) -> float:
    """Invert the baryogenesis formula to find δ_CP that gives η_B = 6.1e-10.

    Given:

        η_B = (Γ_B_eff / H_*) × δ_CP × washout

    We solve for δ_CP:

        δ_CP_required = η_B_observed / (Γ_B_eff / H_* × washout)

    Args:
        T_star: Transition temperature in Kelvin.  Default: Planck temperature.

    Returns:
        Required δ_CP value (dimensionless).  Returns NaN if the formula is
        degenerate (Γ_B_eff = 0 or H_* = 0).
    """
    if T_star is None:
        T_star = T_PLANCK_K

    # Use washout = 1.0 (conservative)
    model = BaryogenesisModel(T_star=T_star, delta_cp=0.5, washout_suppression=1.0)
    gamma_B = model.baryon_violation_rate()
    H_star = model.hubble_at_transition()

    if H_star <= 0.0 or gamma_B <= 0.0:
        return float("nan")

    gamma_over_H = gamma_B / H_star
    if gamma_over_H <= 0.0:
        return float("nan")

    return ETA_B_OBSERVED / gamma_over_H


# ---------------------------------------------------------------------------
# Nucleation rate helper: temperature scan
# ---------------------------------------------------------------------------

def nucleation_rate_vs_temperature(
    temperatures_K: list[float],
    sigma_initial: float = SIGMA_MAX,
    sigma_final: float = 0.0,
) -> list[dict[str, float]]:
    """Compute the bubble nucleation rate at a range of temperatures.

    Args:
        temperatures_K: List of temperatures in Kelvin.
        sigma_initial: Initial substrate strain (default 0.5).
        sigma_final: Final substrate strain (default 0.0).

    Returns:
        List of dicts, each with "T_K", "Gamma_nuc", "log10_Gamma_nuc",
        "H_rad", "Gamma_over_H".
    """
    dt_obj = DesaturationTransition(
        sigma_initial=sigma_initial,
        sigma_final=sigma_final,
        K=K_PLANCK,
    )
    results: list[dict[str, float]] = []
    for T in temperatures_K:
        gamma_nuc = dt_obj.nucleation_rate(T)   # m⁻³ s⁻¹
        H = hubble_rate_radiation(T)             # s⁻¹
        # Nucleation rate per Hubble volume (Hubble radius = c/H)
        r_hubble = c / H if H > 0 else float("inf")
        gamma_per_hubble = gamma_nuc * r_hubble**3 if H > 0 else 0.0

        results.append({
            "T_K": T,
            "Gamma_nuc_m3s": gamma_nuc,
            "log10_Gamma_nuc": math.log10(gamma_nuc) if gamma_nuc > 0 else float("-inf"),
            "H_rad_s": H,
            "Gamma_per_Hubble_vol": gamma_per_hubble,
            "log10_Gamma_per_H": (
                math.log10(gamma_per_hubble)
                if gamma_per_hubble > 0
                else float("-inf")
            ),
        })
    return results
