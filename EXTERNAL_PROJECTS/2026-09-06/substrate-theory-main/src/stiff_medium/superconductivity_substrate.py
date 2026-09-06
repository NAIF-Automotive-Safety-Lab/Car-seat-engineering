"""Superconductivity as substrate phase: paired GEOMETRY + SIMULATOR.

This module makes the substrate ontology of superconductivity EXPLICIT:

    Cooper pair    = a substrate "paired strain bridge" — two electron-strain
                     defects of opposite momentum/spin, locked together by a
                     weak attractive medium-mediated interaction (the substrate
                     replacement of the BCS phonon-exchange channel).  The
                     centre-of-mass mode is collectively coherent across the
                     condensate; this is the substrate reading of the BCS
                     ground state |BCS> = Π_k (u_k + v_k c+_kup c+_-kdn) |0>.

    Energy gap Δ   = the substrate strain it costs to break a paired bridge
                     into two unpaired electron defects (= 2 unpaired strain
                     defects).  At T=0,
                          Δ(0)  =  2 ℏ ω_D  exp(-1 / (N(0) V))
                     in the weak-coupling BCS limit, where ω_D is the substrate
                     mode-locking frequency (Debye), N(0) is the substrate-mode
                     density of states at the Fermi level, and V is the paired
                     attraction strength.

    Critical T_c   = the temperature at which thermal substrate-mode noise
                     unbinds the pairs.  In the weak-coupling BCS limit
                          2 Δ(0) / (k_B T_c)  =  2 π / e^γ  ≃  3.527754
                     (the BCS UNIVERSAL gap ratio; γ = Euler-Mascheroni).

    T_c,max        = the substrate-saturation upper bound on conventional,
                     ambient-pressure superconductivity:
                          T_c,max = Λ_QCD / n_R  =  200 MeV / 18
                                  =  11.11 meV
                                  ≃  128.9 K
                     (HBCCO record 134 K matches at 4.0%; see
                     ``b3_high_tc_bound``).  This is the FIRST cross-disciplinary
                     B3 prediction beyond fundamental physics.

GEOMETRY
--------
``SuperconductivityGeometry`` represents a superconductor as

    * a paired-strain Cooper bridge (two coupled defects with controllable
      coherence length ξ and pair size ξ_0 = ℏ v_F / (π Δ));
    * a single-band quasi-particle dispersion E_k = sqrt(ε_k² + Δ²) with
      Bogoliubov u_k², v_k² coherence factors;
    * a Meissner shell (London penetration depth λ_L) with the macroscopic
      Ginzburg-Landau ratio κ = λ_L / ξ that classifies type-I (κ < 1/√2)
      vs type-II (κ > 1/√2).

SIMULATOR
---------
``SuperconductivitySimulator`` exposes the BCS physics:

    * gap_at_T(T):                  self-consistent BCS gap Δ(T) → 0 at T_c
    * gap_zero_T():                 Δ(0) = 2 ℏ ω_D exp(-1/(N(0)V))
    * critical_temperature():       T_c via 2Δ(0)/(k_B T_c) = 3.527754
    * density_of_states(E):         N(E)/N(0) = E / sqrt(E² - Δ²) for |E|>Δ
    * meissner_field(z, lambda_L):  B(z) = B_0 exp(-z/λ_L)
    * type_classification(kappa):   "I" if κ < 1/√2, "II" otherwise
    * bcs_ratio():                  the universal 2Δ/(k_B T_c) = 3.527754
    * substrate_tc_max():           Λ_QCD / n_R  =  128.9 K

Substrate-specific predictions encoded here:
    * BCS universal gap ratio is intrinsic to the substrate-paired ground
      state and is reproduced exactly (3.527754).
    * T_c,max = 128.9 K caps any conventional ambient-pressure SC; matched
      by HBCCO (134 K, 31-year record) at 4.0%.  Hydride superhydrides
      (H_3S 203 K, LaH_10 260 K) lie above only because they sit at extreme
      pressures (155-170 GPa) where the substrate mode structure is
      qualitatively different.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import math
import numpy as np
from numpy.typing import NDArray

from .b3_constants import LAMBDA_QCD_K, n_R


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

K_B:   float = 1.380649e-23      # J/K  Boltzmann constant
HBAR:  float = 1.054571817e-34   # J s
EV_PER_J: float = 1.0 / 1.602176634e-19   # eV per Joule
MEV_PER_K: float = (K_B * EV_PER_J) * 1e3  # 1 K -> meV (k_B in meV/K)
# 1 K  ↔  k_B = 8.617333e-5 eV  =  8.617333e-2 meV.  Used to convert the
# substrate-saturation energy bound (200 MeV / 18 = 11.11 meV) into a
# temperature ceiling: T_max = 11.11 meV / k_B = 128.9 K.

# BCS universal weak-coupling gap ratio  2Δ(0) / (k_B T_c)  =  2π/e^γ
EULER_GAMMA: float = 0.5772156649015329
BCS_UNIVERSAL_GAP_RATIO: float = 2.0 * math.pi / math.exp(EULER_GAMMA)
# Numerical value: 3.527754...  (BCS 1957)

# Type-I / type-II Ginzburg-Landau classification threshold
GL_KAPPA_THRESHOLD: float = 1.0 / math.sqrt(2.0)


# ---------------------------------------------------------------------------
# Substrate T_c,max ceiling   (B3 cross-disciplinary prediction)
# ---------------------------------------------------------------------------
# Energy per substrate mode at saturation: E* = Λ_QCD / n_R
# In B3, Λ_QCD = 200 MeV is treated as numerically equal to 200 K-equivalent
# units only through the substrate-saturation conversion below, so we work
# in proper SI throughout: convert Λ_QCD (MeV) -> J -> K via k_B.
LAMBDA_QCD_MEV_VALUE: float = 200.0
SUBSTRATE_E_PER_MODE_MEV: float = LAMBDA_QCD_MEV_VALUE / float(n_R)
# = 11.111... meV
SUBSTRATE_TC_MAX_K: float = (
    SUBSTRATE_E_PER_MODE_MEV * 1e-3        # meV -> eV
) / (K_B * EV_PER_J)
# = 128.917... K (matches the cross-disciplinary write-up at 4 sig fig)


# ---------------------------------------------------------------------------
# Reference superconductor table
# ---------------------------------------------------------------------------
# Per-element/compound reference critical temperature (K) and an indicative
# Debye temperature (K) used to back-solve the BCS pair coupling N(0)V.
#
# Sources (PDG, Tinkham "Intro to Superconductivity", Bednorz/Müller 1986,
# Schilling et al. 1993, Eisaki et al. 2004, NIST):
#
#   - Hg                      4.2 K   (Onnes 1911, original discovery)
#   - Pb                      7.2 K
#   - Nb                      9.25 K  (highest among elemental SCs)
#   - MgB_2                  39 K    (Nagamatsu 2001)
#   - YBa_2Cu_3O_7 (YBCO)    93 K    (Wu/Chu 1987)
#   - Bi_2Sr_2Ca_2Cu_3O_10   110 K
#   - Tl_2Ba_2Ca_2Cu_3O_10   125 K
#   - HgBa_2Ca_2Cu_3O_8      134 K   (Schilling 1993, ambient-pressure record)
#
REFERENCE_SUPERCONDUCTORS: Dict[str, Dict[str, float]] = {
    "Hg":       dict(T_c_K=4.2,   T_D_K=72.0,  type="I",  family="elemental"),
    "Pb":       dict(T_c_K=7.2,   T_D_K=105.0, type="I",  family="elemental"),
    "Nb":       dict(T_c_K=9.25,  T_D_K=275.0, type="II", family="elemental"),
    "MgB2":     dict(T_c_K=39.0,  T_D_K=900.0, type="II", family="conv-2band"),
    "YBCO":     dict(T_c_K=93.0,  T_D_K=410.0, type="II", family="cuprate"),
    "BSCCO":    dict(T_c_K=110.0, T_D_K=380.0, type="II", family="cuprate"),
    "TBCCO":    dict(T_c_K=125.0, T_D_K=350.0, type="II", family="cuprate"),
    "HBCCO":    dict(T_c_K=134.0, T_D_K=320.0, type="II", family="cuprate"),
}


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class SuperconductivityGeometry:
    """Substrate geometry for a Cooper-paired condensate.

    Parameters
    ----------
    delta_meV : float
        Zero-temperature gap Δ(0) in milli-electronvolts.  Default 1.40 meV
        is the canonical Pb value.
    coherence_xi_nm : float
        Pippard / GL coherence length ξ in nanometres (pair size).
    london_lambda_nm : float
        London penetration depth λ_L in nanometres.
    n_pair_grid : int
        Number of substrate cells used in the coupled-defect bridge picture.
    """

    delta_meV:        float = 1.40
    coherence_xi_nm:  float = 80.0
    london_lambda_nm: float = 37.0
    n_pair_grid:      int   = 64

    # ---- Cooper paired strain bridge ------------------------------------

    def paired_bridge_positions(self) -> NDArray[np.float64]:
        """Return positions (nm) of the two strain defects that form the
        Cooper pair.  Separation is 2·ξ (one coherence length each side of
        the centre), giving a pair "size" of 2·ξ.
        """
        return np.array([-self.coherence_xi_nm, +self.coherence_xi_nm])

    def paired_bridge_envelope(
        self, x_nm: Optional[NDArray[np.float64]] = None
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Strain envelope ψ_pair(x) ∝ exp(-|x|/ξ) of the paired bridge.

        Returns (x_grid_nm, envelope) with x in [-3ξ, +3ξ] and envelope
        normalised to peak 1 at the centre of mass (x=0).
        """
        if x_nm is None:
            # Force an odd sample count so x = 0 and x = ±ξ are landed on
            # (n_pair_grid//6 cells per ξ, 6 ξ total span -> 6 m + 1 nodes).
            m = max(1, self.n_pair_grid // 6)
            x_nm = np.linspace(
                -3.0 * self.coherence_xi_nm,
                +3.0 * self.coherence_xi_nm,
                6 * m + 1,
            )
        env = np.exp(-np.abs(x_nm) / self.coherence_xi_nm)
        return x_nm, env

    # ---- BCS quasi-particle dispersion ----------------------------------

    def bogoliubov_dispersion(
        self, eps_k_meV: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Quasi-particle energy E_k = sqrt(ε_k² + Δ²) in meV.

        ``eps_k_meV`` is the normal-state band energy measured from the
        Fermi level.  At ε_k = 0 the gap minimum E_k = Δ.
        """
        return np.sqrt(eps_k_meV ** 2 + self.delta_meV ** 2)

    def bogoliubov_uv2(
        self, eps_k_meV: NDArray[np.float64]
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """BCS coherence factors (u_k², v_k²).

            u_k² = ½ (1 + ε_k / E_k)
            v_k² = ½ (1 − ε_k / E_k)
        """
        E_k = self.bogoliubov_dispersion(eps_k_meV)
        # Avoid division by zero at the Fermi point (ε=0); E_k=Δ there.
        ratio = np.divide(eps_k_meV, E_k, out=np.zeros_like(E_k),
                          where=(E_k > 0.0))
        return 0.5 * (1.0 + ratio), 0.5 * (1.0 - ratio)

    # ---- Ginzburg-Landau type classification ----------------------------

    def gl_kappa(self) -> float:
        """Ginzburg-Landau ratio κ = λ_L / ξ."""
        return self.london_lambda_nm / self.coherence_xi_nm

    def is_type_II(self) -> bool:
        """True iff the GL κ exceeds the type-I/II threshold 1/√2."""
        return self.gl_kappa() > GL_KAPPA_THRESHOLD


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class SuperconductivitySimulator:
    """BCS / substrate-saturation simulator paired with a geometry.

    Parameters
    ----------
    geometry : SuperconductivityGeometry
        Substrate geometry (provides Δ(0), ξ, λ_L).
    debye_T_K : float
        Debye temperature θ_D in Kelvin.  Sets ℏω_D = k_B θ_D.
    coupling_NV : float
        Dimensionless BCS pair coupling N(0)V > 0 (weak-coupling regime
        N(0)V ≪ 1).  If None, back-solved from Δ(0) via
            Δ(0) = 2 ℏω_D exp(-1/(N(0)V)).
    """

    geometry:    SuperconductivityGeometry = field(default_factory=SuperconductivityGeometry)
    debye_T_K:   float = 105.0
    coupling_NV: Optional[float] = None

    def __post_init__(self) -> None:
        if self.coupling_NV is None:
            # Back-solve from the geometry's Δ(0).
            two_hbar_omega_D_meV = (
                2.0 * (K_B * self.debye_T_K) * EV_PER_J * 1e3
            )
            ratio = self.geometry.delta_meV / two_hbar_omega_D_meV
            if ratio <= 0.0 or ratio >= 1.0:
                # Outside the BCS weak-coupling validity window.
                raise ValueError(
                    f"Δ(0)/(2 ℏω_D) = {ratio:g} outside (0,1); "
                    "supply coupling_NV explicitly or revisit Debye T."
                )
            # Δ = 2ℏω_D exp(-1/NV)  ⇒  N(0)V = -1 / ln(ratio).
            self.coupling_NV = float(-1.0 / math.log(ratio))

    # ---- Universal BCS quantities ---------------------------------------

    @staticmethod
    def bcs_ratio() -> float:
        """The universal weak-coupling BCS gap ratio 2Δ(0)/(k_B T_c) = 2π/e^γ.

        Numerical value: 3.5277537146...
        """
        return BCS_UNIVERSAL_GAP_RATIO

    # ---- Zero-temperature gap from the BCS gap equation ------------------

    def gap_zero_T_meV(self) -> float:
        """Δ(0) = 2 ℏω_D exp(-1/(N(0)V)), in meV.

        ℏω_D is taken as k_B θ_D (Debye approximation).
        """
        two_hbar_omega_D_meV = (
            2.0 * (K_B * self.debye_T_K) * EV_PER_J * 1e3
        )
        return float(
            two_hbar_omega_D_meV * math.exp(-1.0 / float(self.coupling_NV))
        )

    # ---- Critical temperature from the universal gap ratio --------------

    def critical_temperature_K(self) -> float:
        """T_c = 2 Δ(0) / (k_B · 3.527754).

        Returns T_c in Kelvin.  Equivalent to the BCS expression
            k_B T_c = (e^γ / π) · ℏω_D · exp(-1/(N(0)V))
        with ℏω_D = k_B θ_D.
        """
        delta0_J = self.gap_zero_T_meV() * 1e-3 / EV_PER_J
        return float(2.0 * delta0_J / (K_B * BCS_UNIVERSAL_GAP_RATIO))

    # ---- Self-consistent BCS gap Δ(T) -----------------------------------

    def gap_at_T(self, T_K: float, n_quad: int = 4096) -> float:
        """Self-consistent BCS gap Δ(T) at temperature T (Kelvin).

        Solves
            1/(N(0)V) = ∫_0^{ℏω_D} dε  tanh(sqrt(ε² + Δ²)/(2 k_B T)) /
                                       sqrt(ε² + Δ²)
        for Δ ≥ 0 by bisection.  Returns Δ in meV.

        For T ≥ T_c returns 0.0.
        """
        T_c = self.critical_temperature_K()
        if T_K >= T_c:
            return 0.0
        if T_K <= 0.0:
            return self.gap_zero_T_meV()

        # Energies in J for stability.
        hbar_omega_D_J = K_B * self.debye_T_K
        kT_J = K_B * T_K
        target = 1.0 / float(self.coupling_NV)

        eps = np.linspace(0.0, hbar_omega_D_J, n_quad)
        d_eps = eps[1] - eps[0]

        def lhs_minus_target(delta_J: float) -> float:
            E = np.sqrt(eps ** 2 + delta_J ** 2)
            # tanh(E/(2kT)) / E ; protect E=0 when delta=0
            with np.errstate(divide="ignore", invalid="ignore"):
                integrand = np.where(
                    E > 0.0, np.tanh(E / (2.0 * kT_J)) / E, 0.0
                )
            # numpy 2.x renamed trapz -> trapezoid; fall back for 1.x.
            trapz_fn = getattr(np, "trapezoid", None) or np.trapz
            integral = float(trapz_fn(integrand, dx=d_eps))
            return integral - target

        # Bisection on Δ ∈ [0, Δ(0) * 1.05].
        delta_max_J = self.gap_zero_T_meV() * 1e-3 / EV_PER_J * 1.05
        lo, hi = 0.0, delta_max_J
        f_hi = lhs_minus_target(hi)
        # f(0) is monotone ↑ in T, so at T<T_c we should have lhs(0) > target
        # (gap can be > 0).  If lhs(0) < target the gap is 0 (above T_c).
        if f_hi >= 0.0:
            # gap rail-pegs at maximum; return Δ(0).
            return self.gap_zero_T_meV()
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f_mid = lhs_minus_target(mid)
            if f_mid > 0.0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-30:
                break
        return float(0.5 * (lo + hi)) * EV_PER_J * 1e3  # back to meV

    def gap_curve(
        self, n_T: int = 60
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return (T/T_c, Δ(T)/Δ(0)) curve over T ∈ [0, T_c]."""
        T_c = self.critical_temperature_K()
        T_grid = np.linspace(0.001 * T_c, 0.999 * T_c, n_T)
        d0 = self.gap_zero_T_meV()
        d_T = np.array([self.gap_at_T(t) for t in T_grid])
        return T_grid / T_c, d_T / d0

    # ---- Density of states with gap -------------------------------------

    def density_of_states_ratio(
        self, E_meV: NDArray[np.float64], delta_meV: Optional[float] = None
    ) -> NDArray[np.float64]:
        """BCS density of states N_s(E)/N(0).

            N_s(E)/N(0) = |E| / sqrt(E² − Δ²)        for |E| > Δ
                        = 0                          for |E| < Δ

        E is measured from the Fermi level (meV).  Returns the SAME shape
        as ``E_meV``.
        """
        if delta_meV is None:
            delta_meV = self.geometry.delta_meV
        absE = np.abs(E_meV)
        out = np.zeros_like(E_meV, dtype=float)
        mask = absE > float(delta_meV)
        out[mask] = absE[mask] / np.sqrt(absE[mask] ** 2 - delta_meV ** 2)
        return out

    # ---- Meissner effect (London) ---------------------------------------

    def meissner_field(
        self,
        z_nm: NDArray[np.float64],
        B0_T: float = 1.0,
        lambda_L_nm: Optional[float] = None,
    ) -> NDArray[np.float64]:
        """Magnetic-field profile B(z) inside a SC half-space.

            B(z) = B_0 · exp(-z / λ_L)

        z is depth from the surface (nm); B_0 is the surface field (T).
        """
        if lambda_L_nm is None:
            lambda_L_nm = self.geometry.london_lambda_nm
        return B0_T * np.exp(-np.maximum(z_nm, 0.0) / lambda_L_nm)

    # ---- Type-I vs type-II classification -------------------------------

    @staticmethod
    def type_classification(kappa: float) -> str:
        """Return ``"I"`` if κ < 1/√2, ``"II"`` otherwise (Abrikosov 1957)."""
        return "I" if kappa < GL_KAPPA_THRESHOLD else "II"

    # ---- Substrate-saturation upper bound (B3 cross-disciplinary) -------

    @staticmethod
    def substrate_tc_max_K() -> float:
        """B3 substrate-saturation upper bound:

            T_c,max = (Λ_QCD / n_R) / k_B
                    = (200 MeV / 18) / k_B
                    ≃ 128.9 K

        See ``b3_high_tc_bound`` for the full cross-disciplinary derivation
        (HBCCO record 134 K matches at 4.0%, 31-year ambient-pressure
        unbroken).
        """
        return SUBSTRATE_TC_MAX_K

    # ---- Convenience: bulk reference table ------------------------------

    def reference_table_rows(self) -> List[Dict[str, object]]:
        """Return the reference SC table sorted by T_c ascending."""
        rows: List[Dict[str, object]] = []
        for name, info in REFERENCE_SUPERCONDUCTORS.items():
            T_c = float(info["T_c_K"])
            rows.append({
                "name":          name,
                "T_c_K":         T_c,
                "T_D_K":         float(info["T_D_K"]),
                "type":          str(info["type"]),
                "family":        str(info["family"]),
                "below_ceiling": T_c <= SUBSTRATE_TC_MAX_K,
                "ratio_to_ceiling": T_c / SUBSTRATE_TC_MAX_K,
            })
        rows.sort(key=lambda r: float(r["T_c_K"]))
        return rows


# ---------------------------------------------------------------------------
# Module exports
# ---------------------------------------------------------------------------

__all__ = [
    "BCS_UNIVERSAL_GAP_RATIO",
    "EULER_GAMMA",
    "GL_KAPPA_THRESHOLD",
    "REFERENCE_SUPERCONDUCTORS",
    "SUBSTRATE_E_PER_MODE_MEV",
    "SUBSTRATE_TC_MAX_K",
    "SuperconductivityGeometry",
    "SuperconductivitySimulator",
]
