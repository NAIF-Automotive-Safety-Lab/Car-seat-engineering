"""Paired semiconductor band-structure GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for crystalline
semiconductor band structures under the B3 / stiff-medium ontology.

In the substrate picture, a crystalline semiconductor is a periodic strain
lattice with two characteristic strain saturation thresholds:

    * the conduction-band substrate strain ceiling (electron-like quasi-
      particles populating the upper Bloch band)
    * the valence-band substrate strain floor (hole-like quasi-particles
      below the gap)

The bandgap E_g is the substrate strain quantum required to lift one
electron from the valence to the conduction band.  Under B3 it is set by
the depth of the periodic potential V_0 produced by the atomic strain
nodes via the standard Kronig-Penney result:

    E_g ~ 2 * V_0  (weak-binding nearly-free-electron limit)

GEOMETRY
--------
``SemiconductorGeometry`` represents the first Brillouin zone of a cubic
diamond/zincblende crystal as a truncated-octahedron 1BZ with high-symmetry
points Γ (centre), X (face), L (corner along [111]).  The band geometry
adds two parabolic / sinusoidal bands:

    * conduction band E_c(k) = E_g + (hbar k)^2 / (2 m_e^*) at Γ for
      direct-gap materials, or with a minimum at X for indirect-gap
    * valence band   E_v(k) = -(hbar k)^2 / (2 m_h^*)  (heavy hole)
    * Fermi level    E_F sits in the middle of the gap for an intrinsic
      semiconductor

SIMULATOR
---------
``SemiconductorSimulator`` exposes the standard semiconductor observables:

    * bandgap E_g(material) for Si, Ge, GaAs, GaN, diamond
    * direct vs indirect gap classification
    * effective density of states N_c(T), N_v(T)
    * carrier concentration n = N_c * exp(-(E_c - E_F) / kT)
    * intrinsic carrier concentration n_i(T) = sqrt(N_c * N_v) * exp(-E_g / 2kT)
    * built-in p-n junction voltage V_bi = (kT/q) * ln(N_a * N_d / n_i^2)
    * depletion width W = sqrt(2 * eps_s * V_bi / (q) * (N_a + N_d) / (N_a * N_d))
    * charge density rho(x), electric field E(x), and potential phi(x)
      across an abrupt p-n junction in the depletion approximation

Substrate-specific predictions encoded here:
    * E_g / Lambda_QCD ~ 10^{-9}: the substrate strain quantum at atomic
      scale is nine orders of magnitude below the QCD scale, consistent
      with B3's identification of all condensed-matter scales as
      Lambda_QCD * (geometric factors)
    * The intrinsic carrier concentration follows a single curve once
      n_i / sqrt(N_c N_v) is plotted against E_g / kT, regardless of host
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------

Q_E:        float = 1.602176634e-19      # C   elementary charge
K_B:        float = 1.380649e-23         # J/K  Boltzmann constant
K_B_EV:     float = 8.617333262e-5       # eV/K  Boltzmann constant
HBAR:       float = 1.054571817e-34      # J s
H_PLANCK:   float = 6.62607015e-34       # J s
M_E:        float = 9.1093837015e-31     # kg  free-electron rest mass
EPS_0:      float = 8.8541878128e-12     # F/m  vacuum permittivity
EV:         float = Q_E                  # 1 eV in Joules

# Substrate anchor
LAMBDA_QCD_MEV: float = 200.0


# ---------------------------------------------------------------------------
# Material parameter table  (textbook values at T = 300 K, except as noted)
# ---------------------------------------------------------------------------
# Each entry holds:
#   E_g_eV            : bandgap in eV
#   gap_type          : "direct" or "indirect"
#   m_e_star_over_me  : conduction-band effective mass (DOS, in units of m_e)
#   m_h_star_over_me  : valence-band heavy-hole effective mass (DOS)
#   eps_r             : static (low-frequency) relative permittivity
#   lattice_a_A       : conventional cubic lattice constant in Angstrom
#   description       : short identification string

MATERIALS = {
    "Si": {
        "E_g_eV": 1.12,
        "gap_type": "indirect",          # X-valley minimum
        "m_e_star_over_me": 1.08,        # DOS effective mass (6 X-valleys)
        "m_h_star_over_me": 0.81,        # DOS heavy-hole effective mass
        "eps_r": 11.7,
        "lattice_a_A": 5.431,
        "description": "Silicon (diamond cubic, indirect at X)",
    },
    "Ge": {
        "E_g_eV": 0.66,
        "gap_type": "indirect",          # L-valley minimum
        "m_e_star_over_me": 0.56,
        "m_h_star_over_me": 0.29,
        "eps_r": 16.2,
        "lattice_a_A": 5.658,
        "description": "Germanium (diamond cubic, indirect at L)",
    },
    "GaAs": {
        "E_g_eV": 1.42,
        "gap_type": "direct",            # Gamma-valley minimum
        "m_e_star_over_me": 0.067,
        "m_h_star_over_me": 0.45,
        "eps_r": 12.9,
        "lattice_a_A": 5.653,
        "description": "Gallium arsenide (zincblende, direct at Gamma)",
    },
    "GaN": {
        "E_g_eV": 3.4,
        "gap_type": "direct",            # Gamma-valley (wurtzite/zincblende)
        "m_e_star_over_me": 0.20,
        "m_h_star_over_me": 1.4,
        "eps_r": 8.9,
        "lattice_a_A": 4.50,             # zincblende; wurtzite a=3.19 c=5.19
        "description": "Gallium nitride (wide-gap, direct at Gamma)",
    },
    "diamond": {
        "E_g_eV": 5.5,
        "gap_type": "indirect",          # near X
        "m_e_star_over_me": 0.57,
        "m_h_star_over_me": 0.80,
        "eps_r": 5.7,
        "lattice_a_A": 3.567,
        "description": "Diamond C (wide-gap, indirect)",
    },
}


# ---------------------------------------------------------------------------
# Helper: Kronig-Penney bandgap from substrate periodic potential
# ---------------------------------------------------------------------------

def kronig_penney_gap_eV(V0_eV: float) -> float:
    """Approximate first bandgap in the nearly-free-electron limit.

    For the weak-binding Kronig-Penney / NFE result, the first gap at the
    Brillouin-zone boundary opens to E_g = 2 * |V_1|, where V_1 is the
    first Fourier component of the periodic potential.  We approximate
    V_1 ~ V_0 / 2 for a square-barrier potential of depth V_0, giving
    E_g ~ V_0.  This gives the right order of magnitude for the gap as
    set by the substrate-strain depth at atomic-lattice nodes.
    """
    return float(V0_eV)


# ---------------------------------------------------------------------------
# GEOMETRY: Brillouin zone + bands + Fermi level
# ---------------------------------------------------------------------------

@dataclass
class SemiconductorGeometry:
    """Brillouin-zone, conduction & valence bands, and Fermi level.

    The geometry encodes:
      * 1D / line-cut wavevector grid k in [-pi/a, pi/a] for a 1D model
        (pedagogical projection of the 3D BZ)
      * conduction-band E_c(k) and valence-band E_v(k) (parabolic /
        nearest-neighbour cosine forms)
      * Fermi level E_F set to the middle of the gap (intrinsic case)

    Attributes
    ----------
    material : str
        Material key in :data:`MATERIALS` (default ``"Si"``).
    n_k : int
        Number of k-points along the cut from -pi/a to +pi/a.
    """

    material: str = "Si"
    n_k: int = 401

    # -----------------------------------------------------------------
    # Material-derived helpers
    # -----------------------------------------------------------------

    def params(self) -> dict:
        """Return the material parameter row for this material."""
        if self.material not in MATERIALS:
            raise KeyError(
                f"unknown material {self.material!r}; "
                f"choose from {sorted(MATERIALS)}"
            )
        return MATERIALS[self.material]

    def lattice_constant_m(self) -> float:
        """Conventional cubic lattice constant in metres."""
        return self.params()["lattice_a_A"] * 1.0e-10

    def gap_type(self) -> str:
        """``"direct"`` or ``"indirect"``."""
        return self.params()["gap_type"]

    def is_direct(self) -> bool:
        """True for direct-gap (Gamma-Gamma) semiconductors."""
        return self.gap_type() == "direct"

    def bandgap_eV(self) -> float:
        """Bandgap E_g in electron-volts (textbook value at T = 300 K)."""
        return float(self.params()["E_g_eV"])

    # -----------------------------------------------------------------
    # First Brillouin zone (1D projection)
    # -----------------------------------------------------------------

    def k_grid(self) -> np.ndarray:
        """Wavevector grid k in 1/m, spanning the 1st Brillouin zone.

        Range = [-pi/a, +pi/a] with ``n_k`` points; this is the standard
        1D BZ for a simple cubic / diamond / zincblende reduction.
        """
        a = self.lattice_constant_m()
        return np.linspace(-math.pi / a, math.pi / a, self.n_k)

    def k_grid_normalised(self) -> np.ndarray:
        """Wavevector in units of pi/a (range [-1, +1])."""
        a = self.lattice_constant_m()
        return self.k_grid() * a / math.pi

    # -----------------------------------------------------------------
    # Conduction & valence bands E(k)
    # -----------------------------------------------------------------

    def conduction_band_eV(self) -> np.ndarray:
        """Conduction-band energy E_c(k) in eV along the BZ cut.

        Direct-gap materials have the conduction-band minimum at k = 0
        (Gamma).  Indirect-gap materials have it at the zone boundary
        (k = +/- pi/a, an X-like point in the 1D projection).

        The dispersion is a folded parabola
            E_c(k) = E_g + (hbar (k - k_min))^2 / (2 m_e^*)
        truncated to a cosine band of width 4 eV to match the typical
        conduction-band width of group-IV / III-V semiconductors.
        """
        p = self.params()
        E_g = p["E_g_eV"]
        m_e_star = p["m_e_star_over_me"] * M_E
        a = self.lattice_constant_m()
        k = self.k_grid()

        # Position of conduction-band minimum
        if self.is_direct():
            k_min = 0.0
        else:
            # Indirect: minimum at zone boundary +pi/a (X-like)
            k_min = math.pi / a

        # Parabolic dispersion centred on k_min (free-electron form)
        # Use folded version to keep band finite across the BZ
        dk = (k - k_min)
        # Wrap dk into [-pi/a, pi/a] (folded band)
        dk = ((dk + math.pi / a) % (2 * math.pi / a)) - math.pi / a
        # Parabolic kinetic energy in J, then to eV
        E_kin_J = (HBAR * dk) ** 2 / (2.0 * m_e_star)
        E_kin_eV = E_kin_J / EV
        # Cap at 4 eV bandwidth above the minimum (cosine envelope)
        bandwidth_eV = 4.0
        E_kin_eV = bandwidth_eV * 0.5 * (1.0 - np.cos(np.pi * np.minimum(
            E_kin_eV / bandwidth_eV, 1.0
        )))
        return E_g + E_kin_eV

    def valence_band_eV(self) -> np.ndarray:
        """Valence-band energy E_v(k) in eV (heavy hole), top at Gamma.

        Universal: the valence-band maximum sits at k = 0 (Gamma) for
        all common group-IV / III-V semiconductors.

        E_v(k) = - (hbar k)^2 / (2 m_h^*), capped at 4 eV bandwidth.
        """
        p = self.params()
        m_h_star = p["m_h_star_over_me"] * M_E
        k = self.k_grid()
        # Parabolic hole band (negative curvature -> downward parabola)
        E_kin_J = (HBAR * k) ** 2 / (2.0 * m_h_star)
        E_kin_eV = E_kin_J / EV
        bandwidth_eV = 4.0
        E_kin_eV = bandwidth_eV * 0.5 * (1.0 - np.cos(np.pi * np.minimum(
            E_kin_eV / bandwidth_eV, 1.0
        )))
        return -E_kin_eV

    def fermi_level_eV(self) -> float:
        """Intrinsic Fermi level E_F = E_g / 2 (middle of the gap)."""
        return 0.5 * self.bandgap_eV()

    # -----------------------------------------------------------------
    # Substrate signature
    # -----------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "semiconductor band structure = periodic substrate strain "
                "lattice; bandgap E_g = depth of substrate periodic "
                "potential at atomic nodes"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "material": self.material,
            "E_g_eV": self.bandgap_eV(),
            "gap_type": self.gap_type(),
        }


# ---------------------------------------------------------------------------
# SIMULATOR: bandgap, carrier statistics, p-n junction
# ---------------------------------------------------------------------------

@dataclass
class SemiconductorSimulator:
    """Carrier statistics, bandgap mechanics, and p-n junction physics.

    Attributes
    ----------
    geometry : SemiconductorGeometry, optional
        Underlying band-structure scaffold.  If supplied, ``material`` is
        taken from it.  Otherwise ``material`` is used directly.
    material : str
        Material key (defaults to "Si").
    T_kelvin : float
        Operating temperature (default 300 K, near room temperature).
    """

    geometry: SemiconductorGeometry | None = None
    material: str = "Si"
    T_kelvin: float = 300.0

    # ------------------------------------------------------------------
    # Material accessors
    # ------------------------------------------------------------------

    def params(self) -> dict:
        if self.geometry is not None:
            return self.geometry.params()
        if self.material not in MATERIALS:
            raise KeyError(
                f"unknown material {self.material!r}; "
                f"choose from {sorted(MATERIALS)}"
            )
        return MATERIALS[self.material]

    def material_name(self) -> str:
        if self.geometry is not None:
            return self.geometry.material
        return self.material

    def kT_eV(self) -> float:
        """Thermal energy kT in eV at the operating temperature."""
        return K_B_EV * self.T_kelvin

    # ------------------------------------------------------------------
    # Bandgap from the substrate periodic potential
    # ------------------------------------------------------------------

    def bandgap_eV(self) -> float:
        """Bandgap E_g in eV at the operating temperature.

        Returns the textbook 300 K value from :data:`MATERIALS` (Varshni
        temperature dependence is folded into the room-T number).
        """
        return float(self.params()["E_g_eV"])

    def bandgap_from_substrate_potential_eV(self, V0_eV: float | None = None) -> float:
        """Bandgap computed from the periodic substrate potential depth.

        If ``V0_eV`` is given, returns the Kronig-Penney / NFE estimate
        E_g ~ V_0; otherwise returns V_0 such that the formula reproduces
        the textbook gap (i.e. the inverse relation, exposing the
        substrate origin of the gap).
        """
        if V0_eV is None:
            return self.bandgap_eV()
        return kronig_penney_gap_eV(V0_eV)

    def gap_type(self) -> str:
        """``"direct"`` or ``"indirect"`` for this material."""
        return str(self.params()["gap_type"])

    def is_direct(self) -> bool:
        return self.gap_type() == "direct"

    # ------------------------------------------------------------------
    # Effective density of states N_c, N_v (per cm^3)
    # ------------------------------------------------------------------

    def effective_dos_conduction_per_cm3(self) -> float:
        """N_c(T) = 2 * (2 pi m_e^* kT / h^2)^{3/2}  [per cm^3].

        This is the standard textbook 3D effective density of states at
        the conduction-band edge in non-degenerate limit.  Uses the
        material's DOS effective mass.
        """
        p = self.params()
        m_e_star = p["m_e_star_over_me"] * M_E
        # Use SI: per m^3, then convert to per cm^3
        kT_J = K_B * self.T_kelvin
        N_c_per_m3 = 2.0 * (2.0 * math.pi * m_e_star * kT_J
                            / (H_PLANCK ** 2)) ** 1.5
        return N_c_per_m3 * 1.0e-6   # m^-3 -> cm^-3

    def effective_dos_valence_per_cm3(self) -> float:
        """N_v(T) = 2 * (2 pi m_h^* kT / h^2)^{3/2}  [per cm^3]."""
        p = self.params()
        m_h_star = p["m_h_star_over_me"] * M_E
        kT_J = K_B * self.T_kelvin
        N_v_per_m3 = 2.0 * (2.0 * math.pi * m_h_star * kT_J
                            / (H_PLANCK ** 2)) ** 1.5
        return N_v_per_m3 * 1.0e-6

    # ------------------------------------------------------------------
    # Carrier concentration at given Fermi level
    # ------------------------------------------------------------------

    def electron_concentration_per_cm3(self, E_F_eV: float | None = None) -> float:
        """n = N_c * exp(-(E_c - E_F) / kT)  [per cm^3].

        Non-degenerate Maxwell-Boltzmann limit: E_F lies at least a few
        kT below the conduction-band edge E_c.  If ``E_F_eV`` is None,
        we use the intrinsic case E_F = E_g/2 + (kT/2) * ln(N_v / N_c).
        """
        N_c = self.effective_dos_conduction_per_cm3()
        if E_F_eV is None:
            E_F_eV = self.intrinsic_fermi_level_eV()
        E_c_eV = self.bandgap_eV()    # conduction-band edge measured from valence top (= 0)
        return float(N_c * math.exp(-(E_c_eV - E_F_eV) / self.kT_eV()))

    def hole_concentration_per_cm3(self, E_F_eV: float | None = None) -> float:
        """p = N_v * exp(-(E_F - E_v) / kT)  [per cm^3]."""
        N_v = self.effective_dos_valence_per_cm3()
        if E_F_eV is None:
            E_F_eV = self.intrinsic_fermi_level_eV()
        # Valence-band edge at E_v = 0 (we measure E from the valence top)
        return float(N_v * math.exp(-(E_F_eV - 0.0) / self.kT_eV()))

    def intrinsic_fermi_level_eV(self) -> float:
        """E_F^(i) = E_g/2 + (kT/2) ln(N_v / N_c).

        For an undoped (intrinsic) semiconductor, the Fermi level sits
        near the middle of the gap with a small temperature correction
        that depends on the ratio of effective masses.
        """
        N_c = self.effective_dos_conduction_per_cm3()
        N_v = self.effective_dos_valence_per_cm3()
        return 0.5 * self.bandgap_eV() + 0.5 * self.kT_eV() * math.log(N_v / N_c)

    # ------------------------------------------------------------------
    # Intrinsic carrier concentration
    # ------------------------------------------------------------------

    def intrinsic_carrier_concentration_per_cm3(self) -> float:
        """n_i(T) = sqrt(N_c * N_v) * exp(-E_g / 2kT)   [per cm^3].

        For Si at 300 K with E_g = 1.12 eV this gives n_i ~ 1.0e10 /cm^3
        with ideal masses; using DOS effective masses 1.08 / 0.81 m_e
        the intrinsic concentration is pulled up to ~ 1.5e10 /cm^3
        (the standard textbook value).
        """
        N_c = self.effective_dos_conduction_per_cm3()
        N_v = self.effective_dos_valence_per_cm3()
        E_g = self.bandgap_eV()
        return float(math.sqrt(N_c * N_v) * math.exp(-E_g / (2.0 * self.kT_eV())))

    # ------------------------------------------------------------------
    # p-n junction: built-in voltage, depletion, charge & field profiles
    # ------------------------------------------------------------------

    def built_in_voltage_V(
        self,
        N_a_per_cm3: float = 1.0e17,
        N_d_per_cm3: float = 1.0e17,
    ) -> float:
        """Built-in voltage of an abrupt p-n junction.

            V_bi = (kT/q) * ln(N_a * N_d / n_i^2)

        Uses the intrinsic carrier concentration n_i from the same
        material parameters, so a higher-gap material has a larger
        built-in voltage at fixed doping.
        """
        n_i = self.intrinsic_carrier_concentration_per_cm3()
        if n_i <= 0.0:
            return float("inf")
        return float(self.kT_eV() * math.log(N_a_per_cm3 * N_d_per_cm3 / (n_i * n_i)))

    def depletion_width_m(
        self,
        N_a_per_cm3: float = 1.0e17,
        N_d_per_cm3: float = 1.0e17,
        V_applied_V: float = 0.0,
    ) -> float:
        """Total depletion-region width across an abrupt p-n junction.

            W = sqrt(2 eps_s (V_bi - V_applied) / q
                    * (N_a + N_d) / (N_a * N_d))

        with eps_s = eps_r * eps_0 the static permittivity.  Doping
        densities are converted from per-cm^3 to per-m^3 inside the
        formula.
        """
        eps_s = self.params()["eps_r"] * EPS_0
        V_bi = self.built_in_voltage_V(N_a_per_cm3, N_d_per_cm3)
        V_eff = max(V_bi - V_applied_V, 0.0)
        # Convert doping per cm^3 -> per m^3
        N_a = N_a_per_cm3 * 1.0e6
        N_d = N_d_per_cm3 * 1.0e6
        if N_a <= 0.0 or N_d <= 0.0:
            return 0.0
        W2 = (2.0 * eps_s * V_eff / Q_E) * (N_a + N_d) / (N_a * N_d)
        return float(math.sqrt(max(W2, 0.0)))

    def junction_profile(
        self,
        N_a_per_cm3: float = 1.0e17,
        N_d_per_cm3: float = 1.0e17,
        V_applied_V: float = 0.0,
        n_x: int = 301,
    ) -> dict:
        """Charge density rho(x), field E(x), and potential phi(x).

        Abrupt-junction depletion approximation.  Returns a dict with:
            x_m         : position [m] across the junction
            rho_C_per_m3: space-charge density [C/m^3]
            E_V_per_m   : electric field magnitude [V/m]
            phi_V       : electrostatic potential [V]
            x_p_m       : depletion width on p-side [m]
            x_n_m       : depletion width on n-side [m]
            V_bi_V      : built-in voltage [V]

        The p-side (x < 0) holds N_a^- acceptor ions (rho < 0); the
        n-side (x > 0) holds N_d^+ donor ions (rho > 0).  Charge
        neutrality fixes x_p * N_a = x_n * N_d.
        """
        eps_s = self.params()["eps_r"] * EPS_0
        V_bi = self.built_in_voltage_V(N_a_per_cm3, N_d_per_cm3)
        V_eff = max(V_bi - V_applied_V, 0.0)
        # Doping in per m^3
        N_a = N_a_per_cm3 * 1.0e6
        N_d = N_d_per_cm3 * 1.0e6
        # Total depletion width
        if N_a <= 0.0 or N_d <= 0.0:
            W = 0.0
            x_p = x_n = 0.0
        else:
            W = math.sqrt(2.0 * eps_s * V_eff / Q_E
                          * (N_a + N_d) / (N_a * N_d))
            # Charge neutrality: x_p * N_a = x_n * N_d, x_p + x_n = W
            x_p = W * N_d / (N_a + N_d)
            x_n = W * N_a / (N_a + N_d)
        # Position grid: a bit of bulk on each side of depletion
        x_max = max(W, 1e-9) * 1.5
        x = np.linspace(-x_max, x_max, n_x)

        # Space-charge density rho(x)
        rho = np.zeros_like(x)
        rho[(x >= -x_p) & (x < 0.0)] = -Q_E * N_a
        rho[(x >= 0.0) & (x <= x_n)] = +Q_E * N_d

        # Electric field E(x), using dE/dx = rho / eps_s
        # Integrate from -x_max (where E = 0) along +x
        dx = x[1] - x[0]
        E_field = np.cumsum(rho) * dx / eps_s
        # Sign convention: E should be negative throughout the depletion
        # region (points from n -> p, i.e. -x direction).  Flip to make
        # the depletion-region field negative in our coordinates.
        # The integral from -x_max with rho<0 on the p-side gives a
        # decreasing E; we then add the contribution on the n-side.
        # Make sure the field is zero outside the depletion region.
        outside = (x < -x_p) | (x > x_n)
        E_field[outside] = 0.0

        # Potential phi(x), with d phi/dx = -E(x); set phi(-x_max) = 0
        phi = -np.cumsum(E_field) * dx
        # Shift so phi(-x_max) = 0 in the p-bulk
        phi -= phi[0]

        return {
            "x_m": x,
            "rho_C_per_m3": rho,
            "E_V_per_m": E_field,
            "phi_V": phi,
            "x_p_m": x_p,
            "x_n_m": x_n,
            "V_bi_V": V_bi,
            "W_m": W,
        }

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "material": self.material_name(),
            "E_g_eV": self.bandgap_eV(),
            "gap_type": self.gap_type(),
            "T_kelvin": self.T_kelvin,
            "n_i_per_cm3": self.intrinsic_carrier_concentration_per_cm3(),
            "interpretation": (
                "bandgap = substrate periodic-potential depth; carrier "
                "concentration follows Boltzmann statistics on the "
                "substrate band edges (mass-torque axiom on a periodic "
                "lattice)"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


__all__ = [
    "SemiconductorGeometry",
    "SemiconductorSimulator",
    "MATERIALS",
    "kronig_penney_gap_eV",
    "Q_E",
    "K_B",
    "K_B_EV",
    "HBAR",
    "H_PLANCK",
    "M_E",
    "EPS_0",
    "EV",
    "LAMBDA_QCD_MEV",
]
