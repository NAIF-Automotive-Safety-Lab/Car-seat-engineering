"""Chemical reaction kinetics from the substrate framework — paired GEOMETRY + SIM.

Substrate ontology of "reaction"
--------------------------------
In the B3 substrate picture every covalent bond is a coherent strain knot
in the stiff medium.  A *chemical reaction* is the substrate's traversal
from one knot configuration (the reactants) over a cumulative-torque
saddle (the transition state, TS) into another knot configuration (the
products).  The activation barrier E_a is the cumulative torque needed
to pull the substrate strain pattern through the saddle, and the
pre-exponential factor A is the substrate "attempt frequency"
ν ~ k_B T / h on the bond-breaking degree of freedom.

This module wraps the textbook Arrhenius / Eyring forms in two paired
classes:

  * ``ReactionGeometry``       — substrate strain landscape from
    reactants → TS → products.  Stores the energies of the three
    stationary points and produces a smooth one-dimensional reaction
    coordinate ξ ∈ [0, 1] with the canonical V(ξ) potential
    (low quartic with a bias term so the heights match exactly).

  * ``ReactionKineticsSimulator`` — Arrhenius and Eyring rate
    constants, catalysis (lowered E_a) and helper utilities (half-life,
    rate-ratio between catalysed and uncatalysed pathways, ln k vs 1/T
    Arrhenius plot, derived ΔH‡, ΔS‡).

Public helpers
--------------
``arrhenius_rate(A, E_a, T)``            — k = A · exp(−E_a / R T)
``eyring_rate(deltaG_dagger, T, kappa)``  — k = κ · (k_B T / h) · exp(−ΔG‡/R T)
``catalysed_rate_ratio(E_a, E_a_cat, T)`` — k_cat / k = exp(ΔE_a / R T)

References
----------
* Arrhenius (1889)              — temperature dependence of k.
* Eyring (1935)                 — transition-state theory.
* Atkins, Physical Chemistry, 11th ed., chs 21–22 (collision and TST).
* IUPAC compendium of chemical kinetics (Diels-Alder, SN2, H₂+I₂).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import math

import numpy as np
from numpy.typing import NDArray
from scipy import constants as C


# ---------------------------------------------------------------------------
# Module-level reference values (literature E_a for canonical reactions)
# ---------------------------------------------------------------------------

# Activation energies (kJ/mol) and pre-factors (1/s or M^-1 s^-1) collated
# from physical-chemistry textbooks and IUPAC kinetic databases.  Used both
# in the tests and as defaults for the visualization driver.
LITERATURE: Dict[str, Dict[str, float]] = {
    # Bimolecular SN2: CH3Br + OH- -> CH3OH + Br- (gas/aq phase).
    # IUPAC / Atkins: E_a ~ 90 kJ/mol, A ~ 2.0e11 M^-1 s^-1.
    "SN2_methyl_bromide": {
        "E_a_kJ_per_mol": 90.0,
        "A": 2.0e11,
        "delta_E_kJ_per_mol": -75.0,   # exothermic (bond-energy bookkeeping)
        "label": "CH3Br + OH- -> CH3OH + Br-",
    },
    # H2 + I2 -> 2 HI (Bodenstein 1894): E_a ~ 165 kJ/mol, A ~ 1.0e11 M^-1 s^-1.
    "H2_I2_to_2HI": {
        "E_a_kJ_per_mol": 165.0,
        "A": 1.0e11,
        "delta_E_kJ_per_mol": -9.5,
        "label": "H2 + I2 -> 2 HI",
    },
    # Diels-Alder (butadiene + ethene -> cyclohexene): E_a ~ 115 kJ/mol,
    # A ~ 1.0e6 M^-1 s^-1 (Houk computational, IUPAC).
    "Diels_Alder_butadiene_ethene": {
        "E_a_kJ_per_mol": 115.0,
        "A": 1.0e6,
        "delta_E_kJ_per_mol": -170.0,
        "label": "butadiene + ethene -> cyclohexene (Diels-Alder)",
    },
}


# ---------------------------------------------------------------------------
# Module-level functional helpers (used by tests, viz, and class methods)
# ---------------------------------------------------------------------------

def arrhenius_rate(A: float, E_a_J_per_mol: float, T: float) -> float:
    """Arrhenius rate constant: k = A · exp(−E_a / (R · T)).

    Parameters
    ----------
    A : pre-exponential factor (same units as k, e.g. 1/s or M^-1 s^-1).
    E_a_J_per_mol : activation energy in J/mol.
    T : absolute temperature in K (must be > 0).
    """
    if T <= 0.0:
        raise ValueError("Temperature must be positive (Kelvin).")
    R = float(C.R)
    return float(A * math.exp(-E_a_J_per_mol / (R * T)))


def eyring_rate(
    delta_G_dagger_J_per_mol: float,
    T: float,
    kappa: float = 1.0,
) -> float:
    """Eyring (transition-state) rate constant.

        k = κ · (k_B T / h) · exp(−ΔG‡ / (R T))

    The substrate ontology assigns ν_attempt = k_B T / h to the
    bond-breaking attempt frequency on the saddle; κ is the
    transmission coefficient (κ ≤ 1 is typical; κ = 1 → no recrossing).
    """
    if T <= 0.0:
        raise ValueError("Temperature must be positive (Kelvin).")
    R = float(C.R)
    nu_attempt = (C.Boltzmann * T) / C.Planck   # 1/s  (same for any reaction at fixed T)
    return float(kappa * nu_attempt * math.exp(-delta_G_dagger_J_per_mol / (R * T)))


def catalysed_rate_ratio(E_a_J_per_mol: float,
                         E_a_cat_J_per_mol: float,
                         T: float) -> float:
    """Ratio k_cat / k_uncat assuming the same A.

    k_cat / k = exp((E_a − E_a_cat) / (R T))

    For a 30 kJ/mol drop at 300 K this gives a factor ~ 10^5 — exactly
    the catalytic enhancement that defines "catalysis works".
    """
    if T <= 0.0:
        raise ValueError("Temperature must be positive (Kelvin).")
    R = float(C.R)
    return float(math.exp((E_a_J_per_mol - E_a_cat_J_per_mol) / (R * T)))


def half_life_first_order(k: float) -> float:
    """First-order half-life t_{1/2} = ln(2) / k."""
    if k <= 0.0:
        raise ValueError("Rate constant k must be positive.")
    return float(math.log(2.0) / k)


# ---------------------------------------------------------------------------
# 1. ReactionGeometry — substrate strain landscape: reactants -> TS -> products
# ---------------------------------------------------------------------------

@dataclass
class ReactionGeometry:
    """One-dimensional reaction-coordinate landscape of the substrate.

    The substrate strain energy along the reaction coordinate ξ ∈ [0, 1]
    has three stationary points: reactants at ξ = 0, transition state
    at ξ ≈ 0.5, and products at ξ = 1.  The energies are stored in the
    same units (kJ/mol by convention).

    The smooth potential V(ξ) is reconstructed by a quartic + cubic
    polynomial that is constructed to interpolate the three heights
    exactly:

        V(ξ) = E_R + (E_TS − E_R) · sin²(π ξ / 2 · 2)   (heuristic basis)
              + a smooth bias linear in ξ that reaches E_P − E_R at ξ = 1.

    The quartic basis used here is a normalised double-Morse-type
    barrier; see ``coordinate_grid`` and ``potential`` below for the
    explicit form.  The geometry exposes ``activation_energy``,
    ``reaction_energy`` and ``barrier_position``.

    Attributes
    ----------
    label : human-readable name of the reaction
    energy_reactants_kJ : energy of the reactant well (reference)
    energy_ts_kJ        : energy of the transition state (top of barrier)
    energy_products_kJ  : energy of the product well
    """

    label: str = "generic A -> B"
    energy_reactants_kJ: float = 0.0
    energy_ts_kJ: float = 100.0
    energy_products_kJ: float = -50.0

    # ------------------------------------------------------------------ basic numbers
    def activation_energy_kJ(self) -> float:
        """Forward activation energy: E_a = E_TS − E_R."""
        return float(self.energy_ts_kJ - self.energy_reactants_kJ)

    def activation_energy_reverse_kJ(self) -> float:
        """Reverse activation energy: E_a^rev = E_TS − E_P."""
        return float(self.energy_ts_kJ - self.energy_products_kJ)

    def reaction_energy_kJ(self) -> float:
        """Net reaction energy: ΔE = E_P − E_R (negative ⇒ exothermic)."""
        return float(self.energy_products_kJ - self.energy_reactants_kJ)

    def barrier_position(self) -> float:
        """Reaction-coordinate position of the TS: always 0.5 in this model."""
        return 0.5

    # ------------------------------------------------------------------ smooth potential
    def coordinate_grid(self, n: int = 401) -> NDArray[np.float64]:
        """Reaction coordinate ξ ∈ [0, 1] with n grid points."""
        if n < 3:
            raise ValueError("Need at least 3 grid points.")
        return np.linspace(0.0, 1.0, n)

    def potential(self, xi: NDArray[np.float64] | float) -> NDArray[np.float64]:
        """Substrate strain energy V(ξ) along the reaction coordinate.

        The barrier is an exact-interpolation polynomial:

            V(ξ) = E_R · (1 − ξ)² · (1 + 2ξ)
                  + E_P · ξ² · (3 − 2ξ)
                  + 16 · (E_TS − ½(E_R + E_P)) · ξ² · (1 − ξ)²

        The first two terms are the cubic Hermite blend that connects
        E_R at ξ = 0 to E_P at ξ = 1 with zero slope at both ends.
        The third term is a centred bump that adds exactly
        (E_TS − ½(E_R + E_P)) at ξ = ½ and decays to zero at the wells.

        This guarantees V(0) = E_R, V(1) = E_P, V(½) = E_TS, and that
        the wells are local minima (zero derivative at ξ = 0 and 1).
        """
        x = np.asarray(xi, dtype=np.float64)
        E_R = self.energy_reactants_kJ
        E_TS = self.energy_ts_kJ
        E_P = self.energy_products_kJ

        # Cubic Hermite blend with zero slope at both wells.
        blend = E_R * (1.0 - x) ** 2 * (1.0 + 2.0 * x) \
            + E_P * x ** 2 * (3.0 - 2.0 * x)
        # Centred bump that is zero at the wells, peaks at ξ = 1/2 with
        # value 1, hence multiplier (E_TS - mean of wells).
        bump_height = E_TS - 0.5 * (E_R + E_P)
        bump = 16.0 * bump_height * x ** 2 * (1.0 - x) ** 2
        return blend + bump

    # ------------------------------------------------------------------ derivatives & curvatures
    def force(self, xi: NDArray[np.float64] | float) -> NDArray[np.float64]:
        """−dV/dξ along the reaction coordinate (substrate restoring torque)."""
        x = np.asarray(xi, dtype=np.float64)
        # Numerical derivative via finite differences for simplicity.
        h = 1e-5
        V_plus = self.potential(np.clip(x + h, 0.0, 1.0))
        V_minus = self.potential(np.clip(x - h, 0.0, 1.0))
        return -(V_plus - V_minus) / (2.0 * h)

    def barrier_curvature(self) -> float:
        """Imaginary-frequency curvature ω‡² ∝ d²V/dξ² at ξ = 1/2.

        For the centred-bump term the second derivative at ξ = ½ is
        −8 · bump_height (in kJ/mol per unit ξ²).  A negative number
        signals a saddle (imaginary frequency) — exactly what TST
        requires.
        """
        bump_height = self.energy_ts_kJ - 0.5 * (
            self.energy_reactants_kJ + self.energy_products_kJ
        )
        # Hermite blend has zero second derivative at xi=0.5 by construction.
        # Bump term: 16 · h · ξ²(1-ξ)² has d²/dξ² = 32 · h · (1 − 6ξ + 6ξ²).
        # At ξ = 0.5: 32 · h · (1 − 3 + 1.5) = 32 · h · (-0.5) = −16 h.
        return float(-16.0 * bump_height)

    # ------------------------------------------------------------------ summary
    def summary(self) -> Dict[str, float]:
        return {
            "label": self.label,
            "E_a_forward_kJ": self.activation_energy_kJ(),
            "E_a_reverse_kJ": self.activation_energy_reverse_kJ(),
            "delta_E_kJ": self.reaction_energy_kJ(),
            "barrier_position_xi": self.barrier_position(),
            "barrier_curvature_kJ_per_xi2": self.barrier_curvature(),
        }


# ---------------------------------------------------------------------------
# 2. ReactionKineticsSimulator — Arrhenius / Eyring / catalysis on a geometry
# ---------------------------------------------------------------------------

@dataclass
class ReactionKineticsSimulator:
    """Arrhenius + Eyring rate evaluator pinned to a ``ReactionGeometry``.

    The simulator stores a pre-exponential factor A (in 1/s for
    unimolecular or M^-1 s^-1 for bimolecular reactions) and inherits
    the activation energy from the attached ``ReactionGeometry``.
    Catalysis is modelled as a *new* effective E_a (E_a − ΔE_cat) with
    the same A — exactly the mechanism of an enzyme or a heterogeneous
    catalyst that lowers the saddle without changing the reactant or
    product wells.

    Attributes
    ----------
    geometry : ``ReactionGeometry`` instance with E_a, ΔE built in.
    A_prefactor : Arrhenius pre-exponential (units depend on order).
    transmission_kappa : Eyring transmission coefficient κ ∈ (0, 1].
    catalyst_drop_kJ : amount the catalyst lowers E_a by (kJ/mol).
    """

    geometry: ReactionGeometry = field(default_factory=ReactionGeometry)
    A_prefactor: float = 1.0e11
    transmission_kappa: float = 1.0
    catalyst_drop_kJ: float = 0.0

    # ------------------------------------------------------------------ Arrhenius
    def arrhenius_k(self, T: float, with_catalyst: bool = False) -> float:
        """Arrhenius rate constant at temperature T.

        Optionally apply the catalyst (lower E_a by ``catalyst_drop_kJ``).
        """
        E_a_kJ = self.geometry.activation_energy_kJ()
        if with_catalyst:
            E_a_kJ = max(0.0, E_a_kJ - self.catalyst_drop_kJ)
        E_a_J = E_a_kJ * 1.0e3
        return arrhenius_rate(self.A_prefactor, E_a_J, T)

    # ------------------------------------------------------------------ Eyring (TST)
    def eyring_k(
        self,
        T: float,
        delta_S_dagger_J_per_K_per_mol: float = 0.0,
        with_catalyst: bool = False,
    ) -> float:
        """Eyring rate from the geometry's E_a (treated as ΔH‡).

        ΔG‡ = ΔH‡ − T · ΔS‡; we identify ΔH‡ ≈ E_a (Arrhenius limit
        for condensed-phase or low-T gas-phase reactions).  When
        ``with_catalyst`` is true, the activation enthalpy is lowered
        by ``catalyst_drop_kJ``.
        """
        E_a_kJ = self.geometry.activation_energy_kJ()
        if with_catalyst:
            E_a_kJ = max(0.0, E_a_kJ - self.catalyst_drop_kJ)
        delta_H_dagger_J = E_a_kJ * 1.0e3
        delta_G_dagger_J = (
            delta_H_dagger_J - T * delta_S_dagger_J_per_K_per_mol
        )
        return eyring_rate(delta_G_dagger_J, T, self.transmission_kappa)

    # ------------------------------------------------------------------ catalysis
    def catalysed_speedup(self, T: float) -> float:
        """k_cat / k_uncat = exp(ΔE_a_cat / (R T))."""
        if self.catalyst_drop_kJ <= 0.0:
            return 1.0
        E_a_kJ = self.geometry.activation_energy_kJ()
        E_a_cat_kJ = max(0.0, E_a_kJ - self.catalyst_drop_kJ)
        return catalysed_rate_ratio(E_a_kJ * 1e3, E_a_cat_kJ * 1e3, T)

    # ------------------------------------------------------------------ Arrhenius plot data
    def arrhenius_curve(
        self,
        T_min: float = 250.0,
        T_max: float = 800.0,
        n: int = 80,
        with_catalyst: bool = False,
    ) -> Dict[str, NDArray[np.float64]]:
        """ln k vs 1/T grid suitable for the classic Arrhenius plot.

        Returns a dict with keys ``inv_T`` (1/K), ``ln_k``,
        ``T`` (K) and ``k``.  The slope of ln k vs 1/T is −E_a / R.
        """
        T = np.linspace(T_min, T_max, n)
        k = np.array(
            [self.arrhenius_k(float(t), with_catalyst=with_catalyst) for t in T]
        )
        return {
            "T": T,
            "inv_T": 1.0 / T,
            "k": k,
            "ln_k": np.log(k),
        }

    # ------------------------------------------------------------------ extracted slope
    def extracted_E_a_kJ(
        self,
        T_min: float = 250.0,
        T_max: float = 800.0,
        n: int = 80,
        with_catalyst: bool = False,
    ) -> float:
        """Recover E_a from the slope of ln k vs 1/T.

        slope = −E_a / R → E_a = −R · slope (in J/mol).  The function
        returns the value in kJ/mol so the round-trip with the input
        is exact in the ideal-Arrhenius limit (constant A).
        """
        curve = self.arrhenius_curve(T_min, T_max, n, with_catalyst)
        slope, _intercept = np.polyfit(curve["inv_T"], curve["ln_k"], 1)
        E_a_J = -float(C.R) * float(slope)
        return E_a_J / 1.0e3

    # ------------------------------------------------------------------ summary
    def summary(self, T: float = 298.15) -> Dict[str, float]:
        """One-shot snapshot at temperature T."""
        return {
            "label": self.geometry.label,
            "T_K": T,
            "E_a_kJ": self.geometry.activation_energy_kJ(),
            "k_arrhenius": self.arrhenius_k(T, with_catalyst=False),
            "k_eyring": self.eyring_k(T, with_catalyst=False),
            "k_arrhenius_catalysed": self.arrhenius_k(T, with_catalyst=True),
            "catalyst_speedup": self.catalysed_speedup(T),
            "half_life_s": (
                half_life_first_order(self.arrhenius_k(T))
                if self.arrhenius_k(T) > 0.0 else float("inf")
            ),
        }


# ---------------------------------------------------------------------------
# 3. Convenience constructors for the canonical reactions used in the visuals
# ---------------------------------------------------------------------------

def make_reaction(name: str, catalyst_drop_kJ: float = 0.0) -> ReactionKineticsSimulator:
    """Build a ``ReactionKineticsSimulator`` for one of the LITERATURE entries."""
    if name not in LITERATURE:
        raise KeyError(
            f"Unknown reaction '{name}'. "
            f"Available: {sorted(LITERATURE.keys())}"
        )
    spec = LITERATURE[name]
    geom = ReactionGeometry(
        label=spec["label"],
        energy_reactants_kJ=0.0,
        energy_ts_kJ=spec["E_a_kJ_per_mol"],
        energy_products_kJ=spec["delta_E_kJ_per_mol"],
    )
    return ReactionKineticsSimulator(
        geometry=geom,
        A_prefactor=spec["A"],
        catalyst_drop_kJ=catalyst_drop_kJ,
    )


def all_canonical_reactions(catalyst_drop_kJ: float = 0.0) -> List[ReactionKineticsSimulator]:
    """Return one simulator per LITERATURE reaction (in the order they were defined)."""
    return [make_reaction(n, catalyst_drop_kJ) for n in LITERATURE.keys()]
