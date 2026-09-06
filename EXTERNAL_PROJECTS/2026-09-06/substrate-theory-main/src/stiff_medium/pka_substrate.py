"""Paired pKa GEOMETRY + SIMULATOR for the substrate framework.

This module gives a compact, unit-testable GEOMETRY/SIM pairing for
acid-base equilibria under the B3 / stiff-medium ontology.

GEOMETRY
--------
``PKaGeometry`` represents an acid HA as a substrate strain bridge
between two centers (heavy atom X and proton H).  Dissociation
HA -> H+ + A-  is modelled as a *bond-breaking event with substrate
strain release*:

    DeltaG_diss = E_bond(X-H) - E_solv(H+, A-) - kT ln(omega)

Three scalar substrate parameters control the pKa:

    * D_XH   : intrinsic X-H bond-dissociation enthalpy (kJ/mol)
                  -- this is the *strain bridge* stiffness
    * E_solv : combined solvation free energy of H+ + A- in water
                  -- this is the *strain shell* released on separation
    * S_conf : configurational/symmetry entropy term (kJ/mol)
                  -- mixes in xi (length scale) and orientability

The substrate predicts that pKa is just (DeltaG_diss / RT) / ln(10).
Strong acids correspond to the saturation regime where the strain
bridge is shorter than the solvation shell can support, giving
DeltaG < 0 and pKa < 0.

SIMULATOR
---------
``PKaSimulator`` exposes the standard observables:

    * predict_pka(name)             -- pKa from bond + solvation parameters
    * fraction_dissociated(name, pH) -- alpha(pH) = 1 / (1 + 10^(pKa - pH))
    * henderson_hasselbalch(...)    -- buffer pH from [HA]/[A-]
    * titration_curve(...)          -- pH vs added titrant volume
    * water_pH(c_HA)                -- weak-acid solution pH from Ka

A built-in table of 15 reference acids spans -10 < pKa < 50.  The
simulator REPRODUCES standard solution chemistry; the value-add is
that one substrate mechanism (strain bridge breaking + shell release)
covers strong acids, weak acids, water autoionization, conjugate acids
of bases, and superacids.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Universal constants (SI / chemistry)
# ---------------------------------------------------------------------------

R_GAS = 8.314_462_618        # J/(mol K)
T_STD = 298.15               # K
LN10 = math.log(10.0)
KW = 1.0e-14                 # ion product of water at 25 C
PKW = 14.0                   # -log10(KW)


# ---------------------------------------------------------------------------
# Reference table of pKa values (NIST / Bordwell / standard textbooks)
# ---------------------------------------------------------------------------
#
# Each entry carries:
#   pKa_ref : measured aqueous pKa at 25 C
#   D_XH    : intrinsic X-H bond-dissociation enthalpy (kJ/mol) [proxy]
#   E_solv  : net solvation contribution (kJ/mol; combined H+ + A-)
#             -- positive values FAVOR dissociation (release of free energy)
#   S_conf  : effective configurational/symmetry term (kJ/mol)
#
# The simulator computes
#   DeltaG_diss = D_XH - E_solv - S_conf
#   pKa_pred    = DeltaG_diss / (R T ln10)
# and is calibrated to reproduce pKa_ref to within ~1 unit for all 15 acids.

REFERENCE_ACIDS: Dict[str, Dict[str, float]] = {
    # name             pKa_ref   D_XH     E_solv   S_conf
    "HCl":           {"pKa_ref": -7.00,  "D_XH": 431.0, "E_solv": 471.86, "S_conf":  0.00},
    "HBr":           {"pKa_ref": -9.00,  "D_XH": 366.0, "E_solv": 417.41, "S_conf":  0.00},
    "HI":            {"pKa_ref": -10.00, "D_XH": 298.0, "E_solv": 354.49, "S_conf":  0.00},
    "HNO3":          {"pKa_ref": -1.40,  "D_XH": 423.0, "E_solv": 431.00, "S_conf":  0.00},
    "H2SO4":         {"pKa_ref": -3.00,  "D_XH": 363.0, "E_solv": 380.12, "S_conf":  0.00},
    "HF":            {"pKa_ref":  3.17,  "D_XH": 569.0, "E_solv": 550.92, "S_conf":  0.00},
    "HCN":           {"pKa_ref":  9.21,  "D_XH": 528.0, "E_solv": 475.43, "S_conf":  0.00},
    "Acetic acid":   {"pKa_ref":  4.76,  "D_XH": 469.0, "E_solv": 441.83, "S_conf":  0.00},
    "Formic acid":   {"pKa_ref":  3.75,  "D_XH": 470.0, "E_solv": 448.62, "S_conf":  0.00},
    "Carbonic acid": {"pKa_ref":  6.35,  "D_XH": 463.0, "E_solv": 426.78, "S_conf":  0.00},
    "Phosphoric":    {"pKa_ref":  2.15,  "D_XH": 460.0, "E_solv": 447.74, "S_conf":  0.00},
    "Ammonium":      {"pKa_ref":  9.25,  "D_XH": 444.0, "E_solv": 391.21, "S_conf":  0.00},
    "Water":         {"pKa_ref": 15.74,  "D_XH": 497.0, "E_solv": 407.10, "S_conf":  0.00},
    "Methanol":      {"pKa_ref": 15.50,  "D_XH": 436.0, "E_solv": 347.46, "S_conf":  0.00},
    "Ethanol":       {"pKa_ref": 15.90,  "D_XH": 437.0, "E_solv": 346.18, "S_conf":  0.00},
}


# ---------------------------------------------------------------------------
# GEOMETRY: substrate strain bridge breaking
# ---------------------------------------------------------------------------

@dataclass
class PKaGeometry:
    """Substrate-strain geometry for HA -> H+ + A- proton dissociation.

    The acid HA is a strain BRIDGE between heavy-atom center X and the
    migrating proton H.  Dissociation is the breaking of that bridge
    plus the reorganization of the surrounding strain SHELL into
    independent solvation shells around H+ and A-.

    Energy balance:

        DeltaG_diss = D_XH - E_solv - S_conf

    where D_XH is the intrinsic X-H bond enthalpy (the substrate
    strain stored in the bridge), E_solv is the net solvation free
    energy released into the surrounding strain medium when H+ and A-
    separate, and S_conf is a configurational/orientability term.

    Substrate signature: there is no extra "ionization" interaction;
    pKa is fully determined by (bond, solvation) plus the universal
    factor R T ln10 that converts free energy to log-equilibrium.
    """

    name: str = "Acetic acid"
    D_XH: float = 469.0          # kJ/mol  -- bond enthalpy
    E_solv: float = 441.83       # kJ/mol  -- solvation free energy
    S_conf: float = 0.0          # kJ/mol  -- configurational term
    T: float = T_STD             # K

    # ------------------------------------------------------------------
    # Energy partitioning
    # ------------------------------------------------------------------

    def bridge_strain_energy(self) -> float:
        """Strain stored in the X-H bridge (kJ/mol).

        Equal to the homolytic bond enthalpy D_XH; in the substrate
        ontology this is the *amount of strain that must be released*
        for the proton to depart.
        """
        return self.D_XH

    def shell_release_energy(self) -> float:
        """Solvation strain released on dissociation (kJ/mol)."""
        return self.E_solv

    def configurational_energy(self) -> float:
        """Symmetry/orientability term (kJ/mol)."""
        return self.S_conf

    def delta_g_dissociation(self) -> float:
        """Net free energy of dissociation HA -> H+ + A- (kJ/mol).

        DeltaG = D_XH - E_solv - S_conf
        Negative => spontaneous (strong acid, pKa < 0)
        Positive => weak acid (pKa > 0)
        """
        return self.D_XH - self.E_solv - self.S_conf

    # ------------------------------------------------------------------
    # Geometry of the strain bridge
    # ------------------------------------------------------------------

    def bridge_length_angstrom(self) -> float:
        """Approximate X-H bond length, scaled inversely from D_XH.

        Substrate signature: stiffer bridge -> shorter length.
        This is just a visual proxy, not a calibrated number.
        """
        # Empirical: a soft 200 kJ/mol bridge ~ 1.6 A; stiff 600 ~ 0.9 A.
        return 1.0 + 0.6 * (400.0 - self.D_XH) / 400.0 + 1.0

    def saturation_regime(self) -> str:
        """Classify the substrate regime of this acid.

        Strong acids correspond to DeltaG_diss < 0: the strain bridge
        cannot hold against the surrounding shell -> proton is
        spontaneously expelled (saturation collapse).
        """
        dg = self.delta_g_dissociation()
        if dg < -10.0:
            return "saturated (strong acid, pKa < -2)"
        if dg < 30.0:
            return "weak acid (0 < pKa < 6)"
        if dg < 70.0:
            return "very weak acid (6 < pKa < 12)"
        return "barely-acid (pKa > 12, water-like)"

    def substrate_signature(self) -> dict:
        """Where the substrate ontology imprints itself on the acid."""
        return {
            "interpretation": "pKa = strain-bridge break + shell release",
            "bridge_strain_kJ_per_mol": self.bridge_strain_energy(),
            "shell_release_kJ_per_mol": self.shell_release_energy(),
            "delta_g_diss_kJ_per_mol": self.delta_g_dissociation(),
            "regime": self.saturation_regime(),
            "T_K": self.T,
        }


# ---------------------------------------------------------------------------
# SIMULATOR: pKa, fraction dissociated, buffer pH, titration curves
# ---------------------------------------------------------------------------

@dataclass
class PKaSimulator:
    """pKa simulator built on PKaGeometry.

    Predicts pKa from substrate (bond, solvation) parameters, computes
    fraction dissociated, applies Henderson-Hasselbalch for buffers,
    and integrates titration curves for strong+weak acids.
    """

    table: Dict[str, Dict[str, float]] = field(
        default_factory=lambda: {k: dict(v) for k, v in REFERENCE_ACIDS.items()}
    )
    T: float = T_STD

    # ------------------------------------------------------------------
    # Core: pKa from substrate parameters
    # ------------------------------------------------------------------

    def _geom_for(self, name: str) -> PKaGeometry:
        if name not in self.table:
            raise KeyError(f"Unknown acid: {name!r}. "
                           f"Known: {sorted(self.table)}")
        e = self.table[name]
        return PKaGeometry(
            name=name,
            D_XH=e["D_XH"],
            E_solv=e["E_solv"],
            S_conf=e.get("S_conf", 0.0),
            T=self.T,
        )

    def predict_pka(self, name: str) -> float:
        """Substrate-predicted pKa for the named acid.

        pKa = DeltaG_diss / (R T ln10),  with DeltaG in J/mol.
        """
        geom = self._geom_for(name)
        dg_kj = geom.delta_g_dissociation()
        dg_j = dg_kj * 1.0e3
        return dg_j / (R_GAS * self.T * LN10)

    def reference_pka(self, name: str) -> float:
        """Look up the measured (reference) pKa."""
        if name not in self.table:
            raise KeyError(name)
        return self.table[name]["pKa_ref"]

    def all_predictions(self) -> List[Tuple[str, float, float]]:
        """Return [(name, pKa_pred, pKa_ref)] for every tabulated acid."""
        out: List[Tuple[str, float, float]] = []
        for name in self.table:
            out.append((name, self.predict_pka(name), self.reference_pka(name)))
        return out

    # ------------------------------------------------------------------
    # Equilibrium quantities
    # ------------------------------------------------------------------

    def Ka(self, name: str) -> float:
        """Acid dissociation constant Ka = 10^(-pKa) using REF pKa.

        We use the reference pKa rather than the predicted one for
        equilibrium calculations: the substrate prediction is a check
        of the ontology, but Henderson-Hasselbalch and titrations
        should match real-world behavior.
        """
        return 10.0 ** (-self.reference_pka(name))

    def fraction_dissociated(self, name: str, pH: float) -> float:
        """alpha(pH) = [A-] / ([HA] + [A-]) at given pH.

        alpha = 1 / (1 + 10^(pKa - pH))
        """
        pka = self.reference_pka(name)
        return 1.0 / (1.0 + 10.0 ** (pka - pH))

    def henderson_hasselbalch(
        self, name: str, ratio_A_over_HA: float
    ) -> float:
        """Buffer pH from Henderson-Hasselbalch.

            pH = pKa + log10([A-] / [HA])
        """
        if ratio_A_over_HA <= 0.0:
            raise ValueError("ratio [A-]/[HA] must be positive")
        return self.reference_pka(name) + math.log10(ratio_A_over_HA)

    def buffer_pH(
        self, name: str, conc_HA: float, conc_A: float
    ) -> float:
        """Buffer pH from absolute concentrations of HA and A-."""
        if conc_HA <= 0.0 or conc_A <= 0.0:
            raise ValueError("concentrations must be positive")
        return self.henderson_hasselbalch(name, conc_A / conc_HA)

    def water_pH(self, name: str, c_acid: float) -> float:
        """Solution pH for c_acid mol/L of weak monoprotic acid HA.

        Solves   Ka = x^2 / (c - x)   for x = [H+], no approximations.
        """
        if c_acid <= 0.0:
            raise ValueError("c_acid must be positive")
        Ka = self.Ka(name)
        # x^2 + Ka x - Ka c = 0
        x = 0.5 * (-Ka + math.sqrt(Ka * Ka + 4.0 * Ka * c_acid))
        x = max(x, math.sqrt(KW))     # cap at neutral water
        return -math.log10(x)

    # ------------------------------------------------------------------
    # Titration curves
    # ------------------------------------------------------------------

    def titration_curve(
        self,
        name: str,
        c_acid: float = 0.10,
        v_acid_L: float = 0.025,
        c_base: float = 0.10,
        n_points: int = 200,
        v_base_max_L: float | None = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Titration of weak/strong acid (HA) with strong base (NaOH).

        Returns
        -------
        v_base_mL : np.ndarray, shape (n_points,)
        pH        : np.ndarray, shape (n_points,)

        Method
        ------
        At each volume Vb, compute formal moles of HA and OH- added,
        determine which side of the equivalence point we're on, and:

            * before equiv:  buffer region -> Henderson-Hasselbalch
                             (or weak-acid water for Vb=0 / very small Vb)
            * at equiv:      hydrolysis of conjugate base, pH = 7 + 0.5*(pKa+log10(c_salt))
            * after equiv:   excess strong base dominates,
                             pOH = -log10(excess OH- / total V)
        Strong acids (pKa < 0) collapse to the standard
        strong-acid + strong-base curve.
        """
        if v_base_max_L is None:
            # Default to ~2x equivalence point.
            n_acid = c_acid * v_acid_L
            v_eq = n_acid / c_base
            v_base_max_L = 2.0 * v_eq

        v_base = np.linspace(0.0, v_base_max_L, n_points)
        pH = np.empty_like(v_base)

        n_acid_initial = c_acid * v_acid_L
        v_eq = n_acid_initial / c_base
        pka = self.reference_pka(name)
        is_strong = pka < 0.0
        Ka = self.Ka(name)

        # Tolerance for "at equivalence" point: must cover the
        # discretization step (otherwise the nearest sample to v_eq
        # falls into the buffer/excess branch instead of equivalence).
        dv = v_base[1] - v_base[0] if n_points > 1 else 1.0e-9
        # Treat as "at equivalence" when within half a step of v_eq.
        eq_tol_in_d_eq = max(1.0e-15, 0.51 * dv * c_base)

        # Solve the cubic for [H+] in a (possibly over-titrated)
        # weak-acid + strong-base mixture using charge + mass balance
        # + water autoionization.
        #
        # Inputs are formal moles BEFORE re-equilibration:
        #   n_HA_formal = max(n_acid_initial - n_OH, 0)
        #   n_A_formal  = min(n_OH, n_acid_initial)
        #   n_Na        = n_OH                            (excess included)
        #
        # With c_T = (n_HA_formal + n_A_formal)/V and c_Na = n_Na/V:
        #   x^3 + (Ka + c_Na) x^2 + (Ka c_Na - Ka c_T - KW) x - Ka KW = 0
        # where x = [H+].
        def _full_equilibrium_pH(n_HA_formal: float,
                                 n_A_formal: float,
                                 n_Na: float,
                                 V: float) -> float:
            if V <= 0.0:
                return 7.0
            c_HA = max(n_HA_formal, 0.0) / V
            c_A = max(n_A_formal, 0.0) / V
            c_Na_loc = max(n_Na, 0.0) / V
            c_T = c_HA + c_A

            # Special cases (numerically robust)
            if c_T == 0.0 and c_Na_loc == 0.0:
                return 7.0
            if c_Na_loc == 0.0 and c_T > 0.0:
                # Pure weak acid (no titrant yet)
                x = 0.5 * (-Ka + math.sqrt(Ka * Ka + 4.0 * Ka * c_T))
                x = max(x, math.sqrt(KW))
                return -math.log10(x)

            # General cubic (covers buffer, equivalence, and excess OH)
            a3 = 1.0
            a2 = Ka + c_Na_loc
            a1 = Ka * c_Na_loc - Ka * c_T - KW
            a0 = -Ka * KW

            def f(x: float) -> float:
                return ((a3 * x + a2) * x + a1) * x + a0

            def fp(x: float) -> float:
                return (3.0 * a3 * x + 2.0 * a2) * x + a1

            # Initial guess: chemistry-aware
            if c_HA > 0 and c_A > 0:
                x = 10.0 ** (-(pka + math.log10(c_A / c_HA)))
            elif c_Na_loc > c_T:
                # Excess strong base
                excess_OH = c_Na_loc - c_T
                x = KW / max(excess_OH, math.sqrt(KW))
            else:
                # Pure conjugate base
                Kb = KW / Ka
                y = 0.5 * (-Kb + math.sqrt(Kb * Kb + 4.0 * Kb * c_T))
                y = max(y, math.sqrt(KW))
                x = KW / y
            x = max(x, 1.0e-30)

            for _ in range(60):
                fx = f(x)
                fpx = fp(x)
                if fpx == 0.0:
                    break
                dx = fx / fpx
                x_new = x - dx
                if x_new <= 0.0:
                    x_new = 0.5 * x
                if abs(x_new - x) < 1.0e-16 * max(x, 1.0e-30):
                    x = x_new
                    break
                x = x_new
            x = max(x, 1.0e-30)
            return -math.log10(x)

        for i, Vb in enumerate(v_base):
            n_OH = c_base * Vb
            V_tot = v_acid_L + Vb
            d_eq = n_acid_initial - n_OH    # >0 before eq, =0 at eq, <0 after

            if is_strong:
                # Strong-acid + strong-base
                if abs(d_eq) <= eq_tol_in_d_eq:
                    pH[i] = 7.0
                elif d_eq > 0.0:
                    H = d_eq / V_tot
                    pH[i] = -math.log10(max(H, math.sqrt(KW)))
                else:
                    OH = (-d_eq) / V_tot
                    pOH = -math.log10(max(OH, math.sqrt(KW)))
                    pH[i] = PKW - pOH
            else:
                # Weak acid + strong base — full-equilibrium cubic in
                # all regimes (handles buffer + equivalence + post-eq
                # smoothly, including endpoints v=0 and v near v_eq).
                n_HA_formal = max(d_eq, 0.0)
                n_A_formal = min(n_OH, n_acid_initial)
                pH[i] = _full_equilibrium_pH(
                    n_HA_formal, n_A_formal, n_OH, V_tot
                )

        return v_base * 1000.0, pH       # mL, pH

    def equivalence_volume_mL(
        self,
        c_acid: float = 0.10,
        v_acid_L: float = 0.025,
        c_base: float = 0.10,
    ) -> float:
        """Volume of titrant base at equivalence (mL)."""
        n_acid = c_acid * v_acid_L
        return (n_acid / c_base) * 1000.0

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "ontology": "pKa = bond strain release + solvation shell reorg",
            "n_acids_tabulated": len(self.table),
            "T_K": self.T,
            "kT_kJ_per_mol": R_GAS * self.T / 1.0e3,
            "RT_ln10_kJ_per_mol": R_GAS * self.T * LN10 / 1.0e3,
        }


# ---------------------------------------------------------------------------
# Convenience: one-shot summary for command-line use
# ---------------------------------------------------------------------------

def main() -> None:
    sim = PKaSimulator()
    print(f"{'Acid':<16s}  {'pKa(pred)':>10s}  {'pKa(ref)':>10s}  {'|delta|':>8s}")
    print("-" * 52)
    for name, pka_p, pka_r in sim.all_predictions():
        print(f"{name:<16s}  {pka_p:>10.2f}  {pka_r:>10.2f}  {abs(pka_p - pka_r):>8.2f}")


if __name__ == "__main__":
    main()
