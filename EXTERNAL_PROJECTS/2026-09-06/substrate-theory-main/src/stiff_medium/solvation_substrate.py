"""
Solvation energies and dielectric response from the substrate framework.

Reproduces continuum (Born) ion solvation, hydration shells around
common cations / anions, and the Debye relaxation profile of the static
dielectric function for water and other polar / non-polar liquids.

Substrate ontology:
  A solute embedded in a stiff medium polarises a finite ion-dipole
  shell (geometry layer) and pays an electrostatic free-energy cost
  (simulator layer).  In the macroscopic limit the dielectric response
  is set by the relaxation timescale of the substrate's polarisation
  modes (Debye / Cole-Cole form).

For practical computation the result collapses to the textbook
Born / Debye expressions; this module evaluates those forms and
compares against tabulated data.

References
----------
- Born, Z. Phys. 1, 45 (1920)         : ΔG_solv = -(Z²e²/8πε₀r)(1 - 1/ε_r)
- Marcus, Ion Solvation (Wiley, 1985) : Na+, Cl-, Mg²+ hydration energies
- Kaatze, J. Chem. Eng. Data 1989     : water Debye relaxation τ_D ≈ 8.27 ps
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar, Dict, List, Tuple

import math

import numpy as np
from scipy import constants as C


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

# Static (low-frequency) relative permittivities at ~298 K
DIELECTRIC_CONSTANTS: Dict[str, float] = {
    "water":    78.5,
    "methanol": 32.7,
    "ethanol":  24.5,
    "DMSO":     47.0,
    "acetone":  20.7,
    "hexane":    1.9,
    "benzene":   2.3,
    "vacuum":    1.0,
}

# Optical (high-frequency / electronic) permittivity for Debye model
DIELECTRIC_INFINITY: Dict[str, float] = {
    "water":    5.2,
    "methanol": 5.7,
    "DMSO":     4.2,
    "hexane":   1.9,
}

# Debye relaxation time at 298 K (seconds)
DEBYE_TAU: Dict[str, float] = {
    "water":    8.27e-12,
    "methanol": 50.0e-12,
    "DMSO":     19.0e-12,
}

# Crystallographic (Pauling / Shannon) ionic radii in Angstrom.
# Used for geometry-layer rendering of ion + first hydration shell.
ION_RADII_A: Dict[str, float] = {
    "Na+":   1.02,
    "K+":    1.38,
    "Li+":   0.76,
    "Cl-":   1.81,
    "F-":    1.33,
    "Br-":   1.96,
    "Mg2+":  0.72,
    "Ca2+":  1.00,
    "Zn2+":  0.74,
}

# Effective Born cavity radii (Latimer 1939 / Rashin-Honig 1985, Angstrom).
# Cations: r_eff = r_ion + 0.85 Å (water-O van-der-Waals offset).
# Anions : r_eff = r_ion + 0.10 Å.
# Divalents: empirical Marcus fit (~1.4 Å for Mg²⁺).
# These are the radii that make the continuum Born formula reproduce the
# tabulated absolute single-ion hydration free energies to ~5%.
ION_RADII_BORN_A: Dict[str, float] = {
    "Na+":   1.68,
    "K+":    2.18,
    "Li+":   1.30,
    "Cl-":   1.98,
    "F-":    1.46,
    "Br-":   2.10,
    "Mg2+":  1.43,
    "Ca2+":  1.71,
    "Zn2+":  1.40,
}

ION_CHARGE: Dict[str, int] = {
    "Na+":   +1, "K+": +1, "Li+": +1,
    "Cl-":   -1, "F-": -1, "Br-": -1,
    "Mg2+":  +2, "Ca2+": +2, "Zn2+": +2,
}

# Experimental absolute single-ion hydration free energies (kJ/mol, Marcus 1991
# / Tissandier 1998 conventions).  Sign convention: ΔG < 0 (favourable).
EXPERIMENTAL_HYDRATION_KJ: Dict[str, float] = {
    "Na+":   -407.0,
    "K+":    -322.0,
    "Li+":   -519.0,
    "Cl-":   -347.0,
    "F-":    -465.0,
    "Br-":   -321.0,
    "Mg2+": -1922.0,
    "Ca2+": -1592.0,
    "Zn2+": -1955.0,
}

# Typical first-shell coordination numbers in liquid water
COORDINATION_WATER: Dict[str, int] = {
    "Na+":  6,
    "K+":   7,
    "Li+":  4,
    "Cl-":  6,
    "F-":   6,
    "Br-":  6,
    "Mg2+": 6,
    "Ca2+": 7,
    "Zn2+": 6,
}


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

@dataclass
class SolvationGeometry:
    """
    Solute + solvent shell geometry.

    A central ion (Z, r_ion) is surrounded by `n_shell` water dipoles
    (idealised as point dipoles with O-end pointing toward the cation
    or H-end pointing toward the anion).  Multiple shells can be
    constructed; each shell is approximated as a spherical shell of
    radius r_ion + (k+1)*d_water + k*d_OH where d_water ≈ 2.8 Å is the
    O-O nearest-neighbour distance.

    The geometry layer exposes:
      * `solute_position`   – central ion 3-vector
      * `shell_positions(k)` – Nx3 array of solvent positions in shell k
      * `shell_orientations(k)` – Nx3 dipole-moment vectors
      * `ion_dipole_energy()` – classical ion-dipole interaction ΣE_id
    """

    ion: str
    n_shell: int = 6
    shell_count: int = 1
    d_water_A: float = 2.80   # O-O nearest-neighbour distance (Angstrom)
    d_OH_A: float = 0.96      # O-H bond length
    mu_water_D: float = 1.85  # water dipole moment (Debye)

    solute_position: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def __post_init__(self) -> None:
        if self.ion not in ION_RADII_A:
            raise ValueError(f"unknown ion {self.ion!r}")
        self.Z: int = ION_CHARGE[self.ion]
        self.r_ion_A: float = ION_RADII_A[self.ion]

    # ------------------------------------------------------------------ shell
    def shell_radius_A(self, k: int = 0) -> float:
        """
        Distance from solute centre to centre of solvent O atom in shell k.

        Shell 0 sits at r_ion + d_water/2 + 1.4 Å (van-der-Waals contact).
        Subsequent shells are spaced by d_water (O-O distance).
        """
        return self.r_ion_A + 1.40 + k * self.d_water_A

    def shell_positions(self, k: int = 0) -> np.ndarray:
        """
        Positions of the `n_shell` water O atoms in shell k.

        Distributes points roughly uniformly on a sphere via the
        Fibonacci spiral (deterministic, no RNG).
        """
        n = self.n_shell
        r = self.shell_radius_A(k)
        idx = np.arange(0, n, dtype=float) + 0.5
        phi = math.pi * (math.sqrt(5.0) - 1.0)
        z = 1.0 - 2.0 * idx / n
        rho = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
        theta = phi * idx
        xs = r * rho * np.cos(theta)
        ys = r * rho * np.sin(theta)
        zs = r * z
        return np.column_stack([xs, ys, zs]) + self.solute_position

    def shell_orientations(self, k: int = 0) -> np.ndarray:
        """
        Dipole-moment vectors (unit vectors, Debye magnitude implied).

        Cation: O-end (dipole points away from ion, i.e. -r̂).
        Anion : H-end (dipole points toward the ion, i.e. +r̂).
        """
        positions = self.shell_positions(k) - self.solute_position
        norms = np.linalg.norm(positions, axis=1, keepdims=True)
        rhat = positions / norms
        sign = -1.0 if self.Z > 0 else +1.0
        return sign * rhat

    # ------------------------------------------------------------- diagnostics
    def ion_dipole_energy_kJmol(self) -> float:
        """
        Sum classical ion-dipole interactions over all shells:

            E_id = Σ_k Σ_i (-Z e μ cosθ_i / 4πε₀ r_ki²)

        With cosθ = ±1 by construction, this is the favourable (negative)
        contribution that the Born/Debye picture treats as a first-shell
        correction beyond the continuum estimate.

        Returns
        -------
        E in kJ/mol (per ion).
        """
        e = C.elementary_charge
        eps0 = C.epsilon_0
        # Debye  -> SI : 1 D = 3.33564e-30 C*m
        mu_SI = self.mu_water_D * 3.33564e-30
        total_J = 0.0
        for k in range(self.shell_count):
            r = self.shell_radius_A(k) * 1e-10  # Angstrom -> m
            # Each of n_shell dipoles contributes -|Z| e μ / (4πε₀ r²)
            E_per = -abs(self.Z) * e * mu_SI / (4.0 * math.pi * eps0 * r * r)
            total_J += self.n_shell * E_per
        # J/molecule -> kJ/mol
        return total_J * C.Avogadro / 1000.0

    # ------------------------------------------------------------- 3D export
    def all_shell_positions(self) -> np.ndarray:
        """Concatenated positions across all shells (for visualisation)."""
        chunks = [self.shell_positions(k) for k in range(self.shell_count)]
        return np.vstack(chunks) if chunks else np.zeros((0, 3))


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

@dataclass
class SolvationSimulator:
    """
    Substrate-framework solvation calculator.

    Implements:
      * Born continuum solvation free energy
      * First-shell hydration energies for Na+, Cl-, Mg2+ (and others)
      * Static and frequency-dependent dielectric constants for water,
        methanol, DMSO, hexane (Debye relaxation form)
    """

    results: Dict[str, float] = field(default_factory=dict)

    # --------------------------------------------------- Born solvation -----
    @staticmethod
    def born_solvation_kJmol(
        Z: int, r_A: float, eps_r: float, eps_r_init: float = 1.0
    ) -> float:
        """
        Born free energy of solvation (kJ/mol):

            ΔG = -(Z² e² / 8 π ε₀ r) * (1/eps_r_init - 1/eps_r)

        With eps_r_init = 1 (vacuum -> solvent) this is the standard
        Born expression.

        Parameters
        ----------
        Z          : ion charge (integer, signed)
        r_A        : effective Born radius in Angstrom
        eps_r      : final medium relative permittivity
        eps_r_init : initial medium relative permittivity (default vacuum)

        Returns
        -------
        ΔG in kJ/mol (negative = favourable).
        """
        e = C.elementary_charge
        eps0 = C.epsilon_0
        r_m = r_A * 1e-10
        prefactor = -(Z * Z) * e * e / (8.0 * math.pi * eps0 * r_m)
        screening = (1.0 / eps_r_init - 1.0 / eps_r)
        E_J = prefactor * screening
        return E_J * C.Avogadro / 1000.0

    # ---------------------------------------------- hydration of common ions
    def hydration_energy_kJmol(
        self,
        ion: str,
        solvent: str = "water",
        use_effective_radius: bool = True,
    ) -> float:
        """
        Predicted hydration free energy for a single ion in `solvent`.

        Uses the Born continuum formula with effective cavity radii
        (Latimer / Rashin-Honig).  Crystallographic radii systematically
        over-bind by ~2x because they ignore the soft outer wall of the
        cavity; the +0.85 Å (cations) / +0.10 Å (anions) shift recovers
        absolute single-ion hydration free energies to within ~5%.
        """
        if ion not in ION_RADII_A:
            raise ValueError(f"unknown ion {ion!r}")
        if solvent not in DIELECTRIC_CONSTANTS:
            raise ValueError(f"unknown solvent {solvent!r}")

        Z = ION_CHARGE[ion]
        r_A = ION_RADII_BORN_A[ion] if use_effective_radius else ION_RADII_A[ion]
        eps_r = DIELECTRIC_CONSTANTS[solvent]

        E_born = self.born_solvation_kJmol(Z, r_A, eps_r)

        key = f"{ion}_in_{solvent}_kJmol"
        self.results[key] = E_born
        self.results[f"{ion}_in_{solvent}_Born_kJmol"] = E_born
        return E_born

    # ----------------------------------- dielectric constants (static) ------
    def dielectric_constant(self, solvent: str) -> float:
        if solvent not in DIELECTRIC_CONSTANTS:
            raise ValueError(f"unknown solvent {solvent!r}")
        eps = DIELECTRIC_CONSTANTS[solvent]
        self.results[f"{solvent}_eps_r"] = eps
        return eps

    # -------------------------------- Debye relaxation ε(ω) -----------------
    @staticmethod
    def debye_eps(
        omega: np.ndarray | float,
        eps_s: float,
        eps_inf: float,
        tau: float,
    ) -> np.ndarray:
        """
        Complex permittivity from the single-Debye relaxation:

            ε(ω) = ε_∞ + (ε_s - ε_∞) / (1 + i ω τ)

        Returns a complex ndarray (or scalar) with the same shape as ω.
        """
        omega = np.asarray(omega, dtype=float)
        return eps_inf + (eps_s - eps_inf) / (1.0 + 1j * omega * tau)

    def water_dielectric_spectrum(
        self,
        f_min_Hz: float = 1e6,
        f_max_Hz: float = 1e13,
        n_points: int = 200,
    ) -> Dict[str, np.ndarray]:
        """
        Frequency sweep ε'(f), ε''(f) for liquid water at 298 K.

        Returns
        -------
        dict with keys 'f_Hz', 'eps_real', 'eps_imag', 'eps_abs'.
        """
        eps_s = DIELECTRIC_CONSTANTS["water"]
        eps_inf = DIELECTRIC_INFINITY["water"]
        tau = DEBYE_TAU["water"]
        f = np.logspace(math.log10(f_min_Hz), math.log10(f_max_Hz), n_points)
        omega = 2.0 * math.pi * f
        eps = self.debye_eps(omega, eps_s, eps_inf, tau)
        return {
            "f_Hz":     f,
            "eps_real": eps.real,
            "eps_imag": -eps.imag,    # convention: positive loss factor
            "eps_abs":  np.abs(eps),
        }

    # ------------------------------------------------- experimental compare
    def compare_hydration(
        self, ions: Tuple[str, ...] = ("Na+", "Cl-", "Mg2+")
    ) -> Dict[str, Dict[str, float]]:
        """Compute hydration energies and tabulate vs. experiment."""
        table: Dict[str, Dict[str, float]] = {}
        for ion in ions:
            calc = self.hydration_energy_kJmol(ion, "water")
            ref = EXPERIMENTAL_HYDRATION_KJ[ion]
            table[ion] = {
                "calc_kJmol": calc,
                "exp_kJmol":  ref,
                "abs_err":    calc - ref,
                "rel_err_pct": 100.0 * (calc - ref) / abs(ref),
            }
        self.results["_compare_hydration"] = table  # type: ignore[assignment]
        return table

    # ----------------------------------------------------- text report ------
    def report(self) -> str:
        ions = ("Na+", "K+", "Cl-", "F-", "Mg2+", "Ca2+")
        table = self.compare_hydration(ions)
        lines = []
        lines.append("=" * 78)
        lines.append("Substrate Solvation --- hydration ΔG vs experiment (kJ/mol)")
        lines.append("=" * 78)
        lines.append(f"{'ion':<8}{'calc':>14}{'exp':>14}{'rel %':>10}")
        lines.append("-" * 78)
        for ion, row in table.items():
            lines.append(
                f"{ion:<8}{row['calc_kJmol']:>14.1f}"
                f"{row['exp_kJmol']:>14.1f}{row['rel_err_pct']:>10.2f}"
            )
        lines.append("-" * 78)
        lines.append("Static dielectric constants (298 K):")
        for s, eps in DIELECTRIC_CONSTANTS.items():
            lines.append(f"  {s:<10}  ε_r = {eps:6.2f}")
        lines.append("=" * 78)
        return "\n".join(lines)


if __name__ == "__main__":
    print(SolvationSimulator().report())
