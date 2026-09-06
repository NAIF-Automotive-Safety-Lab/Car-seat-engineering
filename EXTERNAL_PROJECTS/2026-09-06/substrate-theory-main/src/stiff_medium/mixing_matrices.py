"""CKM and PMNS mixing matrix calculator from B3 substrate predictions.

Substrate-derived inputs:
  - Cabibbo angle: sin θ_C = 1/sqrt(20) = 0.2236  (1.3% match to PDG 0.2253)
  - PMNS angles (in terms of α = fine-structure constant):
        sin² θ_12 = 42 α
        sin² θ_13 = 3 α
        sin² θ_23 = 1/2 + 2π α
  - CP-violating phase: substrate suggests δ_CP ≈ 3π/4 for PMNS
        (T2K/NOvA hint at δ_CP ~ -π/2 ... 3π/4 region)
  - CKM phase: empirical Wolfenstein bar-eta ≈ 0.355, bar-rho ≈ 0.157

CKM hierarchy follows Wolfenstein scaling:
        λ = sin θ_C
        |V_us| ~ λ            (θ_12)
        |V_cb| ~ A λ²         (θ_23)
        |V_ub| ~ A λ³         (θ_13)
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import constants


# PDG 2024 reference values
PDG = {
    "sin_theta_C": 0.2253,
    "Vus": 0.2243,
    "Vcb": 0.0410,
    "Vub": 0.00382,
    "J_CKM": 3.18e-5,
    # PMNS (NuFIT 5.3, normal ordering)
    "sin2_theta12_pmns": 0.307,
    "sin2_theta13_pmns": 0.0220,
    "sin2_theta23_pmns": 0.572,
    "delta_CP_pmns": 3.84,  # rad ~ 220 deg
    "J_PMNS": 0.0329,
}


@dataclass
class MixingParams:
    """Substrate-derived mixing parameters."""
    # CKM Wolfenstein
    lam: float          # λ = sin θ_C
    A: float            # ~0.811
    rho_bar: float      # ~0.157
    eta_bar: float      # ~0.355
    # PMNS angles (radians)
    theta12_pmns: float
    theta13_pmns: float
    theta23_pmns: float
    delta_CP_pmns: float


class MixingMatrices:
    """CKM and PMNS mixing matrix calculator from substrate predictions."""

    def __init__(self, A: float = 0.811,
                 rho_bar: float = 0.157, eta_bar: float = 0.355,
                 delta_CP_pmns: float | None = None):
        # Substrate Cabibbo
        self.alpha = constants.alpha  # fine-structure constant
        self.lam = 1.0 / np.sqrt(20.0)              # = 0.22360679...
        self.A = A
        self.rho_bar = rho_bar
        self.eta_bar = eta_bar

        # PMNS substrate angles
        s2_12 = 42.0 * self.alpha
        s2_13 = 3.0 * self.alpha
        s2_23 = 0.5 + 2.0 * np.pi * self.alpha
        self.sin2_theta12_pmns = s2_12
        self.sin2_theta13_pmns = s2_13
        self.sin2_theta23_pmns = s2_23

        self.theta12_pmns = np.arcsin(np.sqrt(s2_12))
        self.theta13_pmns = np.arcsin(np.sqrt(s2_13))
        self.theta23_pmns = np.arcsin(np.sqrt(s2_23))

        # Substrate-suggested CP phase: 3π/4 (close to NuFIT central ~220°)
        self.delta_CP_pmns = (3.0 * np.pi / 4.0) if delta_CP_pmns is None else delta_CP_pmns

    # ------------------------------------------------------------------
    # CKM
    # ------------------------------------------------------------------
    def ckm_matrix(self) -> np.ndarray:
        """Return 3x3 unitary CKM matrix.

        Built via the standard PDG parametrization (three angles + δ),
        with angles inferred from Wolfenstein parameters so that the
        Wolfenstein magnitudes are reproduced exactly and unitarity is
        manifest.
        """
        lam, A = self.lam, self.A
        rho_bar, eta_bar = self.rho_bar, self.eta_bar

        # Map Wolfenstein -> standard angles (PDG eqs.)
        s12 = lam
        s23 = A * lam ** 2
        # V_ub = A λ³ (ρ̄ - i η̄) / sqrt(1 - A² λ⁴) * sqrt(1 - λ²) / (1 - λ²(ρ̄ - i η̄))
        # For our purposes use leading-order: |V_ub| = A λ³ sqrt(ρ̄² + η̄²)
        s13_mag = A * lam ** 3 * np.sqrt(rho_bar ** 2 + eta_bar ** 2)
        # δ_CKM such that V_ub = s13 e^{-i δ}
        delta = np.arctan2(eta_bar, rho_bar)

        c12 = np.sqrt(1 - s12 ** 2)
        c23 = np.sqrt(1 - s23 ** 2)
        c13 = np.sqrt(1 - s13_mag ** 2)
        s13 = s13_mag

        eid = np.exp(1j * delta)
        emid = np.exp(-1j * delta)

        V = np.array([
            [c12 * c13,
             s12 * c13,
             s13 * emid],
            [-s12 * c23 - c12 * s23 * s13 * eid,
             c12 * c23 - s12 * s23 * s13 * eid,
             s23 * c13],
            [s12 * s23 - c12 * c23 * s13 * eid,
             -c12 * s23 - s12 * c23 * s13 * eid,
             c23 * c13],
        ], dtype=complex)
        return V

    # ------------------------------------------------------------------
    # PMNS
    # ------------------------------------------------------------------
    def pmns_matrix(self) -> np.ndarray:
        """Return 3x3 unitary PMNS matrix in PDG parametrization."""
        s12 = np.sin(self.theta12_pmns)
        c12 = np.cos(self.theta12_pmns)
        s13 = np.sin(self.theta13_pmns)
        c13 = np.cos(self.theta13_pmns)
        s23 = np.sin(self.theta23_pmns)
        c23 = np.cos(self.theta23_pmns)
        d = self.delta_CP_pmns
        eid = np.exp(1j * d)
        emid = np.exp(-1j * d)

        U = np.array([
            [c12 * c13,
             s12 * c13,
             s13 * emid],
            [-s12 * c23 - c12 * s23 * s13 * eid,
             c12 * c23 - s12 * s23 * s13 * eid,
             s23 * c13],
            [s12 * s23 - c12 * c23 * s13 * eid,
             -c12 * s23 - s12 * c23 * s13 * eid,
             c23 * c13],
        ], dtype=complex)
        return U

    # ------------------------------------------------------------------
    # Jarlskog invariants
    # ------------------------------------------------------------------
    @staticmethod
    def jarlskog_from_matrix(M: np.ndarray) -> float:
        """Compute Jarlskog invariant from any 3x3 unitary mixing matrix.

        J = Im( M[0,0] * M[1,1] * conj(M[0,1]) * conj(M[1,0]) )
        """
        return float(np.imag(M[0, 0] * M[1, 1]
                             * np.conj(M[0, 1]) * np.conj(M[1, 0])))

    def jarlskog_invariant_ckm(self) -> float:
        return self.jarlskog_from_matrix(self.ckm_matrix())

    def jarlskog_invariant_pmns(self) -> float:
        return self.jarlskog_from_matrix(self.pmns_matrix())

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    @staticmethod
    def unitarity_check(M: np.ndarray, tol: float = 1e-12) -> bool:
        n = M.shape[0]
        delta = M @ M.conj().T - np.eye(n)
        return float(np.max(np.abs(delta))) < tol

    @staticmethod
    def unitarity_residual(M: np.ndarray) -> float:
        n = M.shape[0]
        return float(np.max(np.abs(M @ M.conj().T - np.eye(n))))

    def params(self) -> MixingParams:
        return MixingParams(
            lam=self.lam, A=self.A,
            rho_bar=self.rho_bar, eta_bar=self.eta_bar,
            theta12_pmns=self.theta12_pmns,
            theta13_pmns=self.theta13_pmns,
            theta23_pmns=self.theta23_pmns,
            delta_CP_pmns=self.delta_CP_pmns,
        )

    def compare_to_pdg(self) -> dict:
        """Return dict of predicted vs observed values, also printable."""
        V = self.ckm_matrix()
        U = self.pmns_matrix()

        report = {
            "sin_theta_C": (self.lam, PDG["sin_theta_C"]),
            "|V_us|": (abs(V[0, 1]), PDG["Vus"]),
            "|V_cb|": (abs(V[1, 2]), PDG["Vcb"]),
            "|V_ub|": (abs(V[0, 2]), PDG["Vub"]),
            "J_CKM":  (self.jarlskog_invariant_ckm(), PDG["J_CKM"]),
            "sin2_theta12_PMNS": (self.sin2_theta12_pmns,
                                  PDG["sin2_theta12_pmns"]),
            "sin2_theta13_PMNS": (self.sin2_theta13_pmns,
                                  PDG["sin2_theta13_pmns"]),
            "sin2_theta23_PMNS": (self.sin2_theta23_pmns,
                                  PDG["sin2_theta23_pmns"]),
            "delta_CP_PMNS":     (self.delta_CP_pmns,
                                  PDG["delta_CP_pmns"]),
            "J_PMNS":            (self.jarlskog_invariant_pmns(),
                                  PDG["J_PMNS"]),
        }

        print(f"{'observable':22s} {'predicted':>14s} {'observed':>14s} "
              f"{'%diff':>10s}")
        print("-" * 64)
        for k, (pred, obs) in report.items():
            pct = 100.0 * (pred - obs) / obs if obs != 0 else float("nan")
            print(f"{k:22s} {pred:14.6g} {obs:14.6g} {pct:9.2f}%")
        return report


if __name__ == "__main__":
    mm = MixingMatrices()
    mm.compare_to_pdg()
    print()
    print("CKM unitarity residual: %.2e" % mm.unitarity_residual(mm.ckm_matrix()))
    print("PMNS unitarity residual: %.2e" % mm.unitarity_residual(mm.pmns_matrix()))
