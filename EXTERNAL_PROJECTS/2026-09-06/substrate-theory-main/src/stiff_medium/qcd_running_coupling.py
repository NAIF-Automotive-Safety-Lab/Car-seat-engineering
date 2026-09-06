"""QCD running coupling alpha_s(mu) with 1-loop and 2-loop RG running.

Compares the standard QCD prediction (1- and 2-loop beta function) against
the substrate framework anchor alpha_s = 16 * alpha_em at low energy
(B3 prediction).  Uses scipy.integrate.solve_ivp to integrate the RG
equation in t = ln(mu/mu_ref).

Key references
--------------
* PDG 2024: alpha_s(M_Z) = 0.1179 +/- 0.0010
* M_Z = 91.1876 GeV
* Lambda_QCD (n_f=5, MSbar) ~ 210 MeV (1-loop), ~ 209 MeV (2-loop)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp


# Constants
M_Z = 91.1876  # GeV
ALPHA_S_MZ_PDG = 0.1179
ALPHA_S_MZ_PDG_ERR = 0.0010
ALPHA_EM = 1.0 / 137.035999084  # fine-structure constant at low energy
PI = np.pi


def _beta0(nf: int) -> float:
    """1-loop beta function coefficient.  beta_0 = (33 - 2 n_f) / (12 pi)."""
    return (33.0 - 2.0 * nf) / (12.0 * PI)


def _beta1(nf: int) -> float:
    """2-loop beta function coefficient.  beta_1 = (153 - 19 n_f) / (24 pi^2)."""
    return (153.0 - 19.0 * nf) / (24.0 * PI ** 2)


def _active_flavors(mu_GeV: float) -> int:
    """Return number of active quark flavors at scale mu."""
    if mu_GeV < 1.275:
        return 3  # u,d,s
    if mu_GeV < 4.18:
        return 4  # +c
    if mu_GeV < 173.0:
        return 5  # +b
    return 6  # +t


@dataclass
class RunningResult:
    mu_GeV: float
    alpha_s: float
    loop_order: int
    nf: int


class QCDRunningCoupling:
    """Run alpha_s(mu) under 1- or 2-loop QCD beta function.

    Boundary condition is alpha_s(M_Z) = alpha_s_MZ.  The user may either
    integrate the ODE (more accurate, threshold-aware via active flavors)
    or use the closed-form Lambda_QCD expression.

    Parameters
    ----------
    alpha_s_MZ : float
        Boundary value at the Z mass.  Default = PDG 0.1179.
    """

    def __init__(self, alpha_s_MZ: float = ALPHA_S_MZ_PDG):
        self.alpha_s_MZ = float(alpha_s_MZ)
        self.M_Z = M_Z

    # ------------------------------------------------------------------
    # core RG running
    # ------------------------------------------------------------------
    def _rhs(self, t: float, y: np.ndarray, loop_order: int) -> np.ndarray:
        """RG equation d alpha_s / d ln(mu) = -2 (beta_0 alpha_s^2 + beta_1 alpha_s^3 + ...)."""
        alpha = float(y[0])
        mu = np.exp(t)
        nf = _active_flavors(mu)
        b0 = _beta0(nf)
        d = -2.0 * b0 * alpha ** 2
        if loop_order >= 2:
            b1 = _beta1(nf)
            d += -2.0 * b1 * alpha ** 3
        return np.array([d])

    def alpha_s_at_scale(self, mu_GeV: float, loop_order: int = 2) -> float:
        """Return alpha_s(mu) by integrating the RG equation from M_Z.

        Threshold matching is implicit through _active_flavors.  For the
        accuracy needed in this module we do not perform explicit MSbar
        threshold matching; the discontinuities are O(alpha_s^2) and small.
        """
        if mu_GeV <= 0:
            raise ValueError("mu_GeV must be positive")
        if loop_order not in (1, 2):
            raise ValueError("loop_order must be 1 or 2")

        t0 = np.log(self.M_Z)
        t1 = np.log(mu_GeV)
        if np.isclose(t0, t1):
            return self.alpha_s_MZ

        sol = solve_ivp(
            fun=lambda t, y: self._rhs(t, y, loop_order),
            t_span=(t0, t1),
            y0=np.array([self.alpha_s_MZ]),
            method="RK45",
            rtol=1e-9,
            atol=1e-12,
            max_step=0.05,  # small steps to handle thresholds smoothly
        )
        if not sol.success:
            raise RuntimeError(f"RG integration failed: {sol.message}")
        return float(sol.y[0, -1])

    def alpha_s_at_M_Z(self) -> float:
        """Return alpha_s at the Z pole (the anchor)."""
        return self.alpha_s_MZ

    # ------------------------------------------------------------------
    # substrate prediction
    # ------------------------------------------------------------------
    def alpha_s_substrate_anchor(self) -> float:
        """B3 substrate prediction: alpha_s(low E) = 16 * alpha_em.

        This is the framework-level anchor at the confinement scale.
        16 * alpha ~ 0.1168 -- numerically close to alpha_s(M_Z) by
        coincidence in the substrate ontology where running maps the
        low-energy substrate ratio onto the Z-scale measurement.
        """
        return 16.0 * ALPHA_EM

    # ------------------------------------------------------------------
    # Lambda_QCD extraction
    # ------------------------------------------------------------------
    def lambda_qcd_from_running(
        self, nf: int = 5, loop_order: int = 1, mu_ref_GeV: float = M_Z,
        alpha_ref: Optional[float] = None,
    ) -> float:
        """Extract Lambda_QCD (in GeV) from alpha_s at a reference scale.

        1-loop: Lambda = mu * exp( -1 / (2 beta_0 alpha_s(mu)) )
        2-loop: Lambda = mu * (b0 alpha_s)^(-b1/(2 b0^2)) * exp(-1/(2 b0 alpha_s))
        """
        a = float(alpha_ref if alpha_ref is not None else self.alpha_s_MZ)
        b0 = _beta0(nf)
        if loop_order == 1:
            return float(mu_ref_GeV * np.exp(-1.0 / (2.0 * b0 * a)))
        if loop_order == 2:
            b1 = _beta1(nf)
            pref = (b0 * a) ** (-b1 / (2.0 * b0 ** 2))
            return float(mu_ref_GeV * pref * np.exp(-1.0 / (2.0 * b0 * a)))
        raise ValueError("loop_order must be 1 or 2")

    # ------------------------------------------------------------------
    # asymptotic / confinement diagnostics
    # ------------------------------------------------------------------
    def asymptotic_freedom(self, mu_max_GeV: float = 1.0e6, n: int = 6) -> dict:
        """Verify alpha_s -> 0 as mu -> infinity.

        Returns dict with sampled (mu, alpha_s) values at logarithmically
        spaced scales up to mu_max, plus a monotonicity flag.
        """
        mus = np.logspace(np.log10(self.M_Z), np.log10(mu_max_GeV), n)
        vals = [self.alpha_s_at_scale(float(mu), loop_order=2) for mu in mus]
        monotone = all(vals[i + 1] < vals[i] for i in range(len(vals) - 1))
        return {
            "mu_GeV": mus.tolist(),
            "alpha_s": vals,
            "monotonically_decreasing": bool(monotone),
            "alpha_s_at_max": vals[-1],
        }

    def confinement_scale(self, alpha_s_target: float = 1.0,
                          loop_order: int = 1) -> float:
        """Return the scale mu (in GeV) at which alpha_s ~ alpha_s_target.

        Solved analytically from Lambda_QCD:  alpha_s(mu) = 1 / (2 beta_0 ln(mu/Lambda)).
        For target = 1, mu = Lambda * exp(1 / (2 beta_0)).
        """
        nf = 3  # low energy
        b0 = _beta0(nf)
        Lambda = self.lambda_qcd_from_running(nf=5, loop_order=1)
        # Solve 1 / (2 b0 ln(mu/Lambda)) = alpha_s_target
        return float(Lambda * np.exp(1.0 / (2.0 * b0 * alpha_s_target)))

    # ------------------------------------------------------------------
    # comparison + report
    # ------------------------------------------------------------------
    def compare_to_PDG(self) -> dict:
        """Compare model alpha_s(M_Z) and substrate anchor to PDG."""
        a_model = self.alpha_s_at_M_Z()
        a_sub = self.alpha_s_substrate_anchor()
        return {
            "PDG_alpha_s_MZ": ALPHA_S_MZ_PDG,
            "PDG_uncertainty": ALPHA_S_MZ_PDG_ERR,
            "model_alpha_s_MZ": a_model,
            "model_pct_diff": 100.0 * (a_model - ALPHA_S_MZ_PDG) / ALPHA_S_MZ_PDG,
            "substrate_alpha_s_lowE": a_sub,
            "substrate_pct_diff_vs_MZ": 100.0 * (a_sub - ALPHA_S_MZ_PDG) / ALPHA_S_MZ_PDG,
            "within_5pct_at_MZ": abs(a_model - ALPHA_S_MZ_PDG) / ALPHA_S_MZ_PDG < 0.05,
        }

    def report(self) -> str:
        cmp = self.compare_to_PDG()
        L1 = self.lambda_qcd_from_running(nf=5, loop_order=1)
        L2 = self.lambda_qcd_from_running(nf=5, loop_order=2)
        af = self.asymptotic_freedom()
        conf = self.confinement_scale(1.0)
        a_1 = self.alpha_s_at_scale(1.0, loop_order=2)
        a_10 = self.alpha_s_at_scale(10.0, loop_order=2)
        a_100 = self.alpha_s_at_scale(100.0, loop_order=2)

        lines = [
            "=" * 60,
            "QCD Running Coupling Report",
            "=" * 60,
            f"PDG alpha_s(M_Z)        = {ALPHA_S_MZ_PDG:.4f} +/- {ALPHA_S_MZ_PDG_ERR:.4f}",
            f"Model alpha_s(M_Z)      = {cmp['model_alpha_s_MZ']:.4f}"
            f"   ({cmp['model_pct_diff']:+.2f}%)",
            f"Substrate 16*alpha_em   = {cmp['substrate_alpha_s_lowE']:.4f}"
            f"   ({cmp['substrate_pct_diff_vs_MZ']:+.2f}% vs PDG MZ)",
            "-" * 60,
            f"Lambda_QCD (1-loop, nf=5) = {1000 * L1:.1f} MeV",
            f"Lambda_QCD (2-loop, nf=5) = {1000 * L2:.1f} MeV",
            "-" * 60,
            "Running alpha_s (2-loop):",
            f"  alpha_s(  1 GeV) = {a_1:.4f}",
            f"  alpha_s( 10 GeV) = {a_10:.4f}",
            f"  alpha_s(100 GeV) = {a_100:.4f}",
            "-" * 60,
            f"Asymptotic freedom:  monotone-decreasing = {af['monotonically_decreasing']}",
            f"  alpha_s({af['mu_GeV'][-1]:.1e} GeV) = {af['alpha_s_at_max']:.4f}",
            f"Confinement scale (alpha_s=1)  = {1000 * conf:.0f} MeV",
            "=" * 60,
        ]
        return "\n".join(lines)


if __name__ == "__main__":
    qcd = QCDRunningCoupling()
    print(qcd.report())
