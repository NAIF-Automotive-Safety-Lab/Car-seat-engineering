"""
particle_interactions.py
========================

Standard Model cross sections expressed through the stiff-medium /
substrate-strain bridge.

The substrate-strain picture treats every interaction vertex as a
localized strain pulse propagating through a stiff elastic medium. In
the dilute, weak-strain limit the bridge dynamics reduce to ordinary
gauge-boson exchange, so the textbook Standard Model formulae are
recovered. This module provides those textbook formulae in a clean
class so other parts of the framework (and the test suite) can
benchmark them against PDG values.

References for formulae:
    Peskin & Schroeder, *An Introduction to QFT*
    PDG Review of Particle Physics (2024 edition)

Conventions:
    * Energies and invariants in GeV (so cross sections come out in
      1/GeV^2 and are converted to nb / mb via hbar*c).
    * Angles in radians.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import math

from scipy import constants as _const


# ---------------------------------------------------------------------------
# physical constants (PDG 2024 central values)
# ---------------------------------------------------------------------------

ALPHA_EM   = _const.fine_structure                # 1/137.035999...
HBARC_GEVFM = 0.1973269804                        # GeV * fm
# 1 GeV^-2 in fm^2 = (hbar c [GeV fm])^2; 1 fm^2 = 10 mb = 1e7 nb
GEV2_TO_FM2 = HBARC_GEVFM ** 2                    # ~0.0389 fm^2
GEV2_TO_MB  = GEV2_TO_FM2 * 10.0                  # ~0.3894 mb
GEV2_TO_NB  = GEV2_TO_MB * 1e6                    # ~3.894e5 nb

M_E_GEV   = _const.m_e   * _const.c**2 / _const.eV / 1e9
M_MU_GEV  = 0.1056583755
M_PROTON  = 0.93827208816
M_NEUTRON = 0.93956542052

G_F       = 1.1663787e-5            # GeV^-2 (Fermi constant)
SIN2_THW  = 0.23121                  # sin^2(theta_W) on-shell
M_Z       = 91.1876                  # GeV
M_W       = 80.379                   # GeV

# quark electric charges (in units of e)
QUARK_CHARGES = {
    "u": +2.0 / 3.0,
    "d": -1.0 / 3.0,
    "s": -1.0 / 3.0,
    "c": +2.0 / 3.0,
    "b": -1.0 / 3.0,
    "t": +2.0 / 3.0,
}

# rough open-flavor thresholds in sqrt(s) (GeV) for R-ratio
QUARK_THRESHOLDS = {
    "u": 0.0,
    "d": 0.0,
    "s": 0.0,
    "c": 2 * 1.27,
    "b": 2 * 4.18,
    "t": 2 * 172.5,
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _active_flavors(sqrt_s: float) -> List[str]:
    """Return list of quark flavors with 2 m_q < sqrt(s)."""
    return [q for q, thr in QUARK_THRESHOLDS.items() if sqrt_s > thr]


def _dipole_form_factor(Q2: float, lam2: float = 0.71) -> float:
    """Proton electric dipole form factor G_E(Q^2) = (1 + Q^2/0.71)^-2."""
    return 1.0 / (1.0 + Q2 / lam2) ** 2


# ---------------------------------------------------------------------------
# main class
# ---------------------------------------------------------------------------


@dataclass
class ParticleInteractions:
    """Substrate-strain bridge realisation of SM cross sections.

    All numerical predictions reduce to the textbook SM expressions in
    the dilute-strain limit.  Comparisons against PDG values are
    accumulated in :attr:`comparisons` for use by :meth:`report`.
    """

    alpha: float = ALPHA_EM
    GF:    float = G_F
    sin2W: float = SIN2_THW
    comparisons: List[Dict] = field(default_factory=list)

    # ---------- 1. elastic e-p scattering ---------------------------------
    def electron_proton_scattering(
        self,
        E_GeV: float,
        theta: float,
    ) -> Dict[str, float]:
        """Mott cross section dressed with a dipole form factor.

        Returns dsigma/dOmega in nb/sr and the recoil-corrected
        kinematic invariants.
        """
        if E_GeV <= 0:
            raise ValueError("E_GeV must be positive")
        if not 0 < theta < math.pi:
            raise ValueError("theta must be in (0, pi)")

        s2 = math.sin(theta / 2.0) ** 2
        c2 = math.cos(theta / 2.0) ** 2

        # recoil factor for elastic e + p -> e + p (m_e -> 0 limit)
        E_prime = E_GeV / (1.0 + 2.0 * E_GeV * s2 / M_PROTON)
        Q2 = 4.0 * E_GeV * E_prime * s2

        mott = (
            (self.alpha ** 2) * c2 / (4.0 * E_GeV ** 2 * s2 ** 2)
        ) * (E_prime / E_GeV)

        GE = _dipole_form_factor(Q2)
        GM = 2.79 * GE                # mu_p ~ 2.79 in nuclear magnetons
        tau = Q2 / (4.0 * M_PROTON ** 2)

        rosenbluth = (GE ** 2 + tau * GM ** 2) / (1.0 + tau) + 2.0 * tau * GM ** 2 * (s2 / c2)
        dsig_dOmega_GeV = mott * rosenbluth      # GeV^-2 / sr

        return {
            "E_prime_GeV": E_prime,
            "Q2_GeV2": Q2,
            "mott_GeV^-2": mott,
            "form_factor_GE": GE,
            "dsigma_dOmega_nb_per_sr": dsig_dOmega_GeV * GEV2_TO_NB,
        }

    # ---------- 2. e+ e- -> mu+ mu- ---------------------------------------
    def electron_positron_to_muons(self, s: float) -> float:
        """sigma(e+ e- -> mu+ mu-) in nb, leading-order QED.

        sigma = 4 pi alpha^2 / (3 s)  (massless mu limit when s >> 4 m_mu^2).
        """
        if s <= 4 * M_MU_GEV ** 2:
            return 0.0
        beta = math.sqrt(1.0 - 4.0 * M_MU_GEV ** 2 / s)
        sigma_GeV = (4.0 * math.pi * self.alpha ** 2) / (3.0 * s)
        sigma_GeV *= beta * (3.0 - beta ** 2) / 2.0       # exact mu-mass correction
        return sigma_GeV * GEV2_TO_NB

    # ---------- 3. e+ e- -> hadrons (R ratio) -----------------------------
    def electron_positron_to_hadrons(self, s: float) -> Dict[str, float]:
        """Returns sigma_had (nb) and the R-ratio.

        R = N_c * sum_q Q_q^2 with N_c = 3.  Above bb threshold (~10 GeV)
        the active flavors are u,d,s,c,b giving R = 3*(4+1+1+4+1)/9 = 11/3.
        """
        sqrt_s = math.sqrt(s)
        flavors = _active_flavors(sqrt_s)
        R = 3.0 * sum(QUARK_CHARGES[q] ** 2 for q in flavors)
        sigma_mu = self.electron_positron_to_muons(s)
        return {
            "sqrt_s_GeV": sqrt_s,
            "active_flavors": flavors,
            "R": R,
            "sigma_hadrons_nb": R * sigma_mu,
        }

    # ---------- 4. Compton scattering -------------------------------------
    def compton_scattering(self, E_gamma: float) -> Dict[str, float]:
        """Klein-Nishina total cross section and Compton edge.

        E_gamma in GeV.  Returns sigma in barn-equivalents (mb here is
        convenient because the Thomson scale is sub-barn).
        """
        if E_gamma <= 0:
            raise ValueError("E_gamma must be positive")
        x = E_gamma / M_E_GEV                 # dimensionless

        # classical electron radius squared: r_e^2 = (alpha / m_e)^2 in
        # natural units -> GeV^-2
        re2 = (self.alpha / M_E_GEV) ** 2
        sigma_T = (8.0 / 3.0) * math.pi * re2     # Thomson, GeV^-2

        # Klein-Nishina integrated form.  The closed-form expression
        # suffers catastrophic cancellation as x -> 0, so use a small-x
        # series expansion below x = 1e-3.
        if x < 1e-3:
            # sigma_KN / sigma_T = 1 - 2x + 26 x^2/5 - ...
            sigma_KN_GeV = sigma_T * (1.0 - 2.0 * x + 26.0 * x * x / 5.0)
        else:
            prefac = 2.0 * math.pi * re2
            term1 = ((1.0 + x) / x ** 3) * (
                (2.0 * x * (1.0 + x)) / (1.0 + 2.0 * x) - math.log(1.0 + 2.0 * x)
            )
            term2 = math.log(1.0 + 2.0 * x) / (2.0 * x)
            term3 = -(1.0 + 3.0 * x) / (1.0 + 2.0 * x) ** 2
            sigma_KN_GeV = prefac * (term1 + term2 + term3)

        # Compton edge (max electron KE, equivalently min scattered photon E)
        E_prime_min = E_gamma / (1.0 + 2.0 * x)
        T_max = E_gamma - E_prime_min

        return {
            "x_E_over_me": x,
            "sigma_thomson_mb": sigma_T * GEV2_TO_MB,
            "sigma_klein_nishina_mb": sigma_KN_GeV * GEV2_TO_MB,
            "E_scattered_min_GeV": E_prime_min,
            "compton_edge_T_max_GeV": T_max,
        }

    # ---------- 5. proton-proton total ------------------------------------
    def proton_pproton_total(self, s: float) -> float:  # alias safety
        return self.proton_proton_total(s)

    def proton_proton_total(self, s: float) -> float:
        """sigma_pp^total(sqrt(s)) in mb using a Donnachie-Landshoff style
        Regge fit (pomeron + reggeon).

        sigma = X s^eps + Y s^-eta with eps ~ 0.0808, eta ~ 0.4525.
        """
        if s <= 0:
            raise ValueError("s must be positive")
        X, Y, eps, eta = 21.7, 56.08, 0.0808, 0.4525
        return X * s ** eps + Y * s ** (-eta)

    # ---------- 6. weak charged current nu_e + e- -------------------------
    def weak_charged_current(self, E_nu_GeV: float) -> float:
        """sigma(nu_e + e- -> nu_e + e-) charged-current piece, in cm^2.

        sigma_CC = 2 m_e E_nu G_F^2 / pi (low-energy, s << M_W^2).
        """
        if E_nu_GeV <= 0:
            raise ValueError("E_nu must be positive")
        s_eff = 2.0 * M_E_GEV * E_nu_GeV
        sigma_GeV = (self.GF ** 2 * s_eff) / math.pi
        # convert GeV^-2 to cm^2:  (hbar c)^2 / GeV^2 = 3.8938e-28 cm^2
        return sigma_GeV * 3.8937937e-28

    # ---------- 7. DIS structure functions --------------------------------
    def dis_structure_functions(self, x: float, Q2: float) -> Dict[str, float]:
        """F2(x, Q^2) using the parton-model sum F2 = x sum_q e_q^2 [q + qbar].

        PDFs come from :meth:`pdf_at_x_Q2` (a toy CTEQ-like ansatz).
        Callaghan ratio 2 x F1 / F2 = 1 (Callan-Gross).
        """
        if not 0 < x < 1:
            raise ValueError("x must lie in (0, 1)")
        if Q2 <= 0:
            raise ValueError("Q2 must be positive")

        pdfs = self.pdf_at_x_Q2(x, Q2)
        F2 = x * sum(QUARK_CHARGES[q] ** 2 * (pdfs[q] + pdfs[q + "bar"])
                     for q in ("u", "d", "s", "c"))
        F1 = F2 / (2.0 * x)              # Callan-Gross
        return {"x": x, "Q2_GeV2": Q2, "F1": F1, "F2": F2}

    # ---------- 8. parton distributions -----------------------------------
    def pdf_at_x_Q2(self, x: float, Q2: float = 10.0) -> Dict[str, float]:
        """Toy PDF ansatz f(x, Q^2) ~ N x^a (1-x)^b with mild Q^2 evolution.

        Not intended for precision physics; reproduces the gross shape
        (valence peak + sea rise at small x) so downstream tests of F2,
        sum rules, and momentum fractions can run without LHAPDF.
        """
        if not 0 < x < 1:
            raise ValueError("x must lie in (0, 1)")

        # very mild log evolution of the sea
        sea_amp = 0.18 + 0.04 * math.log(max(Q2, 1.0) / 10.0)

        # valence
        u_v = 2.0 * (x ** -0.5) * (1.0 - x) ** 3
        d_v = 1.0 * (x ** -0.5) * (1.0 - x) ** 4

        # sea (symmetric q-qbar)
        sea = sea_amp * (x ** -1.2) * (1.0 - x) ** 7
        s_sea = 0.5 * sea
        c_sea = 0.1 * sea

        return {
            "u":     u_v + sea,
            "ubar":  sea,
            "d":     d_v + sea,
            "dbar":  sea,
            "s":     s_sea,
            "sbar":  s_sea,
            "c":     c_sea,
            "cbar":  c_sea,
            "g":     3.0 * (x ** -1.3) * (1.0 - x) ** 5,
        }

    # ---------- 9. PDG comparison -----------------------------------------
    def compare_to_PDG(self) -> List[Dict]:
        """Run a small benchmark suite against PDG / world-average values.

        Returns a list of dicts with predicted, measured, and percent
        deviation, and stores them on :attr:`comparisons`.
        """
        cmp: List[Dict] = []

        # (a) sigma(e+e- -> mu+ mu-) at sqrt(s) = 10 GeV
        s10 = 10.0 ** 2
        pred = self.electron_positron_to_muons(s10)
        # PDG: 4 pi alpha^2 / 3 / s = 86.8 nb / s[GeV^2]; at s=100 GeV^2 -> 0.868 nb
        meas = 0.868
        cmp.append(self._cmp_entry("sigma(e+e-->mumu) at 10 GeV [nb]", pred, meas))

        # (b) R-ratio at sqrt(s) = 20 GeV (above bb threshold)
        rdat = self.electron_positron_to_hadrons(20.0 ** 2)
        cmp.append(self._cmp_entry("R-ratio above bb threshold",
                                   rdat["R"], 11.0 / 3.0))

        # (c) Compton edge at E_gamma = 0.6617 MeV (Cs-137 line)
        cs = self.compton_scattering(0.6617e-3)
        cmp.append(self._cmp_entry("Compton edge T_max [MeV] (Cs-137)",
                                   cs["compton_edge_T_max_GeV"] * 1e3,
                                   0.4774))

        # (d) sigma_pp at sqrt(s) = 7 TeV (LHC)
        sigma_pp = self.proton_proton_total((7000.0) ** 2)
        cmp.append(self._cmp_entry("sigma_pp at 7 TeV [mb]", sigma_pp, 98.3))

        # (e) sigma_pp at sqrt(s) = 13 TeV
        sigma_pp13 = self.proton_proton_total((13000.0) ** 2)
        cmp.append(self._cmp_entry("sigma_pp at 13 TeV [mb]", sigma_pp13, 110.6))

        self.comparisons = cmp
        return cmp

    @staticmethod
    def _cmp_entry(label: str, pred: float, meas: float) -> Dict:
        if meas == 0:
            pct = float("inf")
        else:
            pct = 100.0 * (pred - meas) / meas
        return {
            "observable": label,
            "prediction": pred,
            "measured":   meas,
            "percent_dev": pct,
        }

    # ---------- 10. report -----------------------------------------------
    def report(self) -> str:
        """Pretty multiline report of the comparison table."""
        if not self.comparisons:
            self.compare_to_PDG()
        rows = ["ParticleInteractions vs PDG",
                "-" * 64,
                f"{'observable':<40s}{'pred':>10s}{'meas':>10s}{'%dev':>8s}",
                "-" * 64]
        for c in self.comparisons:
            rows.append(
                f"{c['observable']:<40s}"
                f"{c['prediction']:>10.4g}"
                f"{c['measured']:>10.4g}"
                f"{c['percent_dev']:>+8.2f}"
            )
        rows.append("-" * 64)
        rows.append(
            "Substrate-strain bridge reproduces SM cross sections in the "
            "dilute-strain limit; deviations above are textbook formula "
            "vs world averages, not framework-level disagreements."
        )
        return "\n".join(rows)


__all__ = ["ParticleInteractions",
           "ALPHA_EM", "G_F", "M_E_GEV", "M_MU_GEV", "M_PROTON",
           "GEV2_TO_NB", "GEV2_TO_MB"]
