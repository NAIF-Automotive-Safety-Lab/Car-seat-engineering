"""Neutrino oscillation calculator from B3 substrate predictions.

Substrate-derived inputs (see mixing_matrices.py):
  - PMNS angles:
        sin² θ_12 = 42 α
        sin² θ_13 = 3 α
        sin² θ_23 = 1/2 + 2π α
  - δ_CP ≈ 3π/4 (substrate suggestion)
  - Lightest neutrino mass: m_1 = 2.26 meV (NH, B3 prediction)
  - Σm_ν = 60.5 meV (B3, DESI-consistent)

Mass splittings (PDG / NuFIT 5.3 inputs, since B3 supplies the absolute
scale via m_1 + Σm_ν but not the splittings themselves):
        Δm²_21 = 7.42e-5 eV²    (solar)
        |Δm²_31| = 2.515e-3 eV² (atmospheric, NH)

Probabilities use the standard vacuum oscillation formula:

    P(ν_α → ν_β) = δ_αβ
        - 4 Σ_{i<j} Re(U*_αi U_βi U_αj U*_βj) sin²(Δm²_ij L / 4E)
        + 2 Σ_{i<j} Im(U*_αi U_βi U_αj U*_βj) sin(Δm²_ij L / 2E)

with Δm²_ij L / (4E) [dimensionless] = 1.267 · Δm²_ij[eV²] · L[km] / E[GeV].

For solar neutrinos, the matter (MSW) effect dominates and the high-E
survival probability is given by the LMA-MSW limit:
    P(ν_e → ν_e)  →  sin²θ_12 + cos²θ_12 · sin²(2θ_12)/2 · ⟨...⟩
which in the adiabatic limit simplifies to  P_ee ≈ sin²θ_12 · cos²θ_13
≈ 0.30 (at high E).  At low E, vacuum averaging gives ≈ 0.55.
We expose both regimes.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .mixing_matrices import MixingMatrices


# Mass-splitting inputs (PDG / NuFIT 5.3; B3 inherits these for now).
DELTA_M2_21 = 7.42e-5      # eV², solar
DELTA_M2_31_NH = 2.515e-3  # eV², atmospheric, normal hierarchy
DELTA_M2_32_IH = -2.498e-3  # eV², inverted hierarchy

# Substrate predictions.
M1_SUBSTRATE_EV = 2.26e-3      # m_1 in eV (NH)
SIGMA_M_NU_EV = 60.5e-3        # Σm_ν in eV

# Conversion factor: Δm²[eV²]·L[km]/(4·E[GeV]) → dimensionless argument
# of sin².  Standard textbook constant.
KIN_FACTOR = 1.26693281       # = (1e9·1e3) / (4·hbar·c in eV·m), GeV·km units


@dataclass
class OscillationConfig:
    """Hierarchy and CP-phase configuration."""
    hierarchy: str = "NH"      # "NH" or "IH"
    delta_CP: float | None = None  # rad; None → substrate default 3π/4


class NeutrinoOscillation:
    """Three-flavour neutrino oscillation, vacuum + simple MSW.

    Parameters are pulled from the B3 substrate via MixingMatrices,
    with mass-splitting inputs from global fits.
    """

    flavours = ("e", "mu", "tau")

    def __init__(self, config: OscillationConfig | None = None):
        self.config = config or OscillationConfig()
        self.mixing = MixingMatrices(delta_CP_pmns=self.config.delta_CP)
        self._U = self.mixing.pmns_matrix()
        # Mass eigenvalues set lazily.
        self._masses_NH: np.ndarray | None = None
        self._masses_IH: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Mass eigenvalues
    # ------------------------------------------------------------------
    def mass_eigenvalues_NH(self) -> np.ndarray:
        """Return (m1, m2, m3) [eV] in normal hierarchy.

        Anchored by B3 prediction m_1 = 2.26 meV; m_2, m_3 follow from
        Δm²_21 and |Δm²_31|.
        """
        if self._masses_NH is None:
            m1 = M1_SUBSTRATE_EV
            m2 = np.sqrt(m1 ** 2 + DELTA_M2_21)
            m3 = np.sqrt(m1 ** 2 + DELTA_M2_31_NH)
            self._masses_NH = np.array([m1, m2, m3])
        return self._masses_NH

    def mass_eigenvalues_IH(self, m3_eV: float = 0.0) -> np.ndarray:
        """Return (m1, m2, m3) [eV] in inverted hierarchy.

        IH has m_3 < m_1 < m_2 with |Δm²_32| ≈ 2.5e-3 eV² as the large gap.
        Substrate prediction Σm_ν = 60.5 meV is *incompatible* with IH
        (which has a hard floor Σ ≳ 100 meV), so we report the minimum-Σ
        IH spectrum (m_3 = 0 by default) and let callers vary m_3.
        """
        d32 = abs(DELTA_M2_32_IH)
        d21 = DELTA_M2_21
        m3 = float(m3_eV)
        m1 = np.sqrt(m3 ** 2 + d32 - d21)
        m2 = np.sqrt(m3 ** 2 + d32)
        if m3_eV == 0.0 and self._masses_IH is None:
            self._masses_IH = np.array([m1, m2, m3])
            return self._masses_IH
        return np.array([m1, m2, m3])

    def sigma_m_nu(self, hierarchy: str | None = None) -> float:
        h = hierarchy or self.config.hierarchy
        m = self.mass_eigenvalues_NH() if h == "NH" else self.mass_eigenvalues_IH()
        return float(np.sum(m))

    # ------------------------------------------------------------------
    # PMNS matrix
    # ------------------------------------------------------------------
    def pmns_matrix(self) -> np.ndarray:
        """Return the PMNS matrix from substrate angles + δ_CP."""
        return self._U

    # ------------------------------------------------------------------
    # Oscillation probability (vacuum)
    # ------------------------------------------------------------------
    def _flavour_index(self, name: str) -> int:
        return self.flavours.index(name)

    def oscillation_probability(self, alpha: str, beta: str,
                                E_GeV: float, L_km: float,
                                hierarchy: str | None = None) -> float:
        """Vacuum oscillation probability P(ν_α → ν_β)."""
        a = self._flavour_index(alpha)
        b = self._flavour_index(beta)
        U = self._U
        h = hierarchy or self.config.hierarchy
        m = self.mass_eigenvalues_NH() if h == "NH" else self.mass_eigenvalues_IH()
        # Δm²_ij = m_i² - m_j²
        dm2 = np.array([
            [m[i] ** 2 - m[j] ** 2 for j in range(3)] for i in range(3)
        ])

        P = 1.0 if a == b else 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                prod = (np.conj(U[a, i]) * U[b, i]
                        * U[a, j] * np.conj(U[b, j]))
                phase = KIN_FACTOR * dm2[j, i] * L_km / E_GeV  # j>i so positive
                P -= 4.0 * np.real(prod) * np.sin(phase) ** 2
                P += 2.0 * np.imag(prod) * np.sin(2.0 * phase)
        # Numerical clip
        return float(max(0.0, min(1.0, P)))

    # ------------------------------------------------------------------
    # Specific physics regimes
    # ------------------------------------------------------------------
    def solar_neutrino_problem(self, E_MeV: float = 10.0) -> dict:
        """Adiabatic LMA-MSW survival probability for solar ν_e.

        High-E (>5 MeV, dominated by ⁸B): matter-dominated,
            P_ee ≈ sin²θ_12 · cos⁴θ_13 ≈ 0.30
        Low-E (<1 MeV, pp neutrinos): vacuum-averaged,
            P_ee ≈ 1 − ½·sin²(2θ_12)·cos⁴θ_13 + sin⁴θ_13 ≈ 0.55
        We interpolate smoothly and report both limits + value at E_MeV.
        """
        s12 = np.sin(self.mixing.theta12_pmns)
        c12 = np.cos(self.mixing.theta12_pmns)
        s13 = np.sin(self.mixing.theta13_pmns)
        c13 = np.cos(self.mixing.theta13_pmns)
        c13_4 = c13 ** 4
        P_high = (s12 ** 2) * c13_4 + s13 ** 4
        P_low = (1.0 - 0.5 * (2.0 * s12 * c12) ** 2 * c13_4
                 - 2.0 * s12 ** 2 * c12 ** 2 * c13_4 + (s12 ** 2) * c13_4 * 0.0
                 + s13 ** 4)
        # Cleaner low-E vacuum-averaged form:
        P_low = (c13 ** 4) * (1.0 - 0.5 * (2 * s12 * c12) ** 2) + s13 ** 4
        # Smooth interpolation, transition near E ≈ 2 MeV.
        x = (E_MeV - 2.0) / 1.5
        w = 0.5 * (1.0 + np.tanh(x))
        P_at = (1.0 - w) * P_low + w * P_high
        return {
            "E_MeV": E_MeV,
            "P_ee_low_E": float(P_low),     # ≈ 0.55, pp-region
            "P_ee_high_E": float(P_high),   # ≈ 0.30, ⁸B region
            "P_ee_at_E": float(P_at),
            "transition_E_MeV": 2.0,
        }

    def atmospheric_neutrino_oscillation(self,
                                         E_GeV: float = 1.0,
                                         L_km: float = 500.0) -> dict:
        """ν_μ disappearance for atmospheric neutrinos.

        Two-flavour leading-order:
            P_μμ ≈ 1 − sin²(2θ_23) cos⁴θ_13 sin²(1.267 Δm²_31 L/E)
        First minimum (oscillation peak in disappearance) at
            1.267 · Δm²_31 · L/E = π/2  →  L/E ≈ 500 km/GeV.
        """
        P_mm = self.oscillation_probability("mu", "mu", E_GeV, L_km)
        # Locate first minimum of P_μμ vs L/E.
        ratios = np.linspace(50.0, 2000.0, 4000)
        Pvals = np.array([
            self.oscillation_probability("mu", "mu", 1.0, r)
            for r in ratios
        ])
        idx = int(np.argmin(Pvals))
        return {
            "E_GeV": E_GeV,
            "L_km": L_km,
            "P_mu_mu": P_mm,
            "first_min_L_over_E_km_per_GeV": float(ratios[idx]),
            "P_mu_mu_at_min": float(Pvals[idx]),
            "sin2_2theta23": float(np.sin(2 * self.mixing.theta23_pmns) ** 2),
        }

    def T2K_results(self) -> dict:
        """Long-baseline ν_μ → ν_e and ν_μ disappearance at T2K.

        Beam: ⟨E⟩ ≈ 0.6 GeV, baseline L = 295 km.
        """
        E, L = 0.6, 295.0
        P_me = self.oscillation_probability("mu", "e", E, L)
        P_mm = self.oscillation_probability("mu", "mu", E, L)
        # CP-conjugate: flip δ → −δ ⇒ Im part flips sign; we re-build.
        cfg = OscillationConfig(hierarchy=self.config.hierarchy,
                                delta_CP=-self.mixing.delta_CP_pmns)
        anti = NeutrinoOscillation(cfg)
        P_me_bar = anti.oscillation_probability("mu", "e", E, L)
        asym = (P_me - P_me_bar) / (P_me + P_me_bar) if (P_me + P_me_bar) > 0 else 0.0
        return {
            "E_GeV": E,
            "L_km": L,
            "P_mu_e": P_me,
            "P_mu_mu": P_mm,
            "P_mu_e_antinu": P_me_bar,
            "CP_asymmetry": float(asym),
            "delta_CP_rad": self.mixing.delta_CP_pmns,
        }

    def dune_predictions(self) -> dict:
        """DUNE: 1300 km baseline, broadband beam ~1−5 GeV.

        Returns appearance/disappearance over the spectrum and the
        predicted CP asymmetry curve, parameterized by δ_CP.
        """
        L = 1300.0
        Es = np.linspace(0.5, 5.0, 25)
        P_me = np.array([self.oscillation_probability("mu", "e", E, L)
                         for E in Es])
        P_mm = np.array([self.oscillation_probability("mu", "mu", E, L)
                         for E in Es])
        cfg = OscillationConfig(hierarchy=self.config.hierarchy,
                                delta_CP=-self.mixing.delta_CP_pmns)
        anti = NeutrinoOscillation(cfg)
        P_me_bar = np.array([anti.oscillation_probability("mu", "e", E, L)
                             for E in Es])
        asym = np.where(P_me + P_me_bar > 0,
                        (P_me - P_me_bar) / (P_me + P_me_bar), 0.0)
        # Δm²_31 sign discrimination figure of merit
        cfg_ih = OscillationConfig(hierarchy="IH",
                                   delta_CP=self.mixing.delta_CP_pmns)
        ih = NeutrinoOscillation(cfg_ih)
        P_me_IH = np.array([ih.oscillation_probability("mu", "e", E, L)
                            for E in Es])
        return {
            "L_km": L,
            "E_GeV": Es.tolist(),
            "P_mu_e_NH": P_me.tolist(),
            "P_mu_mu_NH": P_mm.tolist(),
            "P_mu_e_antinu_NH": P_me_bar.tolist(),
            "CP_asymmetry": asym.tolist(),
            "P_mu_e_IH": P_me_IH.tolist(),
            "delta_CP_rad": self.mixing.delta_CP_pmns,
            "max_CP_asymmetry": float(np.max(np.abs(asym))),
        }

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    def report(self) -> dict:
        """Aggregate numerical summary."""
        m_NH = self.mass_eigenvalues_NH()
        m_IH = self.mass_eigenvalues_IH()
        atm = self.atmospheric_neutrino_oscillation()
        sol = self.solar_neutrino_problem(E_MeV=10.0)
        sol_low = self.solar_neutrino_problem(E_MeV=0.3)
        t2k = self.T2K_results()
        U = self.pmns_matrix()
        unitarity = self.mixing.unitarity_residual(U)

        out = {
            "hierarchy": self.config.hierarchy,
            "delta_CP_rad": self.mixing.delta_CP_pmns,
            "mass_eigenvalues_NH_eV": m_NH.tolist(),
            "mass_eigenvalues_IH_eV": m_IH.tolist(),
            "Sigma_m_nu_NH_eV": float(np.sum(m_NH)),
            "Sigma_m_nu_IH_eV": float(np.sum(m_IH)),
            "PMNS_unitarity_residual": unitarity,
            "atm_first_min_L_over_E": atm["first_min_L_over_E_km_per_GeV"],
            "atm_P_mu_mu_at_min": atm["P_mu_mu_at_min"],
            "solar_P_ee_high_E": sol["P_ee_high_E"],
            "solar_P_ee_low_E": sol_low["P_ee_low_E"],
            "T2K_P_mu_e": t2k["P_mu_e"],
            "T2K_P_mu_mu": t2k["P_mu_mu"],
            "T2K_CP_asymmetry": t2k["CP_asymmetry"],
        }

        print("=" * 60)
        print("B3 Neutrino Oscillation Report")
        print("=" * 60)
        print(f"Hierarchy:                {out['hierarchy']}")
        print(f"δ_CP [rad]:               {out['delta_CP_rad']:.4f}")
        print(f"m1, m2, m3 (NH) [meV]:    "
              f"{1e3*m_NH[0]:.3f}, {1e3*m_NH[1]:.3f}, {1e3*m_NH[2]:.3f}")
        print(f"m1, m2, m3 (IH) [meV]:    "
              f"{1e3*m_IH[0]:.3f}, {1e3*m_IH[1]:.3f}, {1e3*m_IH[2]:.3f}")
        print(f"Σm_ν (NH) [meV]:          {1e3*out['Sigma_m_nu_NH_eV']:.2f}")
        print(f"Σm_ν (IH) [meV]:          {1e3*out['Sigma_m_nu_IH_eV']:.2f}")
        print(f"PMNS unitarity residual:  {unitarity:.2e}")
        print(f"Atm. first min at L/E:    "
              f"{out['atm_first_min_L_over_E']:.1f} km/GeV  "
              f"(P_μμ_min = {out['atm_P_mu_mu_at_min']:.3f})")
        print(f"Solar P_ee (high E,⁸B):   {out['solar_P_ee_high_E']:.3f}")
        print(f"Solar P_ee (low E, pp):   {out['solar_P_ee_low_E']:.3f}")
        print(f"T2K P(ν_μ→ν_e):           {out['T2K_P_mu_e']:.4f}")
        print(f"T2K P(ν_μ→ν_μ):           {out['T2K_P_mu_mu']:.4f}")
        print(f"T2K CP asymmetry:         {out['T2K_CP_asymmetry']:+.3f}")
        return out


if __name__ == "__main__":
    osc = NeutrinoOscillation()
    osc.report()
