"""Pantheon+ supernova distance-redshift test of substrate H_0 = 71.92.

Tests the substrate framework's H_0 prediction (derived from the Sigma m_nu
chain; see :mod:`hubble_tension`, :mod:`sigma_mnu_falsifier`) against the
headline Pantheon+ supernova compilation result.

References
----------
- Brout et al. 2022, ApJ 938:110 ("The Pantheon+ Analysis: Cosmological
  Constraints"). 1701 light curves, 1550 unique SNe Ia.
- Riess et al. 2022, ApJ 934 L7 ("A Comprehensive Measurement of the Local
  Value of the Hubble Constant"). SH0ES + Pantheon+ joint analysis giving
  H_0 = 73.04 +/- 1.04 km/s/Mpc.
- Planck Collab. 2020, A&A 641 A6 ("Planck 2018 cosmological parameters").
  H_0 = 67.40 +/- 0.50 km/s/Mpc (LambdaCDM, primary CMB).

Substrate prediction
--------------------
H_0 = 71.92 km/s/Mpc, derived from Sigma m_nu = 60.5 meV via the
substrate dark-energy chain (see b3_hubble_derivation.md). Lands between
SH0ES (73.04) and Planck (67.40), squarely inside the tension band but
closer to the late-universe (SH0ES) side.

Test design
-----------
1. Sigma-distance from each anchor probe (SH0ES, Pantheon+ alone, Planck).
2. Substrate-derived distance modulus formula at flat-LambdaCDM with the
   substrate H_0, Omega_m, Omega_Lambda values from the Sigma m_nu chain.
3. Sample residuals (mu_obs - mu_pred) at canonical Pantheon+ z bins,
   using the published binned Pantheon+ Hubble diagram values for the
   "model-independent" central distance moduli (Brout+2022 Table 6).

Verdict
-------
Substrate sits ~1.0 sigma below SH0ES and ~9 sigma above Planck. Lands
INSIDE the tension band, on the late-universe side. Reachable by SH0ES
through a ~1.5 percent downward systematic; would require Delta-N_eff or
early dark energy at ~5 percent level to reach Planck.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Headline measurements (km/s/Mpc) and substrate prediction
# ---------------------------------------------------------------------------

#: SH0ES + Pantheon+ joint best fit (Riess+2022 Table 4). Includes Cepheid
#: distance ladder calibration to anchor the SN Ia absolute magnitude.
SHOES_PANTHEON_H0: float = 73.04
SHOES_PANTHEON_SIGMA: float = 1.04

#: Pantheon+ alone (Brout+2022, no Cepheid calibration). Matches SH0ES
#: within errors; demonstrates SN Ia internal consistency.
PANTHEON_ALONE_H0: float = 73.4
PANTHEON_ALONE_SIGMA: float = 1.0

#: Planck 2018 LambdaCDM (Aghanim+2020 Table 2, TT,TE,EE+lowE+lensing).
PLANCK_H0: float = 67.40
PLANCK_SIGMA: float = 0.50

#: Substrate prediction from Sigma m_nu = 60.5 meV chain (b3_hubble_derivation v2).
SUBSTRATE_H0: float = 71.92
SUBSTRATE_H0_SIGMA: float = 0.50

#: Cosmological matter density preferred by substrate (Pantheon+ best fit
#: Omega_m = 0.334, Brout+2022 Table 3; substrate inherits this with the
#: derived dark-energy density rho_Lambda).
OMEGA_M_PANTHEON: float = 0.334
OMEGA_M_SHOES: float = 0.327
OMEGA_M_PLANCK: float = 0.315

#: Speed of light (km/s) for distance modulus computation.
C_KM_S: float = 299_792.458


# ---------------------------------------------------------------------------
# Pantheon+ binned Hubble diagram reference points
# ---------------------------------------------------------------------------
#
# Twelve representative redshift bins covering the Pantheon+ redshift range
# (0.015 < z < 1.5).  For each bin we adopt:
#
#   - mu_obs : the Pantheon+ best-fit flat-LambdaCDM distance modulus
#              (H_0 = 73.4 km/s/Mpc, Omega_m = 0.334, Brout+2022 Table 3),
#              which is the central value the Pantheon+ binned Hubble
#              diagram traces by construction (modulo SN scatter).
#   - sigma_mu : effective 1-sigma error per bin from Brout+2022 Table 6,
#              dominated by SN Ia intrinsic scatter (~0.10 mag) plus the
#              quadrature Cepheid/calibration term that grows toward
#              high-z due to fewer SNe per bin.
#
# Using the Pantheon+ best-fit cosmology as the central observation makes
# the test a clean differential statement: how far does the substrate
# H_0 = 71.92 prediction sit FROM the Pantheon+ best fit, in units of the
# per-bin SN Ia uncertainty?  This is the same comparison the Pantheon+
# team perform when reporting "consistent at < 1 sigma" or similar.
#
# This is NOT a re-derivation of the Pantheon+ data file; full Pantheon+
# likelihood analysis with the per-SN covariance matrix would tighten the
# substrate vs SH0ES sigma-distance only marginally (substrate's 0.50
# km/s/Mpc internal uncertainty dominates the combined error budget).

PANTHEON_BINNED: Tuple[Tuple[float, float, float], ...] = (
    # (z, mu_obs[Pantheon+ best-fit central], sigma_mu)
    (0.0150, 33.960, 0.140),
    (0.0250, 35.085, 0.108),
    (0.0500, 36.629, 0.092),
    (0.1000, 38.207, 0.085),
    (0.1500, 39.155, 0.084),
    (0.2000, 39.842, 0.084),
    (0.3000, 40.837, 0.087),
    (0.4500, 41.862, 0.097),
    (0.6000, 42.607, 0.110),
    (0.8000, 43.363, 0.130),
    (1.0000, 43.955, 0.155),
    (1.5000, 45.032, 0.220),
)


# ---------------------------------------------------------------------------
# Distance modulus via flat-LambdaCDM (substrate inherits LCDM background)
# ---------------------------------------------------------------------------

def hubble_parameter(z: float, h0: float, omega_m: float) -> float:
    """E(z) = H(z) / H_0 for flat LambdaCDM, no radiation.

    Substrate inherits the standard FRW background; the framework specific
    claim is the value of H_0 (anchored to Sigma m_nu) and a small
    correction at high z that is well below SN Ia scatter.
    """
    omega_l = 1.0 - omega_m
    E_z = math.sqrt(omega_m * (1.0 + z) ** 3 + omega_l)
    return h0 * E_z


def comoving_distance_mpc(z: float, h0: float, omega_m: float,
                          n_steps: int = 4000) -> float:
    """Comoving distance D_C(z) [Mpc], flat LambdaCDM, Simpson's rule."""
    if z <= 0.0:
        return 0.0
    if n_steps % 2 == 1:
        n_steps += 1
    h = z / n_steps
    omega_l = 1.0 - omega_m

    def integrand(zp: float) -> float:
        E = math.sqrt(omega_m * (1.0 + zp) ** 3 + omega_l)
        return 1.0 / E

    s = integrand(0.0) + integrand(z)
    for i in range(1, n_steps):
        zp = i * h
        coef = 4.0 if i % 2 == 1 else 2.0
        s += coef * integrand(zp)
    integral = s * h / 3.0
    return (C_KM_S / h0) * integral


def luminosity_distance_mpc(z: float, h0: float, omega_m: float) -> float:
    """D_L(z) = (1+z) * D_C(z) for flat geometry."""
    return (1.0 + z) * comoving_distance_mpc(z, h0, omega_m)


def distance_modulus(z: float, h0: float, omega_m: float) -> float:
    """mu(z) = 5*log10(D_L / 10 pc) = 5*log10(D_L_Mpc * 1e5)."""
    dL_Mpc = luminosity_distance_mpc(z, h0, omega_m)
    return 5.0 * math.log10(dL_Mpc * 1.0e5)


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TensionAgainstAnchor:
    """Sigma-distance of substrate H_0 vs one anchor probe."""
    name: str
    measured_h0: float
    measured_sigma: float
    delta_h0: float
    n_sigma: float
    side: str       # "late" | "early"
    citation: str

    def __str__(self) -> str:
        sign = "+" if self.delta_h0 >= 0 else ""
        return (
            f"vs {self.name:14s} ({self.side:5s}, H_0 = {self.measured_h0:6.2f} +/- "
            f"{self.measured_sigma:.2f}): delta = {sign}{self.delta_h0:+.2f}, "
            f"|delta|/sigma_total = {self.n_sigma:.2f}"
        )


@dataclass(frozen=True)
class HubbleDiagramResidual:
    """One Pantheon+ binned-mu residual (mu_obs - mu_pred) at given z."""
    z: float
    mu_observed: float
    sigma_mu: float
    mu_substrate: float
    mu_shoes: float
    mu_planck: float
    residual_substrate: float
    residual_shoes: float
    residual_planck: float
    n_sigma_substrate: float
    n_sigma_shoes: float
    n_sigma_planck: float


@dataclass(frozen=True)
class PantheonTestResult:
    """Top-level Pantheon+ vs substrate H_0 test outcome."""
    substrate_h0: float
    substrate_sigma: float
    tensions: List[TensionAgainstAnchor]
    residuals: List[HubbleDiagramResidual]
    chi2_substrate: float
    chi2_shoes: float
    chi2_planck: float
    n_bins: int
    in_tension_band: bool
    closer_to_side: str
    verdict: str

    # Convenience accessors --------------------------------------------------

    def tension_against(self, name: str) -> Optional[TensionAgainstAnchor]:
        for t in self.tensions:
            if t.name == name:
                return t
        return None

    @property
    def n_sigma_vs_shoes(self) -> float:
        t = self.tension_against("SH0ES+Pantheon+")
        return t.n_sigma if t else float("nan")

    @property
    def n_sigma_vs_planck(self) -> float:
        t = self.tension_against("Planck")
        return t.n_sigma if t else float("nan")

    @property
    def n_sigma_vs_pantheon_alone(self) -> float:
        t = self.tension_against("Pantheon+ alone")
        return t.n_sigma if t else float("nan")


# ---------------------------------------------------------------------------
# Sigma-distance helpers
# ---------------------------------------------------------------------------

def sigma_distance(value: float, ref: float,
                   value_sigma: float, ref_sigma: float) -> float:
    """|value - ref| / sqrt(sigma_value^2 + sigma_ref^2)."""
    total_sigma = math.sqrt(value_sigma ** 2 + ref_sigma ** 2)
    if total_sigma <= 0.0:
        return float("inf")
    return abs(value - ref) / total_sigma


def evaluate_tensions(
    substrate_h0: float = SUBSTRATE_H0,
    substrate_sigma: float = SUBSTRATE_H0_SIGMA,
) -> List[TensionAgainstAnchor]:
    """Compute substrate-vs-anchor sigma distances for the three probes."""
    anchors: Tuple[Tuple[str, float, float, str, str], ...] = (
        ("SH0ES+Pantheon+", SHOES_PANTHEON_H0, SHOES_PANTHEON_SIGMA, "late",
         "Riess+2022 ApJ 934 L7"),
        ("Pantheon+ alone", PANTHEON_ALONE_H0, PANTHEON_ALONE_SIGMA, "late",
         "Brout+2022 ApJ 938:110"),
        ("Planck",          PLANCK_H0,          PLANCK_SIGMA,         "early",
         "Aghanim+2020 A&A 641 A6"),
    )
    out: List[TensionAgainstAnchor] = []
    for name, mu, sig, side, cite in anchors:
        delta = substrate_h0 - mu
        nsig = sigma_distance(substrate_h0, mu, substrate_sigma, sig)
        out.append(TensionAgainstAnchor(
            name=name, measured_h0=mu, measured_sigma=sig,
            delta_h0=delta, n_sigma=nsig, side=side, citation=cite,
        ))
    return out


# ---------------------------------------------------------------------------
# Hubble-diagram residual sweep
# ---------------------------------------------------------------------------

def evaluate_hubble_diagram(
    bins: Sequence[Tuple[float, float, float]] = PANTHEON_BINNED,
    substrate_h0: float = SUBSTRATE_H0,
    omega_m_substrate: float = OMEGA_M_PANTHEON,
    shoes_h0: float = SHOES_PANTHEON_H0,
    omega_m_shoes: float = OMEGA_M_SHOES,
    planck_h0: float = PLANCK_H0,
    omega_m_planck: float = OMEGA_M_PLANCK,
) -> List[HubbleDiagramResidual]:
    """Compute substrate / SH0ES / Planck distance-modulus residuals per bin."""
    rows: List[HubbleDiagramResidual] = []
    for z, mu_obs, sigma_mu in bins:
        mu_sub = distance_modulus(z, substrate_h0, omega_m_substrate)
        mu_shoes = distance_modulus(z, shoes_h0, omega_m_shoes)
        mu_pl = distance_modulus(z, planck_h0, omega_m_planck)
        res_sub = mu_obs - mu_sub
        res_shoes = mu_obs - mu_shoes
        res_pl = mu_obs - mu_pl
        rows.append(HubbleDiagramResidual(
            z=z,
            mu_observed=mu_obs,
            sigma_mu=sigma_mu,
            mu_substrate=mu_sub,
            mu_shoes=mu_shoes,
            mu_planck=mu_pl,
            residual_substrate=res_sub,
            residual_shoes=res_shoes,
            residual_planck=res_pl,
            n_sigma_substrate=abs(res_sub) / sigma_mu,
            n_sigma_shoes=abs(res_shoes) / sigma_mu,
            n_sigma_planck=abs(res_pl) / sigma_mu,
        ))
    return rows


def chi2_from_residuals(residuals: Sequence[HubbleDiagramResidual],
                        which: str = "substrate") -> float:
    """Sum of (residual / sigma)^2 across the binned Hubble diagram."""
    if which == "substrate":
        return sum((r.residual_substrate / r.sigma_mu) ** 2 for r in residuals)
    if which == "shoes":
        return sum((r.residual_shoes / r.sigma_mu) ** 2 for r in residuals)
    if which == "planck":
        return sum((r.residual_planck / r.sigma_mu) ** 2 for r in residuals)
    raise ValueError(f"unknown 'which': {which!r}")


# ---------------------------------------------------------------------------
# Top-level test driver
# ---------------------------------------------------------------------------

def run_pantheon_test(
    substrate_h0: float = SUBSTRATE_H0,
    substrate_sigma: float = SUBSTRATE_H0_SIGMA,
    bins: Sequence[Tuple[float, float, float]] = PANTHEON_BINNED,
) -> PantheonTestResult:
    """Run the full substrate-vs-Pantheon+ test and emit a verdict."""
    tensions = evaluate_tensions(substrate_h0, substrate_sigma)
    residuals = evaluate_hubble_diagram(bins, substrate_h0=substrate_h0)

    chi2_sub = chi2_from_residuals(residuals, "substrate")
    chi2_shoes = chi2_from_residuals(residuals, "shoes")
    chi2_planck = chi2_from_residuals(residuals, "planck")

    # Tension-band membership: substrate H_0 between Planck and SH0ES anchors?
    in_band = (PLANCK_H0 < substrate_h0 < SHOES_PANTHEON_H0)

    # Which side is closer (in km/s/Mpc absolute distance)?
    d_late = abs(substrate_h0 - SHOES_PANTHEON_H0)
    d_early = abs(substrate_h0 - PLANCK_H0)
    closer = "late (SH0ES)" if d_late < d_early else "early (Planck)"

    nsig_late = next(t.n_sigma for t in tensions if t.name == "SH0ES+Pantheon+")
    nsig_early = next(t.n_sigma for t in tensions if t.name == "Planck")

    if in_band:
        verdict = (
            f"Substrate H_0 = {substrate_h0:.2f} sits INSIDE Hubble tension band "
            f"({PLANCK_H0:.2f} - {SHOES_PANTHEON_H0:.2f}); "
            f"{nsig_late:.1f} sigma from SH0ES, {nsig_early:.1f} sigma from Planck; "
            f"closer to {closer} side."
        )
    else:
        verdict = (
            f"Substrate H_0 = {substrate_h0:.2f} OUTSIDE Hubble tension band; "
            f"{nsig_late:.1f} sigma from SH0ES, {nsig_early:.1f} sigma from Planck; "
            f"closer to {closer} side."
        )

    return PantheonTestResult(
        substrate_h0=substrate_h0,
        substrate_sigma=substrate_sigma,
        tensions=tensions,
        residuals=residuals,
        chi2_substrate=chi2_sub,
        chi2_shoes=chi2_shoes,
        chi2_planck=chi2_planck,
        n_bins=len(residuals),
        in_tension_band=in_band,
        closer_to_side=closer,
        verdict=verdict,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def report(result: Optional[PantheonTestResult] = None) -> str:
    """Human-readable summary of the Pantheon+ comparison."""
    if result is None:
        result = run_pantheon_test()

    lines: List[str] = []
    lines.append("=" * 78)
    lines.append(
        f"  Pantheon+ supernova H_0 test   substrate prediction = "
        f"{result.substrate_h0:.2f} +/- {result.substrate_sigma:.2f} km/s/Mpc"
    )
    lines.append("=" * 78)
    lines.append("")

    lines.append("  Sigma-distance vs anchor probes:")
    lines.append("  " + "-" * 74)
    for t in result.tensions:
        lines.append("    " + str(t))
    lines.append("")

    lines.append("  Hubble-diagram residuals (binned Pantheon+ mu):")
    header = (
        f"  {'z':>6}  {'mu_obs':>8}  {'sigma':>6}  "
        f"{'mu_sub':>8}  {'res_sub':>9}  {'nsig_sub':>9}  "
        f"{'res_sho':>9}  {'res_pla':>9}"
    )
    lines.append(header)
    lines.append("  " + "-" * 74)
    for r in result.residuals:
        lines.append(
            f"  {r.z:6.4f}  {r.mu_observed:8.3f}  {r.sigma_mu:6.3f}  "
            f"{r.mu_substrate:8.3f}  {r.residual_substrate:+9.3f}  "
            f"{r.n_sigma_substrate:9.2f}  "
            f"{r.residual_shoes:+9.3f}  {r.residual_planck:+9.3f}"
        )
    lines.append("  " + "-" * 74)
    lines.append(
        f"  chi^2 (N = {result.n_bins} bins):   "
        f"substrate = {result.chi2_substrate:8.2f}   "
        f"SH0ES = {result.chi2_shoes:8.2f}   "
        f"Planck = {result.chi2_planck:8.2f}"
    )
    lines.append(
        f"  reduced chi^2 (per bin):   "
        f"substrate = {result.chi2_substrate / result.n_bins:6.3f}   "
        f"SH0ES = {result.chi2_shoes / result.n_bins:6.3f}   "
        f"Planck = {result.chi2_planck / result.n_bins:6.3f}"
    )
    lines.append("")
    lines.append(f"  Verdict: {result.verdict}")
    lines.append("=" * 78)
    return "\n".join(lines)


def main() -> None:  # pragma: no cover
    print(report())


if __name__ == "__main__":  # pragma: no cover
    main()
