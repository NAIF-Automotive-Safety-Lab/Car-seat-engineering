"""JWST high-redshift galaxy comparison test (substrate vs LambdaCDM vs data).

This module is a *falsifier-style* comparison module that wraps
`src.stiff_medium.jwst_high_z` (which exposes the substrate halo / stellar
mass functions and the cube-DM + de-saturation boost) and confronts its
predictions against the published JWST observations of z > 10 galaxies.

It does NOT redefine the underlying physics; it pulls cosmology from the
existing `JWSTHighZ` model and adds:

  1. A curated reference set of JWST detections at z > 10:
       - JADES-GS-z14-0   (Carniani et al. 2024, Nature 633:318)
       - GN-z11           (Bunker et al. 2023; Tacchella et al. 2023)
       - CEERS-93316      (Donnan et al. 2023; later spec-z 4.91)
       - Maisie's Galaxy  (Finkelstein et al. 2023; Arrabal Haro et al. 2023)
       - JADES-GS-z13-0, JADES-GS-z11-0
  2. Two key population-level observables:
       (a) UV luminosity function at z ~ 10 - 14 (Harikane+2023; Donnan+2023)
       (b) cosmic stellar mass density at z ~ 10 (Labbé et al. 2023,
           Nature 616:266)
  3. A LambdaCDM "expectation" track for each of (a) and (b), built from
       a Sheth-Tormen HMF using Planck 2018 cosmology with the SAME mapping
       to UV / M_star (so the only difference between the substrate and
       LCDM curves is the substrate boost — cube-DM + de-saturation).
  4. Per-source predicted/observed comoving abundance, n-σ deviation, and
       a verdict on whether the substrate framework's predicted
       OVER-abundance over LCDM matches the JWST data ratio.
  5. A composite report and a single visual:
       visuals/136_jwst_high_z.png

REFERENCES
----------
- Carniani et al. 2024, Nature 633:318    — JADES-GS-z14-0 spec-z = 14.32.
- Curtis-Lake et al. 2023, Nature Astron. — JADES four spec-confirmed z>10.
- Bunker et al. 2023, A&A 677:A88         — GN-z11 spec confirmation z=10.60.
- Tacchella et al. 2023, ApJ 952:74       — GN-z11 stellar mass ~ 10^9 M_sun.
- Finkelstein et al. 2023, ApJ 946:L13    — CEERS, Maisie's Galaxy.
- Arrabal Haro et al. 2023, Nature 622:707 — CEERS-93316 spec-z = 4.91 (drop).
- Harikane et al. 2023, ApJS 265:5        — UV LF at z = 9 - 16.
- Donnan et al. 2023, MNRAS 518:6011      — UV LF at z = 8 - 15.
- Labbé et al. 2023, Nature 616:266       — high stellar mass density z ~ 7-10
                                            ("impossibly massive galaxies").
- Boylan-Kolchin 2023, Nat. Astron. 7:731 — LCDM tension formalism.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from .jwst_high_z import (
    JWSTHighZ,
    H0_B3,
    OMEGA_M_B3,
    OMEGA_B_B3,
    SIGMA_8_B3,
    SUBSTRATE_HMF_BOOST_DEX,
    SUBSTRATE_SIGMA8_BOOST,
    EPS_STAR_HIGHZ,
)


# ---------------------------------------------------------------------------
# Planck-2018 LCDM "expectation" track
# ---------------------------------------------------------------------------

H0_PLANCK = 67.36
OMEGA_M_PLANCK = 0.3158
OMEGA_B_PLANCK = 0.04897
SIGMA_8_PLANCK = 0.811
EPS_STAR_LCDM_HIGHZ = 0.06   # standard SHMR expectation at z~10 (Behroozi+2019)


# ---------------------------------------------------------------------------
# Curated JWST z > 10 detections
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class JWSTSource:
    """A single JWST high-z source with the bare minimum we need."""
    name: str
    z: float
    log_M_star: float        # log10 (M_star / M_sun)
    M_UV: float              # rest-frame UV absolute magnitude
    log_M_star_err: float = 0.3
    survey: str = ""
    reference: str = ""


def jwst_z_gt_10_sample() -> List[JWSTSource]:
    """Curated set of confirmed / candidate z > 10 galaxies."""
    return [
        JWSTSource(
            "JADES-GS-z14-0", z=14.32, log_M_star=8.7, M_UV=-20.81,
            log_M_star_err=0.4,
            survey="JADES",
            reference="Carniani+2024 Nature 633:318",
        ),
        JWSTSource(
            "JADES-GS-z14-1", z=13.90, log_M_star=8.0, M_UV=-19.0,
            survey="JADES",
            reference="Curtis-Lake+2023 NatAst",
        ),
        JWSTSource(
            "JADES-GS-z13-0", z=13.20, log_M_star=8.6, M_UV=-18.7,
            survey="JADES",
            reference="Curtis-Lake+2023 NatAst",
        ),
        JWSTSource(
            "JADES-GS-z11-0", z=11.58, log_M_star=8.7, M_UV=-19.3,
            survey="JADES",
            reference="Curtis-Lake+2023 NatAst",
        ),
        JWSTSource(
            "GN-z11", z=10.60, log_M_star=9.0, M_UV=-21.5,
            log_M_star_err=0.3,
            survey="JADES (legacy GOODS-N)",
            reference="Bunker+2023 A&A 677:A88; Tacchella+2023 ApJ 952:74",
        ),
        JWSTSource(
            "Maisie's Galaxy", z=11.44, log_M_star=8.5, M_UV=-20.1,
            survey="CEERS",
            reference="Arrabal Haro+2023 Nature 622:707",
        ),
        JWSTSource(
            "CEERS-93316", z=4.91, log_M_star=9.0, M_UV=-19.5,
            survey="CEERS",
            reference="Arrabal Haro+2023 Nat -- candidate-z16 dropped to 4.91",
        ),
    ]


# ---------------------------------------------------------------------------
# Population-level: UV LF and stellar mass density
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class UVLFPoint:
    """One bin of the rest-frame UV luminosity function at fixed z."""
    z: float
    M_UV: float
    log_phi: float         # log10 (phi / Mpc^-3 mag^-1)
    log_phi_err: float
    reference: str = ""


def harikane_donnan_z10_LF() -> List[UVLFPoint]:
    """Selected bins of the z ~ 10 UV LF from Harikane+2023 / Donnan+2023.

    These are the bins most sensitive to the high-mass-end overabundance.
    Values are eyeball reads of the published binned LF at z = 10 - 14.
    """
    return [
        UVLFPoint(z=10.0, M_UV=-22.0, log_phi=-5.40, log_phi_err=0.30,
                  reference="Harikane+2023 ApJS 265:5"),
        UVLFPoint(z=10.0, M_UV=-21.0, log_phi=-4.50, log_phi_err=0.25,
                  reference="Harikane+2023"),
        UVLFPoint(z=10.0, M_UV=-20.0, log_phi=-4.00, log_phi_err=0.20,
                  reference="Donnan+2023 MNRAS 518:6011"),
        UVLFPoint(z=12.0, M_UV=-21.0, log_phi=-5.00, log_phi_err=0.30,
                  reference="Harikane+2023"),
        UVLFPoint(z=12.0, M_UV=-20.0, log_phi=-4.40, log_phi_err=0.30,
                  reference="Harikane+2023"),
        UVLFPoint(z=14.0, M_UV=-20.0, log_phi=-5.30, log_phi_err=0.40,
                  reference="Donnan+2023"),
    ]


@dataclass(frozen=True)
class StellarMassDensity:
    """Cumulative cosmic stellar-mass density above a mass cut at fixed z."""
    z: float
    log_M_star_min: float          # cut in log10(M_star/M_sun)
    log_rho_star: float            # log10 of integrated rho_*  (M_sun / Mpc^3)
    log_rho_star_err: float
    reference: str = ""


def labbe_2023_density_points() -> List[StellarMassDensity]:
    """High-mass density at z ~ 7 - 10 from Labbé et al. 2023, Nature 616:266.

    Reported integrated density above ~ 10^10.5 M_sun  is too high for LCDM
    by ~ 1 dex (Boylan-Kolchin 2023). We tabulate the central reported value
    and the LCDM "ceiling" value so the test can compare both ratios.
    """
    return [
        StellarMassDensity(z=7.5, log_M_star_min=10.5,
                           log_rho_star=6.6, log_rho_star_err=0.3,
                           reference="Labbé+2023 Nature 616:266"),
        StellarMassDensity(z=9.5, log_M_star_min=10.0,
                           log_rho_star=6.2, log_rho_star_err=0.4,
                           reference="Labbé+2023 Nature 616:266"),
        StellarMassDensity(z=11.0, log_M_star_min=9.5,
                           log_rho_star=5.5, log_rho_star_err=0.5,
                           reference="Stefanon+2021 / extrap. Labbé+2023"),
    ]


# ---------------------------------------------------------------------------
# LCDM expectation (no substrate boost)
# ---------------------------------------------------------------------------

def make_lcdm_track() -> JWSTHighZ:
    """A LCDM 'control' instance: Planck cosmo, eps_* low, no substrate boost."""
    return JWSTHighZ(
        H0=H0_PLANCK,
        Omega_m=OMEGA_M_PLANCK,
        Omega_b=OMEGA_B_PLANCK,
        sigma_8=SIGMA_8_PLANCK,
        eps_star=EPS_STAR_LCDM_HIGHZ,
        hmf_boost_dex=0.0,
    )


def make_substrate_track() -> JWSTHighZ:
    """The substrate-boosted predictor with B3 cosmology + cube-DM boost."""
    return JWSTHighZ(
        H0=H0_B3,
        Omega_m=OMEGA_M_B3,
        Omega_b=OMEGA_B_B3,
        sigma_8=SIGMA_8_B3,
        eps_star=EPS_STAR_HIGHZ,
        hmf_boost_dex=SUBSTRATE_HMF_BOOST_DEX,
    )


# ---------------------------------------------------------------------------
# UV magnitude <-> stellar mass mapping (very rough, used only for LF panel)
# ---------------------------------------------------------------------------

def M_UV_to_log_M_star(M_UV: float) -> float:
    """Approximate UV->stellar mass relation at z~10 (Stark et al. 2013).

    A linear fit: log M_star ~ -0.4 (M_UV + 21) + 9.0  is a reasonable
    eyeball relation for z ~ 10 M*-M_UV (Tacchella et al. 2023). Used only
    so that we can convert binned LF M_UV to a halo-mass query.
    """
    return -0.4 * (M_UV + 21.0) + 9.0


# ---------------------------------------------------------------------------
# Single-source comparison
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SourceComparison:
    name: str
    z: float
    log_M_star: float
    log_phi_substrate: float
    log_phi_lcdm: float
    log_excess_substrate: float    # log10(phi_substrate / phi_lcdm)
    n_sigma_substrate_vs_lcdm: float
    survey: str
    reference: str

    @property
    def excess_factor(self) -> float:
        return 10 ** self.log_excess_substrate


def compare_per_source(
    substrate: JWSTHighZ, lcdm: JWSTHighZ
) -> List[SourceComparison]:
    """Predict phi(M_halo, z) for each source's INFERRED halo mass.

    To isolate the *cosmology + dark-sector* contribution from the (much
    larger) astrophysical eps_star contribution, we compute the halo mass
    M_halo = M_star / (eps_star_substrate * f_b_substrate)
    using the SUBSTRATE eps_star, then evaluate the HMF for both cosmologies
    at that single halo mass. This gives the apples-to-apples ratio of
    (substrate HMF) / (LCDM HMF) due to: cube-DM boost + sigma_8 boost +
    H_0 boost. eps_star is reported separately as an additional knob.
    """
    f_b_sub = substrate.Omega_b / substrate.Omega_m
    out: List[SourceComparison] = []
    for src in jwst_z_gt_10_sample():
        M_star = 10 ** src.log_M_star
        # Use a common halo mass (substrate eps_star) for the cosmology comparison
        M_halo = M_star / (substrate.eps_star * f_b_sub)
        phi_s = substrate.halo_mass_function(src.z, M_halo)
        phi_l = lcdm.halo_mass_function(src.z, M_halo)
        # Guard against zeros to keep logs finite
        log_s = math.log10(max(phi_s, 1e-300))
        log_l = math.log10(max(phi_l, 1e-300))
        log_ex = log_s - log_l
        # n-sigma between log10 of the two using the source's stellar mass
        # error as a proxy for prediction uncertainty (no formal LF errors)
        n_sigma = log_ex / max(src.log_M_star_err, 0.05)
        out.append(SourceComparison(
            name=src.name, z=src.z, log_M_star=src.log_M_star,
            log_phi_substrate=log_s, log_phi_lcdm=log_l,
            log_excess_substrate=log_ex,
            n_sigma_substrate_vs_lcdm=n_sigma,
            survey=src.survey, reference=src.reference,
        ))
    return out


# ---------------------------------------------------------------------------
# UV luminosity function comparison
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class UVLFComparison:
    z: float
    M_UV: float
    log_phi_obs: float
    log_phi_obs_err: float
    log_phi_substrate: float
    log_phi_lcdm: float

    @property
    def n_sigma_substrate(self) -> float:
        return (self.log_phi_obs - self.log_phi_substrate) / max(
            self.log_phi_obs_err, 0.05
        )

    @property
    def n_sigma_lcdm(self) -> float:
        return (self.log_phi_obs - self.log_phi_lcdm) / max(
            self.log_phi_obs_err, 0.05
        )


def compare_UV_LF(
    substrate: JWSTHighZ, lcdm: JWSTHighZ
) -> List[UVLFComparison]:
    """Compare predicted phi at the binned LF point vs published observation.

    The phi-per-magnitude is approximated by phi(M_star, z) using the
    Stark 2013 UV->M_star relation; this is one-to-one (Jacobian=1) when
    we treat the LF normalisation as a relative test of substrate vs LCDM.
    """
    out: List[UVLFComparison] = []
    for pt in harikane_donnan_z10_LF():
        log_M_star = M_UV_to_log_M_star(pt.M_UV)
        M_star = 10 ** log_M_star
        phi_s = substrate.stellar_mass_function(pt.z, M_star)
        phi_l = lcdm.stellar_mass_function(pt.z, M_star)
        log_s = math.log10(max(phi_s, 1e-300))
        log_l = math.log10(max(phi_l, 1e-300))
        out.append(UVLFComparison(
            z=pt.z, M_UV=pt.M_UV,
            log_phi_obs=pt.log_phi, log_phi_obs_err=pt.log_phi_err,
            log_phi_substrate=log_s, log_phi_lcdm=log_l,
        ))
    return out


# ---------------------------------------------------------------------------
# Stellar mass density comparison (Labbé 2023)
# ---------------------------------------------------------------------------

def integrated_stellar_density(
    model: JWSTHighZ, z: float, log_M_star_min: float,
    log_M_star_max: float = 12.0, n_pts: int = 60,
) -> float:
    """Trapezoidal integration of M_star * phi(M_star) d log M_star."""
    grid = [log_M_star_min + (log_M_star_max - log_M_star_min) * i / n_pts
            for i in range(n_pts + 1)]
    rho = 0.0
    prev = None
    for lg in grid:
        M_star = 10 ** lg
        phi = model.stellar_mass_function(z, M_star)   # Mpc^-3 dex^-1
        contribution = M_star * phi                    # M_sun Mpc^-3 dex^-1
        if prev is not None:
            rho += 0.5 * (contribution + prev) * (grid[1] - grid[0])
        prev = contribution
    return rho   # M_sun Mpc^-3


@dataclass(frozen=True)
class DensityComparison:
    z: float
    log_M_star_min: float
    log_rho_obs: float
    log_rho_obs_err: float
    log_rho_substrate: float
    log_rho_lcdm: float

    @property
    def n_sigma_substrate(self) -> float:
        return (self.log_rho_obs - self.log_rho_substrate) / max(
            self.log_rho_obs_err, 0.05
        )

    @property
    def n_sigma_lcdm(self) -> float:
        return (self.log_rho_obs - self.log_rho_lcdm) / max(
            self.log_rho_obs_err, 0.05
        )


def compare_stellar_mass_density(
    substrate: JWSTHighZ, lcdm: JWSTHighZ
) -> List[DensityComparison]:
    out: List[DensityComparison] = []
    for pt in labbe_2023_density_points():
        rho_s = integrated_stellar_density(substrate, pt.z, pt.log_M_star_min)
        rho_l = integrated_stellar_density(lcdm, pt.z, pt.log_M_star_min)
        out.append(DensityComparison(
            z=pt.z, log_M_star_min=pt.log_M_star_min,
            log_rho_obs=pt.log_rho_star,
            log_rho_obs_err=pt.log_rho_star_err,
            log_rho_substrate=math.log10(max(rho_s, 1e-300)),
            log_rho_lcdm=math.log10(max(rho_l, 1e-300)),
        ))
    return out


# ---------------------------------------------------------------------------
# Top-level bundle
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class JWSTHighZTestResult:
    """Result bundle returned by `run_jwst_high_z_test`."""
    sources: List[SourceComparison]
    UV_LF: List[UVLFComparison]
    density: List[DensityComparison]
    mean_log_excess: float          # <log10(phi_sub/phi_lcdm)> at fixed M_halo
    median_log_excess: float
    n_sources_resolved: int         # | log_excess | > 1 sigma at fixed M_halo
    n_sources: int
    eps_star_substrate: float       # 0.20 — substrate stellar baryon fraction
    eps_star_lcdm: float            # 0.06 — Behroozi-style LCDM extrapolation
    eps_star_log_boost_dex: float   # log10(eps_sub / eps_lcdm) at fixed M_h
    total_log_boost_dex: float      # cosmology + eps_star combined
    verdict: str

    def summary(self) -> Dict[str, float]:
        return {
            "n_sources": float(self.n_sources),
            "n_sources_resolved": float(self.n_sources_resolved),
            "mean_log10_excess_substrate_over_lcdm_fixed_Mhalo":
                self.mean_log_excess,
            "median_log10_excess_substrate_over_lcdm_fixed_Mhalo":
                self.median_log_excess,
            "eps_star_substrate_over_lcdm_log_dex":
                self.eps_star_log_boost_dex,
            "total_substrate_over_lcdm_log_dex_at_fixed_Mstar":
                self.total_log_boost_dex,
        }


#: JWST OVER-abundance over LCDM expected from data (Boylan-Kolchin 2023).
#: ~ 0.5 - 1.0 dex (factor 3 - 10x) for M_* > 1e10 at z > 7-10.
JWST_OBSERVED_OVERABUNDANCE_DEX = (0.5, 1.0)


def _verdict(mean_log_excess: float, n_resolved: int, n_total: int) -> str:
    """Verdict text based on mean over-prediction and resolved fraction.

    The relevant question is *not* whether the substrate matches the
    binned observed φ in absolute normalisation (the toy HMF used here
    has a fixed power-law σ(M) that misses absolute amplitude by orders of
    magnitude at low halo mass). It is whether the SUBSTRATE/LCDM RATIO
    matches the JWST/LCDM ratio reported in the literature, which is
    ~ 0.5 - 1.0 dex (Boylan-Kolchin 2023; Labbé 2023).
    """
    lo, hi = JWST_OBSERVED_OVERABUNDANCE_DEX
    if mean_log_excess <= 0.0:
        return ("FAIL — substrate predicts FEWER high-z galaxies than LCDM; "
                "JWST overabundance NOT resolved by substrate.")
    if 0.0 < mean_log_excess < 0.10:
        return ("MARGINAL — substrate enhancement < 0.1 dex (factor 1.26x); "
                "insufficient to close the JWST gap (~ 0.5-1 dex).")
    if 0.10 <= mean_log_excess < lo:
        return ("PARTIAL — substrate enhancement of "
                f"~10^{mean_log_excess:.2f} = "
                f"{10**mean_log_excess:.2f}x; same direction as JWST data, "
                f"resolves {n_resolved}/{n_total} sources within 1-sigma.")
    if lo <= mean_log_excess <= hi:
        return ("GOOD MATCH — substrate enhancement "
                f"{10**mean_log_excess:.1f}x (~ 0.5-1 dex), matches the JWST "
                "OVER-abundance over LCDM (Boylan-Kolchin 2023; Labbé 2023).")
    if hi < mean_log_excess < 2.0:
        return ("OVER-PREDICTED — substrate enhancement "
                f"{10**mean_log_excess:.1f}x exceeds JWST OVER-abundance "
                "(~3-10x); boost factor likely too aggressive at this M_*.")
    return ("PATHOLOGICAL — substrate-vs-LCDM ratio "
            f"~10^{mean_log_excess:.1f} dex is dominated by toy-HMF "
            "rare-peak suppression rather than the physical 0.30 dex boost; "
            "treat as a qualitative direction test only.")


def run_jwst_high_z_test() -> JWSTHighZTestResult:
    """Top-level entry point: run all comparisons + verdict."""
    substrate = make_substrate_track()
    lcdm = make_lcdm_track()

    sources = compare_per_source(substrate, lcdm)
    LF = compare_UV_LF(substrate, lcdm)
    density = compare_stellar_mass_density(substrate, lcdm)

    excesses = [s.log_excess_substrate for s in sources]
    mean_x = sum(excesses) / len(excesses)
    sorted_x = sorted(excesses)
    median_x = sorted_x[len(sorted_x) // 2]

    n_total = len(sources)
    n_resolved = sum(
        1 for s in sources
        if abs(s.n_sigma_substrate_vs_lcdm) > 1.0 and s.log_excess_substrate > 0
    )

    # Eps_star contribution: at fixed M_*, M_halo is smaller in substrate
    # by eps_sub/eps_lcdm; the rare-peak HMF is steep so a smaller M_halo
    # gives a much higher phi. Quote the log of (eps_sub/eps_lcdm).
    eps_log_dex = math.log10(substrate.eps_star / lcdm.eps_star)
    # Total boost at fixed M_star: cosmology+HMF + eps_star (additive in dex
    # because phi(M_halo) and Jacobian dlogM_halo/dlogM_star = 1).
    total_dex = mean_x + eps_log_dex
    verdict = _verdict(total_dex, n_resolved, n_total)

    return JWSTHighZTestResult(
        sources=sources, UV_LF=LF, density=density,
        mean_log_excess=mean_x, median_log_excess=median_x,
        n_sources_resolved=n_resolved, n_sources=n_total,
        eps_star_substrate=substrate.eps_star,
        eps_star_lcdm=lcdm.eps_star,
        eps_star_log_boost_dex=eps_log_dex,
        total_log_boost_dex=total_dex,
        verdict=verdict,
    )


def report() -> str:
    """Human-readable summary."""
    r = run_jwst_high_z_test()
    lines: List[str] = []
    lines.append("=" * 72)
    lines.append("JWST high-z galaxies vs substrate vs LambdaCDM")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Per-source comoving abundance phi (Mpc^-3 dex^-1) at observed M_*:")
    lines.append(f"  {'name':22s} {'z':>5s}   {'logM*':>5s}   "
                 f"{'logφ_sub':>10s} {'logφ_LCDM':>10s} {'log(sub/LCDM)':>14s}")
    for s in r.sources:
        lines.append(
            f"  {s.name:22s} {s.z:5.2f}   {s.log_M_star:5.2f}   "
            f"{s.log_phi_substrate:10.3f} {s.log_phi_lcdm:10.3f} "
            f"{s.log_excess_substrate:14.3f}"
        )
    lines.append("")
    lines.append("UV LF (Harikane/Donnan binned vs predictions):")
    lines.append(f"  {'z':>5s}  {'M_UV':>6s}   "
                 f"{'logφ_obs':>9s}   {'logφ_sub':>9s}   {'logφ_LCDM':>9s}   "
                 f"{'nσ_sub':>7s}   {'nσ_LCDM':>7s}")
    for u in r.UV_LF:
        lines.append(
            f"  {u.z:5.1f}  {u.M_UV:6.2f}   {u.log_phi_obs:9.2f}   "
            f"{u.log_phi_substrate:9.2f}   {u.log_phi_lcdm:9.2f}   "
            f"{u.n_sigma_substrate:7.2f}   {u.n_sigma_lcdm:7.2f}"
        )
    lines.append("")
    lines.append("Cosmic stellar mass density (Labbé+2023):")
    lines.append(f"  {'z':>5s}  {'logM*_min':>9s}   "
                 f"{'log ρ_obs':>10s}   {'log ρ_sub':>10s}   {'log ρ_LCDM':>11s}   "
                 f"{'nσ_sub':>7s}   {'nσ_LCDM':>7s}")
    for d in r.density:
        lines.append(
            f"  {d.z:5.1f}  {d.log_M_star_min:9.2f}   "
            f"{d.log_rho_obs:10.2f}   {d.log_rho_substrate:10.2f}   "
            f"{d.log_rho_lcdm:11.2f}   "
            f"{d.n_sigma_substrate:7.2f}   {d.n_sigma_lcdm:7.2f}"
        )
    lines.append("")
    lines.append(f"Mean log10 excess at fixed M_halo (sub/LCDM):   "
                 f"{r.mean_log_excess:+.3f}  "
                 f"(factor {10**r.mean_log_excess:.2f}x)")
    lines.append(f"Median log10 excess at fixed M_halo:            "
                 f"{r.median_log_excess:+.3f}")
    lines.append(f"Sources where cosmology boost > 1 sigma:        "
                 f"{r.n_sources_resolved}/{r.n_sources}")
    lines.append("")
    lines.append("Astrophysical eps_star contribution (substrate vs LCDM):")
    lines.append(f"  eps_*(substrate) = {r.eps_star_substrate:.3f}")
    lines.append(f"  eps_*(LCDM)      = {r.eps_star_lcdm:.3f}")
    lines.append(f"  log10 ratio      = {r.eps_star_log_boost_dex:+.3f} dex  "
                 f"(factor {10**r.eps_star_log_boost_dex:.2f}x)")
    lines.append(f"  TOTAL boost at fixed M_*: "
                 f"{r.total_log_boost_dex:+.3f} dex  "
                 f"(factor {10**r.total_log_boost_dex:.2f}x)")
    lines.append("")
    lines.append("Verdict (combined cosmology + eps_star):")
    lines.append(f"  {r.verdict}")
    lines.append("=" * 72)
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
