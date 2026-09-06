"""Nucleon magnetic moments μ_p, μ_n from substrate Möbius spin structure.

Physics summary
---------------
The proton and neutron magnetic moments are amongst the cleanest QCD
predictions of the constituent (additive) quark model.  In the SU(2)-flavour
ground-state baryon wavefunction (J=1/2, total spin from three S=1/2
constituents), Wigner's spin-flavour algebra gives the famous formulas

    μ_p = (4/3) μ_u − (1/3) μ_d                                          (1)
    μ_n = (4/3) μ_d − (1/3) μ_u                                          (2)

where μ_q = e_q ℏ /(2 m_q c) is a Dirac magnetic moment with charge
e_u = +2/3 e and e_d = −1/3 e.  In QCD the empirical fit requires
m_constituent ≈ 336 MeV, giving μ_p = 2.79 μ_N, μ_n = −1.86 μ_N (vs −1.91
empirical) — a 3% prediction with one input.

In the stiff-medium substrate model we are NOT allowed to import the QCD
constituent mass 336 MeV.  Three kinks live on Möbius half-flux U(1)
bundles inside the baryon and carry Dirac magnetic moments

    μ_K = e_K ℏ / (2 m_K_eff c)                                          (3)

where:
  * e_K is the kink's electric charge in units of e — set by the half-flux
    bundle topology (§18.10) and the colour pair-up that produces u/d
    quark charges (+2/3, −1/3) in the SU(3)-inherited substrate (§18.49).
  * m_K_eff is the constituent (in-baryon) kink mass — the substrate
    analogue of QCD's constituent quark mass.

The same SU(2)-flavour algebra (1)–(2) carries through, since it depends
only on the spin-flavour structure of the J=1/2 baryon wavefunction.  The
substrate inherits this algebra from §18.49's SU(3) gauge inheritance plus
the spin-½ statistics of the half-flux bundle (§18.10).

Substrate-derived candidate effective masses
--------------------------------------------
The N-Δ splitting fixes m_K_eff = √σ ≈ 424 MeV (chromomagnetic_substrate.py
§18.61.1).  With this mass the proton magnetic moment comes out at 2.21 μ_N
(−21%), as the larger-than-QCD constituent mass shrinks each μ_q.

Several other substrate-natural mass candidates are tested:

  (a) m_K_eff = √σ ≈ 424 MeV         — N-Δ splitting (geometric ℏc/R₀)
  (b) m_K = m_p / 3 ≈ 313 MeV         — equal-share proton-third
  (c) m_K = √σ × 3 / 4 ≈ 318 MeV       — 3/4 of the geometric mass
                                         (Y-junction kinematic factor)
  (d) m_K from the linear potential ⟨r⟩_(σ r) virial ≈ 343 MeV
  (e) m_K from the Foot calibration M = 313.86 MeV (lepton-vertex scale)

The cleanest substrate-only choice that reproduces both μ_p and μ_n at
the few-percent level is then identified.

Isoscalar sum
-------------
The SU(2) prediction (1)+(2) gives

    μ_p + μ_n = (4/3 − 1/3)(μ_u + μ_d) = (μ_u + μ_d)
              = (2/3 − 1/3) μ_N × (m_p / m_K_eff)
              = (1/3) × m_p / m_K_eff × μ_N                              (4)

so μ_p + μ_n is set entirely by m_p / m_K_eff, with no spin-flavour mixing.
Empirical: μ_p + μ_n = 0.880 μ_N.  Inverting (4) gives

    m_K_eff (isoscalar) = m_p / (3 × 0.880) = 938 / 2.64 ≈ 355 MeV.

This is the substrate-natural target for the constituent kink mass when
matched to the isoscalar moment alone — independent of any spin algebra
beyond the SU(2) ground-state structure.

Isovector difference
--------------------
Subtracting (2) from (1):

    μ_p − μ_n = (4/3 + 1/3)(μ_u − μ_d) = (5/3)(μ_u − μ_d)
              = (5/3)(+2/3 − (−1/3)) μ_N × (m_p / m_K_eff)
              = (5/3) × m_p / m_K_eff × μ_N                              (5)

Empirical: μ_p − μ_n = 4.706 μ_N → m_K_eff (isovector) = 5 m_p / (3 × 4.706)
                                                       ≈ 332 MeV.

The isoscalar and isovector channels prefer slightly different masses
(355 vs 332 MeV); the empirical "best-fit" QCD value 336 MeV sits between
them.  The substrate's job is to predict where m_K_eff actually lies.

References
----------
spec §18.4 (half-flux), §18.10 (Möbius connection), §18.30 (vertex), §18.49
(SU(3) inheritance), §18.49.5 (string tension), §18.60.2 (R₀ = 1/√σ),
§18.61.1 (N-Δ splitting and constituent kink mass).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Substrate primitives at the QCD scale (inherited from chromomagnetic_substrate)
# ---------------------------------------------------------------------------

SIGMA_QCD_GEV2: Final[float] = 0.180   # GeV²  (string tension, §18.49.5)
XI_QCD_FM: Final[float] = 0.2          # fm    (coherence length at QCD scale)
HBARC_GEVFM: Final[float] = 0.197326980
XI_QCD_INV_GEV: Final[float] = XI_QCD_FM / HBARC_GEVFM   # ≈ 1.0135 GeV⁻¹

# Foot calibration mass (lepton-vertex scale, used elsewhere in the model)
M_FOOT_MEV: Final[float] = 313.86

# ---------------------------------------------------------------------------
# Empirical anchors
# ---------------------------------------------------------------------------

M_PROTON_MEV: Final[float] = 938.272
MU_P_EMPIRICAL: Final[float] = +2.7928473       # in units of μ_N
MU_N_EMPIRICAL: Final[float] = -1.9130427       # in units of μ_N

# Charges of the kink-pairs that produce u, d quarks (§18.49 SU(3) inheritance)
E_U: Final[float] = +2.0 / 3.0    # up-kink charge in units of e
E_D: Final[float] = -1.0 / 3.0    # down-kink charge in units of e


# ---------------------------------------------------------------------------
# Substrate-natural candidate effective masses for the kink
# ---------------------------------------------------------------------------

def m_kink_geometric_MeV(sigma_GeV2: float = SIGMA_QCD_GEV2) -> float:
    """In-baryon constituent kink mass m_K_eff = √σ.

    Geometric: ℏc / R₀ with R₀ = 1/√σ (Y-junction equilibrium length, §18.60.2).
    This is the constituent kink mass that fixes Δm_NΔ in §18.61.1.

    For σ = 0.180 GeV² this returns 424.3 MeV.
    """
    return 1000.0 * math.sqrt(sigma_GeV2)


def m_kink_proton_third_MeV(m_proton_MeV: float = M_PROTON_MEV) -> float:
    """Equal-share proton-third constituent: m_K = m_p / 3 ≈ 312.8 MeV.

    The simplest substrate parametrisation: each of the three kinks carries
    1/3 of the proton's mass.  Comes from the constraint m_p = 3 m_K_eff
    if the kink kinetic energy and binding cancel — the leading-order
    bag/constituent picture.
    """
    return m_proton_MeV / 3.0


def m_kink_y_junction_kinematic_MeV(sigma_GeV2: float = SIGMA_QCD_GEV2) -> float:
    """Y-junction kinematic mass m_K = (3/4) √σ.

    The Y-junction baryon has 3 string arms of length R₀ = 1/√σ meeting at
    a vertex.  In the rotating Y-junction limit, the effective constituent
    mass acquires a kinematic factor of 3/4 from the relative motion of
    each kink-end relative to the centre-of-mass (1/2 from the standard
    rotating-string and 3/2 from the three-arm geometry, giving 3/4 net).

    Returns 318 MeV for σ = 0.180 GeV².
    """
    return 1000.0 * (3.0 / 4.0) * math.sqrt(sigma_GeV2)


def m_kink_linear_potential_MeV(sigma_GeV2: float = SIGMA_QCD_GEV2) -> float:
    """Constituent mass from the Airy-ground-state expectation in V(r)=σr.

    For a single kink in V(r)=σr with bare-mass-zero approximation (relativistic
    light kink), the ground-state energy is

        E_0 = a₁ × (σ²)^(1/3) × (1 / (2 m_bare))^{−1/3}

    For the substrate we set m_bare = √σ (the geometric scale), giving

        E_0 ≈ |a₁| × σ^(1/2) ≈ 2.338 × σ^(1/2)

    where a₁ ≈ −2.338 is the first zero of the Airy function.  But this is
    a binding energy, not an effective mass.  More usefully, the radial
    kinetic energy ⟨T⟩ = E_0 / 3 = 0.78 √σ ≈ 330 MeV plays the role of the
    in-baryon kinetic constituent mass.

    Returns 330 MeV for σ = 0.180 GeV².
    """
    a1 = 2.338107410459767     # |first Airy zero|
    return 1000.0 * (a1 / 3.0) * math.sqrt(sigma_GeV2)


def m_kink_foot_MeV() -> float:
    """Foot calibration mass M_FOOT = 313.86 MeV (lepton-vertex scale).

    This is the same calibration mass used in the lepton sector
    (foot_phase_derivation.py, mobius_dirac_vertex.py).  Including it here
    asks: does the lepton vertex mass also work for the baryon kink?

    Returns 313.86 MeV.
    """
    return M_FOOT_MEV


def m_kink_isoscalar_target_MeV(
    m_proton_MeV: float = M_PROTON_MEV,
    mu_p_plus_n: float = MU_P_EMPIRICAL + MU_N_EMPIRICAL,
) -> float:
    """Mass that exactly reproduces the empirical isoscalar moment μ_p + μ_n.

    From eq. (4):  μ_p + μ_n = (1/3) m_p / m_K_eff μ_N
    =>  m_K_eff = m_p / (3 × (μ_p + μ_n)).

    Returns 355 MeV for the empirical sum 0.880 μ_N — purely diagnostic,
    NOT a substrate prediction.  Reported so we can compare which substrate
    candidate is closest.
    """
    return m_proton_MeV / (3.0 * mu_p_plus_n)


def m_kink_isovector_target_MeV(
    m_proton_MeV: float = M_PROTON_MEV,
    mu_p_minus_n: float = MU_P_EMPIRICAL - MU_N_EMPIRICAL,
) -> float:
    """Mass that exactly reproduces the empirical isovector moment μ_p − μ_n.

    From eq. (5):  μ_p − μ_n = (5/3) m_p / m_K_eff μ_N
    =>  m_K_eff = 5 m_p / (3 × (μ_p − μ_n)).

    Returns ~332 MeV for the empirical difference 4.706 μ_N — diagnostic only.
    """
    return 5.0 * m_proton_MeV / (3.0 * mu_p_minus_n)


# ---------------------------------------------------------------------------
# Magnetic-moment formulas (SU(2)-flavour ground-state baryon wavefunction)
# ---------------------------------------------------------------------------

def quark_moment_in_muN(charge_in_e: float, m_q_MeV: float) -> float:
    """Single-kink Dirac magnetic moment μ_q in units of the nuclear magneton μ_N.

        μ_q = e_q ℏ / (2 m_q c)
            = (e_q / m_q) × m_p μ_N      [since μ_N ≡ eℏ/(2 m_p c)]
            = e_q × m_p / m_q  (in units of μ_N)

    Args:
        charge_in_e: Kink charge in units of the elementary charge e
            (e.g. +2/3 for u-pair, −1/3 for d-pair).
        m_q_MeV: Constituent kink mass [MeV].

    Returns:
        μ_q expressed in units of μ_N.
    """
    return charge_in_e * (M_PROTON_MEV / m_q_MeV)


def proton_moment_SU2(m_q_MeV: float) -> float:
    """Proton magnetic moment μ_p in units of μ_N from eq. (1) (SU(2) symmetric).

    μ_p = (4/3) μ_u − (1/3) μ_d
        = (4/3)(2/3) m_p/m_q + (1/3)(1/3) m_p/m_q
        = (8/9 + 1/9) (m_p / m_q)
        = (m_p / m_q)  μ_N

    The famous result: in the SU(2)-symmetric limit, μ_p = m_p / m_q μ_N.
    """
    mu_u = quark_moment_in_muN(E_U, m_q_MeV)
    mu_d = quark_moment_in_muN(E_D, m_q_MeV)
    return (4.0 / 3.0) * mu_u - (1.0 / 3.0) * mu_d


def neutron_moment_SU2(m_q_MeV: float) -> float:
    """Neutron magnetic moment μ_n in units of μ_N from eq. (2) (SU(2) symmetric).

    μ_n = (4/3) μ_d − (1/3) μ_u
        = (4/3)(−1/3) m_p/m_q − (1/3)(2/3) m_p/m_q
        = (−4/9 − 2/9) (m_p / m_q)
        = −(2/3) (m_p / m_q)  μ_N

    The SU(2) prediction μ_p / μ_n = −3/2 ≈ −1.500 (empirical: −1.460).
    """
    mu_u = quark_moment_in_muN(E_U, m_q_MeV)
    mu_d = quark_moment_in_muN(E_D, m_q_MeV)
    return (4.0 / 3.0) * mu_d - (1.0 / 3.0) * mu_u


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NucleonMomentResult:
    """One mass-candidate evaluation of (μ_p, μ_n).

    Attributes:
        label: Human-readable identifier of the mass candidate.
        m_q_MeV: Effective constituent kink mass [MeV].
        mu_p: Predicted proton moment [μ_N].
        mu_n: Predicted neutron moment [μ_N].
        mu_sum: μ_p + μ_n (isoscalar) [μ_N].
        mu_diff: μ_p − μ_n (isovector) [μ_N].
        ratio: μ_p / μ_n (sign + magnitude diagnostic; SU(2) gives −3/2).
        err_p: (μ_p_pred − μ_p_obs) / |μ_p_obs|.
        err_n: (μ_n_pred − μ_n_obs) / |μ_n_obs|.
        err_sum: (sum_pred − sum_obs) / |sum_obs|.
        err_diff: (diff_pred − diff_obs) / |diff_obs|.
        max_abs_err: max(|err_p|, |err_n|).
    """
    label: str
    m_q_MeV: float
    mu_p: float
    mu_n: float
    mu_sum: float
    mu_diff: float
    ratio: float
    err_p: float
    err_n: float
    err_sum: float
    err_diff: float
    max_abs_err: float


def evaluate_candidate(label: str, m_q_MeV: float) -> NucleonMomentResult:
    """Compute (μ_p, μ_n) and errors for one m_q candidate."""
    mu_p = proton_moment_SU2(m_q_MeV)
    mu_n = neutron_moment_SU2(m_q_MeV)
    mu_sum = mu_p + mu_n
    mu_diff = mu_p - mu_n
    ratio = mu_p / mu_n if abs(mu_n) > 1e-12 else float("nan")

    err_p = (mu_p - MU_P_EMPIRICAL) / abs(MU_P_EMPIRICAL)
    err_n = (mu_n - MU_N_EMPIRICAL) / abs(MU_N_EMPIRICAL)
    sum_obs = MU_P_EMPIRICAL + MU_N_EMPIRICAL
    diff_obs = MU_P_EMPIRICAL - MU_N_EMPIRICAL
    err_sum = (mu_sum - sum_obs) / abs(sum_obs)
    err_diff = (mu_diff - diff_obs) / abs(diff_obs)

    return NucleonMomentResult(
        label=label,
        m_q_MeV=m_q_MeV,
        mu_p=mu_p,
        mu_n=mu_n,
        mu_sum=mu_sum,
        mu_diff=mu_diff,
        ratio=ratio,
        err_p=err_p,
        err_n=err_n,
        err_sum=err_sum,
        err_diff=err_diff,
        max_abs_err=max(abs(err_p), abs(err_n)),
    )


# ---------------------------------------------------------------------------
# Sweep across all substrate-natural mass candidates
# ---------------------------------------------------------------------------

def sweep_mass_candidates(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    m_proton_MeV: float = M_PROTON_MEV,
) -> list[NucleonMomentResult]:
    """Evaluate (μ_p, μ_n) for each substrate-natural mass candidate.

    Returns a list ordered by max(|err_p|, |err_n|) ascending — best fit first.
    """
    candidates: list[tuple[str, float]] = [
        ("(a) m_K = √σ  [N-Δ geometric]",
            m_kink_geometric_MeV(sigma_GeV2)),
        ("(b) m_K = m_p/3  [proton-third]",
            m_kink_proton_third_MeV(m_proton_MeV)),
        ("(c) m_K = (3/4)√σ  [Y-junction kinematic]",
            m_kink_y_junction_kinematic_MeV(sigma_GeV2)),
        ("(d) m_K = (a₁/3)√σ  [Airy linear-V kinetic]",
            m_kink_linear_potential_MeV(sigma_GeV2)),
        ("(e) m_K = M_Foot  [lepton-vertex scale]",
            m_kink_foot_MeV()),
        # Diagnostic targets (NOT substrate predictions — exact fits to data)
        ("(diag-isoscalar) m_K matches μ_p+μ_n",
            m_kink_isoscalar_target_MeV(m_proton_MeV)),
        ("(diag-isovector) m_K matches μ_p−μ_n",
            m_kink_isovector_target_MeV(m_proton_MeV)),
    ]

    results = [evaluate_candidate(label, m) for label, m in candidates]
    return results


# ---------------------------------------------------------------------------
# Substrate prediction: pick the mass-channel with both nucleon moments < 5%
# ---------------------------------------------------------------------------

def best_substrate_candidate(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    m_proton_MeV: float = M_PROTON_MEV,
    threshold: float = 0.05,
) -> NucleonMomentResult | None:
    """Return the lowest-error substrate candidate that satisfies |err| < threshold.

    "Substrate" here excludes the (diag-…) targets, which are not predictions
    but exact fits to the data.

    Args:
        sigma_GeV2: String tension [GeV²].
        m_proton_MeV: Proton mass [MeV].
        threshold: Fractional-error threshold (0.05 = 5%).

    Returns:
        :class:`NucleonMomentResult` with the smallest max(|err_p|, |err_n|)
        among substrate candidates whose max error is below ``threshold``.
        ``None`` if no substrate candidate passes.
    """
    results = sweep_mass_candidates(sigma_GeV2, m_proton_MeV)
    substrate_only = [r for r in results if not r.label.startswith("(diag")]

    passing = [r for r in substrate_only if r.max_abs_err < threshold]
    if not passing:
        return None
    return min(passing, key=lambda r: r.max_abs_err)


# ---------------------------------------------------------------------------
# Substrate prediction: closed-form isoscalar moment from σ and m_p
# ---------------------------------------------------------------------------

def isoscalar_moment_substrate(
    m_q_MeV: float,
    m_proton_MeV: float = M_PROTON_MEV,
) -> float:
    """μ_p + μ_n = (1/3) m_p / m_K_eff  μ_N    (eq. 4).

    For m_K_eff = √σ ≈ 424 MeV → μ_p + μ_n ≈ 0.737 μ_N (vs empirical 0.880,
    error −16%).  For m_K_eff = m_p/3 ≈ 313 MeV → exactly 1.000 μ_N.

    Args:
        m_q_MeV: Constituent kink mass [MeV].
        m_proton_MeV: Proton mass [MeV].

    Returns:
        Isoscalar nucleon moment in units of μ_N.
    """
    return (m_proton_MeV / m_q_MeV) / 3.0


def isovector_moment_substrate(
    m_q_MeV: float,
    m_proton_MeV: float = M_PROTON_MEV,
) -> float:
    """μ_p − μ_n = (5/3) m_p / m_K_eff  μ_N    (eq. 5).

    Args:
        m_q_MeV: Constituent kink mass [MeV].
        m_proton_MeV: Proton mass [MeV].

    Returns:
        Isovector nucleon moment in units of μ_N.
    """
    return (5.0 / 3.0) * (m_proton_MeV / m_q_MeV)


# ---------------------------------------------------------------------------
# Top-level summary container
# ---------------------------------------------------------------------------

@dataclass
class NucleonMomentsSummary:
    """All mass-candidate sweeps plus the chosen best substrate candidate."""
    sigma_GeV2: float
    m_proton_MeV: float
    candidates: list[NucleonMomentResult]
    best_substrate: NucleonMomentResult | None
    threshold: float


def compute_nucleon_moments_summary(
    sigma_GeV2: float = SIGMA_QCD_GEV2,
    m_proton_MeV: float = M_PROTON_MEV,
    threshold: float = 0.05,
) -> NucleonMomentsSummary:
    """Compute the full mass-candidate sweep and identify the best substrate fit.

    Args:
        sigma_GeV2: String tension [GeV²].
        m_proton_MeV: Proton mass [MeV].
        threshold: Fractional-error threshold for "best substrate" qualification.

    Returns:
        :class:`NucleonMomentsSummary` containing all candidates and the
        identified best (or None if none passes the threshold).
    """
    return NucleonMomentsSummary(
        sigma_GeV2=sigma_GeV2,
        m_proton_MeV=m_proton_MeV,
        candidates=sweep_mass_candidates(sigma_GeV2, m_proton_MeV),
        best_substrate=best_substrate_candidate(sigma_GeV2, m_proton_MeV, threshold),
        threshold=threshold,
    )


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

def format_candidate_row(r: NucleonMomentResult) -> str:
    """One-line summary for a single mass candidate."""
    return (
        f"  {r.label:<48}  m_K={r.m_q_MeV:>7.2f} MeV  "
        f"μ_p={r.mu_p:+7.4f}  μ_n={r.mu_n:+7.4f}  "
        f"err(p)={100.0*r.err_p:+7.2f}%  err(n)={100.0*r.err_n:+7.2f}%"
    )


def print_full_report(summary: NucleonMomentsSummary) -> None:
    """Print the full nucleon-moment analysis (sweep + best substrate)."""
    print("=" * 96)
    print(" Substrate prediction of nucleon magnetic moments μ_p, μ_n")
    print(" (SU(2)-flavour additive-quark formulas, m_K_eff from substrate)")
    print("=" * 96)
    print()
    print(f"  Substrate inputs:  σ = {summary.sigma_GeV2:.4f} GeV²,  "
          f"m_p = {summary.m_proton_MeV:.3f} MeV")
    print(f"  Empirical anchors: μ_p = {MU_P_EMPIRICAL:+.4f} μ_N,  "
          f"μ_n = {MU_N_EMPIRICAL:+.4f} μ_N,")
    print(f"                     μ_p+μ_n = {MU_P_EMPIRICAL + MU_N_EMPIRICAL:+.4f} μ_N "
          f"(isoscalar),")
    print(f"                     μ_p−μ_n = {MU_P_EMPIRICAL - MU_N_EMPIRICAL:+.4f} μ_N "
          f"(isovector),")
    print(f"                     μ_p/μ_n = {MU_P_EMPIRICAL / MU_N_EMPIRICAL:+.4f}  "
          f"(SU(2): −3/2 = −1.500)")
    print()
    print("-" * 96)
    print(" Mass-candidate sweep")
    print("-" * 96)
    print()

    # Sort by max abs error so the best one comes first
    ordered = sorted(summary.candidates, key=lambda r: r.max_abs_err)
    for r in ordered:
        print(format_candidate_row(r))
    print()

    print("-" * 96)
    print(" Isoscalar / isovector channel decomposition")
    print("-" * 96)
    print()
    print(f"  {'candidate':<48}  {'m_K [MeV]':>10}  "
          f"{'μ_p+μ_n':>10}  {'err':>9}  {'μ_p-μ_n':>10}  {'err':>9}")
    for r in ordered:
        print(
            f"  {r.label:<48}  {r.m_q_MeV:>10.2f}  "
            f"{r.mu_sum:>+10.4f}  {100.0*r.err_sum:>+8.2f}%  "
            f"{r.mu_diff:>+10.4f}  {100.0*r.err_diff:>+8.2f}%"
        )
    print()

    print("-" * 96)
    print(" Best substrate prediction (excluding diagnostic exact-fit targets)")
    print("-" * 96)
    print()
    if summary.best_substrate is None:
        print(f"  No substrate candidate satisfies max|err| < "
              f"{100.0 * summary.threshold:.1f}%.")
        # Report the closest substrate candidate anyway
        substrate_only = [
            r for r in summary.candidates if not r.label.startswith("(diag")
        ]
        closest = min(substrate_only, key=lambda r: r.max_abs_err)
        print(f"  Closest substrate candidate:  {closest.label}")
        print(f"    m_K_eff           = {closest.m_q_MeV:.2f} MeV")
        print(f"    μ_p (pred / obs)  = {closest.mu_p:+.4f} / {MU_P_EMPIRICAL:+.4f}  "
              f"(err {100.0*closest.err_p:+.2f}%)")
        print(f"    μ_n (pred / obs)  = {closest.mu_n:+.4f} / {MU_N_EMPIRICAL:+.4f}  "
              f"(err {100.0*closest.err_n:+.2f}%)")
        print(f"    μ_p/μ_n           = {closest.ratio:+.4f}  (SU(2): −1.500)")
    else:
        b = summary.best_substrate
        print(f"  WINNER:  {b.label}")
        print(f"    m_K_eff           = {b.m_q_MeV:.2f} MeV")
        print(f"    μ_p (pred / obs)  = {b.mu_p:+.4f} / {MU_P_EMPIRICAL:+.4f}  "
              f"(err {100.0*b.err_p:+.2f}%)")
        print(f"    μ_n (pred / obs)  = {b.mu_n:+.4f} / {MU_N_EMPIRICAL:+.4f}  "
              f"(err {100.0*b.err_n:+.2f}%)")
        print(f"    μ_p/μ_n           = {b.ratio:+.4f}  (SU(2): −1.500)")
        print(f"    threshold         = {100.0 * summary.threshold:.1f}%")
    print()

    print("-" * 96)
    print(" Sign structure (cleanest test — purely from spin-flavour algebra)")
    print("-" * 96)
    print()
    print("  All substrate candidates predict μ_p > 0 and μ_n < 0 with μ_p/μ_n = −3/2.")
    print("  The sign and the −3/2 ratio are FORCED by the SU(2) spin-flavour structure")
    print("  and the kink-pair charges (+2/3, −1/3) — no fit, no QCD α_s, no constituent")
    print("  mass input.  The substrate inherits these from the §18.49 SU(3) gauge")
    print("  inheritance plus the half-flux Möbius spin-½ statistics (§18.10).")
    print()
    print("  Empirical μ_p/μ_n = −1.460 (vs SU(2) prediction −1.500): the small")
    print("  4% deviation is the well-known SU(6) → SU(3)-flavour symmetry breaking,")
    print("  beyond the SU(2)-symmetric ground-state baryon wavefunction used here.")
    print()


if __name__ == "__main__":
    print_full_report(compute_nucleon_moments_summary())
