"""Does the substrate K(ξ) running reproduce the QCD α_s(Q²) running?

Physics question
----------------
The substrate stiffness runs with scale.  Anchored to two physical points

    K(ξ_e   = 3.86e-13 m) = ℏc/ξ_e⁴           = 1.42e24 Pa     (electron scale)
    K(ξ_QCD = 2.0e-16 m) = σ_QCD_SI / σ_lat  = 2.87e5  Pa     (QCD scale)

a pure-power-law fit gives

    K(ξ) = K_e × (ξ_e/ξ)^a    with    a = -5.69                            (1)

(Sign convention: a < 0 because K is SMALLER at the QCD scale than at the
electron scale.  K reaches its physical Planck value K_P = c⁷/(ℏG²) ≈
4.6e113 Pa only at ξ = l_Planck, far below ξ_QCD; the run from ξ_e to
ξ_QCD passes through a minimum.)

At the QCD scale the §18.61.1 dimensionless Möbius coupling

    α_M(ξ) ≡ σ(ξ)[GeV²] × ξ²[GeV⁻²]                                       (2)

evaluates to ≈ 0.185, which played the role of α_s in the +6.8% N-Δ
chromomagnetic prediction (chromomagnetic_substrate.py).

The natural follow-up: does α_M evaluated at OTHER scales reproduce the
QCD α_s(Q²) running?  In particular, does it match α_s(M_Z) = 0.118?

QCD reference (PDG 2024)
------------------------
    α_s(M_Z = 91.2 GeV) = 0.118
    α_s(10 GeV)         ≈ 0.18
    α_s(1 GeV)          ≈ 0.45
    α_s(M_τ = 1.78 GeV) ≈ 0.32

Asymptotic freedom: α_s decreases monotonically as Q increases.

What this module computes
-------------------------
We test FIVE candidate dimensionless substrate couplings against α_s(Q):

  A. α_M    = σ × ξ²        [§18.61.1 Möbius coupling, in (GeV²)(GeV⁻²)]
  B. α_lat  = σ_SI × ξ² / ℏc [substrate-natural fully dimensionless]
  C. α_C    = ℏc / (K × ξ⁴)  [action-quantum residual; = 1 if ℏ = Kξ⁴/c
                              holds at every scale]
  D. α_D    = α_M(QCD) × (Q_QCD/Q)^(a-2)
                              [a power-law DECREASE forced through QCD anchor]
  E. α_E    = QCD-style log: 1/α_E(Q) = 1/α_E(Q_0) + b ln(Q/Q_0), b chosen
              from the substrate's own α_M values at Q=1 GeV and Q=M_Z.

Findings (see verdict at end of run)
------------------------------------
• Candidate A (the §18.61.1 coupling) DECREASES with Q (since a < 0, σ also
  shrinks below the electron scale, faster than 1/ξ² grows).  But it
  decreases as a STEEP POWER LAW, not the gentle logarithm of QCD.

• Candidate B is identically σ_lat ≈ 0.51 at every scale (scale-INVARIANT).

• Candidate C grows with Q as a power law — the WRONG direction from
  asymptotic freedom.

• Candidates D and E are constrained fits, not derivations, and demonstrate
  that no LOGARITHMIC running emerges from a power-law K(ξ).

The honest result: the substrate's α_M ≈ 0.185 at the QCD scale is a
MAGNITUDE COINCIDENCE.  A power-law K(ξ) cannot reproduce QCD's logarithmic
β-function, so the running of α_M is not the running of α_s — and at M_Z
α_M is many orders of magnitude smaller than α_s.

References
----------
spec §18.61.1 (α_M = σ ξ² ≈ 0.185 at QCD scale);
src/stiff_medium/substrate_rg_running.py (K(ξ) power-law, a ≈ -5.69);
src/stiff_medium/regge_spectrum.py (σ_QCD = 0.180 GeV²);
src/stiff_medium/chromomagnetic_substrate.py (Möbius-coupling N-Δ usage);
PDG 2024 (α_s(Q²) reference data).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

# Reuse substrate constants and the fitted RG running
from stiff_medium.substrate_rg_running import (
    HBAR_SI,
    C_SI,
    GEV_TO_J,
    FM_TO_M,
    XI_ELECTRON,
    K_ELECTRON,
    SIGMA_LATTICE_NATURAL,
    build_anchors,
    fit_triple_constrained_n,
)


# ---------------------------------------------------------------------------
# Constants for unit conversion
# ---------------------------------------------------------------------------

HBARC_GEVFM: Final[float] = 0.197326980          # GeV·fm
HBARC_J_M: Final[float] = HBAR_SI * C_SI         # J·m

# Substrate-derived running parameters (set lazily; from triple anchor fit)
_TRIPLE_FIT_CACHE: dict[str, float] | None = None


def _triple_fit_params() -> tuple[float, float]:
    """Run the triple-anchor fit once and cache (a, n) of the K-running.

    Returns:
        (a_substrate, n_substrate) — the substrate's fitted RG exponent and
        power-law index.  In the convention K(ξ)^(1-n) = K_e^(1-n) +
        a(1-n) ln(ξ/ξ_e), n=1 case reduces to K = K_e exp(-a ln(ξ/ξ_e)) =
        K_e (ξ/ξ_e)^(-a) = K_e (ξ_e/ξ)^a.
    """
    global _TRIPLE_FIT_CACHE
    if _TRIPLE_FIT_CACHE is None:
        anchor_e, anchor_qcd, anchor_planck = build_anchors()
        n_best, a_best, _ = fit_triple_constrained_n(
            anchor_e, anchor_qcd, anchor_planck
        )
        _TRIPLE_FIT_CACHE = {"n": n_best, "a": a_best}
    return _TRIPLE_FIT_CACHE["a"], _TRIPLE_FIT_CACHE["n"]


# Anchor: σ_QCD = 0.180 GeV² at ξ_QCD = 0.2 fm gives K_QCD via σ_SI = σ_lat K
_XI_QCD: Final[float] = 0.2 * FM_TO_M
_SIGMA_QCD_GEV2: Final[float] = 0.180
_SIGMA_QCD_SI: Final[float] = _SIGMA_QCD_GEV2 * GEV_TO_J ** 2 / HBARC_J_M
_K_QCD_TARGET: Final[float] = _SIGMA_QCD_SI / SIGMA_LATTICE_NATURAL  # ≈ 2.87e5 Pa


def _running_exponent_n1() -> float:
    """Compute the exponent of K(ξ) = K_e (ξ_e/ξ)^a for n=1.

    Uses the (electron, QCD) anchor pair.  Returns a as a SIGNED number;
    a < 0 means K decreases as ξ decreases below ξ_e.
    """
    # K_QCD = K_e * (xi_e/xi_QCD)^a  =>  a = ln(K_QCD/K_e) / ln(xi_e/xi_QCD)
    return math.log(_K_QCD_TARGET / K_ELECTRON) / math.log(XI_ELECTRON / _XI_QCD)


# ---------------------------------------------------------------------------
# Substrate K(ξ) and σ(ξ) running
# ---------------------------------------------------------------------------

def K_at(xi_m: float) -> float:
    """Substrate stiffness K [Pa] at scale ξ [m].

    Pure n=1 power law fitted between (electron, QCD) anchors:

        K(ξ) = K_e × (ξ_e / ξ)^a    with    a = ln(K_QCD/K_e)/ln(ξ_e/ξ_QCD)

    With K_e = 1.42e24 Pa and K_QCD = 2.87e5 Pa, a ≈ -5.69, so K decreases
    as ξ shrinks below the electron scale.

    Args:
        xi_m: Length scale [m].

    Returns:
        K [Pa].
    """
    a = _running_exponent_n1()
    return K_ELECTRON * (XI_ELECTRON / xi_m) ** a


def K_at_triple(xi_m: float) -> float:
    """Substrate K [Pa] from the triple-anchor (a, n) fit (alternative to power-law).

    Args:
        xi_m: Length scale [m].

    Returns:
        K [Pa] using n_triple, a_triple from the simultaneous (electron, QCD,
        Planck) anchor fit.  In practice the triple-anchor n is very close
        to 1 so this gives nearly the same K as :func:`K_at`.
    """
    a, n = _triple_fit_params()
    t = math.log(xi_m / XI_ELECTRON)
    if abs(n - 1.0) < 1e-10:
        return K_ELECTRON * math.exp(-a * t)
    val = K_ELECTRON ** (1.0 - n) + a * (1.0 - n) * t
    if val <= 0.0:
        return float("inf")
    return val ** (1.0 / (1.0 - n))


def sigma_GeV2_at(xi_m: float) -> float:
    """Substrate string tension σ [GeV²] at scale ξ [m].

    σ_SI(ξ) = σ_lat × K(ξ) [J/m] is converted to GeV² by σ × ℏc / GeV²
    (using the natural-units identity 1 GeV² ↔ 1 GeV² × (1/ℏc) J/m).

    Args:
        xi_m: Length scale [m].

    Returns:
        σ in GeV².
    """
    K_xi = K_at(xi_m)
    sigma_SI = SIGMA_LATTICE_NATURAL * K_xi              # J/m
    sigma_GeV2 = sigma_SI * HBARC_J_M / GEV_TO_J ** 2    # → GeV²
    return sigma_GeV2


# ---------------------------------------------------------------------------
# Candidate dimensionless substrate couplings
# ---------------------------------------------------------------------------

def alpha_M_naive(xi_m: float) -> float:
    """Candidate A: α_M = σ × ξ² (σ in GeV², ξ in GeV⁻¹).

    This is the §18.61.1 "Möbius coupling".  At ξ = ξ_QCD = 0.2 fm with
    σ = 0.180 GeV² it equals 0.185.

    Scale-dependence: σ(ξ) ∝ K(ξ) ∝ ξ^(-a) and ξ² ∝ ξ², so

        α_M ∝ ξ^(2 - a)    with    a ≈ -5.69
             ∝ ξ^(2 - (-5.69))
             ∝ ξ^(+7.69)
             ∝ Q^(-7.69)

    So α_M DECREASES as Q increases.  Direction matches QCD asymptotic
    freedom; magnitude / power-law shape does not.

    Args:
        xi_m: Length scale [m].

    Returns:
        Dimensionless σ × ξ².
    """
    sigma = sigma_GeV2_at(xi_m)                # GeV²
    xi_invGeV = xi_m / HBARC_GEVFM / FM_TO_M   # ξ [GeV⁻¹]
    return sigma * xi_invGeV ** 2


def alpha_lattice_invariant() -> float:
    """Candidate B (idealised): α_lat ≡ σ_lat = 0.51, scale-INVARIANT by construction.

    If the substrate identity ℏ = K ξ⁴ / c held at EVERY scale (i.e. if K(ξ)
    ran exactly as ℏc/ξ⁴ rather than as K_e (ξ_e/ξ)^a with a ≈ -5.69), then

        σ × ξ² / (ℏc) = σ_lat × K × ξ⁴ / (ℏc × ξ²) × 1/ξ²|_{ℏ-id} = σ_lat,

    i.e. the dimensionless string tension σ_lat ≈ 0.51 would be the true
    scale-INVARIANT substrate coupling.  In the actual substrate this
    identity holds only at the electron scale; the running K(ξ) breaks it
    elsewhere.  See :func:`alpha_M_naive` (Candidate A): mathematically

        α_lat_running(ξ) = σ_SI(ξ) × ξ² / (ℏc) ≡ α_M(ξ)

    so the "running" version of B is identically Candidate A.

    Returns:
        σ_lat = 0.51 (dimensionless, scale-invariant).
    """
    return SIGMA_LATTICE_NATURAL


def alpha_lat_running(xi_m: float) -> float:
    """Candidate B (running version): σ_SI × ξ² / (ℏc).

    Mathematically IDENTICAL to :func:`alpha_M_naive` because the units
    work out the same way (the "GeV² × GeV⁻²" of A is exactly "J·m / J·m"
    of B once you trace through the conversions).  Provided here to make
    the equality explicit in the comparison table and to confirm that
    there is no extra information in B beyond A.

    Args:
        xi_m: Length scale [m].

    Returns:
        Dimensionless σ_SI × ξ² / (ℏc) — equals α_M_naive(ξ) up to
        floating-point precision.
    """
    K_xi = K_at(xi_m)
    sigma_SI = SIGMA_LATTICE_NATURAL * K_xi
    return sigma_SI * xi_m ** 2 / HBARC_J_M


def alpha_inverse(xi_m: float) -> float:
    """Candidate C: α_C = ℏc / (K × ξ⁴) — action-quantum residual.

    The substrate identity ℏ = K ξ⁴ / c (action quantum, §18.46.1) means
    ℏc / (K ξ⁴) = 1 — but this is only exact at the electron scale.
    With running K(ξ) = K_e (ξ_e/ξ)^a:

        ℏc / (K(ξ) ξ⁴) = ℏc / (K_e ξ_e^a × ξ^(4-a))
                       ∝ ξ^(a-4)
                       ∝ Q^(4-a) = Q^(4-(-5.69)) = Q^(+9.69)

    α_C GROWS with Q — the OPPOSITE direction from asymptotic freedom.

    Args:
        xi_m: Length scale [m].

    Returns:
        Dimensionless ℏc / (K ξ⁴).
    """
    K_xi = K_at(xi_m)
    return HBARC_J_M / (K_xi * xi_m ** 4)


def alpha_inverse_running(xi_m: float) -> float:
    """Candidate D: power-law decrease forced through the QCD anchor.

    α_D(Q) = α_M(Q_QCD) × (Q_QCD/Q)^(2-a)
             where 2 - a > 0 for a ≈ -5.69, so α_D decreases as Q grows.

    Args:
        xi_m: Length scale [m].

    Returns:
        α_D [dimensionless].
    """
    a = _running_exponent_n1()
    Q_GeV = HBARC_GEVFM * 1e-15 / xi_m
    Q_QCD = HBARC_GEVFM * 1e-15 / _XI_QCD
    alpha_M_at_QCD = alpha_M_naive(_XI_QCD)
    return alpha_M_at_QCD * (Q_QCD / Q_GeV) ** (2.0 - a)


def alpha_QCD_log(Q_GeV: float) -> float:
    """Candidate E: QCD-style logarithmic running fit to substrate at two scales.

        1/α_E(Q) = 1/α_E(Q_0) + b × ln(Q/Q_0)

    We fit α_E to match α_M(Q) at TWO points (Q=1 GeV and Q=M_Z) — using the
    SUBSTRATE's α_M values as targets, not PDG values.  This isolates the
    question "is the substrate's running shape logarithmic?" from "does its
    magnitude match α_s?".

    Args:
        Q_GeV: Energy scale [GeV].

    Returns:
        α_E(Q) — the QCD-style log fit through substrate points.
    """
    Q0 = 1.0
    Q1 = 91.2
    alpha0 = alpha_M_naive(Q_to_xi_m(Q0))
    alpha1 = alpha_M_naive(Q_to_xi_m(Q1))
    if alpha0 <= 0.0 or alpha1 <= 0.0:
        return float("nan")
    b = (1.0 / alpha1 - 1.0 / alpha0) / math.log(Q1 / Q0)
    inv_alpha = 1.0 / alpha0 + b * math.log(Q_GeV / Q0)
    if inv_alpha <= 0.0:
        return float("inf")
    return 1.0 / inv_alpha


# ---------------------------------------------------------------------------
# Q ↔ ξ conversion
# ---------------------------------------------------------------------------

def xi_m_to_Q_GeV(xi_m: float) -> float:
    """Convert ξ [m] to Q [GeV] via Q = ℏc/ξ.

    Args:
        xi_m: Length [m].

    Returns:
        Q [GeV].
    """
    return HBARC_GEVFM * 1e-15 / xi_m


def Q_to_xi_m(Q_GeV: float) -> float:
    """Convert Q [GeV] to ξ [m] via ξ = ℏc/Q.

    Args:
        Q_GeV: Energy [GeV].

    Returns:
        ξ [m].
    """
    return HBARC_GEVFM * 1e-15 / Q_GeV


# ---------------------------------------------------------------------------
# QCD α_s reference (PDG 2024) at selected scales
# ---------------------------------------------------------------------------

ALPHA_S_PDG: Final[dict[float, float]] = {
    # Q [GeV]: α_s(Q) reference value (PDG 2024 / world-average)
    1.0:    0.45,
    1.78:   0.32,    # M_τ
    2.0:    0.30,
    5.0:    0.21,
    10.0:   0.18,
    30.0:   0.135,
    91.1876: 0.1180,  # M_Z (PDG world average 2024)
    100.0:  0.116,
    1000.0: 0.082,
}


# ---------------------------------------------------------------------------
# Comparison machinery
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CouplingPoint:
    """One Q-point with substrate predictions and the PDG reference.

    Attributes:
        xi_m:                Length scale [m].
        Q_GeV:               Equivalent energy scale [GeV] = ℏc/ξ.
        K_Pa:                Substrate stiffness K(ξ) [Pa] from the n=1 fit.
        sigma_GeV2:          Substrate string tension σ(ξ) [GeV²].
        alpha_M_naive:       Candidate A: σ × ξ² (dimensionless).
        alpha_lat_inv:       Candidate B: σ_lat = 0.51 (scale-invariant).
        alpha_inverse:       Candidate C: ℏc / (K ξ⁴).
        alpha_inverse_run:   Candidate D: power-law decrease normalised at QCD.
        alpha_QCD_log:       Candidate E: QCD-style log running fit to A.
        alpha_s_pdg:         PDG reference α_s(Q) (or NaN if not tabulated).
    """
    xi_m: float
    Q_GeV: float
    K_Pa: float
    sigma_GeV2: float
    alpha_M_naive: float
    alpha_lat_inv: float
    alpha_inverse: float
    alpha_inverse_run: float
    alpha_QCD_log: float
    alpha_s_pdg: float


def compute_at_xi(xi_m: float, alpha_s_ref: float = float("nan")) -> CouplingPoint:
    """Compute all candidate substrate couplings at a single ξ.

    Args:
        xi_m:        Length scale [m].
        alpha_s_ref: Reference α_s(Q) (PDG) at this scale; NaN if unknown.

    Returns:
        :class:`CouplingPoint` with all candidate couplings populated.
    """
    Q_GeV = xi_m_to_Q_GeV(xi_m)
    return CouplingPoint(
        xi_m=xi_m,
        Q_GeV=Q_GeV,
        K_Pa=K_at(xi_m),
        sigma_GeV2=sigma_GeV2_at(xi_m),
        alpha_M_naive=alpha_M_naive(xi_m),
        alpha_lat_inv=alpha_lattice_invariant(),
        alpha_inverse=alpha_inverse(xi_m),
        alpha_inverse_run=alpha_inverse_running(xi_m),
        alpha_QCD_log=alpha_QCD_log(Q_GeV),
        alpha_s_pdg=alpha_s_ref,
    )


def standard_scale_grid() -> list[float]:
    """Return a list of physically meaningful ξ scales [m] for comparison.

    Spans Q ∈ [0.2, 1000] GeV in roughly logarithmic steps, plus the
    actual ξ_QCD = 0.2 fm anchor where α_M = 0.185 by construction.

    Returns:
        Sorted list of ξ values [m] (largest ξ / smallest Q first).
    """
    Q_GeV_values = [
        0.2, 0.5, 0.9866,           # 0.9866 GeV ⇔ ξ = 0.2 fm exactly
        1.0, 1.78, 2.0, 5.0, 10.0,
        30.0, 91.1876, 100.0, 1000.0,
    ]
    return [Q_to_xi_m(Q) for Q in Q_GeV_values]


@dataclass
class AlphaRunningSummary:
    """Aggregate summary of substrate α_s-running test.

    Attributes:
        a_substrate:   Substrate K(ξ) running exponent (n=1 form, signed).
        n_triple:      n from the triple-anchor fit (alt RGE).
        a_triple:      a from the triple-anchor fit.
        points:        List of CouplingPoints across the standard grid.
        verdict:       Human-readable conclusion about which candidate (if
                       any) reproduces the QCD α_s running.
        summary_lines: The verdict as a list of lines (for programmatic use).
    """
    a_substrate: float
    n_triple: float
    a_triple: float
    points: list[CouplingPoint]
    verdict: str = ""
    summary_lines: list[str] = field(default_factory=list)


def run_full_comparison() -> AlphaRunningSummary:
    """Run the full α_M(Q) vs α_s(Q) comparison and produce a summary.

    Returns:
        :class:`AlphaRunningSummary` with all candidate evaluations and a
        verdict on whether any substrate combination reproduces α_s(Q).
    """
    a_n1 = _running_exponent_n1()
    a_triple, n_triple = _triple_fit_params()

    points: list[CouplingPoint] = []
    for xi_m in standard_scale_grid():
        Q_GeV = xi_m_to_Q_GeV(xi_m)
        ref = float("nan")
        for Q_key, val in ALPHA_S_PDG.items():
            # ±0.1% tolerance on the key match
            if abs(Q_key - Q_GeV) / Q_key < 1e-3:
                ref = val
                break
        points.append(compute_at_xi(xi_m, alpha_s_ref=ref))

    verdict_lines = _build_verdict(points, a_n1)

    return AlphaRunningSummary(
        a_substrate=a_n1,
        n_triple=n_triple,
        a_triple=a_triple,
        points=points,
        verdict="\n".join(verdict_lines),
        summary_lines=verdict_lines,
    )


def _pick(points: list[CouplingPoint], target_Q: float) -> CouplingPoint:
    """Return the CouplingPoint nearest to target_Q [GeV]."""
    return min(points, key=lambda p: abs(p.Q_GeV - target_Q))


def _build_verdict(points: list[CouplingPoint], a_substrate: float) -> list[str]:
    """Produce a structured human-readable verdict on the comparison.

    Args:
        points:        Computed CouplingPoints.
        a_substrate:   K(ξ) running exponent (n=1 form, signed).

    Returns:
        List of report lines.
    """
    lines: list[str] = []
    lines.append("VERDICT: substrate K(ξ) running vs QCD α_s(Q) running")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Substrate running exponent (n=1 fit):  a = {a_substrate:+.4f}")
    lines.append("  K(ξ) = K_e × (ξ_e/ξ)^a, so K decreases as ξ shrinks below ξ_e.")
    lines.append("  σ(ξ) = σ_lat × K(ξ), so σ also decreases below ξ_e.")
    lines.append("")

    qcd_anchor = _pick(points, 0.9866)   # closest to ξ = 0.2 fm
    one_GeV = _pick(points, 1.0)
    mz = _pick(points, 91.1876)
    high = points[-1]
    low = points[0]

    # Candidate A: σ ξ²
    lines.append("Candidate A:  α_M = σ × ξ²  [§18.61.1 Möbius coupling]")
    lines.append(f"  α_M(Q≈0.987 GeV ⇔ ξ_QCD)  = {qcd_anchor.alpha_M_naive:.4f}     "
                 f"(spec target ≈ 0.185 ✓)")
    lines.append(f"  α_M(Q=1 GeV)              = {one_GeV.alpha_M_naive:.4e}")
    lines.append(f"  α_M(M_Z = 91.2 GeV)       = {mz.alpha_M_naive:.4e}")
    lines.append(f"  α_M(Q=1 TeV)              = {high.alpha_M_naive:.4e}")
    direction = "DECREASES" if mz.alpha_M_naive < qcd_anchor.alpha_M_naive else "INCREASES"
    arrow = "matches QCD direction" if direction == "DECREASES" else "OPPOSITE of QCD"
    lines.append(f"  Sign of running: {direction} with Q  ({arrow})")
    lines.append(f"  Power: α_M ∝ ξ^({2.0 - a_substrate:+.3f}) = "
                 f"Q^({a_substrate - 2.0:+.3f}) — POWER LAW, not log")
    lines.append("")

    # Candidate B: scale-invariant σ_lat = 0.51
    lines.append("Candidate B:  α_lat ≡ σ_lat = 0.51  [idealised, scale-INVARIANT]")
    lines.append("  σ_lat is the dimensionless lattice string tension; if the")
    lines.append("  substrate identity ℏ = K ξ⁴/c held at every scale, σ × ξ²/(ℏc)")
    lines.append("  would equal 0.51 always.  Trivially constant ⇒ cannot reproduce")
    lines.append("  α_s(Q) running.  (The 'running σ × ξ²/ℏc' is mathematically")
    lines.append("  identical to Candidate A — same value, same direction, same")
    lines.append("  failure mode at high Q.)")
    lines.append("")

    # Candidate C: ℏc / (K ξ⁴)
    lines.append("Candidate C:  α_C = ℏc / (K × ξ⁴)  [action-quantum residual]")
    lines.append(f"  α_C(Q≈0.987 GeV ⇔ ξ_QCD) = {qcd_anchor.alpha_inverse:.4e}")
    lines.append(f"  α_C(M_Z)                 = {mz.alpha_inverse:.4e}")
    lines.append(f"  α_C(Q=1 TeV)             = {high.alpha_inverse:.4e}")
    direction_c = ("DECREASES"
                   if mz.alpha_inverse < qcd_anchor.alpha_inverse
                   else "INCREASES")
    arrow_c = ("matches QCD direction"
               if direction_c == "DECREASES"
               else "OPPOSITE of QCD")
    lines.append(f"  Sign of running: {direction_c} with Q ({arrow_c})")
    lines.append(f"  Power: α_C ∝ ξ^({a_substrate - 4:+.3f}) = "
                 f"Q^({4 - a_substrate:+.3f})")
    lines.append("")

    # Candidate D: power-law decrease normalised at QCD
    lines.append("Candidate D:  α_D = α_M(QCD) × (Q_QCD/Q)^(2-a)  [forced power-law decrease]")
    lines.append(f"  α_D(Q≈0.987 GeV ⇔ ξ_QCD) = {qcd_anchor.alpha_inverse_run:.4e}")
    lines.append(f"  α_D(M_Z)                 = {mz.alpha_inverse_run:.4e}")
    lines.append(f"  α_D(Q=1 TeV)             = {high.alpha_inverse_run:.4e}")
    if not math.isnan(mz.alpha_s_pdg):
        ratio = mz.alpha_inverse_run / mz.alpha_s_pdg
        lines.append(f"  α_D(M_Z) / α_s_PDG(M_Z) = {ratio:.4e}")
    lines.append("  This is a POWER LAW, not the LOGARITHMIC running of QCD.")
    lines.append("")

    # Candidate E: QCD-style log running, fit to A at Q=1 GeV and Q=M_Z
    lines.append("Candidate E:  1/α_E(Q) = 1/α_E(Q_0) + b ln(Q/Q_0), fit to A at 1 GeV and M_Z")
    lines.append(f"  α_E(Q=1 GeV)   = {one_GeV.alpha_QCD_log:.4e}")
    lines.append(f"  α_E(M_Z)       = {mz.alpha_QCD_log:.4e}")
    lines.append(f"  α_E(Q=1 TeV)   = {high.alpha_QCD_log:.4e}")
    lines.append("  This ansatz matches A's two values BY CONSTRUCTION; the fact that")
    lines.append("  it fits A's intermediate points poorly demonstrates that A is NOT")
    lines.append("  logarithmic in Q — its true shape is a steep power law.")
    lines.append("")

    # Summary verdict
    lines.append("-" * 70)
    lines.append("CONCLUSION")
    lines.append("-" * 70)
    lines.append("")
    lines.append("The substrate K(ξ) running with a ≈ -5.69 is a POWER LAW in ξ.")
    lines.append("It does NOT reproduce the LOGARITHMIC QCD running of α_s(Q).")
    lines.append("")
    lines.append("Detailed comparison at standard scales:")
    lines.append("  At Q ≈ 1 GeV (ξ near QCD anchor):")
    lines.append(f"    α_M = {qcd_anchor.alpha_M_naive:.4f}  vs PDG α_s(1 GeV) ≈ 0.45  "
                 f"({qcd_anchor.alpha_M_naive / 0.45:.2%})")
    lines.append("  At Q = M_Z = 91.2 GeV:")
    if not math.isnan(mz.alpha_s_pdg):
        ratio_A_mz = mz.alpha_M_naive / mz.alpha_s_pdg
        lines.append(f"    α_M = {mz.alpha_M_naive:.3e}  vs PDG α_s(M_Z) = "
                     f"{mz.alpha_s_pdg:.4f}  (ratio {ratio_A_mz:.2e})")
    lines.append("")
    lines.append("Sign analysis of α_M = σ × ξ² (the §18.61.1 coupling):")
    lines.append(f"  σ ∝ K ∝ ξ^(-a) with a < 0 ⇒ σ ∝ ξ^(|a|) GROWS with ξ.")
    lines.append(f"  So as Q increases (ξ decreases), σ shrinks fast and α_M = σ×ξ²")
    lines.append(f"  shrinks even faster (ξ^(2+|a|) ≈ ξ^7.69 ∝ Q^(-7.69)).")
    lines.append("  Direction matches QCD asymptotic freedom, but power law is far")
    lines.append("  steeper than the log running — at M_Z it is below 10⁻¹⁵ instead")
    lines.append("  of 0.118.")
    lines.append("")
    lines.append("HONEST RESULT")
    lines.append("─" * 14)
    lines.append("  No simple combination of substrate primitives K(ξ), σ(ξ), ξ runs")
    lines.append("  with the LOGARITHMIC structure of the QCD β-function:")
    lines.append("    α_M = σ ξ²       — power-law DECREASE, far too steep.")
    lines.append("    α_lat = σ_SI ξ²/ℏc — same power-law as α_M (same sign of slope).")
    lines.append("    α_C = ℏc/(K ξ⁴)  — power-law INCREASE, wrong direction.")
    lines.append("    α_D, α_E         — constrained fits that demonstrate, by their")
    lines.append("                       failure to interpolate, that K-running is")
    lines.append("                       not logarithmic.")
    lines.append("")
    lines.append("  The §18.61.1 result α_M(QCD) ≈ 0.185 ≈ α_s(QCD) is therefore a")
    lines.append("  COINCIDENCE OF MAGNITUDE at one scale.  It supports the chromo-")
    lines.append("  magnetic N-Δ calculation at that scale, but it does NOT constitute")
    lines.append("  a derivation of QCD running from substrate physics.")
    lines.append("")
    lines.append("  To produce QCD-like running in the substrate, one would need a")
    lines.append("  separate mechanism (e.g. dynamical generation of ln Q in K-flow")
    lines.append("  via vacuum-loop screening of σ, or an effective α emerging from")
    lines.append("  Möbius-flux self-interactions) that is not present in the current")
    lines.append("  power-law K(ξ).")
    return lines


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

def format_table(points: list[CouplingPoint]) -> str:
    """Format a comparison table of all candidates vs PDG.

    Args:
        points: Sorted list of CouplingPoints (low Q to high Q).

    Returns:
        Multi-line string showing Q, all candidates, and the PDG reference.
    """
    lines: list[str] = []
    lines.append("=" * 116)
    lines.append(
        f"{'Q [GeV]':>10} {'σ [GeV²]':>12} {'α_M=σξ² (A)':>14} "
        f"{'σ_lat (B)':>10} {'α_C=ℏc/Kξ⁴ (C)':>16} "
        f"{'α_D pow':>12} {'α_E log':>12} {'α_s PDG':>10} {'A/PDG':>10}"
    )
    lines.append("=" * 116)
    for p in points:
        ratio_str = "—"
        if not math.isnan(p.alpha_s_pdg) and p.alpha_s_pdg > 0:
            ratio_str = f"{p.alpha_M_naive / p.alpha_s_pdg:.2e}"
        pdg_str = f"{p.alpha_s_pdg:.4f}" if not math.isnan(p.alpha_s_pdg) else "—"
        lines.append(
            f"{p.Q_GeV:>10.4f} {p.sigma_GeV2:>12.4e} {p.alpha_M_naive:>14.4e} "
            f"{p.alpha_lat_inv:>10.4f} {p.alpha_inverse:>16.4e} "
            f"{p.alpha_inverse_run:>12.4e} {p.alpha_QCD_log:>12.4e} "
            f"{pdg_str:>10} {ratio_str:>10}"
        )
    lines.append("=" * 116)
    return "\n".join(lines)


def print_full_report() -> None:
    """Print the full comparison report (table + verdict)."""
    summary = run_full_comparison()
    print("=" * 78)
    print(" SUBSTRATE α_M(Q) vs QCD α_s(Q) — running comparison")
    print("=" * 78)
    print()
    print(f"  Substrate K-running exponent (n=1):  a = {summary.a_substrate:+.4f}")
    print(f"  Triple-anchor fit:                    n = {summary.n_triple:.4f}, "
          f"a = {summary.a_triple:+.4e}")
    print()
    print(format_table(summary.points))
    print()
    print(summary.verdict)


if __name__ == "__main__":
    print_full_report()
