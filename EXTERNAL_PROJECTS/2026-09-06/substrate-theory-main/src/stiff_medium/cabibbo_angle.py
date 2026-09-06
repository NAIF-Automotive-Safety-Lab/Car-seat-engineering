"""Cabibbo angle θ_C from substrate topology — multi-hypothesis test.

Spec context (§§18.10, 18.30, 18.49):
- §18.10: Möbius half-flux U(1) bundle, A = (1/2) dθ, holonomy = -1
- §18.30: vertex closure → exactly 3 generations
- §18.49: SU(3) extension; CKM matrix has 4 free parameters (3 angles + 1 phase)
- §18.49.7: "CKM mixing comes from rotations between different stress-loaded
            vertex states. The Möbius bundle topology constrains which mixings
            are allowed. Specific values are empirical input."

The Cabibbo angle parameterizes V_us in the CKM matrix:
    V_us = sin(θ_C) cos(θ_13) ≈ sin(θ_C) ≈ 0.2255
    θ_C ≈ 13.04° ≈ 0.2275 rad

Empirical relation (the "GIM identity"):
    sin(θ_C) ≈ √(m_d/m_s)   with  m_d ≈ 4.7 MeV, m_s ≈ 95 MeV
    → √(0.0495) = 0.2222

GOAL OF THIS MODULE:
    Test multiple substrate-topological hypotheses for θ_C with NO hand-tuning.
    Catalogue every reasonable candidate; report numerical errors; identify
    the cleanest substrate match (if any).

CANDIDATE FAMILIES:
  (A) GIM/dimensional: sin(θ_C) = √(m_d/m_s).  Reduces θ_C to a single quark-
      mass ratio.  In the substrate this becomes a Yukawa-ratio question
      (still empirical for individual Yukawas).
  (B) Topological angles: π/N for integer N near the empirical value of
      ~13.8 (since 180/13.04 ≈ 13.80).  Variants: 2π/N, π/(N·M), etc.
  (C) Z₃ vertex / Möbius half-flux combinations: (1/2)π / k, (2π/3)/k, etc.
  (D) Foot-phase residual: the empirical Foot phase is δ_emp ≈ 1.404π.
      Differences (δ - 7π/5), (1.5π - δ), and Foot-related angles are tested.
  (E) Substrate-scale ratios: ξ_QCD / λ_strange and similar dimensional
      constructions.
  (F) Geometric ratios from the 45° cone and the 4π fermion period.

CRITICAL RULES (from the task):
- DO NOT hand-tune the angle.
- Report ALL reasonable substrate candidates with errors.
- The Cabibbo angle is an EXTRA dimensionless number beyond what we have
  derived (σ, ξ_QCD set hadron masses; θ_C is a separate flavor-mixing
  parameter).

OUTPUT INTERPRETATION:
- If a clean substrate-derived angle gives 13.04 ± 0.3° → MAJOR result.
- If close but not exact → identify the missing structure.
- If no clean substrate match → report as open (matches §18.49.7 status).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import math
import numpy as np


# ---------------------------------------------------------------------------
# Empirical landmarks (PDG 2024)
# ---------------------------------------------------------------------------

# Cabibbo angle (PDG values)
SIN_THETA_C: float = 0.2255       # sin(θ_C) = |V_us|/cos(θ_13) ≈ |V_us|
THETA_C_RAD: float = 0.2275       # θ_C in radians (≈ asin(0.2255))
THETA_C_DEG: float = 13.04        # θ_C in degrees
TAN_THETA_C: float = 0.2316       # tan(θ_C) ≈ √(m_d/m_s)

# Light quark masses (PDG MS-bar at 2 GeV)
M_U_MEV: float = 2.16             # up quark
M_D_MEV: float = 4.67             # down quark
M_S_MEV: float = 93.4             # strange quark
M_C_MEV: float = 1.27e3           # charm quark
M_B_MEV: float = 4.18e3           # bottom quark
M_T_MEV: float = 172.7e3          # top quark

# Substrate primitives (from spec / src/)
XI_QCD_FM: float = 0.2            # QCD coherence length [fm]
SIGMA_QCD_GEV2: float = 0.180     # string tension at ξ_QCD [GeV²]

# Hadron mass anchors (PDG)
M_PI_MEV: float = 139.57          # charged pion
M_K_MEV: float = 493.68           # charged kaon
M_PROTON_MEV: float = 938.272

# Foot phase (from foot_phase_derivation.py): δ_emp ≈ 1.404π
DELTA_FOOT_PI: float = 1.404
DELTA_FOOT_RAD: float = DELTA_FOOT_PI * math.pi


# ---------------------------------------------------------------------------
# Candidate dataclass
# ---------------------------------------------------------------------------


@dataclass
class CabibboCandidate:
    """A candidate substrate-derived prediction for the Cabibbo angle.

    Attributes:
        name: Human-readable label.
        family: Strategy family (A=GIM, B=topo π/N, C=Z₃/Möbius, D=Foot, ...).
        formula: Symbolic / verbal formula.
        sin_theta: Predicted sin(θ_C).
        theta_deg: Predicted θ_C [degrees].
        theta_rad: Predicted θ_C [radians].
        substrate_basis: Connection to spec ingredients.
        pct_error_sin: % error in sin(θ_C).
        pct_error_deg: % error in θ_C[deg].
        verdict: Honest assessment.
    """

    name: str
    family: str
    formula: str
    sin_theta: float
    theta_deg: float
    theta_rad: float
    substrate_basis: str
    pct_error_sin: float
    pct_error_deg: float
    verdict: str


def _make_candidate(
    name: str,
    family: str,
    formula: str,
    sin_theta: float = float("nan"),
    theta_rad: float | None = None,
    theta_deg: float | None = None,
    substrate_basis: str = "",
    verdict: str = "",
) -> CabibboCandidate:
    """Build a candidate from any of (sin_theta, theta_rad, theta_deg).

    Exactly one of the three values defines the angle; the others are
    derived consistently.
    """
    if theta_deg is not None:
        deg = float(theta_deg)
        rad = math.radians(deg)
        s = math.sin(rad)
    elif theta_rad is not None:
        rad = float(theta_rad)
        deg = math.degrees(rad)
        s = math.sin(rad)
    elif not math.isnan(sin_theta):
        s = float(sin_theta)
        # asin handles |s| <= 1; otherwise NaN
        if abs(s) > 1.0:
            rad = float("nan")
            deg = float("nan")
        else:
            rad = math.asin(s)
            deg = math.degrees(rad)
    else:
        raise ValueError("must specify one of sin_theta, theta_rad, theta_deg")

    pct_s = abs(s - SIN_THETA_C) / SIN_THETA_C * 100.0 if math.isfinite(s) else float("inf")
    pct_d = (
        abs(deg - THETA_C_DEG) / THETA_C_DEG * 100.0
        if math.isfinite(deg) else float("inf")
    )

    return CabibboCandidate(
        name=name,
        family=family,
        formula=formula,
        sin_theta=float(s),
        theta_deg=float(deg),
        theta_rad=float(rad),
        substrate_basis=substrate_basis,
        pct_error_sin=float(pct_s),
        pct_error_deg=float(pct_d),
        verdict=verdict,
    )


# ---------------------------------------------------------------------------
# Family A: GIM identity (mass-ratio based)
# ---------------------------------------------------------------------------


def gim_identity_pdg() -> CabibboCandidate:
    """sin(θ_C) = √(m_d/m_s) using PDG quark masses.

    The "GIM identity" is the empirical observation that, for the SM CKM
    structure, sin(θ_C) is well approximated by the geometric mean of the
    down/strange quark mass ratio.

    With PDG m_d ≈ 4.67 MeV, m_s ≈ 93.4 MeV:
        √(m_d/m_s) = √(0.0500) = 0.2236   (vs 0.2255 empirical)
    """
    s = math.sqrt(M_D_MEV / M_S_MEV)
    return _make_candidate(
        name="GIM identity (PDG quark masses)",
        family="A_GIM",
        formula="sin(θ_C) = √(m_d/m_s)",
        sin_theta=s,
        substrate_basis=(
            "Reduces θ_C to the d/s quark Yukawa-mass ratio.  In our model "
            "the quark Yukawas are still 6 free parameters per §18.49.6 "
            "— the GIM identity is therefore not yet a substrate-first "
            "PREDICTION, only an empirical reduction."
        ),
        verdict=(
            "Excellent numerical match (~1% error) but not derived: just "
            "shifts the question from θ_C to m_d/m_s."
        ),
    )


def gim_identity_with_charm_correction() -> CabibboCandidate:
    """sin(θ_C) = √(m_d/m_s) × (1 - m_d/m_s)  (subleading correction).

    In the diagonalisation of a 2×2 mass matrix, the rotation angle
    satisfies tan(2θ) = 2|M_12| / (M_22 - M_11) and for a "democratic"
    ansatz one recovers tan(θ) = √(m_d/m_s) at leading order.  The next
    correction is multiplicative ∼ (1 - m_d/m_s).
    """
    r = M_D_MEV / M_S_MEV
    s_lead = math.sqrt(r)
    s = s_lead * (1.0 - r)
    return _make_candidate(
        name="GIM with subleading mass-matrix correction",
        family="A_GIM",
        formula="sin(θ_C) = √(m_d/m_s) × (1 - m_d/m_s)",
        sin_theta=s,
        substrate_basis="2×2 mass-matrix diagonalisation, NLO in m_d/m_s.",
        verdict="Marginal change vs LO; subleading correction reduces match.",
    )


def gim_with_tan() -> CabibboCandidate:
    """tan(θ_C) = √(m_d/m_s)  (alternative form).

    Some textbooks give tan(θ_C) rather than sin(θ_C) ≈ √(m_d/m_s).
    For small angles the two coincide.
    """
    t = math.sqrt(M_D_MEV / M_S_MEV)
    rad = math.atan(t)
    return _make_candidate(
        name="GIM identity using tan(θ_C)",
        family="A_GIM",
        formula="tan(θ_C) = √(m_d/m_s)",
        theta_rad=rad,
        substrate_basis="Alt form of GIM (tan rather than sin); same content small-angle.",
        verdict="Tiny difference for ~13°; same status as sin(θ_C) form.",
    )


# ---------------------------------------------------------------------------
# Family B: Pure topological angles π/N and rational multiples
# ---------------------------------------------------------------------------


def topological_pi_over_N(N: int) -> CabibboCandidate:
    """θ_C = π/N for integer N, in radians.

    The empirical θ_C ≈ 0.2275 rad ⇒ N = π/0.2275 ≈ 13.80.
    So N = 14 gives π/14 ≈ 12.857° (close).
    N = 13 gives 13.846° (also close).
    Neither matches exactly — but report all to be thorough.
    """
    rad = math.pi / N
    return _make_candidate(
        name=f"θ_C = π/{N}",
        family="B_topo_piN",
        formula=f"θ_C = π/{N}",
        theta_rad=rad,
        substrate_basis=f"N-fold partitioning of half-flux holonomy π.",
        verdict="",  # filled by pct_error
    )


def topological_2pi_over_N(N: int) -> CabibboCandidate:
    """θ_C = 2π/N for integer N.  N = 27 (= 3³) is special."""
    rad = 2.0 * math.pi / N
    return _make_candidate(
        name=f"θ_C = 2π/{N}",
        family="B_topo_2piN",
        formula=f"θ_C = 2π/{N}",
        theta_rad=rad,
        substrate_basis=(
            f"N-fold partitioning of full holonomy 2π.  N=27 = 3³ from "
            "vertex closure cubed (3 generations × 3 colours × 3 spatial)."
            if N == 27 else f"N-fold partitioning of full holonomy 2π."
        ),
        verdict="",
    )


def topological_pi_over_NM(N: int, M: int) -> CabibboCandidate:
    """θ_C = π/(N·M) for integers N, M (composite denominators)."""
    rad = math.pi / (N * M)
    return _make_candidate(
        name=f"θ_C = π/({N}·{M}) = π/{N*M}",
        family="B_topo_composite",
        formula=f"θ_C = π/({N}·{M})",
        theta_rad=rad,
        substrate_basis=f"Composite denominator: {N} × {M}.",
        verdict="",
    )


# ---------------------------------------------------------------------------
# Family C: Möbius half-flux × Z₃ vertex closure combinations
# ---------------------------------------------------------------------------


def half_flux_over_3sq() -> CabibboCandidate:
    """θ_C = (π/2) / 3² = π/18  (half-flux divided by 3² = 9 from §18.30 cubed).

    Möbius half-flux holonomy = π.  Half of that = π/2 (Chern × full).
    Divide by 3² (vertex-closure squared = 9 stress states).
    """
    rad = math.pi / 18.0
    return _make_candidate(
        name="θ_C = π/18 (half-flux/9)",
        family="C_mobius_z3",
        formula="θ_C = (π/2) / 9 = π/18",
        theta_rad=rad,
        substrate_basis=(
            "Half-flux π/2 divided by 9 = 3² stress states in vertex "
            "closure; per §18.30 each generation has its own Z₃ phase."
        ),
        verdict="",
    )


def vertex_closure_over_N(N: int) -> CabibboCandidate:
    """θ_C = (2π/3) / N  — Z₃ sector subdivided."""
    rad = (2.0 * math.pi / 3.0) / N
    return _make_candidate(
        name=f"θ_C = (2π/3)/{N}",
        family="C_mobius_z3",
        formula=f"θ_C = (2π/3) / {N}",
        theta_rad=rad,
        substrate_basis=f"Z₃ sector angle 2π/3 subdivided by {N}.",
        verdict="",
    )


def chern_winding_combinations() -> list[CabibboCandidate]:
    """Various combinations of c₁ (Chern) and winding numbers."""
    cands = []
    # c1 × winding × 1/N for various N
    for N in [12, 13, 14, 15, 16]:
        rad = math.pi / (2 * N)  # = (1/2 chern) × (π) × (1/N)
        cands.append(_make_candidate(
            name=f"θ_C = c₁·π/{N} = π/{2*N}",
            family="C_mobius_z3",
            formula=f"θ_C = c₁ × π / {N} = π/{2*N}",
            theta_rad=rad,
            substrate_basis=f"Möbius Chern (1/2) × half-flux (π) / {N}.",
            verdict="",
        ))
    return cands


# ---------------------------------------------------------------------------
# Family D: Foot-phase residual
# ---------------------------------------------------------------------------


def foot_phase_residual() -> list[CabibboCandidate]:
    """Try angles related to the Foot phase δ_emp ≈ 1.404π.

    Hypothesis: the Cabibbo angle is the residual mismatch between the Foot
    phase and the nearest "topologically clean" rational multiple of π.

    Candidates:
    - δ_emp - 7π/5 ≈ 0.004π × π = 0.013 rad ≈ 0.72°  (way too small)
    - δ_emp - π = 0.404π = 1.270 rad ≈ 72.7°  (way too large)
    - π/2 - δ_emp/something
    - The angle of the Möbius "tilt" between adjacent Foot generations
    """
    cands = []

    # D1: δ - 7π/5 (residual to nearest rational)
    delta_rad = DELTA_FOOT_RAD
    res = delta_rad - 1.4 * math.pi
    cands.append(_make_candidate(
        name="θ_C = δ_emp − 7π/5  (Foot-rational residual)",
        family="D_foot",
        formula="θ_C = δ_emp − 7π/5",
        theta_rad=abs(res),
        substrate_basis="Foot phase residual to nearest simple rational 7π/5.",
        verdict="",
    ))

    # D2: π/3 - δ_emp/something  (try generation-phase difference)
    # Actually try: angle between adjacent generations in the Foot
    # parametrization at the empirical δ.  The angles are 2πn/3 + δ for
    # n=0,1,2; the "tilt" between adjacent ones is just 2π/3 (geometric).
    cands.append(_make_candidate(
        name="θ_C = (3π/2 − δ_emp)/2",
        family="D_foot",
        formula="θ_C = (3π/2 − δ_emp)/2",
        theta_rad=(3.0 * math.pi / 2.0 - delta_rad) / 2.0,
        substrate_basis=(
            "Half-residual to 3π/2 (full Foot inversion).  No deep "
            "topological motivation."
        ),
        verdict="",
    ))

    # D3: empirical δ_emp − 1.4π in degrees
    res_deg = abs(delta_rad - 1.4 * math.pi) * 180.0 / math.pi
    # If we interpret this as DEGREES rather than radians
    cands.append(_make_candidate(
        name="θ_C [deg] = (δ_emp − 7π/5) × 180/π × scale",
        family="D_foot",
        formula="ratio × Foot residual",
        theta_deg=res_deg * 1000.0,  # exaggerated check
        substrate_basis="Foot residual interpreted in degrees (way too small).",
        verdict="",
    ))

    return cands


# ---------------------------------------------------------------------------
# Family E: Substrate-scale dimensional ratios
# ---------------------------------------------------------------------------


def substrate_scale_ratios() -> list[CabibboCandidate]:
    """Dimensional ratios from substrate scales.

    Candidates:
    - sin(θ_C) = ξ_QCD / λ_strange where λ_strange = ℏc/m_s
    - sin(θ_C) = m_pi / m_proton  (pion-to-proton mass ratio)
    - sin(θ_C) = (m_pi / m_K)²
    """
    cands = []

    # E1: sin(θ_C) = ξ_QCD / λ_strange
    # λ_strange = ℏc / m_s = 197 MeV·fm / 93.4 MeV = 2.110 fm
    HBARC_MEV_FM = 197.327
    lam_strange_fm = HBARC_MEV_FM / M_S_MEV
    s = XI_QCD_FM / lam_strange_fm
    cands.append(_make_candidate(
        name="sin(θ_C) = ξ_QCD / λ_strange",
        family="E_dim_substrate",
        formula="sin(θ_C) = ξ_QCD × m_s / (ℏc)",
        sin_theta=s,
        substrate_basis="Dimensional ratio of QCD coherence length to strange Compton wavelength.",
        verdict="",
    ))

    # E2: sin(θ_C) = m_d / something with substrate dim
    # m_d / m_pi = 4.67 / 139.6 = 0.0334 (too small)
    s_dpi = M_D_MEV / M_PI_MEV
    cands.append(_make_candidate(
        name="sin(θ_C) = m_d / m_π",
        family="E_dim_substrate",
        formula="sin(θ_C) = m_d / m_π",
        sin_theta=s_dpi,
        substrate_basis="Down-to-pion mass ratio.",
        verdict="",
    ))

    # E3: sin(θ_C) = m_π / m_p
    s_pip = M_PI_MEV / M_PROTON_MEV
    cands.append(_make_candidate(
        name="sin(θ_C) = m_π / m_p",
        family="E_dim_substrate",
        formula="sin(θ_C) = m_π / m_proton",
        sin_theta=s_pip,
        substrate_basis="Pion-to-proton ratio (Goldstone vs nucleon scales).",
        verdict="",
    ))

    # E4: sin(θ_C) = (m_π / m_K)²
    s_pik = (M_PI_MEV / M_K_MEV) ** 2
    cands.append(_make_candidate(
        name="sin(θ_C) = (m_π/m_K)²",
        family="E_dim_substrate",
        formula="sin(θ_C) = (m_π/m_K)²",
        sin_theta=s_pik,
        substrate_basis="Squared pion/kaon ratio (light/strange Goldstone hierarchy).",
        verdict="",
    ))

    # E5: sin(θ_C) = m_π / m_K  (light/strange Goldstone)
    s_lk = M_PI_MEV / M_K_MEV
    cands.append(_make_candidate(
        name="sin(θ_C) = m_π / m_K",
        family="E_dim_substrate",
        formula="sin(θ_C) = m_π / m_K",
        sin_theta=s_lk,
        substrate_basis="Direct π-to-K mass ratio.",
        verdict="",
    ))

    return cands


# ---------------------------------------------------------------------------
# Family F: Geometric / 45° cone constructions
# ---------------------------------------------------------------------------


def geometric_45_cone() -> list[CabibboCandidate]:
    """Geometric constructions from the 45° cone (§18.3) and 4π fermion period."""
    cands = []

    # F1: 45° / 3 (Z₃ sector partitioning of cone half-angle)
    cands.append(_make_candidate(
        name="θ_C = 45°/3 = 15°",
        family="F_geom",
        formula="θ_C = (45°)/3 = 15°",
        theta_deg=15.0,
        substrate_basis="Cone half-angle 45° divided into 3 generation sectors.",
        verdict="",
    ))

    # F2: 4π / N for N = 100 (~13.0°)... actually 4π/N works for generation period
    # 4π in degrees = 720; 720/55 ≈ 13.09°  (N=55 is non-special)
    cands.append(_make_candidate(
        name="θ_C = 4π/55 (fermion period / 55)",
        family="F_geom",
        formula="θ_C = 720°/55",
        theta_deg=720.0 / 55.0,
        substrate_basis=(
            "4π fermion period (Möbius double-cover) divided by 55 — "
            "no obvious substrate origin for 55; reported for completeness."
        ),
        verdict="",
    ))

    # F3: 360° / 3³ = 360°/27 = 13.333°
    cands.append(_make_candidate(
        name="θ_C = 360°/27 (full circle / 3³)",
        family="F_geom",
        formula="θ_C = 360°/27 = 360°/3³",
        theta_deg=360.0 / 27.0,
        substrate_basis=(
            "Full azimuth 360° divided by 3³ = 27 (vertex-closure cubed: "
            "3 generations × 3 colours × 3 spatial)."
        ),
        verdict="",
    ))

    # F4: 360° / 28 = 12.857° = π/14
    cands.append(_make_candidate(
        name="θ_C = 360°/28 = 180°/14",
        family="F_geom",
        formula="θ_C = 360°/28",
        theta_deg=360.0 / 28.0,
        substrate_basis="Full azimuth divided by 28 = 4×7 — no clean substrate basis.",
        verdict="",
    ))

    # F5: arctan(1/2 × √(m_d/m_s)) — small-angle Möbius half-shift
    rad = math.atan(0.5 * math.sqrt(M_D_MEV / M_S_MEV))
    cands.append(_make_candidate(
        name="θ_C = arctan(½ √(m_d/m_s)) — half-flux GIM",
        family="F_geom",
        formula="θ_C = arctan(c₁ × √(m_d/m_s)),  c₁=1/2",
        theta_rad=rad,
        substrate_basis="Möbius Chern factor 1/2 multiplied into the GIM identity.",
        verdict="",
    ))

    # F6: arcsin(2 × √(m_d/m_s) × sin(π/14))
    # purely combinatorial: tests whether topology + GIM combine
    rad = math.asin(2.0 * math.sqrt(M_D_MEV / M_S_MEV) * math.sin(math.pi / 14.0))
    cands.append(_make_candidate(
        name="θ_C = arcsin(2 × √(m_d/m_s) × sin(π/14))",
        family="F_geom",
        formula="topo·GIM combination",
        theta_rad=rad,
        substrate_basis="Combinatorial test (no first-principles motivation).",
        verdict="",
    ))

    return cands


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------


def all_candidates() -> list[CabibboCandidate]:
    """Build the complete list of substrate candidate predictions.

    Returns:
        List of CabibboCandidate, sorted by absolute pct_error_deg ascending.
    """
    cands: list[CabibboCandidate] = []

    # Family A: GIM
    cands.append(gim_identity_pdg())
    cands.append(gim_identity_with_charm_correction())
    cands.append(gim_with_tan())

    # Family B: π/N
    for N in [10, 11, 12, 13, 14, 15, 16, 17, 18]:
        cands.append(topological_pi_over_N(N))

    # Family B: 2π/N (special: N=27)
    for N in [24, 25, 26, 27, 28, 29, 30]:
        cands.append(topological_2pi_over_N(N))

    # Family B: composite π/(N·M)
    for N, M in [(2, 7), (3, 5), (3, 4), (4, 4), (2, 8)]:
        cands.append(topological_pi_over_NM(N, M))

    # Family C: Möbius × Z₃
    cands.append(half_flux_over_3sq())
    for N in [4, 5, 6, 7]:
        cands.append(vertex_closure_over_N(N))
    cands.extend(chern_winding_combinations())

    # Family D: Foot residuals
    cands.extend(foot_phase_residual())

    # Family E: substrate-scale ratios
    cands.extend(substrate_scale_ratios())

    # Family F: geometric / 45° cone
    cands.extend(geometric_45_cone())

    # Sort by absolute % error in θ_C (degrees)
    cands.sort(key=lambda c: c.pct_error_deg)
    return cands


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def report_candidates(cands: Sequence[CabibboCandidate],
                      max_rows: int = 999) -> str:
    """Format candidates as a textual report.

    Args:
        cands: Iterable of CabibboCandidate.
        max_rows: Limit number of rows printed.

    Returns:
        Multi-line string.
    """
    lines = []
    lines.append("=" * 95)
    lines.append("CABIBBO ANGLE CANDIDATES — substrate hypothesis test")
    lines.append("=" * 95)
    lines.append(f"Empirical: sin θ_C = {SIN_THETA_C:.4f},  θ_C = {THETA_C_DEG:.2f}°,  "
                 f"θ_C = {THETA_C_RAD:.4f} rad")
    lines.append("-" * 95)
    head = (
        f"{'#':>3} {'family':<14} {'name':<43} {'sin θ':>7} {'θ[deg]':>7} "
        f"{'%err sin':>8} {'%err deg':>8}"
    )
    lines.append(head)
    lines.append("-" * 95)
    for i, c in enumerate(cands[:max_rows], start=1):
        s = c.sin_theta
        d = c.theta_deg
        s_str = f"{s:.4f}" if math.isfinite(s) else "  NaN "
        d_str = f"{d:.2f}" if math.isfinite(d) else "  NaN "
        lines.append(
            f"{i:>3} {c.family:<14} {c.name[:43]:<43} {s_str:>7} {d_str:>7} "
            f"{c.pct_error_sin:>7.2f}% {c.pct_error_deg:>7.2f}%"
        )
    lines.append("-" * 95)
    return "\n".join(lines)


def best_n(cands: Sequence[CabibboCandidate], n: int = 5) -> list[CabibboCandidate]:
    """Top-n candidates by % error in degrees."""
    return sorted(cands, key=lambda c: c.pct_error_deg)[:n]


# ---------------------------------------------------------------------------
# Honest assessment
# ---------------------------------------------------------------------------


def honest_assessment(cands: Sequence[CabibboCandidate]) -> dict:
    """Apply the result-classification rules from the task.

    Returns:
        Dict with keys: best_candidate, classification, message.
    """
    best = sorted(cands, key=lambda c: c.pct_error_deg)[0]

    # Threshold: 13.04 ± 0.3° in degrees ⇒ |Δθ| ≤ 0.3 ⇒ ~2.3% relative error
    # Rule from task:
    #  - clean substrate-derived & ≤ 2.3% (≤ 0.3°) → MAJOR result
    #  - close but not exact → identify missing structure
    #  - no clean match → open

    msg_lines = []
    msg_lines.append("=" * 95)
    msg_lines.append("HONEST ASSESSMENT")
    msg_lines.append("=" * 95)
    msg_lines.append(f"Best candidate: {best.name}")
    msg_lines.append(f"  formula:          {best.formula}")
    msg_lines.append(f"  family:           {best.family}")
    msg_lines.append(f"  predicted θ_C:    {best.theta_deg:.4f}°  (sin = {best.sin_theta:.4f})")
    msg_lines.append(f"  empirical θ_C:    {THETA_C_DEG:.4f}°    (sin = {SIN_THETA_C:.4f})")
    msg_lines.append(f"  % error (deg):    {best.pct_error_deg:.3f}%")
    msg_lines.append(f"  substrate basis:  {best.substrate_basis}")

    basis_lower = best.substrate_basis.lower()
    has_clean_basis = not any(
        marker in basis_lower
        for marker in (
            "no obvious substrate origin",
            "reported for completeness",
            "not yet a substrate-first",
            "no first-principles motivation",
            "no deep topological motivation",
        )
    )

    if best.pct_error_deg <= 2.3:
        # Could be a major result — but check that it's substrate-derived,
        # not a tautology like the GIM identity (which uses input quark masses
        # that are themselves free Yukawas in our model per §18.49.6).
        if best.family == "A_GIM":
            classification = "GOOD numerical match but NOT first-principles"
            msg_lines.append("")
            msg_lines.append(
                f"VERDICT: {classification}.  The GIM identity reduces θ_C to the "
                f"d/s quark mass ratio, but in our model those quark masses are "
                f"6 free Yukawas (§18.49.6) — same status as the SM.  This is "
                f"a SHIFT of the open question, not a derivation."
            )
        elif not has_clean_basis:
            classification = "NUMERICAL NEAR-MISS — no clean substrate denominator"
            msg_lines.append("")
            msg_lines.append(
                f"VERDICT: {classification}.  The angle lands within 0.3°, "
                f"but the denominator/origin is not selected by the current "
                f"Möbius + Z₃ substrate rules.  This is useful evidence for "
                f"the flavor-mixing scale, not a first-principles CKM "
                f"derivation.  Per §18.49.7, CKM angles remain empirical "
                f"until the flavor-mixing operator is derived."
            )
        else:
            classification = "MAJOR — clean substrate match within 0.3°"
            msg_lines.append("")
            msg_lines.append(
                f"VERDICT: {classification}.  This is a topology-only candidate "
                f"with no empirical input.  If the substrate basis is rigorous, "
                f"this is a genuine first-principles prediction."
            )
    elif best.pct_error_deg <= 5.0:
        classification = "CLOSE — within 5%, missing structure"
        msg_lines.append("")
        msg_lines.append(
            f"VERDICT: {classification}.  The closest substrate angle is "
            f"~{best.theta_deg:.2f}° vs empirical {THETA_C_DEG:.2f}° (off by "
            f"{abs(best.theta_deg - THETA_C_DEG):.2f}°).  Missing structure: "
            f"likely the second-order corrections to the Möbius bundle's "
            f"Dirac equation (an open Path B computation per §18.11)."
        )
    else:
        classification = "OPEN — no clean substrate match"
        msg_lines.append("")
        msg_lines.append(
            f"VERDICT: {classification}.  No purely topological combination "
            f"of {{half-flux π, vertex-closure 2π/3, Chern 1/2}} reproduces "
            f"θ_C ≈ 13.04° to better than {best.pct_error_deg:.1f}%.  Per "
            f"§18.49.7, the CKM angles remain empirical input in our model."
        )

    return {
        "best_candidate": best,
        "classification": classification,
        "report": "\n".join(msg_lines),
    }


# ---------------------------------------------------------------------------
# Module entry point
# ---------------------------------------------------------------------------


def run_full_analysis(verbose: bool = True) -> dict:
    """Run the complete Cabibbo-angle multi-hypothesis test.

    Returns:
        Dict with all candidates, top-N, and honest assessment.
    """
    cands = all_candidates()

    full_report = report_candidates(cands)
    best5 = best_n(cands, n=5)
    assess = honest_assessment(cands)

    if verbose:
        print(full_report)
        print()
        print(assess["report"])

    return {
        "all_candidates": cands,
        "top5": best5,
        "assessment": assess,
        "report_string": full_report,
    }


if __name__ == "__main__":
    run_full_analysis(verbose=True)
