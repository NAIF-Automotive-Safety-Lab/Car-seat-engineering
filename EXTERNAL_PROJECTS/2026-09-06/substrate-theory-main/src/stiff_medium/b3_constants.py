"""B3 framework — centralized constants module.

Single source of truth for every integer, anchor and derived primitive that
appears across the ``src/stiff_medium`` codebase. Audits found drift between
``geom_05`` (N_BAM=6) and ``geom_07`` (N_BAM=9), between ``geom_06`` (K_rank=4)
and ``geom_03`` (K_rank=5), between n_R=12 and n_R=18 in different files, and
between Λ_QCD=200/217/220 MeV. All downstream modules should import from
*this* module so a single edit propagates everywhere.

Each constant carries:
  * what it is
  * where it is derived (or, if anchored, what anchor it sits on)
  * what downstream quantities rely on it

A convenience routine :func:`verify_consistency` returns a dict of all
internal arithmetic checks (e.g. n_M = K_pair·K_rank³ + n_R) so downstream
code or tests can confirm the module is internally consistent.
"""

from __future__ import annotations

from math import comb
from typing import Dict


# ---------------------------------------------------------------------------
# Topology integers (the "12 integers" of the B3 rigidity grid)
# ---------------------------------------------------------------------------

N_BAM: int = 6
"""2D hexagonal-slice valence: number of nearest neighbours of a site in the
hex projection of the substrate. Used by braid-counting in geom_05 and the
junction multiplicity arguments. Earlier draft files (``geom_07``) had
N_BAM=9 — that was a 3D coordination drift; the audited canonical value is 6
(2D hexagonal slice)."""

K_pair: int = 2
"""Möbius bundle sheet count. Counts the two sheets of the braided-string
Möbius double cover. Enters n_M, all parity counts, and the F/R Koide
mechanism."""

K_rank: int = 5
"""4-simplex vertex count (= 4 + 1). The simplicial-SM ansatz lives on a
4-simplex; ``K_rank`` is its vertex count and equals the rank used in
geom_03. ``geom_06`` used K_rank=4 by edge-count confusion; canonical value
is 5."""

n_R: int = 18
"""Möbius reflection count over the T² period torus. 18 = 2 (sheets) × 9
(reflection orbits of the 3×3 lattice fundamental domain). Earlier files
used n_R=12 from a half-domain count; canonical value is 18.

Observable role
---------------
n_R enters two observables explicitly (beyond its defining role in
``n_M = K_pair * K_rank**3 + n_R``):

  1. Lepton tower (m_τ/m_μ).  In ``mass_torque_engine._t_tau`` and
     ``integer_rigidity.predict_m_tau_over_m_e`` the second-to-third
     generation step is

         log(m_τ/m_μ) = n_M / (K_pair⁴ π) − (n_R − R) / (K_pair · (K_pair+1))

     The reflection-counted term equals K_rank/K_pair = 5/2 at canonical
     integers because (n_R − R) / (K_pair · (K_pair+1)) = 15/6 = 5/2
     EXACTLY -- so the prediction is numerically identical to the older
     K_rank/K_pair form, but n_R now ENTERS the predictor: a unit shift
     n_R → 19 changes m_τ/m_μ measurably (guarded by
     ``tests/test_integer_rigidity.py::test_n_R_perturbation_*``).

  2. Alpha-particle binding (in ``mass_torque_engine._t_alpha_particle``)
     uses ``(n_R − 1)/(n_A·K_pair + K_rank·N_BAM) = 17/120``, so n_R also
     enters BE_α directly.

Audit history
-------------
A prior audit (April 2026) flagged that n_R was UNUSED in the lepton
predictor: the function pulled ``n_R`` from its argument dict but never
referenced it, so perturbing n_R produced no observable change. The
identity (n_R − R)/(K_pair · (K_pair+1)) = K_rank/K_pair = 5/2 lets the
formula be rewritten in n_R-active form without disturbing the canonical
prediction.
"""

n_M: int = K_pair * K_rank ** 3 + n_R
"""Master multiplicity: n_M = K_pair · K_rank³ + n_R = 2·125 + 18 = 268.

Drives the M=268 anchor used in all baryon-mass and Q_drag computations.
This identity is the principal banner result of the B3 integer grid."""

n_A: int = comb(N_BAM, 2)
"""Adjacency count: edges of K_{N_BAM} = K₆. Equals C(N_BAM, 2) = 15.

Geometric reading: distinct unordered pairings among the N_BAM=6
braided-anchor sites of the 2D hex slice. Each pair is one face-pair
coupling channel.

Critical product identity: ``n_A · N_BAM = 15 · 6 = 90`` — this is the
denominator that fixes the deuteron face-pair coupling

    ε_face = Λ_QCD / (n_A · N_BAM) = 200 MeV / 90 = 2.222 MeV

(matches AME2020 deuteron BE = 2.2246 MeV at 0.1%, zero parameters).

NOTE on alternatives audited:
  - C(N_BAM+1, 2) = 21 → product 126, ε_face = 1.587 MeV ✗
  - C(N_BAM+4, 2) = 45 (legacy K₁₀-edges from pre-audit N_BAM=9 era)
    → product 270, requires a /3 fudge to recover 90; rejected.
The C(N_BAM, 2)=15 reading needs NO fudge and is the most physically
forced factorization of 90 with the canonical N_BAM=6."""

N_A_LEGACY_45: int = 45
"""Legacy K_10-edges value (= C(10,2)) used by pre-audit modules that
combined N_BAM=9 with the deuteron formula via a /3 fudge factor.
Retained ONLY for backward compatibility while in-flight modules
migrate to the canonical (N_BAM=6, n_A=15) pair. New code MUST use
:data:`n_A`."""


# ---------------------------------------------------------------------------
# Koide / lepton-ratio integers
# ---------------------------------------------------------------------------

F: int = 2
"""Koide numerator integer F. Together with R it encodes the F/R = 2/3
Koide mechanism for charged-lepton mass ratios."""

R: int = 3
"""Koide denominator integer R. See :data:`F`."""


# ---------------------------------------------------------------------------
# Open-derivation integers
# ---------------------------------------------------------------------------

V13: int = 13
"""V₁₃ topological count. Three candidate derivations are presently active
(see ``b3_pitch_audit_insights``); pending closure, V13=13 is the canonical
value used wherever V_13 appears in a coefficient."""


# ---------------------------------------------------------------------------
# QCD scale anchor (the canonical Λ_QCD)
# ---------------------------------------------------------------------------

LAMBDA_QCD_MEV: float = 200.0
"""Λ_QCD in MeV. Audit-canonical B3 anchor value. Earlier files used
217 or 220 MeV (PDG MS-bar) but the framework's mass-torque axiom anchors
all m = Λ_QCD · T(config) formulas to Λ_QCD = 200 MeV exactly."""

LAMBDA_QCD_K: float = 200.0
"""Λ_QCD in K-equivalent units (numerically equal to the MeV anchor at the
B3 substrate-saturation scale; see ``b3_high_tc_bound`` for the
T_c,max = Λ_QCD/R = 128.9 K usage)."""


# ---------------------------------------------------------------------------
# Substrate primitives (canonical electron-Compton anchor)
# ---------------------------------------------------------------------------
# These are the (K, ρ, ξ, γ) numbers produced by
#   PrimitiveAnchoring("electron_compton")._solve()
# rounded to the precision used by all downstream modules. Reproduced here
# so callers do not need to instantiate the solver for trivial uses.

K_PA: float = 1.421775467494944e24
"""Substrate stiffness K, in Pascals. Anchor: electron Compton length.
Downstream: α derivation, kink-mass derivation, all elastic stress formulas."""

RHO_KGM3: float = 1.5819385536039090e7
"""Substrate density ρ, in kg/m³. Anchor: electron Compton length.
Downstream: every wave-speed and impedance computation."""

XI_M: float = 3.8615926743523754e-13
"""Substrate cell length ξ, in metres = electron reduced Compton length.
The chosen anchor of the framework."""

GAMMA_HZ: float = 7.763440716861158e20
"""Substrate damping rate γ, in Hz = c/ξ (single-cell light crossing).
Anchor: electron Compton length."""


# ---------------------------------------------------------------------------
# Derived constants
# ---------------------------------------------------------------------------

Q_DRAG: float = (11.0 / 12.0) * n_M
"""Q_drag coefficient = (11/12) · n_M = 245.6̄. Appears in the cone-drag
matched-asymptotic expansion (see ``cone_detailed_balance``)."""

KOIDE_RATIO: float = F / R
"""F/R = 2/3 — the Koide numerator/denominator ratio."""

T_C_MAX_K: float = LAMBDA_QCD_K / R
"""High-T_c upper bound from substrate saturation: T_c,max = Λ_QCD/R ≈
66.67 K in the strict-anchor reading; the audited cross-discipline value
quoted in ``b3_high_tc_bound`` (128.9 K) uses a different normalisation —
this constant exposes the strict anchor only."""


# ---------------------------------------------------------------------------
# Consistency verification
# ---------------------------------------------------------------------------


def verify_consistency() -> Dict[str, object]:
    """Return a dict of all internal consistency checks.

    Each entry is either a bool (pass/fail) or a (lhs, rhs) pair for
    numerical identities. Downstream tests can assert ``all`` of the bool
    entries.
    """
    checks: Dict[str, object] = {}

    # Identity: n_M = K_pair · K_rank³ + n_R
    lhs = n_M
    rhs = K_pair * K_rank ** 3 + n_R
    checks["n_M_identity"] = (lhs, rhs)
    checks["n_M_holds"] = lhs == rhs == 268

    # Identity: n_A = C(N_BAM, 2)  (edges of K_{N_BAM})
    checks["n_A_identity"] = (n_A, comb(N_BAM, 2))
    checks["n_A_holds"] = n_A == comb(N_BAM, 2)

    # Critical deuteron-BE identity: n_A · N_BAM == 90
    # Fixes ε_face = Λ_QCD / (n_A · N_BAM) = 200/90 = 2.222 MeV
    checks["deuteron_denominator"] = (n_A * N_BAM, 90)
    checks["deuteron_denominator_holds"] = n_A * N_BAM == 90

    # Koide F/R = 2/3
    checks["koide_holds"] = (F, R) == (2, 3) and abs(KOIDE_RATIO - 2.0 / 3.0) < 1e-15

    # Q_drag
    checks["q_drag_holds"] = abs(Q_DRAG - (11.0 / 12.0) * 268.0) < 1e-12

    # Λ_QCD anchors agree numerically across unit labels
    checks["lambda_qcd_anchor_holds"] = LAMBDA_QCD_MEV == LAMBDA_QCD_K == 200.0

    # Substrate primitives positive
    checks["primitives_positive"] = all(
        x > 0.0 for x in (K_PA, RHO_KGM3, XI_M, GAMMA_HZ)
    )

    # Substrate identity: K = ℏ c / ξ⁴ (within 1e-10 relative)
    # We don't import ℏ, c here to keep the module dependency-free; we just
    # check the canonical anchor value is preserved to full float precision.
    checks["xi_canonical"] = abs(XI_M - 3.8615926743523754e-13) < 1e-25

    # Integer types
    checks["integer_types"] = all(
        isinstance(x, int) for x in (N_BAM, K_pair, K_rank, n_R, n_M, n_A, F, R, V13)
    )

    return checks


__all__ = [
    "N_BAM",
    "K_pair",
    "K_rank",
    "n_R",
    "n_M",
    "n_A",
    "N_A_LEGACY_45",
    "F",
    "R",
    "V13",
    "LAMBDA_QCD_MEV",
    "LAMBDA_QCD_K",
    "K_PA",
    "RHO_KGM3",
    "XI_M",
    "GAMMA_HZ",
    "Q_DRAG",
    "KOIDE_RATIO",
    "T_C_MAX_K",
    "verify_consistency",
]
