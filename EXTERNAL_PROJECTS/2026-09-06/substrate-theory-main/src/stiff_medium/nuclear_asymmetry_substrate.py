"""Substrate-derived nuclear asymmetry coefficient a_sym.

Derives the SEMF asymmetry coefficient a_sym ≈ 23 MeV (the empirical
Bohr–Mottelson / Wapstra value) from substrate primitives, *without*
the textbook Fermi-gas + isovector-potential calibration.

DERIVATION CHAIN (Category A: substrate-derived)
------------------------------------------------
The K_4 cell has K_rank = 5 vertices (4 quark faces + 1 apex closure)
and carries K_pair = 2 isospin sheets (Möbius bundle u/d).  When an
isobar has N ≠ Z the cell stack cannot match all faces u↔d-symmetrically:
some face pairs end up forced into same-isospin (n-n or p-p) couplings
which carry the full per-cell torque cost rather than the deuteron-style
mixed-isospin coupling.

Per unmatched-pair cost = (face-pair binding) × (Möbius sheets) × (vertex count)

    a_sym = ε_face · K_pair · K_rank
          = (Λ_QCD / (n_A · N_BAM)) · K_pair · K_rank
          = (200 / 90) · 2 · 5
          = 200 / 9
          = 22.222… MeV

Compare with empirical Bohr–Mottelson / Wapstra a_sym = 23 MeV:
    relative error  =  (22.222 − 23) / 23  =  −3.4 %.

The (N − Z)² / A scaling that multiplies a_sym in the SEMF arises from
the standard combinatorial argument (fraction of forced same-isospin
pairs in a stack with given isospin imbalance).  This module assumes the
SEMF combinatorial form and supplies only the coefficient.

EQUIVALENT EXPRESSIONS (numerically identical, all = Λ_QCD / 9)
---------------------------------------------------------------
1.  ε_face · K_pair · K_rank                        — primitive product
2.  Λ_QCD · (K_pair · K_rank) / (n_A · N_BAM)       — fully Λ-anchored
3.  Λ_QCD / (K_pair · K_rank − 1)                   — direct integer denominator
4.  Λ_QCD / R²                                       — Koide-denominator squared

Form 1 is the cleanest physical reading: it expresses a_sym as the
deuteron face-pair binding (ε_face) multiplied by the two integer
counts that distinguish "isospin-mixed" pairs (K_pair = u/d sheets)
from "rank-closed" cells (K_rank = 4-simplex vertices).

CATEGORY UPGRADE
----------------
This module promotes a_sym from Category B (in-principle substrate-
derivable) to Category A (substrate-derived from primitives) in
``analysis/substrate_completion_roadmap.md``, sector 1.2.

HONEST CAVEAT
-------------
The combinatorial (N-Z)²/A form is *imported* from the standard SEMF
derivation; the substrate work here only fixes the COEFFICIENT.  A
fully substrate-native re-derivation of the (N-Z)²/A form (rather
than just the prefactor) would require an explicit count of forced
same-isospin face-pairs in an N ≠ Z stack — a topological calculation
that has not been done.  The 3.4 % match to the empirical value is
striking but does not in itself force the (K_pair · K_rank) factor:
several other inventory products (e.g. ε_face × N_BAM × K_pair = 26.67 MeV
at +15.9 %, ε_face × (K_pair·K_rank+1) = 24.44 MeV at +6.3 %) also land
in the right ballpark.  The K_pair · K_rank form is preferred because
of its direct physical reading (Möbius sheets × simplex vertices), not
because alternatives have been ruled out.

REFERENCES
----------
* :data:`src.stiff_medium.b3_constants.LAMBDA_QCD_MEV` — 200 MeV anchor.
* :data:`src.stiff_medium.b3_constants.K_pair`, K_rank — integer primitives.
* :data:`src.stiff_medium.k4_face_pair_geometry.EPS_FACE_MEV` — 2.222 MeV.
* `analysis/substrate_completion_roadmap.md` — sector 1.2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from src.stiff_medium.b3_constants import (
    K_pair,
    K_rank,
    LAMBDA_QCD_MEV,
    N_BAM,
    R,
    n_A,
)
from src.stiff_medium.k4_face_pair_geometry import EPS_FACE_MEV


# ---------------------------------------------------------------------------
# Empirical reference (NOT used as input — only for error reporting)
# ---------------------------------------------------------------------------

A_SYM_EMPIRICAL_MEV: Final[float] = 23.0
"""Bohr–Mottelson / Wapstra empirical SEMF asymmetry coefficient (MeV).

Literature range across major SEMF fits:

    Bohr–Mottelson (1969)  :  23.0 MeV
    Wapstra (1971)         :  23.7 MeV
    Möller–Nix (1995)      :  25.2 MeV  (with surface-asymmetry split)
    Myers–Swiatecki (1966) :  21.0 MeV  (with separate surface term)

The substrate prediction 22.22 MeV sits inside the literature range.
"""


# ---------------------------------------------------------------------------
# Substrate-derived a_sym
# ---------------------------------------------------------------------------

def a_sym_substrate() -> float:
    """Substrate-derived SEMF asymmetry coefficient (MeV).

    Returns:
        a_sym = ε_face · K_pair · K_rank
              = (Λ_QCD / (n_A · N_BAM)) · K_pair · K_rank
              = 200 / 9
              ≈ 22.222 MeV

    Physical reading: each unmatched p-n face-pair in an N ≠ Z stack
    costs the deuteron face-pair binding (ε_face = 2.222 MeV) times the
    Möbius sheet count (K_pair = 2) times the 4-simplex vertex count
    (K_rank = 5).  Compare empirical Bohr–Mottelson 23 MeV (3.4 % match).
    """
    return EPS_FACE_MEV * K_pair * K_rank


def a_sym_lambda_form() -> float:
    """Equivalent Λ-anchored form of a_sym (MeV).

    a_sym = Λ_QCD · (K_pair · K_rank) / (n_A · N_BAM)
          = 200 · 10 / 90  =  Λ_QCD / 9

    Numerically identical to :func:`a_sym_substrate`; provided for the
    Λ-anchored audit trail.
    """
    return LAMBDA_QCD_MEV * (K_pair * K_rank) / (n_A * N_BAM)


def a_sym_koide_form() -> float:
    """Equivalent Koide-denominator form (MeV).

    a_sym = Λ_QCD / R²  with R = 3 the Koide F/R denominator.

    Numerically identical to :func:`a_sym_substrate` because
    K_pair · K_rank − 1 = 9 = R²; provided for the Koide-mechanism
    audit trail.
    """
    return LAMBDA_QCD_MEV / (R * R)


def a_sym_relative_error() -> float:
    """Fractional error vs empirical 23 MeV: (a_sym_substrate − 23) / 23."""
    return (a_sym_substrate() - A_SYM_EMPIRICAL_MEV) / A_SYM_EMPIRICAL_MEV


# ---------------------------------------------------------------------------
# Comparison record
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AsymmetryDerivation:
    """Substrate-derived a_sym alongside literature variants and error."""

    a_sym_substrate_MeV: float
    a_sym_empirical_MeV: float
    derivation: str

    @property
    def err_pct(self) -> float:
        """100 · (substrate − empirical) / empirical (signed, percent)."""
        return 100.0 * (self.a_sym_substrate_MeV - self.a_sym_empirical_MeV) / self.a_sym_empirical_MeV

    @property
    def err_abs_MeV(self) -> float:
        """|substrate − empirical| (MeV, unsigned)."""
        return abs(self.a_sym_substrate_MeV - self.a_sym_empirical_MeV)


def report() -> AsymmetryDerivation:
    """Build the headline AsymmetryDerivation record."""
    return AsymmetryDerivation(
        a_sym_substrate_MeV=a_sym_substrate(),
        a_sym_empirical_MeV=A_SYM_EMPIRICAL_MEV,
        derivation="ε_face · K_pair · K_rank = 2.222 · 2 · 5",
    )


# ---------------------------------------------------------------------------
# Convenience: print the substrate-derived value with context
# ---------------------------------------------------------------------------

def demo() -> None:
    """Print the headline derivation."""
    rep = report()
    print("Substrate-derived SEMF asymmetry coefficient")
    print("=" * 60)
    print(f"  ε_face          = {EPS_FACE_MEV:.4f} MeV")
    print(f"  K_pair          = {K_pair}")
    print(f"  K_rank          = {K_rank}")
    print()
    print(f"  a_sym (primitive product) = ε · K_p · K_r")
    print(f"                            = {EPS_FACE_MEV:.4f} · {K_pair} · {K_rank}")
    print(f"                            = {a_sym_substrate():.4f} MeV")
    print()
    print(f"  Equivalent forms:")
    print(f"    Λ_QCD · K_pair·K_rank / (n_A·N_BAM) = {a_sym_lambda_form():.4f} MeV")
    print(f"    Λ_QCD / (K_pair·K_rank − 1)         = {LAMBDA_QCD_MEV / (K_pair*K_rank - 1):.4f} MeV")
    print(f"    Λ_QCD / R² (Koide denominator²)     = {a_sym_koide_form():.4f} MeV")
    print()
    print(f"  Empirical (Bohr–Mottelson)            = {rep.a_sym_empirical_MeV:.4f} MeV")
    print(f"  Relative error                        = {rep.err_pct:+.2f} %")
    print(f"  Absolute residual                     = {rep.err_abs_MeV:.4f} MeV")


if __name__ == "__main__":
    demo()
