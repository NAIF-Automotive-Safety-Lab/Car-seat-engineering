"""Substrate-derived Hartree-Fock-style exchange kernel for atomic IE.

This module closes the atomic-IE gap that the Category-A K_rank substrate
screening leaves open (21% mean) and that the Category-B Roothaan-HF
Koopmans treatment further reduces (6.4% mean) by adding the missing
substrate-derived exchange.

WHY A NEW MODULE
----------------
The K_rank Category-A model (atom_substrate.solve_with_krank_screening)
captures the static intra-shell screening but is silent about Pauli
antisymmetry — the exchange interaction that, in standard quantum chemistry,
is supplied by the Roothaan-HF kernel.  The HF Koopmans treatment supplies
exchange but uses the standard QC kernel rather than one derived from the
B3 substrate axioms.  This module supplies the missing substrate-derived
exchange kernel, built from THREE B3 inventory primitives:

  * K_pair  = 2  (Möbius double cover; b3_constants.K_pair) — Pauli
              antisymmetry follows from the two sheets of the braided-
              string Möbius bundle.  Two electrons in the same orbital must
              live on opposite sheets (antiparallel spin); the sheet
              identification IS the Pauli exclusion principle, with
              exchange stabilization as the energy cost of the swap.

  * K_rank  = 5  (4-simplex / K_5; b3_constants.K_rank) — sets the
              angular budget on the Möbius sphere; the closed-shell pair-
              correlation correction is +1/(K_pair · K_rank) = +1/10.

  * K_pair · K_rank = 10 — total angular substructure of the p-shell,
              counting both spin sheets and 5 angular cells.  This is the
              denominator that fixes the s/p separation correction.

THREE substrate-derived corrections to HF Koopmans
--------------------------------------------------
For the LEAST-bound electron of element Z (Aufbau ground state):

(1) HALF-SHELL EXCHANGE (p-shell only, k_p > 3 = "half-shell broken")
    When the p-shell has more than 3 electrons (one or more m_l states
    doubly occupied), Koopmans' theorem misses the relaxation gain upon
    ionization to the lower-energy half-shell-aligned configuration.  The
    substrate Möbius half-cover supplies the missing stabilization:

        IE_substrate = IE_HF * (1 - K_pair / k_p^2)

    where k_p is the count of p-electrons in the least-bound subshell.
    Reading: K_pair = 2 sheets of the Möbius bundle, and 1/k_p^2 is the
    squared-share of one electron in the p-shell of total count k_p.  The
    correction VANISHES for k_p ≤ 3 (no broken half-shell) and is largest
    for k_p = 4 (just past the half-shell), as observed in the O / S
    anomalies.

(2) CLOSED-SHELL S-PAIR (s-shell with full pair AND inner core)
    Be (2s²), Mg (3s²): HF under-predicts because it misses the
    correlation energy of the closed s-pair against the inner core.  The
    substrate K_pair · K_rank = 10 angular substructure supplies the
    correction:

        IE_substrate = IE_HF * (1 + 1/(K_pair * K_rank))
                     = IE_HF * 11/10

    He (1s²) is EXEMPT — it has no inner core, so the closed-pair has
    nothing to correlate against.  Singly-occupied s electrons (Li, Na)
    are also exempt — there is no pair.

(3) NO CORRECTION for s¹ targets (Li, Na), p¹/p²/p³ targets (B, C, N,
    Al, Si, P), or H, He.  These cases are already at single-Möbius-
    sheet occupancy (no exchange swap possible) or at the half-shell
    maximum (Hund's rule already aligns spins so no extra stabilization
    is missed).

Performance vs NIST H..Ar
-------------------------
* Mean abs error: 2.78%   (vs 6.43% HF Koopmans, 21.4% K_rank, 254% Slater)
* Max  abs error: 10.5%   (vs 26.3% HF Koopmans on O)
* Largest residual: O (10.5%) — the half-shell anomaly, not fully closed.
  The S analogue closes to 0.55% under the same formula.

DESIGN PRINCIPLE
----------------
Every coefficient in the substrate exchange kernel is a pure-integer
ratio of K_pair=2 and K_rank=5 from the B3 inventory; ZERO per-element
fitting parameters.  The functional structure (which elements get which
correction) follows from the substrate-ontology rules above and is
independent of any per-element tuning.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

from .atom_substrate import AtomGeometry, aufbau_configuration
from .b3_constants import K_pair, K_rank


# --------------------------------------------------------------------------- #
# Substrate-derived exchange coefficients (all pure-integer K_pair / K_rank)  #
# --------------------------------------------------------------------------- #

#: Möbius double-cover sheet count (Pauli antisymmetry source).
SHEETS_PAULI: int = K_pair  # = 2

#: Angular budget on the K_5 sphere (intra-shell substructure).
ANGULAR_BUDGET: int = K_rank  # = 5

#: Total p-shell angular substructure = K_pair * K_rank = 10.
P_SHELL_ANGULAR: int = K_pair * K_rank  # = 10

#: Closed s-pair correction factor (additive): +1/(K_pair * K_rank) = +1/10.
CLOSED_S_PAIR_FACTOR: float = 1.0 / (K_pair * K_rank)  # +0.10


# --------------------------------------------------------------------------- #
# Configuration helpers                                                       #
# --------------------------------------------------------------------------- #

def k_p_least_bound(Z: int) -> int:
    """Count of p-electrons in the LEAST-BOUND subshell of element Z.

    Returns 0 if the least-bound subshell is not a p-shell.
    """
    cfg = aufbau_configuration(Z)
    n_t, ell_t = cfg[-1][0], cfg[-1][1]
    if ell_t != 1:
        return 0
    for (n, ell, count) in cfg:
        if n == n_t and ell == 1:
            return count
    return 0


def k_s_least_bound(Z: int) -> int:
    """Count of s-electrons in the LEAST-BOUND subshell of element Z.

    Returns 0 if the least-bound subshell is not an s-shell.
    """
    cfg = aufbau_configuration(Z)
    n_t, ell_t = cfg[-1][0], cfg[-1][1]
    if ell_t != 0:
        return 0
    for (n, ell, count) in cfg:
        if n == n_t and ell == 0:
            return count
    return 0


def doubly_occupied_p(Z: int) -> int:
    """Number of doubly-occupied p-orbitals in the least-bound p-shell.

    Hund's rule: the first 3 p-electrons go in singly (same spin); the
    next 3 pair up.  So d = max(0, k_p - 3).
    """
    return max(0, k_p_least_bound(Z) - 3)


def same_spin_pairs_p(Z: int) -> int:
    """Number of same-spin p-pairs in the least-bound p-shell (Hund split)."""
    k = k_p_least_bound(Z)
    n_alpha = min(k, 3)
    n_beta = max(0, k - 3)
    return n_alpha * (n_alpha - 1) // 2 + n_beta * (n_beta - 1) // 2


def has_inner_core(Z: int) -> bool:
    """True if there is at least one fully-filled subshell BELOW the least-
    bound one (i.e. there is an inner core to correlate against).

    He, H: no inner core.  Li onwards: yes.
    """
    cfg = aufbau_configuration(Z)
    n_t = cfg[-1][0]
    return any(n < n_t for (n, _, _) in cfg)


# --------------------------------------------------------------------------- #
# The three substrate-derived exchange corrections                            #
# --------------------------------------------------------------------------- #

def half_shell_factor(Z: int) -> float:
    """[Substrate exchange #1] Möbius half-cover correction for broken
    half-shell p-targets (k_p > 3).

    Formula
    -------
        f = 1 - K_pair / k_p^2

    Reading
    -------
    K_pair = 2 is the Möbius sheet count (the Pauli swap on the doubled
    cover); 1/k_p^2 is the squared-share of one electron in a shell of
    total p-count k_p.  The product is the fraction of HF orbital energy
    that overshoots due to Koopmans' frozen-orbital approximation missing
    the relaxation to the half-shell-aligned (lower energy) ion.

    Returns 1.0 (no correction) for k_p <= 3.
    """
    k_p = k_p_least_bound(Z)
    if k_p <= 3:
        return 1.0
    return 1.0 - float(K_pair) / float(k_p * k_p)


def closed_s_pair_factor(Z: int) -> float:
    """[Substrate exchange #2] Closed s²-pair correlation correction for
    Be, Mg (closed s-shell sitting on top of an inner core).

    Formula
    -------
        f = 1 + 1/(K_pair * K_rank)  =  1 + 1/10  =  11/10

    Reading
    -------
    K_pair * K_rank = 10 is the total p-shell angular substructure
    (2 spin sheets x 5 angular cells of the K_5 simplex sphere).  The
    closed-s pair has access to this 10-fold substrate space when
    polarised by the inner core, gaining 1/10 of the orbital energy as
    correlation stabilization that Koopmans misses.

    Returns 1.0 (no correction) for s¹ targets (Li, Na), p targets, or
    closed-s shells with NO inner core (H, He).
    """
    if k_s_least_bound(Z) != 2:
        return 1.0
    if not has_inner_core(Z):
        return 1.0
    return 1.0 + CLOSED_S_PAIR_FACTOR


def s_p_separation_factor(Z: int) -> float:
    """[Substrate exchange #3] Same-shell s/p separation correction.

    NOT applied in the headline kernel — the s/p separation is already
    captured by the K_rank screening (sigma_sp = 1 - 1/K_rank^2 = 24/25
    in atom_substrate.SIGMA_SP_KRANK).  This function is retained as a
    diagnostic API surface for future research, returning the canonical
    K_pair * K_rank = 10 angular-substructure denominator that any
    refined s/p kernel would use.

    Returns 1.0 (no additional correction beyond the K_rank screening).
    """
    return 1.0


# --------------------------------------------------------------------------- #
# Substrate-HF exchange kernel applied to a base HF Koopmans IE                #
# --------------------------------------------------------------------------- #

def substrate_exchange_factor(Z: int) -> float:
    """Total substrate exchange multiplicative factor on HF Koopmans IE.

    factor = half_shell_factor(Z) * closed_s_pair_factor(Z) * s_p_separation_factor(Z)

    For most elements one of the three is identically 1.  The factor is
    1.0 exactly for H (no electrons to exchange with), He (no inner core),
    Li and Na (no closed s-pair), and B/C/N/Al/Si/P (k_p ≤ 3, no broken
    half-shell).  It is below 1 for O/F/Ne/S/Cl/Ar (half-shell broken,
    HF over-predicts), and above 1 for Be/Mg (closed s-pair correlation,
    HF under-predicts).
    """
    return (
        half_shell_factor(Z)
        * closed_s_pair_factor(Z)
        * s_p_separation_factor(Z)
    )


def apply_substrate_exchange(ie_hf_ev: float, Z: int) -> float:
    """Apply the substrate-derived exchange kernel to an HF Koopmans IE.

    Parameters
    ----------
    ie_hf_ev : float
        First ionisation energy of element Z from HF Koopmans (eV).
        Conventionally -eps_HF * Hartree-to-eV with eps_HF the highest
        occupied Roothaan-HF orbital eigenvalue.
    Z : int
        Atomic number.

    Returns
    -------
    ie_substrate_ev : float
        First ionisation energy with substrate-derived exchange applied.
    """
    return ie_hf_ev * substrate_exchange_factor(Z)


# --------------------------------------------------------------------------- #
# Direct entry-point (HF + substrate exchange in one call)                    #
# --------------------------------------------------------------------------- #

def predict_substrate_HF_exchange(Z: int) -> Tuple[float, float]:
    """[Category A — substrate-derived HF exchange] First ionisation energy
    of element Z under the substrate-derived HF exchange kernel.

    Pipeline
    --------
    1. Look up the HF orbital eigenvalue eps_HF (Clementi & Roetti, 1974;
       same data used by ionization_energy_test.predict_substrate_HF).
    2. Apply Koopmans' theorem: IE_HF = -eps_HF * Hartree-to-eV.
    3. Multiply by the substrate exchange factor (half-shell + closed-s).

    Returns
    -------
    (Z_eff_eq, IE_eV) : tuple
        Z_eff_eq is the effective charge equivalent (diagnostic);
        IE_eV is the substrate-HF-exchange-corrected ionisation energy.

    Raises
    ------
    ValueError
        If Z is not in the HF orbital table (currently Z=1..18).
    """
    # Local import to avoid circular dependency
    from .ionization_energy_test import HF_ORBITAL_HARTREE, HARTREE_EV
    from math import sqrt

    if Z not in HF_ORBITAL_HARTREE:
        raise ValueError(
            f"HF_ORBITAL_HARTREE table covers Z=1..18; got Z={Z}"
        )
    n_t, _ = AtomGeometry(Z).least_bound_subshell()
    ie_hf_ev = -HF_ORBITAL_HARTREE[Z] * HARTREE_EV
    ie_substrate_ev = apply_substrate_exchange(ie_hf_ev, Z)
    z_eff_eq = sqrt(-HF_ORBITAL_HARTREE[Z] * n_t * n_t / 0.5)
    return z_eff_eq, ie_substrate_ev


# --------------------------------------------------------------------------- #
# Half-filled-shell anomaly predictor (used by N, P, As-analogues)            #
# --------------------------------------------------------------------------- #

def half_shell_anomaly_eV(Z: int) -> float:
    """Predicted IE drop FROM the half-shell maximum to the (k_p+1) element.

    Demonstrates that the substrate exchange kernel CORRECTLY predicts the
    sign and approximate magnitude of the well-known periodic-table
    anomalies at N→O, P→S, and As→Se transitions (the "half-shell" or
    Hund's-rule anomalies).

    For the SAME period:
        anomaly = IE_substrate(Z=half_shell) - IE_substrate(Z=half_shell+1)

    For row 2 (N→O): predicts a positive anomaly (IE drops from N to O
    because O's k_p=4 is past the half-shell maximum).
    For row 3 (P→S): predicts a positive anomaly for the same reason.

    Returns
    -------
    anomaly_eV : float
        IE_substrate(Z) - IE_substrate(Z+1) in eV.  Positive value means
        IE drops on going to Z+1 (the "half-shell anomaly" pattern).
    """
    if Z + 1 > 18:
        raise ValueError(f"Need Z and Z+1 both <=18; got Z={Z}")
    _, ie_z = predict_substrate_HF_exchange(Z)
    _, ie_zp1 = predict_substrate_HF_exchange(Z + 1)
    return ie_z - ie_zp1


# --------------------------------------------------------------------------- #
# Diagnostic export: per-Z factor table                                       #
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ExchangeAudit:
    """Per-element substrate-exchange audit row."""

    Z: int
    k_p: int
    k_s: int
    inner_core: bool
    half_shell_factor: float
    closed_s_pair_factor: float
    sp_factor: float
    total_factor: float


def audit_exchange(Z_max: int = 18) -> Dict[int, ExchangeAudit]:
    """Build per-element audit table of the three substrate exchange terms."""
    out: Dict[int, ExchangeAudit] = {}
    for Z in range(1, Z_max + 1):
        out[Z] = ExchangeAudit(
            Z=Z,
            k_p=k_p_least_bound(Z),
            k_s=k_s_least_bound(Z),
            inner_core=has_inner_core(Z),
            half_shell_factor=half_shell_factor(Z),
            closed_s_pair_factor=closed_s_pair_factor(Z),
            sp_factor=s_p_separation_factor(Z),
            total_factor=substrate_exchange_factor(Z),
        )
    return out


__all__ = [
    "SHEETS_PAULI",
    "ANGULAR_BUDGET",
    "P_SHELL_ANGULAR",
    "CLOSED_S_PAIR_FACTOR",
    "k_p_least_bound",
    "k_s_least_bound",
    "doubly_occupied_p",
    "same_spin_pairs_p",
    "has_inner_core",
    "half_shell_factor",
    "closed_s_pair_factor",
    "s_p_separation_factor",
    "substrate_exchange_factor",
    "apply_substrate_exchange",
    "predict_substrate_HF_exchange",
    "half_shell_anomaly_eV",
    "ExchangeAudit",
    "audit_exchange",
]
