"""Baryon Regge slope from a relativistic Y-junction of three strings.

Physics summary
---------------
A meson is an open string of length 2L with two quark endpoints moving at
the speed of light.  A baryon is modelled here as a Y-junction: three
strings of length L, each running from a central vertex to a quark at
distance L, with the whole configuration rotating rigidly about the centre.

Both configurations are pure substrate strings — the same tension σ from
the running stiffness K(ξ) — but the geometry differs.  The Regge slope
α' = dJ/dM² depends on geometry, so the meson slope cannot be assumed to
apply to baryons.  This module derives α'_Y from the relativistic Y-junction
classical motion.

Derivation
----------
For a single rotating Nambu-Goto open string (the meson) of total length
2L with c=1 and ω = 1/L (so the endpoints move at c):

    E_meson = 2 ∫_0^L σ dr / sqrt(1 - (r/L)²) = π σ L
    J_meson = 2 ∫_0^L σ (r²/L) dr / sqrt(1 - (r/L)²) = π σ L² / 2

eliminating L:

    M² = π² σ² L²,   J = π σ L² / 2
    ⇒ α'_meson = J / M² = 1 / (2 π σ).

For a Y-junction of THREE strings (each of length L) all rotating with
the same ω = 1/L about the central vertex, the integrals are over each
string separately, with the integrand identical to the half-string of
the meson (one half of the meson is a string from r=0 at v=0 to r=L at
v=c — exactly one Y-junction arm):

    E_Y = 3 ∫_0^L σ dr / sqrt(1 - (r/L)²) = 3 σ L · (π/2) = (3π/2) σ L
    J_Y = 3 ∫_0^L σ (r²/L) dr / sqrt(1 - (r/L)²) = 3 σ L² · (π/4) = (3π/4) σ L²

eliminating L:

    M² = E² = (3π/2)² σ² L² = (9π²/4) σ² L²
    J = (3π/4) σ L²

    ⇒ α'_Y = J / M² = (3π/4 σ L²) / (9π²/4 σ² L²) = 1 / (3 π σ).

So the Y-junction Regge slope is

    α'_Y = (2/3) · α'_meson.

This is a SMALLER slope than the meson — the trajectory is steeper, i.e.
states with the same J are MORE massive than the meson trajectory would
predict.

Honest assessment of what this captures and what it misses
----------------------------------------------------------
What this captures: pure relativistic-string angular momentum from rigid
rotation of a planar Y-junction.  The geometric factor 2/3 is exact for
this idealised configuration.

What this misses, and why the N-to-Δ jump is NOT a string prediction:
The transition N(939, J=1/2) → Δ(1232, J=3/2) is NOT an orbital
excitation.  Both states have L=0; the J change is from a quark spin flip
(symmetric vs antisymmetric SU(6) state).  In QCD this is the
chromomagnetic spin-spin interaction, mediated by gluon exchange, which
is NOT contained in the rigid-string picture.  It is a SEPARATE physical
mechanism that the substrate model has not yet derived.

What the Y-junction slope SHOULD predict cleanly: orbital excitations
along a fixed-flavour, fixed-spin trajectory.  E.g.  the Δ trajectory
from Δ(1232, J=3/2) up to Δ(1950, J=7/2) is a pure ΔL=2 orbital
excitation, with no spin flip.  Likewise the N → N(1680) trajectory has
ΔL=2.  These are the right comparisons for a string-only theory.

This module therefore (i) reports the Y-junction slope, (ii) uses it to
predict the Δ(1232) → Δ(1950) and N(939) → N(1680) ΔJ=2 orbital steps,
and (iii) flags the residual at the N → Δ jump as the chromomagnetic
spin-spin physics that the substrate has not yet derived.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# Substrate input
# ---------------------------------------------------------------------------

GEV_TO_MEV = 1e3
SIGMA_QCD_GEV2 = 0.180  # from K(ξ_QCD) running

# Meson Regge slope
ALPHA_PRIME_MESON = 1.0 / (2.0 * math.pi * SIGMA_QCD_GEV2)   # ≈ 0.8842 GeV⁻²

# Y-junction Regge slope, derived above
ALPHA_PRIME_Y_JUNCTION = 1.0 / (3.0 * math.pi * SIGMA_QCD_GEV2)  # ≈ 0.5895 GeV⁻²
Y_OVER_MESON = ALPHA_PRIME_Y_JUNCTION / ALPHA_PRIME_MESON         # = 2/3 exactly


# ---------------------------------------------------------------------------
# Y-junction relativistic computation, as a sanity-check function
# ---------------------------------------------------------------------------

def y_junction_E_and_J(sigma_GeV2: float, L: float) -> tuple[float, float]:
    """Return total energy E and angular momentum J of a rotating Y-junction.

    Three rigid-rotor strings of length L, tension sigma_GeV2, ω = 1/L
    (endpoints at speed of light), rotating about the central vertex.

    Args:
        sigma_GeV2: Substrate string tension [GeV²].
        L:          Length of each arm in natural units.

    Returns:
        (E_GeV, J_natural) where J is the dimensionless angular momentum
        in natural units (it has dimensions of L*E in physical units).
    """
    E = 1.5 * math.pi * sigma_GeV2 * L          # (3π/2) σ L
    J = 0.75 * math.pi * sigma_GeV2 * L ** 2    # (3π/4) σ L²
    return E, J


def alpha_prime_from_Y(sigma_GeV2: float) -> float:
    """Return α'_Y derived from the Y-junction.

    Args:
        sigma_GeV2: Substrate string tension [GeV²].

    Returns:
        α'_Y = 1 / (3 π σ) in GeV⁻².
    """
    return 1.0 / (3.0 * math.pi * sigma_GeV2)


# ---------------------------------------------------------------------------
# PDG anchor data for baryon trajectories
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BaryonState:
    """One baryon resonance with its PDG mass and quantum numbers.

    Attributes:
        name:     PDG name.
        J:        Total angular momentum.
        L:        Orbital angular momentum (for spin-spin diagnostics).
        S:        Total quark spin.
        mass_pdg: Observed mass [MeV] from PDG.
    """
    name: str
    J: float
    L: int
    S: float
    mass_pdg: float


# Light-quark baryons. L,S noted so we can separate orbital vs spin-spin.
BARYON_STATES: list[BaryonState] = [
    BaryonState("N(939)",    J=0.5, L=0, S=0.5, mass_pdg=938.92),
    BaryonState("Δ(1232)",   J=1.5, L=0, S=1.5, mass_pdg=1232.0),
    BaryonState("N(1680)",   J=2.5, L=2, S=0.5, mass_pdg=1685.0),
    BaryonState("Δ(1950)",   J=3.5, L=2, S=1.5, mass_pdg=1930.0),
]


# ---------------------------------------------------------------------------
# Trajectory predictions under a hypothesised slope
# ---------------------------------------------------------------------------

@dataclass
class SlopeHypothesis:
    """A candidate Regge slope and the predictions it makes.

    Attributes:
        label:        Human-readable label for the slope.
        alpha_prime:  Slope α' [GeV⁻²].
        anchor_name:  Name of the baryon used to fix α₀.
        alpha_zero:   Trajectory intercept (computed from anchor).
        predictions:  {name: (M_pred_MeV, M_obs_MeV, frac_err)}.
    """
    label: str
    alpha_prime: float
    anchor_name: str
    alpha_zero: float = 0.0
    predictions: dict[str, tuple[float, float, float]] = field(default_factory=dict)


def regge_mass_GeV(J: float, alpha_prime: float, alpha_zero: float) -> float:
    """Compute mass on a Regge trajectory: M² = (J - α₀)/α'.

    Args:
        J:           Total angular momentum.
        alpha_prime: Slope α' [GeV⁻²].
        alpha_zero:  Intercept (dimensionless).

    Returns:
        Mass in GeV; 0 if (J - α₀) ≤ 0.
    """
    diff = J - alpha_zero
    if diff <= 0.0:
        return 0.0
    return math.sqrt(diff / alpha_prime)


def evaluate_slope(
    label: str,
    alpha_prime: float,
    states: list[BaryonState],
    anchor: BaryonState,
) -> SlopeHypothesis:
    """Predict baryon masses for a trial slope, anchored on one state.

    Args:
        label:       Label for the hypothesis.
        alpha_prime: Trial slope [GeV⁻²].
        states:      All baryon states to predict.
        anchor:      Baryon used to fix α₀.

    Returns:
        SlopeHypothesis with per-state (M_pred, M_obs, frac_err).
    """
    M_anchor = anchor.mass_pdg / GEV_TO_MEV
    alpha_zero = anchor.J - alpha_prime * M_anchor ** 2
    preds: dict[str, tuple[float, float, float]] = {}
    for s in states:
        M_pred_GeV = regge_mass_GeV(s.J, alpha_prime, alpha_zero)
        M_pred_MeV = M_pred_GeV * GEV_TO_MEV
        frac_err = (M_pred_MeV - s.mass_pdg) / s.mass_pdg
        preds[s.name] = (M_pred_MeV, s.mass_pdg, frac_err)
    return SlopeHypothesis(
        label=label,
        alpha_prime=alpha_prime,
        anchor_name=anchor.name,
        alpha_zero=alpha_zero,
        predictions=preds,
    )


# ---------------------------------------------------------------------------
# Spin-spin diagnostic
# ---------------------------------------------------------------------------

def spin_spin_residual(
    n_state: BaryonState, delta_state: BaryonState, alpha_prime: float
) -> dict[str, float]:
    """Diagnose how much of the N→Δ mass-squared splitting is non-string.

    Both states have L=0 in the constituent-quark picture; their J differs
    only by the symmetric-vs-antisymmetric SU(6) spin coupling.  The pure
    Regge contribution to ΔM² for a J change at fixed L would be (ΔJ)/α'.
    Anything beyond that is chromomagnetic.

    Args:
        n_state:     N(939).
        delta_state: Δ(1232).
        alpha_prime: Trial slope [GeV⁻²].

    Returns:
        Dict with the observed ΔM², the pure-string ΔM² for the same ΔJ,
        the residual ΔM², and the fraction of ΔM² that is non-string.
    """
    M_N = n_state.mass_pdg / GEV_TO_MEV
    M_D = delta_state.mass_pdg / GEV_TO_MEV
    dM2_obs = M_D ** 2 - M_N ** 2
    dJ = delta_state.J - n_state.J
    dM2_string = dJ / alpha_prime
    dM2_residual = dM2_obs - dM2_string
    fraction_non_string = dM2_residual / dM2_obs if dM2_obs != 0 else 0.0
    return {
        "dM2_obs_GeV2": dM2_obs,
        "dJ": dJ,
        "dM2_string_GeV2": dM2_string,
        "dM2_residual_GeV2": dM2_residual,
        "fraction_non_string": fraction_non_string,
    }


# ---------------------------------------------------------------------------
# Top-level summary
# ---------------------------------------------------------------------------

@dataclass
class YJunctionReport:
    """Aggregate result for the Y-junction baryon analysis."""
    sigma_GeV2: float
    alpha_prime_meson: float
    alpha_prime_y: float
    hypotheses: list[SlopeHypothesis]
    spin_spin: dict[str, float]
    delta_anchored_meson_slope: SlopeHypothesis
    delta_anchored_y_slope: SlopeHypothesis


def compute_report(sigma_GeV2: float = SIGMA_QCD_GEV2) -> YJunctionReport:
    """Build the full report comparing slope hypotheses for baryons.

    Args:
        sigma_GeV2: Substrate string tension [GeV²].

    Returns:
        Populated YJunctionReport.
    """
    states = BARYON_STATES
    nucleon = next(s for s in states if s.name == "N(939)")
    delta = next(s for s in states if s.name == "Δ(1232)")

    a_meson = 1.0 / (2.0 * math.pi * sigma_GeV2)
    a_y = 1.0 / (3.0 * math.pi * sigma_GeV2)
    a_three_halves = a_meson * 1.5

    # All anchored at the nucleon (the standard convention).
    hypotheses = [
        evaluate_slope(
            "α'_meson  = 1/(2πσ)        (current implementation)",
            a_meson, states, nucleon,
        ),
        evaluate_slope(
            "α'_Y       = (2/3)·α'_meson  (rigid Y-junction derivation)",
            a_y, states, nucleon,
        ),
        evaluate_slope(
            "α' × 3/2   = 1.5·α'_meson  (hypothesised inverted geometry)",
            a_three_halves, states, nucleon,
        ),
    ]

    # Spin-spin diagnostic with the meson slope (the standard reference).
    spin_spin = spin_spin_residual(nucleon, delta, a_meson)

    # Anchor on the Δ instead of the N — isolates the orbital trajectory
    # by side-stepping the chromomagnetic spin-spin shift in the anchor.
    delta_anchored_meson = evaluate_slope(
        "α'_meson, anchored at Δ(1232)",
        a_meson, states, delta,
    )
    delta_anchored_y = evaluate_slope(
        "α'_Y, anchored at Δ(1232)",
        a_y, states, delta,
    )

    return YJunctionReport(
        sigma_GeV2=sigma_GeV2,
        alpha_prime_meson=a_meson,
        alpha_prime_y=a_y,
        hypotheses=hypotheses,
        spin_spin=spin_spin,
        delta_anchored_meson_slope=delta_anchored_meson,
        delta_anchored_y_slope=delta_anchored_y,
    )


def print_report(report: YJunctionReport) -> None:
    """Print the Y-junction analysis as a clean text table.

    Args:
        report: Result of compute_report.
    """
    print("=" * 78)
    print("Y-JUNCTION BARYON REGGE SLOPE ANALYSIS")
    print("=" * 78)
    print(f"  σ        = {report.sigma_GeV2:.4f} GeV²")
    print(f"  α'_meson = {report.alpha_prime_meson:.4f} GeV⁻²  = 1/(2πσ)")
    print(f"  α'_Y     = {report.alpha_prime_y:.4f} GeV⁻²  = 1/(3πσ) = (2/3)·α'_meson")
    print()
    print("  Derivation: rigid Y-junction of three strings of length L,")
    print("  rotating with ω=1/L so the outer ends move at c.  Integrating the")
    print("  three relativistic arms gives E = (3π/2)σL and J = (3π/4)σL²,")
    print("  hence α'_Y = J/M² = 1/(3πσ).")
    print()

    print("-" * 78)
    print("  CANDIDATE SLOPES, anchored at N(939) (J = 1/2)")
    print("-" * 78)
    for h in report.hypotheses:
        print(f"\n  {h.label}")
        print(f"    α'    = {h.alpha_prime:.4f} GeV⁻²")
        print(f"    α₀    = {h.alpha_zero:+.4f}")
        print(f"    {'state':<14}{'M_pred [MeV]':>15}{'M_PDG [MeV]':>15}{'err':>10}")
        for name, (Mp, Mo, err) in h.predictions.items():
            tag = "  ← anchor" if name == h.anchor_name else ""
            if Mp > 0.0:
                print(f"    {name:<14}{Mp:>15.4g}{Mo:>15.4g}{100*err:>+9.2f}%{tag}")
            else:
                print(f"    {name:<14}{'—':>15}{Mo:>15.4g}{'(below)':>10}")

    print()
    print("-" * 78)
    print("  SPIN-SPIN DIAGNOSTIC: how much of N→Δ is NOT string physics?")
    print("-" * 78)
    ss = report.spin_spin
    print(f"    observed ΔM²        = {ss['dM2_obs_GeV2']:.4f} GeV²")
    print(f"    string ΔM² for ΔJ=1 = {ss['dM2_string_GeV2']:.4f} GeV²"
          f"   (using α'_meson)")
    print(f"    residual ΔM²        = {ss['dM2_residual_GeV2']:.4f} GeV²")
    print(f"    non-string fraction = {100.0 * ss['fraction_non_string']:.2f}%")
    print()
    print("  Interpretation: in the constituent-quark picture both N(939) and")
    print("  Δ(1232) have L=0 — they differ only in quark-spin coupling, so the")
    print("  ENTIRE ΔM² should be chromomagnetic spin-spin physics, NOT Regge.")
    print("  Forcing an orbital trajectory through them therefore mis-attributes")
    print("  spin-spin energy to the slope.")
    print()

    print("-" * 78)
    print("  CLEAN ORBITAL TRAJECTORY: anchor on Δ(1232), test Δ(1950) only")
    print("-" * 78)
    for h in (report.delta_anchored_meson_slope, report.delta_anchored_y_slope):
        print(f"\n  {h.label}")
        print(f"    α'    = {h.alpha_prime:.4f} GeV⁻²")
        print(f"    α₀    = {h.alpha_zero:+.4f}")
        print(f"    {'state':<14}{'M_pred [MeV]':>15}{'M_PDG [MeV]':>15}{'err':>10}")
        for name, (Mp, Mo, err) in h.predictions.items():
            tag = "  ← anchor" if name == h.anchor_name else ""
            if Mp > 0.0:
                print(f"    {name:<14}{Mp:>15.4g}{Mo:>15.4g}{100*err:>+9.2f}%{tag}")
            else:
                print(f"    {name:<14}{'—':>15}{Mo:>15.4g}{'(below)':>10}")
    print()
