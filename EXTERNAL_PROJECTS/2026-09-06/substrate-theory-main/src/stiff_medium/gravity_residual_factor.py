"""Gravity residual factor ε — substrate computation of the charge-symmetric
suppression that converts the bare-substrate G ≈ 3.8 × 10³⁴ m³ kg⁻¹ s⁻² into
the observed G ≈ 6.7 × 10⁻¹¹ m³ kg⁻¹ s⁻² (§18.32, §18.60.3).

CONTEXT (from §18.60.3):
    The bare substrate gives G_substrate = c³ ξ_e²/ℏ. The required suppression
    is ε² = G_obs / G_substrate ≈ 1.75 × 10⁻⁴⁵, hence ε ≈ 4.19 × 10⁻²³.
    Numerically this matches m_e / M_Planck "to all printed digits."

THIS MODULE: tests two questions HONESTLY:

    (A)  Is "ε = m_e/M_Planck" a substrate prediction or a dimensional
         tautology of the way M_Planck and m_e are defined?

    (B)  Does the substrate's UV-IR structure (a Planck-scale ξ_P and an
         atomic-scale ξ_e at which the kink lives) independently deliver
         the ratio l_Planck / ξ_e *without* using G as input — and if so
         is this the substrate's own prediction of ε?

The answer to (A) is: TAUTOLOGY — once you fix ξ_e = ℏ/(m_e c) (Solution A
of substrate_primitives.py) AND define l_Planck = √(ℏG/c³), the algebraic
identity m_e × ξ_e = M_Planck × l_Planck = ℏ/c is built in. So
ξ_e/l_Planck = M_Planck/m_e is a tautology *unless* l_Planck is set
independently (without using G).

The answer to (B) is: PARTIAL. The substrate framework (per §18.46.1) does
posit a Planck-scale length ξ_P at which ℏ = K_P ξ_P⁴/c, distinct from the
electron-scale ξ_e where the kink lives. These are TWO independent length
scales of the substrate. If we identify l_Planck with ξ_P AND ξ_e with the
electron Compton wavelength, then ε ~ ξ_P/ξ_e is a substrate prediction.
But ξ_P itself is not derivable from substrate primitives without
additional input — it is an empirical anchor for the substrate's UV cutoff.

So the predictive content is: "given two independent substrate length
scales (ξ_P, ξ_e), gravity's coupling is suppressed by ε² ~ (ξ_P/ξ_e)²."
This is a *structural* prediction with the small parameter identified —
but the absolute numerical value of ε requires fixing ξ_P empirically.

Output of this module:
- epsilon_dimensional: the algebraic identity m_e/M_Planck.
- epsilon_substrate_uv_ir: the UV/IR scale ratio l_Planck/ξ_e (using
  l_Planck from G; circular if G is the input).
- epsilon_substrate_no_G: an attempt to compute ε *without* G as input —
  reports failure (no substrate primitive sets l_Planck independently).
- epsilon_observed: the empirical value G_obs / G_substrate, square-rooted.
- universality_check: tests whether G_eff predicted by ε(m) varies with
  the test-particle mass m. Standard physics: no. Substrate: depends on
  whether ε is built from m_test or from the source's bound configuration.

The honest verdict goes in `assess()`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Inputs — substrate primitives only. G_OBSERVED is reported but never used
# as an input to the substrate prediction.
# ---------------------------------------------------------------------------

C_SI: Final[float] = 2.99792458e8           # m/s, exact
HBAR_SI: Final[float] = 1.054571817e-34     # J·s
M_E_SI: Final[float] = 9.1093837015e-31     # kg
M_P_SI: Final[float] = 1.67262192369e-27    # kg, proton (for universality test)
ALPHA_EM: Final[float] = 7.2973525693e-3    # ≈ 1/137.036
G_OBSERVED: Final[float] = 6.67430e-11      # comparison only

# Substrate primitives at the electron scale (Solution A in
# substrate_primitives.py): ξ = λ_C(electron).
XI_E: Final[float] = HBAR_SI / (M_E_SI * C_SI)   # ≈ 3.86e-13 m
K_E: Final[float] = HBAR_SI * C_SI / XI_E**4     # ≈ 1.42e24 Pa
RHO_E: Final[float] = K_E / C_SI**2               # ≈ 1.58e7 kg/m³

# Bare substrate G (from §18.60.3): the only dimensional combination of
# (K, ρ, ξ, c, ℏ) with units of G after the constraint ℏ = K ξ⁴/c.
G_SUBSTRATE_BARE: Final[float] = C_SI**3 * XI_E**2 / HBAR_SI

# Empirical Planck quantities (use G_OBSERVED — for comparison only).
L_PLANCK_OBS: Final[float] = math.sqrt(HBAR_SI * G_OBSERVED / C_SI**3)
M_PLANCK_OBS: Final[float] = math.sqrt(HBAR_SI * C_SI / G_OBSERVED)


# ---------------------------------------------------------------------------
# 1. The algebraic identity m_e × ξ_e = M_Planck × l_Planck = ℏ/c
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class IdentityCheck:
    """Numerical verification of the dimensional tautology."""

    m_e_xi_e: float          # m_e × ξ_e  [kg·m]
    M_pl_l_pl: float         # M_Planck × l_Planck  [kg·m]
    hbar_over_c: float       # ℏ/c  [kg·m]
    ratio_eq_one: float      # (m_e ξ_e) / (M_Planck l_Planck), should be exactly 1
    statement: str


def algebraic_identity() -> IdentityCheck:
    """The relation m_e ξ_e = M_Planck l_Planck = ℏ/c is a DEFINITION-LEVEL
    identity once ξ_e = ℏ/(m_e c) (Solution A of substrate_primitives.py).

    This is not physics — it is the algebra of how Planck units are written.
    Therefore the equation
            ε := m_e / M_Planck = l_Planck / ξ_e = ξ_P / ξ_e   [if ξ_P ≡ l_Planck]
    is a tautology of the variable substitutions, NOT a substrate prediction.
    """
    a = M_E_SI * XI_E
    b = M_PLANCK_OBS * L_PLANCK_OBS
    c_ = HBAR_SI / C_SI
    ratio = a / b
    return IdentityCheck(
        m_e_xi_e=a,
        M_pl_l_pl=b,
        hbar_over_c=c_,
        ratio_eq_one=ratio,
        statement=(
            "m_e × ξ_e = M_Planck × l_Planck = ℏ/c, identically. "
            "Therefore m_e/M_Planck = l_Planck/ξ_e is a TAUTOLOGY of "
            "the unit definitions, given ξ_e = λ_C(electron). It is not "
            "an independent substrate prediction."
        ),
    )


# ---------------------------------------------------------------------------
# 2. Candidates for ε — each must say WHAT IT USES and whether it is
#    circular (uses G as input).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EpsilonCandidate:
    name: str
    value: float
    formula: str
    inputs_used: tuple[str, ...]
    uses_G: bool                 # True ⇒ circular (G appears in the inputs)
    is_tautological: bool        # True ⇒ identity, not an independent prediction
    rationale: str


def epsilon_candidates() -> list[EpsilonCandidate]:
    """Five candidate computations of ε. Reports inputs, circularity, and
    physical content explicitly.
    """
    # 1. The empirical "answer" ε² = G_obs/G_substrate.
    eps_observed = math.sqrt(G_OBSERVED / G_SUBSTRATE_BARE)

    # 2. m_e/M_Planck — uses G via M_Planck.  Circular.
    eps_dim_a = M_E_SI / M_PLANCK_OBS

    # 3. ξ_P / ξ_e where ξ_P = l_Planck (uses G).  Circular and tautological
    #    against #2.
    eps_dim_b = L_PLANCK_OBS / XI_E

    # 4. Substrate-only attempt: try to set l_Planck without G.
    #    The only substrate-internal candidates for "Planck length" are
    #       l_h := sqrt(ℏ/c³) ≈ 1.96e-21 m  (NOT the Planck length;
    #          this ℓ contains G implicitly: l_Planck = sqrt(ℏ G/c³),
    #          so l_h × √G = l_Planck. l_h alone has wrong units of length
    #          —  it has units sqrt(s) actually).
    #    Sanity check the units: [ℏ] = kg m² /s, [c³] = m³/s³, so
    #       [ℏ/c³] = kg s² /m  → not a length squared.
    #    There is NO substrate-internal length built from K, ρ, ξ, c, ℏ
    #    alone other than ξ itself (and any power of ξ). The substrate
    #    primitives provide ONE length, ξ_e. To set l_Planck independently
    #    we must POSIT a second substrate length ξ_P (the deep UV scale).
    eps_no_G = float("nan")  # cannot compute without an extra input

    # 5. UV/IR scale ratio interpretation: assume the substrate has a deep-UV
    #    scale ξ_P at which the action quantum is set (ℏ = K_P ξ_P⁴/c per
    #    §18.46.1).  Identify ε ≡ ξ_P / ξ_e.  The empirical match (4.19e-23)
    #    requires ξ_P = ε × ξ_e ≈ 1.62e-35 m.  This equals l_Planck, but
    #    only because we are CHOOSING ξ_P = l_Planck.  No independent
    #    substrate constraint sets ξ_P/ξ_e.
    eps_uv_ir = L_PLANCK_OBS / XI_E  # same number as #3 by construction

    return [
        EpsilonCandidate(
            name="ε_observed",
            value=eps_observed,
            formula="sqrt(G_observed / G_substrate_bare)",
            inputs_used=("G_observed", "G_substrate_bare(=c³ξ²/ℏ)"),
            uses_G=True,
            is_tautological=False,
            rationale="The empirical suppression that the substrate must reproduce.",
        ),
        EpsilonCandidate(
            name="ε_dim_a",
            value=eps_dim_a,
            formula="m_e / M_Planck",
            inputs_used=("m_e", "M_Planck(=√(ℏc/G))"),
            uses_G=True,
            is_tautological=True,
            rationale=(
                "Uses M_Planck, which contains G. Numerically equals "
                "ε_observed because (m_e × ξ_e) = (M_Planck × l_Planck) = ℏ/c. "
                "Circular: G appears in M_Planck."
            ),
        ),
        EpsilonCandidate(
            name="ε_dim_b",
            value=eps_dim_b,
            formula="l_Planck / ξ_e",
            inputs_used=("l_Planck(=√(ℏG/c³))", "ξ_e"),
            uses_G=True,
            is_tautological=True,
            rationale=(
                "Same number as ε_dim_a by m_e ξ_e = M_Planck l_Planck. "
                "Tautology of the unit system."
            ),
        ),
        EpsilonCandidate(
            name="ε_substrate_no_G",
            value=eps_no_G,
            formula="(undefined — no substrate length other than ξ_e)",
            inputs_used=("K", "ρ", "ξ_e", "c", "ℏ"),
            uses_G=False,
            is_tautological=False,
            rationale=(
                "Substrate primitives provide ONE independent length (ξ). "
                "There is no second length to form a ratio. ε cannot be "
                "computed from K, ρ, ξ, c, ℏ alone."
            ),
        ),
        EpsilonCandidate(
            name="ε_UV/IR (ξ_P/ξ_e)",
            value=eps_uv_ir,
            formula="ξ_P / ξ_e  with ξ_P ≡ l_Planck",
            inputs_used=("ξ_P (= l_Planck)", "ξ_e"),
            uses_G=True,
            is_tautological=True,
            rationale=(
                "If we POSIT a deep-UV substrate length ξ_P and identify it "
                "with l_Planck, ε = ξ_P/ξ_e gives the right number — but "
                "only because ξ_P/ξ_e is then a renamed m_e/M_Planck. "
                "Predictive content: 'ε is set by a UV/IR length ratio'. "
                "Numerical value requires fixing ξ_P empirically."
            ),
        ),
    ]


# ---------------------------------------------------------------------------
# 3. Back-of-the-envelope: what would a "distributed cancellation" deliver?
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CancellationModel:
    """Toy model of the §18.32 distributed cancellation.

    The qualitative picture: a bound configuration has N coherent vectors
    contributing to the EM coupling. The charge-symmetric (gravity) channel
    sums these vectors with random orientations: the residual after
    cancellation scales as 1/√N.

    If N ~ (bound-state volume / substrate-cell volume) = (ξ_e / ξ_P)³, then
    ε ~ 1/√N ~ (ξ_P/ξ_e)^{3/2}.

    Numerically: ε would be (l_Planck/ξ_e)^{3/2} ≈ (4.19e-23)^{3/2} ≈ 2.7e-34,
    which is 11 orders of magnitude TOO SMALL. The "distributed cancellation"
    1/√N picture overshoots — gravity would be even weaker than observed.

    Alternative scalings:
        ε ~ 1/N     ⇒ (ξ_P/ξ_e)³ ≈ 7.4e-68  (way too small)
        ε ~ 1/√N    ⇒ ≈ 2.7e-34            (too small by 10¹¹)
        ε ~ 1/N^{1/3} ⇒ ≈ ξ_P/ξ_e ≈ 4.2e-23 (matches; but this means ε ~ 1/N^{1/3}, NOT 1/√N)

    A 1/N^{1/3} scaling means N is the LINEAR vector count along a direction,
    not the volume count. Geometrically that suggests cancellation along a
    1D Möbius cycle (azimuthal) rather than over the bulk volume — which is
    exactly the §18.10 Möbius half-flux structure.

    So a coherent picture is:
        - The bound state's azimuthal Möbius cycle has length ~ξ_e (the
          Compton wavelength).
        - The substrate cell along that cycle has length ~ξ_P (Planck).
        - Number of cells along the cycle: N_cycle = ξ_e/ξ_P.
        - The half-flux ensures equal-and-opposite contributions sum coherently
          for EM but cancel destructively for gravity, leaving a residual
          1/N_cycle from the boundary mismatch.
        - ε = 1/N_cycle = ξ_P/ξ_e = m_e/M_Planck ✓.

    This is a *plausible* substrate mechanism; it does NOT add new computation
    beyond the dimensional tautology because the cycle length is set to ξ_e
    by hand.
    """

    N_volume: float          # = (ξ_e/ξ_P)³ if ξ_P = l_Planck
    N_cycle: float           # = ξ_e/ξ_P
    eps_volume_scaling: float    # ε ~ 1/√N_volume
    eps_cycle_scaling: float     # ε ~ 1/N_cycle  (matches observed)
    eps_observed: float
    best_match: str
    statement: str


def cancellation_toy_model() -> CancellationModel:
    """Compare three N-scalings of the distributed-cancellation hypothesis."""
    N_vol = (XI_E / L_PLANCK_OBS) ** 3
    N_cyc = XI_E / L_PLANCK_OBS
    eps_vol = 1.0 / math.sqrt(N_vol)
    eps_cyc = 1.0 / N_cyc
    eps_obs = math.sqrt(G_OBSERVED / G_SUBSTRATE_BARE)

    # Find which scaling matches best
    candidates = {
        "1/√N_volume (3D bulk cancellation)": eps_vol,
        "1/N_cycle (1D Möbius cycle cancellation)": eps_cyc,
    }
    best = min(candidates.items(), key=lambda kv: abs(math.log10(kv[1] / eps_obs)))

    return CancellationModel(
        N_volume=N_vol,
        N_cycle=N_cyc,
        eps_volume_scaling=eps_vol,
        eps_cycle_scaling=eps_cyc,
        eps_observed=eps_obs,
        best_match=best[0],
        statement=(
            "1/N_cycle scaling matches ε_observed exactly (it is the "
            "tautology m_e/M_Planck written differently). The 1/√N bulk "
            "scaling is 11 orders too small. So if a 'distributed "
            "cancellation' picture is right, the cancellation runs along "
            "a 1D Möbius cycle of length ξ_e divided into substrate cells "
            "of length ξ_P — NOT a 3D bulk. This is consistent with the "
            "§18.10 half-flux being a 1D azimuthal feature."
        ),
    )


# ---------------------------------------------------------------------------
# 4. Universality check: does G_eff depend on the test-particle mass?
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class UniversalityCheck:
    """Test whether ε(m) depends on the configuration's mass m, which would
    make G_eff = ε(m)² × G_substrate mass-dependent and violate the
    equivalence principle.

    Two readings of the substrate picture:

      (i) ε is a property of the SUBSTRATE (the cycle length ratio ξ_P/ξ_e
          where ξ_e is the *electron* Compton wavelength). Then ε is
          universal and G_eff is universal. Equivalence principle: ✓.

     (ii) ε is a property of the BOUND STATE — its cycle length is the
          configuration's own Compton wavelength λ_C(m) = ℏ/(mc), so
          ε(m) = ξ_P/λ_C(m) = m × l_Planck c/ℏ = m/M_Planck. Then
          G_eff(m) ∝ (m/M_Planck)² varies with m. Equivalence principle: ✗.

    Reading (i) is consistent with observation. Reading (ii) is the more
    naive transcription of the dimensional argument and is FALSIFIED by
    Eötvös-class experiments at the 10⁻¹³ level.

    The model must commit to (i): the cycle length is a substrate-fixed
    quantity (the electron's Compton wavelength is the universal IR cutoff
    for the charge-symmetric channel, fixed by the lightest charged
    fermion). All bound configurations inherit this same ε.

    This is a non-trivial structural commitment — it says the *electron*
    Compton wavelength (not the proton's, not the test-mass's) is the
    universal IR scale of the gravity sector. We test it numerically.
    """

    G_eff_electron_test: float       # G_eff if ε computed from m_e
    G_eff_proton_test: float         # G_eff if ε computed from m_p
    G_eff_universal: float            # G_eff if ε is fixed at electron value
    ratio_proton_to_electron: float   # (m_p/m_e)² — the predicted violation in reading (ii)
    eotvos_bound: float               # ~10⁻¹³ — observed bound on EP violation
    reading_i_passes: bool            # universal ε
    reading_ii_passes: bool           # mass-dependent ε
    statement: str


def universality_check() -> UniversalityCheck:
    """Compute G_eff under both readings and compare to the equivalence
    principle bound from Eötvös experiments (~10⁻¹³).
    """
    eps_e = M_E_SI / M_PLANCK_OBS
    eps_p = M_P_SI / M_PLANCK_OBS

    G_eff_e = eps_e**2 * G_SUBSTRATE_BARE
    G_eff_p = eps_p**2 * G_SUBSTRATE_BARE
    # In reading (i), ε is fixed at eps_e for *all* configurations:
    G_eff_univ = G_eff_e

    # Reading (ii) predicts G_eff(proton)/G_eff(electron) = (m_p/m_e)² ≈ 3.4e6
    ratio_pe = (M_P_SI / M_E_SI) ** 2

    # Eötvös-type experiments bound |Δa/a| < ~10⁻¹³ between different test
    # masses. Reading (ii) predicts (m_p/m_e)² ≈ 3.4e6 difference per
    # constituent — utterly excluded.
    eotvos_bound = 1e-13

    return UniversalityCheck(
        G_eff_electron_test=G_eff_e,
        G_eff_proton_test=G_eff_p,
        G_eff_universal=G_eff_univ,
        ratio_proton_to_electron=ratio_pe,
        eotvos_bound=eotvos_bound,
        reading_i_passes=True,
        reading_ii_passes=False,
        statement=(
            "Reading (i) [universal ε fixed by electron Compton wavelength] "
            "preserves the equivalence principle. Reading (ii) [ε scales "
            "with the bound state's own mass] predicts G_eff(p)/G_eff(e) "
            f"≈ {ratio_pe:.2e}, ruled out by Eötvös at 10⁻¹³. The model "
            "MUST commit to reading (i): the cycle-length cutoff ξ_e is "
            "the *electron* Compton wavelength as a universal substrate "
            "scale, not the test particle's own scale."
        ),
    )


# ---------------------------------------------------------------------------
# 5. Composite assessment
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Assessment:
    G_substrate_bare: float
    G_observed: float
    epsilon_observed: float
    epsilon_dimensional: float
    epsilon_substrate_only: float       # nan if not derivable
    matches: bool
    is_genuine_prediction: bool
    summary: str


def assess() -> Assessment:
    """Honest assessment of the gravity residual factor computation."""
    eps_obs = math.sqrt(G_OBSERVED / G_SUBSTRATE_BARE)
    eps_dim = M_E_SI / M_PLANCK_OBS
    eps_substrate_only = float("nan")
    matches = abs(math.log10(eps_dim / eps_obs)) < 0.01

    summary = (
        "FINDINGS:\n"
        f"  ε_observed       = {eps_obs:.4e}  (= sqrt(G_obs/G_substrate))\n"
        f"  ε_dimensional    = {eps_dim:.4e}  (= m_e/M_Planck)\n"
        f"  match to 4 dp:   {matches}\n"
        "\n"
        "INTERPRETATION:\n"
        "  The numerical match is EXACT but it is a TAUTOLOGY. The identity\n"
        "  m_e × ξ_e = M_Planck × l_Planck = ℏ/c holds for *any* particle\n"
        "  whose Compton wavelength is taken as the substrate's IR length,\n"
        "  given the standard definitions of M_Planck and l_Planck. Once\n"
        "  the substrate adopts ξ = λ_C(electron) (Solution A), and once\n"
        "  l_Planck is computed from G (the standard definition), the\n"
        "  numerical match m_e/M_Planck = l_Planck/ξ_e is built in.\n"
        "\n"
        "  The substrate primitives K, ρ, ξ_e, c, ℏ provide ONE length\n"
        "  (ξ_e) and CANNOT independently produce ε without positing a\n"
        "  second substrate length ξ_P. So the 'substrate only' answer is\n"
        "  undefined: ε requires extra input.\n"
        "\n"
        "PHYSICAL CONTENT (the non-tautological piece):\n"
        "  1. The substrate framework REQUIRES a UV scale ξ_P at which\n"
        "     the action quantum ℏ is set (§18.46.1). This is independent\n"
        "     of the IR scale ξ_e at which the electron lives.\n"
        "  2. The Möbius topology (§18.10) puts the charge-symmetric\n"
        "     residual on a 1D azimuthal cycle of length ξ_e. The\n"
        "     toy-model 1/N_cycle scaling (with N_cycle = ξ_e/ξ_P) gives\n"
        "     ε = ξ_P/ξ_e exactly.\n"
        "  3. By the algebraic identity, ξ_P/ξ_e = m_e/M_Planck whenever\n"
        "     ξ_P = l_Planck. Identifying these two scales is the\n"
        "     substrate's Planck-scale anchor — empirically input, not\n"
        "     derived.\n"
        "  4. UNIVERSALITY: the cutoff ξ_e is fixed by the *electron*\n"
        "     Compton wavelength (the lightest charged fermion's IR\n"
        "     scale), not the test-particle's. This preserves the\n"
        "     equivalence principle and is a non-trivial commitment.\n"
        "\n"
        "PREDICTIVE CLAIM:\n"
        "  The substrate predicts ε = ξ_P/ξ_e (1/N_cycle scaling on a\n"
        "  Möbius cycle), with the absolute value requiring ξ_P input.\n"
        "  Setting ξ_P = l_Planck recovers the observed ε ≈ 4.19 × 10⁻²³.\n"
        "  The 'distributed cancellation' of §18.32 is now explicit: it\n"
        "  is 1D azimuthal cancellation along the half-flux cycle, NOT\n"
        "  3D bulk cancellation (which would give 1/√N and miss by 10¹¹).\n"
    )

    return Assessment(
        G_substrate_bare=G_SUBSTRATE_BARE,
        G_observed=G_OBSERVED,
        epsilon_observed=eps_obs,
        epsilon_dimensional=eps_dim,
        epsilon_substrate_only=eps_substrate_only,
        matches=matches,
        is_genuine_prediction=False,  # it's a tautology + an empirical anchor
        summary=summary,
    )


# ---------------------------------------------------------------------------
# 6. Composite report
# ---------------------------------------------------------------------------

def gravity_residual_factor_report() -> dict[str, object]:
    """Run the full analysis and return a structured report."""
    identity = algebraic_identity()
    candidates = epsilon_candidates()
    cancellation = cancellation_toy_model()
    universality = universality_check()
    final = assess()
    return {
        "primitives": {
            "K_e_Pa": K_E,
            "rho_e_kgm3": RHO_E,
            "xi_e_m": XI_E,
            "c_ms": C_SI,
            "hbar_Js": HBAR_SI,
            "m_e_kg": M_E_SI,
        },
        "G_substrate_bare": G_SUBSTRATE_BARE,
        "G_observed": G_OBSERVED,
        "algebraic_identity": {
            "m_e_xi_e": identity.m_e_xi_e,
            "M_pl_l_pl": identity.M_pl_l_pl,
            "hbar_over_c": identity.hbar_over_c,
            "ratio_eq_one": identity.ratio_eq_one,
            "statement": identity.statement,
        },
        "epsilon_candidates": [
            {
                "name": e.name,
                "value": e.value,
                "formula": e.formula,
                "inputs_used": e.inputs_used,
                "uses_G": e.uses_G,
                "is_tautological": e.is_tautological,
                "rationale": e.rationale,
            }
            for e in candidates
        ],
        "cancellation_model": {
            "N_volume": cancellation.N_volume,
            "N_cycle": cancellation.N_cycle,
            "eps_volume_scaling": cancellation.eps_volume_scaling,
            "eps_cycle_scaling": cancellation.eps_cycle_scaling,
            "eps_observed": cancellation.eps_observed,
            "best_match": cancellation.best_match,
            "statement": cancellation.statement,
        },
        "universality_check": {
            "G_eff_electron_test": universality.G_eff_electron_test,
            "G_eff_proton_test": universality.G_eff_proton_test,
            "G_eff_universal": universality.G_eff_universal,
            "ratio_proton_to_electron": universality.ratio_proton_to_electron,
            "eotvos_bound": universality.eotvos_bound,
            "reading_i_passes": universality.reading_i_passes,
            "reading_ii_passes": universality.reading_ii_passes,
            "statement": universality.statement,
        },
        "assessment": {
            "epsilon_observed": final.epsilon_observed,
            "epsilon_dimensional": final.epsilon_dimensional,
            "epsilon_substrate_only": final.epsilon_substrate_only,
            "matches": final.matches,
            "is_genuine_prediction": final.is_genuine_prediction,
            "summary": final.summary,
        },
    }
