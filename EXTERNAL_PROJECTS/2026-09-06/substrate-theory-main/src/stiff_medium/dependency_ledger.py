"""Dependency ledger for stiff-medium model claims.

The ledger is deliberately conservative.  It separates claims that are
directly derived from substrate assumptions from claims that are inherited
from standard effective theories, use anchors, are calibrated, or have failed.

This is not a physics engine.  It is bookkeeping that keeps the model honest.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from typing import Iterable


class ClaimTag(StrEnum):
    """Evidence/dependency status for a model claim."""

    PRIMITIVE = "primitive"
    STRUCTURAL = "structural"
    DERIVED = "derived"
    INHERITED = "inherited"
    ANCHORED = "anchored"
    CALIBRATED = "calibrated"
    FAILED = "failed"
    OPEN = "open"


class Risk(StrEnum):
    """Residual model risk for a ledger entry."""

    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    BLOCKER = "blocker"


@dataclass(frozen=True)
class LedgerEntry:
    """One model claim and its dependency status."""

    claim_id: str
    sector: str
    claim: str
    tag: ClaimTag
    risk: Risk
    evidence: str
    dependencies: tuple[str, ...] = ()
    next_step: str = ""

    @property
    def counts_as_independent_prediction(self) -> bool:
        """True only for direct derived substrate results."""
        return self.tag == ClaimTag.DERIVED

    @property
    def needs_work(self) -> bool:
        """True for entries that should stay on the active work queue."""
        return self.tag in {
            ClaimTag.PRIMITIVE,
            ClaimTag.ANCHORED,
            ClaimTag.CALIBRATED,
            ClaimTag.FAILED,
            ClaimTag.OPEN,
        } or self.risk in {Risk.HIGH, Risk.BLOCKER}


LEDGER: tuple[LedgerEntry, ...] = (
    LedgerEntry(
        "P001",
        "substrate",
        "Core medium primitives K, rho, xi",
        ClaimTag.PRIMITIVE,
        Risk.BLOCKER,
        "Required inputs; not derived from deeper substrate dynamics.",
        next_step="Derive primitives or declare minimal axiom set.",
    ),
    LedgerEntry(
        "P002",
        "substrate",
        "45 degree cone propagation constraint",
        ClaimTag.PRIMITIVE,
        Risk.HIGH,
        "Geometric compaction replaces lambda bookkeeping with Q-null cone form; "
        "equal-partition elastic penalty selects 45 degrees; lattice-invariant audit "
        "shows the quartic is first only if self-dual longitudinal/transverse "
        "exchange removes the allowed quadratic bias. Paired dual-branch exchange "
        "cancels the bias and gives beta>0 only when branch weights are exactly equal. "
        "Swap-symmetric local detailed balance gives exact 50/50 weights if the "
        "dual exchange generator is degenerate. A branch-swap elastic-cell "
        "automorphism J^T H J = H is sufficient to force that generator; a "
        "saturated diamond spring cell supplies that automorphism conditionally. "
        "Graph enumeration shows the diamond is uniquely minimal only when two "
        "saturated anchors and direct branch exchange are required. "
        "Schur-complement elimination of shared finite-compliance saturated "
        "anchors induces that direct L-T exchange, so it need not be a separate "
        "primitive. A neutral saturated phase-slip segment conditionally selects "
        "paired endpoint anchors, and the saturation barrier supplies finite "
        "anchor compliance below the exact cap. A discrete lattice boundary "
        "realizes that segment as the minimal nonzero saturated 1-chain; the "
        "cone angle is robust for any finite symmetric anchor stiffness ratio, "
        "but topology does not fix the ratio. A pure saturation-barrier energy "
        "delocalizes an imposed phase slip instead of selecting one bond; a "
        "Peierls/core localization cost or equivalent loaded saddle is still "
        "needed.",
        dependencies=("P001",),
        next_step=(
            "Derive the Peierls/core localization term, anchor/branch stiffness "
            "ratio, or loaded saturated-bond saddle from the actual substrate "
            "stiffness tensor or saturation microstructure."
        ),
    ),
    LedgerEntry(
        "P003",
        "topology",
        "Mobius half-flux topology for spin and Pauli structure",
        ClaimTag.STRUCTURAL,
        Risk.MEDIUM,
        "Implemented dynamically; geometric compaction packages it as Mobius connection holonomy.",
        dependencies=("P001", "P002"),
        next_step="Derive half-flux as lowest-energy substrate closure.",
    ),
    LedgerEntry(
        "P004",
        "substrate",
        "Saturation cap sigma_max = 1/2",
        ClaimTag.PRIMITIVE,
        Risk.HIGH,
        "Central to black holes and cosmology; value is committed, not derived.",
        dependencies=("P001",),
        next_step="Derive from elastic stability or deeper lattice geometry.",
    ),
    LedgerEntry(
        "D001",
        "kinematics",
        "Wave speed c = sqrt(K/rho)",
        ClaimTag.DERIVED,
        Risk.LOW,
        "Direct linear-elastic wave result.",
        dependencies=("P001",),
    ),
    LedgerEntry(
        "D002",
        "kinematics",
        "Planck action identity hbar = K xi^4 / c",
        ClaimTag.DERIVED,
        Risk.MEDIUM,
        "Dimensional substrate identity used throughout; depends on xi primitive.",
        dependencies=("P001", "D001"),
        next_step="Clarify relation between K_action and coarse-grained K.",
    ),
    LedgerEntry(
        "D003",
        "kinematics",
        "Mass-energy identity E = m c^2 from internal vector energy",
        ClaimTag.DERIVED,
        Risk.LOW,
        "Follows from cone/internal energy bookkeeping in current model.",
        dependencies=("P001", "P002", "D001"),
    ),
    LedgerEntry(
        "D004",
        "atomic",
        "Hydrogen isotope shifts from back-reaction and measured mass ratios",
        ClaimTag.DERIVED,
        Risk.MEDIUM,
        "Sub-ppm agreement in scripts; uses real isotope masses, not fitted shifts.",
        dependencies=("P001", "P002"),
    ),
    LedgerEntry(
        "I001",
        "atomic",
        "Hydrogen/helium/chemistry benchmarks using standard Coulomb quantum mechanics",
        ClaimTag.INHERITED,
        Risk.MEDIUM,
        "Valid if substrate reduces to Coulomb + Pauli; many computations are standard QM.",
        dependencies=("P003", "D001"),
        next_step="Separate direct substrate predictions from standard QM inheritance.",
    ),
    LedgerEntry(
        "I002",
        "qed",
        "QED precision results: g-2, Lamb shift, Casimir, 21cm inheritance",
        ClaimTag.INHERITED,
        Risk.MEDIUM,
        "Inherited through low-energy QED matching, not independent substrate derivations.",
        dependencies=("P003", "D002"),
        next_step="Compute alpha and loop structure from substrate Lagrangian.",
    ),
    LedgerEntry(
        "D005",
        "gravity",
        "Gravity/EM force hierarchy numerical ratio",
        ClaimTag.DERIVED,
        Risk.MEDIUM,
        "Matches 0.06% in current formulation but depends on Planck-scale anchor.",
        dependencies=("P001", "D002"),
        next_step="Remove UV/Planck circularity.",
    ),
    LedgerEntry(
        "S001",
        "black_holes",
        "Black-hole saturation/no-singularity ontology",
        ClaimTag.STRUCTURAL,
        Risk.MEDIUM,
        "Coherent consequence of sigma cap; precision deviations not computed.",
        dependencies=("P004", "D005"),
        next_step="Derive ringdown/echo/entropy deviations.",
    ),
    LedgerEntry(
        "D006",
        "qcd",
        "QCD string scale from substrate K(xi) and xi_QCD",
        ClaimTag.DERIVED,
        Risk.MEDIUM,
        "Feeds Regge, proton radius, f_pi, and magnetic moment checks.",
        dependencies=("P001", "D002"),
        next_step="Clarify K_action vs K_coarse.",
    ),
    LedgerEntry(
        "D007",
        "qcd",
        "Proton radius v3 from 3D bound-state wavefunction",
        ClaimTag.DERIVED,
        Risk.LOW,
        "0.808 fm vs 0.841 fm (-3.85%).",
        dependencies=("D006",),
    ),
    LedgerEntry(
        "D008",
        "qcd",
        "Pion decay constant f_pi = 1/2 sigma xi_QCD",
        ClaimTag.DERIVED,
        Risk.LOW,
        "91.22 MeV vs 92.4 MeV (-1.3%).",
        dependencies=("D006",),
    ),
    LedgerEntry(
        "A001",
        "qcd",
        "Hyperon spectrum after one strange-flavour anchor",
        ClaimTag.ANCHORED,
        Risk.MEDIUM,
        "Non-anchor hyperons within 1.68%, but m_s_struct is fixed from Lambda.",
        dependencies=("D006",),
        next_step="Derive strange structural mass instead of anchoring it.",
    ),
    LedgerEntry(
        "A002",
        "qcd",
        "f_K and strange constituent scaling conditional on m_K/m_pi",
        ClaimTag.ANCHORED,
        Risk.MEDIUM,
        "Excellent f_K result, but empirical endpoint mass ratio is supplied.",
        dependencies=("D008",),
        next_step="Derive pseudoscalar endpoint ratio from SU(3) breaking.",
    ),
    LedgerEntry(
        "A003",
        "qcd",
        "p-n mass splitting with substrate isospin candidate",
        ClaimTag.ANCHORED,
        Risk.MEDIUM,
        "Improved to no empirical isospin anchor, but Airy selector is not unique.",
        dependencies=("D006", "P003"),
        next_step="Derive why orientation selects Airy kinetic mass.",
    ),
    LedgerEntry(
        "F001",
        "leptons",
        "Lepton hierarchy from extended Mobius topology",
        ClaimTag.FAILED,
        Risk.BLOCKER,
        "Extended topology gives O(1) spreads, not x207/x3477.",
        dependencies=("P003",),
        next_step="Stop adding Mobius variants; derive vertex eigenvalues kappa_n.",
    ),
    LedgerEntry(
        "O001",
        "leptons",
        "Vertex eigenvalues for charged lepton masses",
        ClaimTag.OPEN,
        Risk.BLOCKER,
        "Needs kappa_mu/kappa_e ~= 4.28e4, kappa_tau/kappa_e ~= 1.21e7, "
        "and positive-root Foot phase; current boundary-loop trial gives ~0.2% errors.",
        dependencies=("P002", "P003"),
        next_step="Derive O_vertex giving delta/pi = 7/12 + (8pi^2)^-1[1-(16pi^2)^-1].",
    ),
    LedgerEntry(
        "O002",
        "flavour",
        "CKM/PMNS flavour mixing operator",
        ClaimTag.OPEN,
        Risk.HIGH,
        "Half-flux overlap trial sin(theta_C)=1/(pi*sqrt(2)) is -0.187%, but no H_mix derivation.",
        dependencies=("P003", "D006"),
        next_step="Derive H_mix overlap integral and test full CKM/PMNS.",
    ),
    LedgerEntry(
        "F002",
        "uv",
        "Planck length from saturation alone",
        ClaimTag.FAILED,
        Risk.BLOCKER,
        "Substrate-only length misses by ~10^22; UV target is chi_UV ~= 4.2e-23 "
        "or action S_UV ~= 51.53; phase-slip trials hit action within 0.1%.",
        dependencies=("P001", "P004", "D002"),
        next_step="Derive closed saturated phase-slip saddle or declare xi_P primitive.",
    ),
    LedgerEntry(
        "O003",
        "cosmology",
        "Pre-CMB proto-matter transfer windows",
        ClaimTag.OPEN,
        Risk.BLOCKER,
        "Needs f_vis <= 4e-4; biharmonic opacity with k_c=0.1 Mpc^-1 passes numerically.",
        dependencies=("P004", "D006"),
        next_step="Derive W_m(k), W_gamma(k), sigma_m, sigma_gamma, and k_c.",
    ),
    LedgerEntry(
        "C001",
        "cosmology",
        "Hubble shift via imposed 7% sound-horizon suppression",
        ClaimTag.CALIBRATED,
        Risk.HIGH,
        "Gives H0 ~= 72.47, but shift is imposed and acoustic peak consistency is open.",
        dependencies=("O003",),
        next_step="Run full CMB/BAO fit from derived P_substrate(k).",
    ),
    LedgerEntry(
        "C002",
        "cosmology",
        "JWST high-z excess factor",
        ClaimTag.CALIBRATED,
        Risk.HIGH,
        "182x at z=14 is calibrated by construction in current module.",
        dependencies=("O003",),
        next_step="Predict luminosity/mass functions without calibration.",
    ),
    LedgerEntry(
        "O004",
        "dark_matter",
        "Dark substrate-stress sector",
        ClaimTag.OPEN,
        Risk.HIGH,
        "Closure candidates: Omega_dark/Omega_b=5.350 (-0.185%), "
        "f_mobile=0.8408, sigma/m~=0.27, ell_pol=alpha^3(c/H0)/sqrt(3), "
        "v_dark=alpha*c/sqrt(5), tau_pol~=48.77 Myr (-0.239%); factor scan "
        "shows tau-only degeneracy but one physical subpercent candidate; "
        "cluster audit requires mobile stress to carry the 150 kpc offset; "
        "finite-speed 1D transport keeps lensing peak offset; sqrt(5) is "
        "rank-5 symmetric-traceless mode count; alpha speed factor equals "
        "second-order neutral stiffness K_eff/K=alpha^2; EM-darkness gate "
        "requires no emission/absorption/reflection channel; geometric action "
        "candidate packages dark stress as -1/2 K alpha^2 Tr(ST(strain)^2).",
        dependencies=("D006", "P003"),
        next_step=(
            "Derive phase-space measure, zero-mode mobile split, alpha^3/sqrt(3) "
            "coherence filtering, second-order neutral stiffness K_eff/K=alpha^2, "
            "and coupled rho_kink/rho_pol dynamics."
        ),
    ),
    LedgerEntry(
        "F003",
        "baryogenesis",
        "One-shot Big-Bang baryogenesis",
        ClaimTag.FAILED,
        Risk.LOW,
        "Retired framing; bubble/sphaleron calculations fail or miss target.",
        dependencies=("P003", "P004"),
        next_step="Replace with orientation inheritance/selection calculation.",
    ),
    LedgerEntry(
        "O005",
        "matter_orientation",
        "Matter-sector orientation selection / inheritance",
        ClaimTag.OPEN,
        Risk.HIGH,
        "Needed to explain stable macroscopic orientation; f_anti<1e-18 requires "
        "DeltaE/T_eff >= 41.45; orientation-vortex trials give S ~= 41.3.",
        dependencies=("P003", "P004"),
        next_step="Derive de-saturation orientation vortex action and determinant.",
    ),
)


def all_entries() -> tuple[LedgerEntry, ...]:
    """Return the immutable ledger entries."""
    return LEDGER


def entries_by_tag(tag: ClaimTag) -> tuple[LedgerEntry, ...]:
    """Return all entries with a given tag."""
    return tuple(entry for entry in LEDGER if entry.tag == tag)


def entries_by_sector(sector: str) -> tuple[LedgerEntry, ...]:
    """Return all entries in a sector."""
    return tuple(entry for entry in LEDGER if entry.sector == sector)


def active_work_queue() -> tuple[LedgerEntry, ...]:
    """Return entries that need active theoretical or numerical work."""
    priority = {
        Risk.BLOCKER: 0,
        Risk.HIGH: 1,
        Risk.MEDIUM: 2,
        Risk.LOW: 3,
    }
    return tuple(
        sorted(
            (entry for entry in LEDGER if entry.needs_work),
            key=lambda entry: (priority[entry.risk], entry.sector, entry.claim_id),
        )
    )


def tag_summary(entries: Iterable[LedgerEntry] = LEDGER) -> Counter[str]:
    """Count entries by dependency tag."""
    return Counter(entry.tag.value for entry in entries)


def risk_summary(entries: Iterable[LedgerEntry] = LEDGER) -> Counter[str]:
    """Count entries by risk level."""
    return Counter(entry.risk.value for entry in entries)


def independent_prediction_count(entries: Iterable[LedgerEntry] = LEDGER) -> int:
    """Count direct derived entries only."""
    return sum(1 for entry in entries if entry.counts_as_independent_prediction)


def format_markdown(entries: Iterable[LedgerEntry] = LEDGER) -> str:
    """Render a compact Markdown dependency ledger."""
    rows = [
        "| ID | Sector | Tag | Risk | Claim | Next step |",
        "|---|---|---|---|---|---|",
    ]
    for entry in entries:
        rows.append(
            "| {id} | {sector} | {tag} | {risk} | {claim} | {next_step} |".format(
                id=entry.claim_id,
                sector=entry.sector,
                tag=entry.tag.value,
                risk=entry.risk.value,
                claim=entry.claim.replace("|", "/"),
                next_step=(entry.next_step or "-").replace("|", "/"),
            )
        )
    return "\n".join(rows)


def format_summary() -> str:
    """Render a human-readable ledger summary."""
    tag_counts = tag_summary()
    risk_counts = risk_summary()
    independent = independent_prediction_count()
    total = len(LEDGER)
    lines = [
        "Dependency ledger summary",
        "=========================",
        f"Total entries: {total}",
        f"Independent derived entries: {independent}",
        "",
        "By tag:",
    ]
    for tag in ClaimTag:
        lines.append(f"  {tag.value:<11} {tag_counts.get(tag.value, 0):>3}")
    lines.extend(["", "By risk:"])
    for risk in Risk:
        lines.append(f"  {risk.value:<8} {risk_counts.get(risk.value, 0):>3}")
    lines.extend(["", "Active blocker/high-risk queue:"])
    for entry in active_work_queue():
        if entry.risk in {Risk.BLOCKER, Risk.HIGH}:
            lines.append(
                f"  {entry.claim_id} [{entry.sector}/{entry.tag.value}/{entry.risk.value}] "
                f"{entry.claim} -> {entry.next_step}"
            )
    return "\n".join(lines)
