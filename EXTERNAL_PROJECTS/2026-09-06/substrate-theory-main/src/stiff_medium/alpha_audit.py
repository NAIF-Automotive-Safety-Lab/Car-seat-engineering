"""Reusable fine-structure-constant audit summary."""

from __future__ import annotations

from dataclasses import dataclass

from stiff_medium.alpha_derivation import (
    AlphaDerivationAudit,
    run_alpha_derivation_audit,
)
from stiff_medium.alpha_two_loop import (
    AlphaTwoLoopAudit,
    TwoLoopResult,
    run_alpha_two_loop_audit,
)
from stiff_medium.physical_constants import (
    ALPHA_THOMSON_CODATA_2022,
    INV_ALPHA_THOMSON_CODATA_2022,
)


@dataclass(frozen=True)
class AlphaAuditSummary:
    """Single machine-checkable conclusion for the current alpha work."""

    derivation: AlphaDerivationAudit
    two_loop: AlphaTwoLoopAudit
    fixed_delta_inv_alpha: float
    fixed_abs_pct_inv_alpha: float
    fixed_delta_alpha: float
    conclusion: str


def signed_delta_inv(result: TwoLoopResult, target_inv: float) -> float:
    return result.inv_alpha_2loop - target_inv


def signed_delta_alpha(result: TwoLoopResult, target_alpha: float) -> float:
    return result.alpha_2loop - target_alpha


def pct_inv(result: TwoLoopResult, target_inv: float) -> float:
    return abs(signed_delta_inv(result, target_inv)) / target_inv * 100.0


def validate_target_constants(*audits: AlphaDerivationAudit | AlphaTwoLoopAudit) -> None:
    """Fail fast if an audit section is not using the shared CODATA target."""
    for audit in audits:
        if audit.target_inv_alpha != INV_ALPHA_THOMSON_CODATA_2022:
            raise RuntimeError(
                "audit target inverse alpha drifted from "
                "INV_ALPHA_THOMSON_CODATA_2022"
            )
        if audit.target_alpha != ALPHA_THOMSON_CODATA_2022:
            raise RuntimeError("audit target alpha drifted from ALPHA_THOMSON_CODATA_2022")


def build_alpha_audit_summary() -> AlphaAuditSummary:
    """Run both alpha audit sections and return the compact conclusion."""
    derivation = run_alpha_derivation_audit()
    two_loop = run_alpha_two_loop_audit()
    validate_target_constants(derivation, two_loop)

    fixed = two_loop.fixed_beta
    fixed_delta_inv = signed_delta_inv(fixed, two_loop.target_inv_alpha)
    fixed_pct_inv = pct_inv(fixed, two_loop.target_inv_alpha)
    fixed_delta_alpha = signed_delta_alpha(fixed, two_loop.target_alpha)

    return AlphaAuditSummary(
        derivation=derivation,
        two_loop=two_loop,
        fixed_delta_inv_alpha=fixed_delta_inv,
        fixed_abs_pct_inv_alpha=fixed_pct_inv,
        fixed_delta_alpha=fixed_delta_alpha,
        conclusion=(
            "No first-principles alpha derivation: best fixed two-loop "
            f"1/alpha={fixed.inv_alpha_2loop:.9f}, "
            f"delta={fixed_delta_inv:+.9f} ({fixed_pct_inv:.9f}%), "
            "and exact agreement requires fitting beta."
        ),
    )
