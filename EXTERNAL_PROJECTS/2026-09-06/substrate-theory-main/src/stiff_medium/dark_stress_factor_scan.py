"""Small-integer factor scan for dark-stress scale closure.

The scale-closure candidate

    ell_pol = alpha^3 (c/H0) / sqrt(3)
    v_dark = alpha c / sqrt(5)

matches the required cluster memory time. This module checks whether that match
is isolated or part of a wider family of small-integer fits.
"""

from __future__ import annotations

from dataclasses import dataclass

from .dark_stress_hybrid import hybrid_cluster_check
from .dark_stress_memory_clock import relaxation_count
from .dark_stress_scale_closure import ALPHA_EM, MYR_S, hubble_length_kpc
from .substrate_polarization_dm import C_SI, KPC_M


@dataclass(frozen=True)
class FactorCandidate:
    """One dark-stress scale factor candidate."""

    ell_power: int
    ell_projection: int
    speed_power: int
    shear_modes: int
    ell_pol_kpc: float
    v_dark_km_s: float
    tau_pol_myr: float
    tau_error_pct: float

    @property
    def key(self) -> tuple[int, int, int, int]:
        return (
            self.ell_power,
            self.ell_projection,
            self.speed_power,
            self.shear_modes,
        )


@dataclass(frozen=True)
class FactorScanAssessment:
    """Summary of the small-integer scan."""

    total_candidates: int
    physical_candidates: int
    subpercent_tau_candidates: int
    physical_subpercent_tau_candidates: int
    best_tau_candidate: FactorCandidate
    best_physical_candidate: FactorCandidate
    verdict: str


def tau_pol_myr(ell_pol_kpc: float, v_dark_km_s: float) -> float:
    """Return polarization-memory time for a coherence crossing clock."""

    tau_clock_myr = ell_pol_kpc * KPC_M / (v_dark_km_s * 1000.0) / MYR_S
    return relaxation_count() * tau_clock_myr


def scan_factor_candidates(
    max_ell_power: int = 5,
    max_speed_power: int = 3,
    max_projection: int = 9,
    max_shear_modes: int = 9,
) -> tuple[FactorCandidate, ...]:
    """Scan alpha powers and small square-root projection factors."""

    required = hybrid_cluster_check().required_memory_myr
    hubble_kpc = hubble_length_kpc()
    candidates: list[FactorCandidate] = []

    for ell_power in range(1, max_ell_power + 1):
        for ell_projection in range(1, max_projection + 1):
            ell_pol_kpc = ALPHA_EM**ell_power * hubble_kpc / ell_projection**0.5
            for speed_power in range(0, max_speed_power + 1):
                for shear_modes in range(1, max_shear_modes + 1):
                    v_dark_km_s = (
                        ALPHA_EM**speed_power * C_SI / shear_modes**0.5 / 1000.0
                    )
                    tau = tau_pol_myr(ell_pol_kpc, v_dark_km_s)
                    err = (tau / required - 1.0) * 100.0
                    candidates.append(
                        FactorCandidate(
                            ell_power=ell_power,
                            ell_projection=ell_projection,
                            speed_power=speed_power,
                            shear_modes=shear_modes,
                            ell_pol_kpc=ell_pol_kpc,
                            v_dark_km_s=v_dark_km_s,
                            tau_pol_myr=tau,
                            tau_error_pct=err,
                        )
                    )

    return tuple(sorted(candidates, key=lambda item: abs(item.tau_error_pct)))


def is_physical_candidate(
    candidate: FactorCandidate,
    ell_bounds_kpc: tuple[float, float] = (0.1, 10.0),
    speed_bounds_km_s: tuple[float, float] = (300.0, 3000.0),
) -> bool:
    """Return whether the candidate has halo-scale length and cluster-scale speed."""

    return (
        ell_bounds_kpc[0] <= candidate.ell_pol_kpc <= ell_bounds_kpc[1]
        and speed_bounds_km_s[0] <= candidate.v_dark_km_s <= speed_bounds_km_s[1]
    )


def assess_factor_scan() -> FactorScanAssessment:
    """Assess whether the scale-closure factors are unique enough to keep."""

    candidates = scan_factor_candidates()
    physical = tuple(item for item in candidates if is_physical_candidate(item))
    subpercent = tuple(item for item in candidates if abs(item.tau_error_pct) < 1.0)
    physical_subpercent = tuple(
        item for item in physical if abs(item.tau_error_pct) < 1.0
    )
    best_tau = candidates[0]
    best_physical = physical[0]

    if physical_subpercent and physical_subpercent[0].key == (3, 3, 1, 5):
        verdict = (
            "combined tau alone is degenerate, but halo-scale length plus "
            "cluster-speed filters select alpha^3/sqrt(3) and alpha/sqrt(5); "
            "derive these filters from the neutral-stress field equations"
        )
    else:
        verdict = (
            "factor scan does not uniquely support the proposed dark-stress "
            "scale closure"
        )

    return FactorScanAssessment(
        total_candidates=len(candidates),
        physical_candidates=len(physical),
        subpercent_tau_candidates=len(subpercent),
        physical_subpercent_tau_candidates=len(physical_subpercent),
        best_tau_candidate=best_tau,
        best_physical_candidate=best_physical,
        verdict=verdict,
    )
