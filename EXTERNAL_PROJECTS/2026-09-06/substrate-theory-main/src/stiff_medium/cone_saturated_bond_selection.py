"""Energetic selection test for a one-bond saturated phase slip.

The lattice boundary audit shows that a single open bond has the right
two-anchor topology.  This module checks a harder question: does the saturation
barrier alone energetically localize an imposed phase slip onto one bond?

It does not.  For a convex barrier, spreading a fixed slip over more bonds
lowers the barrier energy.  A one-bond saturated segment therefore needs an
additional localization mechanism, such as a Peierls/core cost for activating
multiple bonds or a solved external load/saddle that pins the slip.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


def saturation_barrier_energy(strain_fraction: float) -> float:
    """Return the dimensionless barrier energy with the zero-strain value removed."""

    f = abs(float(strain_fraction))
    if not 0.0 <= f < 1.0:
        raise ValueError("strain_fraction must satisfy |f| < 1")
    return 1.0 / math.sqrt(1.0 - f * f) - 1.0


def distributed_slip_energy(
    total_strain_fraction: float,
    bond_count: int,
    *,
    core_cost_per_active_bond: float = 0.0,
) -> float:
    """Return energy for distributing a fixed slip uniformly over active bonds."""

    if bond_count <= 0:
        raise ValueError("bond_count must be positive")
    if core_cost_per_active_bond < 0.0:
        raise ValueError("core_cost_per_active_bond must be non-negative")
    per_bond_fraction = total_strain_fraction / bond_count
    return (
        bond_count * saturation_barrier_energy(per_bond_fraction)
        + bond_count * core_cost_per_active_bond
    )


def slip_energy_table(
    total_strain_fraction: float,
    bond_counts: tuple[int, ...],
    *,
    core_cost_per_active_bond: float = 0.0,
) -> tuple[tuple[int, float], ...]:
    """Return (bond_count, energy) rows for an imposed phase slip."""

    return tuple(
        (
            count,
            distributed_slip_energy(
                total_strain_fraction,
                count,
                core_cost_per_active_bond=core_cost_per_active_bond,
            ),
        )
        for count in bond_counts
    )


def selected_bond_count(table: tuple[tuple[int, float], ...]) -> int:
    """Return the bond count with minimal energy."""

    if not table:
        raise ValueError("energy table must be non-empty")
    return min(table, key=lambda row: row[1])[0]


def critical_core_cost_for_single_bond(
    total_strain_fraction: float,
    bond_counts: tuple[int, ...],
) -> float:
    """Return the minimum per-bond core cost that makes one bond competitive."""

    if 1 not in bond_counts:
        raise ValueError("bond_counts must include 1")
    one_bond_energy = distributed_slip_energy(total_strain_fraction, 1)
    required = 0.0
    for count in bond_counts:
        if count == 1:
            continue
        spread_energy = distributed_slip_energy(total_strain_fraction, count)
        required = max(required, (one_bond_energy - spread_energy) / (count - 1))
    return required


@dataclass(frozen=True)
class SaturatedBondSelectionAssessment:
    """Assessment of saturated-bond energetic selection."""

    total_strain_fraction: float
    pure_selected_bond_count: int
    pure_one_bond_energy: float
    pure_widest_energy: float
    pure_barrier_selects_single_bond: bool
    critical_core_cost: float
    trial_core_cost: float
    core_selected_bond_count: int
    core_cost_fixed_by_substrate: bool
    fully_derived: bool
    verdict: str


def assess_saturated_bond_selection(
    *,
    total_strain_fraction: float = 0.9,
    bond_counts: tuple[int, ...] = (1, 2, 4, 8, 16, 32, 64),
) -> SaturatedBondSelectionAssessment:
    """Assess whether the barrier alone selects the one-bond phase-slip segment."""

    pure_table = slip_energy_table(total_strain_fraction, bond_counts)
    pure_selected = selected_bond_count(pure_table)
    one_energy = dict(pure_table)[1]
    widest_count = max(bond_counts)
    widest_energy = dict(pure_table)[widest_count]
    critical_core = critical_core_cost_for_single_bond(
        total_strain_fraction,
        bond_counts,
    )
    trial_core = 1.01 * critical_core
    core_table = slip_energy_table(
        total_strain_fraction,
        bond_counts,
        core_cost_per_active_bond=trial_core,
    )
    core_selected = selected_bond_count(core_table)
    barrier_selects = pure_selected == 1

    if not barrier_selects and core_selected == 1:
        verdict = (
            "the saturation barrier alone delocalizes an imposed phase slip; "
            "a positive Peierls/core localization cost can select one bond, "
            "but that cost is not yet derived from the substrate"
        )
    elif barrier_selects:
        verdict = "the saturation barrier alone selects the one-bond phase slip"
    else:
        verdict = "one-bond phase-slip selection remains energetically unclosed"

    return SaturatedBondSelectionAssessment(
        total_strain_fraction=total_strain_fraction,
        pure_selected_bond_count=pure_selected,
        pure_one_bond_energy=one_energy,
        pure_widest_energy=widest_energy,
        pure_barrier_selects_single_bond=barrier_selects,
        critical_core_cost=critical_core,
        trial_core_cost=trial_core,
        core_selected_bond_count=core_selected,
        core_cost_fixed_by_substrate=False,
        fully_derived=False,
        verdict=verdict,
    )
