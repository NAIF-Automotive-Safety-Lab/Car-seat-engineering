"""Neutral-stress tensor mode count for the dark-stress speed factor.

The dark-stress scale closure uses

    v_dark = alpha c / sqrt(5).

This module checks the structural part of that claim.  In three spatial
dimensions, a symmetric stress tensor has six independent components. Removing
the scalar trace leaves a rank-5 symmetric-traceless sector.  If neutral stress
energy is distributed isotropically over those five components, the RMS group
speed per mode is alpha*c/sqrt(5).

This is still not the full alpha derivation.  It only turns the denominator
sqrt(5) from a fit into a tensor-mode count.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .dark_stress_scale_closure import ALPHA_EM
from .substrate_polarization_dm import C_SI


def symmetric_tensor_trace_projector() -> np.ndarray:
    """Return the 6D projector onto 3D symmetric-traceless tensors.

    Basis order: xx, yy, zz, xy, xz, yz.  The trace direction is proportional
    to (1, 1, 1, 0, 0, 0).
    """

    identity = np.eye(6)
    trace_direction = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    trace_direction /= np.linalg.norm(trace_direction)
    return identity - np.outer(trace_direction, trace_direction)


def projector_rank(projector: np.ndarray, tolerance: float = 1.0e-10) -> int:
    """Return numerical rank of a projector-like matrix."""

    eigenvalues = np.linalg.eigvalsh(projector)
    return int(np.count_nonzero(eigenvalues > tolerance))


@dataclass(frozen=True)
class NeutralStressModeAssessment:
    """Assessment of the symmetric-traceless neutral-stress sector."""

    projector_rank: int
    trace_eigenvalue: float
    shear_eigenvalues: tuple[float, ...]
    idempotence_error: float
    speed_formula: str
    v_dark_km_s: float
    verdict: str


def assess_neutral_stress_modes() -> NeutralStressModeAssessment:
    """Assess whether the neutral-stress mode count supports sqrt(5)."""

    projector = symmetric_tensor_trace_projector()
    eigenvalues = tuple(float(value) for value in np.linalg.eigvalsh(projector))
    rank = projector_rank(projector)
    trace_eigenvalue = min(eigenvalues)
    shear_eigenvalues = tuple(value for value in eigenvalues if value > 0.5)
    idempotence_error = float(np.linalg.norm(projector @ projector - projector))
    v_dark_km_s = ALPHA_EM * C_SI / math.sqrt(rank) / 1000.0

    if rank == 5 and idempotence_error < 1.0e-12:
        verdict = (
            "sqrt(5) is the symmetric-traceless stress mode count; alpha "
            "suppression still needs a neutral-coupling derivation"
        )
    else:
        verdict = "neutral-stress projector does not support the sqrt(5) factor"

    return NeutralStressModeAssessment(
        projector_rank=rank,
        trace_eigenvalue=trace_eigenvalue,
        shear_eigenvalues=shear_eigenvalues,
        idempotence_error=idempotence_error,
        speed_formula="alpha*c/sqrt(rank_ST), rank_ST=5",
        v_dark_km_s=v_dark_km_s,
        verdict=verdict,
    )
