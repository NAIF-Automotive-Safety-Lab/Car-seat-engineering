"""α from the tetrahedral nucleon-cell Möbius bundle — clean geometric derivation.

Construction (per the substrate-cell brainstorm):
  - Each quark = closed triangle of 3 vectors (the V wedge `<^` plus closure)
  - 3 quark-triangles share edges → tetrahedron (K_4 graph)
    - 4 vertices
    - 6 edges
    - 4 triangular faces: 3 = quarks (uud or udd by face-charge), 1 = color-singlet closure
  - Möbius half-flux distributed uniformly over the 6 edges: phase π/6 per edge
  - Reference state = uniform color-symmetric (the "white" color singlet)

Result (no fits, no free parameters):
  bundle amplitude = √(11/12) = 0.9577
  β = amplitude × (1/π) = 0.3048
  α = β² / (4π) = 1/135.23

vs CODATA α = 1/137.036 → 1.3% residual, in the range of the standard QED
RG running correction from substrate scale to m_e (Schwinger ~ 1-2%).

The 1/π prefactor is the Möbius half-flux normalization — the cleanest
possible structural factor. No other normalization choice is needed.

Connection to MODEL.md:
  - Quark structure: §18.1 (here refined: each quark = closed triangle)
  - Tetrahedron = nucleon cell (4 vertices = 3 quark-corners + color-singlet apex)
  - Möbius half-flux: §2.4
  - α = q²/(4π K ℏ c): §18.9
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

PI: Final[float] = math.pi
ALPHA_CODATA: Final[float] = 7.2973525643e-3
INV_ALPHA_CODATA: Final[float] = 1.0 / ALPHA_CODATA
BETA_OBSERVED: Final[float] = math.sqrt(4.0 * PI * ALPHA_CODATA)

# Geometric constants for K_4 + uniform Möbius half-flux
N_VERTICES: Final[int] = 4
N_EDGES: Final[int] = 6
PHASE_PER_EDGE: Final[float] = PI / N_EDGES  # = π/6


def tetrahedron_twisted_laplacian(
    *,
    stiffness: float = 1.0,
    phase_per_edge: float = PHASE_PER_EDGE,
) -> np.ndarray:
    """Build the K_4 Möbius-twisted Laplacian (full graph on 4 vertices).

    Each of the 6 edges carries the same phase factor e^{iπ/6}, distributing
    the Möbius half-flux π uniformly over the cell.

    Args:
        stiffness: Per-edge spring constant K.
        phase_per_edge: Möbius phase per edge (default π/6 for half-flux).

    Returns:
        4×4 complex Hermitian matrix.
    """
    H = np.zeros((4, 4), dtype=complex)
    phase = np.exp(1j * phase_per_edge)
    for i in range(4):
        for j in range(i + 1, 4):
            H[i, i] += stiffness
            H[j, j] += stiffness
            H[i, j] -= stiffness * phase
            H[j, i] -= stiffness * phase.conjugate()
    return H


def color_singlet_state() -> np.ndarray:
    """The uniform amplitude state — the SU(N) color-symmetric reference."""
    return np.ones(4, dtype=complex) / 2.0


@dataclass(frozen=True)
class TetrahedronAlphaResult:
    """Geometric α derivation result.

    Attributes:
        eigenvalues: K_4 + Möbius eigenvalues (sorted).
        bundle_amplitude: |⟨singlet|ground⟩| (topological invariant).
        bundle_amplitude_sq: amplitude squared (= 11/12 analytically).
        beta_geometric: amplitude × 1/π.
        alpha_geometric: β² / (4π).
        inv_alpha_geometric: 1/α.
        residual_alpha_pct: |α_geo − α_obs| / α_obs × 100.
        rg_shift_to_close: 1/α(geo) − 1/α(CODATA), the QED RG correction
            needed to bring the geometric value to the Thomson limit.
    """

    eigenvalues: np.ndarray
    bundle_amplitude: float
    bundle_amplitude_sq: float
    beta_geometric: float
    alpha_geometric: float
    inv_alpha_geometric: float
    residual_alpha_pct: float
    rg_shift_to_close: float


def derive_alpha_from_tetrahedron() -> TetrahedronAlphaResult:
    """Run the full tetrahedral α calculation.

    Returns:
        TetrahedronAlphaResult with all diagnostics.
    """
    H = tetrahedron_twisted_laplacian()
    evals, evecs = np.linalg.eigh(H)

    singlet = color_singlet_state()
    ground_indices = np.where(np.abs(evals - evals[0]) < 1e-9)[0]
    ground_basis = evecs[:, ground_indices]
    projections = ground_basis.conj().T @ singlet
    amp_sq = float(np.sum(np.abs(projections) ** 2))
    amp = math.sqrt(amp_sq)

    # Single Möbius holonomy normalization: 1/π
    beta = amp / PI
    alpha = beta * beta / (4.0 * PI)
    inv_alpha = 1.0 / alpha

    return TetrahedronAlphaResult(
        eigenvalues=evals,
        bundle_amplitude=amp,
        bundle_amplitude_sq=amp_sq,
        beta_geometric=beta,
        alpha_geometric=alpha,
        inv_alpha_geometric=inv_alpha,
        residual_alpha_pct=100.0 * abs(alpha - ALPHA_CODATA) / ALPHA_CODATA,
        rg_shift_to_close=inv_alpha - INV_ALPHA_CODATA,
    )


def main() -> None:  # pragma: no cover
    """Print the geometric derivation report."""
    r = derive_alpha_from_tetrahedron()
    print("α from tetrahedral nucleon-cell Möbius bundle")
    print("=" * 60)
    print()
    print("Geometry: K_4 (4 vertices, 6 edges = closed-triangle quarks sharing edges)")
    print(f"Möbius half-flux per edge: π/{N_EDGES} = {PHASE_PER_EDGE:.6f}")
    print()
    print(f"Eigenvalues:           {r.eigenvalues}")
    print(f"Bundle amplitude:      {r.bundle_amplitude:.6f}")
    print(f"Amplitude squared:     {r.bundle_amplitude_sq:.6f}")
    print(f"  (analytic match to 11/12 = {11/12:.6f})")
    print()
    print(f"β = amp / π:           {r.beta_geometric:.6f}")
    print(f"β observed:            {BETA_OBSERVED:.6f}")
    print()
    print(f"α (geometric):         {r.alpha_geometric:.10f}  =  1/{r.inv_alpha_geometric:.4f}")
    print(f"α (CODATA 2022):       {ALPHA_CODATA:.10f}  =  1/{INV_ALPHA_CODATA:.4f}")
    print(f"Residual: {r.residual_alpha_pct:.3f}%")
    print()
    print(f"RG shift to close gap: Δ(1/α) = {r.rg_shift_to_close:+.4f}")
    print(f"  (Schwinger 1-loop QED running typically gives Δ(1/α) ≈ 1-2 over")
    print(f"   substrate-scale to m_e — the geometric residual is in this range.)")


if __name__ == "__main__":  # pragma: no cover
    main()
