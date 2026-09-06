"""1D propagating EM field — coupling sources to distant absorbers.

Implements spec §18.20: charges couple to the medium's wave field,
disturbances propagate at c, distant resonant masses absorb the
energy. Wave field φ(x, t) on a grid evolves per ∂²φ/∂t² = c² ∂²φ/∂x².
Charged particles act as sources (proportional to acceleration);
distant particles feel the field as a force (proportional to gradient).

Energy budget tracked through:
- Source kinetic energy → field energy → absorber kinetic energy
"""

import numpy as np


class EMField1D:
    """1D scalar wave field on a uniform grid.

    Discretized wave equation:
        φ(x, t+dt) = 2 φ(x, t) − φ(x, t−dt) + c²·dt²/dx² · ∇²φ

    Sources add to ∂²φ/∂t² at their positions (proportional to charge × acceleration).
    Field evaluated at particle positions provides the EM force.
    """

    def __init__(
        self,
        x_min: float,
        x_max: float,
        n_points: int,
        c: float = 1.0,
    ):
        self.x_min = x_min
        self.x_max = x_max
        self.n = n_points
        self.dx = (x_max - x_min) / (n_points - 1)
        self.c = c
        self.x = np.linspace(x_min, x_max, n_points)
        self.phi_prev = np.zeros(n_points)
        self.phi_curr = np.zeros(n_points)

    def step(self, dt: float, sources: dict[float, float] | None = None) -> None:
        """One wave-equation timestep.

        sources: dict of {x_position: source_amplitude} — adds to ∂²φ/∂t² at those points.
        """
        # CFL stability: c·dt/dx < 1
        if self.c * dt / self.dx > 1.0:
            raise ValueError(f"CFL violated: c·dt/dx = {self.c * dt / self.dx} > 1")

        coeff = (self.c * dt / self.dx) ** 2
        # Compute Laplacian on interior points
        new_phi = np.zeros_like(self.phi_curr)
        new_phi[1:-1] = (
            2 * self.phi_curr[1:-1]
            - self.phi_prev[1:-1]
            + coeff * (self.phi_curr[2:] - 2 * self.phi_curr[1:-1] + self.phi_curr[:-2])
        )
        # Boundary conditions: absorbing (1st-order Mur) so waves don't reflect
        new_phi[0] = self.phi_curr[1] + (self.c * dt - self.dx) / (self.c * dt + self.dx) * (
            new_phi[1] - self.phi_curr[0]
        )
        new_phi[-1] = self.phi_curr[-2] + (self.c * dt - self.dx) / (self.c * dt + self.dx) * (
            new_phi[-2] - self.phi_curr[-1]
        )

        # Apply sources: add to ∂²φ/∂t² at source positions
        if sources is not None:
            for x_src, amplitude in sources.items():
                idx = int(round((x_src - self.x_min) / self.dx))
                if 0 <= idx < self.n:
                    new_phi[idx] += amplitude * dt * dt

        # Roll buffers
        self.phi_prev = self.phi_curr
        self.phi_curr = new_phi

    def gradient_at(self, x: float) -> float:
        """Field gradient at position x (gives force on a charge)."""
        idx = int(round((x - self.x_min) / self.dx))
        if idx <= 0 or idx >= self.n - 1:
            return 0.0
        return (self.phi_curr[idx + 1] - self.phi_curr[idx - 1]) / (2 * self.dx)

    def value_at(self, x: float) -> float:
        idx = int(round((x - self.x_min) / self.dx))
        if 0 <= idx < self.n:
            return float(self.phi_curr[idx])
        return 0.0

    def total_energy(self) -> float:
        """Field energy: ∫ ½(∂_t φ)² + ½ c² (∂_x φ)² dx (in arbitrary units)."""
        # Time derivative approximated by (φ_curr - φ_prev) / dt — caller
        # would need to pass dt; for static evaluation, just spatial part:
        if len(self.phi_curr) < 2:
            return 0.0
        spatial = 0.5 * self.c ** 2 * np.sum(np.diff(self.phi_curr) ** 2) / self.dx
        return float(spatial)
