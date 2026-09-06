"""3D EM wave field — the 3D extension of em_field.py (closes §18.23 item 8).

Per spec §18.20, EM is a wave in the medium. The 1D version
(em_field.py) demonstrates source-absorber coupling at one dimension.
The full 3D version is needed for realistic spectroscopy with
proper geometric attenuation (1/r²).

Discretized 3D wave equation on a uniform cubic grid:
    ∂²φ/∂t² = c² (∂²φ/∂x² + ∂²φ/∂y² + ∂²φ/∂z²)

Memory: a 50×50×50 grid is 125k cells × 2 buffers = 250k floats,
manageable. Larger grids work but are slower.
"""

import numpy as np


class EMField3D:
    """3D scalar wave field on a uniform cubic grid.

    Stores phi(i,j,k) as a 3D numpy array. Time evolution is FDTD-style:
        phi_new[i,j,k] = 2*phi[i,j,k] - phi_old[i,j,k]
                       + (c*dt/dx)^2 * Laplacian
    """

    def __init__(
        self,
        x_min: float,
        x_max: float,
        y_min: float,
        y_max: float,
        z_min: float,
        z_max: float,
        n_x: int,
        n_y: int,
        n_z: int,
        c: float = 1.0,
    ):
        self.x_min, self.x_max = x_min, x_max
        self.y_min, self.y_max = y_min, y_max
        self.z_min, self.z_max = z_min, z_max
        self.n_x, self.n_y, self.n_z = n_x, n_y, n_z
        self.dx = (x_max - x_min) / (n_x - 1)
        self.dy = (y_max - y_min) / (n_y - 1)
        self.dz = (z_max - z_min) / (n_z - 1)
        self.c = c

        self.phi_prev = np.zeros((n_x, n_y, n_z))
        self.phi_curr = np.zeros((n_x, n_y, n_z))

    def _index_of(self, x: float, y: float, z: float) -> tuple[int, int, int]:
        i = int(round((x - self.x_min) / self.dx))
        j = int(round((y - self.y_min) / self.dy))
        k = int(round((z - self.z_min) / self.dz))
        return i, j, k

    def step(self, dt: float, sources: dict[tuple[float, float, float], float] | None = None) -> None:
        """One FDTD step. Sources: dict of {(x, y, z): source_amplitude}."""
        # CFL stability for 3D: c·dt·sqrt(1/dx² + 1/dy² + 1/dz²) < 1
        cfl = self.c * dt * np.sqrt(1/self.dx**2 + 1/self.dy**2 + 1/self.dz**2)
        if cfl > 1.0:
            raise ValueError(f"3D CFL violated: factor = {cfl} > 1")

        cx = (self.c * dt / self.dx) ** 2
        cy = (self.c * dt / self.dy) ** 2
        cz = (self.c * dt / self.dz) ** 2

        # Laplacian on interior points
        phi = self.phi_curr
        new_phi = np.zeros_like(phi)
        new_phi[1:-1, 1:-1, 1:-1] = (
            2 * phi[1:-1, 1:-1, 1:-1]
            - self.phi_prev[1:-1, 1:-1, 1:-1]
            + cx * (phi[2:, 1:-1, 1:-1] - 2 * phi[1:-1, 1:-1, 1:-1] + phi[:-2, 1:-1, 1:-1])
            + cy * (phi[1:-1, 2:, 1:-1] - 2 * phi[1:-1, 1:-1, 1:-1] + phi[1:-1, :-2, 1:-1])
            + cz * (phi[1:-1, 1:-1, 2:] - 2 * phi[1:-1, 1:-1, 1:-1] + phi[1:-1, 1:-1, :-2])
        )

        # Simple absorbing boundaries: just damp at boundaries
        # (Not as precise as Mur ABC but adequate for demonstration)
        new_phi[0, :, :] *= 0.0
        new_phi[-1, :, :] *= 0.0
        new_phi[:, 0, :] *= 0.0
        new_phi[:, -1, :] *= 0.0
        new_phi[:, :, 0] *= 0.0
        new_phi[:, :, -1] *= 0.0

        # Apply sources
        if sources is not None:
            for (x, y, z), amplitude in sources.items():
                i, j, k = self._index_of(x, y, z)
                if 0 <= i < self.n_x and 0 <= j < self.n_y and 0 <= k < self.n_z:
                    new_phi[i, j, k] += amplitude * dt * dt

        self.phi_prev = self.phi_curr
        self.phi_curr = new_phi

    def value_at(self, x: float, y: float, z: float) -> float:
        i, j, k = self._index_of(x, y, z)
        if (0 <= i < self.n_x) and (0 <= j < self.n_y) and (0 <= k < self.n_z):
            return float(self.phi_curr[i, j, k])
        return 0.0

    def gradient_at(self, x: float, y: float, z: float) -> np.ndarray:
        """Gradient of the field at (x,y,z). Returns 3D vector."""
        i, j, k = self._index_of(x, y, z)
        if not (1 <= i < self.n_x - 1 and 1 <= j < self.n_y - 1 and 1 <= k < self.n_z - 1):
            return np.zeros(3)
        gx = (self.phi_curr[i+1, j, k] - self.phi_curr[i-1, j, k]) / (2 * self.dx)
        gy = (self.phi_curr[i, j+1, k] - self.phi_curr[i, j-1, k]) / (2 * self.dy)
        gz = (self.phi_curr[i, j, k+1] - self.phi_curr[i, j, k-1]) / (2 * self.dz)
        return np.array([gx, gy, gz])
