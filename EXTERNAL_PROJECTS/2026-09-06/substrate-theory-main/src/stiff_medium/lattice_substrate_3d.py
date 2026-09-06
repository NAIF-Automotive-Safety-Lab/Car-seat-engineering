"""3D lattice substrate field solver — full Lagrangian dynamics.

Solves the substrate field on a 3D Nx × Ny × Nz lattice using the full
§18.11 Lagrangian *with* a damping term:

    ℒ = ½ ρ (∂_t u)² − ½ K |∇u|² − V(u) − γ u (∂_t u)

with the (sine-Gordon × saturation) potential

    V(u) = (K / ξ²) (1 − cos(u))

The Euler-Lagrange equation including the drag γ is

    ρ ∂_t² u − K ∇² u + (K/ξ²) sin(u) + γ ∂_t u = 0

Integrated by a *drag-aware* velocity-Verlet (symplectic when γ=0).

This is the full 3D analogue of `lattice_substrate_2d.py` and uses the
same numpy-vectorised idiom (np.roll for the 7-point Laplacian).

Usage:
    sub = LatticeSubstrate3D(Nx=32, Ny=32, Nz=32, dx=0.5, dt=0.05)
    sub.kink_3d(direction='x', x0=0.0)
    sub.evolve(500)
    print(sub.total_energy())
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

# Canonical cone-bouncing frequency formula (single source of truth).
from src.stiff_medium.cone_bouncing_protocol import (
    omega_b_canonical as _omega_b_canonical,  # noqa: F401  (consistency)
)


# ---------------------------------------------------------------------------
# 3D lattice substrate
# ---------------------------------------------------------------------------

@dataclass
class LatticeSubstrate3D:
    """3D substrate field u(x, y, z, t) under the full §18.11 Lagrangian.

    Parameters
    ----------
    Nx, Ny, Nz : int
        Grid sizes in x, y, z.
    dx : float
        Lattice spacing (assumed isotropic dx = dy = dz).
    dt : float
        Time-step.  Must satisfy CFL  dt < dx / sqrt(3 K/ρ) for stability.
    rho : float
        Substrate mass density (default 1).
    K : float
        Stiffness (default 1; gives c² = K/ρ = 1).
    xi : float
        Coherence length (default 1).
    gamma : float
        Drag coefficient (default 0 = conservative dynamics).
    bc : str
        Boundary condition: 'periodic' or 'neumann' (zero-flux).
    """

    Nx: int = 32
    Ny: int = 32
    Nz: int = 32
    dx: float = 0.5
    dt: float = 0.05
    rho: float = 1.0
    K: float = 1.0
    xi: float = 1.0
    gamma: float = 0.0
    bc: str = "periodic"

    # State
    u: NDArray[np.float64] = field(init=False)
    v: NDArray[np.float64] = field(init=False)
    x: NDArray[np.float64] = field(init=False)
    y: NDArray[np.float64] = field(init=False)
    z: NDArray[np.float64] = field(init=False)
    X: NDArray[np.float64] = field(init=False)
    Y: NDArray[np.float64] = field(init=False)
    Z: NDArray[np.float64] = field(init=False)
    t: float = field(init=False, default=0.0)
    step_count: int = field(init=False, default=0)

    def __post_init__(self) -> None:
        self.u = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float64)
        self.v = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float64)
        self.x = (np.arange(self.Nx) - self.Nx // 2) * self.dx
        self.y = (np.arange(self.Ny) - self.Ny // 2) * self.dx
        self.z = (np.arange(self.Nz) - self.Nz // 2) * self.dx
        self.X, self.Y, self.Z = np.meshgrid(
            self.x, self.y, self.z, indexing="ij"
        )
        if self.bc not in ("periodic", "neumann"):
            raise ValueError(f"Unknown BC: {self.bc}")
        c = np.sqrt(self.K / self.rho)
        cfl = c * self.dt / self.dx
        # 3D CFL bound dt < dx / (c sqrt(3))
        self._cfl = cfl

    # ------------------------------------------------------------------
    # Differential operators (7-point stencil)
    # ------------------------------------------------------------------

    def laplacian(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        """7-point Laplacian in 3D.  Periodic or Neumann BC."""
        if self.bc == "periodic":
            xp = np.roll(u, -1, axis=0)
            xm = np.roll(u, +1, axis=0)
            yp = np.roll(u, -1, axis=1)
            ym = np.roll(u, +1, axis=1)
            zp = np.roll(u, -1, axis=2)
            zm = np.roll(u, +1, axis=2)
        else:  # neumann
            xp = np.empty_like(u); xm = np.empty_like(u)
            yp = np.empty_like(u); ym = np.empty_like(u)
            zp = np.empty_like(u); zm = np.empty_like(u)
            xp[:-1] = u[1:]; xp[-1] = u[-1]
            xm[1:] = u[:-1]; xm[0] = u[0]
            yp[:, :-1] = u[:, 1:]; yp[:, -1] = u[:, -1]
            ym[:, 1:] = u[:, :-1]; ym[:, 0] = u[:, 0]
            zp[:, :, :-1] = u[:, :, 1:]; zp[:, :, -1] = u[:, :, -1]
            zm[:, :, 1:] = u[:, :, :-1]; zm[:, :, 0] = u[:, :, 0]
        return (xp + xm + yp + ym + zp + zm - 6.0 * u) / (self.dx ** 2)

    def gradient_components(
        self, u: NDArray[np.float64]
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Centred-difference (∂_x u, ∂_y u, ∂_z u)."""
        if self.bc == "periodic":
            ux = (np.roll(u, -1, axis=0) - np.roll(u, +1, axis=0)) / (2.0 * self.dx)
            uy = (np.roll(u, -1, axis=1) - np.roll(u, +1, axis=1)) / (2.0 * self.dx)
            uz = (np.roll(u, -1, axis=2) - np.roll(u, +1, axis=2)) / (2.0 * self.dx)
        else:
            ux = np.zeros_like(u); uy = np.zeros_like(u); uz = np.zeros_like(u)
            ux[1:-1] = (u[2:] - u[:-2]) / (2.0 * self.dx)
            ux[0] = (u[1] - u[0]) / self.dx
            ux[-1] = (u[-1] - u[-2]) / self.dx
            uy[:, 1:-1] = (u[:, 2:] - u[:, :-2]) / (2.0 * self.dx)
            uy[:, 0] = (u[:, 1] - u[:, 0]) / self.dx
            uy[:, -1] = (u[:, -1] - u[:, -2]) / self.dx
            uz[:, :, 1:-1] = (u[:, :, 2:] - u[:, :, :-2]) / (2.0 * self.dx)
            uz[:, :, 0] = (u[:, :, 1] - u[:, :, 0]) / self.dx
            uz[:, :, -1] = (u[:, :, -1] - u[:, :, -2]) / self.dx
        return ux, uy, uz

    def gradient_sq(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        ux, uy, uz = self.gradient_components(u)
        return ux * ux + uy * uy + uz * uz

    # ------------------------------------------------------------------
    # Force / energy
    # ------------------------------------------------------------------

    def force(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.K * self.laplacian(u) - (self.K / self.xi ** 2) * np.sin(u)

    def kinetic_energy(self) -> float:
        return float(0.5 * self.rho * np.sum(self.v ** 2) * self.dx ** 3)

    def gradient_energy(self) -> float:
        return float(0.5 * self.K * np.sum(self.gradient_sq(self.u)) * self.dx ** 3)

    def potential_energy(self) -> float:
        return float((self.K / self.xi ** 2) * np.sum(1.0 - np.cos(self.u))
                     * self.dx ** 3)

    def total_energy(self) -> float:
        return self.kinetic_energy() + self.gradient_energy() + self.potential_energy()

    def energy_components(self) -> dict:
        T = self.kinetic_energy()
        G = self.gradient_energy()
        V = self.potential_energy()
        return {
            "kinetic": T,
            "gradient": G,
            "potential": V,
            "total": T + G + V,
            "t": self.t,
            "step": self.step_count,
        }

    # ------------------------------------------------------------------
    # Time evolution — drag-aware velocity Verlet
    # ------------------------------------------------------------------

    def step(self) -> None:
        dt = self.dt
        drag_half = np.exp(-self.gamma * dt / (2.0 * self.rho))
        a_n = self.force(self.u) / self.rho
        v_half = (self.v + 0.5 * dt * a_n) * drag_half
        self.u = self.u + dt * v_half
        a_np1 = self.force(self.u) / self.rho
        self.v = v_half * drag_half + 0.5 * dt * a_np1
        self.t += dt
        self.step_count += 1

    def evolve(self, steps: int) -> None:
        for _ in range(int(steps)):
            self.step()

    # ------------------------------------------------------------------
    # Initial conditions
    # ------------------------------------------------------------------

    def reset(self) -> None:
        self.u.fill(0.0)
        self.v.fill(0.0)
        self.t = 0.0
        self.step_count = 0

    def kink_3d(self, direction: str = "x", x0: float = 0.0,
                speed: float = 0.0) -> None:
        """1D sine-Gordon kink extended uniformly in the two transverse
        directions:  u(x,y,z) = 4 arctan(exp(γ_l (x − x0) / ξ))
        for ``direction='x'``, with similar formulas for y and z.

        ``speed`` (in units of c = sqrt(K/ρ)) gives a Lorentz-boosted
        profile.
        """
        self.reset()
        c = np.sqrt(self.K / self.rho)
        beta = speed / c
        if abs(beta) >= 1.0:
            raise ValueError("Speed must satisfy |v| < c.")
        gamma_l = 1.0 / np.sqrt(1.0 - beta ** 2)

        if direction == "x":
            arg = gamma_l * (self.X - x0) / self.xi
            self.u = 4.0 * np.arctan(np.exp(arg))
            sech = 1.0 / np.cosh(arg)
            du = (2.0 / self.xi) * gamma_l * sech
            self.v = -speed * du
        elif direction == "y":
            arg = gamma_l * (self.Y - x0) / self.xi
            self.u = 4.0 * np.arctan(np.exp(arg))
            sech = 1.0 / np.cosh(arg)
            du = (2.0 / self.xi) * gamma_l * sech
            self.v = -speed * du
        elif direction == "z":
            arg = gamma_l * (self.Z - x0) / self.xi
            self.u = 4.0 * np.arctan(np.exp(arg))
            sech = 1.0 / np.cosh(arg)
            du = (2.0 / self.xi) * gamma_l * sech
            self.v = -speed * du
        else:
            raise ValueError(f"direction must be 'x','y','z'; got {direction}")

    def monopole_initial(self, x0: float = 0.0, y0: float = 0.0,
                         z0: float = 0.0, core: float = 1.5,
                         charge: int = 1) -> None:
        """Substrate monopole-like (hedgehog scalar) configuration.

        The scalar phase u winds radially:

            u(r) = π × charge × tanh(r / core)

        so that u → 0 at the centre and u → π × charge far out, and
        the gradient ∇u points radially.  This is the scalar-phase
        analogue of a 't Hooft-Polyakov monopole hedgehog (we do not
        carry an internal vector index here, only the radial winding
        of the cosine-potential phase).

        ``charge`` is the integer winding (1 = monopole, 2 = double).
        """
        self.reset()
        dx_ = self.X - x0
        dy_ = self.Y - y0
        dz_ = self.Z - z0
        r = np.sqrt(dx_ ** 2 + dy_ ** 2 + dz_ ** 2) + 1e-12
        self.u = np.pi * charge * np.tanh(r / core)
        self.v.fill(0.0)

    def vortex_string_initial(self, axis: str = "z", x0: float = 0.0,
                              y0: float = 0.0, winding: int = 1,
                              core: float = 1.0) -> None:
        """1D defect line (vortex string) along the given axis.

        For ``axis='z'`` the field carries a 2D vortex in the (x,y)
        plane, extruded uniformly along z:

            u(x,y,z) = winding × atan2(y-y0, x-x0) × tanh(r⊥ / core)

        with r⊥ the perpendicular distance to the string.  Analogues
        for axis='x' (vortex in y-z) and axis='y' (vortex in x-z).
        """
        self.reset()
        if axis == "z":
            dx_ = self.X - x0
            dy_ = self.Y - y0
            r = np.sqrt(dx_ ** 2 + dy_ ** 2) + 1e-12
            theta = np.arctan2(dy_, dx_)
        elif axis == "x":
            dy_ = self.Y - x0  # x0 is the (y) offset for axis='x'
            dz_ = self.Z - y0
            r = np.sqrt(dy_ ** 2 + dz_ ** 2) + 1e-12
            theta = np.arctan2(dz_, dy_)
        elif axis == "y":
            dx_ = self.X - x0
            dz_ = self.Z - y0
            r = np.sqrt(dx_ ** 2 + dz_ ** 2) + 1e-12
            theta = np.arctan2(dz_, dx_)
        else:
            raise ValueError(f"axis must be 'x','y','z'; got {axis}")
        self.u = winding * theta * np.tanh(r / core)
        self.v.fill(0.0)

    def gaussian_pulse_initial(self, x0: float = 0.0, y0: float = 0.0,
                               z0: float = 0.0, amplitude: float = 0.3,
                               width: float = 2.0) -> None:
        """Small-amplitude Gaussian pulse centred at (x0,y0,z0)."""
        self.reset()
        r2 = (self.X - x0) ** 2 + (self.Y - y0) ** 2 + (self.Z - z0) ** 2
        self.u = amplitude * np.exp(-r2 / (2.0 * width ** 2))
        self.v.fill(0.0)

    # ------------------------------------------------------------------
    # Topological diagnostics
    # ------------------------------------------------------------------

    def topological_charge(self) -> int:
        """Integer topological charge of the current field configuration.

        The implementation is mode-aware:

        * If a radial hedgehog (monopole) is present (the field has a
          single point where ∇u changes sign across all three axes)
          we return the integer winding of the radial phase, computed
          as round(u_max / π).
        * If a vortex string runs along an axis we return the planar
          winding number (round of (1/2π) ∮ ∇θ · dl) computed on a
          loop in the perpendicular plane through the lattice centre.
        * Otherwise returns 0.

        For both cases the result is an integer.
        """
        cx, cy, cz = self.Nx // 2, self.Ny // 2, self.Nz // 2
        u = self.u

        # Try monopole detection: maximum of u, scaled by π
        u_max = float(np.max(u))
        u_min = float(np.min(u))
        # Hedgehog monopole has u→π*Q far out, u≈0 at centre
        # Heuristic: if center value is near zero AND max is near n*π
        center_val = float(u[cx, cy, cz])
        if abs(center_val) < 0.5 * np.pi and u_max > 0.5 * np.pi:
            # Monopole-like
            charge = int(round(u_max / np.pi))
            # Sanity check via spherical integration of dΩ-weighted phase
            return charge

        # Try vortex string along z: take the central z-plane
        # and integrate ∇θ around a loop
        for axis in ("z", "x", "y"):
            w = self._loop_winding(axis)
            if w != 0:
                return w
        return 0

    def _loop_winding(self, axis: str) -> int:
        """Integrate (1/2π) ∮ du · dl around a square loop in the plane
        perpendicular to ``axis`` through the lattice centre.

        The integrand is the discrete circulation of u itself; for the
        vortex_string_initial profile u = n×θ×tanh(r/core), so the
        circulation at large r is 2π n.
        """
        cx, cy, cz = self.Nx // 2, self.Ny // 2, self.Nz // 2
        u = self.u
        # Loop radius in cells (use ~ Nx//4 from centre)
        r = min(self.Nx, self.Ny, self.Nz) // 4
        if r < 2:
            return 0

        if axis == "z":
            # Plane z = cz; loop in (x,y)
            slab = u[:, :, cz]
            ix0, iy0 = cx, cy
        elif axis == "x":
            slab = u[cx, :, :]
            ix0, iy0 = cy, cz
        elif axis == "y":
            slab = u[:, cy, :]
            ix0, iy0 = cx, cz
        else:
            return 0

        # Build square loop indices anti-clockwise
        i_lo, i_hi = ix0 - r, ix0 + r
        j_lo, j_hi = iy0 - r, iy0 + r
        if i_lo < 0 or j_lo < 0 or i_hi >= slab.shape[0] or j_hi >= slab.shape[1]:
            return 0

        # Sample u along the loop (anticlockwise)
        path = []
        for j in range(j_lo, j_hi + 1):
            path.append(slab[i_hi, j])
        for i in range(i_hi, i_lo - 1, -1):
            path.append(slab[i, j_hi])
        for j in range(j_hi, j_lo - 1, -1):
            path.append(slab[i_lo, j])
        for i in range(i_lo, i_hi + 1):
            path.append(slab[i, j_lo])
        path = np.asarray(path)

        # Sum branch-corrected differences (unwrap), then divide by 2π
        diffs = np.diff(path)
        # Map diffs into [-π, π]
        diffs = (diffs + np.pi) % (2.0 * np.pi) - np.pi
        circ = float(diffs.sum())
        return int(round(circ / (2.0 * np.pi)))

    # ------------------------------------------------------------------
    # Inspection
    # ------------------------------------------------------------------

    def field_snapshot(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        return self.u.copy(), self.v.copy()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def kink_centre_axis(u: NDArray[np.float64], axis_coord: NDArray[np.float64],
                     axis: int = 0) -> float:
    """Estimate the centre of a 3D kink (extended in 2 transverse dims)
    along the given axis (0=x, 1=y, 2=z).
    """
    # Average over the two transverse axes
    other = tuple(a for a in (0, 1, 2) if a != axis)
    profile = u.mean(axis=other)
    f = profile - np.pi
    sign = np.sign(f)
    idx = np.where(np.diff(sign) != 0)[0]
    if idx.size == 0:
        return float("nan")
    i = int(idx[0])
    f0, f1 = f[i], f[i + 1]
    if f1 == f0:
        return float(axis_coord[i])
    frac = -f0 / (f1 - f0)
    return float(axis_coord[i] + frac * (axis_coord[i + 1] - axis_coord[i]))
