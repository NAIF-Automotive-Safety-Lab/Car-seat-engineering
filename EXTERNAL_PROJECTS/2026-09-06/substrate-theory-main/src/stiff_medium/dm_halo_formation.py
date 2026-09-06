"""Dark matter halo formation simulator (substrate cube-DM).

Cube-DM is EM-inert and stable, with mass m ~ 27.5 GeV. It interacts only
gravitationally so the macroscopic dynamics are identical to standard
cold collisionless dark matter. As a consequence, gravitational collapse of
an over-dense region of cube-DM should produce a halo with the standard NFW
density profile

    rho(r) = rho_s / [ (r / r_s) * (1 + r / r_s)**2 ]

with inner slope -1 and outer slope -3, and a concentration c = R_vir / r_s
that follows c(M) ~ M^{-0.1} for CDM-like halos.

This module provides a small, self-contained N-body code (direct summation
with a Plummer softening), a leapfrog integrator, a friends-of-friends-style
halo finder, a radial-profile builder, an NFW fit, and a c-M relation, plus
a virial-theorem diagnostic.

Units: we use code units in which G = 1. Length, mass and time can be mapped
to physical kpc / 1e10 M_sun / Gyr via the standard rescaling used in
galactic N-body simulations -- the validation tests here are scale-free.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple, List

import numpy as np
from scipy.optimize import curve_fit


G_CODE = 1.0  # gravitational constant in code units


# ----------------------------------------------------------------------------
# data containers
# ----------------------------------------------------------------------------


@dataclass
class Halo:
    """A halo identified in the particle distribution."""

    indices: np.ndarray
    center: np.ndarray
    radius: float          # virial-ish radius used by the finder
    mass: float            # sum of particle masses
    n_particles: int


# ----------------------------------------------------------------------------
# main simulator
# ----------------------------------------------------------------------------


class DMHaloFormation:
    """N-body dark-matter halo formation simulator.

    Direct O(N^2) gravity with Plummer softening. Adequate for N ~ 10^2-10^3,
    which is what we need for the unit tests. A `tree_force` stub is provided
    for future extension to larger N.
    """

    def __init__(self, softening: float = 0.05, particle_mass: float = 1.0,
                 seed: Optional[int] = None):
        self.softening = float(softening)
        self.particle_mass = float(particle_mass)
        self.rng = np.random.default_rng(seed)

        self.positions: Optional[np.ndarray] = None   # (N, 3)
        self.velocities: Optional[np.ndarray] = None  # (N, 3)
        self.masses: Optional[np.ndarray] = None      # (N,)
        self.time: float = 0.0
        self.history: List[float] = []

    # ------------------------------------------------------------------
    # initial conditions
    # ------------------------------------------------------------------

    def initialize_cosmological(self, N: int, box_size: float,
                                sigma_v: float) -> None:
        """Gaussian random positions in a centered box, isotropic Gaussian
        velocity field with dispersion sigma_v per axis. The mean velocity
        is removed so the centre of mass is at rest."""
        self.positions = self.rng.normal(0.0, box_size / 4.0, size=(N, 3))
        self.velocities = self.rng.normal(0.0, sigma_v, size=(N, 3))
        self.velocities -= self.velocities.mean(axis=0, keepdims=True)
        self.masses = np.full(N, self.particle_mass)
        self.time = 0.0
        self.history = []

    def initialize_halo(self, N: int, radius: float, sigma_v: float,
                        center: Optional[np.ndarray] = None) -> None:
        """Place N particles roughly inside a sphere of given radius with an
        isotropic velocity dispersion. Used to seed merger tests."""
        if center is None:
            center = np.zeros(3)
        # uniform-in-sphere via rejection
        pts = []
        while len(pts) < N:
            cand = self.rng.uniform(-radius, radius, size=(N, 3))
            mask = np.sum(cand ** 2, axis=1) < radius ** 2
            pts.extend(cand[mask].tolist())
        pos = np.asarray(pts[:N]) + center
        vel = self.rng.normal(0.0, sigma_v, size=(N, 3))
        vel -= vel.mean(axis=0, keepdims=True)
        if self.positions is None:
            self.positions = pos
            self.velocities = vel
            self.masses = np.full(N, self.particle_mass)
        else:
            self.positions = np.vstack([self.positions, pos])
            self.velocities = np.vstack([self.velocities, vel])
            self.masses = np.concatenate([self.masses,
                                          np.full(N, self.particle_mass)])
        self.time = 0.0

    # ------------------------------------------------------------------
    # forces
    # ------------------------------------------------------------------

    def compute_forces(self) -> np.ndarray:
        """Direct-sum gravity with Plummer softening. Returns acceleration
        array of shape (N, 3)."""
        pos = self.positions
        m = self.masses
        # pairwise displacement r_ij = r_j - r_i
        diff = pos[None, :, :] - pos[:, None, :]            # (N, N, 3)
        r2 = np.sum(diff ** 2, axis=-1) + self.softening ** 2
        inv_r3 = r2 ** (-1.5)
        np.fill_diagonal(inv_r3, 0.0)
        # acceleration on i = G * sum_j m_j (r_j - r_i) / |r_ij|^3
        acc = G_CODE * np.einsum('ijk,ij,j->ik', diff, inv_r3, m)
        return acc

    def tree_force(self) -> np.ndarray:
        """Stub for a tree (Barnes-Hut) force calculation. For the small N
        used by the tests we just dispatch to direct summation."""
        return self.compute_forces()

    # ------------------------------------------------------------------
    # integration
    # ------------------------------------------------------------------

    def evolve(self, t_max: float, dt: float = 0.01) -> None:
        """Leapfrog (kick-drift-kick) integration to time t_max."""
        if self.positions is None:
            raise RuntimeError("initial conditions not set")
        n_steps = max(1, int(np.ceil(t_max / dt)))
        dt = t_max / n_steps
        acc = self.compute_forces()
        for _ in range(n_steps):
            self.velocities += 0.5 * dt * acc
            self.positions += dt * self.velocities
            acc = self.compute_forces()
            self.velocities += 0.5 * dt * acc
            self.time += dt
            self.history.append(self.total_energy())

    # ------------------------------------------------------------------
    # diagnostics
    # ------------------------------------------------------------------

    def kinetic_energy(self, indices: Optional[np.ndarray] = None) -> float:
        if indices is None:
            v = self.velocities
            m = self.masses
        else:
            v = self.velocities[indices]
            m = self.masses[indices]
        return 0.5 * np.sum(m * np.sum(v ** 2, axis=1))

    def potential_energy(self, indices: Optional[np.ndarray] = None) -> float:
        if indices is None:
            pos = self.positions
            m = self.masses
        else:
            pos = self.positions[indices]
            m = self.masses[indices]
        diff = pos[None, :, :] - pos[:, None, :]
        r = np.sqrt(np.sum(diff ** 2, axis=-1) + self.softening ** 2)
        np.fill_diagonal(r, np.inf)
        u = -G_CODE * np.sum(np.triu(np.outer(m, m) / r, k=1))
        return u

    def total_energy(self) -> float:
        return self.kinetic_energy() + self.potential_energy()

    def velocity_dispersion(self, halo: Halo) -> float:
        """1-D velocity dispersion (averaged over the three axes) of the
        particles belonging to `halo`, in the halo's rest frame."""
        v = self.velocities[halo.indices]
        v = v - v.mean(axis=0, keepdims=True)
        return float(np.sqrt(np.mean(np.sum(v ** 2, axis=1)) / 3.0))

    # ------------------------------------------------------------------
    # halo finding
    # ------------------------------------------------------------------

    def halo_finder(self, sigma_threshold: float = 2.0,
                    min_particles: int = 20) -> List[Halo]:
        """Friends-of-friends-style finder. The mean inter-particle spacing is
        used to set a linking length, and `sigma_threshold` modulates how
        aggressively we group. Returns a list of Halo objects sorted by mass.
        """
        N = self.positions.shape[0]
        # rough mean separation from bounding box
        extent = self.positions.max(axis=0) - self.positions.min(axis=0)
        volume = max(np.prod(extent), 1e-12)
        mean_spacing = (volume / N) ** (1.0 / 3.0)
        link = mean_spacing / max(sigma_threshold, 1e-3)

        # union-find on pairs within `link`
        parent = list(range(N))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i, j):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj

        # vectorized neighbour search; OK at N <= a few thousand
        for i in range(N):
            d2 = np.sum((self.positions[i + 1:] - self.positions[i]) ** 2,
                        axis=1)
            neighbours = np.where(d2 < link ** 2)[0]
            for n in neighbours:
                union(i, i + 1 + int(n))

        groups: dict = {}
        for i in range(N):
            r = find(i)
            groups.setdefault(r, []).append(i)

        halos: List[Halo] = []
        for idxs in groups.values():
            if len(idxs) < min_particles:
                continue
            arr = np.asarray(idxs)
            pos = self.positions[arr]
            mass = float(np.sum(self.masses[arr]))
            center = (self.masses[arr][:, None] * pos).sum(axis=0) / mass
            r_max = float(np.max(np.linalg.norm(pos - center, axis=1)))
            halos.append(Halo(indices=arr, center=center, radius=r_max,
                              mass=mass, n_particles=len(arr)))
        halos.sort(key=lambda h: -h.mass)
        return halos

    # ------------------------------------------------------------------
    # profiles
    # ------------------------------------------------------------------

    def radial_profile(self, halo_center: np.ndarray, max_r: float,
                       n_bins: int = 20,
                       indices: Optional[np.ndarray] = None
                       ) -> Tuple[np.ndarray, np.ndarray]:
        """Build a logarithmic-radius density profile.

        Returns (r_centers, rho) where rho is the mass density in each shell.
        """
        if indices is None:
            pos = self.positions
            m = self.masses
        else:
            pos = self.positions[indices]
            m = self.masses[indices]

        r = np.linalg.norm(pos - halo_center, axis=1)
        r_min = max(self.softening, max_r / 1000.0)
        edges = np.logspace(np.log10(r_min), np.log10(max_r), n_bins + 1)
        rho = np.zeros(n_bins)
        centers = np.sqrt(edges[:-1] * edges[1:])
        for i in range(n_bins):
            mask = (r >= edges[i]) & (r < edges[i + 1])
            shell_mass = float(np.sum(m[mask]))
            shell_vol = (4.0 / 3.0) * np.pi * (edges[i + 1] ** 3 - edges[i] ** 3)
            rho[i] = shell_mass / shell_vol
        return centers, rho

    # ------------------------------------------------------------------
    # NFW fitting and c-M relation
    # ------------------------------------------------------------------

    @staticmethod
    def nfw_density(r: np.ndarray, rho_s: float, r_s: float) -> np.ndarray:
        x = r / r_s
        return rho_s / (x * (1.0 + x) ** 2)

    def nfw_fit(self, profile: Tuple[np.ndarray, np.ndarray]
                ) -> Tuple[float, float]:
        """Fit (rho_s, r_s) to a (r, rho) profile in log-log space."""
        r, rho = profile
        mask = (rho > 0) & np.isfinite(rho)
        r = r[mask]
        rho = rho[mask]
        if len(r) < 4:
            raise RuntimeError("not enough non-empty bins to fit NFW")

        log_rho = np.log(rho)

        def model(r_, log_rho_s, log_r_s):
            return np.log(self.nfw_density(r_, np.exp(log_rho_s),
                                           np.exp(log_r_s)))

        p0 = (np.log(np.median(rho)), np.log(np.median(r)))
        popt, _ = curve_fit(model, r, log_rho, p0=p0, maxfev=20000)
        rho_s = float(np.exp(popt[0]))
        r_s = float(np.exp(popt[1]))
        return rho_s, r_s

    @staticmethod
    def concentration(M_halo: float, M_pivot: float = 1e12,
                      c_pivot: float = 10.0, slope: float = -0.1) -> float:
        """Standard CDM c-M relation: c(M) = c_pivot * (M / M_pivot)^slope.

        Default values reproduce c ~ 10 for a Milky-Way-like 1e12 M_sun halo.
        """
        return float(c_pivot * (M_halo / M_pivot) ** slope)

    # ------------------------------------------------------------------
    # virial theorem
    # ------------------------------------------------------------------

    def virial_ratio(self, halo: Halo) -> float:
        """Return |2T + U| / |U|. For a virialized system this is small."""
        T = self.kinetic_energy(halo.indices)
        # subtract bulk motion of the halo before computing T
        v = self.velocities[halo.indices]
        m = self.masses[halo.indices]
        v_cm = (m[:, None] * v).sum(axis=0) / m.sum()
        v_rel = v - v_cm
        T = 0.5 * np.sum(m * np.sum(v_rel ** 2, axis=1))
        U = self.potential_energy(halo.indices)
        return float(abs(2.0 * T + U) / max(abs(U), 1e-30))


# ----------------------------------------------------------------------------
# convenience: build a relaxed halo for tests
# ----------------------------------------------------------------------------


def build_relaxed_halo(N: int = 200, radius: float = 1.0, M: float = 1.0,
                       seed: int = 0,
                       relax_time: float = 5.0,
                       softening: float = 0.05) -> DMHaloFormation:
    """Make a roughly virialised, NFW-like halo by initialising N particles
    drawn from an approximate NFW density and giving them the velocity
    dispersion implied by the virial theorem, then evolving briefly so that
    they relax."""
    sim = DMHaloFormation(softening=softening,
                          particle_mass=M / N, seed=seed)
    rng = sim.rng

    # sample radii from an approximate NFW M(<r) ~ ln(1+x) - x/(1+x) cdf
    # (we use a pivot scale-radius r_s = 0.3 * radius so c ~ 3.3, fine for tests)
    r_s = 0.3 * radius
    # rejection sample
    pts = []
    while len(pts) < N:
        r_try = rng.uniform(softening, radius, size=N)
        x = r_try / r_s
        # rho * 4 pi r^2 propto r / (1 + x)^2
        weight = r_try / (1.0 + x) ** 2
        weight /= weight.max()
        keep = rng.uniform(0.0, 1.0, size=N) < weight
        pts.extend(r_try[keep].tolist())
    r_samples = np.asarray(pts[:N])

    # random angles
    cos_th = rng.uniform(-1.0, 1.0, size=N)
    sin_th = np.sqrt(1.0 - cos_th ** 2)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=N)
    pos = np.column_stack([
        r_samples * sin_th * np.cos(phi),
        r_samples * sin_th * np.sin(phi),
        r_samples * cos_th,
    ])

    # virial velocity dispersion: sigma^2 ~ G M / (2 R)
    sigma = np.sqrt(G_CODE * M / (2.0 * radius))
    vel = rng.normal(0.0, sigma, size=(N, 3))
    vel -= vel.mean(axis=0, keepdims=True)

    sim.positions = pos
    sim.velocities = vel
    sim.masses = np.full(N, M / N)
    sim.evolve(relax_time, dt=0.05)
    return sim
