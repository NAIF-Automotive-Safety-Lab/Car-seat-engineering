"""substrate_3d_movie.py
========================

Full-evolution 3D substrate-field rendering pipeline.

This module is the *movie* counterpart to ``bound_state_3d_extractor``: it
evolves the §18.11 sine-Gordon Lagrangian on a 32³ cubic lattice,
collects snapshots at uniformly spaced times across one full breathing
period, and exposes those snapshots as 3D iso-surface point-clouds for
matplotlib rendering.

Where ``bound_state_3d_extractor`` extracts the rest-mass identity
m c² = ℏ ω_b from the breathing frequency, ``substrate_3d_movie`` shows
the *spatial* unfolding: u(x, y, z, t) is plotted at each of N_FRAMES
times so the bound state can be SEEN to breathe in place.

Components
----------

* ``Substrate3DGeometry`` — 32³ lattice geometry, bound-state seed
  (Gaussian bump matched to the ground mode of ``bound_state_3d_extractor``),
  CFL-safe time step.

* ``Substrate3DSimulator`` — evolves the full §18.11 Lagrangian with
  drag-aware velocity-Verlet, samples u(x,y,z,t) at N_FRAMES uniformly
  across one breathing period, and on demand extracts iso-surface
  point-clouds at any threshold |u| = level.

* Diagnostics — total energy E(t), centre-cell amplitude u(0,0,0,t),
  iso-surface volume V_iso(t) at the canonical breathing iso-level
  (counts cells whose |u| ≥ level · |u|_max).

The defaults are tuned so the test-suite finishes well under 60 s:
N = 32, L = 12 ξ, n_frames = 8, n_periods = 1.0.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from src.stiff_medium.lattice_substrate_3d import LatticeSubstrate3D
from src.stiff_medium.cone_bouncing_protocol import (
    omega_b_canonical,
    mass_from_omega,
    HBAR,
    C,
)


# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------


@dataclass
class Substrate3DGeometry:
    """32³ cubic-lattice geometry seeded with a bound-state bump.

    Parameters
    ----------
    K, rho, xi, gamma : float
        Substrate primitives.  Defaults are the natural baseline
        K = ρ = ξ = 1, γ = 0 so c_eff = 1, ω_b = 1, T = 2π.
    sigma_max : float
        Saturation cap (default 0.5 per §2.6).
    L : float
        Half-width of the cubic simulation box.  Default 12 ξ — wide
        enough that the Gaussian-bump bound state vanishes (<1%) at the
        boundary so periodic BC do not contaminate the breathing mode.
    N : int
        Grid points along each axis.  Default 32 (32 768 cells, fits in
        memory and finishes one period in ~10 s).
    """

    K: float = 1.0
    rho: float = 1.0
    xi: float = 1.0
    gamma: float = 0.0
    sigma_max: float = 0.5
    L: float = 12.0
    N: int = 32

    # ---- derived ----
    @property
    def c_eff(self) -> float:
        """Substrate wave speed √(K/ρ)."""
        return float(np.sqrt(self.K / self.rho))

    @property
    def dx(self) -> float:
        """Grid spacing 2L / N."""
        return float(2.0 * self.L / self.N)

    @property
    def dt(self) -> float:
        """3D-CFL-safe time step (safety factor 0.4)."""
        return float(0.4 * self.dx / (self.c_eff * np.sqrt(3.0)))

    @property
    def omega_b(self) -> float:
        """Canonical cone-bouncing angular frequency."""
        return omega_b_canonical(self.K, self.rho, self.xi, self.gamma)

    @property
    def period(self) -> float:
        """Breathing period T = 2π / ω_b."""
        omega = self.omega_b
        return float("inf") if omega <= 0 else float(2.0 * np.pi / omega)

    @property
    def envelope_width(self) -> float:
        """Width σ_env of the Gaussian bound-state envelope.  4 ξ is
        wide enough that |k| ≪ 1/ξ so the small-amplitude expansion of
        the cosine well linearises cleanly to ω₀ = c/ξ."""
        return float(4.0 * self.xi)

    # ---- seed ----

    def bound_state_seed(
        self,
        X: NDArray[np.float64],
        Y: NDArray[np.float64],
        Z: NDArray[np.float64],
        amplitude_frac: float = 0.05,
        breath_frac: float = 0.10,
    ) -> NDArray[np.float64]:
        """Breathing-mode Gaussian bound-state seed.

            u(r) = (1 + δ) · A · exp(−r² / (2 σ_env²))

        with A = amplitude_frac · σ_max and δ = breath_frac.  The
        breathing factor pushes the bump slightly off equilibrium so its
        centre-cell amplitude oscillates at ω_b without introducing
        wavevector content that would shift the FFT peak above ω₀.
        """
        if not 0.0 < amplitude_frac < 1.0:
            raise ValueError("amplitude_frac must lie in (0, 1).")
        if abs(breath_frac) >= 1.0:
            raise ValueError("|breath_frac| must be < 1.")
        A = amplitude_frac * self.sigma_max * (1.0 + breath_frac)
        r2 = X * X + Y * Y + Z * Z
        return A * np.exp(-r2 / (2.0 * self.envelope_width ** 2))


# ---------------------------------------------------------------------------
# 2. Simulator
# ---------------------------------------------------------------------------


@dataclass
class Substrate3DSimulator:
    """Evolve u(x,y,z,t) under the full §18.11 Lagrangian and sample
    snapshots for movie / iso-surface rendering.

    Parameters
    ----------
    geometry : Substrate3DGeometry
        Source of K, ρ, ξ, γ, σ_max + bound-state seed.
    n_frames : int
        Number of snapshots to record across the evolution window.
        Default 8 (matches the 4-snapshot strip + 6-panel iso threshold
        layout used by the renderer).
    n_periods : float
        Number of analytic breathing periods to evolve.  Default 1.0 so
        snapshots span exactly one full breathing cycle.
    amplitude_frac : float
        Initial bump amplitude as a fraction of σ_max.  Default 0.05 →
        A = 0.025, deeply inside the linear cosine well.
    breath_frac : float
        Initial breathing offset, default 0.10.
    """

    geometry: Substrate3DGeometry = field(default_factory=Substrate3DGeometry)
    n_frames: int = 8
    n_periods: float = 1.0
    amplitude_frac: float = 0.05
    breath_frac: float = 0.10

    # populated by run()
    _result: Optional[Dict[str, NDArray[np.float64]]] = field(
        default=None, init=False, repr=False)

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def _make_lattice(self) -> LatticeSubstrate3D:
        g = self.geometry
        return LatticeSubstrate3D(
            Nx=g.N, Ny=g.N, Nz=g.N,
            dx=g.dx, dt=g.dt,
            rho=g.rho, K=g.K, xi=g.xi, gamma=g.gamma,
            bc="periodic",
        )

    def seed(self, sub: LatticeSubstrate3D) -> None:
        """Apply the bound-state seed in place."""
        sub.u = self.geometry.bound_state_seed(
            sub.X, sub.Y, sub.Z,
            amplitude_frac=self.amplitude_frac,
            breath_frac=self.breath_frac,
        )
        sub.v.fill(0.0)

    # ------------------------------------------------------------------
    # Time evolution
    # ------------------------------------------------------------------

    def run(self) -> Dict[str, NDArray[np.float64]]:
        """Evolve and snapshot u(x,y,z,t) at n_frames uniformly spaced
        times over [0, n_periods · T_b].

        Returns a dict with keys

            t          : time samples [period units when natural]
            u_frames   : (n_frames, N, N, N) array of snapshots
            u_center   : centre-cell strain trace
            energy     : total Lagrangian energy at each frame
            iso_vol    : iso-surface volume (cell count) at each frame
                         using the canonical breathing iso-level
                         (35% of |u|_max over the run)
            x, y, z    : axis coordinates [units of ξ]
            T_period   : analytic breathing period
            iso_level  : the level used for ``iso_vol``
        """
        g = self.geometry
        sub = self._make_lattice()
        self.seed(sub)

        T_anal = g.period
        T_total = self.n_periods * T_anal
        n_samples = int(self.n_frames)
        if n_samples < 2:
            raise ValueError("n_frames must be >= 2.")
        dt_sample = T_total / (n_samples - 1) if n_samples > 1 else T_total
        steps_between = max(1, int(round(dt_sample / sub.dt)))

        t_arr = np.empty(n_samples, dtype=float)
        u_center = np.empty(n_samples, dtype=float)
        energy = np.empty(n_samples, dtype=float)
        u_frames = np.empty((n_samples, g.N, g.N, g.N), dtype=float)

        cx, cy, cz = sub.Nx // 2, sub.Ny // 2, sub.Nz // 2

        for k in range(n_samples):
            t_arr[k] = sub.t
            u_center[k] = float(sub.u[cx, cy, cz])
            energy[k] = float(sub.total_energy())
            u_frames[k] = sub.u.copy()
            if k < n_samples - 1:
                sub.evolve(steps_between)

        # Iso-volume series at the canonical breathing iso-level
        abs_max = float(np.abs(u_frames).max())
        if abs_max <= 0.0:
            abs_max = 1.0
        iso_level = 0.35 * abs_max
        iso_vol = np.array([
            int(np.sum(np.abs(u_frames[k]) >= iso_level))
            for k in range(n_samples)
        ], dtype=float)

        self._result = {
            "t": t_arr,
            "u_frames": u_frames,
            "u_center": u_center,
            "energy": energy,
            "iso_vol": iso_vol,
            "x": sub.x,
            "y": sub.y,
            "z": sub.z,
            "dt": np.array([sub.dt]),
            "dx": np.array([sub.dx]),
            "T_period": np.array([T_anal]),
            "iso_level": np.array([iso_level]),
            "abs_max": np.array([abs_max]),
        }
        return self._result

    # ------------------------------------------------------------------
    # Iso-surface extractor
    # ------------------------------------------------------------------

    @property
    def result(self) -> Dict[str, NDArray[np.float64]]:
        if self._result is None:
            self.run()
        assert self._result is not None
        return self._result

    def iso_points(
        self,
        frame_index: int,
        level_frac: float = 0.35,
        max_points: int = 4000,
        rng_seed: int = 0,
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64],
               NDArray[np.float64]]:
        """Extract an iso-surface as a coloured scatter point-cloud.

        Returns ``(xs, ys, zs, vals)`` where each array has length
        ≤ ``max_points`` and represents lattice cells whose |u| sits
        within a thin shell around ``level_frac * |u|_max``.

        Marching-cubes-free so it has no skimage dependency.

        Parameters
        ----------
        frame_index : int
            Which frame in ``u_frames`` to extract.
        level_frac : float
            Iso-level as a fraction of the run's max |u|.  Default 0.35.
        max_points : int
            Maximum number of scatter points to return.  If the raw
            mask exceeds this, points are subsampled deterministically
            using ``rng_seed``.
        rng_seed : int
            Subsampling seed (for reproducible visuals).
        """
        res = self.result
        u = res["u_frames"][int(frame_index)]
        x = res["x"]; y = res["y"]; z = res["z"]
        abs_max = float(res["abs_max"][0])
        if abs_max <= 0.0:
            abs_max = 1.0
        level = level_frac * abs_max
        # Thin shell: level <= |u| <= 1.5*level
        absu = np.abs(u)
        mask = (absu >= level) & (absu <= 1.5 * level)
        if mask.sum() == 0:
            # Relax to all cells above the level
            mask = absu >= level
        if mask.sum() > max_points:
            idx = np.where(mask.ravel())[0]
            chosen = np.random.default_rng(rng_seed).choice(
                idx, size=max_points, replace=False
            )
            mask = np.zeros(mask.size, dtype=bool)
            mask[chosen] = True
            mask = mask.reshape(absu.shape)
        ix, iy, iz = np.where(mask)
        xs = x[ix]
        ys = y[iy]
        zs = z[iz]
        vals = absu[mask]
        return xs, ys, zs, vals

    def diagnostics(self) -> Dict[str, float]:
        """Summary numbers used by the test suite + render labels."""
        res = self.result
        e = res["energy"]
        u_c = res["u_center"]
        iv = res["iso_vol"]
        # Energy conservation ratio
        e_init = float(e[0])
        e_drift = float(np.abs(e - e_init).max())
        e_drift_pct = 100.0 * e_drift / max(abs(e_init), 1e-30)
        # Iso volume oscillation amplitude
        iv_mean = float(iv.mean())
        iv_amp = float(iv.max() - iv.min())
        iv_amp_frac = iv_amp / max(iv_mean, 1e-30)
        return {
            "n_frames": int(len(res["t"])),
            "T_period": float(res["T_period"][0]),
            "iso_level": float(res["iso_level"][0]),
            "energy_init": e_init,
            "energy_max_drift": e_drift,
            "energy_drift_pct": e_drift_pct,
            "u_center_init": float(u_c[0]),
            "u_center_amp": float(np.std(u_c)),
            "iso_vol_mean": iv_mean,
            "iso_vol_amp": iv_amp,
            "iso_vol_amp_frac": iv_amp_frac,
            "omega_pred": float(self.geometry.omega_b),
            "rest_E_pred": float(HBAR * self.geometry.omega_b),
            "mass_pred_kg": float(mass_from_omega(self.geometry.omega_b)),
            "abs_max": float(res["abs_max"][0]),
        }
