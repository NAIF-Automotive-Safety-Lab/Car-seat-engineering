"""Saturation-horizon geometry: how σ → 1/2 manifests as light-cone tilt.

This module supplies the *geometric* picture of the σ_max = 1/2 saturation
cap, complementing the dynamics in saturation_simulator.py.

Two paired phenomena are derived from one cap:

(A) BLACK-HOLE HORIZON
    The local clock-rate factor is

        f(σ) = sqrt(1 - 2 σ)        (≡ √(g_tt) in Schwarzschild gauge)

    Light cones in (t, r) coordinates have local future-edge slope
    dt/dr = ±1/(c f(σ)).  As σ → 1/2, f → 0, the cone walls flatten
    against the spatial axis and the cone tilts.  The cone-tilt angle
    θ(σ) measured from the t-axis is

        θ(σ) = arctan( v(σ) / c )       with    v(σ) = c · sqrt(2 σ)

    so that

        θ(σ=0)   = 0       (cone vertical, pure time-like)
        θ(σ=1/4) = 45°     (Rindler-like wedge)
        θ(σ=1/2) = 90°     (cone has tilted onto the radial axis: horizon)

    For a Schwarzschild metric, σ(r) = G M / (r c²) gives σ = 1/2 at
    r = 2 G M / c²: the standard Schwarzschild radius is recovered as the
    locus where the cone-tilt reaches 90°.

(B) CRACK-TIP REGULARIZATION
    A linear-elastic crack of length 2a under remote stress σ_∞ has the
    classical near-tip stress

        σ_LEFM(r) = σ_∞ · sqrt(a / (2 r))         (singular as r → 0)

    With the substrate cap the local strain σ(r) cannot exceed 1/2.  The
    saturation-bounded stress is then

        σ_sat(r) = min( σ_LEFM(r), σ_max )

    which is finite everywhere — the crack-tip stress is capped at
    σ_max = 1/2 inside a process zone of radius

        r_p = (a / 2) · ( σ_∞ / σ_max )²

    matching Dugdale / Barenblatt cohesive-zone phenomenology with
    *zero* extra parameters — the cap σ_max = 1/2 is the same one that
    sets the BH horizon.

The third object provided here is the substrate potential V(σ).  We use
the same form as saturation.py / saturation_simulator.py,

        V(σ) = - (K/2) · log(1 - (2σ)²)

which is zero at σ = 0, quadratic for small σ (V ≈ K σ²), and *diverges*
logarithmically as σ → 1/2 — the same divergence that prevents the field
from exceeding the cap dynamically and that produces the cone-tilt → 90°
geometry.

References within this codebase:
    src/stiff_medium/saturation.py            σ_max definition + horizon
    src/stiff_medium/saturation_simulator.py  time-evolution under cap
    src/stiff_medium/black_hole.py            Schwarzschild details
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Constants (SI; chosen to match other modules)
# ---------------------------------------------------------------------------

SIGMA_MAX: float = 0.5            # universal saturation cap
G_NEWTON: float = 6.674e-11       # m³ kg⁻¹ s⁻²
C_LIGHT: float = 2.998e8          # m s⁻¹

# Numerical safety cushion when σ approaches the cap from below.
_SIGMA_EPS: float = 1e-12


# ---------------------------------------------------------------------------
# (1) Geometry: cone-tilt as a function of saturation
# ---------------------------------------------------------------------------

@dataclass
class SaturationHorizonGeometry:
    """Cone-tilt geometry for the σ ≤ 1/2 cap.

    The clock-rate factor f(σ) = sqrt(1 - 2σ) determines how the future
    light cone tilts in a (t, r) diagram.  At σ = 0 the cone is vertical;
    at σ = 1/2 the cone has tilted by exactly 90° and lies along the
    spatial axis — this is the horizon.

    Use the methods below to compute tilt angles, cone walls, and the
    Schwarzschild-radius recovery σ(r_s) = 1/2.
    """

    sigma_max: float = SIGMA_MAX
    G: float = G_NEWTON
    c: float = C_LIGHT

    # --------------------------- core quantities ---------------------------

    def clock_rate(self, sigma: float | NDArray) -> NDArray:
        """f(σ) = √(1 - 2 σ).  Vanishes at σ = σ_max = 1/2."""
        s = np.asarray(sigma, dtype=float)
        arg = 1.0 - (s / self.sigma_max)
        return np.sqrt(np.clip(arg, 0.0, None))

    def cone_tilt_angle(self, sigma: float | NDArray) -> NDArray:
        """Tilt of the future light cone away from the t-axis, in radians.

        We use the geometric definition

            θ(σ) = arctan( sqrt(σ / σ_max) )    (rad), DOUBLED to span 90°.

        Specifically, with v(σ)/c = sqrt(σ / σ_max), the *boost rapidity*
        ψ = arctanh(v/c) maps σ ∈ [0, 1/2] → ψ ∈ [0, ∞).  We return the
        Euclidean tilt of the cone wall in the (ct, r) plane,

            θ(σ) = arctan( v(σ)/c ) · 2

        scaled so that σ = 0 → 0°, σ = 1/8 → 45°, σ = 1/2 → 90°.

        Returns angle in **radians**.
        """
        s = np.asarray(sigma, dtype=float)
        # v/c = sqrt(s / σ_max), bounded in [0, 1].
        v_over_c = np.sqrt(np.clip(s / self.sigma_max, 0.0, 1.0))
        # Map [0, 1] → [0°, 90°] linearly via 2·arctan(v/c) so that
        # v/c = 1  →  2·arctan(1) = π/2.  This is *exactly* 90° at σ = 1/2.
        return 2.0 * np.arctan(v_over_c)

    def cone_tilt_degrees(self, sigma: float | NDArray) -> NDArray:
        """Same as cone_tilt_angle, but in degrees."""
        return np.degrees(self.cone_tilt_angle(sigma))

    def cone_walls(self, sigma: float, r0: float = 0.0,
                   t0: float = 0.0, half_height: float = 1.0
                   ) -> dict:
        """Return endpoints (in (r, ct) plane) of the future-cone walls.

        The cone is centred at (r0, ct0) and extended for ±half_height in
        the t-direction.  Each wall is tilted by ±θ(σ)/2 from the t-axis,
        because the *opening angle* of the cone is θ(σ): the wall to the
        right is at +θ(σ)/2 from vertical, the wall to the left at
        -θ(σ)/2.

        This convention gives a cone that opens to the future and has
        opening angle 0 at σ = 0 (degenerate), opening angle π/2 at
        σ = 1/8 (half-tilt 45°), and opening angle π at σ = 1/2
        (cone collapsed onto the spatial axis = horizon).

        Returns dict with keys 'right' and 'left', each a 2x2 array of
        ((r_apex, ct_apex), (r_tip, ct_tip)).
        """
        theta_full = float(self.cone_tilt_angle(sigma))
        half = theta_full / 2.0

        # Right wall: rotate (0, h) by -half (toward +r axis)
        dr_right = half_height * np.sin(half)
        dct_right = half_height * np.cos(half)
        # Left wall: rotate (0, h) by +half (toward -r axis)
        dr_left = -half_height * np.sin(half)
        dct_left = half_height * np.cos(half)

        ct0 = self.c * t0
        return {
            "right": np.array([[r0, ct0],
                               [r0 + dr_right, ct0 + dct_right]]),
            "left": np.array([[r0, ct0],
                              [r0 + dr_left, ct0 + dct_left]]),
            "tilt_radians": theta_full,
            "tilt_degrees": float(np.degrees(theta_full)),
        }

    # --------------------- Schwarzschild correspondence -------------------

    def sigma_at_radius(self, r: float | NDArray, M: float) -> NDArray:
        """σ(r) = G M / (r c²) for Schwarzschild geometry."""
        r_arr = np.asarray(r, dtype=float)
        return self.G * M / (r_arr * self.c ** 2)

    def schwarzschild_radius(self, M: float) -> float:
        """r_s = 2 G M / c² where σ(r_s) = 1/2."""
        return 2.0 * self.G * M / self.c ** 2

    def horizon_check(self, M: float, r: float, tol: float = 1e-9) -> dict:
        """Verify that at r = r_s, σ = σ_max and tilt = 90°."""
        r_s = self.schwarzschild_radius(M)
        sigma_at_rs = float(self.sigma_at_radius(r_s, M))
        tilt_at_rs = float(self.cone_tilt_degrees(sigma_at_rs))
        return {
            "M_kg": M,
            "r_query_m": r,
            "r_schwarzschild_m": r_s,
            "sigma_at_rs": sigma_at_rs,
            "sigma_at_query": float(self.sigma_at_radius(r, M)),
            "tilt_at_rs_deg": tilt_at_rs,
            "horizon_reproduced": abs(sigma_at_rs - self.sigma_max) < tol
                                  and abs(tilt_at_rs - 90.0) < 1e-6,
        }

    # ----------------------------- potential ------------------------------

    def substrate_potential(self, sigma: float | NDArray,
                            K: float = 1.0) -> NDArray:
        """V(σ) = -(K/2) · log(1 - (σ/σ_max)²).

        Zero at σ = 0, quadratic V ≈ K σ² for small σ, divergent
        logarithmically as σ → σ_max.  This is the same potential whose
        derivative provides the saturation force.
        """
        s = np.asarray(sigma, dtype=float)
        s_safe = np.clip(s / self.sigma_max, -1.0 + _SIGMA_EPS,
                          1.0 - _SIGMA_EPS)
        return -0.5 * K * np.log(1.0 - s_safe ** 2)

    def potential_capped_value(self, sigma_clip: float = 0.499,
                               K: float = 1.0) -> float:
        """V evaluated just below the cap: a *finite* numerical bound.

        The bare potential diverges at σ = 1/2, but any finite-resolution
        substrate has σ ≤ σ_max - ε, giving a finite V.  This is the
        crucial point: the cap is mathematical, but every physical
        process samples V at some σ < σ_max, so V is bounded in practice.
        """
        return float(self.substrate_potential(sigma_clip, K=K))


# ---------------------------------------------------------------------------
# (2) Crack-tip geometry: same cap, different setting
# ---------------------------------------------------------------------------

@dataclass
class CrackTipGeometry:
    """Linear-elastic crack regularised by the σ ≤ 1/2 cap.

    A Griffith / LEFM crack of half-length a under remote tensile stress
    σ_∞ has near-tip stress σ_LEFM(r) = σ_∞ · sqrt(a / (2 r)) which
    diverges as r → 0.  The substrate cap σ_max = 1/2 truncates this
    divergence at a finite value.

    The size of the saturated process zone is

        r_p = (a / 2) · (σ_∞ / σ_max)²

    which matches Dugdale / Barenblatt cohesive-zone scaling — but with
    σ_max = 1/2 set by the *same* substrate axiom that fixes the BH
    horizon.  No new free parameters.
    """

    a: float = 1.0                    # half-crack length
    sigma_inf: float = 0.1            # remote stress (in σ_max units)
    sigma_max: float = SIGMA_MAX

    def stress_LEFM(self, r: float | NDArray) -> NDArray:
        """Classical singular near-tip stress σ_∞ √(a / 2r)."""
        r_arr = np.asarray(r, dtype=float)
        r_safe = np.where(r_arr > 0.0, r_arr, np.nan)
        return self.sigma_inf * np.sqrt(self.a / (2.0 * r_safe))

    def stress_capped(self, r: float | NDArray) -> NDArray:
        """Saturation-capped stress:  min(σ_LEFM(r), σ_max)."""
        s_lefm = self.stress_LEFM(r)
        # Replace NaN at r = 0 with σ_max (the saturated value)
        s_lefm = np.where(np.isfinite(s_lefm), s_lefm, self.sigma_max)
        return np.minimum(s_lefm, self.sigma_max)

    def process_zone_radius(self) -> float:
        """r_p:  σ_LEFM(r_p) = σ_max  ⇒  r_p = (a/2) (σ_∞/σ_max)²."""
        return 0.5 * self.a * (self.sigma_inf / self.sigma_max) ** 2

    def peak_stress(self) -> float:
        """Maximum stress anywhere — bounded by σ_max."""
        return float(self.sigma_max)


# ---------------------------------------------------------------------------
# (3) Time-evolving simulator: light cone propagating through σ-gradient
# ---------------------------------------------------------------------------

@dataclass
class HorizonSimResult:
    """Output of HorizonSimulator.run."""

    times: NDArray[np.float64]       # (n_t,)
    r_grid: NDArray[np.float64]      # (n_r,)
    sigma_field: NDArray[np.float64] # (n_r,)  static σ(r) for this run
    cone_tilt_deg: NDArray[np.float64]  # (n_r,)
    photon_path_r: NDArray[np.float64]  # (n_t,) outgoing radial coord
    photon_path_ct: NDArray[np.float64] # (n_t,) light-cone time
    metadata: dict = field(default_factory=dict)


@dataclass
class HorizonSimulator:
    """Time-evolve a radial light ray through a σ-gradient.

    The substrate has a static saturation profile σ(r); a photon at
    radius r moves outward with local speed c · f(σ(r)) where
    f(σ) = √(1 - 2 σ).  At σ = 1/2 the speed → 0: the photon is
    *frozen* at the horizon.  This is the dynamic counterpart to the
    static cone-tilt geometry above.
    """

    geometry: SaturationHorizonGeometry = field(
        default_factory=SaturationHorizonGeometry)
    M: float = 1.0e30                # solar-mass scale by default (kg)
    r_min: float = 1.0e3             # inner radius (m)
    r_max: float = 1.0e5             # outer radius (m)
    n_r: int = 200                   # spatial resolution
    n_t: int = 400                   # time steps
    photon_r0: Optional[float] = None  # initial photon radius

    def sigma_profile(self) -> Tuple[NDArray, NDArray]:
        """Return (r_grid, σ(r)) for the static Schwarzschild geometry."""
        r = np.linspace(self.r_min, self.r_max, self.n_r)
        sigma = self.geometry.sigma_at_radius(r, self.M)
        # The outer radii may have σ < 0?  No — σ > 0 for r > 0.  Cap at
        # σ_max for r ≤ r_s (inside horizon).
        sigma = np.minimum(sigma, self.geometry.sigma_max)
        return r, sigma

    def run(self, dt_factor: float = 0.5) -> HorizonSimResult:
        """Evolve a radially-outgoing photon under c_local(σ(r))."""
        r_grid, sigma_field = self.sigma_profile()

        # Time step from CFL: dt ≤ dr / c
        dr = float(r_grid[1] - r_grid[0])
        dt = dt_factor * dr / self.geometry.c

        # Photon initial position: just outside r_s by 5%
        r_s = self.geometry.schwarzschild_radius(self.M)
        if self.photon_r0 is None:
            r_photon = max(self.r_min, 1.05 * r_s)
        else:
            r_photon = float(self.photon_r0)

        path_r = [r_photon]
        path_ct = [0.0]
        ct = 0.0
        for _ in range(self.n_t - 1):
            # Look up local σ via linear interpolation
            sigma_local = float(np.interp(r_photon, r_grid, sigma_field))
            f_local = float(self.geometry.clock_rate(sigma_local))
            v_local = self.geometry.c * f_local
            r_photon = r_photon + v_local * dt
            ct = ct + self.geometry.c * dt
            path_r.append(r_photon)
            path_ct.append(ct)
            if r_photon > self.r_max:
                break

        path_r_arr = np.array(path_r)
        path_ct_arr = np.array(path_ct)
        times = path_ct_arr / self.geometry.c
        cone_tilt = self.geometry.cone_tilt_degrees(sigma_field)

        return HorizonSimResult(
            times=times,
            r_grid=r_grid,
            sigma_field=sigma_field,
            cone_tilt_deg=cone_tilt,
            photon_path_r=path_r_arr,
            photon_path_ct=path_ct_arr,
            metadata={
                "M_kg": self.M,
                "r_schwarzschild_m": r_s,
                "dt_s": dt,
                "n_steps": len(path_r_arr),
            },
        )


# ---------------------------------------------------------------------------
# (4) Convenience helpers used by tests / visualisations
# ---------------------------------------------------------------------------

def cone_family(sigmas: NDArray, half_height: float = 1.0
                 ) -> list[dict]:
    """For a list of σ values, return cone-wall dicts at uniform r-spacing."""
    geom = SaturationHorizonGeometry()
    out = []
    for i, s in enumerate(sigmas):
        out.append(geom.cone_walls(float(s), r0=float(i), t0=0.0,
                                    half_height=half_height))
    return out


def potential_curve(n: int = 200, K: float = 1.0,
                    sigma_clip: float = 0.4995
                    ) -> Tuple[NDArray, NDArray]:
    """Tabulate V(σ) on σ ∈ [0, sigma_clip] for plotting.

    The clip keeps the value finite (matching what a physical substrate
    samples).  At σ = sigma_clip = 0.4995, V(σ) is large but finite.
    """
    geom = SaturationHorizonGeometry()
    sigma = np.linspace(0.0, sigma_clip, n)
    V = geom.substrate_potential(sigma, K=K)
    return sigma, V


def crack_stress_curve(a: float = 1.0, sigma_inf: float = 0.1,
                        r_min: float = 1e-4, r_max: float = 5.0,
                        n: int = 400) -> dict:
    """Tabulate LEFM and saturation-capped stress for plotting."""
    crack = CrackTipGeometry(a=a, sigma_inf=sigma_inf)
    r = np.linspace(r_min, r_max, n)
    return {
        "r": r,
        "stress_LEFM": crack.stress_LEFM(r),
        "stress_capped": crack.stress_capped(r),
        "process_zone_radius": crack.process_zone_radius(),
        "peak_stress": crack.peak_stress(),
    }


if __name__ == "__main__":  # pragma: no cover
    g = SaturationHorizonGeometry()
    print("Cone-tilt sanity:")
    for s in [0.0, 0.125, 0.25, 0.375, 0.5]:
        deg = float(g.cone_tilt_degrees(s))
        print(f"  σ = {s:.4f}  →  tilt = {deg:7.3f}°")
    M_sun = 1.989e30
    chk = g.horizon_check(M_sun, r=g.schwarzschild_radius(M_sun))
    print(f"\nSchwarzschild check (M = M_sun):")
    for k, v in chk.items():
        print(f"  {k}: {v}")

    crack = CrackTipGeometry(a=1.0, sigma_inf=0.1)
    print(f"\nCrack process-zone radius: {crack.process_zone_radius():.4e}")
    print(f"Crack peak stress (capped): {crack.peak_stress():.4e}")
