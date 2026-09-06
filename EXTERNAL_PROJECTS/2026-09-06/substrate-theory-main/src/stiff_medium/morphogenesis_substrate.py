"""Morphogenesis / Turing patterns from the substrate framework — paired GEOMETRY + SIM.

Substrate ontology of "morphogenesis"
-------------------------------------
In the B3 substrate picture a developing tissue is a 2D sheet of the
stiff medium where each cell is a stationary strain pocket. Two
substrate "morphogen" fields propagate as long-wavelength strain waves
on top of the sheet:

  * an *activator* u (compressive strain — autocatalytic, short-range
    diffusion D_u),
  * an *inhibitor* v (tensile strain — counter-autocatalytic,
    long-range diffusion D_v).

The local nonlinearity (Gierer-Meinhardt 1972) is

    f(u, v) = ρ · (u² / v − μ · u)
    g(u, v) = ρ · u² − ν · v

i.e. the substrate's nonlinear response ("reaction") feeds the
activator off itself and limits it via the inhibitor; diffusion
("substrate wave propagation") moves the two fields at very different
speeds. Turing's 1952 theorem: when D_v / D_u exceeds a critical
ratio, the homogeneous steady state is unstable and a stationary
spatial pattern with intrinsic wavelength

    λ ~ 2π · √(D_v / k)        (k = curvature of the linearised reaction term)

emerges spontaneously — exactly the leopard spots / zebra stripes /
Belousov-Zhabotinsky maze that the experimental record shows.

The same picture also covers Wolpert's *positional information*: a
1D substrate strain gradient set up between a source (left edge) and
a sink (right edge) divides the tissue into stripes via concentration
*thresholds* — the canonical "French flag" model.

This module wraps the textbook reaction-diffusion equations in two
paired classes:

  * ``MorphogenesisGeometry`` — 2D tissue (Lx, Ly, n_x, n_y) plus
    morphogen concentration field initialisers (random, gradient,
    central-spike) and the predicted Turing wavelength
    λ ≈ 2π √(D_v / k).

  * ``MorphogenesisSimulator`` — explicit-Euler integrator for the
    coupled Gierer-Meinhardt + diffusion PDE on a periodic grid with
    the 5-point Laplacian; pattern classifier (stripe / spot / maze)
    based on the Fourier-amplitude axial asymmetry; French-flag
    gradient solver and three-band cell-fate readout (blue / white /
    red bands of the Wolpert flag).

Public helpers
--------------
``turing_wavelength(D_v, k)``               — λ = 2π · √(D_v / k)
``turing_instability_satisfied(D_u, D_v, ratio_crit)`` — D_v/D_u > critical
``french_flag_fates(c, theta_low, theta_high)`` — band assignment

References
----------
* Turing (1952)              — "The chemical basis of morphogenesis", Phil. Trans. R. Soc. B 237, 37.
* Gierer & Meinhardt (1972)  — "A theory of biological pattern formation", Kybernetik 12, 30.
* Wolpert (1969)             — "Positional information and the spatial pattern of cellular differentiation", J. Theor. Biol. 25, 1.
* Murray, Mathematical Biology II (Springer 2003) chs 2-3 (RD textbook).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import math

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Module-level reference parameters (literature Gierer-Meinhardt regimes)
# ---------------------------------------------------------------------------

# Canonical reaction-diffusion parameter sets that produce the three
# textbook patterns: stripes (anisotropic), spots (bumpy), and labyrinths
# / mazes (isotropic, high gain). Values from Murray (2003) and the
# Meinhardt review (Models of Biological Pattern Formation, 1982).
LITERATURE: Dict[str, Dict[str, float]] = {
    # Stripes: D_v/D_u ≈ 40, moderate gain — leopard markings, cuttlefish
    # display, finger ridges. Activator scale << inhibitor scale.
    "stripes_zebra": {
        "D_u": 1.0e-5,
        "D_v": 4.0e-4,
        "rho": 0.05,
        "mu": 0.10,
        "nu": 0.10,
        "label": "stripes (zebra / cuttlefish / fingerprints)",
    },
    # Spots: high inhibitor diffusion, lower gain — leopard spots,
    # giraffe patches. Same kinetics, different D ratio.
    "spots_leopard": {
        "D_u": 1.0e-5,
        "D_v": 1.0e-3,
        "rho": 0.05,
        "mu": 0.10,
        "nu": 0.10,
        "label": "spots (leopard / giraffe)",
    },
    # Maze / labyrinth: even higher gain produces the BZ-style
    # convoluted bands seen in cerebral cortex folding & BZ chemistry.
    "maze_belousov": {
        "D_u": 1.0e-5,
        "D_v": 2.5e-4,
        "rho": 0.10,
        "mu": 0.10,
        "nu": 0.10,
        "label": "maze / labyrinth (cortex / BZ chemistry)",
    },
}


# ---------------------------------------------------------------------------
# Module-level functional helpers
# ---------------------------------------------------------------------------

def turing_wavelength(D_v: float, k: float) -> float:
    """Predicted Turing pattern wavelength λ = 2π · √(D_v / k).

    Parameters
    ----------
    D_v : inhibitor diffusion coefficient (length²/time).
    k : effective reaction-rate curvature at the homogeneous fixed
        point (1/time). For Gierer-Meinhardt with ρ ~ k_B T-scale
        nonlinearity, k ≈ ρ at first order.

    The expression is the standard linear-stability result: the
    fastest-growing Turing mode has wavenumber q* = (k / D_v)^(1/4)
    in the symmetric activator/inhibitor case; the corresponding
    wavelength is λ = 2π / q* ∝ √(D_v / k) up to an O(1) prefactor.
    """
    if D_v <= 0.0:
        raise ValueError("Inhibitor diffusion D_v must be positive.")
    if k <= 0.0:
        raise ValueError("Reaction-rate curvature k must be positive.")
    return float(2.0 * math.pi * math.sqrt(D_v / k))


def turing_instability_satisfied(
    D_u: float, D_v: float, ratio_critical: float = 10.0,
) -> bool:
    """Return True iff D_v / D_u exceeds the Turing instability threshold.

    For the classical Gierer-Meinhardt activator-inhibitor system the
    bifurcation occurs when the inhibitor diffuses faster than the
    activator by a factor of order 10-30 (the exact value depends on
    ρ, μ, ν but is always > 1). The default ``ratio_critical = 10``
    captures the textbook "short-range activation / long-range
    inhibition" regime.
    """
    if D_u <= 0.0:
        raise ValueError("Activator diffusion D_u must be positive.")
    if D_v <= 0.0:
        raise ValueError("Inhibitor diffusion D_v must be positive.")
    if ratio_critical <= 1.0:
        raise ValueError("Critical ratio must exceed 1 (short < long range).")
    return bool((D_v / D_u) > ratio_critical)


def french_flag_fates(
    c: NDArray[np.float64],
    theta_low: float = 0.33,
    theta_high: float = 0.66,
) -> NDArray[np.int_]:
    """Assign cell fates from morphogen concentration via two thresholds.

    Wolpert's French-flag model: a 1D morphogen gradient running from
    high (source) to low (sink) divides the tissue into three bands:

        c >= theta_high          → fate 2 (red — high morphogen)
        theta_low <= c < theta_high → fate 1 (white — middle band)
        c < theta_low            → fate 0 (blue — low morphogen)

    Parameters
    ----------
    c : concentration array (any shape).
    theta_low, theta_high : threshold values; must satisfy 0 ≤ low < high ≤ max.

    Returns
    -------
    Integer array of the same shape as c with fate codes 0, 1, 2.
    """
    if theta_low >= theta_high:
        raise ValueError("Lower threshold must be strictly less than upper.")
    arr = np.asarray(c, dtype=np.float64)
    fates = np.zeros(arr.shape, dtype=np.int_)
    fates = np.where(arr >= theta_high, 2, fates)
    fates = np.where((arr >= theta_low) & (arr < theta_high), 1, fates)
    fates = np.where(arr < theta_low, 0, fates)
    return fates


# ---------------------------------------------------------------------------
# 1. MorphogenesisGeometry — 2D tissue + morphogen field
# ---------------------------------------------------------------------------

@dataclass
class MorphogenesisGeometry:
    """2D tissue substrate with a morphogen concentration field.

    The tissue is a rectangle (Lx, Ly) discretised on an n_x × n_y
    grid. The morphogen "concentration" is the substrate strain
    amplitude — high u corresponds to a compressive strain pocket
    (e.g. nascent feather follicle), high v to a tensile relief field
    that travels far and damps the activator.

    Attributes
    ----------
    Lx, Ly : tissue extents (length units, e.g. mm).
    n_x, n_y : number of grid points along each axis.
    D_u, D_v : activator / inhibitor diffusion coefficients.
    rho : autocatalysis amplitude (sets ρ in Gierer-Meinhardt).
    mu : activator linear-decay rate.
    nu : inhibitor linear-decay rate.
    label : human-readable name of the regime.
    """

    Lx: float = 1.0
    Ly: float = 1.0
    n_x: int = 64
    n_y: int = 64
    D_u: float = 1.0e-5
    D_v: float = 4.0e-4
    rho: float = 0.05
    mu: float = 0.10
    nu: float = 0.10
    label: str = "generic activator-inhibitor tissue"

    # ------------------------------------------------------------------ grid helpers
    def dx(self) -> float:
        """Cell spacing along x: Δx = Lx / n_x."""
        return float(self.Lx / self.n_x)

    def dy(self) -> float:
        """Cell spacing along y: Δy = Ly / n_y."""
        return float(self.Ly / self.n_y)

    def coordinate_grids(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return (X, Y) meshgrids of cell centre coordinates."""
        x = (np.arange(self.n_x) + 0.5) * self.dx()
        y = (np.arange(self.n_y) + 0.5) * self.dy()
        return np.meshgrid(x, y, indexing="xy")

    # ------------------------------------------------------------------ field initialisers
    def homogeneous_field(self, value: float = 1.0) -> NDArray[np.float64]:
        """Uniform field of given value over the tissue."""
        return np.full((self.n_y, self.n_x), float(value), dtype=np.float64)

    def random_perturbation_field(
        self,
        base: float = 1.0,
        amplitude: float = 0.05,
        rng_seed: Optional[int] = 0,
    ) -> NDArray[np.float64]:
        """Uniform field + small white-noise perturbation (Turing seed).

        Turing patterns require an arbitrarily small perturbation off
        the homogeneous fixed point to start growing; a default
        amplitude of 5% is well inside the linear-instability regime.
        """
        rng = np.random.default_rng(rng_seed)
        noise = rng.standard_normal(size=(self.n_y, self.n_x))
        return base + amplitude * noise

    def gradient_field_1d(
        self,
        c_source: float = 1.0,
        c_sink: float = 0.0,
    ) -> NDArray[np.float64]:
        """Linear 1D morphogen gradient from x = 0 (source) to x = Lx (sink).

        Constant along y; this is the standard Wolpert geometry that
        generates positional information for the French-flag readout.
        """
        x = np.linspace(0.0, 1.0, self.n_x)
        line = c_source + (c_sink - c_source) * x
        return np.tile(line, (self.n_y, 1)).astype(np.float64)

    def central_spike_field(
        self,
        base: float = 1.0,
        spike: float = 5.0,
        radius_fraction: float = 0.05,
    ) -> NDArray[np.float64]:
        """Uniform field + a Gaussian spike at the centre (organiser cell).

        Useful for showing the spreading of an inhibitor wave from a
        localised source — substrate "wave propagation" from a single
        compressive seed.
        """
        X, Y = self.coordinate_grids()
        x0 = 0.5 * self.Lx
        y0 = 0.5 * self.Ly
        sigma = radius_fraction * 0.5 * (self.Lx + self.Ly)
        bump = (spike - base) * np.exp(
            -((X - x0) ** 2 + (Y - y0) ** 2) / (2.0 * sigma ** 2)
        )
        return base + bump.astype(np.float64)

    # ------------------------------------------------------------------ Turing prediction
    def predicted_wavelength(self) -> float:
        """λ ≈ 2π · √(D_v / ρ) — the dominant Turing mode wavelength.

        We identify the reaction-rate curvature k with ρ (the
        autocatalysis amplitude), which is the leading nonlinearity
        scale at the homogeneous fixed point.
        """
        return turing_wavelength(self.D_v, self.rho)

    def diffusion_ratio(self) -> float:
        """D_v / D_u — must exceed ~10 for Turing instability."""
        if self.D_u <= 0.0:
            raise ValueError("D_u must be positive.")
        return float(self.D_v / self.D_u)

    def turing_unstable(self, ratio_critical: float = 10.0) -> bool:
        """True iff the parameters lie in the Turing-unstable regime."""
        return turing_instability_satisfied(self.D_u, self.D_v, ratio_critical)

    # ------------------------------------------------------------------ summary
    def summary(self) -> Dict[str, float]:
        return {
            "label": self.label,
            "Lx": self.Lx,
            "Ly": self.Ly,
            "n_x": self.n_x,
            "n_y": self.n_y,
            "D_u": self.D_u,
            "D_v": self.D_v,
            "D_v_over_D_u": self.diffusion_ratio(),
            "rho": self.rho,
            "mu": self.mu,
            "nu": self.nu,
            "predicted_wavelength": self.predicted_wavelength(),
            "turing_unstable": self.turing_unstable(),
        }


# ---------------------------------------------------------------------------
# 2. MorphogenesisSimulator — Gierer-Meinhardt RD on the geometry
# ---------------------------------------------------------------------------

@dataclass
class MorphogenesisSimulator:
    """Gierer-Meinhardt activator-inhibitor RD integrator on a periodic grid.

    Discretises

        ∂u/∂t = D_u ∇²u + ρ (u²/v − μ u)
        ∂v/∂t = D_v ∇²v + ρ u² − ν v

    with the 5-point Laplacian on the geometry grid (periodic
    boundaries — the tissue is a torus). Time step dt is chosen
    automatically from the diffusion CFL condition unless overridden.

    Attributes
    ----------
    geometry : MorphogenesisGeometry providing Lx/Ly/n_x/n_y and parameters.
    dt : explicit-Euler time step (auto-set if 0).
    n_steps : number of integration steps.
    rng_seed : seed for the initial random perturbation.
    """

    geometry: MorphogenesisGeometry = field(default_factory=MorphogenesisGeometry)
    dt: float = 0.0
    n_steps: int = 4000
    rng_seed: int = 0

    # ------------------------------------------------------------------ time-step helper
    def stable_dt(self, safety: float = 0.20) -> float:
        """CFL-stable explicit-Euler dt for the diffusion operator.

        For an explicit 5-point Laplacian, stability requires
            dt < 0.5 · min(Δx, Δy)² / (2 · max(D_u, D_v))
        We multiply by ``safety`` to leave headroom for the reaction
        terms, which can also impose a stiffness constraint.
        """
        h = min(self.geometry.dx(), self.geometry.dy())
        D_max = max(self.geometry.D_u, self.geometry.D_v)
        if D_max <= 0.0:
            raise ValueError("Need at least one positive diffusion coefficient.")
        return float(safety * h ** 2 / (4.0 * D_max))

    # ------------------------------------------------------------------ Laplacian
    @staticmethod
    def laplacian5(field: NDArray[np.float64], dx: float, dy: float) -> NDArray[np.float64]:
        """5-point periodic Laplacian on a 2D grid (∂²/∂x² + ∂²/∂y²)."""
        lap = (
            (np.roll(field, +1, axis=1) - 2.0 * field + np.roll(field, -1, axis=1)) / (dx ** 2)
            + (np.roll(field, +1, axis=0) - 2.0 * field + np.roll(field, -1, axis=0)) / (dy ** 2)
        )
        return lap

    # ------------------------------------------------------------------ integrator
    def evolve(
        self,
        u_init: Optional[NDArray[np.float64]] = None,
        v_init: Optional[NDArray[np.float64]] = None,
        record_every: int = 200,
    ) -> Dict[str, NDArray[np.float64]]:
        """Integrate Gierer-Meinhardt forward in time.

        Initial conditions default to a small white-noise perturbation
        about the homogeneous fixed point (u₀ = v₀ = 1).

        Returns
        -------
        Dict with keys
          ``u_final`` (n_y, n_x), ``v_final`` (n_y, n_x),
          ``u_history`` (n_records, n_y, n_x), ``t_history`` (n_records,),
          ``dt`` (float), ``n_steps`` (int).
        """
        geom = self.geometry
        dx, dy = geom.dx(), geom.dy()
        dt = self.dt if self.dt > 0.0 else self.stable_dt()
        D_u, D_v = geom.D_u, geom.D_v
        rho, mu, nu = geom.rho, geom.mu, geom.nu

        if u_init is None:
            u = geom.random_perturbation_field(
                base=1.0, amplitude=0.05, rng_seed=self.rng_seed
            )
        else:
            u = np.array(u_init, dtype=np.float64, copy=True)
        if v_init is None:
            v = geom.random_perturbation_field(
                base=1.0, amplitude=0.05, rng_seed=self.rng_seed + 1
            )
        else:
            v = np.array(v_init, dtype=np.float64, copy=True)

        # Keep v away from zero to avoid divide-by-zero in u²/v.
        v = np.maximum(v, 1.0e-6)

        history_u: List[NDArray[np.float64]] = []
        history_t: List[float] = []
        for step in range(self.n_steps):
            lap_u = self.laplacian5(u, dx, dy)
            lap_v = self.laplacian5(v, dx, dy)
            f = rho * (u ** 2 / np.maximum(v, 1.0e-6) - mu * u)
            g = rho * u ** 2 - nu * v
            u = u + dt * (D_u * lap_u + f)
            v = v + dt * (D_v * lap_v + g)
            v = np.maximum(v, 1.0e-6)
            if step % record_every == 0:
                history_u.append(u.copy())
                history_t.append(step * dt)

        history_u.append(u.copy())
        history_t.append(self.n_steps * dt)
        return {
            "u_final": u,
            "v_final": v,
            "u_history": np.array(history_u),
            "t_history": np.array(history_t),
            "dt": np.array([dt]),
            "n_steps": np.array([self.n_steps]),
        }

    # ------------------------------------------------------------------ pattern wavelength readout
    def measured_wavelength(self, u_field: NDArray[np.float64]) -> float:
        """Estimate the dominant wavelength from the radial Fourier spectrum.

        Computes |FFT(u − ⟨u⟩)|², bins by radial wavenumber, and returns
        2π / q* where q* is the bin centre with the largest power.
        """
        geom = self.geometry
        u_centred = u_field - float(np.mean(u_field))
        spec = np.abs(np.fft.fft2(u_centred)) ** 2
        n_y, n_x = u_field.shape
        kx = 2.0 * np.pi * np.fft.fftfreq(n_x, d=geom.dx())
        ky = 2.0 * np.pi * np.fft.fftfreq(n_y, d=geom.dy())
        KX, KY = np.meshgrid(kx, ky, indexing="xy")
        K = np.sqrt(KX ** 2 + KY ** 2)
        # Radial bin
        k_max = float(np.max(K))
        n_bins = max(8, min(n_x, n_y) // 2)
        bins = np.linspace(0.0, k_max, n_bins + 1)
        which = np.digitize(K.ravel(), bins) - 1
        which = np.clip(which, 0, n_bins - 1)
        power = np.zeros(n_bins, dtype=np.float64)
        counts = np.zeros(n_bins, dtype=np.int_)
        for idx in range(n_bins):
            mask = (which == idx)
            power[idx] = float(np.sum(spec.ravel()[mask]))
            counts[idx] = int(np.sum(mask))
        # Avoid the DC bin (idx 0) and any zero-count bins
        valid = (counts > 0) & (np.arange(n_bins) > 0)
        if not np.any(valid):
            return float("inf")
        idx_star = int(np.argmax(np.where(valid, power, -np.inf)))
        q_star = 0.5 * (bins[idx_star] + bins[idx_star + 1])
        if q_star <= 0.0:
            return float("inf")
        return float(2.0 * np.pi / q_star)

    # ------------------------------------------------------------------ pattern classifier
    def classify_pattern(
        self,
        u_field: NDArray[np.float64],
        anisotropy_threshold: float = 1.5,
        spot_excess_kurtosis: float = 1.0,
    ) -> str:
        """Heuristic stripe / spot / maze classifier from the FFT.

        Spots are isotropic + the field has high (positive) excess
        kurtosis (sharp peaks above a calm baseline). Stripes are
        anisotropic in the FFT (one direction dominates). Mazes are
        isotropic + low kurtosis (band-like, no sharp peaks).
        """
        u_centred = u_field - float(np.mean(u_field))
        spec = np.abs(np.fft.fft2(u_centred)) ** 2
        n_y, n_x = u_field.shape
        # Compare the power along the kx and ky axes (excluding DC).
        kx_power = float(np.sum(spec[0, 1:]))
        ky_power = float(np.sum(spec[1:, 0]))
        ratio = max(kx_power, ky_power) / max(min(kx_power, ky_power), 1.0e-20)
        # Excess kurtosis of the field
        x = u_field.ravel()
        m = float(np.mean(x))
        s = float(np.std(x))
        if s <= 0.0:
            return "uniform"
        kurt = float(np.mean(((x - m) / s) ** 4)) - 3.0
        if ratio > anisotropy_threshold:
            return "stripes"
        if kurt > spot_excess_kurtosis:
            return "spots"
        return "maze"

    # ------------------------------------------------------------------ French-flag readout
    def french_flag_gradient(
        self,
        c_source: float = 1.0,
        c_sink: float = 0.0,
        D: float = 1.0e-3,
        decay: float = 1.0,
        n_x: Optional[int] = None,
    ) -> Dict[str, NDArray[np.float64]]:
        """Steady-state 1D morphogen gradient with diffusion + decay.

        Solves D · ∂²c/∂x² − k · c = 0 with c(0) = c_source,
        c(L) = c_sink. The closed-form solution is exponentially
        decaying:

            c(x) = c_sink + (c_source − c_sink) · sinh(α(L − x)) / sinh(α L)

        where α = √(decay / D) is the inverse "morphogen length scale".

        Returns a dict with ``x``, ``c``, ``alpha`` (1/length).
        """
        if D <= 0.0:
            raise ValueError("Diffusion D must be positive.")
        if decay < 0.0:
            raise ValueError("Decay must be non-negative.")
        L = float(self.geometry.Lx)
        if n_x is None:
            n_x = self.geometry.n_x
        x = np.linspace(0.0, L, n_x)
        if decay == 0.0:
            # Pure diffusion ⇒ linear gradient
            c = c_source + (c_sink - c_source) * (x / L)
            alpha = 0.0
        else:
            alpha = float(math.sqrt(decay / D))
            c = c_sink + (c_source - c_sink) * np.sinh(alpha * (L - x)) / math.sinh(alpha * L)
        return {"x": x, "c": c, "alpha": np.array([alpha])}

    def cell_fates(
        self,
        c: NDArray[np.float64],
        theta_low: float = 0.33,
        theta_high: float = 0.66,
    ) -> NDArray[np.int_]:
        """Wolpert French-flag fate readout (0 / 1 / 2 = blue / white / red)."""
        return french_flag_fates(c, theta_low, theta_high)

    # ------------------------------------------------------------------ summary
    def summary(self) -> Dict[str, float]:
        geom = self.geometry
        dt = self.dt if self.dt > 0.0 else self.stable_dt()
        return {
            "label": geom.label,
            "dt": dt,
            "n_steps": self.n_steps,
            "total_time": dt * self.n_steps,
            "predicted_wavelength": geom.predicted_wavelength(),
            "D_v_over_D_u": geom.diffusion_ratio(),
            "turing_unstable": geom.turing_unstable(),
        }


# ---------------------------------------------------------------------------
# 3. Convenience constructors for the canonical regimes
# ---------------------------------------------------------------------------

def make_morphogenesis(
    name: str,
    Lx: float = 1.0,
    Ly: float = 1.0,
    n_x: int = 64,
    n_y: int = 64,
    n_steps: int = 4000,
    rng_seed: int = 0,
) -> MorphogenesisSimulator:
    """Build a MorphogenesisSimulator for one of the LITERATURE regimes."""
    if name not in LITERATURE:
        raise KeyError(
            f"Unknown morphogenesis regime '{name}'. "
            f"Available: {sorted(LITERATURE.keys())}"
        )
    spec = LITERATURE[name]
    geom = MorphogenesisGeometry(
        Lx=Lx, Ly=Ly, n_x=n_x, n_y=n_y,
        D_u=spec["D_u"], D_v=spec["D_v"],
        rho=spec["rho"], mu=spec["mu"], nu=spec["nu"],
        label=spec["label"],
    )
    return MorphogenesisSimulator(
        geometry=geom, n_steps=n_steps, rng_seed=rng_seed,
    )


def all_canonical_regimes(
    Lx: float = 1.0, Ly: float = 1.0,
    n_x: int = 64, n_y: int = 64,
    n_steps: int = 4000, rng_seed: int = 0,
) -> List[MorphogenesisSimulator]:
    """Return one simulator per LITERATURE regime."""
    return [
        make_morphogenesis(n, Lx=Lx, Ly=Ly, n_x=n_x, n_y=n_y,
                           n_steps=n_steps, rng_seed=rng_seed)
        for n in LITERATURE.keys()
    ]
