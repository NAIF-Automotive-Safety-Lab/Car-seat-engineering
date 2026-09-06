"""Paired chaos / dynamical-systems GEOMETRY + SIMULATOR.

Connects the B3 / stiff-medium ontology to deterministic chaos by treating
trajectories in phase space as substrate strain orbits, attractors as
invariant strain manifolds, and Lyapunov exponents as the local rate at
which two infinitesimally separated strain configurations diverge.

GEOMETRY
--------
``ChaosGeometry`` is the substrate-side description of a dissipative
dynamical system.  It carries:

    * ``dim``                 -- phase-space dimension n
    * ``state_bounds``        -- bounding box used for plotting / sampling
    * ``attractor_type``      -- "fixed_point", "limit_cycle", "torus" or
                                 "strange" (chaotic)
    * ``hausdorff_estimate``  -- estimated correlation/box-counting
                                  dimension D_2 of the attractor
    * ``lyapunov_spectrum``   -- list (lambda_1, lambda_2, ...) sorted in
                                  descending order

The geometry knows that for a strange attractor:
    * sum(lambda_i) < 0   (dissipative -- volumes shrink)
    * lambda_1 > 0        (sensitive dependence on initial conditions)
    * D_KY = j + sum_{i<=j}(lambda_i)/|lambda_{j+1}|  (Kaplan-Yorke
                                                       dimension)

SIMULATOR
---------
``ChaosSimulator`` exposes the canonical demonstration systems:

    * ``simulate_lorenz``     -- Lorenz 1963 with classical chaotic
                                  parameters sigma=10, beta=8/3, rho=28
                                  via a fixed-step RK4 integrator.
    * ``lorenz_lyapunov``     -- largest Lyapunov exponent of Lorenz
                                  estimated by integrating a normalised
                                  tangent vector along the trajectory
                                  (Benettin et al., 1980).  Expected
                                  numerical value ~0.9 (the textbook value
                                  is approximately 0.9056).
    * ``logistic_map``        -- iteration of x_{n+1} = r * x_n * (1 - x_n).
    * ``logistic_bifurcation``-- bifurcation diagram (set of asymptotic
                                  iterates at each r) for the logistic map.
    * ``logistic_lyapunov``   -- closed-form Lyapunov exponent
                                  lambda(r) = (1/N) sum log|r * (1 - 2x_n)|.
    * ``feigenbaum_constants``-- universal constants of period-doubling
                                  cascades:
                                      delta = lim (r_n - r_{n-1})/(r_{n+1} - r_n)
                                            ~ 4.669 201 609
                                      alpha = lim d_n / d_{n+1}
                                            ~ 2.502 907 875

Substrate interpretation:
    * a state x in R^n is the strain configuration of a coarse-grained
      substrate cell tracked over time;
    * the vector field x_dot = F(x) is the local stress-relaxation law;
    * a strange attractor is the fractal set of strain configurations
      that the substrate visits indefinitely;
    * positive Lyapunov exponent lambda_1 measures the rate at which an
      infinitesimal strain perturbation amplifies under the substrate
      dynamics, i.e. the rate of substrate-information loss into
      coarse-grained scales;
    * the Feigenbaum constants are universal because the period-doubling
      cascade is governed by a single universal renormalisation-group
      fixed-point operator, which is metric-free and hence does not
      depend on the microscopic substrate parameters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Substrate / framework constants
# ---------------------------------------------------------------------------

LAMBDA_QCD_MEV: float = 200.0

# Classical Lorenz parameters (Lorenz 1963).  rho > rho_H ~ 24.74 puts the
# system in the chaotic regime.
LORENZ_SIGMA: float = 10.0
LORENZ_BETA: float = 8.0 / 3.0
LORENZ_RHO: float = 28.0

# Universal Feigenbaum constants (high-precision reference values).
FEIGENBAUM_DELTA: float = 4.669201609102990
FEIGENBAUM_ALPHA: float = 2.502907875095892

# Logistic-map period-doubling onset r-values (Strogatz, table 10.6.1).
LOGISTIC_R_PERIOD_DOUBLING: Tuple[float, ...] = (
    3.000000000,    # 1 -> 2
    3.449490000,    # 2 -> 4
    3.544090000,    # 4 -> 8
    3.564407000,    # 8 -> 16
    3.568759000,    # 16 -> 32
    3.569692000,    # 32 -> 64
    3.569891000,    # 64 -> 128
)
LOGISTIC_R_INFINITY: float = 3.569945672  # accumulation point


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class ChaosGeometry:
    """Substrate-side description of a (possibly chaotic) attractor."""

    dim: int
    state_bounds: Tuple[Tuple[float, float], ...]
    attractor_type: str = "strange"
    lyapunov_spectrum: Optional[Sequence[float]] = None
    hausdorff_estimate: Optional[float] = None
    name: str = "anonymous"

    def __post_init__(self) -> None:
        if self.dim <= 0:
            raise ValueError(f"dim must be positive, got {self.dim}")
        if len(self.state_bounds) != self.dim:
            raise ValueError(
                f"state_bounds must have {self.dim} entries, "
                f"got {len(self.state_bounds)}"
            )
        for lo, hi in self.state_bounds:
            if not hi > lo:
                raise ValueError(
                    f"state_bounds entries must satisfy hi > lo, got ({lo}, {hi})"
                )
        if self.attractor_type not in {
            "fixed_point", "limit_cycle", "torus", "strange",
        }:
            raise ValueError(
                f"attractor_type must be one of "
                f"fixed_point/limit_cycle/torus/strange, "
                f"got {self.attractor_type!r}"
            )
        if self.lyapunov_spectrum is not None:
            spec = list(float(x) for x in self.lyapunov_spectrum)
            if len(spec) != self.dim:
                raise ValueError(
                    f"lyapunov_spectrum length {len(spec)} != dim {self.dim}"
                )
            # ensure descending
            if any(spec[i] < spec[i + 1] - 1e-12 for i in range(len(spec) - 1)):
                raise ValueError("lyapunov_spectrum must be sorted descending")
            self.lyapunov_spectrum = spec

    # -- Sampling and validation helpers ------------------------------------

    def in_bounds(self, x: np.ndarray) -> bool:
        x = np.asarray(x, dtype=float).ravel()
        if x.shape != (self.dim,):
            return False
        for xi, (lo, hi) in zip(x, self.state_bounds):
            if xi < lo or xi > hi:
                return False
        return True

    # -- Diagnostics --------------------------------------------------------

    def is_chaotic(self) -> bool:
        """Operational definition: attractor_type == 'strange' and
        the leading Lyapunov exponent (if known) is strictly positive."""
        if self.attractor_type != "strange":
            return False
        if self.lyapunov_spectrum is None:
            return True  # take attractor_type as the source of truth
        return self.lyapunov_spectrum[0] > 0.0

    def is_dissipative(self) -> bool:
        """sum(lambda_i) < 0 means phase-space volume contracts."""
        if self.lyapunov_spectrum is None:
            return False
        return sum(self.lyapunov_spectrum) < 0.0

    def kaplan_yorke_dimension(self) -> Optional[float]:
        """Kaplan-Yorke (Lyapunov) dimension:

            D_KY = j + (sum_{i<=j} lambda_i) / |lambda_{j+1}|

        where j is the largest index with sum_{i<=j} lambda_i >= 0.
        Returns None when the spectrum is missing or the formula does not
        apply (all-negative spectrum).
        """
        if self.lyapunov_spectrum is None:
            return None
        spec = list(self.lyapunov_spectrum)
        running = 0.0
        j = -1
        for i, lam in enumerate(spec):
            if running + lam < 0.0:
                break
            running += lam
            j = i
        if j < 0:
            # No positive partial sum -> stable fixed point
            return 0.0
        if j == len(spec) - 1:
            # Whole spectrum is non-negative; D_KY = n by convention
            return float(self.dim)
        denom = abs(spec[j + 1])
        if denom < 1e-15:
            return float(j + 1)
        return float(j + 1) + running / denom


# Pre-built geometries for the canonical demos
def lorenz_geometry() -> ChaosGeometry:
    """Classical Lorenz attractor at sigma=10, beta=8/3, rho=28.

    Numerical Lyapunov spectrum is approximately
        lambda_1 = +0.9056
        lambda_2 =  0
        lambda_3 = -14.572
    (Sprott 2003, table 4.3).  D_KY = 2 + 0.9056/14.572 ~ 2.062.
    """
    return ChaosGeometry(
        dim=3,
        state_bounds=((-30.0, 30.0), (-30.0, 30.0), (0.0, 60.0)),
        attractor_type="strange",
        lyapunov_spectrum=(0.9056, 0.0, -14.572),
        hausdorff_estimate=2.062,
        name="Lorenz",
    )


def logistic_geometry(r: float = 3.9) -> ChaosGeometry:
    """1D logistic map.  For r = 3.9 the leading exponent is ~0.494 (Strogatz)."""
    return ChaosGeometry(
        dim=1,
        state_bounds=((0.0, 1.0),),
        attractor_type="strange" if r > LOGISTIC_R_INFINITY else "limit_cycle",
        lyapunov_spectrum=None,
        hausdorff_estimate=None,
        name=f"Logistic(r={r:.4f})",
    )


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class ChaosSimulator:
    """Numerical integrators for the canonical chaotic systems."""

    sigma: float = LORENZ_SIGMA
    beta: float = LORENZ_BETA
    rho: float = LORENZ_RHO

    # ---- Lorenz vector field and its Jacobian -----------------------------

    def lorenz_rhs(self, state: np.ndarray) -> np.ndarray:
        """Lorenz 1963 vector field F(x, y, z)."""
        x, y, z = state
        return np.array([
            self.sigma * (y - x),
            x * (self.rho - z) - y,
            x * y - self.beta * z,
        ])

    def lorenz_jacobian(self, state: np.ndarray) -> np.ndarray:
        x, y, z = state
        return np.array([
            [-self.sigma, self.sigma, 0.0],
            [self.rho - z, -1.0, -x],
            [y, x, -self.beta],
        ])

    # ---- Generic RK4 integrator -------------------------------------------

    def rk4_step(self, rhs: Callable[[np.ndarray], np.ndarray],
                 state: np.ndarray, dt: float) -> np.ndarray:
        k1 = rhs(state)
        k2 = rhs(state + 0.5 * dt * k1)
        k3 = rhs(state + 0.5 * dt * k2)
        k4 = rhs(state + dt * k3)
        return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # ---- Lorenz attractor trajectory --------------------------------------

    def simulate_lorenz(
        self,
        x0: Tuple[float, float, float] = (1.0, 1.0, 1.0),
        t_max: float = 40.0,
        dt: float = 0.01,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Integrate the Lorenz system with fixed-step RK4.

        Returns (t, X) where X has shape (n_steps, 3).
        """
        if t_max <= 0:
            raise ValueError("t_max must be positive")
        if dt <= 0:
            raise ValueError("dt must be positive")
        n_steps = int(np.ceil(t_max / dt)) + 1
        t = np.arange(n_steps, dtype=float) * dt
        X = np.empty((n_steps, 3), dtype=float)
        X[0] = np.asarray(x0, dtype=float)
        for k in range(n_steps - 1):
            X[k + 1] = self.rk4_step(self.lorenz_rhs, X[k], dt)
        return t, X

    # ---- Lorenz Lyapunov exponent (Benettin et al., 1980) ----------------

    def lorenz_lyapunov(
        self,
        x0: Tuple[float, float, float] = (1.0, 1.0, 1.0),
        t_max: float = 100.0,
        dt: float = 0.01,
        transient: float = 5.0,
    ) -> float:
        """Largest Lyapunov exponent of the Lorenz attractor.

        The textbook value is 0.9056 (Sprott 2003).  With t_max=100 the
        Benettin estimator typically lands within a few percent of that.
        """
        if t_max <= transient:
            raise ValueError("t_max must exceed transient")
        n_steps = int(np.ceil(t_max / dt))
        n_trans = int(np.ceil(transient / dt))
        x = np.asarray(x0, dtype=float).copy()
        # tangent vector starts as a unit perturbation
        d = np.array([1.0, 0.0, 0.0], dtype=float)
        log_sum = 0.0
        n_count = 0
        for k in range(n_steps):
            # propagate the trajectory
            x = self.rk4_step(self.lorenz_rhs, x, dt)
            # propagate the tangent linearly using the Jacobian (Euler is
            # sufficient and stable when combined with periodic
            # renormalisation)
            J = self.lorenz_jacobian(x)
            d = d + dt * J @ d
            norm = float(np.linalg.norm(d))
            if not math.isfinite(norm) or norm == 0.0:
                # numerical pathology -- bail out
                break
            if k >= n_trans:
                log_sum += math.log(norm)
                n_count += 1
            d = d / norm
        if n_count == 0:
            raise RuntimeError("no post-transient steps recorded")
        return log_sum / (n_count * dt)

    # ---- Logistic map -----------------------------------------------------

    def logistic_map(self, x: float, r: float) -> float:
        """One step of x_{n+1} = r * x_n * (1 - x_n)."""
        if not 0.0 <= x <= 1.0:
            raise ValueError(f"x must lie in [0, 1], got {x}")
        if r < 0.0:
            raise ValueError(f"r must be non-negative, got {r}")
        return r * x * (1.0 - x)

    def iterate_logistic(
        self,
        x0: float,
        r: float,
        n_iter: int,
        n_skip: int = 0,
    ) -> np.ndarray:
        """Return the orbit (x_{n_skip}, ..., x_{n_skip + n_iter - 1})."""
        if n_iter <= 0:
            raise ValueError("n_iter must be positive")
        x = float(x0)
        for _ in range(n_skip):
            x = r * x * (1.0 - x)
        out = np.empty(n_iter, dtype=float)
        for k in range(n_iter):
            x = r * x * (1.0 - x)
            out[k] = x
        return out

    def logistic_bifurcation(
        self,
        r_values: np.ndarray,
        x0: float = 0.5,
        n_skip: int = 1000,
        n_keep: int = 200,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Sample the asymptotic iterates of the logistic map for many r.

        Returns two flat arrays (R, X) of equal length suitable for a
        scatter-plot of the bifurcation diagram.
        """
        r_values = np.asarray(r_values, dtype=float).ravel()
        n_r = r_values.size
        R_out = np.empty(n_r * n_keep, dtype=float)
        X_out = np.empty(n_r * n_keep, dtype=float)
        idx = 0
        for r in r_values:
            x = float(x0)
            for _ in range(n_skip):
                x = r * x * (1.0 - x)
            for _ in range(n_keep):
                x = r * x * (1.0 - x)
                R_out[idx] = r
                X_out[idx] = x
                idx += 1
        return R_out, X_out

    def logistic_lyapunov(
        self,
        r: float,
        x0: float = 0.5,
        n_iter: int = 5000,
        n_skip: int = 1000,
    ) -> float:
        """Closed-form Lyapunov exponent for the logistic map.

            lambda(r) = (1/N) sum_{n=1..N} log|r * (1 - 2 x_n)|

        Positive iff the orbit is chaotic.
        """
        if not 0.0 <= x0 <= 1.0:
            raise ValueError("x0 must lie in [0, 1]")
        if n_iter <= 0:
            raise ValueError("n_iter must be positive")
        x = float(x0)
        for _ in range(n_skip):
            x = r * x * (1.0 - x)
        log_sum = 0.0
        n_used = 0
        for _ in range(n_iter):
            d = abs(r * (1.0 - 2.0 * x))
            if d > 0.0:
                log_sum += math.log(d)
                n_used += 1
            x = r * x * (1.0 - x)
        if n_used == 0:
            return float("-inf")
        return log_sum / n_used

    def logistic_lyapunov_curve(
        self,
        r_values: np.ndarray,
        x0: float = 0.5,
        n_iter: int = 2000,
        n_skip: int = 500,
    ) -> np.ndarray:
        r_values = np.asarray(r_values, dtype=float).ravel()
        out = np.empty_like(r_values)
        for i, r in enumerate(r_values):
            out[i] = self.logistic_lyapunov(r, x0=x0,
                                            n_iter=n_iter, n_skip=n_skip)
        return out

    # ---- Feigenbaum constants --------------------------------------------

    def feigenbaum_constants(self) -> Tuple[float, float]:
        """Return the universal Feigenbaum constants (delta, alpha).

        These are reference values; the simulator does not re-derive them
        from scratch (period-doubling cascade detection is provided
        separately, see ``feigenbaum_delta_estimate``).
        """
        return FEIGENBAUM_DELTA, FEIGENBAUM_ALPHA

    def feigenbaum_delta_estimate(
        self,
        r_doublings: Sequence[float] = LOGISTIC_R_PERIOD_DOUBLING,
    ) -> float:
        """Estimate Feigenbaum delta from period-doubling onset r-values.

            delta_n = (r_{n+1} - r_n) / (r_{n+2} - r_{n+1})

        We average the last few estimates (which are closest to the
        asymptotic value 4.6692...).
        """
        rs = list(r_doublings)
        if len(rs) < 4:
            raise ValueError(
                "need at least 4 period-doubling onsets to estimate delta")
        gaps = [rs[i + 1] - rs[i] for i in range(len(rs) - 1)]
        ratios = [gaps[i] / gaps[i + 1] for i in range(len(gaps) - 1)]
        # last ratio is most converged
        return float(ratios[-1])


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------

def detect_period_doubling_cascade(
    r_doublings: Sequence[float] = LOGISTIC_R_PERIOD_DOUBLING,
) -> bool:
    """True iff successive gaps r_{n+1} - r_n are strictly decreasing."""
    rs = list(r_doublings)
    if len(rs) < 3:
        return False
    gaps = [rs[i + 1] - rs[i] for i in range(len(rs) - 1)]
    return all(gaps[i] > gaps[i + 1] > 0.0 for i in range(len(gaps) - 1))


__all__ = [
    "LAMBDA_QCD_MEV",
    "LORENZ_SIGMA",
    "LORENZ_BETA",
    "LORENZ_RHO",
    "FEIGENBAUM_DELTA",
    "FEIGENBAUM_ALPHA",
    "LOGISTIC_R_PERIOD_DOUBLING",
    "LOGISTIC_R_INFINITY",
    "ChaosGeometry",
    "ChaosSimulator",
    "lorenz_geometry",
    "logistic_geometry",
    "detect_period_doubling_cascade",
]
