"""Cone-bouncing oscillator — the §18.35 mass mechanism.

A vector traveling at speed c on the 45° cone has a "preferred direction"
it would like to travel in. The medium's directional stiffness κ pulls it
back when it strays. The vector oscillates around its preferred direction.

The momentum of this oscillation IS the rest mass:
    m c² = ℏ × ω_bounce
    ω_bounce = √(κ / I)

Particles with no preferred direction (photons): κ = 0 → m = 0
Light neutrinos: small κ → small mass
Charged leptons: moderate κ
Heavy carriers (kinks): large κ

3D extension (§18.35 full 3D wobble):
    The preferred direction in 3D defines two perpendicular wobble planes.
    Each plane carries its own stiffness (κ_x, κ_y) and inertia (I_x, I_y),
    giving two independent harmonic modes with frequencies ω_x and ω_y.
    The combined rest mass energy is the quadrature sum:
        (m c²)² = (ℏ ω_x)² + (ℏ ω_y)²
    Per §18.53.7 the minimum wobble of each mode enforces Δx·Δp ≥ ℏ/2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import numpy as np


# ---------------------------------------------------------------------------
# Existing 1-D wobble class — kept for full backward compatibility
# ---------------------------------------------------------------------------

@dataclass
class ConeBouncingParticle:
    """A particle described as a wobbling vector on the 45° cone.

    Attributes:
        preferred_direction: 3D unit vector defining the particle's "intended" motion
        kappa: directional stiffness (medium's restoring strength)
        I: effective inertia (vector mass × geometric factors)
        theta: current wobble angle (small = aligned with preferred direction)
        theta_dot: angular velocity of wobble
    """
    preferred_direction: np.ndarray
    kappa: float
    I: float = 1.0
    theta: float = 0.1
    theta_dot: float = 0.0

    def __post_init__(self):
        # Normalize preferred direction
        norm = np.linalg.norm(self.preferred_direction)
        if norm > 0:
            self.preferred_direction = self.preferred_direction / norm

    @property
    def omega_bounce(self) -> float:
        """Bouncing frequency ω = √(κ/I)."""
        if self.kappa <= 0:
            return 0.0
        return np.sqrt(self.kappa / self.I)

    @property
    def rest_mass(self) -> float:
        """Rest mass m c² = ℏ × ω (with ℏ = 1 in natural units)."""
        return self.omega_bounce  # natural units

    def step(self, dt: float) -> None:
        """Time-step the wobble dynamics: θ̈ = -(κ/I) × θ."""
        # Velocity-Verlet for harmonic oscillator
        if self.kappa <= 0:
            # Massless: no oscillation, just propagate freely
            return

        omega_sq = self.kappa / self.I
        theta_ddot = -omega_sq * self.theta
        theta_dot_half = self.theta_dot + 0.5 * theta_ddot * dt
        new_theta = self.theta + theta_dot_half * dt
        new_theta_ddot = -omega_sq * new_theta
        self.theta_dot = theta_dot_half + 0.5 * new_theta_ddot * dt
        self.theta = new_theta

    def instantaneous_velocity(self, c: float = 1.0) -> np.ndarray:
        """Current velocity vector (tilted by θ from preferred direction).

        Returns 3D velocity at speed c, tilted by self.theta from the
        preferred direction in some perpendicular plane.
        """
        # Build orthonormal frame around preferred direction
        z_hat = self.preferred_direction
        # Choose any vector not parallel to z_hat
        if abs(z_hat[2]) < 0.9:
            ref = np.array([0.0, 0.0, 1.0])
        else:
            ref = np.array([1.0, 0.0, 0.0])
        x_hat = np.cross(z_hat, ref)
        x_hat = x_hat / np.linalg.norm(x_hat)

        # Velocity tilted by θ
        return c * (np.cos(self.theta) * z_hat + np.sin(self.theta) * x_hat)

    def average_translational_speed(self) -> float:
        """⟨v_z⟩ = c × ⟨cos θ⟩ averaged over wobble period."""
        # For small wobble: ⟨cos θ⟩ ≈ 1 - ⟨θ²⟩/2 = 1 - amplitude²/4
        amplitude_sq = self.theta**2 + (self.theta_dot / self.omega_bounce)**2 if self.omega_bounce > 0 else 0
        return 1.0 - amplitude_sq / 4


# ---------------------------------------------------------------------------
# New 3-D wobble class — §18.35 with two independent angular DoF
# ---------------------------------------------------------------------------

def _build_perp_frame(z_hat: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return two unit vectors (x_hat, y_hat) perpendicular to z_hat.

    The frame (x_hat, y_hat, z_hat) is right-handed.

    Args:
        z_hat: Unit vector defining the preferred direction.

    Returns:
        Tuple of (x_hat, y_hat) both orthogonal to z_hat and to each other.
    """
    if abs(z_hat[2]) < 0.9:
        ref = np.array([0.0, 0.0, 1.0])
    else:
        ref = np.array([1.0, 0.0, 0.0])
    x_hat = np.cross(z_hat, ref)
    x_hat = x_hat / np.linalg.norm(x_hat)
    y_hat = np.cross(z_hat, x_hat)
    y_hat = y_hat / np.linalg.norm(y_hat)
    return x_hat, y_hat


@dataclass
class ConeBouncingParticle3D:
    """A particle with full 3-D wobble dynamics around its preferred axis.

    Two independent harmonic modes capture wobble in the x̂ and ŷ planes
    perpendicular to the preferred direction (ẑ).  Per §18.35 the rest-mass
    energy follows from the quadrature combination of both modes:

        (m c²)² = (ℏ ω_x)² + (ℏ ω_y)²

    so that for the isotropic case (κ_x = κ_y, I_x = I_y) this simplifies to

        m c² = ℏ ω √2.

    Attributes:
        preferred_direction: 3-D unit vector — the particle's preferred axis.
        kappa_x: Directional stiffness in the x̂ wobble plane.
        kappa_y: Directional stiffness in the ŷ wobble plane.
        I_x: Effective inertia for the x̂ mode.
        I_y: Effective inertia for the ŷ mode.
        theta_x: Current wobble angle away from ẑ toward x̂.
        theta_y: Current wobble angle away from ẑ toward ŷ.
        theta_x_dot: Angular velocity of the θ_x mode.
        theta_y_dot: Angular velocity of the θ_y mode.
        hbar: Reduced Planck constant (natural units: 1.0).
    """

    preferred_direction: np.ndarray
    kappa_x: float
    kappa_y: float
    I_x: float = 1.0
    I_y: float = 1.0
    theta_x: float = 0.1
    theta_y: float = 0.0
    theta_x_dot: float = 0.0
    theta_y_dot: float = 0.0
    hbar: float = 1.0

    # Derived frame vectors — populated in __post_init__
    _x_hat: np.ndarray = field(init=False, repr=False)
    _y_hat: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Normalize preferred direction and build the local perpendicular frame."""
        norm = np.linalg.norm(self.preferred_direction)
        if norm > 0:
            self.preferred_direction = self.preferred_direction / norm
        self._x_hat, self._y_hat = _build_perp_frame(self.preferred_direction)

    # ------------------------------------------------------------------
    # Frequency properties
    # ------------------------------------------------------------------

    @property
    def omega_x(self) -> float:
        """Bouncing frequency for the θ_x mode: ω_x = √(κ_x / I_x).

        Returns:
            Angular frequency in rad / time_unit.  Zero when κ_x ≤ 0.
        """
        if self.kappa_x <= 0:
            return 0.0
        return float(np.sqrt(self.kappa_x / self.I_x))

    @property
    def omega_y(self) -> float:
        """Bouncing frequency for the θ_y mode: ω_y = √(κ_y / I_y).

        Returns:
            Angular frequency in rad / time_unit.  Zero when κ_y ≤ 0.
        """
        if self.kappa_y <= 0:
            return 0.0
        return float(np.sqrt(self.kappa_y / self.I_y))

    # ------------------------------------------------------------------
    # Mass
    # ------------------------------------------------------------------

    @property
    def total_rest_mass(self) -> float:
        """Combined rest-mass energy from both wobble modes (natural units ℏ = 1).

        Per §18.35 the two independent angular modes each contribute
        ℏ ω to the rest-mass energy in quadrature:

            m c² = √( (ℏ ω_x)² + (ℏ ω_y)² )

        Returns:
            Rest-mass energy m c² in natural units where ℏ = self.hbar.
        """
        return float(np.sqrt((self.hbar * self.omega_x) ** 2 + (self.hbar * self.omega_y) ** 2))

    # ------------------------------------------------------------------
    # Dynamics
    # ------------------------------------------------------------------

    def step(self, dt: float) -> None:
        """Advance both wobble angles by one Velocity-Verlet step of size dt.

        Each mode evolves independently as a simple harmonic oscillator:

            θ̈_i = -(κ_i / I_i) θ_i

        Args:
            dt: Time step size.
        """
        self._verlet_step_x(dt)
        self._verlet_step_y(dt)

    def _verlet_step_x(self, dt: float) -> None:
        """Velocity-Verlet integration for the θ_x harmonic mode."""
        if self.kappa_x <= 0:
            return
        omega_sq = self.kappa_x / self.I_x
        ddot = -omega_sq * self.theta_x
        dot_half = self.theta_x_dot + 0.5 * ddot * dt
        new_theta = self.theta_x + dot_half * dt
        new_ddot = -omega_sq * new_theta
        self.theta_x_dot = dot_half + 0.5 * new_ddot * dt
        self.theta_x = new_theta

    def _verlet_step_y(self, dt: float) -> None:
        """Velocity-Verlet integration for the θ_y harmonic mode."""
        if self.kappa_y <= 0:
            return
        omega_sq = self.kappa_y / self.I_y
        ddot = -omega_sq * self.theta_y
        dot_half = self.theta_y_dot + 0.5 * ddot * dt
        new_theta = self.theta_y + dot_half * dt
        new_ddot = -omega_sq * new_theta
        self.theta_y_dot = dot_half + 0.5 * new_ddot * dt
        self.theta_y = new_theta

    # ------------------------------------------------------------------
    # Kinematics
    # ------------------------------------------------------------------

    def instantaneous_velocity_3d(self, c: float = 1.0) -> np.ndarray:
        """Current 3-D velocity vector incorporating both wobble angles.

        The velocity is expressed as a small-angle tilt of the preferred
        direction ẑ simultaneously toward x̂ (by θ_x) and ŷ (by θ_y):

            v = c × (sin θ_x x̂ + sin θ_y ŷ + cos θ_total ẑ) / |...|

        For small angles this approximates to
            v ≈ c × (θ_x x̂ + θ_y ŷ + (1 - (θ_x² + θ_y²)/2) ẑ)

        and is renormalised to maintain |v| = c.

        Args:
            c: Speed (default 1 in natural units).

        Returns:
            3-D numpy array with |v| = c.
        """
        # Raw tilt vector (not yet unit length)
        v_raw = (
            np.sin(self.theta_x) * self._x_hat
            + np.sin(self.theta_y) * self._y_hat
            + np.cos(self.wobble_radius()) * self.preferred_direction
        )
        norm = np.linalg.norm(v_raw)
        if norm < 1e-15:
            return c * self.preferred_direction.copy()
        return c * v_raw / norm

    def wobble_radius(self) -> float:
        """Total angular displacement from the preferred axis.

        For small angles: r ≈ √(θ_x² + θ_y²).

        Returns:
            Angular displacement in radians.
        """
        return float(np.sqrt(self.theta_x ** 2 + self.theta_y ** 2))

    # ------------------------------------------------------------------
    # Amplitude helpers (useful for uncertainty calculations)
    # ------------------------------------------------------------------

    def _amplitude_x(self) -> float:
        """Oscillation amplitude of the θ_x mode (energy-based).

        A_x² = θ_x² + (θ̇_x / ω_x)²
        """
        if self.omega_x <= 0:
            return abs(self.theta_x)
        return float(np.sqrt(self.theta_x ** 2 + (self.theta_x_dot / self.omega_x) ** 2))

    def _amplitude_y(self) -> float:
        """Oscillation amplitude of the θ_y mode (energy-based).

        A_y² = θ_y² + (θ̇_y / ω_y)²
        """
        if self.omega_y <= 0:
            return abs(self.theta_y)
        return float(np.sqrt(self.theta_y ** 2 + (self.theta_y_dot / self.omega_y) ** 2))


# ---------------------------------------------------------------------------
# §18.53.7 — Heisenberg uncertainty verification
# ---------------------------------------------------------------------------

def verify_uncertainty_principle(
    particle: ConeBouncingParticle3D,
    length_scale: float = 1.0,
) -> Dict[str, float]:
    """Verify Δx · Δp ≥ ℏ/2 for each wobble mode.

    Per §18.53.7 the cone-bouncing geometry imposes a minimum wobble whose
    amplitude-momentum product satisfies the Heisenberg bound structurally.

    The position uncertainty for mode i is estimated as the oscillation
    amplitude scaled by a characteristic length L (Compton wavelength proxy):

        Δx_i = A_i × L

    The momentum uncertainty is the oscillation momentum of that mode:

        Δp_i = I_i × ω_i × A_i

    Their product:

        Δx_i · Δp_i = I_i × ω_i × A_i² × L

    For a particle sitting exactly at its zero-point (A² = ℏ / (2 I ω)):

        Δx · Δp = ℏ/2  (at the minimum)

    Args:
        particle: A ``ConeBouncingParticle3D`` instance.
        length_scale: Characteristic length scale L (default 1 in natural
            units; set to Compton wavelength ℏ / (m c) for physical units).

    Returns:
        Dictionary with keys:
            ``delta_x_x``, ``delta_p_x``, ``product_x``: mode-x values.
            ``delta_x_y``, ``delta_p_y``, ``product_y``: mode-y values.
            ``hbar_half``: ℏ/2 for comparison.
            ``satisfies_x``, ``satisfies_y``: bool whether ≥ ℏ/2.
            ``satisfies_both``: combined flag.
    """
    hbar = particle.hbar
    hbar_half = hbar / 2.0

    A_x = particle._amplitude_x()
    A_y = particle._amplitude_y()

    # Position uncertainty: amplitude × length scale
    delta_x_x = A_x * length_scale
    delta_x_y = A_y * length_scale

    # Momentum uncertainty: I × ω × A  (mv for the oscillator)
    delta_p_x = particle.I_x * particle.omega_x * A_x if particle.omega_x > 0 else 0.0
    delta_p_y = particle.I_y * particle.omega_y * A_y if particle.omega_y > 0 else 0.0

    product_x = delta_x_x * delta_p_x
    product_y = delta_x_y * delta_p_y

    satisfies_x = bool(product_x >= hbar_half - 1e-12)
    satisfies_y = bool(product_y >= hbar_half - 1e-12)

    return {
        "delta_x_x": delta_x_x,
        "delta_p_x": delta_p_x,
        "product_x": product_x,
        "delta_x_y": delta_x_y,
        "delta_p_y": delta_p_y,
        "product_y": product_y,
        "hbar_half": hbar_half,
        "satisfies_x": satisfies_x,
        "satisfies_y": satisfies_y,
        "satisfies_both": bool(satisfies_x and satisfies_y),
    }


# ---------------------------------------------------------------------------
# Ensemble simulation
# ---------------------------------------------------------------------------

def simulate_ensemble(
    N_particles: int,
    kappa_distribution: np.ndarray | List[float],
    dt: float,
    n_steps: int,
    I_x: float = 1.0,
    I_y: float = 1.0,
    initial_theta: float = 0.1,
    hbar: float = 1.0,
    rng: np.random.Generator | None = None,
) -> Dict[str, float | np.ndarray]:
    """Create an ensemble of 3-D particles and evolve them, returning statistics.

    Each particle is drawn with κ_x = κ_y sampled from ``kappa_distribution``.
    If ``kappa_distribution`` has fewer than ``N_particles`` entries it is
    sampled with replacement; if it has exactly ``N_particles`` entries it is
    used directly.

    Initial conditions: θ_x = initial_theta, θ_y = 0, both velocities = 0,
    preferred direction along +z̃ axis.

    Args:
        N_particles: Number of particles in the ensemble.
        kappa_distribution: 1-D array of κ values to draw from.
        dt: Integration time step.
        n_steps: Number of Verlet steps to run.
        I_x: Inertia for the x̂ mode (same for all particles).
        I_y: Inertia for the ŷ mode (same for all particles).
        initial_theta: Initial θ_x displacement (same for all particles).
        hbar: Reduced Planck constant.
        rng: Optional numpy random generator for reproducibility.

    Returns:
        Dictionary with keys:
            ``kappas``: array of sampled κ values.
            ``omega_x``: array of ω_x per particle.
            ``omega_y``: array of ω_y per particle.
            ``masses``: array of total rest-mass energies.
            ``wobble_amplitudes``: array of final wobble radii.
            ``mean_mass``: float mean of masses.
            ``std_mass``: float std-dev of masses.
            ``mean_wobble``: float mean of wobble radii.
            ``std_wobble``: float std-dev of wobble radii.
            ``uncertainty_products_x``: Δx·Δp for each particle (x-mode).
            ``uncertainty_products_y``: Δx·Δp for each particle (y-mode).
            ``all_satisfy_uncertainty``: bool — every particle satisfies ≥ ℏ/2.
    """
    if rng is None:
        rng = np.random.default_rng()

    kappa_distribution = np.asarray(kappa_distribution, dtype=float)
    indices = rng.choice(len(kappa_distribution), size=N_particles, replace=True)
    kappas = kappa_distribution[indices]

    preferred = np.array([0.0, 0.0, 1.0])

    particles: List[ConeBouncingParticle3D] = [
        ConeBouncingParticle3D(
            preferred_direction=preferred.copy(),
            kappa_x=float(k),
            kappa_y=float(k),
            I_x=I_x,
            I_y=I_y,
            theta_x=initial_theta,
            theta_y=0.0,
            theta_x_dot=0.0,
            theta_y_dot=0.0,
            hbar=hbar,
        )
        for k in kappas
    ]

    # Evolve
    for _ in range(n_steps):
        for p in particles:
            p.step(dt)

    # Collect final statistics
    omega_x_vals = np.array([p.omega_x for p in particles])
    omega_y_vals = np.array([p.omega_y for p in particles])
    masses = np.array([p.total_rest_mass for p in particles])
    wobble_amps = np.array([p.wobble_radius() for p in particles])

    uncertainty_x = np.array([verify_uncertainty_principle(p)["product_x"] for p in particles])
    uncertainty_y = np.array([verify_uncertainty_principle(p)["product_y"] for p in particles])
    hbar_half = hbar / 2.0
    all_satisfy = bool(np.all(uncertainty_x >= hbar_half - 1e-12) and np.all(uncertainty_y >= hbar_half - 1e-12))

    return {
        "kappas": kappas,
        "omega_x": omega_x_vals,
        "omega_y": omega_y_vals,
        "masses": masses,
        "wobble_amplitudes": wobble_amps,
        "mean_mass": float(np.mean(masses)),
        "std_mass": float(np.std(masses)),
        "mean_wobble": float(np.mean(wobble_amps)),
        "std_wobble": float(np.std(wobble_amps)),
        "uncertainty_products_x": uncertainty_x,
        "uncertainty_products_y": uncertainty_y,
        "all_satisfy_uncertainty": all_satisfy,
    }


# ---------------------------------------------------------------------------
# Unchanged helper functions — backward compatible
# ---------------------------------------------------------------------------

def categorize_particle_by_kappa(kappa: float) -> str:
    """Classify particle type by directional stiffness κ.

    Per §18.35:
    - κ = 0: photon (no preferred direction)
    - κ very small: light neutrino
    - κ moderate: charged lepton
    - κ large: heavy carrier (kink, W/Z analog)
    """
    if kappa < 1e-30:
        return "photon (m=0)"
    elif kappa < 1e-20:
        return "light neutrino"
    elif kappa < 1e-10:
        return "charged lepton"
    elif kappa < 1.0:
        return "intermediate boson"
    else:
        return "heavy kink / GUT-scale"


def mass_from_kappa(kappa: float, I: float = 1.0, hbar: float = 1.0) -> float:
    """Compute m c² from directional stiffness κ.

    m c² = ℏ × √(κ/I)
    """
    if kappa <= 0:
        return 0.0
    return hbar * np.sqrt(kappa / I)


def kappa_for_mass(m_c2: float, I: float = 1.0, hbar: float = 1.0) -> float:
    """Inverse: required κ to produce given mass-energy."""
    return I * (m_c2 / hbar)**2
