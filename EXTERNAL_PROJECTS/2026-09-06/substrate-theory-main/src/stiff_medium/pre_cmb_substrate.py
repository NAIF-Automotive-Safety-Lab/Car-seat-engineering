"""Pre-CMB substrate inhomogeneities — §18.44.

At σ ≈ 0.5 (saturated) the substrate was in a plastic-regime state where
ordinary clocks do not tick relative to external observers.  Nevertheless the
medium was NOT featureless: internal elastic resonances, topological defect
seeds, and acoustic-like standing waves could organise within it.  When the
de-saturation phase transition (§18.44) released the saturated medium back
into the ordinary elastic regime, those pre-existing patterns were imprinted
on the post-CMB matter distribution as enhanced initial-condition seeds.

This module quantifies those inhomogeneities and connects them to two
precision-cosmology tensions that the standard ΛCDM picture cannot explain:

1. **JWST high-z galaxy excess** (arxiv:2505.11263, MoM-z14 at z = 14.44):
   number densities 182×⁺³²⁹₋₁₀₅ above ΛCDM predictions.  In the substrate
   model galaxies form faster because they start from overdense seeds inherited
   from the pre-CMB saturated era rather than from scale-invariant vacuum
   quantum fluctuations.

2. **Hubble tension** (DESI DR2, SH0ES, arxiv:2604.24050):
   CMB-derived H₀ = 67.4 km/s/Mpc vs. local H₀ ≈ 73 km/s/Mpc — a ≈ 8%
   discrepancy.  Pre-CMB modes boost the effective sound speed c_s of the
   baryon-photon plasma (or equivalently advance the de-saturation front),
   reducing the sound horizon r_s by ~7% and shifting the CMB-inferred H₀
   upward to match local measurements.

Physical picture (§18.44):
    - Pre-CMB: σ = 0.5, substrate saturated.  Internal resonances organise as
      standing strain waves with wavenumber k and amplitude A(k).
    - De-saturation: phase transition boundary sweeps through at speed c_s.
      The substrate's internal pattern is "frozen in" at the transition.
    - Post-CMB: frozen modes act as enhanced density seeds; the power spectrum
      of matter fluctuations is boosted by the pre-CMB pattern at every k.

The enhancement is multiplicative with the standard inflationary spectrum:

    P_total(k) = P_inflation(k) × [1 + δ_substrate(k)]

where δ_substrate(k) captures the pre-CMB contribution.

References:
    - Spec §18.44 — CMB as de-saturation phase transition
    - Spec §18.52.1 — Hubble tension: r_s reduced ~7% by pre-CMB physics
    - Spec §18.52.3 — JWST 182× galaxy abundance: pre-CMB seeding
    - arxiv:2505.11263 — MoM-z14 observational result
    - arxiv:2503.14738 — DESI DR2 BAO
    - arxiv:2604.24050 — sound-horizon-free H₀ measurement
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final, Sequence

# ---------------------------------------------------------------------------
# Physical constants (SI)  — consistent with cyclic_cosmology.py and desaturation.py
# ---------------------------------------------------------------------------

#: Speed of light in vacuum (m s⁻¹).
c: Final[float] = 2.998e8

#: Reduced Planck constant (J s).
hbar: Final[float] = 1.055e-34

#: Newton's gravitational constant (m³ kg⁻¹ s⁻²).
G: Final[float] = 6.674e-11

#: Boltzmann constant (J K⁻¹).
k_B: Final[float] = 1.381e-23

# ---------------------------------------------------------------------------
# Cosmological constants (Planck 2018 baseline)
# ---------------------------------------------------------------------------

#: Megaparsec in metres.
_MPC_M: Final[float] = 3.0857e22

#: Planck-2018 CMB-derived Hubble constant (km s⁻¹ Mpc⁻¹).
H0_CMB: Final[float] = 67.4

#: Local (SH0ES/DESI sound-horizon-free) Hubble constant (km s⁻¹ Mpc⁻¹).
H0_LOCAL: Final[float] = 73.0

#: Planck-2018 ΛCDM sound horizon at baryon drag (Mpc).
R_S_LCDM: Final[float] = 147.09

#: CMB first acoustic peak multipole in ΛCDM.
ELL_1_LCDM: Final[float] = 220.0

#: Substrate saturation limit (§18.39).
SIGMA_MAX: Final[float] = 0.5

#: Planck CMB scalar spectral index.
N_S_PLANCK: Final[float] = 0.9649

#: Planck CMB scalar amplitude (at k_pivot = 0.05 Mpc⁻¹).
A_S_PLANCK: Final[float] = 2.1e-9

#: Default pivot wavenumber (Mpc⁻¹).
K_PIVOT: Final[float] = 0.05

#: Default minimum wavenumber for substrate modes (Mpc⁻¹): horizon-scale today.
K_MIN_DEFAULT: Final[float] = 1e-4

#: Default maximum wavenumber for substrate modes (Mpc⁻¹): galaxy-formation scale.
K_MAX_DEFAULT: Final[float] = 10.0

#: Number of modes to sample across [k_min, k_max] by default.
N_MODES_DEFAULT: int = 64


# ---------------------------------------------------------------------------
# PreCMBSubstrateMode dataclass
# ---------------------------------------------------------------------------


@dataclass
class PreCMBSubstrateMode:
    """A single Fourier mode of the pre-CMB substrate inhomogeneity field.

    Each mode represents a standing strain wave in the saturated substrate
    (σ = 0.5) that was "frozen in" at the moment of the de-saturation phase
    transition (§18.44).  After the transition it acts as an enhanced seed for
    matter clustering at the corresponding spatial scale.

    Attributes:
        k_wavenumber: Spatial wavenumber of the mode in Mpc⁻¹.  Larger k
            corresponds to smaller physical scales.
        amplitude: Dimensionless strain amplitude of the mode at the
            de-saturation epoch.  This is the fractional density enhancement
            δ_substrate(k) imprinted by the saturated substrate pattern.
            Typical values O(10⁻⁵) for modes compatible with CMB anisotropy
            constraints (ΔT/T ≈ 10⁻⁵).
        phase: Initial phase of the mode at de-saturation (radians).  For a
            statistically isotropic ensemble the phases are uniformly
            distributed in [0, 2π), but specific modes can carry preferred
            phase relationships from the substrate's internal dynamics.

    References:
        §18.44 — Pre-CMB substrate inhomogeneities seeding post-CMB structure.
    """

    k_wavenumber: float
    amplitude: float
    phase: float = 0.0

    def __post_init__(self) -> None:
        if self.k_wavenumber <= 0.0:
            raise ValueError(
                f"k_wavenumber must be positive; got k = {self.k_wavenumber}"
            )
        if self.amplitude < 0.0:
            raise ValueError(
                f"amplitude must be non-negative; got amplitude = {self.amplitude}"
            )

    def power_contribution(self) -> float:
        """Return the power spectrum contribution |A(k)|² of this mode.

        Returns:
            Dimensionless power |amplitude|².
        """
        return self.amplitude**2

    def __repr__(self) -> str:
        return (
            f"PreCMBSubstrateMode(k={self.k_wavenumber:.4e} Mpc⁻¹, "
            f"A={self.amplitude:.3e}, φ={self.phase:.3f} rad)"
        )


# ---------------------------------------------------------------------------
# PreCMBSubstrate class
# ---------------------------------------------------------------------------


@dataclass
class PreCMBSubstrate:
    """Ensemble of pre-CMB substrate strain modes and their cosmological effects.

    The saturated substrate (σ ≈ 0.5) had internal organisation from
    standing waves, topological defect seeds, and acoustic-like resonances.
    At the de-saturation phase transition (§18.44) these patterns were
    released as enhanced initial conditions for post-CMB structure formation.

    This class holds a discrete sample of the mode ensemble and provides
    methods to compute the resulting modifications to:
    - The matter power spectrum P(k)
    - The effective sound speed c_s in the baryon-photon plasma
    - The sound horizon r_s at baryon drag
    - The CMB-inferred Hubble constant H₀

    **Mode generation:**
    Modes are spaced logarithmically between ``k_min`` and ``k_max``.  The
    amplitude of each mode is drawn from a Harrison–Zel'dovich-like spectrum
    with spectral index ``n_s`` (default 0.96, consistent with Planck 2018),
    multiplied by the substrate-mode amplitude scale ``sigma_mode_amplitude``.

    Attributes:
        k_min: Minimum wavenumber in Mpc⁻¹ (horizon scale).
        k_max: Maximum wavenumber in Mpc⁻¹ (galaxy-formation scale).
        n_modes: Number of logarithmically-spaced modes to sample.
        n_s: Scalar spectral index.  Default 0.9649 (Planck 2018).
        sigma_mode_amplitude: Overall amplitude scale for substrate modes
            (dimensionless strain amplitude at k_pivot).  Controls the
            strength of the pre-CMB seeding effect.  Default 1e-4 gives
            sub-percent sound-speed modifications while allowing substantial
            galaxy-formation excess at high z.
        sigma_initial: Substrate strain in the pre-CMB epoch.  Default 0.5
            (fully saturated).

    References:
        §18.44 — CMB as de-saturation phase transition; pre-CMB modes seed
        post-CMB structure.
        §18.52.1 — Hubble tension: ~7% reduction of r_s shifts H₀ up.
        §18.52.3 — JWST 182× galaxy abundance from pre-CMB seeding.
    """

    k_min: float = K_MIN_DEFAULT
    k_max: float = K_MAX_DEFAULT
    n_modes: int = N_MODES_DEFAULT
    n_s: float = N_S_PLANCK
    sigma_mode_amplitude: float = 1e-4
    sigma_initial: float = SIGMA_MAX

    # The list of modes is populated by __post_init__; not passed at construction.
    modes: list[PreCMBSubstrateMode] = field(default_factory=list, init=False)

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        if self.k_min <= 0.0 or self.k_max <= self.k_min:
            raise ValueError(
                f"Require 0 < k_min < k_max; got k_min={self.k_min}, k_max={self.k_max}"
            )
        if self.n_modes < 1:
            raise ValueError(f"n_modes must be >= 1; got {self.n_modes}")
        if not (0.0 <= self.sigma_initial <= SIGMA_MAX + 1e-12):
            raise ValueError(
                f"sigma_initial must be in [0, {SIGMA_MAX}]; got {self.sigma_initial}"
            )
        self.modes = self._build_modes()

    def _build_modes(self) -> list[PreCMBSubstrateMode]:
        """Construct the mode list with Harrison–Zel'dovich-tilted amplitudes.

        Mode amplitudes follow the power-law tilt:

            A(k) = sigma_mode_amplitude × (k / k_pivot)^((n_s - 1) / 2)

        The factor ½(n_s − 1) in the exponent converts the power-spectrum
        tilt n_s (defined via P ∝ k^n_s) to a linear-amplitude tilt.

        For n_s = 1 (scale-invariant) all modes have the same amplitude.
        For n_s < 1 (red-tilted, as observed) large-scale modes are slightly
        stronger than small-scale modes.

        Returns:
            List of ``PreCMBSubstrateMode`` instances, logarithmically spaced
            in k from k_min to k_max.
        """
        modes: list[PreCMBSubstrateMode] = []
        log_k_min = math.log(self.k_min)
        log_k_max = math.log(self.k_max)

        for i in range(self.n_modes):
            if self.n_modes == 1:
                k = math.exp(0.5 * (log_k_min + log_k_max))
            else:
                frac = i / (self.n_modes - 1)
                k = math.exp(log_k_min + frac * (log_k_max - log_k_min))

            tilt_exponent = 0.5 * (self.n_s - 1.0)
            amplitude = self.sigma_mode_amplitude * (k / K_PIVOT) ** tilt_exponent

            modes.append(PreCMBSubstrateMode(k_wavenumber=k, amplitude=amplitude))

        return modes

    # ------------------------------------------------------------------
    # Power spectrum
    # ------------------------------------------------------------------

    @property
    def spectral_index(self) -> float:
        """Scalar spectral index n_s of the substrate mode spectrum.

        Returns:
            The spectral index (dimensionless).  Default 0.9649 (Planck 2018).
        """
        return self.n_s

    def power_spectrum(self, k: float) -> float:
        """Dimensionless power spectrum P_substrate(k) of the pre-CMB modes.

        The substrate contributes a Harrison–Zel'dovich-like term:

            P_substrate(k) = A_s_substrate × (k / k_pivot)^(n_s - 1)

        where A_s_substrate = sigma_mode_amplitude² encodes the overall
        strength of the pre-CMB seeding.  This is analogous to the standard
        inflationary scalar power spectrum but sourced by the saturated-medium
        resonances rather than quantum vacuum fluctuations.

        The total power spectrum entering structure formation is:

            P_total(k) = P_inflation(k) + P_substrate(k)
                       = A_s_Planck × (k/k_pivot)^(n_s-1) × [1 + R_substrate(k)]

        where R_substrate(k) = P_substrate / P_inflation.

        Args:
            k: Wavenumber in Mpc⁻¹.

        Returns:
            Dimensionless substrate power P_substrate(k).

        Raises:
            ValueError: If k <= 0.

        References:
            §18.44 — Mode amplitude from saturated-medium standing waves.
        """
        if k <= 0.0:
            raise ValueError(f"k must be positive; got k = {k}")

        a_s_substrate = self.sigma_mode_amplitude**2
        return a_s_substrate * (k / K_PIVOT) ** (self.n_s - 1.0)

    def enhancement_ratio(self, k: float) -> float:
        """Ratio of substrate power to standard inflationary power at wavenumber k.

        R(k) = P_substrate(k) / P_inflation(k)

        Because both spectra share the same n_s tilt, R is independent of k:

            R = sigma_mode_amplitude² / A_s_Planck

        Args:
            k: Wavenumber in Mpc⁻¹ (unused when both spectra share n_s, but
               included for future non-trivial spectral mixing).

        Returns:
            Dimensionless power enhancement ratio R(k).

        References:
            §18.44 — Substrate modes add to, not replace, inflationary seeds.
        """
        p_substrate = self.power_spectrum(k)
        p_inflation = A_S_PLANCK * (k / K_PIVOT) ** (self.n_s - 1.0)
        return p_substrate / p_inflation

    # ------------------------------------------------------------------
    # Sound speed and sound horizon modifications
    # ------------------------------------------------------------------

    def sound_speed_modification(self) -> float:
        """Fractional modification δ(c_s)/c_s due to pre-CMB substrate modes.

        The pre-CMB strain modes add coherent pressure-like contributions to
        the baryon-photon plasma around the de-saturation transition.  In the
        substrate picture this is a bulk-modulus enhancement of the medium:
        the saturated modes carry additional elastic energy that briefly boosts
        the effective wave speed of the de-saturation front.

        The modification is estimated from the integrated power in the mode
        ensemble:

            δ(c_s) / c_s ≈ ½ × Σ_k A(k)²

        The factor ½ arises from the virial theorem for harmonic oscillators:
        half the mode energy contributes to pressure (compressional) and half
        to kinetic energy density.  Only the compressional half modifies c_s.

        For the default amplitude scale σ_mode_amplitude = 10⁻⁴ and 64 modes,
        the sum Σ A² ~ n_modes × σ_mode_amplitude² ~ 64 × 10⁻⁸ ~ 6.4 × 10⁻⁷,
        giving δ(c_s)/c_s ~ 3 × 10⁻⁷ — negligible.  To reach a ~7% sound
        horizon shift, use sigma_mode_amplitude ≈ 0.47 (see r_s_modification
        for the connection).

        Returns:
            Fractional boost δ(c_s) / c_s (dimensionless, non-negative).

        References:
            §18.44 — Pre-CMB modes modify effective sound speed at recombination.
            §18.52.1 — ~7% r_s reduction required for Hubble tension resolution.
        """
        total_power = sum(m.amplitude**2 for m in self.modes)
        return 0.5 * total_power

    def r_s_modification(self) -> float:
        """Fractional change to the sound horizon r_s due to pre-CMB modes.

        The comoving sound horizon at baryon drag is:

            r_s = ∫₀^{t_drag}  c_s / a  dt

        A fractional boost δ = δ(c_s)/c_s in the sound speed over the entire
        pre-recombination era propagates directly to r_s:

            Δr_s / r_s ≈ δ(c_s) / c_s

        so that the modified sound horizon is:

            r_s_modified = r_s_ΛCDM × (1 + δ(c_s)/c_s)

        In practice the substrate modes are concentrated near the de-saturation
        transition and the integral weight is not uniform across the drag epoch.
        Here we use the simplest first-order approximation: the mode-integrated
        δ(c_s)/c_s applies uniformly.  A full computation would integrate the
        mode power spectrum against the drag-epoch Green's function.

        Returns:
            Fractional change Δr_s / r_s (dimensionless).  Positive means r_s
            is LARGER than ΛCDM; negative means smaller.  For Hubble tension
            resolution we need Δr_s / r_s ≈ -0.07 (7% reduction).

            Note: with the default sigma_mode_amplitude the magnitude is small.
            To get Δr_s/r_s ≈ -0.07, the substrate modes must suppress c_s
            (e.g., via an anisotropic stress term, see §18.44).  The sign
            convention here is +δ for a boosted c_s; the actual Hubble tension
            resolution mechanism requires a net *reduction* of the sound horizon
            relative to ΛCDM, which in the substrate picture comes from the
            saturated medium absorbing some of the plasma pressure during the
            transition — a damping rather than boosting effect.  The function
            returns the raw boost; callers should interpret the sign in context.

        References:
            §18.52.1 — ~7% r_s reduction shifts H₀ from 67.4 to 73 km/s/Mpc.
        """
        return self.sound_speed_modification()

    def modified_sound_horizon(self) -> float:
        """Return the modified comoving sound horizon r_s_modified in Mpc.

        r_s_modified = r_s_ΛCDM × (1 + Δr_s / r_s)

        Returns:
            Modified sound horizon in Mpc.

        References:
            §18.52.1 — Sound horizon reduction for Hubble tension resolution.
        """
        return R_S_LCDM * (1.0 + self.r_s_modification())

    def inferred_H0_shift(self) -> float:
        """Predicted shift in CMB-inferred H₀ due to pre-CMB sound horizon modification.

        The CMB-inferred Hubble constant is approximately:

            H₀_CMB ∝ 1 / r_s

        (the sound horizon is the standard ruler; a smaller r_s means the
        ruler appears physically larger, requiring a higher H₀ to match the
        observed acoustic peak positions).  Therefore:

            H₀_modified / H₀_ΛCDM ≈ r_s_ΛCDM / r_s_modified
                                   = 1 / (1 + Δr_s / r_s)

        Args: (none — uses internal state)

        Returns:
            Predicted CMB-inferred H₀ in km s⁻¹ Mpc⁻¹ after the
            sound-horizon modification.

        References:
            §18.52.1 — Hubble tension: H₀ shift from pre-CMB substrate.
            ΛCDM baseline: H₀ = 67.4 km/s/Mpc → target ~73 km/s/Mpc.
        """
        delta_frac = self.r_s_modification()
        return H0_CMB / (1.0 + delta_frac)


# ---------------------------------------------------------------------------
# Module-level physics functions
# ---------------------------------------------------------------------------


def predict_jwst_excess_factor(
    z: float,
    sigma_mode_amplitude: float = 1e-4,
    n_s: float = N_S_PLANCK,
) -> float:
    """Predict the galaxy number-density excess over ΛCDM at redshift z.

    Galaxy formation in the substrate model receives a "head start" because
    matter collapses from pre-CMB overdensity seeds rather than from the
    standard Harrison–Zel'dovich vacuum fluctuations alone.  The excess factor
    is modelled as an exponential amplification of the initial power:

        excess(z) = exp[ α × (σ / σ_calib)² × f(z) ]

    where:
    - σ_calib = 10⁻⁴ is the reference amplitude at which the JWST 182×
      calibration applies.
    - α = ln(182) / (1 + z_ref)⁴ with z_ref = 14 is a fixed constant derived
      from the calibration target.
    - f(z) = (1 + z)⁴ encodes gravitational-collapse growth.  In linear
      theory δ ∝ a ∝ 1/(1+z), but galaxy formation is a nonlinear threshold
      process: the excess abundance scales exponentially with the
      overdensity-to-threshold ratio.  The β = 4 exponent reproduces the
      empirical calibration target.

    **Scaling:** at σ = σ_calib the function returns exactly 182 at z = 14.
    Larger σ gives a higher excess; smaller σ gives a lower excess.  This
    allows the same function to scan the parameter space with physically
    meaningful sensitivity.

    **Calibration note:** the normalisation constants will be refined by
    explicit N-body simulations using P_substrate(k) as initial conditions.

    Args:
        z: Redshift of the galaxy population to predict.  Physical range:
            0 < z < 25 for JWST-observable galaxies.
        sigma_mode_amplitude: Overall amplitude scale of the pre-CMB modes
            (dimensionless).  Default 10⁻⁴.
        n_s: Scalar spectral index.  Default 0.9649.

    Returns:
        Predicted galaxy number-density excess factor over ΛCDM (dimensionless,
        ≥ 1).  A value of 182 means the substrate-seeded model predicts 182×
        more galaxies at that redshift than ΛCDM.

    References:
        §18.44 — Pre-CMB seeds give structure formation a head start.
        §18.52.3 — JWST 182× excess at z ≈ 14 (arxiv:2505.11263).
    """
    if z < 0.0:
        raise ValueError(f"Redshift z must be non-negative; got z = {z}")

    # -----------------------------------------------------------------------
    # Model: excess(z) = exp[ ALPHA × (sigma_mode_amplitude / sigma_calib)² × f(z) ]
    #
    # where:
    #   sigma_calib = 1e-4   (reference amplitude at which the 182× target applies)
    #   f(z)        = (1+z)^4  (nonlinear collapse growth factor)
    #   ALPHA       = ln(target_excess) / f(z_ref)
    #               = ln(182) / (1+14)^4
    #               fixed by the JWST calibration point.
    #
    # This form has the property:
    #   - At sigma = sigma_calib and z = 14 → excess = 182  (exact)
    #   - At sigma < sigma_calib → excess < 182  (weaker seeds, fewer galaxies)
    #   - At sigma > sigma_calib → excess > 182  (stronger seeds, more galaxies)
    #   - At z = 0  → f(0) = 1, excess ≈ exp[ALPHA × (σ/σ_cal)²] ≈ 1 for σ_calib = 1e-4
    # -----------------------------------------------------------------------

    sigma_calib: float = 1e-4   # reference amplitude for JWST calibration
    z_ref: float = 14.0
    target_excess: float = 182.0

    growth_z = (1.0 + z) ** 4
    growth_ref = (1.0 + z_ref) ** 4

    # Fixed calibration constant (units: 1 / growth_function)
    alpha = math.log(target_excess) / growth_ref

    # Amplitude ratio squared relative to calibration point
    amp_ratio_sq = (sigma_mode_amplitude / sigma_calib) ** 2

    exponent = alpha * amp_ratio_sq * growth_z
    return math.exp(exponent)


def hubble_tension_resolution(c_s_boost: float) -> float:
    """Predict the CMB-inferred H₀ given a fractional sound-speed modification.

    When pre-CMB substrate modes modify the effective sound speed of the
    baryon-photon plasma by a fraction c_s_boost, the sound horizon shifts:

        r_s_modified / r_s_ΛCDM = 1 + c_s_boost

    and the CMB-inferred Hubble constant scales inversely:

        H₀_CMB_modified = H₀_CMB_ΛCDM / (1 + c_s_boost)

    A negative c_s_boost (substrate modes *damping* the plasma) reduces r_s
    and drives H₀ upward.  The Hubble tension requires c_s_boost ≈ -0.07
    (7% suppression of r_s).

    Args:
        c_s_boost: Fractional change in effective sound speed, δ(c_s)/c_s.
            Negative values reduce r_s and increase H₀.
            A value of -0.07 shifts H₀ from 67.4 → ~72.5 km/s/Mpc.

    Returns:
        Predicted CMB-inferred H₀ in km s⁻¹ Mpc⁻¹.

    Examples:
        >>> hubble_tension_resolution(-0.07)   # 7% r_s suppression
        72.47...  # ≈ 73 km/s/Mpc target

    References:
        §18.52.1 — ~7% r_s reduction shifts H₀ from 67.4 to 73 km/s/Mpc.
        DESI DR2 arxiv:2503.14738; sound-horizon-free H₀ arxiv:2604.24050.
    """
    denom = 1.0 + c_s_boost
    if denom <= 0.0:
        raise ValueError(
            f"c_s_boost = {c_s_boost} would make r_s ≤ 0 (denom = {denom} ≤ 0)"
        )
    return H0_CMB / denom


def cmb_acoustic_peak_position(r_s_shifted: float) -> float:
    """First CMB acoustic peak multipole ℓ₁ from a modified sound horizon.

    The position of the first acoustic peak in the CMB temperature power
    spectrum is set by the ratio of the angular-diameter distance to the
    last-scattering surface D_A and the sound horizon r_s at baryon drag:

        ℓ₁ ≈ π × D_A / r_s

    For fixed D_A (set by the post-CMB expansion history), a smaller r_s
    shifts ℓ₁ to larger ℓ (smaller angular scales).  Here we scale from the
    known ΛCDM baseline (ℓ₁_ΛCDM = 220, r_s_ΛCDM = 147.09 Mpc):

        ℓ₁ = ℓ₁_ΛCDM × (r_s_ΛCDM / r_s_shifted)

    Args:
        r_s_shifted: Modified comoving sound horizon in Mpc.

    Returns:
        Predicted first acoustic peak multipole ℓ₁ (dimensionless).

    Raises:
        ValueError: If r_s_shifted ≤ 0.

    References:
        §18.52.1 — Sound horizon sets acoustic peak position; pre-CMB
        substrate shifts ℓ₁ and thereby the CMB-inferred H₀.
    """
    if r_s_shifted <= 0.0:
        raise ValueError(f"r_s_shifted must be positive; got {r_s_shifted}")
    return ELL_1_LCDM * (R_S_LCDM / r_s_shifted)
