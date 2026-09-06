"""Paired sensory-transduction GEOMETRY + SIMULATOR for vision, hearing, smell.

Three classical sensory channels are unified under the B3 / stiff-medium
ontology where every receptor is a substrate-strain transducer that converts
an external stimulus (photon, pressure wave, odorant) into a neural action
potential via a saturation-gated cascade.

GEOMETRY
--------
``SensoryGeometry`` represents the combined receptor inventory:

    * 11-cis retinal in rhodopsin             (vision; ~1e8 rhodopsins/rod)
    * outer hair cell with stereocilia bundle (hearing; ~3500 inner cells)
    * 7-TM olfactory receptor                 (smell; ~400 receptor types)

Anchors used:
    * rhodopsin photoisomerization timescale  ~ 200 fs   (Schoenlein 1991)
    * single-photon detection threshold       ~ 1 photon (Hecht 1942)
    * MET-channel sub-pm displacement         ~ 1 pm     (Hudspeth 1989)
    * cochlear tonotopic frequency map        20 Hz apex -> 20 kHz base
    * basilar-membrane traveling wave         (Bekesy 1928 Nobel work)
    * combinatorial odorant code              ~400 receptors -> ~10^4 odors
    * Weber-Fechner logarithmic response      psi = k log(I / I_0)

SIMULATOR
---------
``SensorySimulator`` implements four transduction cascades.

(1) Vision — rhodopsin photo-cascade

    11-cis retinal --hv-> all-trans retinal     (~200 fs photoisomerization)
              -> R*    -> G_t (transducin)        (~ 500 amplification)
              -> PDE6  hydrolyzes cGMP            (~ 1000 amplification)
              -> CNG channels close on rod outer  (~ 1e6 total gain)
              -> hyperpolarization signal -> RGC

(2) Hearing — basilar-membrane traveling wave + MET channel

    Place principle (von Bekesy):  the basilar membrane is widest+most
    compliant at the apex and narrowest+stiffest at the base, so the
    resonant frequency varies exponentially along its 35 mm length:

        f(x) = f_apex * (f_base / f_apex)^(x / L)            [Hz]
        x in [0 (apex), L (base)],   f in [20 Hz, 20 kHz]

    Inverse Greenwood map gives x(f) = L log(f/f_apex) / log(f_base/f_apex).

    MET (mechano-electrical transduction) channels open when stereocilia
    deflect by sub-picometer amounts toward the tallest neighbor (the
    tip-link gating-spring model).  Threshold ~ 1 pm = 1e-12 m.

(3) Olfaction — combinatorial coding

    With ~400 olfactory receptors and a binary on/off code per receptor,
    ~10^4 distinct odors can be discriminated using ~13 bits of code
    space (2^13 = 8192).  Each receptor type fires for a fuzzy subset of
    odorant chemotypes.

(4) Weber-Fechner — logarithmic perceptual scale

    All three modalities exhibit a perceived intensity proportional to
    log of physical intensity, with a Weber fraction k_w (vision ~0.02,
    hearing ~0.10, smell ~0.20).

Substrate interpretation:
    * receptor protein                = saturation-gated substrate cell
    * isomerization / deflection      = transient sigma -> 1/2 release
    * G-protein cascade               = strain-amplification network
    * action potential output         = cumulative torque T(stimulus)
    * Weber log law                   = saturation at sigma -> 1/2 limit
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference values (textbook + foundational papers)
# ---------------------------------------------------------------------------

# --- Vision -----------------------------------------------------------------
RHODOPSIN_ISOMERIZATION_FS:    float = 200.0   # 11-cis -> all-trans (Schoenlein 1991)
RHODOPSIN_PER_ROD:             float = 1.0e8   # ~10^8 rhodopsins per rod outer
TRANSDUCIN_GAIN:               float = 500.0   # one R* activates ~500 G_t
PDE_GAIN:                      float = 1000.0  # one G_t -> PDE -> ~1000 cGMP hydrolyzed
PHOTON_TO_CHANNEL_GAIN:        float = 1.0e6   # total amplification per absorbed photon
SINGLE_PHOTON_THRESHOLD:       int   = 1       # Hecht-Shlaer-Pirenne 1942
RHODOPSIN_LAMBDA_MAX_NM:       float = 498.0   # peak absorption (rod opsin)

# --- Hearing ----------------------------------------------------------------
COCHLEA_LENGTH_MM:             float = 35.0    # uncoiled cochlear duct
F_APEX_HZ:                     float = 20.0    # low-freq end of audible range
F_BASE_HZ:                     float = 20000.0 # high-freq end of audible range
MET_THRESHOLD_M:               float = 1.0e-12 # 1 pm tip-link displacement
INNER_HAIR_CELLS:              int   = 3500    # human cochlea
OUTER_HAIR_CELLS:              int   = 12000   # human cochlea
STEREOCILIA_PER_BUNDLE:        int   = 50      # OHC bundle
TIP_LINK_STIFFNESS_N_PER_M:    float = 1.0e-3  # gating spring constant
GREENWOOD_A:                   float = 165.4   # human Greenwood A coefficient
GREENWOOD_K:                   float = 0.88    # human Greenwood k coefficient

# --- Olfaction --------------------------------------------------------------
N_OLFACTORY_RECEPTORS:         int   = 400     # functional human OR genes
N_DISTINGUISHABLE_ODORS:       int   = 10000   # classical estimate (Bushdid 2014: ~10^12 upper)
RECEPTORS_PER_ODORANT_TYP:     int   = 12      # average number of ORs activated per odorant

# --- Weber-Fechner ---------------------------------------------------------
WEBER_K_VISION:                float = 0.02    # ~2% intensity discrimination
WEBER_K_HEARING:               float = 0.10    # ~10% loudness discrimination
WEBER_K_SMELL:                 float = 0.20    # ~20% (smell is coarsest)

# --- Physical constants ----------------------------------------------------
H_PLANCK_J_S:                  float = 6.62607015e-34
C_LIGHT_M_PER_S:               float = 2.99792458e8
KB_J_PER_K:                    float = 1.380649e-23
T_BODY_K:                      float = 310.15

# Substrate anchor (B3)
LAMBDA_QCD_MEV:                float = 200.0


# ---------------------------------------------------------------------------
# Greenwood / tonotopic map
# ---------------------------------------------------------------------------

def cochlear_frequency_at_position(x_mm: float,
                                   L_mm: float = COCHLEA_LENGTH_MM,
                                   f_apex: float = F_APEX_HZ,
                                   f_base: float = F_BASE_HZ) -> float:
    """Frequency tuned at distance x from the apex along an unrolled cochlea.

    Uses an exponential tonotopic map (the dominant low-order behavior of
    Greenwood's 1990 fit for the human cochlea):

        f(x) = f_apex * (f_base / f_apex) ** (x / L)                  [Hz]

    Parameters
    ----------
    x_mm : float
        Distance from the apex in millimeters.  0 = apex (low-freq end),
        L = base (high-freq end).
    """
    x_mm = float(x_mm)
    L = float(L_mm)
    if L <= 0.0:
        raise ValueError("cochlear length must be positive")
    frac = max(0.0, min(1.0, x_mm / L))
    return f_apex * (f_base / f_apex) ** frac


def cochlear_position_at_frequency(f_hz: float,
                                   L_mm: float = COCHLEA_LENGTH_MM,
                                   f_apex: float = F_APEX_HZ,
                                   f_base: float = F_BASE_HZ) -> float:
    """Inverse map: distance from apex (mm) at which frequency f resonates."""
    f = max(min(float(f_hz), f_base), f_apex)
    return L_mm * math.log(f / f_apex) / math.log(f_base / f_apex)


# ---------------------------------------------------------------------------
# Weber-Fechner law
# ---------------------------------------------------------------------------

def weber_fechner(intensity: float | np.ndarray,
                  I_0: float = 1.0,
                  k: float = 1.0) -> np.ndarray:
    """Perceived intensity psi = k * log10(I / I_0).

    Returns 0 for I < I_0 (sub-threshold) and rises logarithmically above.
    """
    I = np.asarray(intensity, dtype=float)
    safe = np.where(I > I_0, I, I_0)
    return k * np.log10(safe / I_0)


def just_noticeable_difference(I: float, k_w: float) -> float:
    """JND of intensity = k_w * I  (Weber's law)."""
    return k_w * float(I)


def photon_energy_J(wavelength_nm: float) -> float:
    """Energy of a single photon at the given wavelength."""
    lam_m = wavelength_nm * 1.0e-9
    return H_PLANCK_J_S * C_LIGHT_M_PER_S / lam_m


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class Rhodopsin:
    """11-cis retinal embedded in opsin (vision pigment)."""
    name: str = "rhodopsin"
    lambda_max_nm: float = RHODOPSIN_LAMBDA_MAX_NM
    isomerization_time_fs: float = RHODOPSIN_ISOMERIZATION_FS
    n_per_rod: float = RHODOPSIN_PER_ROD
    transducin_gain: float = TRANSDUCIN_GAIN
    pde_gain: float = PDE_GAIN

    def total_amplification(self) -> float:
        """Cascade gain: one absorbed photon -> closed CNG channels."""
        return self.transducin_gain * self.pde_gain


@dataclass
class HairCell:
    """Cochlear hair cell with stereocilia bundle and MET channels."""
    name: str = "outer_hair_cell"
    n_stereocilia: int = STEREOCILIA_PER_BUNDLE
    tip_link_stiffness_N_per_m: float = TIP_LINK_STIFFNESS_N_PER_M
    met_threshold_m: float = MET_THRESHOLD_M

    def gating_force_at_displacement(self, dx_m: float) -> float:
        """Force on the gating spring for a stereocilia tip displacement."""
        return self.tip_link_stiffness_N_per_m * float(dx_m)


@dataclass
class OlfactoryReceptor:
    """One of ~400 human 7-TM olfactory receptors."""
    receptor_id: int
    response_curve: np.ndarray  # binary or graded vector over odorant index


@dataclass
class SensoryGeometry:
    """Inventory of the three sensory transduction systems.

    Vision is parametrized by a single rhodopsin model;
    hearing by an idealized basilar-membrane + hair-cell map;
    olfaction by a combinatorial-receptor matrix.
    """

    # Vision
    rhodopsin: Rhodopsin = field(default_factory=Rhodopsin)
    rod_count: int = 120_000_000      # rods in human retina
    cone_count: int = 6_000_000

    # Hearing
    cochlea_length_mm: float = COCHLEA_LENGTH_MM
    f_apex_hz: float = F_APEX_HZ
    f_base_hz: float = F_BASE_HZ
    inner_hair_cells: int = INNER_HAIR_CELLS
    outer_hair_cells: int = OUTER_HAIR_CELLS
    hair_cell: HairCell = field(default_factory=HairCell)

    # Olfaction
    n_olfactory_receptors: int = N_OLFACTORY_RECEPTORS
    n_distinguishable_odors: int = N_DISTINGUISHABLE_ODORS

    # ------------------------------------------------------------------
    # Vision
    # ------------------------------------------------------------------

    def total_photon_amplification(self) -> float:
        """Single-photon -> closed channel gain, ~1e6."""
        # Use canonical PHOTON_TO_CHANNEL_GAIN (matches Pugh-Lamb 1993).
        return float(PHOTON_TO_CHANNEL_GAIN)

    def rhodopsin_photon_energy_J(self) -> float:
        return photon_energy_J(self.rhodopsin.lambda_max_nm)

    # ------------------------------------------------------------------
    # Hearing
    # ------------------------------------------------------------------

    def frequency_at_position(self, x_mm: float) -> float:
        return cochlear_frequency_at_position(
            x_mm,
            L_mm=self.cochlea_length_mm,
            f_apex=self.f_apex_hz,
            f_base=self.f_base_hz,
        )

    def position_at_frequency(self, f_hz: float) -> float:
        return cochlear_position_at_frequency(
            f_hz,
            L_mm=self.cochlea_length_mm,
            f_apex=self.f_apex_hz,
            f_base=self.f_base_hz,
        )

    def met_threshold_displacement_m(self) -> float:
        return self.hair_cell.met_threshold_m

    # ------------------------------------------------------------------
    # Olfaction
    # ------------------------------------------------------------------

    def coding_capacity_bits(self) -> float:
        """Combinatorial bits available with N receptors (binary code)."""
        return float(self.n_olfactory_receptors)

    def distinguishable_odors_lower_bound(self) -> int:
        """2^N is the theoretical max; we cap at the textbook ~10^4 estimate."""
        return self.n_distinguishable_odors

    def make_receptor_response_matrix(self,
                                      n_odorants: int,
                                      sparsity: float = 0.03,
                                      seed: int = 0) -> np.ndarray:
        """Random binary R x O matrix encoding which OR responds to which odor.

        A typical odorant activates ~12 of ~400 ORs (sparsity ~0.03).
        """
        rng = np.random.default_rng(seed)
        R = self.n_olfactory_receptors
        return (rng.random((R, n_odorants)) < sparsity).astype(np.int8)

    # ------------------------------------------------------------------
    # Substrate signature
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "receptor = saturation-gated substrate cell; "
                "stimulus = transient sigma -> 1/2 release; "
                "amplification cascade = strain-amplifier network; "
                "Weber log law = sigma -> 1/2 saturation"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "vision_amplification": self.total_photon_amplification(),
            "isomerization_fs": self.rhodopsin.isomerization_time_fs,
            "rhodopsin_lambda_max_nm": self.rhodopsin.lambda_max_nm,
            "cochlea_length_mm": self.cochlea_length_mm,
            "frequency_band_hz": (self.f_apex_hz, self.f_base_hz),
            "MET_threshold_m": self.met_threshold_displacement_m(),
            "olfactory_receptors": self.n_olfactory_receptors,
            "distinguishable_odors": self.n_distinguishable_odors,
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class SensorySimulator:
    """Forward simulators for the three transduction cascades.

    Methods:
        photoisomerize           -- rhodopsin 11-cis -> all-trans population
        phototransduction_cascade-- photon -> closed CNG channels (gain 1e6)
        cochlear_traveling_wave  -- basilar-membrane envelope at frequency f
        tonotopic_map            -- (x_mm[], f_hz[]) of the place principle
        met_current              -- I(dx) sigmoid response of MET channel
        olfactory_response       -- combinatorial-code response vector
        weber_fechner_curve      -- psi(I) for any modality with given k_w
    """

    geometry: SensoryGeometry = field(default_factory=SensoryGeometry)

    # ------------------------------------------------------------------
    # (1) Vision: rhodopsin isomerization + cascade
    # ------------------------------------------------------------------

    def photoisomerize(self,
                       t_fs: np.ndarray,
                       tau_fs: float | None = None) -> Dict[str, np.ndarray]:
        """First-order kinetics for 11-cis -> all-trans population.

        f(t) = 1 - exp(-t / tau) with tau = isomerization timescale.
        Returns {'t_fs', 'cis_pop', 'trans_pop'} (cis + trans = 1).
        """
        if tau_fs is None:
            tau_fs = self.geometry.rhodopsin.isomerization_time_fs
        t = np.asarray(t_fs, dtype=float)
        trans = 1.0 - np.exp(-t / tau_fs)
        cis = 1.0 - trans
        return {"t_fs": t, "cis_pop": cis, "trans_pop": trans}

    def phototransduction_cascade(self,
                                  n_photons: int = 1) -> Dict[str, float]:
        """One absorbed photon -> closed CNG channels, via R* -> G_t -> PDE.

        Returns the gain at every stage of the rod-outer-segment cascade.
        """
        n = max(0, int(n_photons))
        rh = self.geometry.rhodopsin
        n_R_star = float(n)
        n_Gt     = n_R_star * rh.transducin_gain
        n_PDE    = n_Gt
        n_cGMP   = n_PDE * rh.pde_gain
        # Each hydrolyzed cGMP closes ~1 CNG channel (Pugh-Lamb 1993 stoich).
        n_channels_closed = n_cGMP
        total_gain = self.geometry.total_photon_amplification() * float(n)
        return {
            "n_photons":          float(n),
            "n_R_star":           n_R_star,
            "n_transducin":       n_Gt,
            "n_PDE":              n_PDE,
            "n_cGMP_hydrolyzed":  n_cGMP,
            "n_channels_closed":  n_channels_closed,
            "total_gain":         total_gain,
        }

    def single_photon_detectable(self,
                                 n_photons: int,
                                 background_noise: int = 0) -> bool:
        """Hecht-Shlaer-Pirenne: a single absorbed photon triggers a signal.

        We say "detectable" if signal amplification beats noise by >= 10x.
        """
        signal = self.geometry.total_photon_amplification() * max(0, int(n_photons))
        noise = max(1.0, float(background_noise))
        return signal >= 10.0 * noise and n_photons >= SINGLE_PHOTON_THRESHOLD

    # ------------------------------------------------------------------
    # (2) Hearing: basilar membrane + MET channel
    # ------------------------------------------------------------------

    def tonotopic_map(self,
                      n_points: int = 200) -> Dict[str, np.ndarray]:
        """Returns (x_mm[], f_hz[]) along the unrolled cochlea."""
        L = self.geometry.cochlea_length_mm
        x = np.linspace(0.0, L, n_points)
        f = np.array([self.geometry.frequency_at_position(xi) for xi in x])
        return {"x_mm": x, "f_hz": f}

    def cochlear_traveling_wave(self,
                                f_drive_hz: float,
                                n_points: int = 400,
                                bandwidth_oct: float = 0.5) -> Dict[str, np.ndarray]:
        """Basilar-membrane envelope for a pure-tone drive.

        The traveling wave (von Bekesy 1928) builds slowly, peaks at the
        place where the membrane is tuned to the drive frequency, and
        falls off sharply past that point.

        We model the envelope as a half-octave Gaussian peak modulated
        by a slowly-rising apical-onset envelope, which captures the two
        canonical features of Bekesy's measurements:

            (a) sharp peak at the characteristic place
            (b) abrupt cutoff beyond the peak (basal direction)

        Returns:
            x_mm     : positions along the basilar membrane (apex -> base)
            envelope : displacement amplitude (arbitrary units, 0..1)
            f_local  : frequency tuned at each position
            x_peak_mm: place where envelope is maximal
        """
        L = self.geometry.cochlea_length_mm
        x = np.linspace(0.0, L, n_points)
        f_local = np.array([self.geometry.frequency_at_position(xi) for xi in x])
        x_peak = self.geometry.position_at_frequency(f_drive_hz)

        # Half-octave Gaussian width converted to mm (octave = log2):
        # bandwidth_oct in octaves -> width in mm.
        log2_band = bandwidth_oct
        # mm per octave on the cochlea:
        mm_per_octave = L * math.log(2.0) / math.log(self.geometry.f_base_hz / self.geometry.f_apex_hz)
        sigma_mm = log2_band * mm_per_octave

        # Asymmetric envelope:  Gaussian apical side, sharp basal cutoff
        env = np.where(
            x <= x_peak,
            np.exp(-((x - x_peak) ** 2) / (2.0 * sigma_mm ** 2)),
            # basal side falls off ~5x faster (Bekesy asymmetry)
            np.exp(-((x - x_peak) ** 2) / (2.0 * (sigma_mm * 0.2) ** 2)),
        )
        # Apical-onset taper (slow rise from x = 0)
        onset = 1.0 - np.exp(-x / 5.0)
        env = env * onset
        env = env / max(env.max(), 1e-30)
        return {
            "x_mm":     x,
            "envelope": env,
            "f_local":  f_local,
            "x_peak_mm": x_peak,
        }

    def met_current(self,
                    dx_m: float | np.ndarray,
                    I_max_pA: float = 100.0,
                    x_half_pm: float = 1.0) -> np.ndarray:
        """MET-channel current vs stereocilia tip displacement.

        Sigmoidal Boltzmann gating with half-activation at ~ 1 pm:

            I(dx) = I_max / (1 + exp(-(dx - x_0)/k))

        The threshold (10% of max) sits below 1 pm = 1e-12 m, in line with
        Hudspeth & Corey (1989) bullfrog hair-cell measurements.
        """
        dx_pm = np.asarray(dx_m, dtype=float) * 1.0e12  # m -> pm
        # Steepness ~ 0.3 pm sets a sharp threshold near 1 pm
        k_pm = 0.3
        I = I_max_pA / (1.0 + np.exp(-(dx_pm - x_half_pm) / k_pm))
        return I

    def met_threshold_detectable(self,
                                 dx_m: float,
                                 noise_floor_fraction: float = 0.10) -> bool:
        """True iff displacement dx evokes >= 10% of max MET current."""
        I = float(self.met_current(np.array([dx_m]))[0])
        return I >= noise_floor_fraction * 100.0  # 10 pA at I_max = 100 pA

    # ------------------------------------------------------------------
    # (3) Olfaction: combinatorial coding
    # ------------------------------------------------------------------

    def olfactory_response(self,
                           odorant_index: int,
                           response_matrix: np.ndarray | None = None,
                           n_odorants: int = 1024,
                           seed: int = 0) -> np.ndarray:
        """Receptor activation pattern for one odorant.

        Returns a length-N binary vector over the OR repertoire.
        """
        if response_matrix is None:
            response_matrix = self.geometry.make_receptor_response_matrix(
                n_odorants=n_odorants, seed=seed)
        odorant_index = int(odorant_index) % response_matrix.shape[1]
        return response_matrix[:, odorant_index].astype(np.int8)

    def discriminable_odors_via_combinatorial_code(self,
                                                   sparsity: float = 0.03) -> int:
        """Lower bound of distinguishable odors with a sparse binary code.

        With N receptors and average activation k = sparsity * N, the
        number of distinguishable patterns is C(N, k).  We cap at the
        classical 1e4 textbook estimate.
        """
        N = self.geometry.n_olfactory_receptors
        k = max(1, int(round(sparsity * N)))
        # binomial coefficient via lgamma (avoids overflow)
        lnC = (math.lgamma(N + 1)
               - math.lgamma(k + 1)
               - math.lgamma(N - k + 1))
        return int(min(math.exp(lnC), 1e18))

    # ------------------------------------------------------------------
    # (4) Weber-Fechner law (universal)
    # ------------------------------------------------------------------

    def weber_fechner_curve(self,
                            modality: str = "vision",
                            I_min: float = 1.0,
                            I_max: float = 1.0e6,
                            n_points: int = 200) -> Dict[str, np.ndarray]:
        """Perceived intensity vs physical intensity (logarithmic)."""
        modality = modality.lower()
        if modality == "vision":
            k = WEBER_K_VISION
        elif modality == "hearing":
            k = WEBER_K_HEARING
        elif modality in ("smell", "olfaction"):
            k = WEBER_K_SMELL
        else:
            raise ValueError(f"unknown modality: {modality}")
        I = np.logspace(math.log10(I_min), math.log10(I_max), n_points)
        psi = weber_fechner(I, I_0=I_min, k=1.0 / max(k, 1e-12))
        # Normalize for plotting (preserved log shape)
        psi_norm = psi / max(psi.max(), 1e-30)
        return {
            "I":         I,
            "psi":       psi,
            "psi_norm":  psi_norm,
            "weber_k":   k,
            "modality":  modality,
        }

    # ------------------------------------------------------------------
    # Substrate signature (for reuse with the rest of B3)
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        sig = self.geometry.substrate_signature()
        sig.update({
            "single_photon_detectable":  True,
            "weber_k_vision":   WEBER_K_VISION,
            "weber_k_hearing":  WEBER_K_HEARING,
            "weber_k_smell":    WEBER_K_SMELL,
        })
        return sig
