"""Mitochondrial bioenergetics — paired GEOMETRY + SIM + VIZ.

Substrate ontology of mitochondria
----------------------------------
The B3 / stiff-medium framework treats mitochondrial bioenergetics as a
substrate-strain pump: electrons flowing down a redox potential ladder
(NADH → O₂) are converted into a proton chemical potential across the
inner membrane (Δμ_H+), which is then converted back into mechanical
torque on the F0F1 ATP synthase rotor — and the rotor torque is finally
converted into ATP-bond strain energy.

Three coupled energy currencies, three substrate transduction stages:

    e⁻ flow            ⇌  ΔE_redox         (Complex I, II, III, IV)
    H⁺ flow            ⇌  Δμ_H+ (PMF)      (chemiosmosis, Mitchell 1961)
    rotor torque       ⇌  ATP synthesis    (F0F1, Boyer 1979 / Walker 1997)

This module supplies a paired GEOMETRY + SIM that reproduces the textbook
quantitative chemistry of these three stages plus the cristae topology.

Topology
~~~~~~~~
Inner membrane is folded into cristae (sheet-like invaginations).  Surface
area ratio cristae/IM_smooth ≈ 5 — the substrate uses topology to multiply
its working surface for ETC complex packing.

ETC complexes (4 of them, embedded in inner membrane):
    Complex I   (NADH dehydrogenase)        4 H⁺/2 e⁻ pumped, NADH → Q
    Complex II  (succinate dehydrogenase)   0 H⁺/2 e⁻ pumped, FADH₂ → Q
    Complex III (cyt bc1)                   4 H⁺/2 e⁻ pumped, Q → cyt c
    Complex IV  (cyt c oxidase)             2 H⁺/2 e⁻ pumped, cyt c → O₂

Mobile carriers:
    ubiquinone Q  — lipid-soluble, shuttles 2 e⁻ from I/II to III
    cyt c         — water-soluble, shuttles 1 e⁻ from III to IV (×2)

Redox potentials (E°' at pH 7, V):
    NADH / NAD⁺           = -0.32 V
    FADH₂ / FAD           = -0.22 V (succinate / fumarate)
    Q / QH₂               = +0.04 V (lipid-soluble pool)
    cyt c (Fe³⁺/Fe²⁺)     = +0.25 V
    ½ O₂ / H₂O            = +0.82 V
Total drop NADH→O₂ = 1.14 V → ΔG = −nFΔE = −2 × 96485 × 1.14
                              = -220 kJ/mol per 2 e⁻

Proton motive force (PMF / Δμ_H+)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Δμ_H+ = ΔΨ + (2.303 RT/F) · (pH_out − pH_in)
          = -0.16 V + 0.06 V × (-1)         (typical mito: ΔpH=1, ΔΨ=-160 mV)
          = -0.22 V → quoted positive as Δμ_H+ ≈ 200 mV (Mitchell 1961)
ΔG per H⁺ = F · |Δμ_H+| = 96485 × 0.22 = 21 kJ/mol per H⁺

F0F1 ATP synthase
~~~~~~~~~~~~~~~~~
F0 rotor: c-ring with c=8 (mammals), 10 (yeast), 14 (chloroplasts)
F1 head: 3 catalytic β-subunits, each cycling through O→L→T states
Stoichiometry: c protons drive 360° rotation = 3 ATP synthesized
    → 3 H⁺/ATP (yeast c10),  2.7 H⁺/ATP (mammal c8),  4.7 H⁺/ATP (plants c14)
Boyer's binding-change mechanism: rotor torque transduces PMF energy
    into ATP bond strain (substrate axiom: mass = cumulative torque).

P/O ratio
~~~~~~~~~
    P/O = (H⁺ pumped per 2 e⁻) / (H⁺ per ATP)
    NADH path:    P/O = 10 / 4 = 2.5  ATP per NADH (mammal)
    FADH₂ path:   P/O =  6 / 4 = 1.5  ATP per FADH₂

Glucose oxidation total
~~~~~~~~~~~~~~~~~~~~~~~
Glycolysis:        2 ATP + 2 NADH (cyto)
Pyruvate→AcCoA:    2 NADH
TCA cycle (×2):    6 NADH + 2 FADH₂ + 2 GTP
Substrate-level:   4 ATP equivalents
Oxidative phosph:  10 NADH × 2.5 + 2 FADH₂ × 1.5 = 25 + 3 = 28 ATP
TOTAL:             ≈ 30 ATP per glucose (mammal, cyto NADH penalised by shuttle)
                   ≈ 32 ATP per glucose (idealised, malate-aspartate shuttle)

WHAT THIS MODULE DOES
---------------------
1. ``MitoGeometry``      — inner membrane + cristae + 4 ETC complexes +
                           ATP synthase placement, proton gradient profile.
2. ``MitoSimulator``     — ETC electron transport (NADH→Q→cyt c→O₂),
                           PMF ledger, F0F1 rotor kinetics, ATP synthesis
                           stoichiometry, P/O ratio, glucose accounting.
3. Standalone helpers    — proton_motive_force(), atp_per_nadh(),
                           glucose_to_atp(), c_ring_atp_ratio().

VISUALS PRODUCED BY ``render_mito``
-----------------------------------
106_etc_complexes.png   — 4 ETC complexes (I/II/III/IV) with redox arrows,
                          Q and cyt c shuttle paths, H⁺ pumping count.
107_atp_synthase.png    — F0F1 rotor structure (c-ring + γ-shaft + α₃β₃
                          head) + proton flux → ATP synthesis curve.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray


# ---------------------------------------------------------------------------
# Physical constants (SI; biochemistry units)
# ---------------------------------------------------------------------------

R_GAS = 8.314_462_618          # J / (mol · K)
F_FARADAY = 96_485.332_12      # C / mol
T_BODY = 310.15                # K (37°C, mammalian body temperature)
RT_F = R_GAS * T_BODY / F_FARADAY      # 0.0267 V at 310 K
LN10_RT_F = 2.302_585 * RT_F           # 0.0615 V at 310 K (Nernst slope)


# ---------------------------------------------------------------------------
# Reference values (textbook biochemistry)
# ---------------------------------------------------------------------------

# Redox potentials E°' (V at pH 7)
E_NADH_NAD     = -0.320
E_FADH2_FAD    = -0.220       # succinate / fumarate couple
E_Q_QH2        = +0.045       # ubiquinone pool
E_CYTC         = +0.254       # cytochrome c (Fe³⁺ / Fe²⁺)
E_O2_H2O       = +0.815

# Proton stoichiometry per 2 e⁻ traversal (canonical mammalian values)
H_PER_2E_COMPLEX_I   = 4
H_PER_2E_COMPLEX_II  = 0
H_PER_2E_COMPLEX_III = 4
H_PER_2E_COMPLEX_IV  = 2

# c-ring sizes for the F0 rotor across organisms
C_RING_MAMMAL   = 8           # bovine, human (Watt 2010)
C_RING_YEAST    = 10          # S. cerevisiae (Stock 1999)
C_RING_BACTERIA = 10          # E. coli (most bacteria)
C_RING_PLANT    = 14          # spinach chloroplast

# ATP per c-ring rotation = 3 (3 catalytic β subunits)
ATP_PER_REV = 3

# Glucose accounting (mammalian, with malate-aspartate shuttle)
NADH_PER_GLUCOSE = 10         # 2 cyto + 2 pyruvate→AcCoA + 6 TCA
FADH2_PER_GLUCOSE = 2         # 2 from TCA succinate→fumarate
SUBSTRATE_LEVEL_ATP = 4       # 2 from glycolysis + 2 from TCA (GTP→ATP)

# Mitchell PMF target (chemiosmotic theory, Mitchell 1961)
DELTA_PSI_DEFAULT_V = -0.160          # membrane potential, matrix-negative
DELTA_PH_DEFAULT    = +1.0            # matrix more alkaline by ~1 unit


# ---------------------------------------------------------------------------
# Standalone helper functions
# ---------------------------------------------------------------------------

def proton_motive_force(delta_psi_v: float = DELTA_PSI_DEFAULT_V,
                        delta_ph: float = DELTA_PH_DEFAULT,
                        T: float = T_BODY) -> float:
    """Mitchell PMF in volts (positive convention).

        |Δμ_H+| / F  =  |ΔΨ|  +  (2.303 RT/F) · ΔpH

    With the standard mammalian numbers ΔΨ = -160 mV, ΔpH = 1:
        |PMF|  =  0.160 + 0.0615 × 1  ≈ 0.22 V (≈ 200 mV)

    Returns the absolute value so it reads as a positive driving force.
    """
    rt_f = R_GAS * T / F_FARADAY
    ln10_rt_f = 2.302_585 * rt_f
    pmf = abs(delta_psi_v) + ln10_rt_f * delta_ph
    return float(pmf)


def delta_g_per_proton(pmf_v: Optional[float] = None,
                       T: float = T_BODY) -> float:
    """ΔG per H⁺ traversing the membrane = F · PMF (in J/mol).

    For the default 200 mV PMF this is 21 kJ/mol per proton.
    """
    if pmf_v is None:
        pmf_v = proton_motive_force(T=T)
    return float(F_FARADAY * pmf_v)


def delta_g_redox(delta_e_v: float, n_electrons: int = 2) -> float:
    """ΔG = -n F ΔE (J/mol) for a redox half-cell drop ΔE.

    NADH → O₂ has ΔE = +1.14 V, n = 2 ⇒ ΔG = −220 kJ/mol per 2 e⁻.
    Sign convention: ΔG < 0 means the reaction is exergonic.
    """
    return float(-n_electrons * F_FARADAY * delta_e_v)


def c_ring_atp_ratio(c: int) -> float:
    """H⁺ per ATP = c / 3 (one full 360° rotation = 3 ATP).

    c = 8  (mammal):  2.667 H⁺/ATP
    c = 10 (yeast):   3.333 H⁺/ATP
    c = 14 (plant):   4.667 H⁺/ATP
    """
    if c <= 0:
        raise ValueError("c-ring size must be positive")
    return float(c) / float(ATP_PER_REV)


def atp_per_nadh(h_per_atp: float = c_ring_atp_ratio(C_RING_MAMMAL),
                 h_pumped_per_2e: int = (H_PER_2E_COMPLEX_I +
                                          H_PER_2E_COMPLEX_III +
                                          H_PER_2E_COMPLEX_IV)) -> float:
    """ATP yielded per NADH oxidised = (H⁺ pumped) / (H⁺ per ATP).

    Mammal default: 10 H⁺ pumped (4+4+2) / 2.667 H⁺ per ATP ≈ 3.75 ATP
    With the textbook leak-corrected P/O = 2.5 the effective yield is 2.5.
    The 10/4 = 2.5 textbook number assumes 4 H⁺/ATP including the
    phosphate antiporter cost; we expose both.
    """
    return float(h_pumped_per_2e) / float(h_per_atp)


def atp_per_fadh2(h_per_atp: float = c_ring_atp_ratio(C_RING_MAMMAL),
                  h_pumped_per_2e: int = (H_PER_2E_COMPLEX_III +
                                           H_PER_2E_COMPLEX_IV)) -> float:
    """ATP yielded per FADH₂ = (H⁺ pumped via III+IV only) / (H⁺ per ATP).

    Mammal default: 6 H⁺ / 2.667 ≈ 2.25 ATP (textbook leak-corrected: 1.5).
    """
    return float(h_pumped_per_2e) / float(h_per_atp)


def glucose_to_atp(p_o_nadh: float = 2.5,
                   p_o_fadh2: float = 1.5,
                   nadh_per_glucose: int = NADH_PER_GLUCOSE,
                   fadh2_per_glucose: int = FADH2_PER_GLUCOSE,
                   substrate_level_atp: int = SUBSTRATE_LEVEL_ATP) -> float:
    """Total ATP per glucose = SLP + NADH×P/O_NADH + FADH₂×P/O_FADH₂.

    Default yields the textbook ~30 ATP per glucose.
    """
    return float(substrate_level_atp +
                 nadh_per_glucose * p_o_nadh +
                 fadh2_per_glucose * p_o_fadh2)


# ---------------------------------------------------------------------------
# 1. MitoGeometry — inner membrane, cristae, complex placement
# ---------------------------------------------------------------------------

@dataclass
class ETCComplex:
    """One inner-membrane ETC complex.

    Attributes
    ----------
    name :
        "I", "II", "III", or "IV".
    e_in_v, e_out_v :
        Redox potentials of the donor and acceptor (V at pH 7).
    h_pumped_per_2e :
        Stoichiometric protons translocated per 2-electron passage.
    position_x :
        Lateral position along the inner membrane (arbitrary cartoon).
    """

    name: str
    e_in_v: float
    e_out_v: float
    h_pumped_per_2e: int
    position_x: float

    def delta_e(self) -> float:
        """Redox drop ΔE = E_out − E_in (positive for downhill)."""
        return float(self.e_out_v - self.e_in_v)

    def delta_g_per_2e(self) -> float:
        """ΔG (J/mol) per 2 e⁻ traversal = -nF(E_out − E_in)."""
        return delta_g_redox(self.delta_e(), n_electrons=2)


@dataclass
class MitoGeometry:
    """Inner membrane + cristae + ETC complexes + ATP synthase layout.

    Attributes
    ----------
    delta_psi_v :
        Membrane potential ΔΨ in volts (matrix-negative, default -160 mV).
    delta_ph :
        ΔpH = pH_matrix − pH_IMS (default +1, matrix more alkaline).
    cristae_fold_factor :
        Surface-area amplification from cristae folding (default 5).
        Inner membrane area ≈ cristae_fold_factor × outer membrane area.
    c_ring_size :
        Number of c-subunits in the F0 rotor (default 8 = mammalian).
    temperature :
        Bath temperature (K), default 310.15 K (37°C).
    """

    delta_psi_v: float = DELTA_PSI_DEFAULT_V
    delta_ph: float = DELTA_PH_DEFAULT
    cristae_fold_factor: float = 5.0
    c_ring_size: int = C_RING_MAMMAL
    temperature: float = T_BODY

    complexes: List[ETCComplex] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.complexes:
            # canonical cartoon: 4 complexes laid out left→right
            self.complexes = [
                ETCComplex("I",   E_NADH_NAD,  E_Q_QH2, H_PER_2E_COMPLEX_I,   0.15),
                ETCComplex("II",  E_FADH2_FAD, E_Q_QH2, H_PER_2E_COMPLEX_II,  0.30),
                ETCComplex("III", E_Q_QH2,     E_CYTC,  H_PER_2E_COMPLEX_III, 0.55),
                ETCComplex("IV",  E_CYTC,      E_O2_H2O, H_PER_2E_COMPLEX_IV, 0.80),
            ]

    # ------------------------------------------------------------------
    # Bioenergetic accessors
    # ------------------------------------------------------------------
    def proton_motive_force(self) -> float:
        """|Δμ_H+|/F = |ΔΨ| + (2.303 RT/F) · ΔpH   (V)."""
        return proton_motive_force(self.delta_psi_v, self.delta_ph,
                                    self.temperature)

    def delta_g_per_proton(self) -> float:
        """ΔG per H⁺ traversing F0 = F · PMF  (J/mol)."""
        return delta_g_per_proton(self.proton_motive_force(),
                                   self.temperature)

    def h_per_atp(self) -> float:
        """H⁺ per ATP = c-ring size / 3 catalytic β-sites."""
        return c_ring_atp_ratio(self.c_ring_size)

    def total_h_pumped_nadh(self) -> int:
        """Sum of H⁺/2e⁻ across the NADH pathway (I + III + IV)."""
        return int(H_PER_2E_COMPLEX_I + H_PER_2E_COMPLEX_III +
                    H_PER_2E_COMPLEX_IV)

    def total_h_pumped_fadh2(self) -> int:
        """Sum of H⁺/2e⁻ across the FADH₂ pathway (III + IV only)."""
        return int(H_PER_2E_COMPLEX_III + H_PER_2E_COMPLEX_IV)

    # ------------------------------------------------------------------
    # Spatial profiles for the renderer
    # ------------------------------------------------------------------
    def membrane_profile(self, x: NDArray[np.float64]
                          ) -> NDArray[np.float64]:
        """Inner membrane y-profile with cristae folds along x in [0, 1]."""
        x = np.asarray(x, dtype=float)
        # Smooth folded membrane: low-amplitude cristae corrugation
        n_folds = max(int(round(self.cristae_fold_factor)), 1)
        return 0.05 * np.sin(2.0 * np.pi * n_folds * x)

    def proton_gradient(self, y: NDArray[np.float64]
                         ) -> NDArray[np.float64]:
        """Proton concentration profile across the membrane (cartoon).

        y in [-1, 1], membrane at y = 0.  Returns relative [H⁺] with
        intermembrane space (y > 0) high, matrix (y < 0) low.
        """
        y = np.asarray(y, dtype=float)
        # 10× concentration ratio = 1 pH unit
        ratio = 10.0 ** self.delta_ph
        return np.where(y > 0.0, 1.0, 1.0 / ratio)

    # ------------------------------------------------------------------
    def summary(self) -> Dict[str, float]:
        return {
            "delta_psi_v":     float(self.delta_psi_v),
            "delta_ph":        float(self.delta_ph),
            "pmf_v":           float(self.proton_motive_force()),
            "dG_per_H_kjmol":  float(self.delta_g_per_proton() / 1000.0),
            "c_ring":          int(self.c_ring_size),
            "h_per_atp":       float(self.h_per_atp()),
            "h_pumped_nadh":   int(self.total_h_pumped_nadh()),
            "h_pumped_fadh2":  int(self.total_h_pumped_fadh2()),
            "n_complexes":     len(self.complexes),
        }


# ---------------------------------------------------------------------------
# 2. MitoSimulator — ETC + PMF + ATP synthase + glucose accounting
# ---------------------------------------------------------------------------

@dataclass
class MitoSimulator:
    """End-to-end mitochondrial bioenergetics simulator.

    Brings together
      * Electron transport chain (NADH → Q → cyt c → O₂)
      * Proton motive force ledger (Δμ_H+ ≈ 200 mV)
      * F0F1 ATP synthase rotor (c rotations, 3 ATP per 360°)
      * P/O ratio and glucose → ATP accounting
    """

    geometry: MitoGeometry = field(default_factory=MitoGeometry)
    p_o_nadh:  float = 2.5      # textbook leak-corrected
    p_o_fadh2: float = 1.5      # textbook leak-corrected

    # ------------------------------------------------------------------
    # ETC observables
    # ------------------------------------------------------------------
    def etc_redox_ladder(self) -> List[Tuple[str, float]]:
        """Energy-ordered redox levels traversed by the chain.

        Returns [(label, E°' in V), ...] from NADH down to ½ O₂/H₂O.
        """
        return [
            ("NADH/NAD⁺",   E_NADH_NAD),
            ("FADH₂/FAD",   E_FADH2_FAD),
            ("Q / QH₂",     E_Q_QH2),
            ("cyt c³⁺/²⁺",  E_CYTC),
            ("½ O₂ / H₂O",  E_O2_H2O),
        ]

    def etc_total_delta_e_nadh(self) -> float:
        """Total redox drop NADH → O₂ (V) = E_O2 − E_NADH ≈ 1.14 V."""
        return float(E_O2_H2O - E_NADH_NAD)

    def etc_total_delta_g_nadh_kjmol(self) -> float:
        """ΔG (kJ/mol) for NADH + ½O₂ + H⁺ → NAD⁺ + H₂O   (≈ -220 kJ/mol)."""
        return float(delta_g_redox(self.etc_total_delta_e_nadh(),
                                    n_electrons=2) / 1000.0)

    def etc_h_pumped_per_nadh(self) -> int:
        """H⁺ translocated per NADH oxidised (I + III + IV)."""
        return self.geometry.total_h_pumped_nadh()

    def etc_h_pumped_per_fadh2(self) -> int:
        """H⁺ translocated per FADH₂ oxidised (III + IV only)."""
        return self.geometry.total_h_pumped_fadh2()

    # ------------------------------------------------------------------
    # PMF observables
    # ------------------------------------------------------------------
    def pmf_v(self) -> float:
        """Proton motive force in volts (≈ 0.20 — 0.22 V)."""
        return self.geometry.proton_motive_force()

    def pmf_mv(self) -> float:
        """PMF in millivolts (Mitchell quotes ≈ 200 mV)."""
        return 1000.0 * self.pmf_v()

    def delta_g_per_proton_kjmol(self) -> float:
        """ΔG per H⁺ traversing F0 in kJ/mol (≈ 21 kJ/mol)."""
        return self.geometry.delta_g_per_proton() / 1000.0

    # ------------------------------------------------------------------
    # F0F1 ATP synthase observables
    # ------------------------------------------------------------------
    def h_per_atp(self) -> float:
        """H⁺ per ATP = c-ring size / 3 catalytic β-sites."""
        return self.geometry.h_per_atp()

    def atp_per_revolution(self) -> int:
        """ATP synthesized per 360° rotor rotation (= 3, fixed by α₃β₃ head)."""
        return ATP_PER_REV

    def atp_per_c_protons(self, n_protons: float) -> float:
        """ATP synthesised after n_protons have flowed through F0.

        Continuous form: ATP = n_protons / (c/3) = n_protons × 3 / c.
        Discrete bursts of 3 happen each full revolution; this is the
        continuous-flow version useful for the rate plot.
        """
        return float(n_protons) * float(ATP_PER_REV) / float(
            self.geometry.c_ring_size)

    def revolutions_per_atp(self) -> float:
        """Fractional revolution per ATP = 1/3 (since 3 ATP/rev)."""
        return 1.0 / float(ATP_PER_REV)

    def torque_required_nm(self) -> float:
        """Rotor torque (pN·nm) required to drive F0F1 against PMF.

        Each 360° revolution moves c protons across the membrane, doing
        work W = c · q_e · PMF (per single rotor molecule).
        Torque τ = W / (2π rad) = c · q_e · PMF / 2π   (J/rad = N·m)
        For c=8, PMF=0.22 V: τ ≈ 45 pN·nm — matches single-molecule data
        (Yasuda 2001: ~40 pN·nm; Toyabe 2011 single-molecule torque).

        Returns torque in pN·nm (the unit conventional in the F1 literature).
        """
        c = float(self.geometry.c_ring_size)
        pmf = self.pmf_v()
        q_e = 1.602_176_634e-19         # elementary charge (C)
        # Energy per revolution per single rotor molecule (J)
        e_per_rev = c * q_e * pmf
        # Torque (N·m) = E / (2π rad)
        torque_nm = e_per_rev / (2.0 * np.pi)
        # Convert N·m → pN·nm:  1 N·m = 1e12 pN × 1e9 nm = 1e21 pN·nm
        torque_pnm = torque_nm * 1.0e21
        return float(torque_pnm)

    # ------------------------------------------------------------------
    # P/O and glucose accounting
    # ------------------------------------------------------------------
    def p_o_ratio_nadh(self) -> float:
        """P/O for the NADH pathway (default 2.5 ATP per NADH)."""
        return float(self.p_o_nadh)

    def p_o_ratio_fadh2(self) -> float:
        """P/O for the FADH₂ pathway (default 1.5 ATP per FADH₂)."""
        return float(self.p_o_fadh2)

    def atp_per_glucose(self) -> float:
        """Total ATP per glucose ≈ 30 (mammal) or up to 32 (idealised)."""
        return glucose_to_atp(p_o_nadh=self.p_o_nadh,
                              p_o_fadh2=self.p_o_fadh2)

    def atp_per_glucose_breakdown(self) -> Dict[str, float]:
        """Itemised ATP yield per glucose."""
        nadh_atp = NADH_PER_GLUCOSE * self.p_o_nadh
        fadh2_atp = FADH2_PER_GLUCOSE * self.p_o_fadh2
        return {
            "substrate_level":    float(SUBSTRATE_LEVEL_ATP),
            "from_nadh":          float(nadh_atp),
            "from_fadh2":         float(fadh2_atp),
            "total":              float(SUBSTRATE_LEVEL_ATP + nadh_atp +
                                        fadh2_atp),
            "n_nadh":             int(NADH_PER_GLUCOSE),
            "n_fadh2":            int(FADH2_PER_GLUCOSE),
        }

    # ------------------------------------------------------------------
    # Time-series simulation of ETC + PMF + ATP synthesis
    # ------------------------------------------------------------------
    def simulate(self,
                 n_nadh: float = 100.0,
                 n_fadh2: float = 0.0,
                 n_steps: int = 200,
                 dt: float = 1.0) -> Dict[str, NDArray]:
        """Toy first-order kinetics: NADH consumption → H⁺ build-up → ATP.

        Single first-order rate constant k_etc per NADH equivalent;
        proton flux through F0 set by F1 catalytic rate k_atp.
        Not meant to be quantitatively predictive — only for visualising
        coupled kinetics in the renderer.
        """
        n_steps = int(max(n_steps, 1))
        dt = float(dt)

        k_etc = 0.05         # NADH oxidation rate per arbitrary timestep
        k_atp = 0.10         # ATP synthase turnover rate

        nadh = np.zeros(n_steps + 1, dtype=float)
        fadh2_arr = np.zeros(n_steps + 1, dtype=float)
        h_pool = np.zeros(n_steps + 1, dtype=float)
        atp = np.zeros(n_steps + 1, dtype=float)
        t = np.arange(n_steps + 1, dtype=float) * dt

        nadh[0] = float(n_nadh)
        fadh2_arr[0] = float(n_fadh2)
        h_per_atp = self.h_per_atp()
        h_per_nadh = self.etc_h_pumped_per_nadh()
        h_per_fadh2 = self.etc_h_pumped_per_fadh2()

        for i in range(1, n_steps + 1):
            # First-order consumption of NADH and FADH₂
            d_nadh = -k_etc * nadh[i - 1] * dt
            d_fadh2 = -k_etc * fadh2_arr[i - 1] * dt
            nadh[i] = max(nadh[i - 1] + d_nadh, 0.0)
            fadh2_arr[i] = max(fadh2_arr[i - 1] + d_fadh2, 0.0)

            # H⁺ pumped into the pool by ETC
            h_in = -d_nadh * h_per_nadh + (-d_fadh2) * h_per_fadh2

            # H⁺ flux through F0 saturates at k_atp · H_pool
            flux = min(k_atp * h_pool[i - 1] * dt, h_pool[i - 1])
            atp[i] = atp[i - 1] + flux / h_per_atp
            h_pool[i] = max(h_pool[i - 1] + h_in - flux, 0.0)

        return {
            "t":       t,
            "nadh":    nadh,
            "fadh2":   fadh2_arr,
            "h_pool":  h_pool,
            "atp":     atp,
        }

    # ------------------------------------------------------------------
    def summary(self) -> Dict[str, object]:
        gluc = self.atp_per_glucose_breakdown()
        return {
            "pmf_mv":                 float(self.pmf_mv()),
            "dG_per_H_kjmol":         float(
                self.delta_g_per_proton_kjmol()),
            "c_ring":                 int(self.geometry.c_ring_size),
            "h_per_atp":              float(self.h_per_atp()),
            "atp_per_revolution":     int(self.atp_per_revolution()),
            "torque_pNnm":            float(self.torque_required_nm()),
            "etc_total_delta_e_v":    float(self.etc_total_delta_e_nadh()),
            "etc_total_dG_kjmol":     float(
                self.etc_total_delta_g_nadh_kjmol()),
            "h_pumped_per_nadh":      int(self.etc_h_pumped_per_nadh()),
            "h_pumped_per_fadh2":     int(self.etc_h_pumped_per_fadh2()),
            "p_o_nadh":               float(self.p_o_nadh),
            "p_o_fadh2":              float(self.p_o_fadh2),
            "atp_per_glucose":        float(self.atp_per_glucose()),
            "atp_per_glucose_break":  gluc,
        }


# ---------------------------------------------------------------------------
# 3. Convenience helpers used by the renderer
# ---------------------------------------------------------------------------

def etc_complexes_default() -> List[ETCComplex]:
    """Cartoon list of 4 ETC complexes laid out along the inner membrane."""
    return MitoGeometry().complexes


def example_c_rings() -> List[Dict[str, float]]:
    """Catalogue of c-ring sizes across organisms (for the H⁺/ATP table)."""
    rows = [
        ("E. coli (bacterial)",   C_RING_BACTERIA),
        ("S. cerevisiae (yeast)", C_RING_YEAST),
        ("Bovine / human",        C_RING_MAMMAL),
        ("Spinach chloroplast",   C_RING_PLANT),
    ]
    out: List[Dict[str, float]] = []
    for name, c in rows:
        out.append({
            "name": name,
            "c": int(c),
            "h_per_atp": float(c_ring_atp_ratio(c)),
            "atp_per_rev": int(ATP_PER_REV),
        })
    return out
