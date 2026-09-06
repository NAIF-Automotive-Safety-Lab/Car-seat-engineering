"""Paired neural action-potential GEOMETRY + SIMULATOR for the substrate framework.

Models the Hodgkin-Huxley action potential (1952 Nobel work) under the B3 /
stiff-medium ontology, where ion channels are interpreted as locally
saturated cells of the substrate that gate substrate-strain transmission
along the axon membrane.

GEOMETRY
--------
``NeuralGeometry`` represents the axon as a 1D cable bounded by a
phospholipid bilayer membrane (~5 nm thick).  Three substrate "channel"
species pierce the membrane:

    * voltage-gated Na+ channels (fast activation, slow inactivation)
    * voltage-gated K+ channels (delayed rectifier)
    * leak channels (passive Cl-/K+ conductance)

For each species the Nernst (equilibrium) potential is computed from the
intracellular/extracellular concentration ratio at body temperature:

    E_ion = (RT/zF) ln([ion]_out / [ion]_in)

Standard mammalian values (Hodgkin-Huxley squid axon variant):

    E_Na = +50  mV       (Na+ rushes in to depolarize)
    E_K  = -77  mV       (K+ rushes out to repolarize)
    E_L  = -54.4 mV      (leak sets resting near threshold)
    V_rest    = -65  mV
    V_threshold = -55 mV
    V_peak    ≈ +40  mV

SIMULATOR
---------
``NeuralSimulator`` integrates the Hodgkin-Huxley equations:

    C_m dV/dt = -[ I_Na + I_K + I_L ] + I_inj
    I_Na = g_Na_max * m^3 * h * (V - E_Na)
    I_K  = g_K_max  * n^4       * (V - E_K)
    I_L  = g_L      *           (V - E_L)

with the gating variables m, h, n following first-order kinetics

    dx/dt = alpha_x(V) (1 - x) - beta_x(V) x

using the canonical Hodgkin-Huxley 1952 alpha/beta rate functions
(shifted so V_rest = -65 mV in modern convention rather than 0 mV).

Substrate interpretation:
    * V(t) = strain potential across membrane = saturation differential
      (intracellular vs extracellular substrate density)
    * channel opening = cell de-saturation gate triggered by local
      strain crossing the threshold V_th = -55 mV
    * refractory period = h-gate slow re-arming = substrate cell needs
      time to re-saturate after firing
    * conduction velocity scales with (myelin sheath insulation)^(1/2)
      because saltatory conduction skips strain-leak loss between nodes
      of Ranvier:
          v_unmyelinated  ~  1   m/s   (continuous)
          v_myelinated    ~ 100  m/s   (saltatory)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference values (Hodgkin-Huxley 1952, modern V_rest = -65 mV convention)
# ---------------------------------------------------------------------------

# Membrane capacitance per unit area  [microF / cm^2]
C_M_DEFAULT: float = 1.0

# Maximum conductances per unit area  [mS / cm^2]
G_NA_MAX_DEFAULT: float = 120.0
G_K_MAX_DEFAULT:  float = 36.0
G_L_DEFAULT:      float = 0.3

# Reversal (Nernst / leak) potentials  [mV]
E_NA_DEFAULT: float = +50.0
E_K_DEFAULT:  float = -77.0
E_L_DEFAULT:  float = -54.4

# Voltage milestones  [mV]
V_REST:       float = -65.0
V_THRESHOLD:  float = -55.0
V_PEAK_TYP:   float = +40.0

# Refractory bookkeeping  [ms]
ABSOLUTE_REFRACTORY_MS: float = 1.0
RELATIVE_REFRACTORY_MS: float = 4.0

# Conduction velocities  [m / s]
V_UNMYELINATED_M_PER_S: float = 1.0
V_MYELINATED_M_PER_S:   float = 100.0

# Mammalian physiological ion concentrations  [mM]  (intra / extra)
NA_IN_MM:  float = 15.0
NA_OUT_MM: float = 145.0
K_IN_MM:   float = 140.0
K_OUT_MM:  float = 4.0
CL_IN_MM:  float = 10.0
CL_OUT_MM: float = 110.0

# Physical constants for the Nernst equation
R_GAS_CONST_J_PER_MOL_K: float = 8.314_462_618
F_FARADAY_C_PER_MOL:     float = 96_485.332_12
T_BODY_K:                float = 310.15  # 37 deg C

# Substrate anchor (B3 framework)
LAMBDA_QCD_MEV: float = 200.0


# ---------------------------------------------------------------------------
# Hodgkin-Huxley gating-variable rate functions
# ---------------------------------------------------------------------------

def _hh_alpha_m(V: np.ndarray | float) -> np.ndarray | float:
    """Na+ activation forward rate alpha_m(V).  V in mV.  Result in 1/ms."""
    V = np.asarray(V, dtype=float)
    # Shifted Hodgkin-Huxley 1952 form with V_rest = -65 mV
    dv = V + 40.0
    out = np.where(
        np.abs(dv) > 1.0e-6,
        0.1 * dv / (1.0 - np.exp(-dv / 10.0)),
        1.0,  # L'Hopital limit
    )
    return out


def _hh_beta_m(V: np.ndarray | float) -> np.ndarray | float:
    V = np.asarray(V, dtype=float)
    return 4.0 * np.exp(-(V + 65.0) / 18.0)


def _hh_alpha_h(V: np.ndarray | float) -> np.ndarray | float:
    V = np.asarray(V, dtype=float)
    return 0.07 * np.exp(-(V + 65.0) / 20.0)


def _hh_beta_h(V: np.ndarray | float) -> np.ndarray | float:
    V = np.asarray(V, dtype=float)
    return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))


def _hh_alpha_n(V: np.ndarray | float) -> np.ndarray | float:
    V = np.asarray(V, dtype=float)
    dv = V + 55.0
    out = np.where(
        np.abs(dv) > 1.0e-6,
        0.01 * dv / (1.0 - np.exp(-dv / 10.0)),
        0.1,  # L'Hopital limit
    )
    return out


def _hh_beta_n(V: np.ndarray | float) -> np.ndarray | float:
    V = np.asarray(V, dtype=float)
    return 0.125 * np.exp(-(V + 65.0) / 80.0)


def _x_inf(alpha, beta):
    return alpha / (alpha + beta)


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class IonChannel:
    """Single ion-channel species (Na+, K+, leak)."""
    name: str
    g_max: float          # max conductance [mS/cm^2]
    E_rev: float          # reversal potential [mV]
    voltage_gated: bool   # true for Na+/K+, false for leak
    valence: int = 1      # +1 for Na+/K+, -1 for Cl-, 0 for leak


@dataclass
class NeuralGeometry:
    """Axon membrane + ion-channel inventory.

    The axon is a 1D cable of length ``axon_length_m`` and radius
    ``axon_radius_m`` bounded by a phospholipid bilayer of thickness
    ``membrane_thickness_m``.  The membrane carries three channel
    populations whose densities set the per-area conductances.

    Reversal potentials are stored as the Hodgkin-Huxley 1952 calibrated
    values (E_Na = +50, E_K = -77, E_L = -54.4 mV) with which the canonical
    alpha/beta gating-rate functions were determined; these slightly
    differ from textbook Nernst values for mammalian concentrations
    (E_Na ~ +61, E_K ~ -95) which are exposed separately as
    ``E_Na_nernst`` / ``E_K_nernst``.
    """

    axon_length_m:        float = 1.0e-3        # 1 mm
    axon_radius_m:        float = 5.0e-6        # 5 microns (typical)
    membrane_thickness_m: float = 5.0e-9        # 5 nm bilayer
    myelinated:           bool  = False

    # Ion concentrations (mM)
    na_in_mM:  float = NA_IN_MM
    na_out_mM: float = NA_OUT_MM
    k_in_mM:   float = K_IN_MM
    k_out_mM:  float = K_OUT_MM
    cl_in_mM:  float = CL_IN_MM
    cl_out_mM: float = CL_OUT_MM

    # Channel densities (per cm^2 of membrane) and per-area conductance
    g_Na_max: float = G_NA_MAX_DEFAULT
    g_K_max:  float = G_K_MAX_DEFAULT
    g_L:      float = G_L_DEFAULT

    # HH 1952 reversal potentials (used by simulator; calibrated to the
    # canonical alpha/beta rate functions). Stored as fields so they can
    # be overridden if a user wants to switch to Nernst-derived values.
    E_Na_rev: float = E_NA_DEFAULT
    E_K_rev:  float = E_K_DEFAULT
    E_L_rev:  float = E_L_DEFAULT

    # Membrane capacitance per area [microF/cm^2]
    C_m: float = C_M_DEFAULT

    # ------------------------------------------------------------------
    # Channel inventory
    # ------------------------------------------------------------------

    def channels(self) -> Tuple[IonChannel, IonChannel, IonChannel]:
        return (
            IonChannel("Na+",  self.g_Na_max, self.E_Na(),  True,  +1),
            IonChannel("K+",   self.g_K_max,  self.E_K(),   True,  +1),
            IonChannel("leak", self.g_L,      self.E_L(),   False,  0),
        )

    # ------------------------------------------------------------------
    # Nernst (equilibrium) potentials
    # ------------------------------------------------------------------

    @staticmethod
    def nernst_potential_mV(c_in_mM: float, c_out_mM: float,
                            valence: int = 1, T_K: float = T_BODY_K) -> float:
        """Nernst equation: E = (RT/zF) ln(c_out / c_in)  in mV."""
        return (
            (R_GAS_CONST_J_PER_MOL_K * T_K)
            / (valence * F_FARADAY_C_PER_MOL)
            * math.log(c_out_mM / c_in_mM)
            * 1000.0
        )

    def E_Na(self) -> float:
        """HH 1952 reversal potential for Na+ (calibrated)."""
        return self.E_Na_rev

    def E_K(self) -> float:
        """HH 1952 reversal potential for K+ (calibrated)."""
        return self.E_K_rev

    def E_L(self) -> float:
        """HH 1952 leak reversal (calibrated, phenomenological)."""
        return self.E_L_rev

    # Nernst-derived counterparts (textbook physiological values)
    def E_Na_nernst(self) -> float:
        return self.nernst_potential_mV(self.na_in_mM, self.na_out_mM, +1)

    def E_K_nernst(self) -> float:
        return self.nernst_potential_mV(self.k_in_mM, self.k_out_mM, +1)

    def E_Cl(self) -> float:
        return self.nernst_potential_mV(self.cl_in_mM, self.cl_out_mM, -1)

    # ------------------------------------------------------------------
    # Conduction
    # ------------------------------------------------------------------

    def conduction_velocity_m_per_s(self) -> float:
        """Saltatory (myelin) conduction is ~100x faster than continuous."""
        return V_MYELINATED_M_PER_S if self.myelinated else V_UNMYELINATED_M_PER_S

    def membrane_area_cm2(self) -> float:
        return 2.0 * math.pi * self.axon_radius_m * self.axon_length_m * 1.0e4  # m^2 -> cm^2

    def membrane_capacitance_uF(self) -> float:
        """Total membrane capacitance from per-area C_m."""
        return self.C_m * self.membrane_area_cm2()

    # ------------------------------------------------------------------
    # Substrate signature
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "axon = 1D substrate cable; ion channel = local saturation gate; "
                "AP = strain pulse propagating along membrane"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "channel_species": 3,
            "voltage_gated": 2,
            "myelinated": self.myelinated,
            "v_conduction_m_per_s": self.conduction_velocity_m_per_s(),
        }


# ---------------------------------------------------------------------------
# SIMULATOR: Hodgkin-Huxley action-potential integrator
# ---------------------------------------------------------------------------

@dataclass
class NeuralSimulator:
    """Hodgkin-Huxley equations integrated by explicit Euler.

    Solves the coupled ODE system

        C_m dV/dt = -[I_Na + I_K + I_L] + I_inj(t)
        dm/dt    = alpha_m(V)(1 - m) - beta_m(V) m
        dh/dt    = alpha_h(V)(1 - h) - beta_h(V) h
        dn/dt    = alpha_n(V)(1 - n) - beta_n(V) n

    Output arrays: t [ms], V [mV], m, h, n, I_Na, I_K, I_L  [microA/cm^2].
    """

    geometry: NeuralGeometry = field(default_factory=NeuralGeometry)
    dt_ms: float = 0.01
    duration_ms: float = 50.0

    # ------------------------------------------------------------------
    # Steady-state gating values at a given V
    # ------------------------------------------------------------------

    @staticmethod
    def steady_state_gating(V: float) -> dict:
        return {
            "m_inf": float(_x_inf(_hh_alpha_m(V), _hh_beta_m(V))),
            "h_inf": float(_x_inf(_hh_alpha_h(V), _hh_beta_h(V))),
            "n_inf": float(_x_inf(_hh_alpha_n(V), _hh_beta_n(V))),
        }

    # ------------------------------------------------------------------
    # Resting state at V_rest
    # ------------------------------------------------------------------

    def resting_state(self) -> Tuple[float, float, float, float]:
        """Return (V_rest, m0, h0, n0) at resting potential."""
        V0 = V_REST
        ss = self.steady_state_gating(V0)
        return V0, ss["m_inf"], ss["h_inf"], ss["n_inf"]

    # ------------------------------------------------------------------
    # Single-step ODE update
    # ------------------------------------------------------------------

    def _ionic_currents(self, V: float, m: float, h: float, n: float) -> Tuple[float, float, float]:
        g = self.geometry
        I_Na = g.g_Na_max * (m ** 3) * h * (V - g.E_Na())
        I_K  = g.g_K_max  * (n ** 4)     * (V - g.E_K())
        I_L  = g.g_L                     * (V - g.E_L())
        return I_Na, I_K, I_L

    # ------------------------------------------------------------------
    # Full simulation -- explicit Euler with adaptive sub-step on V spike
    # ------------------------------------------------------------------

    def simulate(self, I_inj_uA_per_cm2=10.0,
                 stim_start_ms: float = 5.0,
                 stim_duration_ms: float = 0.5,
                 second_stim_ms: float | None = None,
                 second_stim_duration_ms: float = 0.5) -> dict:
        """Simulate the action potential; returns arrays vs t.

        Parameters
        ----------
        I_inj_uA_per_cm2 : float
            Magnitude of the injected current pulse during the stimulus.
        stim_start_ms : float
            Time at which the first stimulus begins.
        stim_duration_ms : float
            Duration of the first stimulus pulse.
        second_stim_ms : float, optional
            Time at which an optional second stimulus begins (used for
            measuring the refractory period).  If ``None``, no second
            stimulus is delivered.
        """
        n_steps = int(round(self.duration_ms / self.dt_ms)) + 1
        t  = np.linspace(0.0, self.duration_ms, n_steps)
        V  = np.zeros(n_steps)
        m  = np.zeros(n_steps)
        h  = np.zeros(n_steps)
        n_ = np.zeros(n_steps)
        INa = np.zeros(n_steps)
        IK  = np.zeros(n_steps)
        IL  = np.zeros(n_steps)
        I_inj_arr = np.zeros(n_steps)

        V[0], m[0], h[0], n_[0] = self.resting_state()

        C_m = self.geometry.C_m

        for i in range(1, n_steps):
            t_now = t[i - 1]
            in_stim_1 = stim_start_ms <= t_now < stim_start_ms + stim_duration_ms
            in_stim_2 = (
                second_stim_ms is not None
                and second_stim_ms <= t_now < second_stim_ms + second_stim_duration_ms
            )
            I_inj = float(I_inj_uA_per_cm2) if (in_stim_1 or in_stim_2) else 0.0
            I_inj_arr[i - 1] = I_inj

            I_Na_now, I_K_now, I_L_now = self._ionic_currents(V[i - 1], m[i - 1], h[i - 1], n_[i - 1])
            INa[i - 1], IK[i - 1], IL[i - 1] = I_Na_now, I_K_now, I_L_now

            dV = (-(I_Na_now + I_K_now + I_L_now) + I_inj) / C_m
            V_new = V[i - 1] + self.dt_ms * dV

            am, bm = float(_hh_alpha_m(V[i - 1])), float(_hh_beta_m(V[i - 1]))
            ah, bh = float(_hh_alpha_h(V[i - 1])), float(_hh_beta_h(V[i - 1]))
            an, bn = float(_hh_alpha_n(V[i - 1])), float(_hh_beta_n(V[i - 1]))

            m_new = m[i - 1] + self.dt_ms * (am * (1.0 - m[i - 1]) - bm * m[i - 1])
            h_new = h[i - 1] + self.dt_ms * (ah * (1.0 - h[i - 1]) - bh * h[i - 1])
            n_new = n_[i - 1] + self.dt_ms * (an * (1.0 - n_[i - 1]) - bn * n_[i - 1])

            # Clamp gating variables to [0, 1] for numerical safety
            V[i]  = V_new
            m[i]  = max(0.0, min(1.0, m_new))
            h[i]  = max(0.0, min(1.0, h_new))
            n_[i] = max(0.0, min(1.0, n_new))

        # Fill last currents
        INa[-1], IK[-1], IL[-1] = self._ionic_currents(V[-1], m[-1], h[-1], n_[-1])

        return {
            "t_ms":  t,
            "V_mV":  V,
            "m":     m,
            "h":     h,
            "n":     n_,
            "I_Na":  INa,
            "I_K":   IK,
            "I_L":   IL,
            "I_inj": I_inj_arr,
        }

    # ------------------------------------------------------------------
    # Action-potential analysis
    # ------------------------------------------------------------------

    @staticmethod
    def peak_voltage_mV(result: dict) -> float:
        return float(np.max(result["V_mV"]))

    @staticmethod
    def crossed_threshold(result: dict, V_th: float = V_THRESHOLD) -> bool:
        return bool(np.max(result["V_mV"]) >= V_th)

    @staticmethod
    def fired_action_potential(result: dict, V_peak_min: float = 0.0) -> bool:
        """An AP requires the spike to overshoot 0 mV."""
        return bool(np.max(result["V_mV"]) >= V_peak_min)

    @staticmethod
    def spike_times_ms(result: dict, V_threshold: float = 0.0) -> np.ndarray:
        """Return times at which V crosses V_threshold from below."""
        V = result["V_mV"]
        t = result["t_ms"]
        crossings = (V[:-1] < V_threshold) & (V[1:] >= V_threshold)
        return t[1:][crossings]

    def absolute_refractory_ms(self) -> float:
        return ABSOLUTE_REFRACTORY_MS

    def relative_refractory_ms(self) -> float:
        return RELATIVE_REFRACTORY_MS

    def measure_refractory_ms(self,
                              stim_start_ms: float = 5.0,
                              I_inj_uA_per_cm2: float = 10.0,
                              gap_ms: tuple = (1.0, 2.0, 4.0, 6.0, 10.0),
                              stim_duration_ms: float = 0.5) -> dict:
        """Probe whether a second stimulus at +gap fires a 2nd AP.

        Returns a dict mapping each gap (ms) to True/False (fired or not)
        and the inferred minimum-fireable gap (refractory boundary).
        """
        results = {}
        for g in gap_ms:
            sim_out = self.simulate(
                I_inj_uA_per_cm2=I_inj_uA_per_cm2,
                stim_start_ms=stim_start_ms,
                stim_duration_ms=stim_duration_ms,
                second_stim_ms=stim_start_ms + g,
                second_stim_duration_ms=stim_duration_ms,
            )
            spikes = self.spike_times_ms(sim_out, V_threshold=0.0)
            results[g] = (len(spikes) >= 2)
        # Smallest gap that produced 2 spikes
        fired_gaps = [g for g, ok in results.items() if ok]
        min_gap = min(fired_gaps) if fired_gaps else float("inf")
        results["min_fireable_gap_ms"] = min_gap
        return results

    # ------------------------------------------------------------------
    # Conduction velocity (axon-level)
    # ------------------------------------------------------------------

    def conduction_velocity_m_per_s(self) -> float:
        return self.geometry.conduction_velocity_m_per_s()

    def conduction_speedup_myelin(self) -> float:
        return V_MYELINATED_M_PER_S / V_UNMYELINATED_M_PER_S

    # ------------------------------------------------------------------
    # Substrate signature
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "AP = substrate strain pulse; ion channel = local saturation gate; "
                "refractory period = h-gate slow re-arming "
                "(substrate cell needs time to re-saturate after firing)"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "v_conduction_unmyelinated_m_per_s": V_UNMYELINATED_M_PER_S,
            "v_conduction_myelinated_m_per_s":   V_MYELINATED_M_PER_S,
            "myelin_speedup": V_MYELINATED_M_PER_S / V_UNMYELINATED_M_PER_S,
        }


__all__ = [
    "IonChannel",
    "NeuralGeometry",
    "NeuralSimulator",
    "C_M_DEFAULT",
    "G_NA_MAX_DEFAULT",
    "G_K_MAX_DEFAULT",
    "G_L_DEFAULT",
    "E_NA_DEFAULT",
    "E_K_DEFAULT",
    "E_L_DEFAULT",
    "V_REST",
    "V_THRESHOLD",
    "V_PEAK_TYP",
    "ABSOLUTE_REFRACTORY_MS",
    "RELATIVE_REFRACTORY_MS",
    "V_UNMYELINATED_M_PER_S",
    "V_MYELINATED_M_PER_S",
    "NA_IN_MM", "NA_OUT_MM",
    "K_IN_MM",  "K_OUT_MM",
    "CL_IN_MM", "CL_OUT_MM",
    "R_GAS_CONST_J_PER_MOL_K",
    "F_FARADAY_C_PER_MOL",
    "T_BODY_K",
    "LAMBDA_QCD_MEV",
]
