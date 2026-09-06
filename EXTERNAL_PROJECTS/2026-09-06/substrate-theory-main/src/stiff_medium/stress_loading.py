"""Stress-loading vertex mechanism — §18.30 refined lepton mass.

Per user clarification (§18.30 refined): muon and tau are NOT separate
particle species. They are the SAME electron in stable excited states,
formed when collider-supplied energy is loaded onto the vertex.

The vertex absorbs discrete stress-quanta:
- Electron: ground state (0 stress quanta)
- Muon: 1 quantum of vertex stress
- Tau: 2 quanta
- ≥ 3: cannot close; immediate decay (forbidden)

This eliminates 2 distinct Dirac fields and 2 Yukawa couplings vs SM,
leaving 2 excitation energies (Δ_1, Δ_2) of ONE field.

Decay: vertex sheds stress quanta as 2 neutrinos + electron (V-A weak).
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class StressLoadedConfiguration:
    """An electron-like bound configuration with stress quanta loaded."""
    n_quanta: int = 0  # 0=electron, 1=muon, 2=tau
    base_mass: float = 0.000511  # GeV (electron mass)

    # Excitation energies (empirical; from measured masses)
    delta_1: float = 0.105148  # GeV → muon mass = base + delta_1
    delta_2: float = 1.776349  # GeV → tau mass = base + delta_2

    @property
    def mass(self) -> float:
        """Mass of stress-loaded configuration."""
        if self.n_quanta == 0:
            return self.base_mass
        elif self.n_quanta == 1:
            return self.base_mass + self.delta_1
        elif self.n_quanta == 2:
            return self.base_mass + self.delta_2
        else:
            raise ValueError(f"n_quanta = {self.n_quanta} not allowed (max 2)")

    @property
    def name(self) -> str:
        """Particle name based on n_quanta."""
        return ["electron", "muon", "tau"][self.n_quanta]

    @property
    def is_stable(self) -> bool:
        """Only n=0 (electron) is absolutely stable."""
        return self.n_quanta == 0

    def vertex_closure_check(self) -> bool:
        """Per §6.4: vertex closure caps stress quanta at 3.

        n=0,1,2 admit coherent bound states.
        n≥3 fails closure → no bound state (forbidden).
        """
        return self.n_quanta < 3

    def lifetime_factor(self) -> float:
        """Approximate lifetime ratio relative to muon lifetime.

        Decay rate Γ ∝ G_F² × m^5 (V-A weak).
        For muon at base level (n=1), τ_μ = 2.197 μs.
        """
        if self.n_quanta == 0:
            return float('inf')  # stable
        m_n = self.mass
        m_mu = self.base_mass + self.delta_1
        # Γ ∝ m^5; τ ∝ 1/Γ
        return (m_mu / m_n)**5

    def winding_number(self) -> int:
        """Möbius half-flux winding for this state.

        Maps to charge: q = -1 for all charged leptons (single winding).
        Stress quanta don't change winding; they add angular momentum
        at the vertex.
        """
        return 1


def stress_load_electron(electron, n_quanta_to_add: int) -> StressLoadedConfiguration:
    """Add stress quanta to an electron, producing μ or τ excited state."""
    if n_quanta_to_add < 0 or n_quanta_to_add > 2:
        raise ValueError(f"Cannot load {n_quanta_to_add} quanta (max 2)")
    return StressLoadedConfiguration(n_quanta=n_quanta_to_add)


def vertex_stress_quanta_decay_products(state: StressLoadedConfiguration) -> dict:
    """Predict decay products when vertex sheds stress quanta.

    For n=1 (muon): μ → e + ν_μ + ν̄_e (vertex sheds 1 quantum)
    For n=2 (tau): τ → ν_τ + W* (W* → various final states with branching)
    """
    if state.n_quanta == 0:
        return {"products": ["stable"]}
    elif state.n_quanta == 1:
        return {
            "products": ["e⁻", "ν̄_e", "ν_μ"],
            "lifetime_us": 2.1969811,
            "channel": "leptonic μ → e ν̄_e ν_μ",
        }
    elif state.n_quanta == 2:
        return {
            "products": ["multiple channels (e+2ν, μ+2ν, hadrons+ν)"],
            "lifetime_fs": 290.3,
            "branching": {
                "τ → e + ν̄_e + ν_τ": 0.1782,
                "τ → μ + ν̄_μ + ν_τ": 0.1739,
                "τ → hadrons + ν_τ": 0.6478,
            },
        }
    else:
        raise ValueError("n_quanta > 2 forbidden")


def koide_relation_check(states: list) -> dict:
    """Verify Koide's relation Q = (Σm)/(Σ√m)² = 2/3 for given states.

    Should give Q ≈ 2/3 to high precision for the lepton triplet."""
    masses = [s.mass for s in states]
    sum_m = sum(masses)
    sum_sqrt_m = sum(np.sqrt(m) for m in masses)
    Q = sum_m / sum_sqrt_m**2
    return {
        "masses_GeV": masses,
        "Q": Q,
        "Q_target": 2/3,
        "deviation_pct": abs(Q - 2/3) / (2/3) * 100,
    }


def all_three_lepton_states():
    """Return the three lepton states (electron, muon, tau)."""
    return [
        StressLoadedConfiguration(n_quanta=0),
        StressLoadedConfiguration(n_quanta=1),
        StressLoadedConfiguration(n_quanta=2),
    ]
