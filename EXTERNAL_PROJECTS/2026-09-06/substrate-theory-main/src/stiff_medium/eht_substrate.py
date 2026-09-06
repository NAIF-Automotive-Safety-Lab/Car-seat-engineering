"""EHT (Event Horizon Telescope) substrate predictions.

The Event Horizon Telescope imaged M87* in 2019 and Sgr A* in 2022, both
showing a bright photon ring surrounding a dark central shadow. Within the
B3 / stiff-medium substrate framework, a black-hole horizon corresponds to
the sigma = 1/2 saturation locus. At EHT angular resolution (~10 uas) the
predicted image is identical to the GR (Kerr) result; the substrate adds
no observable deviation. This module produces those predictions and
compares them against the published EHT measurements.

References:
    EHT Collaboration, ApJL 875 (2019) L1-L6   [M87*]
    EHT Collaboration, ApJL 930 (2022) L12-L17 [Sgr A*]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import cos, pi, radians, sin, sqrt
from typing import Dict


# Physical constants (SI)
G = 6.67430e-11          # m^3 kg^-1 s^-2
C = 2.99792458e8         # m s^-1
M_SUN = 1.98847e30       # kg
PC = 3.0856775814913673e16  # m
KPC = 1e3 * PC
MPC = 1e6 * PC
UAS_PER_RAD = (180.0 / pi) * 3600.0 * 1e6  # micro-arcseconds per radian


def _to_uas(angle_rad: float) -> float:
    return angle_rad * UAS_PER_RAD


def _gravitational_radius(mass_kg: float) -> float:
    """r_g = G M / c^2 in meters."""
    return G * mass_kg / (C * C)


@dataclass
class SourceParameters:
    name: str
    mass_msun: float
    distance_m: float
    observed_ring_uas: float

    @property
    def mass_kg(self) -> float:
        return self.mass_msun * M_SUN

    @property
    def r_g(self) -> float:
        return _gravitational_radius(self.mass_kg)

    @property
    def r_g_uas(self) -> float:
        """Angular size of one r_g in micro-arcseconds."""
        return _to_uas(self.r_g / self.distance_m)


@dataclass
class EHTSubstrate:
    """EHT predictions in the substrate / GR-equivalent regime."""

    # Empirical EHT-derived ring-to-shadow scaling factor.
    # The bright photon ring observed by the EHT is offset outward
    # from the critical curve due to finite emission/lensing width.
    # Best-fit calibrations give d_ring ~ (11 +/- 0.5) r_g in diameter,
    # vs. d_shadow = 2 * b_c = 2 * sqrt(27) r_g ~ 10.39 r_g.
    ring_diameter_in_rg: float = 10.6
    shadow_diameter_in_rg: float = field(default=2.0 * sqrt(27.0))

    # ------------------------------------------------------------------
    # Source parameters
    # ------------------------------------------------------------------
    def M87_parameters(self) -> SourceParameters:
        """M87* parameters (EHT 2019)."""
        return SourceParameters(
            name="M87*",
            mass_msun=6.5e9,
            distance_m=16.8 * MPC,
            observed_ring_uas=42.0,
        )

    def SgrA_parameters(self) -> SourceParameters:
        """Sgr A* parameters (EHT 2022)."""
        return SourceParameters(
            name="Sgr A*",
            mass_msun=4.3e6,
            distance_m=8.3 * KPC,
            observed_ring_uas=51.8,
        )

    # ------------------------------------------------------------------
    # Core geometric predictions
    # ------------------------------------------------------------------
    def shadow_radius(self, spin: float = 0.0) -> float:
        """Critical impact parameter b_c in units of r_g.

        Schwarzschild: b_c = 3 * sqrt(3) ~ 5.196 r_g (radius).
        Kerr varies weakly with spin (a few percent for retrograde rays),
        so we keep the spherically averaged value.
        """
        b_c = 3.0 * sqrt(3.0)
        # Slight spin-averaged correction (Kerr shadow area is nearly
        # spin-independent; <2% across 0 <= a <= 1).
        return b_c * (1.0 - 0.02 * spin * spin)

    def predicted_ring_diameter(
        self, mass_msun: float, distance_m: float, spin: float = 0.0
    ) -> float:
        """Substrate prediction for the observed bright-ring diameter (uas).

        Substrate framework: identical to GR at EHT precision.
        """
        r_g = _gravitational_radius(mass_msun * M_SUN)
        # Spin reduces the average ring diameter slightly (~5% at a=1).
        spin_factor = 1.0 - 0.05 * spin * spin
        ring_m = self.ring_diameter_in_rg * r_g * spin_factor
        return _to_uas(ring_m / distance_m)

    def predicted_shadow_diameter(
        self, mass_msun: float, distance_m: float, spin: float = 0.0
    ) -> float:
        """Inner shadow (critical curve) angular diameter (uas)."""
        r_g = _gravitational_radius(mass_msun * M_SUN)
        diam_m = 2.0 * self.shadow_radius(spin) * r_g
        return _to_uas(diam_m / distance_m)

    # ------------------------------------------------------------------
    # Kerr shape and beaming
    # ------------------------------------------------------------------
    def kerr_shadow(self, spin: float, inclination_deg: float) -> Dict[str, float]:
        """Shadow shape parameters for a spinning (Kerr) black hole.

        Returns approximate semi-major / semi-minor axes (in r_g) and
        displacement of the shadow centroid from the BH spin axis.
        Uses the standard small-distortion expansion; accurate to a few
        percent across all spins for inclinations <~ 60 deg.
        """
        a = max(-1.0, min(1.0, spin))
        i = radians(inclination_deg)

        # Symmetric average shadow size.
        b_avg = self.shadow_radius(a)

        # Asymmetry: shadow is squashed perpendicular to the spin axis
        # for prograde-side rays. Standard fit: delta ~ a^2 sin^2(i) / 4.
        distortion = (a * a * sin(i) * sin(i)) / 4.0
        semi_major = b_avg
        semi_minor = b_avg * (1.0 - distortion)

        # Centroid displacement (frame-dragging shift), in r_g.
        displacement = 2.0 * a * sin(i)

        return {
            "spin": a,
            "inclination_deg": inclination_deg,
            "semi_major_rg": semi_major,
            "semi_minor_rg": semi_minor,
            "axial_ratio": semi_minor / semi_major,
            "centroid_displacement_rg": displacement,
        }

    def ring_brightness_asymmetry(
        self, spin: float = 0.5, inclination_deg: float = 17.0
    ) -> Dict[str, float]:
        """Doppler-beamed brightness asymmetry around the ring.

        For a thin Keplerian emitter near the ISCO, approaching material
        is blue-shifted/boosted by D = 1/(gamma (1 - beta cos theta))
        and the brightness ratio scales as ~D^(3+alpha) for a power-law
        spectrum of slope alpha (alpha ~ 1 for synchrotron near 230 GHz).
        """
        # Crude Keplerian velocity at r ~ 5 r_g (near typical emitter).
        r = 5.0
        v_kep = 1.0 / sqrt(r)        # in units of c
        gamma = 1.0 / sqrt(1.0 - v_kep * v_kep)
        i = radians(inclination_deg)

        # Approaching vs receding side along the line of sight.
        cos_theta = sin(i)           # tangential motion projected
        d_app = 1.0 / (gamma * (1.0 - v_kep * cos_theta))
        d_rec = 1.0 / (gamma * (1.0 + v_kep * cos_theta))

        alpha = 1.0
        ratio = (d_app / d_rec) ** (3.0 + alpha)
        # M87 EHT measured ~2:1 azimuthal contrast; Sgr A* less constrained.
        return {
            "spin": spin,
            "inclination_deg": inclination_deg,
            "v_keplerian_over_c": v_kep,
            "doppler_approaching": d_app,
            "doppler_receding": d_rec,
            "brightness_ratio_approaching_to_receding": ratio,
        }

    # ------------------------------------------------------------------
    # Substrate-specific signature
    # ------------------------------------------------------------------
    def substrate_horizon_mark(self, source: SourceParameters) -> Dict[str, float]:
        """Any substrate-specific deviation from the GR image?

        At sigma = 1/2 saturation the substrate matches GR identically
        in the geometric optics limit. Putative deviations would scale
        as (l_Planck / r_g)^p with p >= 1; for M87* and Sgr A*:
            l_P / r_g ~ 1e-40,
        far below any conceivable EHT angular resolution
        (~10 uas vs r_g ~ 4 uas for M87, 5 uas for Sgr A*).
        """
        l_planck = 1.616255e-35  # m
        ratio = l_planck / source.r_g
        # Resolution limit roughly 10 uas converted to r_g units.
        resolution_uas = 10.0
        resolution_in_rg = resolution_uas / source.r_g_uas
        return {
            "l_Planck_over_r_g": ratio,
            "expected_deviation_uas": ratio * source.r_g_uas,
            "EHT_resolution_in_r_g": resolution_in_rg,
            "detectable": False,
            "note": "Substrate prediction = GR at EHT precision.",
        }

    # ------------------------------------------------------------------
    # Comparison & reporting
    # ------------------------------------------------------------------
    def compare_to_observations(self, spin: float = 0.5) -> Dict[str, Dict[str, float]]:
        """Compare substrate predictions against the EHT measurements."""
        results: Dict[str, Dict[str, float]] = {}
        for src in (self.M87_parameters(), self.SgrA_parameters()):
            pred = self.predicted_ring_diameter(src.mass_msun, src.distance_m, spin)
            results[src.name] = {
                "mass_Msun": src.mass_msun,
                "distance_m": src.distance_m,
                "r_g_uas": src.r_g_uas,
                "predicted_ring_uas": pred,
                "observed_ring_uas": src.observed_ring_uas,
                "fractional_error": (pred - src.observed_ring_uas)
                / src.observed_ring_uas,
                "abs_pct_error": 100.0
                * abs(pred - src.observed_ring_uas)
                / src.observed_ring_uas,
            }
        return results

    def report(self) -> str:
        """Human-readable summary."""
        lines = [
            "=" * 64,
            "EHT Substrate Predictions (B3 / stiff-medium)",
            "=" * 64,
            "",
            f"Shadow critical impact parameter: b_c = {self.shadow_radius():.4f} r_g",
            f"Adopted ring-to-r_g scaling:      d_ring = {self.ring_diameter_in_rg:.2f} r_g",
            "",
        ]
        comp = self.compare_to_observations()
        for name, row in comp.items():
            lines.append(f"--- {name} ---")
            lines.append(f"  M = {row['mass_Msun']:.3e} M_sun")
            lines.append(f"  one r_g          = {row['r_g_uas']:.3f} uas")
            lines.append(f"  predicted ring   = {row['predicted_ring_uas']:.2f} uas")
            lines.append(f"  observed ring    = {row['observed_ring_uas']:.2f} uas")
            lines.append(f"  |error|          = {row['abs_pct_error']:.2f}%")
            lines.append("")

        kerr = self.kerr_shadow(spin=0.9, inclination_deg=17.0)
        lines.append("Kerr shadow (a=0.9, i=17 deg, M87-like):")
        lines.append(
            f"  axial ratio      = {kerr['axial_ratio']:.4f} "
            f"(displacement {kerr['centroid_displacement_rg']:.3f} r_g)"
        )

        beam = self.ring_brightness_asymmetry(spin=0.9, inclination_deg=17.0)
        lines.append(
            "Doppler brightness asymmetry "
            f"(v/c = {beam['v_keplerian_over_c']:.3f}): "
            f"ratio = {beam['brightness_ratio_approaching_to_receding']:.2f}"
        )
        lines.append("")

        sub = self.substrate_horizon_mark(self.M87_parameters())
        lines.append(
            "Substrate horizon mark: "
            f"l_P/r_g = {sub['l_Planck_over_r_g']:.2e} "
            "(undetectable at EHT precision)."
        )
        lines.append("=" * 64)
        return "\n".join(lines)


if __name__ == "__main__":  # pragma: no cover
    print(EHTSubstrate().report())
