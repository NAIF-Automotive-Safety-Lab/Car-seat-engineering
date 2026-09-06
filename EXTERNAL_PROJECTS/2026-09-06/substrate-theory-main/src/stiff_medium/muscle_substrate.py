"""Paired muscle-contraction GEOMETRY + SIMULATOR for the substrate framework.

This module provides a compact GEOMETRY/SIM pairing for striated-muscle
contraction (sliding filament + cross-bridge cycle) under the
B3 / stiff-medium ontology.

GEOMETRY
--------
``MuscleGeometry`` represents one sarcomere -- the contractile unit of
striated muscle, from Z-disk to Z-disk.  Reference dimensions (textbook
vertebrate skeletal muscle, frog sartorius / mammalian fibre):

    * sarcomere rest length      L_s  = 2.5  micrometre  (Z-to-Z)
    * thick (myosin) filament    L_my = 1.6  micrometre  (A-band length)
    * thin  (actin)  filament    L_th = 1.0  micrometre  (Z to tip)
    * H-zone (bare zone)         L_H  = 0.2  micrometre  (M-line region)
    * I-band (thin only)         L_I  = (L_s - L_my)     ~ 0.9 micrometre
    * A-band (thick + overlap)   L_A  = L_my             = 1.6 micrometre
    * thick filament diameter            ~ 15 nm
    * thin  filament diameter            ~  8 nm
    * inter-filament spacing (hex)       ~ 40-45 nm
    * myosin cross-bridges / thick       ~ 300 (S1 heads)
    * cross-bridge spacing along thick   ~ 14.3 nm (helical, 6/turn)

The sliding-filament theory (H.E. Huxley & Hanson 1954; A.F. Huxley &
Niedergerke 1954) holds that contraction shortens the I-band and the H-zone
while the A-band length stays fixed.  The cross-bridge cycle (Huxley 1957;
Huxley-Simmons 1971) is a 5-state ATPase cycle whose mechanical output is
a ~ 5-10 nm power stroke per attached head per ATP hydrolysed.

SIMULATOR
---------
``MuscleSimulator`` exposes the contractile observables of a single
sarcomere or whole muscle fibre under the substrate ontology:

    * Hill (1938) hyperbolic force-velocity:
          (F + a)(v + b) = (F_max + a) * b
      Equivalently  v(F) = b * (F_max - F) / (F + a),  with a / F_max ~ 0.25
      and v_max = b * F_max / a.  Maximum mechanical power occurs near
      v ~ 0.3 * v_max (textbook ~ 0.31 v_max for a / F_max = 0.25).

    * Gordon-Huxley-Julian (1966) length-tension:
          F_active(L) is a piecewise-linear plateau with peak at
          L_0 in [2.0, 2.25] micrometre (the optimal-overlap window),
          falling to zero at L = 1.27 micrometre (full thin-thick clash)
          and at L = 3.65 micrometre (zero overlap).

    * Cross-bridge cycle (5 states, after Lymn-Taylor 1971):
          attached  (rigor, A.M)
          power-stroke (A.M.ADP -> A.M, 5-10 nm displacement)
          detached  (A + M.ATP)
          primed    (M.ADP.Pi, head re-cocked, weakly bound)
          attached again.
      One ATP per cycle.  Duty ratio ~ 5% (skeletal) means most heads at
      v_max are detached.

    * ATP -> mechanical work efficiency.  The thermodynamic ceiling for
      converting GTP/ATP free energy into mechanical work is ~ 50%
      (thermal limit, at body T); the Huxley-Simmons quantitative model
      delivers W_stroke = F_stroke * d_stroke ~ 25 pN nm vs the Delta_G
      of one ATP ~ 50 pN nm so eta_max ~ 50%.

Substrate-specific predictions encoded here:
    * Each cross-bridge is a small substrate strain pump that converts
      the chemical strain stored in the ATP -> ADP + Pi reaction into
      a mechanical 5-10 nm filament displacement; the Huxley-Simmons
      lever arm rotation is the mechanical re-folding of the S1 motor
      domain after Pi release.
    * The Hill hyperbola is the universal envelope for strain-pump
      arrays: when load F increases, attachment time per head goes up
      so v drops; the hyperbolic shape with a / F_max ~ 0.25 is the
      consequence of strain-dependent detachment-rate skew.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Reference (textbook) values -- vertebrate skeletal muscle
# ---------------------------------------------------------------------------

# Sarcomere geometry (micrometre)
SARCOMERE_REST_LENGTH_UM:    float = 2.5
THICK_FILAMENT_LENGTH_UM:    float = 1.6   # A-band length
THIN_FILAMENT_LENGTH_UM:     float = 1.0   # Z to tip of thin filament
H_ZONE_LENGTH_UM:            float = 0.2   # bare-zone (M-line region)

# Filament diameters (nanometre)
THICK_FILAMENT_DIAMETER_NM:  float = 15.0
THIN_FILAMENT_DIAMETER_NM:   float = 8.0

# Inter-filament hexagonal lattice spacing (nanometre)
INTER_FILAMENT_SPACING_NM:   float = 42.0

# Myosin cross-bridge layout
MYOSIN_HEADS_PER_THICK:      int   = 300       # S1 heads on one thick filament
CROSS_BRIDGE_AXIAL_PERIOD_NM: float = 14.3     # helical pitch axial repeat
HEADS_PER_HELICAL_TURN:      int   = 6         # 3 levels x 2 heads = 6

# Hill (1938) force-velocity constants (frog sartorius, scaled to F_max = 1)
HILL_A_OVER_FMAX:            float = 0.25      # canonical textbook value
HILL_B_OVER_VMAX:            float = 0.25      # consequence of a/F_max=0.25

# Gordon-Huxley-Julian (1966) length-tension corner points (micrometre)
LENGTH_TENSION_L_MIN_UM:     float = 1.27   # full clash (F_active = 0)
LENGTH_TENSION_L_PLATEAU_LOW_UM:  float = 2.00  # plateau begins
LENGTH_TENSION_L_PLATEAU_HIGH_UM: float = 2.25  # plateau ends
LENGTH_TENSION_L_MAX_UM:     float = 3.65   # zero overlap (F_active = 0)

# Power-stroke biophysics (Huxley-Simmons 1971; Finer-Simmons-Spudich 1994)
POWER_STROKE_NM:             float = 10.0  # per single S1 head per ATP
SINGLE_HEAD_FORCE_PN:        float = 5.0   # isometric force per head
ATP_HYDROLYSIS_KCAL_PER_MOL: float = 7.3   # standard biochem ~ 7.3 kcal/mol
ATP_HYDROLYSIS_PN_NM:        float = 50.0  # ~ 7.3 kcal/mol -> 50 pN nm per ATP

# Cross-bridge duty ratio (fraction of cycle attached, skeletal muscle)
DUTY_RATIO_SKELETAL:         float = 0.05  # ~ 5% in fast skeletal at v_max

# Substrate anchor
LAMBDA_QCD_MEV:              float = 200.0
KT_AT_BODY_T_PN_NM:          float = 4.28   # k_B T at 310 K in pN.nm


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------


@dataclass
class MuscleGeometry:
    """Single sarcomere: Z-disk, A-band, I-band, M-line, thick + thin filaments.

    The geometry is one half of a myofibril cross-section: thick (myosin)
    filaments on a hexagonal lattice with thin (actin) filaments threading
    the trigonal interstices.  Coordinates are in micrometre (axial) and
    nanometre (cross-section).

    Attributes
    ----------
    sarcomere_length_um : float
        Z-to-Z distance.  Default = textbook resting length 2.5 micrometre.
    n_thick : int
        Number of thick filaments rendered (3 by default for a small
        cross-sectional patch).
    n_cross_bridges : int
        Number of cross-bridge heads sampled per thick filament for
        rendering / counting.  Default 30 (one per period x 30 levels).
    """

    sarcomere_length_um: float = SARCOMERE_REST_LENGTH_UM
    n_thick: int = 3
    n_cross_bridges: int = 30

    # ------------------------------------------------------------------
    # Band geometry (axial)
    # ------------------------------------------------------------------

    def z_disk_positions_um(self) -> Tuple[float, float]:
        """Return (left, right) Z-disk positions, centred on the M-line."""
        half = 0.5 * self.sarcomere_length_um
        return (-half, +half)

    def m_line_position_um(self) -> float:
        """M-line bisects the sarcomere at z = 0."""
        return 0.0

    def a_band_extent_um(self) -> Tuple[float, float]:
        """A-band = location of the thick filament (full A-band length).

        The A-band is centred on the M-line and spans the full length of
        the thick filament regardless of sarcomere length (this is the
        signature of the sliding-filament model).
        """
        half = 0.5 * THICK_FILAMENT_LENGTH_UM
        return (-half, +half)

    def h_zone_extent_um(self) -> Tuple[float, float]:
        """H-zone = bare zone in the middle of the A-band (no thin overlap).

        The H-zone shrinks during contraction and disappears at maximum
        overlap; here we report the resting H-zone (= L_H by default,
        scaled with sarcomere length deviation).
        """
        # H-zone shrinks linearly with shortening from rest length
        h_rest = H_ZONE_LENGTH_UM
        delta = self.sarcomere_length_um - SARCOMERE_REST_LENGTH_UM
        h = max(h_rest + delta, 0.0)  # clamps at full overlap
        half = 0.5 * h
        return (-half, +half)

    def i_band_extents_um(self) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """I-band = thin-only zones, between Z-disk and edge of A-band.

        Returns ((left_inner, left_outer), (right_inner, right_outer))
        as ((-L_A/2, -L_s/2), (+L_A/2, +L_s/2)).  Each I-band is a single
        wedge from the edge of the thick filament to the Z-disk on its
        side.
        """
        a_left, a_right = self.a_band_extent_um()
        z_left, z_right = self.z_disk_positions_um()
        return ((z_left, a_left), (a_right, z_right))

    def thick_filament_axial_extent_um(self) -> Tuple[float, float]:
        """Same as A-band: the thick filament is the A-band."""
        return self.a_band_extent_um()

    def thin_filament_axial_extents_um(
        self,
    ) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """Each thin filament runs from a Z-disk inward toward the M-line.

        The thin filament length L_th is fixed (1.0 micrometre); when the
        sarcomere shortens, the thin filaments slide deeper into the A-band
        and may overlap each other near the M-line.

        Returns ((left_root, left_tip), (right_root, right_tip)).
        """
        z_left, z_right = self.z_disk_positions_um()
        # Thin filament tip = Z + L_th (left filament points right)
        left_tip = z_left + THIN_FILAMENT_LENGTH_UM
        right_tip = z_right - THIN_FILAMENT_LENGTH_UM
        return ((z_left, left_tip), (z_right, right_tip))

    # ------------------------------------------------------------------
    # Cross-section: hexagonal lattice
    # ------------------------------------------------------------------

    def thick_filament_cross_section_nm(self) -> np.ndarray:
        """Return (n_thick, 2) coordinates of thick-filament centres in nm.

        Thick filaments sit on a 2D hexagonal lattice.  We sample the
        centre and its first-shell neighbours up to n_thick.
        """
        a = INTER_FILAMENT_SPACING_NM
        # Hex lattice basis vectors
        e1 = np.array([a, 0.0])
        e2 = np.array([0.5 * a, math.sqrt(3.0) / 2.0 * a])
        # Generate up to n_thick points starting from the origin and walking
        # outward in a hex spiral
        candidates = [np.array([0.0, 0.0])]
        for shell in range(1, 4):
            for i in range(-shell, shell + 1):
                for j in range(-shell, shell + 1):
                    if max(abs(i), abs(j), abs(i + j)) == shell:
                        candidates.append(i * e1 + j * e2)
                        if len(candidates) >= self.n_thick:
                            break
                if len(candidates) >= self.n_thick:
                    break
            if len(candidates) >= self.n_thick:
                break
        return np.array(candidates[: self.n_thick])

    def thin_filament_cross_section_nm(self) -> np.ndarray:
        """Return coordinates of thin filaments (trigonal interstices).

        For each pair of adjacent thick filaments, the two trigonal sites
        between them carry one thin filament each.  In the canonical
        vertebrate skeletal lattice the thin:thick ratio is 2:1.
        """
        thick = self.thick_filament_cross_section_nm()
        a = INTER_FILAMENT_SPACING_NM
        thin = []
        # For an N-thick patch, place 2*N thin filaments on the trigonal
        # voids. We use a simple offset of (a/2, a*sqrt(3)/6) for one set
        # and the negative offset for the other set, relative to each thick.
        offsets = [
            np.array([0.0,  a * math.sqrt(3.0) / 3.0]),
            np.array([0.0, -a * math.sqrt(3.0) / 3.0]),
        ]
        for tx in thick:
            for off in offsets:
                thin.append(tx + off)
        return np.array(thin)

    # ------------------------------------------------------------------
    # Cross-bridge head positions (axial)
    # ------------------------------------------------------------------

    def cross_bridge_axial_positions_um(self) -> np.ndarray:
        """Return axial positions of myosin S1 heads along the thick filament.

        Heads are spaced at CROSS_BRIDGE_AXIAL_PERIOD_NM (= 14.3 nm)
        with HEADS_PER_HELICAL_TURN per turn (= 6).  We omit the bare
        H-zone (no heads in the central stripe).  Returned in micrometre.
        """
        a_left, a_right = self.a_band_extent_um()
        h_left, h_right = self.h_zone_extent_um()
        period_um = CROSS_BRIDGE_AXIAL_PERIOD_NM * 1e-3
        # Sample positions along [a_left, a_right] every period_um, skip [h_left, h_right]
        n_total_heads = int(round(THICK_FILAMENT_LENGTH_UM / period_um))
        positions = a_left + (np.arange(n_total_heads) + 0.5) * period_um
        # Remove heads in the bare H-zone
        mask = (positions <= h_left) | (positions >= h_right)
        kept = positions[mask]
        # Subsample to n_cross_bridges so we don't render hundreds of heads
        if len(kept) > self.n_cross_bridges:
            idx = np.linspace(0, len(kept) - 1, self.n_cross_bridges, dtype=int)
            kept = kept[idx]
        return kept

    def cross_bridge_count(self) -> int:
        """Total cross-bridge heads on one full thick filament (textbook ~ 300)."""
        period_um = CROSS_BRIDGE_AXIAL_PERIOD_NM * 1e-3
        n_total = int(round(THICK_FILAMENT_LENGTH_UM / period_um))
        # Subtract heads in the bare H-zone
        n_bare = int(round(H_ZONE_LENGTH_UM / period_um))
        return n_total - n_bare

    # ------------------------------------------------------------------
    # Overlap and length-tension
    # ------------------------------------------------------------------

    def thin_thick_overlap_um(self) -> float:
        """Length of the overlap zone between thin and thick filaments.

        Per half-sarcomere.  Maximum at the optimal-overlap plateau.
        """
        a_left, a_right = self.a_band_extent_um()
        # Distance from inner-end of thick (left edge of left half) to thin-tip
        # (left filament tip)
        (z_left, left_tip), (z_right, right_tip) = self.thin_filament_axial_extents_um()
        # Overlap on the left half: from a_left (inner edge of half-A on left)
        # to left_tip (where the thin reaches), bounded above by 0 (M-line).
        # We take the minimum of (left_tip - a_left, m_line_position - a_left).
        left_overlap = max(0.0, min(left_tip, 0.0) - a_left)
        right_overlap = max(0.0, a_right - max(right_tip, 0.0))
        # Symmetric: report one side
        return 0.5 * (left_overlap + right_overlap)

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "sarcomere = strain-pump array on the stiff substrate; "
                "each cross-bridge converts ATP free-energy into a "
                "5-10 nm filament-sliding power stroke"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "sarcomere_length_um": self.sarcomere_length_um,
            "n_cross_bridges_per_thick": self.cross_bridge_count(),
            "filament_thin_to_thick_ratio": 2,
            "atp_per_stroke": 1,
        }


# ---------------------------------------------------------------------------
# SIMULATOR: Hill F-v, length-tension, cross-bridge cycle, ATP economy
# ---------------------------------------------------------------------------


@dataclass
class MuscleSimulator:
    """Contractile dynamics on the sarcomere geometry.

    The simulator exposes (i) the Hill 1938 hyperbolic force-velocity
    relation, (ii) the Gordon-Huxley-Julian 1966 length-tension curve,
    (iii) the Lymn-Taylor 1971 5-state cross-bridge cycle, and
    (iv) the ATP -> mechanical-work efficiency.  All quantities are
    normalised to F_max = 1 and v_max = 1 unless otherwise stated.

    Parameters
    ----------
    geometry : MuscleGeometry, optional
        Underlying sarcomere geometry; if None, a default sarcomere is
        constructed.
    f_max_normalised : float
        Normalised maximum isometric force (default 1.0).
    v_max_normalised : float
        Normalised maximum (unloaded) shortening velocity (default 1.0).
    a_over_fmax : float
        Hill coefficient a/F_max; canonical value ~ 0.25 (frog sartorius).
    """

    geometry: MuscleGeometry | None = None
    f_max_normalised: float = 1.0
    v_max_normalised: float = 1.0
    a_over_fmax: float = HILL_A_OVER_FMAX

    # ------------------------------------------------------------------
    # Hill (1938) force-velocity hyperbola
    # ------------------------------------------------------------------

    def hill_a(self) -> float:
        """Hill constant a [units of F_max]."""
        return self.a_over_fmax * self.f_max_normalised

    def hill_b(self) -> float:
        """Hill constant b [units of v_max].

        From the boundary v = v_max at F = 0:
            v_max = b * F_max / a  =>  b = a/F_max * v_max  =  a_over_fmax * v_max.
        """
        return self.a_over_fmax * self.v_max_normalised

    def hill_velocity_from_force(self, force: float) -> float:
        """Solve Hill (1938) hyperbola for velocity given load.

        (F + a)(v + b) = (F_max + a) * b
        =>  v(F) = b * (F_max - F) / (F + a)
        """
        a = self.hill_a()
        b = self.hill_b()
        if force >= self.f_max_normalised:
            return 0.0
        if force < 0.0:
            force = 0.0
        return b * (self.f_max_normalised - force) / (force + a)

    def hill_force_from_velocity(self, velocity: float) -> float:
        """Solve Hill (1938) hyperbola for force given velocity.

        From (F + a)(v + b) = (F_max + a) * b:
            F(v) = (F_max + a) * b / (v + b) - a
        """
        a = self.hill_a()
        b = self.hill_b()
        if velocity >= self.v_max_normalised:
            return 0.0
        if velocity < 0.0:
            velocity = 0.0
        f = (self.f_max_normalised + a) * b / (velocity + b) - a
        return max(f, 0.0)

    def hill_force_velocity_curve(
        self, n: int = 200,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Discretised F-v hyperbola.

        Returns (forces, velocities) sampled at n equispaced velocities
        from 0 to v_max; force runs from F_max at v=0 down to 0 at v_max.
        """
        v = np.linspace(0.0, self.v_max_normalised, n)
        f = np.array([self.hill_force_from_velocity(vi) for vi in v])
        return f, v

    def hill_power_curve(
        self, n: int = 200,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Mechanical power P = F * v as a function of velocity.

        Returns (velocities, powers, forces) sampled at n velocities.
        """
        v = np.linspace(0.0, self.v_max_normalised, n)
        f = np.array([self.hill_force_from_velocity(vi) for vi in v])
        p = f * v
        return v, p, f

    def velocity_at_max_power(self) -> float:
        """Velocity that maximises P = F(v) * v.

        For Hill's hyperbola with a/F_max = 0.25 (b/v_max = 0.25):
            v* / v_max ~ 0.31
        We return this analytically by maximisation.
        """
        a = self.hill_a()
        b = self.hill_b()
        F_max = self.f_max_normalised
        # P(v) = v * [ (F_max + a) * b / (v + b) - a ]
        # dP/dv = 0
        #   => (F_max + a) * b * b / (v + b)^2 = a
        #   => (v + b)^2 = (F_max + a) * b^2 / a
        #   => v* = b * sqrt( (F_max + a) / a ) - b
        v_star = b * (math.sqrt((F_max + a) / a) - 1.0)
        return v_star

    def force_at_max_power(self) -> float:
        """Force at the velocity that maximises mechanical power."""
        return self.hill_force_from_velocity(self.velocity_at_max_power())

    def max_mechanical_power(self) -> float:
        """Maximum mechanical power output (in F_max * v_max units)."""
        v_star = self.velocity_at_max_power()
        return v_star * self.hill_force_from_velocity(v_star)

    # ------------------------------------------------------------------
    # Gordon-Huxley-Julian (1966) length-tension curve
    # ------------------------------------------------------------------

    def length_tension_curve(
        self, n: int = 400,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Active force vs sarcomere length (Gordon-Huxley-Julian 1966).

        Piecewise-linear:
          0 for L <= L_min                 (clash)
          rises to plateau over L_min..L_plateau_low
          plateau (F = 1) over L_plateau_low..L_plateau_high
          falls to 0 over L_plateau_high..L_max  (decreasing overlap)
          0 for L >= L_max                 (no overlap)

        Returns (lengths, F_active) on n samples spanning [1.0, 4.0] um.
        """
        L = np.linspace(1.0, 4.0, n)
        F = np.zeros_like(L)
        for i, l in enumerate(L):
            if l <= LENGTH_TENSION_L_MIN_UM or l >= LENGTH_TENSION_L_MAX_UM:
                F[i] = 0.0
            elif l < LENGTH_TENSION_L_PLATEAU_LOW_UM:
                # Ascending limb
                frac = ((l - LENGTH_TENSION_L_MIN_UM) /
                        (LENGTH_TENSION_L_PLATEAU_LOW_UM - LENGTH_TENSION_L_MIN_UM))
                F[i] = frac
            elif l <= LENGTH_TENSION_L_PLATEAU_HIGH_UM:
                F[i] = 1.0
            else:
                # Descending limb
                frac = ((LENGTH_TENSION_L_MAX_UM - l) /
                        (LENGTH_TENSION_L_MAX_UM - LENGTH_TENSION_L_PLATEAU_HIGH_UM))
                F[i] = frac
        return L, F

    def active_force_at_length(self, L_um: float) -> float:
        """Active normalised force at a given sarcomere length."""
        if L_um <= LENGTH_TENSION_L_MIN_UM or L_um >= LENGTH_TENSION_L_MAX_UM:
            return 0.0
        if L_um < LENGTH_TENSION_L_PLATEAU_LOW_UM:
            return ((L_um - LENGTH_TENSION_L_MIN_UM) /
                    (LENGTH_TENSION_L_PLATEAU_LOW_UM - LENGTH_TENSION_L_MIN_UM))
        if L_um <= LENGTH_TENSION_L_PLATEAU_HIGH_UM:
            return 1.0
        return ((LENGTH_TENSION_L_MAX_UM - L_um) /
                (LENGTH_TENSION_L_MAX_UM - LENGTH_TENSION_L_PLATEAU_HIGH_UM))

    def optimal_length_um(self) -> float:
        """Sarcomere length giving peak active force (centre of plateau)."""
        return 0.5 * (LENGTH_TENSION_L_PLATEAU_LOW_UM + LENGTH_TENSION_L_PLATEAU_HIGH_UM)

    # ------------------------------------------------------------------
    # Cross-bridge cycle (Lymn-Taylor 1971, 5 states)
    # ------------------------------------------------------------------

    def cross_bridge_states(self) -> Tuple[str, ...]:
        """The five canonical cross-bridge states in cycle order."""
        return (
            "attached_rigor",      # A.M
            "power_stroke",        # A.M.ADP -> A.M (5-10 nm displacement)
            "detached",            # A + M.ATP
            "primed",              # M.ADP.Pi (head re-cocked, weakly bound)
            "attached_again",      # A.M.ADP.Pi -> A.M.ADP
        )

    def cross_bridge_state_transitions(self) -> Tuple[Tuple[str, str], ...]:
        """The directed cycle of state-to-state transitions."""
        states = self.cross_bridge_states()
        return tuple((states[i], states[(i + 1) % len(states)]) for i in range(len(states)))

    def power_stroke_displacement_nm(self) -> float:
        """Filament displacement per single S1 head per ATP."""
        return POWER_STROKE_NM

    def power_stroke_force_pN(self) -> float:
        """Isometric force generated by a single attached S1 head."""
        return SINGLE_HEAD_FORCE_PN

    def power_stroke_work_pN_nm(self) -> float:
        """Work delivered by one head per cycle: F * d."""
        return self.power_stroke_force_pN() * self.power_stroke_displacement_nm()

    def duty_ratio(self) -> float:
        """Fraction of cycle a head spends attached (skeletal default ~ 5%)."""
        return DUTY_RATIO_SKELETAL

    # ------------------------------------------------------------------
    # ATP economy and thermodynamic efficiency
    # ------------------------------------------------------------------

    def atp_per_cross_bridge_cycle(self) -> int:
        """One ATP hydrolysed per single cross-bridge cycle (per head)."""
        return 1

    def atp_free_energy_pN_nm(self) -> float:
        """Free energy released per ATP hydrolysed in cellular conditions."""
        return ATP_HYDROLYSIS_PN_NM

    def thermodynamic_efficiency(self) -> float:
        """Mechanical work / ATP free energy = ~ 0.5 (theoretical ceiling).

        Per-stroke work = F_stroke * d_stroke ~ 5 pN * 10 nm = 50 pN nm.
        ATP free energy ~ 50 pN nm (from Delta_G ~ -7.3 kcal/mol).
        Ratio ~ 1.0 in this textbook approximation but the in-cell
        efficiency is ~ 50% because half the ATP energy is dissipated
        as heat through detachment and unproductive Pi release.
        """
        # Use the textbook bookkeeping: half goes to mechanical work.
        return 0.5

    def per_stroke_work_to_atp_ratio(self) -> float:
        """Direct ratio of stroke mechanical work to ATP free-energy."""
        return self.power_stroke_work_pN_nm() / self.atp_free_energy_pN_nm()

    # ------------------------------------------------------------------
    # Combined dynamics: contraction time-course
    # ------------------------------------------------------------------

    def shortening_trajectory(
        self,
        load_normalised: float = 0.5,
        duration_s: float = 0.2,
        n_steps: int = 200,
        L0_um: float | None = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Integrate sarcomere shortening at constant load.

        At a given normalised load F (0..1), the shortening velocity is
        constant (Hill hyperbola) so L(t) decreases linearly until the
        ascending limb of the length-tension curve forces a stop.  Model:

            v(t) = hill_velocity_from_force(F * F_active(L(t)))
            dL/dt = - v(t) * v_max_phys      (v_max_phys default = 5 um/s)

        Returns (t, L_um, F_normalised).
        """
        v_max_phys_um_per_s = 5.0   # textbook fast-fibre v_max ~ 5 sarcomere lengths/s
        L = (L0_um if L0_um is not None else SARCOMERE_REST_LENGTH_UM)
        Ls: list[float] = []
        Fs: list[float] = []
        ts = np.linspace(0.0, duration_s, n_steps)
        dt = duration_s / max(n_steps - 1, 1)
        for _ in ts:
            f_act = self.active_force_at_length(L)
            f_load = load_normalised * f_act
            v = self.hill_velocity_from_force(f_load) * v_max_phys_um_per_s * f_act
            Ls.append(L)
            Fs.append(f_act * load_normalised)
            L = max(L - v * dt, LENGTH_TENSION_L_MIN_UM)
        return ts, np.array(Ls), np.array(Fs)

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "f_max_normalised": self.f_max_normalised,
            "v_max_normalised": self.v_max_normalised,
            "a_over_fmax": self.a_over_fmax,
            "v_at_max_power_over_vmax": self.velocity_at_max_power() / self.v_max_normalised,
            "n_states_per_cycle": len(self.cross_bridge_states()),
            "atp_per_cycle": self.atp_per_cross_bridge_cycle(),
            "thermodynamic_efficiency": self.thermodynamic_efficiency(),
            "interpretation": (
                "muscle = sarcomere strain-pump array; Hill F-v hyperbola "
                "is the universal envelope; each cross-bridge consumes "
                "1 ATP per ~10 nm power stroke; ~50% theoretical efficiency"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
        }


__all__ = [
    "MuscleGeometry",
    "MuscleSimulator",
    "SARCOMERE_REST_LENGTH_UM",
    "THICK_FILAMENT_LENGTH_UM",
    "THIN_FILAMENT_LENGTH_UM",
    "H_ZONE_LENGTH_UM",
    "THICK_FILAMENT_DIAMETER_NM",
    "THIN_FILAMENT_DIAMETER_NM",
    "INTER_FILAMENT_SPACING_NM",
    "MYOSIN_HEADS_PER_THICK",
    "CROSS_BRIDGE_AXIAL_PERIOD_NM",
    "HEADS_PER_HELICAL_TURN",
    "HILL_A_OVER_FMAX",
    "HILL_B_OVER_VMAX",
    "LENGTH_TENSION_L_MIN_UM",
    "LENGTH_TENSION_L_PLATEAU_LOW_UM",
    "LENGTH_TENSION_L_PLATEAU_HIGH_UM",
    "LENGTH_TENSION_L_MAX_UM",
    "POWER_STROKE_NM",
    "SINGLE_HEAD_FORCE_PN",
    "ATP_HYDROLYSIS_KCAL_PER_MOL",
    "ATP_HYDROLYSIS_PN_NM",
    "DUTY_RATIO_SKELETAL",
    "LAMBDA_QCD_MEV",
    "KT_AT_BODY_T_PN_NM",
]
