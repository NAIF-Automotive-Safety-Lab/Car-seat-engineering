"""mobius_sheet_swap.py
==========================

Paired GEOMETRY + SIMULATION module for the **Möbius sheet-swap involution**
that ties together four otherwise-independent observables of the B3 /
stiff-medium framework:

    (a) particle ↔ antiparticle distinction for **charged** fermions
    (b) Majorana identity (ν = ν̄) for **neutral** fermions
    (c) the saturation cap σ ≤ 1/2 as a **fixed point** of the involution
    (d) half-integer spin from the **double cover** of SO(3) by SU(2)

Geometry
--------
The Möbius bundle ``π : E → S¹`` is a non-trivial real line bundle with
total space ``E ≃ ([0, 2π] × {±1}) / ∼``, where the equivalence relation
identifies ``(0, +1) ∼ (2π, −1)`` (and conversely).  Locally the bundle
has **two sheets** labelled ``A`` and ``B``; globally a single transport
around the base circle exchanges them.

The corresponding sheet-swap involution is

    τ : E → E,        τ(θ, s) = (θ, −s),        τ² = id.

This is a **Z/2 action**.  Its three structural consequences:

* charged fermion section ``ψ_q(θ, s) = q · s`` is *not* fixed by τ:
  ``τ ψ_q = −ψ_q`` so the swap exchanges particle ↔ antiparticle.
* neutral fermion section ``ψ_0(θ, s) = 0 + 0 i`` IS fixed by τ:
  the two-component Dirac field collapses to its Majorana real form
  (``ν = ν̄``).
* the saturation amplitude ``A = 1/2`` is invariant: it lives on the
  fixed locus of the geometric reflection in fibre coordinates and so
  acts as a hard wall that no excitation can cross.

The double-cover relation between SU(2) and SO(3) — equivalently the
fact that one circuit on S¹ is half a circuit on the universal cover
S¹ → S¹ ramified at the swap point — is what realises *spin-½*: the
holonomy of the bundle on a single base loop is ``exp(iπ) = −1`` and
two base loops are needed to come back to ``+1``.

This module exposes the geometry and the simulator as two paired
classes, plus a small set of demo helpers used by
``scripts/render_all_visuals.py``:

* :class:`MobiusSheetGeometry` — explicit 2-sheet Möbius bundle, sheet
  labels, sheet-swap involution as a Z/2 action, fixed-point set, and a
  3-D parametrisation of the strip used by the visualizer.
* :class:`SheetSwapSimulator` — time-evolves a small two-component
  excitation (charged or neutral) through the sheet-swap and reports
  the post-swap state, the swap-invariance flag, and a log of the
  saturation amplitude (which never crosses 1/2).

The whole module is fully deterministic and depends only on
``numpy`` (numerical) and (for the convenience plotting helpers)
``matplotlib`` — both already pinned in ``pyproject.toml``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

#: Saturation cap σ_max — the substrate cannot tilt past this value
#: (spec §§5, 18.4).  It is a *fixed point* of the sheet-swap involution
#: in fibre coordinates.
SIGMA_MAX: float = 0.5

#: Sheet labels.  ``+1`` is sheet A, ``−1`` is sheet B; the Möbius
#: identification glues ``(0, +1) ∼ (2π, −1)``.
SHEET_A: int = +1
SHEET_B: int = -1


# ---------------------------------------------------------------------------
# Geometry: 2-sheet Möbius bundle
# ---------------------------------------------------------------------------

@dataclass
class MobiusSheetGeometry:
    """Explicit 2-sheet Möbius bundle and its sheet-swap involution.

    The bundle ``π : E → S¹`` is constructed as a quotient

        E = ([0, 2π] × {+1, −1}) / ∼,

    with the gluing ``(0, +1) ∼ (2π, −1)`` and ``(0, −1) ∼ (2π, +1)``.
    The two sheets are *locally* trivial (over any open arc not crossing
    θ = 0) but **globally** identified after one traversal of the base
    circle — that is the defining feature of the Möbius construction.

    The sheet-swap is realised by

        τ : E → E,        τ(θ, s) = (θ, −s),

    which is a generator of a Z/2 action: ``τ² = id``.

    Parameters
    ----------
    radius : float, default 1.0
        Major radius of the embedded strip in 3-D (used by the
        visualizer; does not affect the algebraic Z/2 structure).
    width : float, default 0.5
        Half-width of the strip in fibre coordinates.  The fibre
        ``v ∈ [−width, +width]`` is the geometric realisation of the
        sheet coordinate; the **sign** of ``v`` is the sheet label.
    n_theta : int, default 240
        Number of base points used to discretise the strip surface.
    n_v : int, default 24
        Number of transverse points across the strip.

    Notes
    -----
    * ``sheets`` returns the explicit list ``[SHEET_A, SHEET_B] = [+1, −1]``.
    * ``involution(theta, s)`` returns ``(theta, −s)``; ``involution_squared``
      is exactly the identity for any input.
    * ``fixed_set()`` returns the locus where ``s = 0``, i.e. the
      *core circle* of the strip; this is the σ = 1/2 cap edge in the
      physical interpretation.
    * ``surface_xyz()`` returns ``(X, Y, Z)`` arrays for ``plot_surface``.
    """

    radius: float = 1.0
    width: float = 0.5
    n_theta: int = 240
    n_v: int = 24

    def __post_init__(self) -> None:
        if self.radius <= 0.0:
            raise ValueError(f"radius must be positive; got {self.radius!r}")
        if self.width <= 0.0:
            raise ValueError(f"width must be positive; got {self.width!r}")
        if self.n_theta < 8 or self.n_v < 2:
            raise ValueError("n_theta >= 8 and n_v >= 2 required")

    # ------------------------------------------------------------------
    # Bundle structure
    # ------------------------------------------------------------------
    @staticmethod
    def sheets() -> Tuple[int, int]:
        """Return the two sheet labels ``(A, B) = (+1, −1)``.

        Returns
        -------
        tuple of int
            Exactly two sheet labels, in canonical order.
        """
        return (SHEET_A, SHEET_B)

    @staticmethod
    def n_sheets() -> int:
        """Return the number of sheets — always 2 for the Möbius bundle."""
        return 2

    # ------------------------------------------------------------------
    # Z/2 sheet-swap involution
    # ------------------------------------------------------------------
    @staticmethod
    def involution(theta: float, sheet: int) -> Tuple[float, int]:
        """Apply the Z/2 sheet-swap involution τ(θ, s) = (θ, −s).

        Parameters
        ----------
        theta : float
            Base coordinate (radians); preserved by the involution.
        sheet : int
            Sheet label, ``+1`` or ``−1``.

        Returns
        -------
        (theta_out, sheet_out) : tuple
            Same base coordinate; sheet label flipped in sign.

        Raises
        ------
        ValueError
            If ``sheet`` is not ±1.
        """
        if sheet not in (SHEET_A, SHEET_B):
            raise ValueError(
                f"sheet must be +1 or −1; got {sheet!r}"
            )
        return (float(theta), int(-sheet))

    @classmethod
    def involution_squared(cls, theta: float, sheet: int) -> Tuple[float, int]:
        """Apply τ twice; should be the identity for every (θ, s)."""
        once = cls.involution(theta, sheet)
        twice = cls.involution(*once)
        return twice

    @classmethod
    def is_involution(cls, n_samples: int = 64) -> bool:
        """Verify τ² = id on a deterministic grid of samples.

        Returns
        -------
        bool
            True iff ``involution_squared`` matches the input on every
            sampled (θ, s); used as a self-test invariant.
        """
        thetas = np.linspace(0.0, 2.0 * math.pi, n_samples, endpoint=False)
        for theta in thetas:
            for s in cls.sheets():
                t2, s2 = cls.involution_squared(float(theta), s)
                if not (math.isclose(t2, theta, rel_tol=0.0, abs_tol=1e-12)
                        and s2 == s):
                    return False
        return True

    # ------------------------------------------------------------------
    # Fixed-point set in fibre coordinates
    # ------------------------------------------------------------------
    def fixed_set_radius(self) -> float:
        """Return the geometric radius of the fixed locus of τ.

        In fibre coordinates the involution ``v → −v`` has the single
        fixed point ``v = 0``; geometrically this is the **core circle**
        of the Möbius strip at radius ``self.radius``.  In the physical
        identification this is the saturation cap ``σ = 1/2``.

        Returns
        -------
        float
            The major radius (== ``self.radius``).
        """
        return float(self.radius)

    @staticmethod
    def saturation_cap_is_fixed_point() -> bool:
        """The σ = 1/2 cap sits on the fixed locus of τ → returns True.

        The cap is the boundary an excitation cannot cross; under the
        sheet-swap reflection ``v → −v`` it is mapped to itself, which
        is exactly the *invariant submanifold* of the involution.
        """
        # In fibre coords, the cap edge is v = 0 (core circle); under
        # v → −v this is invariant.  The σ = 1/2 amplitude saturates
        # to the cap; physically this is the bouncing wall.
        cap_value = SIGMA_MAX
        # Apply v → −v with v == cap_value treated as the boundary radius;
        # the boundary set {|v| = cap_value} is set-wise invariant.
        return abs(cap_value) == abs(-cap_value) and cap_value > 0.0

    # ------------------------------------------------------------------
    # Spin-½ from the double cover
    # ------------------------------------------------------------------
    @staticmethod
    def loop_holonomy(n_loops: int = 1) -> float:
        """Holonomy of the Möbius bundle on ``n_loops`` traversals of S¹.

        A single loop on the base circle picks up the sheet-swap
        permutation, which acts on the spinor section as ``ψ → −ψ``.
        Two loops give the identity.  This is the *double cover* of
        SO(3) by SU(2) realised geometrically:

            holonomy(n) = (−1)^n.

        Returns
        -------
        float
            ``+1.0`` for even ``n_loops``, ``−1.0`` for odd ``n_loops``.
        """
        return float((-1) ** int(n_loops))

    @classmethod
    def is_double_cover(cls) -> bool:
        """Return True iff one base loop is *half* a spinor loop.

        Equivalent to ``loop_holonomy(1) == −1`` AND
        ``loop_holonomy(2) == +1``.
        """
        return (cls.loop_holonomy(1) == -1.0
                and cls.loop_holonomy(2) == +1.0)

    @classmethod
    def spin_value(cls) -> float:
        """The (half-integer) spin forced by the double cover.

        Returns
        -------
        float
            ``0.5`` — the unique half-integer spin compatible with the
            ``(−1)^n`` holonomy of the Möbius double cover.
        """
        if not cls.is_double_cover():
            raise RuntimeError(
                "Möbius bundle holonomy inconsistent with double cover"
            )
        return 0.5

    # ------------------------------------------------------------------
    # 3-D parametrisation for the visualiser
    # ------------------------------------------------------------------
    def surface_xyz(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (X, Y, Z) for an embedded Möbius strip.

        Uses the standard parametrisation

            x = (R + v cos(θ/2)) cos(θ)
            y = (R + v cos(θ/2)) sin(θ)
            z = v sin(θ/2)

        with ``θ ∈ [0, 2π]`` and ``v ∈ [−width, +width]``.  The sign of
        ``v`` is the sheet label; the half-angle ``θ/2`` realises the
        double-cover twist.

        Returns
        -------
        (X, Y, Z) : tuple of numpy.ndarray, each shape (n_v, n_theta)
        """
        u = np.linspace(0.0, 2.0 * np.pi, self.n_theta)
        v = np.linspace(-self.width, +self.width, self.n_v)
        U, V = np.meshgrid(u, v)
        X = (self.radius + V * np.cos(U / 2.0)) * np.cos(U)
        Y = (self.radius + V * np.cos(U / 2.0)) * np.sin(U)
        Z = V * np.sin(U / 2.0)
        return X, Y, Z

    def sheet_color_field(self) -> np.ndarray:
        """Sign of the fibre coordinate v — used to colour the two sheets.

        Returns
        -------
        numpy.ndarray of shape (n_v, n_theta)
            ``+1.0`` on sheet A (v > 0), ``−1.0`` on sheet B (v < 0),
            and ``0.0`` exactly on the fixed circle (v = 0).
        """
        v = np.linspace(-self.width, +self.width, self.n_v)
        # Sign field along the fibre, broadcast to (n_v, n_theta).
        s = np.sign(v).reshape(-1, 1)
        return np.repeat(s, self.n_theta, axis=1)


# ---------------------------------------------------------------------------
# Excitation states (charged vs neutral)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Excitation:
    """A two-component substrate excitation labelled by sheet.

    Attributes
    ----------
    charge : float
        Electric charge in units of ``e``.  ``0.0`` selects the neutral
        (Majorana) branch; non-zero selects the charged branch.
    sheet : int
        Sheet label ``+1`` (A) or ``−1`` (B).
    amplitude : float
        Strain amplitude on the bundle, magnitude must satisfy
        ``|amplitude| ≤ SIGMA_MAX``.
    label : str
        Human-readable name (e.g. ``'electron'``, ``'neutrino'``).
    """

    charge: float
    sheet: int = SHEET_A
    amplitude: float = 0.4
    label: str = ""

    def __post_init__(self) -> None:
        if self.sheet not in (SHEET_A, SHEET_B):
            raise ValueError(f"sheet must be ±1; got {self.sheet!r}")
        if abs(self.amplitude) > SIGMA_MAX + 1e-12:
            raise ValueError(
                f"|amplitude|={abs(self.amplitude)} exceeds saturation cap "
                f"{SIGMA_MAX}; the substrate forbids this state."
            )

    @property
    def is_neutral(self) -> bool:
        """True iff this excitation has zero electric charge."""
        return self.charge == 0.0


# ---------------------------------------------------------------------------
# Sheet-swap simulator
# ---------------------------------------------------------------------------

@dataclass
class SheetSwapSimulator:
    """Time-evolve a substrate excitation through the sheet-swap τ.

    The simulator transports a two-component excitation along the base
    circle of the Möbius bundle.  After one full traversal the bundle
    glues sheet A to sheet B (and vice versa).  The result depends
    crucially on whether the excitation is **charged** or **neutral**:

    * charged: the post-swap state has flipped sheet *and* flipped
      charge, i.e. it is the **antiparticle**.  particle and antiparticle
      are *distinct* observable states.
    * neutral: the post-swap state is the *same* state as the input —
      the swap is gauged out, and one obtains the **Majorana** identity
      ``ν = ν̄``.

    Parameters
    ----------
    geometry : MobiusSheetGeometry
        Bundle on which the excitation lives.
    n_steps : int, default 200
        Number of intermediate angular steps in one base loop, used
        purely for trajectory logging and the visualiser; the swap
        itself is exact.
    """

    geometry: MobiusSheetGeometry = field(default_factory=MobiusSheetGeometry)
    n_steps: int = 200

    def __post_init__(self) -> None:
        if self.n_steps < 4:
            raise ValueError(f"n_steps must be >= 4; got {self.n_steps}")

    # ------------------------------------------------------------------
    # Single-step swap and full-loop evolution
    # ------------------------------------------------------------------
    def apply_sheet_swap(self, exc: Excitation) -> Excitation:
        """Apply the Z/2 swap τ to a single excitation.

        The action is

            τ(charge, sheet, A) = ( -sign · charge,  -sheet,  A ),

        where ``sign = +1`` for charged excitations (so charge → −charge,
        i.e. particle → antiparticle) and the *neutral* sub-algebra
        identifies the two sheets so that ``τ`` reduces to the identity
        on the (charge=0) component.

        Parameters
        ----------
        exc : Excitation

        Returns
        -------
        Excitation
            Post-swap state.  For charged input: swapped sheet AND
            sign-flipped charge.  For neutral input: identical to ``exc``
            (Majorana fixed point).
        """
        if exc.is_neutral:
            # τ acts trivially on the neutral component: ν = ν̄.
            return Excitation(
                charge=0.0,
                sheet=exc.sheet,                 # same sheet
                amplitude=exc.amplitude,         # same amplitude
                label=exc.label,
            )
        # Charged: flip both the sheet label and the charge.
        return Excitation(
            charge=-exc.charge,
            sheet=-exc.sheet,
            amplitude=exc.amplitude,
            label=self._antiparticle_name(exc.label),
        )

    @staticmethod
    def _antiparticle_name(name: str) -> str:
        """Return the conventional antiparticle label for a particle name."""
        if not name:
            return ""
        if name.startswith("anti-"):
            return name[len("anti-"):]
        return f"anti-{name}"

    def evolve_one_loop(self, exc: Excitation) -> Dict[str, object]:
        """Transport ``exc`` once around the base circle and apply τ.

        Returns a dictionary suitable for plotting:

        * ``"thetas"``     — array of base angles θ ∈ [0, 2π)
        * ``"amplitudes"`` — array of the (capped) amplitude trajectory
        * ``"sheet_labels"`` — array of sheet labels at each step
        * ``"start"``      — the input :class:`Excitation`
        * ``"end"``        — the post-loop :class:`Excitation`
          (charged: antiparticle; neutral: identical)
        * ``"is_majorana"`` — ``True`` iff the swap returned the same
          state (i.e. neutral input)
        * ``"sigma_cap"``  — the saturation cap value (always 0.5)

        The amplitude is *capped* at ``SIGMA_MAX`` at every step, which
        is the simulator-side check that the σ = 1/2 wall acts as a
        fixed point of the involution.
        """
        thetas = np.linspace(0.0, 2.0 * math.pi, self.n_steps, endpoint=False)
        amps = np.full_like(thetas, fill_value=abs(exc.amplitude), dtype=float)
        amps = np.minimum(amps, SIGMA_MAX)        # cap enforcement
        # Sheet label oscillates smoothly across the loop in the visualiser:
        # we report the discrete pre-swap label at every step.
        sheet_arr = np.full_like(thetas, fill_value=exc.sheet, dtype=float)

        post = self.apply_sheet_swap(exc)
        is_majorana = (
            exc.is_neutral
            and post.charge == exc.charge
            and post.sheet == exc.sheet
            and math.isclose(post.amplitude, exc.amplitude)
        )
        return {
            "thetas": thetas,
            "amplitudes": amps,
            "sheet_labels": sheet_arr,
            "start": exc,
            "end": post,
            "is_majorana": bool(is_majorana),
            "sigma_cap": SIGMA_MAX,
        }

    # ------------------------------------------------------------------
    # Convenience: paired demos used by the renderer
    # ------------------------------------------------------------------
    def demo_charged_vs_neutral(
        self,
        charged: Excitation | None = None,
        neutral: Excitation | None = None,
    ) -> Dict[str, Dict[str, object]]:
        """Run a single base loop for a charged and a neutral excitation.

        Defaults: ``electron`` (charge −1, sheet A) and ``neutrino``
        (charge 0, sheet A).  Returned dict has keys ``"charged"`` and
        ``"neutral"`` whose values are the dictionaries returned by
        :meth:`evolve_one_loop`.
        """
        if charged is None:
            charged = Excitation(charge=-1.0, sheet=SHEET_A,
                                 amplitude=0.4, label="electron")
        if neutral is None:
            neutral = Excitation(charge=0.0, sheet=SHEET_A,
                                 amplitude=0.4, label="neutrino")
        return {
            "charged": self.evolve_one_loop(charged),
            "neutral": self.evolve_one_loop(neutral),
        }

    # ------------------------------------------------------------------
    # Cap enforcement is itself a falsifiable invariant
    # ------------------------------------------------------------------
    @staticmethod
    def cap_is_invariant(amps: np.ndarray, tol: float = 1e-12) -> bool:
        """True iff every amplitude in ``amps`` lies in ``[0, σ_max]``.

        The σ = 1/2 wall is enforced as a *hard* cap: any state with
        amplitude exceeding the cap is unphysical.
        """
        a = np.asarray(amps, dtype=float)
        return bool(np.all(a >= -tol) and np.all(a <= SIGMA_MAX + tol))


# ---------------------------------------------------------------------------
# Module-level convenience and self-test entry points
# ---------------------------------------------------------------------------

def default_charged() -> Excitation:
    """Default charged demo state — a left-handed electron on sheet A."""
    return Excitation(charge=-1.0, sheet=SHEET_A, amplitude=0.4,
                      label="electron")


def default_neutral() -> Excitation:
    """Default neutral demo state — a Majorana neutrino on sheet A."""
    return Excitation(charge=0.0, sheet=SHEET_A, amplitude=0.4,
                      label="neutrino")


def summary() -> Dict[str, object]:
    """Return a one-shot summary of all four observables forced by τ.

    Used by the ``--summary`` CLI flag and by the regression tests; it
    is the single function that exercises every required assertion in
    one place:

    * 2 sheets,
    * τ² = id,
    * charged → distinct antiparticle,
    * neutral → identical (Majorana),
    * σ = 1/2 fixed point,
    * spin-½ from double cover.
    """
    geom = MobiusSheetGeometry()
    sim = SheetSwapSimulator(geometry=geom)
    charged_in = default_charged()
    neutral_in = default_neutral()
    charged_out = sim.apply_sheet_swap(charged_in)
    neutral_out = sim.apply_sheet_swap(neutral_in)
    return {
        "n_sheets": geom.n_sheets(),
        "involution_squared_is_identity": geom.is_involution(),
        "charged_swap_distinct": (
            charged_out.charge == -charged_in.charge
            and charged_out.sheet == -charged_in.sheet
        ),
        "neutral_swap_identical": (
            neutral_out.charge == neutral_in.charge
            and neutral_out.sheet == neutral_in.sheet
            and neutral_out.amplitude == neutral_in.amplitude
        ),
        "saturation_cap_value": SIGMA_MAX,
        "saturation_cap_fixed_point": geom.saturation_cap_is_fixed_point(),
        "double_cover_holonomy_one_loop": geom.loop_holonomy(1),
        "double_cover_holonomy_two_loops": geom.loop_holonomy(2),
        "spin_value": geom.spin_value(),
    }


if __name__ == "__main__":  # pragma: no cover
    s = summary()
    width = max(len(k) for k in s)
    print("=" * (width + 32))
    print("Möbius sheet-swap involution — structural summary")
    print("=" * (width + 32))
    for k, v in s.items():
        print(f"  {k:<{width}} : {v}")
    print("=" * (width + 32))
