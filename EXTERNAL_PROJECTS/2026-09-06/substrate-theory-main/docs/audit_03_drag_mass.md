# Audit 03 — Drag-as-Mass-Generation Mechanism

**Scope:** the 6 modules that implement, or claim to use, the drag-as-mass mechanism.
**Question:** do they share a single, consistent formula `m = f(gamma, K, rho, xi)`?
**Verdict (preview):** the *anchor-level identity* `m c^2 = ℏ ω_b` is universal. The
*formula relating ω_b to (γ, K, ρ, ξ)* is **not** universal — at least three
distinct functional forms appear, two of which are mutually contradictory in
their γ-dependence. Each module is internally calibrated to the electron
anchor, so all six reproduce m_e ≈ 0.511 MeV; that masks the disagreement.
The mechanism is **rigorously specified per module but ad-hoc across modules.**

---

## 1. Per-module formulas

| # | Module | ω_b formula | γ-dependence at fixed (K, ρ, ξ) |
|---|---|---|---|
| A | `drag_mass_generator.py` (l. 198–213) | `ω_b = (c/ξ) · √(1 + α_drag · γ̃)` ;  γ̃ = γξ/√(Kρ) | √(1+αγ̃) — sub-linear |
| B | `mass_torque_engine.py` (l. 222–247) | `ω_b = γ · ξ · ρ · √K · f_int` | linear in γ |
| C | `primitive_anchoring.py` (l. 130–173) | `ω_b ≡ γ` (the *identity* `m c² = ℏ γ`) | γ IS ω_b |
| D | `lattice_substrate_2d.py` (l. 360–416) | `ω_b² = ω_0² + (γ/2ρ)²`,  ω_0 = √(K/ρ)/ξ | √(ω_0²+γ²/4ρ²) — quadratic-in-γ correction on top of bare ω_0 |
| E | `bound_state_spectrum.py` (l. 60–77) | `ω_b = γ · √(K/ρ) / ξ` | linear in γ (with prefactor c_s/ξ, not ξρ√K) |
| F | `generation_map.py` (l. 246–316) | `m = γ_gen · f(particle)`; γ_gen carries MeV units, with `γ_2/γ_1 = exp(n_M/(K_pair⁴π))` | γ literally IS the mass; mass-ratio relation only |

Across all six modules the **closing identity** is the same:

    m c² = ℏ ω_b                                    (rest-mass = bouncing action/cycle)

What differs is the **dynamical content** that fixes ω_b in terms of the
substrate primitives.

---

## 2. Pairwise consistency analysis

### 2a. A vs. B — the central contradiction

- **A (drag_mass_generator):** `ω_b = (c/ξ) · √(1 + α γξ/√(Kρ))`.
  At γ → 0 the prediction is the **bare Compton frequency** `c/ξ`; drag is a
  sub-leading dressing. The "Compton oscillator without drag" is the zero-th
  order of the mass.
- **B (mass_torque_engine):** `ω_b = γ ξ ρ √K · f_int`.
  At γ → 0 the prediction is **zero mass**. The Compton oscillator does not
  exist without drag; the entire ω_b scale is multiplicatively proportional
  to γ.

These are not the same formula, not even up to a calibration. They disagree
on the most basic question — *does mass survive at γ = 0?*
- A says **yes** (bare Compton).
- B says **no** (mass IS the drag).

This is precisely the discrepancy flagged by the user. Both forms are
*self-consistent* and both produce 0.511 MeV at the electron anchor (B
absorbs everything into the calibration constant `m_e_natural = Λ_QCD/m_e`),
but they encode different physics. **At least one must be wrong.**

### 2b. C vs. all others

C is the *limit case* of any drag-driven formula in which ω_b ≡ γ. It can
be obtained from B by setting `ξρ√K · f_int = 1` (the unit-baseline
primitives), and from A by going to the *strong-drag* limit where the
dressing factor dominates. C is therefore not in conflict with A or B; it
is a calibration choice. But it is not in conflict only because it has no
independent dynamical content — `m c² = ℏγ` is a definition of γ, not a
prediction.

### 2c. B vs. E — the "linear-in-γ" cousins

- **B:** `ω_b = γ · (ξ ρ √K) · f_int`. Combination has dimensions
  `[γ][ξ][ρ][K]^½` — depends on **ξ ρ √K**.
- **E:** `ω_b = γ · √(K/ρ)/ξ = γ · c_s/ξ`. Depends on **c_s/ξ = √(K/ρ)/ξ**.

Both are linear in γ, but the prefactor differs by a factor of
`(ξρ√K) / (√(K/ρ)/ξ) = ξ²ρ √(Kρ)/√K = ξ²ρ²` (after algebra).
For the canonical baseline `K = ρ = 1` and ξ in natural units this factor
is `ξ²`, which is **not unity** (xi_e ≈ 3.86×10⁻¹³ m). So B and E disagree
quantitatively whenever ξ ≠ 1, despite both being "linear in γ".

E is the form that appears in the original `cone_bouncing.py` and in the
docstring of `drag_mass_generator.py`'s introduction: `ω_b = c/ξ · f(γ̃)`
with f → 1 at zero drag. E with γ chosen to be the electron Compton angular
frequency reproduces the same number A produces at γ=0 (because E is then
ω_b = ω_e = c/ξ_e). So E and A are *compatible at the electron anchor* but
encode different functional dependencies on γ away from it.

### 2d. D — damped-oscillator shift

D's formula `ω_b² = ω_0² + (γ/2ρ)²` is the **standard damped harmonic
oscillator frequency shift**, applied here as a *correction* to a separately
measured (FFT) bare frequency. It is not derived from substrate primitives;
it is a phenomenological loading added so the FFT-measured frequency rises
monotonically with γ in the simulator output (the code comment is explicit
about this: "Add a small drag-dependent loading so ω_b grows monotonically
with γ even when the linear-KG bare frequency dominates the spectrum").

This formula is consistent with **none** of A, B, E. It is a third, distinct
parametric form that lives only in the lattice simulator. It happens to
agree with A in the small-γ expansion only after relabeling
(`(γ/2ρ)² ↔ α γ̃ ω_0²` requires `α γ̃ = γ²/(4ρ²ω_0²) = γ²ξ²/(4 Kρ)`, which
is *quadratic* in γ̃, not *linear* as A specifies). So D and A disagree at
leading order in γ̃: D contributes O(γ̃²), A contributes O(γ̃).

### 2e. F — generation map

F doesn't touch (K, ρ, ξ) at all. It uses `γ_gen` as a stand-in for *mass at
generation g*, and the mass ratios come from the lepton-tower exponentials
`exp(n_M/(K_pair⁴ π))`. The only drag-mechanism *claim* it makes is that
the ratios `γ_2/γ_1` and `γ_3/γ_2` ARE the mass ratios. Numerically this
gives:

    n_M/(K_pair⁴·π) = 268/(16π) = 5.33169
    exp(5.33169)    = 206.787   (PDG m_μ/m_e = 206.768; 9.2×10⁻⁵ relative)

So the **generation scaling is correct in F**, but it is decoupled from any
substrate-primitive formula. F could be removed and re-implemented on top
of any of A/B/C/D/E without loss.

---

## 3. Anchor verification (each module gives m_e = 0.511 MeV)

Numerical check at the electron Compton anchor:

| Module | Formula reduced at anchor | m_e (MeV) | Status |
|---|---|---|---|
| A | γ=0 → ω_b = c/ξ_e | 0.51099895 | ✓ identity |
| B | gamma=K=rho=xi=1, T = 1/(Λ_QCD/m_e_obs) | 0.510999 | ✓ pinned |
| C | γ ≡ m_e c²/ℏ → m = ℏγ/c² | 0.51099895 | ✓ tautology |
| D | bare KG mode at ω_0 = c/ξ_e (γ=0 limit) | 0.51099895 | ✓ identity |
| E | γ = ω_e ξ_e/c_s → ω_b = ω_e | 0.51099895 | ✓ anchor by construction |
| F | γ_1 := m_e (in MeV-equiv units) | 0.51099895 | ✓ definitional |

**All six pass — but each by definition / construction, not as a prediction.**
The electron mass is the *anchor* in every module; no module *predicts* it
from independent substrate parameters. This is the correct B3 calibration
strategy, but it means anchor agreement is not evidence of formula consistency.

---

## 4. The dimensionless ratio m c² / (ℏ γ)

This ratio should be O(1) and ideally identical across modules.

| Module | m c²/(ℏγ) at anchor | Notes |
|---|---|---|
| A | special-cased to 1 when γ→0 (code returns 1.0) | tautological |
| B | m_e/Λ_QCD ≈ 2.55 × 10⁻³ | **NOT O(1)** — γ is in different units |
| C | exactly 1 | identity |
| D | depends on FFT measurement; γ=0 → diverges | undefined at anchor |
| E | exactly 1 at electron anchor | by construction |
| F | exactly 1 (γ_1 := m_e in MeV) | definitional |

**B is the outlier.** Its γ is the *substrate-natural* drag (dimensionless, set
to 1.0 by primitives), not an angular frequency. The dimensionless ratio
m c²/(ℏγ) in B's units is therefore the Compton frequency *in units of Λ_QCD*,
which is m_e/Λ_QCD ~ 0.0026 for the electron — not O(1).

This means the "γ" in module B is **not the same physical quantity** as the γ
in modules A, C, D, E. A,C,D,E use γ with dimensions of [angular frequency]
(SI Hz). B uses γ as a dimensionless substrate primitive that needs an
explicit Λ_QCD prefactor to recover MeV mass. The two conventions are
related by `γ_B = γ_A·ξρ√K · 1/Λ_QCD` after substituting f_int=1 — i.e. by a
nontrivial dimensional reshuffle that no module documents.

---

## 5. Generation scaling check

The generation lift `γ_g+1/γ_g = exp(n_M/(K_pair⁴ π))` should produce the
PDG mass ratios when the same particle type is anchored at gen 1.

    γ_2/γ_1 = exp(268/(16π)) = exp(5.33169) = 206.787
    PDG m_μ/m_e                              = 206.768
    relative error                            = 9.2 × 10⁻⁵

✓ Works in F where it is explicit.
✓ Works in B's `_t_muon` (uses the same exp(n_M/(K_pair⁴ π))).
✗ Not implemented in A, C, D, E. A's `electron_mass`/`quark_mass_spectrum`
   simply re-anchors ξ to each particle's PDG Compton wavelength rather than
   propagating γ through the generation lift. This is the **opposite**
   philosophy from F: A says "different particle = different ξ"; F says
   "different generation = different γ". Both reproduce the PDG values, but
   they put the physical content in different primitives.

---

## 6. Honest assessment

**Is the drag-as-mass mechanism rigorously specified?**

- **Within each module: yes.** Each of the six modules specifies a definite
  formula, calibrates it at the electron anchor, and produces consistent
  internal predictions.
- **Across modules: no.** Three distinct formulas for ω_b coexist:
  1. **A's "Compton + drag dressing":** ω_b = (c/ξ)√(1+αγ̃) — drag is
     sub-leading dressing on a bare Compton oscillator.
  2. **B's "mass IS drag":** ω_b ∝ γ — no mass at γ=0.
  3. **D's "damped oscillator shift":** ω_b² = ω_0² + (γ/2ρ)² — quadratic
     correction.

  E is a special case of the linear family; C and F are anchor-tautologies
  that are compatible with anything.

**Which formula is correct?**

The user's MEMORY note "mass = cumulative torque on substrate" and the
B3 mass-torque axiom `m = Λ_QCD · T(config)` favor B's interpretation —
mass IS drag-derived torque, no Compton oscillator without drag. But A's
docstring is *explicit* that "for zero drag (gamma=0) the bound state is a
pure Compton oscillator (omega_b = c/xi)" — the opposite. So the framework
itself contains the contradiction, not just the implementation.

**Recommendation for unification:**

The cleanest reconciliation is to elevate **C's identity m c² = ℏγ** as the
*definition* and demand all other modules express their ω_b such that
γ ≡ ω_b in the appropriate physical units. Then:

- A's formula becomes a *prediction*: `γ = (c/ξ)√(1+αγ̃)` is an implicit
  equation for γ in terms of (K, ρ, ξ), reducing to γ = c/ξ at α→0. This is
  consistent with C and is the cleanest "drag from primitives" derivation.
- B's `ω_b = γ ξ ρ √K · f_int` should be reinterpreted as a *definition of
  γ in substrate-natural units* (γ_B := ω_b · 1/(ξρ√K · f_int)), so the
  mass formula is `m c² = ℏω_b = ℏ · γ_B · ξρ√K · f_int`. This is what the
  code computes; the docstring just needs to say "γ_B is dimensionless,
  not a Hz drag rate."
- D's damped-oscillator correction should be dropped or relabeled as a
  numerical artifact of the FFT-peak measurement, not as a physics formula.
- E and F are already compatible with C.

Concretely: a single new `drag_mass_protocol.py` that exports
`omega_b_from_primitives(K, rho, xi, gamma) -> float` and is called by all
six modules would eliminate the contradiction. Today each module
re-implements its own.

**Bottom line:** the drag-as-mass *concept* is universal; the *formula* is
not. The framework needs a single canonical ω_b(K,ρ,ξ,γ) statement before
"drag-as-mass" can be called rigorously specified. Until then, calibrating
to the electron anchor in each module is what makes everything *appear* to
agree — and is also what hides the disagreement from external observers.
