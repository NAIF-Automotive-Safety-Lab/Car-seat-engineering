# Path B: Continuum Field Theory Roadmap

> Path B is mostly analytical (derive equations, find soliton solutions, compute numerical predictions). It is not a coded implementation plan in the same sense as Path C. This document is a research roadmap with milestone deliverables, not a TDD task list.

**Goal:** Predict measurable numerical values of the Standard Model from the substrate parameters K (medium stiffness) and ρ (medium effective density), with c derived as the natural wave speed of the substrate. Per spec §2: direct derivation only; no renormalization or perturbative correction loops.

**Pre-requisite (already known from Path C):** Spec §6's 2D orbital cone requires 3D simulation/analysis. Path B should work in 3D from the start. The 45° rule in 3D specifies a *cone* of allowed directions (U(1) family per axis), giving the freedom needed for true orbital motion.

**Status:** Not yet started. This document is the plan.

---

## Phase 0: Pre-work

### 0.1 Spec revision

Update `docs/superpowers/specs/2026-04-29-stiff-medium-theory-design.md` to absorb the Path C findings:

- §5 (displacement rule, vectors preserved) is **validated** as a binding mechanism in 2D — note this explicitly.
- §6 (2D bound orbit sweeping a 3D cone) is **partially validated**: a 1D-relative bound state is achieved with linear COM drift; full 2D circulation requires 3D dimensionality and is the target of Path B.
- Add a short section noting that the 45° constraint, when applied in 2D, has only 4 discrete velocity directions and forces 1D bound states under zero-net-momentum. In 3D, the constraint specifies a U(1) cone and admits 2D orbital motion.

**Deliverable:** revised spec doc, committed.

### 0.2 Choice of continuum model

Pick a specific nonlinear elastic continuum theory whose linearization gives wave speed c and whose nonlinear sector admits soliton solutions. Candidates:

1. **Nonlinear sine-Gordon-like field theory** in 3D, with the substrate's "slope" represented as a phase field. Slopes (kink/antikink) are localized; their bound states have known math.
2. **Skyrme-type model** — the medium's strain is parameterized by a unit vector field; topological solitons (skyrmions) carry conserved winding number; well-studied for nucleon physics.
3. **Custom elastic strain field** — a stress-strain tensor formulation tailored to "stiff isotropic 3D medium" with explicit K and ρ. More work to set up, but most directly aligned with the spec's language.

**Decision criterion:** Pick the model whose neutrino-soliton solution most cleanly matches spec §5 (1D propagating strain pulse with ± slope and 45° propagation). The Skyrme route has the cleanest topology story but the least direct connection to "stiff medium."

**Deliverable:** a brief decision document (`docs/superpowers/specs/2026-XX-XX-path-b-model-choice.md`) naming the chosen model and its linearization.

---

## Phase 1: Substrate dynamics

### 1.1 Linear regime: derive c from K and ρ

Starting from the chosen field theory, write the small-amplitude wave equation. For an elastic continuum with stiffness modulus K and density ρ, the wave speed is c = √(K/ρ).

**Deliverable:** an equation showing c² = K/ρ, with the model's specific K and ρ identified.

### 1.2 Nonlinear regime: soliton ansatz

Find a 1D propagating localized solution that carries a ± slope along its length. This is the *neutrino* in the model. Confirm:

- Propagation speed = c (the natural wave speed).
- Direction at 45° to a chosen axis (in 3D, on a cone around the axis).
- Localized: amplitude → 0 at infinity.
- Carries a finite energy E_ν.

**Deliverable:** the soliton solution as a closed-form expression (or numerical solution if not analytical), and its energy.

---

## Phase 2: Bound states (the electron)

### 2.1 Two-soliton bound state

Find a stationary or slowly-varying solution where two neutrino-solitons orbit a common center. Spec §6 says this is the electron. In the chosen model:

- Look for a two-soliton bound state with zero net momentum (rest electron).
- Confirm the bound state's geometry is a 2D rotation in 3D (a ring, sweeping the 3D cone of velocity directions over time).
- Compute the binding energy E_e (= rest mass × c²).

**Deliverable:** the bound-state solution and its energy.

### 2.2 The first numerical checkpoint: electron rest mass

Per spec §2, this is the moment of truth. Compute m_e = E_e / c² from the bound-state energy. Compare to the measured value (511 keV).

If the prediction matches measured to ~1% accuracy: **the model is in business**. The next phases are worth doing.

If the prediction is off by orders of magnitude with no obvious dimensional issue: **the model needs revision**. Per §2, no correction loops — we go back to Phase 0.2 and pick a different continuum model, or revise the bound-state ansatz, until either we get the right number or conclude the spec needs structural change.

**Deliverable:** numerical value m_e^model, comparison to 511 keV, and a verdict (match/mismatch).

---

## Phase 3: Higher-energy structures

If Phase 2.2 succeeds:

### 3.1 Lepton mass spectrum

Use the spec §6 picture (electron + stress-quanta on the vertex) to find the muon and tau:

- One quantum of vertex stress added → muon. Compute its mass; compare to 105.66 MeV (~207× electron).
- Two quanta → tau. Compute its mass; compare to 1776.86 MeV (~3477× electron).
- Verify the closure condition forbids a fourth quantum.

**Deliverable:** lepton mass spectrum from the model, comparison to measured.

### 3.2 Bi-pyramid → nucleon

Find a four-soliton or topological bound state corresponding to spec §7 (bi-pyramidal nucleon). Compute:

- Vertex count (gives quark count).
- Fractional slope per vertex (gives quark charge fractions).
- Total binding energy = nucleon rest mass.
- Compare to proton (938 MeV) and neutron (939.6 MeV).

**Deliverable:** nucleon mass and quark charge fractions from the model.

### 3.3 Hydrogen Rydberg

For the locked e-p pair (spec §7.1), find the excitation modes (rocking, breathing, twisting) and their energies. Compare the spectrum to Rydberg's 1/n² law and the Rydberg constant (13.6 eV).

**Deliverable:** hydrogen spectrum from the model.

### 3.4 Fine-structure constant

If the medium has a natural dimensionless ratio (e.g., a coupling between the strain field and the photon-equivalent oscillation), it should equal α ≈ 1/137.036. Compute and compare.

**Deliverable:** α from the model.

---

## Phase 4: Tests against established physics

Run the model against high-precision tests:

- Inertial mass = gravitational mass to measured precision (already structurally guaranteed by spec §8).
- Pair production / annihilation thresholds (γ → e⁺e⁻ requires γ energy ≥ 2 m_e c²).
- β-decay endpoint energy (n → p + e + ν̄ has a known endpoint at 0.782 MeV).

**Deliverable:** a comparison table of model predictions vs. measurements.

---

## Phase 5: Predictions beyond the SM

If Phases 1–4 succeed, the model has earned the right to predict beyond what's measured:

- A new stable particle from a yet-uncatalogued geometric closure?
- A specific deviation from QED at some energy scale?
- A relation between coupling constants that the SM treats as independent?

**Deliverable:** at least one falsifiable prediction the SM does not make, with a measurable signature.

---

## Methodology bar (per spec §2)

- No renormalization. Loop corrections are not allowed to close gaps to measurement.
- No fitting of free parameters to data. K, ρ, and the model's structure are fixed once at the start.
- If a derived number doesn't match measurement, the response is: revise the substrate model or the closure rules, *not* introduce a correction term.
- Honest reporting: if predictions disagree, the spec is wrong; we update it and try again, or accept falsification.

---

## Estimated effort

This is research-grade theoretical physics. Conservative estimate, even for someone with a strong background:

- **Phase 0:** days (model choice + spec revision).
- **Phase 1:** weeks (soliton ansatz; finding closed-form solutions in nonlinear field theories is hard).
- **Phase 2:** weeks to months (two-soliton bound states are notoriously difficult; even getting a numerical solution is non-trivial).
- **Phase 2.2 (electron mass checkpoint):** is the gate. If it fails, may require returning to Phase 0.2 multiple times.
- **Phase 3+:** months to years if Phase 2 succeeds.

This is a long-horizon program. Path C v1 + v2 took us hours, end-to-end; Path B's first checkpoint is much further out.

---

## What this document is for

This is the *roadmap*, not the *plan*. A real plan for Phase 0 + 1 (model choice + linearization + soliton) would be the next deliverable, ideally written in a fresh session with full context budget. It would specify the exact field equations, the ansatz to try, the math tools to use (Mathematica, sympy, NumPy), and the success criteria.

The user should consider:

1. **Whether to do Phase 0.1 (spec revision) immediately** — small, high-value, captures Path C learnings.
2. **Whether Phase 0.2 + 1 is the right next chunk** — requires a session dedicated to math, with focused tools (computer algebra preferred).
3. **Whether to first do "Path C 3D"** — upgrade the simulation to 3D and re-test for 2D orbits. This is intermediate effort (a few hours, similar shape to Path C v1 + v2) and would give us empirical confidence in the 3D orbital picture before committing to analytical work.

The user has flexibility on ordering. The key constraint is the §2 methodology: at every step, predictions must be derivable directly from the substrate, with no parameter tuning.
