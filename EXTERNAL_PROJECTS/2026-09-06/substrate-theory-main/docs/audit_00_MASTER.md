# Master Audit — Stiff Medium / B3 Substrate Framework

**Date:** 2026-05-01
**Auditor:** Independent synthesis pass (audit_00)
**Scope:** `MODEL.md`, `RESULTS.md`, `scripts/` (~250 files), `src/stiff_medium/` (~120 modules), `tests/` (60 passing).
**Status of sibling audits:** audit_01 through audit_09 not present at audit time. This master document therefore stands as a primary, not a synthesis. When the sibling audits land they should be folded in via a revision pass.

This audit is deliberately critical-but-fair. It is written assuming a reader who is sympathetic to substrate ontologies but skeptical of any specific numerology, and who wants to know where the framework would actually break.

---

## 1. Overall consistency

### 1.1 What hangs together

The framework has a genuine internal architecture, not a bag of formulas:

- **One ontological commitment** (a 3D stiff elastic medium with four primitives K, ρ, ξ, γ + Möbius half-flux topology) is propagated downward into every sector. Mass = cone-bouncing drag, charge = chirality channel of back-reaction, gravity = symmetric channel of the same back-reaction, spin-½ = Möbius topology, saturation σ ≤ ½ = horizons + cosmology. This is a single mechanism family, not five glued-together mechanisms.
- **The B3 inventory layer** (12 integers + 2 axioms + 4 anchors) is presented as a *derived* structure on top of the substrate primitives, not as a competing ontology. `strain_medium_is_tighter.py` makes the claim explicit: B3's integers are supposed to fall out of substrate packing, Möbius topology, and saturation. The framework is honest that this derivation is partial.
- **A canonical Lagrangian** (§18.45 / §18.84 in MODEL.md §9) is written down end-to-end and unifies the substrate, Dirac, EM bundle, and saturation potential in one expression. Not every observable is computed *from* this Lagrangian (most are computed from the inventory layer), but the existence of the Lagrangian gives the framework a definite mathematical object to attack.

### 1.2 Where it strains

- **Two "engine" languages running in parallel.** The substrate side (K, ρ, ξ, γ, σ-cap, sine-Gordon × saturation potential) and the B3 inventory side (n_M=268, n_R=18, n_A=45, K_pair=2, K_rank=5, etc.) are both used to derive observables. When they overlap (α, lepton ratios, PMNS), the answers agree, which is a real consistency check. But the *map* between them is still partly verbal. The "tighter" claim in `strain_medium_is_tighter.py` is correct in input-counting but several of the bottom-up derivations there are sketches, not proofs (e.g. n_M = 268 from "Möbius bundle index" — the actual integer drop has not been demonstrated from the substrate Lagrangian).
- **Drag coefficient γ does too much work.** Once γ is admitted, it sets electron mass, the Q-factor in the α closed form (Q = (11/12)·n_M = 245.67), the Higgs mass closure, and the entire lepton ladder via cone-bouncing. This is a genuine unification *if* γ is independently determined; it becomes near-tautological if γ is read off α and then re-used. The current scripts read γ off K, ρ, ξ via `mass_torque_engine.py`, which is the right direction but not yet a fully closed loop.
- **The "drag replaces Higgs" claim** is consistent in the new clean-Lagrangian proposal (`clean_lagrangian_proposal.py`) but the encompassing Lagrangian in §9 still carries an explicit Yukawa term `g_Y φ ψ̄ψ`. The framework is mid-transition; documents need to converge.

**Verdict:** Internally consistent at the conceptual level. Mid-transition at the formal-Lagrangian level. The two parallel languages (substrate vs. B3) are not yet a single derivation chain end-to-end, but they are demonstrably compatible everywhere they have been jointly tested.

---

## 2. Strongest results — top 5 most rigorous derivations

These are ranked by *honest derivation depth × experimental match × parameter-freeness*. Curve-fits and integer-recyclings are explicitly excluded.

1. **Gravity / EM force ratio = 8.09 × 10⁻³⁷ (0.06%).** The strongest single number in the framework. Two channels of the same Poisson back-reaction predict the ratio (m_p/M_Pl)² / α structurally. Independently verifiable, no tunable knob, sits on a derivation that uses only the two-channel structure of §2.3. If this number had been *off* by an order of magnitude the framework would be dead.

2. **GW speed = c at 10⁻¹⁵ (GW170817).** Not derived after the fact — it is *forced* by the substrate ontology. Anything propagating in the medium has speed √(K/ρ). EM and gravity are two channels of the same medium; identical c is a structural prediction, not a calculation. This is a genuine ontology-forced result.

3. **CPT (m_p̄ = m_p) at 16 ppt (BASE).** Same status: the substrate has no preferred orientation, so antimatter is the orientation-flipped configuration with identical strain spectrum. Forced, not fit. The fact that BASE saw zero is consistent with a "the universe is the substrate, antimatter is the same substrate run backward" ontology and would have falsified it otherwise.

4. **Mercury perihelion = 42.99″/century, light bending 1.7509″.** These come from the strong-field nonlinear back-reaction of the substrate (§18.39) and reproduce GR at the precision of measurement. They are not as ontology-forced as #1–#3 — they require the nonlinear σ structure to be exactly the right form — but the agreement is at experimental precision and the path through the substrate Lagrangian is explicit.

5. **Closed-form α to 0.004%.** α = (11/(48π³)) · exp(−π/Q), Q = (11/12)·n_M = 245.67. The match is genuinely impressive. *Caveat:* this is partly an integer-drop result (n_M = 268 is not yet derived from the Lagrangian, only motivated by Möbius bundle counting), so it sits between "derivation" and "well-motivated formula." Still belongs in the top 5 because if any of the integer constraints had been wrong the formula would have missed by orders of magnitude.

**Honorable mentions:** Hylleraas-helium at 5×10⁻⁵ (inherited from QED, not a substrate triumph but a consistency check), 5-loop g-2 to 10⁻¹⁰ (same — inheritance), Σm_ν = 60.5 meV vs DESI bound (passes ΛCDM cleanly), Cabibbo angle = 1/(π√2) at 0.035% (clean integer-free formula).

---

## 3. Weakest claims — top 5 most ad-hoc or hand-waved

1. **n_M = 268 has no derivation.** It enters the α formula via Q = (11/12)·n_M and sets the precision of the headline result. It is presented as a B3 inventory integer "from Möbius bundle index over 2-torus" but no script in `scripts/` actually derives 268 from the substrate Lagrangian. This is the single highest-leverage open derivation in the framework.

2. **The 12 B3 integers, collectively, are not yet derived from the substrate.** `strain_medium_is_tighter.py` *claims* the derivations (N_BAM=6 from hexagonal packing, K_pair=2 from Möbius sheets, etc.) but most are one-paragraph sketches. The framework's "6 inputs vs. SM's 30" compression is honest only if those derivations close. They mostly haven't.

3. **Density perturbation amplitude fails by 10¹³.** Substrate-saturation cosmology (per memory note `b3_substrate_saturation_cosmology.md`) resolves horizon, flatness, and singularity but *fails* the perturbation amplitude by 13 orders of magnitude. This is the largest open quantitative gap in the framework, larger than any active falsifier. The papers downplay it; an audit cannot.

4. **PMNS angles as α-multiples (sin²θ_12 = 42α, sin²θ_13 = 3α).** These match cleanly but the integer prefactors (42, 3) are fit to inventory after the fact. The framework calls them "cell-inventory sums" but does not derive *which* sum is required *before* knowing the answer. This is post-hoc integer-matching, not prediction. Same critique applies to many of the "B3 integer recycling" entries in §5.1 of MODEL.md.

5. **The Mexican-hat-vs-sine-Gordon transition is incomplete.** §6.1 of MODEL.md proposes dropping the Higgs quartic in favor of sine-Gordon × saturation, but §9 still writes the encompassing Lagrangian with a Yukawa term. The "cleaner formulation" closes the Higgs mass at 0.23%, but it has not been integrated into the rest of the symbolic chain. Until it is, the claim that the framework eliminates Higgs fine-tuning is provisional.

**Honorable dishonorable mentions:** η_baryon at 24% (called a success in MODEL.md but 24% is large), δ_CP at 1.83% off (fine, but "= -π/2 maximal" is the kind of formula that is right precisely when it is right and gives no further information), Ω_DM/Ω_b = 5.35 (clean formula but the dark-matter cube-cell candidate at 27.5 GeV has zero direct experimental support).

---

## 4. Critical open gaps — what most needs to be closed

Ranked by leverage on the framework's credibility:

1. **Derive n_M = 268 from the substrate Lagrangian.** Without this, the α closed form is "well-motivated numerology." With it, α is a derivation chain from K, ρ, ξ alone. Highest single payoff.

2. **Derive V_13 (equivalently the missing "cell-inventory" integer) from substrate.** Same problem class as n_M. Needed for the PMNS-as-α-multiples claim to mean anything.

3. **Solve the density perturbation amplitude problem.** The 10¹³ mismatch in substrate-saturation cosmology is the biggest *quantitative* failure in the framework. Either a mechanism is missing (likely something about how strain modes seed inhomogeneities at de-saturation) or the cosmology section needs to retreat.

4. **Close the loop on γ (drag).** Demonstrate that γ is fixed by K, ρ, ξ via the substrate Lagrangian alone, with no input from observed masses. Until then, γ is a fifth primitive masquerading as derived.

5. **Derive M_Pl/v_EW = exp(4π² − 1) from the substrate, not from inventory.** The hierarchy formula is striking (0.093% in the exponent) but currently rests on the same "B3 integer happens to land here" pattern. A substrate-Lagrangian derivation would convert this from suspicious to definitive.

6. **Reconcile the parallel languages.** Pick one: either the substrate Lagrangian is the engine and B3 integers are derived consequences, or B3 integers are the engine and substrate is one realization of them. The current "both at once" lets the framework dodge each derivation challenge by pointing at the other side.

7. **Derive the cube-cell dark-matter mass scale (27.5 GeV) from substrate dynamics**, and/or honestly accept that it is currently a postulate. Direct-detection experiments will start probing this regime in the next 5 years.

---

## 5. Falsifiability timeline — experiments testing the substrate in the next 5 years

| Year | Experiment | What it tests | Substrate prediction | Failure mode |
|---|---|---|---|---|
| 2026 | DESI DR3 (full Σm_ν) | Neutrino mass sum | 60.5 meV; passes ΛCDM bound (<64.2), fails strict FC (<53) | Tight FC bound that holds across analyses kills the chain |
| 2026 | KATRIN final m_ν_e | Lightest-state floor | m_1 ≈ 2.26 meV (lightest); Σm_ν chain | Direct laboratory mass above ~0.5 eV breaks the chain |
| 2026–27 | LIGO O5 / Voyager | GW speed, polarization modes | Speed = c exactly; only 2 tensor polarizations | Any scalar/vector polarization mode disconfirms transverse-traceless-only graviton |
| 2026–28 | LHC Run 4 / HL-LHC | New particles in the 1–10 TeV window | NONE in this band; dark matter is 27.5 GeV cube-cell, EM-dark | Discovery of any new charged or colored state would force major revision |
| 2027 | Hyper-K, DUNE first physics | δ_CP, θ_23 octant | δ_CP = −π/2 maximal, sin²θ_23 = ½ + 2πα | Strong non-maximal δ_CP or wrong θ_23 octant disconfirms PMNS-from-α |
| 2027–28 | LZ, XENONnT, PandaX-4T at full exposure | Direct DM detection | NULL (cube-cell has zero charge radius, only quadrupole) | A WIMP-like signal at standard cross-sections kills the cube-cell candidate |
| 2027–29 | CMB-S4 / LiteBIRD | n_s, r, Σm_ν, polarization | n_s = 1 − 1/(8π²) = 0.9873 (vs Planck 0.965); r small | A precise n_s of 0.965 kills the formula; r > 0.01 stresses the saturated-bouncing picture |
| 2028 | Euclid full release | σ_8, structure growth | σ_8 = 0.783 (captures both tensions) | Convergence of σ_8 to ΛCDM Planck value with no tension would be awkward |
| 2028–30 | EHT next-gen | BH horizon structure | σ = ½ saturated interior; no singularity | Sharp deviation from Kerr metric in horizon-scale imaging |
| 2029–30 | LIGO/LISA stochastic GW background | Pre-CMB phase transition | De-saturation transition signature; not inflation | A clean inflationary stochastic background pattern forces retreat |
| 2030 | DUNE atmospheric ν, JUNO mass ordering | Mass hierarchy | Normal ordering compatible | Inverted ordering would force re-examination of m_1 = 2.26 meV |

The framework is *more* falsifiable than typical BSM extensions because it has clean numerical predictions, not loose parameter ranges. The next 5 years will substantially narrow the live-vs-dead window.

---

## 6. Honest scoreboard

### Sector-by-sector

| Sector | Status | Confidence |
|---|---|---|
| Ontology-forced (CPT, GW speed, ALPHA-g antimatter gravity) | PASSED | High — these were forced before measurement |
| Atomic spectroscopy (H, He hydrogenic, Lamb, 21cm) | PASSED inherited | High — but inheritance from QED, not original |
| Strong-field GR (Mercury, light bending, GPS, Pound-Rebka) | PASSED | High — substrate Lagrangian gives these from §18.39 |
| Lepton ratios (m_μ/m_e, m_τ/m_μ, Koide) | MATCHES | Medium — depends on B3 integer assignment, not yet derived |
| α (closed form 0.004%) | MATCHES | Medium — striking, but n_M = 268 not derived |
| PMNS / CKM | MATCHES | Medium-low — α-multiple formulas are post-hoc integer matches |
| QCD spectrum (current quarks, α_s) | MATCHES at 1–7% | Medium — expected accuracy of an inventory matching |
| Cosmology (Σm_ν, H_0, σ_8) | PASSES bounds, captures tensions | Medium — internally consistent, depends on chain |
| Dark matter identity (27.5 GeV cube-cell) | UNTESTED prediction | Low — concrete and falsifiable, but no direct support |
| Density perturbation amplitude | FAILS by 10¹³ | Hard fail — biggest open problem |
| Hierarchy (M_Pl/v_EW = exp(4π²−1)) | MATCHES at 0.093% in exponent | Medium — same caveat as α |
| Baryogenesis η = exp(−21) | MATCHES at 24% | Low — large residual called a success is generous |
| Higgs mass via sine-Gordon × saturation | MATCHES at 0.23% | Medium — clean formula, transition incomplete |
| Strong CP | RESOLVED via UV symmetry | Medium — needs more derivation |

### Aggregate scoreboard

- **Ontology-forced wins:** 3 (CPT, GW speed, antimatter gravity). These are the framework's spine.
- **Sub-1% matches with parameter-freeness defensible:** ~10 (gravity/EM ratio, Mercury, light bending, α, Cabibbo, several PMNS, Higgs mass via clean Lagrangian, m_W, sin²θ_W, m_b, m_t, Ω_DM/Ω_b, Ω_b).
- **Sub-1% matches with parameter-freeness contestable:** ~15 (most of the §5.1 table — they involve integer assignments that are not yet derived from the substrate alone).
- **Hard quantitative failures:** 1 (density perturbation amplitude).
- **Live falsifiers in next 5 years:** ~10 (see §5).
- **Free parameter count:** Honestly somewhere between 6 (substrate primitives + saturation cap + orientability) and 10 (if Δ₁, Δ₂, ε_0 are still in play). MODEL.md §9 shows 9+1 = 10 explicitly. The 4-parameter claim is aspirational, not yet operational.

### One-line verdict

The framework's *ontological core* is unusually rigid for a fundamental-physics proposal — it makes ontology-forced predictions that all three of CPT, GW-speed, and antimatter-gravity have actually confirmed. Its *quantitative scaffolding* (B3 integers, drag closures, PMNS-as-α-multiples) is in mid-construction: many of the matches are real, but the derivation chains from the substrate Lagrangian to the integers are not yet closed. The single largest failure (density perturbations, off by 10¹³) is acknowledged in the memory notes but not in MODEL.md, and that asymmetry of presentation is itself a thing to fix.

---

## 7. Most important next steps

In priority order:

1. **Close the n_M = 268 derivation from the substrate Lagrangian.** Single highest-leverage open problem.
2. **Pick one engine.** Either substrate is foundational and B3 integers fall out, or B3 integers are foundational and substrate is one model. Stop running both in parallel as alternative excuses.
3. **Address the density perturbation 10¹³ gap explicitly in MODEL.md**, even if it stays open. Honesty about a hard failure increases overall credibility.
4. **Close γ as a derived quantity**, not a fifth primitive. If γ comes from K, ρ, ξ via Lagrangian dynamics, the framework's input count claim is real. If it doesn't, say so.
5. **Run pre-registered predictions** against DUNE δ_CP, LZ/XENONnT DM searches, CMB-S4 n_s. The framework is at the stage where pre-registration is its biggest credibility multiplier.
6. **Converge MODEL.md §9 with §6.1.** Either keep the Yukawa or drop it for sine-Gordon × saturation. Currently both exist.
7. **Audit the 12 B3 integers** one-by-one against substrate: which are derivations, which are sketches, which are postulates. Produce a table. The framework cannot claim "tighter than SM" until that table is honest.

---

## 8. What would falsify the framework

A short, sharp list. The framework is dead if:

1. Any propagating field is ever measured at speed ≠ c with high confidence (kills the medium ontology).
2. Antimatter is measured to gravitate differently from matter (kills the symmetric-channel hypothesis).
3. CPT is broken at any level (kills substrate orientation symmetry).
4. δ_CP turns out to be far from −π/2, or sin²θ_13 is far from 3α, with no integer-relabeling rescue available.
5. Direct DM detection finds a WIMP-like signal at standard cross-sections (kills cube-cell EM-dark candidate).
6. CMB n_s converges to 0.965 ± 0.001 with no overlap with 0.9873 (kills the substrate inflation-replacement formula).
7. A new charged or colored state is found at LHC in the 0.1–10 TeV band (kills "no new particles" claim).
8. The density perturbation problem is shown to be unfixable in any saturation-class cosmology (the 10¹³ gap becomes terminal rather than open).

The framework is robust if ≥6 of these 8 keep going its way over the next 5 years.

---

## 9. Final note on style

The framework is being developed in the open by a single researcher with a clear methodological principle ("just do the math, don't hedge"). That style is visible in the codebase: 250+ scripts, most of them under a few hundred lines, each attacking one observable or one mechanism. It is not the style of a finished theory. It is the style of an *exploration program* with a definite ontological commitment, and it should be evaluated as such.

The strongest critique of the framework is not that any one prediction is wrong (most of them are right or close). The strongest critique is that the parallel languages (substrate Lagrangian vs. B3 inventory) make it hard to tell which predictions are *forced* and which are *fit*. Every numerical success in §5.1 deserves a tag: "ontology-forced," "Lagrangian-derived," "inventory-matched," or "post-hoc fit." Until that tagging exists, outsiders will (rightly) discount the aggregate.

Build the tagging. The framework is good enough to deserve it.

---

*audit_00 — written 2026-05-01. To be revised when audit_01 through audit_09 are produced and the synthesis can fold in their independent findings.*
