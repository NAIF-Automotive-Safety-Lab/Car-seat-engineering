# Audit 08 — MODEL.md & geom_*.py vs Simulation Code

**Date:** 2026-05-01
**Scope:** `/Users/hendrixx./Desktop/untitled folder/MODEL.md` + `scripts/geom_01..10_*.py`
versus `src/stiff_medium/` (126 modules) and `tests/` (62 test files).
**Status:** Read-only audit. No files modified.

---

## 0. Top-line verdict

MODEL.md is in **substantially good shape** post the 2026-05-01 substrate-drag
closure pass: drag γ is correctly listed as a primitive (line 32), the Higgs
Mexican-hat is explicitly dropped (line 41, §6.1), and §5.1 reflects the
26+34-observable inventory. The 10 geom scripts cover almost exactly the right
list of mechanisms.

The drift problems are concentrated in five areas:

1. **§8 module list (lines 567–593) is severely stale** — it lists only 9 of
   126 actual modules and omits every key 2026-05-01 module (drag, mass-torque,
   primitive-anchoring, alpha-closed-form, mobius_k4_numerical, the
   geom_*-companion src modules, etc.).
2. **§9 free-parameter table (lines 649–664) is internally inconsistent** with
   §1.2/§4.3** — §1.2 and §4.3 claim 4 continuous primitives (K, ρ, ξ, γ) but
   §9 still lists 9 continuous parameters including g_Y, e, ε_0, Δ₁, Δ₂.
   §6.0 explicitly closed Δ₁, Δ₂ at <0.5%.
3. **§9 Lagrangian is the OLD Mexican-hat-adjacent form** with `g_Y φ ψ̄ψ`
   Yukawa coupling and `−ε_0` vacuum offset — both contradict the §2.5 +
   §6.1 statements that drag replaces Yukawa and the sine-Gordon × saturation
   has no µ²/no offset.
4. **Several MODEL.md topics are completely missing** that the code now
   produces: high-T_c bound (128.9 K), Majorana neutrino prediction, m_1
   lightest neutrino (2.26 meV), the 9+ banner falsifier results from
   `b3_passed_banner_tests`, and the σ=½ Möbius-Z/2 fixed-point derivation.
5. **The Lieb-Oxford retraction is not in MODEL.md** because Lieb-Oxford was
   never in MODEL.md to begin with — but `scripts/lieb_oxford_closure.py` now
   says "double-counted" (lines 92, 124, 210). Not a doc-fix issue for
   MODEL.md, but RESULTS.md may need a check.

---

## 1. OUTDATED CLAIMS (need correction)

### 1.1 §8 Implementation list (lines 565–593) — CRITICAL DRIFT

**Current text (lines 567–578):**
> The core dynamics is implemented in `src/stiff_medium/`:
>   - `neutrino.py`: 45° cone primitive
>   - `three_d.py`: 3D propagation
>   - `dynamics.py`: time evolution
>   - `back_reaction.py`: medium response forces
>   - `mobius_dynamics.py`: half-flux topology
>   - `atomic.py`: multi-electron N-body
>   - `spinor.py`: Möbius spinor
>   - `em_field.py` / `em_field_3d.py`: EM wave propagation
>   - `detector.py`: bound-state tracking

**Reality:** `src/stiff_medium/` has **126 .py modules**. Missing references include
the most important post-2026-05-01 ones:

- `mass_torque_engine.py` (drag torque T_drag, cone-bouncing freq)
- `drag_mass_generator.py`
- `primitive_anchoring.py` (electron_drag anchor)
- `alpha_geometric.py`, `alpha_bundle.py`, `alpha_derivation.py`,
  `alpha_two_loop.py`, `alpha_audit.py`, `alpha_tetrahedron.py`
- `mobius_k4_numerical.py`, `mobius_quantization.py`,
  `mobius_dirac_vertex.py`, `mobius_dirac_vertex_extended.py`
- `cabibbo_angle.py`, `cabibbo_substrate.py`, `chromomagnetic_substrate.py`
- `cone_*` family (12+ scripts on cone variational origin, lattice
  selection, two-anchor origin, swap generator origin, etc.)
- `dark_stress_*` family (8+ scripts)
- `bound_state_spectrum.py`, `confinement_potential.py`,
  `cosmology_simulator.py`, `bbn_predictions.py`, `baryogenesis.py`,
  `back_reaction.py`/`back_reaction_v2.py`
- 343 tests across `tests/` (the user's context says 343, the file count is
  62 test files — likely 343 = total test functions).

**Recommended replacement (insert at line 567):**

> The core dynamics live in `src/stiff_medium/` (126 modules) with 343
> tests in `tests/`. Major sectors:
>
> *Foundations:* `neutrino.py` (45° cone), `three_d.py`, `dynamics.py`,
> `back_reaction.py`, `primitive_anchoring.py` (4-primitive solve incl.
> `electron_drag` anchor).
>
> *Topology / Möbius:* `mobius_bundle.py`, `mobius_dynamics.py`,
> `mobius_quantization.py`, `mobius_k4_numerical.py`,
> `mobius_dirac_vertex[_extended].py`, `spinor.py`.
>
> *Drag-as-mass:* `mass_torque_engine.py`, `drag_mass_generator.py`.
>
> *α derivation:* `alpha_geometric.py`, `alpha_bundle.py`,
> `alpha_tetrahedron.py`, `alpha_derivation.py`, `alpha_two_loop.py`,
> `alpha_audit.py`, `alpha_s_running_from_K.py`.
>
> *Cone-origin program:* `cone_variational_origin.py`,
> `cone_two_anchor_origin.py`, `cone_swap_generator_origin.py`,
> `cone_diamond_cell_*.py`, `cone_phase_slip_lattice.py`,
> `cone_lattice_microgeometry.py`, `cone_self_dual_exchange.py`,
> `cone_detailed_balance.py`, `cone_bouncing.py`, `cone_anchor_induced_exchange.py`.
>
> *Mixing & EW:* `cabibbo_angle.py`, `cabibbo_substrate.py`,
> `chromomagnetic_substrate.py`.
>
> *Hadronic / nuclear:* `bound_state_spectrum.py`, `confinement_potential.py`,
> `axial_cloud_quenching.py`, `baryon_y_junction.py`.
>
> *Cosmology:* `cosmology_simulator.py`, `cyclic_cosmology.py`,
> `bbn_predictions.py`, `bbn_from_substrate_thermal.py`,
> `baryogenesis.py`, `black_hole.py`.
>
> *Dark sector:* `dark_stress_em_darkness.py`, `dark_stress_transport.py`,
> `dark_stress_cluster_dynamics.py`, `dark_stress_memory_clock.py`,
> `dark_stress_factor_scan.py`, `dark_stress_hybrid.py`,
> `dark_stress_parameter_closure.py`, `dark_stress_scale_closure.py`,
> `cube_cell_dm_simulator.py`.

### 1.2 §9 free-parameter table (lines 651–664) — INCONSISTENT WITH §1.2/§4.3

**Current (line 664):** "9 continuous + 1 binary = 10 parameters total."
**§1.2 (line 25) and §4.3 (line 320):** "4 continuous + 1 binary."

The §9 table still lists `g_Y`, `e`, `ε_0`, `Δ₁`, `Δ₂` as separate
parameters — but:

- **g_Y (Yukawa)** — §2.5 (lines 117–148) and §6.1 (lines 462–476) say
  drag γ replaces Yukawa entirely. `mass_torque_engine.py` lines 56–72
  confirm code uses drag γ, not g_Y. **Remove from table.**
- **e (bundle charge)** — §6.0 line 428 closes α to 0.004% from
  K_4 + Möbius + drag, no free e. **Remove from table.**
- **ε_0 (vacuum offset)** — §6.1 line 471 says
  "No cosmological-constant catastrophe (V(0) = 0)" → no ε_0. **Remove.**
- **Δ₁, Δ₂** — §6.0 line 432–435 closed both at <0.5% via
  `n_G(k_r²-k_p)` and `n_A·n_G/(k_e-k_p)`. **Remove from table.**
- **φ_max (saturation)** — §2.6 says σ_max = ½ structurally; not free.
  Could remove, but §9 introduces it as a Lagrangian symbol; mark
  "= ½ structurally."

**Recommended new table:** K, ρ, ξ, γ (drag), Möbius binary = 4+1.

### 1.3 §9 Lagrangian (lines 599–625) — STILL HAS DROPPED TERMS

Lines 603–606 still print:

```
+ ψ̄(iℏγ^μ(∂_μ + ieA_μ) − g_Y φ)ψ           [Dirac fermion + EM + Yukawa]
```

and line 610:

```
V(φ) = (K/ξ²)(1 − cos(φ/ξ))/√(1 − (φ/φ_max)²) − ε_0
```

But §6.1 lines 466–471 explicitly drop the Mexican hat *and* set V(0)=0:

> Use only sine-Gordon × saturation:
>     V(φ) = (K/ξ²)(1 − cos(φ/ξ)) / √(1 − (φ/φ_max)²)
> - V(0) = 0 (no cosmological-constant catastrophe)

And §2.5 says drag replaces Yukawa entirely. Recommendation:

- **Line 604:** replace `− g_Y φ` with `− γ ξ ρ √K · f_int(config)` (the
  drag mass term as written in `mass_torque_engine.py` line 167).
- **Line 610:** drop `− ε_0` and update §9 free-parameter table.
- **Line 643 ("g_Y φ ψ̄ψ Yukawa | Mass generation, gravity"):** rewrite as
  "γ drag coefficient | Mass generation via cone-bouncing ω_b = ℏ/m c²".

### 1.4 §6.3 Open Problems list — line numbering broken (lines 519–530)

The list jumps `13, 14, 4, 5, 6, 7, 8, 9` — items 4–9 should be 15–20.
Pure cosmetic, but worth fixing along with the substantive cleanup.

### 1.5 §10 Status (lines 686–696)

Line 692: "Working parameters: 10–20 depending on effective-sector counting"
**should be** "4 continuous + 1 binary topology" to match §1.2/§4.3/§6.0.

Line 694 lists "α, lepton excitation energies, CKM/PMNS, current quarks" as
open work. **All four are now closed in §6.0 and §6.1.** Should read:

> Open work: full CKM beyond Cabibbo (V_cb, V_ub, V_td), color-confinement
> dynamics and α_s(μ) running, M_Pl UV closure, full CMB/Hubble transfer
> function, dark-sector neutral-stress dynamics (cluster lensing offset),
> nonlinear strong-field GR.

---

## 2. MISSING TOPICS (should be added)

### 2.1 High-T_c bound — `T_c,max = Λ_QCD / R_baryon = 128.9 K`

Not mentioned anywhere in MODEL.md. Fully derived in
`scripts/condensed_matter_substrate.py` lines 39, 220–222, 493–509. This
is a banner cross-disciplinary result: matches the 31-year ambient-pressure
record (HgBaCa₂Cu₃O₈ 134 K) at 4%. **Add to §5.1c table** as:

| Condensed matter | T_c,max ambient SC | Λ_QCD/R_baryon = 128.9 K | 4% |

### 2.2 Lightest-neutrino prediction (m_1 = 2.26 meV) and Σm_ν

§5.1 row "Σm_ν | from cell-inventory | <DESI bound" exists but is vague.
The current state:

- m_1 = 2.26 meV (line 487 mentions in passing for ρ_Λ)
- Σm_ν = 60.5 meV (DESI DR2 falsifier per `b3_critical_falsifier_sigma_mnu`)
- Tension: passes ΛCDM (<64.2) but fails strict FC (<53) by 14%

**Add a sub-section** under §6 or §5.1 explaining this is a *live* falsifier
not yet decided, and that B3 predicts a specific Σm_ν, not a bound.

### 2.3 Majorana neutrino prediction

User context lists this as a B3 prediction, but neither MODEL.md nor any
`scripts/geom_*.py` mentions it explicitly. No dedicated `*neutrino_mass*`,
`*0vbb*`, or `*majorana*` script in `src/`. **Either add to docs** or
**flag as missing implementation** — currently inconsistent with claimed
prediction inventory.

### 2.4 σ = ½ as Möbius Z/2 fixed point

`geom_08_saturation_cap.py` line 101–114, 380 establishes σ=½ as the
"unique fixed point of the Z/2 sheet-swap involution s ↔ 1−s" and marks
it FORCED. But MODEL.md §2.6 (line 153) just *asserts* "σ_max = ½"
without the topological derivation. **Add 2-line derivation** under §2.6:

> The cap σ = ½ is the unique fixed point of the Z/2 sheet-swap
> involution on the Möbius bundle (s ↔ 1−s), so the saturation value
> is forced by the same topology that gives spin-½. See
> `scripts/geom_08_saturation_cap.py`.

### 2.5 Drag-as-mass mechanism — under-prominent in §1/§9

§2.5 explains drag-bouncing well, but §1.2 line 32 only says "dissipation;
generates rest mass via cone-bouncing" without the explicit replacement
of Higgs Yukawa. §9 still leaves Yukawa in the Lagrangian. Add an
explicit subsection §1.2.1 or §2.5.1 stating:

> **Drag γ replaces the Higgs Yukawa coupling.** Particle rest mass is
> set by `m c² = ℏ ω_bounce` where ω_bounce = γ ξ ρ √K · f_int(config).
> The SM's Higgs vev + Yukawa table compresses to a single substrate
> primitive γ plus combinatorial f_int from configuration topology.
> See `src/stiff_medium/mass_torque_engine.py` (lines 56–72, 160–173)
> and `src/stiff_medium/drag_mass_generator.py`.

### 2.6 m_μ/m_e = exp(n_M / 16π) at 0.009% with 16 = K_pair⁴ derived

§5.1 line 347 lists `m_μ/m_e = n_G(k_r²-k_p) = 207` at 0.11%. But
`geom_06_generations.py` line 84–91 shows the *exp(n_M/(16π)) form* at
0.009% with `16 = K_pair⁴ = 2⁴` (4-form on Möbius doubled space). User
context confirms this is the canonical current form. **Replace line 347**
with:

| Lepton | m_μ/m_e | exp(n_M/(K_pair⁴·π)), n_M=268 | **0.009%** |

(or list both, with the exp form as the primary derivation.)

### 2.7 3 generations forced from D=3

User context: "3 generations forced from D=3". `geom_06_generations.py`
lines 76–77, 244 derive this. MODEL.md §3.4 line 226 just says
"3 lepton generations exactly (vertex closure caps stress quanta at 3)" —
the **vertex closure** explanation is *one* derivation but the cleaner
**D=3 spatial dimensions ↔ 3 generations** mapping is missing. Add to
§3.4 or §2:

> **3 generations from D=3.** Each generation = one orthogonal cone-bounce
> harmonic in the substrate's 3 spatial directions. D > 3 substrate would
> give 4+ generations; D = 2 would give 2. Empirical 3-generation
> exhaustion forces D = 3. See `scripts/geom_06_generations.py` lines 76–77.

### 2.8 Banner zero-parameter wins (CPT, GW=c, ALPHA-g)

User memory `b3_passed_banner_tests` lists three zero-param wins not
explicit in MODEL.md table:

- CPT: m_p̄ = m_p at 16 ppt (BASE)
- GW speed = c at 10⁻¹⁵ (GW170817)
- Antimatter gravity (ALPHA-g)

§3.3 line 215 mentions GW speed = c but not the precision; CPT and ALPHA-g
absent. **Add to §5.1 or new §5.1d** as zero-parameter falsifier passes.

---

## 3. INCONSISTENCIES BETWEEN MODEL.md AND geom_*.py SCRIPTS

| Topic | MODEL.md says | geom_*.py / src says | Action |
|---|---|---|---|
| n_R | n_R = 18 (line 39 chain) | `geom_06.py` line 28: `n_R = 12`; `geom_02.py` derives n_R = 18; `geom_06.py` then *adopts* n_M = 268 directly | Reconcile: state n_R = 18 canonically and explain `geom_06` numerical-ansatz override |
| K_rank | implied K_rank = 5 in §5.1 footer text | `geom_03.py` line 36: `K_RANK = 5`; `geom_06.py` line 27: `K_rank = 4` | Inconsistency *within* the geom scripts, propagates to MODEL.md silently. Fix `geom_06.py` to K_rank = 5, document |
| Λ_QCD anchor | not pinned in MODEL.md | `geom_03,05,07,10`: 200 MeV; `geom_04`: 217.9 MeV; `geom_08`: 217 MeV; user MEMORY: 220 MeV | **Add §1.3.1** stating the canonical Λ_QCD = 200 MeV anchor and flagging the ±10% scatter across scripts as a known cleanup |
| Mass mechanism | §2.5: drag γ; §9: g_Y Yukawa | All src/ uses drag γ exclusively (`mass_torque_engine`, `drag_mass_generator`, `primitive_anchoring.electron_drag`) | Strip Yukawa from §9 Lagrangian (covered in 1.3) |
| ε_face | §5.1 lines 351, 354: ε_face used without definition | `geom_07.py` line 39: `ε_face = Λ_QCD/(n_A·N_BAM/4.5) = 200/90 = 2.222 MeV` | Add 1-line definition in §5.1 footer or §3 |
| n_BAM | not in MODEL.md | `geom_05.py` line 29: derives n_BAM = 6; `geom_07.py` line 36: uses N_BAM = 9 | Pick one: 6 or 9. Currently the geom scripts disagree |

---

## 4. CROSS-CHECK: §5.1 RESULTS vs CODE OUTPUT

I did not run the simulations during this audit (read-only). The following
formula-vs-script consistency check on §5.1:

| MODEL.md row | Formula in MODEL.md | Script that produces it | Match |
|---|---|---|---|
| α(0) | 11/(48π³)·exp(-3π/737) | `scripts/alpha_closed_form.py`, `src/.../alpha_geometric.py` | ✓ formula present |
| Cabibbo | 1/(π√2) | `scripts/cabibbo_substrate_test.py`, `src/.../cabibbo_angle.py` | ✓ |
| Higgs m_H | √(4/15)·v_EW | `clean_lagrangian_proposal.py` | ✓ |
| sin²θ_12 = 42α | — | `scripts/pmns_complete.py` | needs verification (script not opened) |
| m_μ/m_e = 207 | n_G(k_r²-k_p) | `geom_06_generations.py` uses exp(n_M/16π) instead | **mismatched form** (see 2.6) |
| Higgs λ_H = 2/15 | K_pair/n_A | needs `clean_lagrangian_proposal.py` re-check | likely OK |
| α_s(M_Z) = 16α | 16 = (Strand+1)² | `src/.../alpha_s_running_from_K.py` | ✓ src module exists |
| Σm_ν | "from cell-inventory" | no dedicated `neutrino_mass.py` in src; lives in `cosmology_simulator.py`? | **needs locator** |
| η baryogenesis | exp(−21) | `scripts/baryogenesis_test.py`, `src/.../baryogenesis.py` | ✓ |

Recommend: add a column `script` to §5.1 table giving the canonical
producer for each row (the table is a claim, but currently un-traceable
to specific code).

---

## 5. RECOMMENDED DOC UPDATES (priority-ordered)

### P0 (false claims / contradictions)
1. **Rewrite §9 free-param table (lines 651–664)** to 4+1 (K, ρ, ξ, γ + Möbius).
2. **Strip g_Y Yukawa from §9 Lagrangian (line 604)** and `−ε_0` from V(φ) (line 610). Update line 643 of "what emerges" table.
3. **Update §10 line 692 + 694** to remove already-closed items.
4. **Reconcile K_rank inconsistency (geom_03 vs geom_06: 5 vs 4).** Likely a pure typo in `geom_06.py` line 27, but it makes n_M derivation print "140 ?" which is sloppy.

### P1 (missing topics)
5. **Add §1.2.1 drag-as-mass subsection** (text in 2.5 above).
6. **Add §2.6 Z/2 fixed-point derivation of σ=½** (text in 2.4).
7. **Add §3.4 or §2 D=3 → 3-generation derivation** (text in 2.7).
8. **Add §5.1c rows for T_c,max = 128.9 K and CPT/GW/ALPHA-g zero-param wins.**
9. **Add Σm_ν / m_1 prediction subsection under §6** with falsifier status.

### P2 (stale module list)
10. **Replace §8 module list (lines 567–593)** with the 12-sector grouping in §1.1 above.

### P3 (cosmetic)
11. Fix §6.3 list numbering (lines 519–530 jump 13,14,4,5,...).
12. Add `script:` column to §5.1 table for traceability.
13. Pin Λ_QCD anchor in MODEL.md and audit the geom scripts for the 200/217/220 MeV scatter.
14. Resolve N_BAM = 6 (geom_05) vs N_BAM = 9 (geom_07) within the geom series before any of these propagate into MODEL.md.

### P4 (Majorana)
15. Either add Majorana neutrino prediction to MODEL.md and create a script
    in `src/`, or remove it from the user's prediction inventory. Currently
    neither doc nor code references it.

---

## 6. WHAT IS ACCURATE (worth preserving)

- §1.2 4-primitive table including drag γ (line 25–34) — correct.
- §2.5 cone-bouncing / drag mass mechanism — correct and complete.
- §5.1 main table — formulas correct as far as spot-checked; coverage
  matches user memory `b3_pitch_audit_insights`.
- §6.0 closures — accurate against current code.
- §6.1 clean Lagrangian without Mexican hat — correct (and matches
  `clean_lagrangian_proposal.py`).
- §6.2 cosmology closures (hierarchy, Λ, n_s) — accurate.
- The geom_*.py series itself is internally well-organized: each script
  has a clear scope, honest failure-mode flagging (`geom_05`, `geom_08`),
  and follows a consistent "rigorous + honest gap" pattern.

---

## 7. Files referenced

Primary docs:
- `/Users/hendrixx./Desktop/untitled folder/MODEL.md`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_01_substrate_foundations.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_02_mobius_bundle.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_03_k4_tetrahedron.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_04_cube_cell.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_05_packing_NBAM6.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_06_generations.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_07_cell_stacking.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_08_saturation_cap.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_09_orientability_indices.py`
- `/Users/hendrixx./Desktop/untitled folder/scripts/geom_10_mass_torque.py`

Key src modules sampled:
- `/Users/hendrixx./Desktop/untitled folder/src/stiff_medium/mass_torque_engine.py`
- `/Users/hendrixx./Desktop/untitled folder/src/stiff_medium/drag_mass_generator.py`
- `/Users/hendrixx./Desktop/untitled folder/src/stiff_medium/primitive_anchoring.py`

Cross-disciplinary (T_c bound):
- `/Users/hendrixx./Desktop/untitled folder/scripts/condensed_matter_substrate.py`
  (lines 39, 220–222, 493–509)

Lieb-Oxford retraction (not in MODEL.md, just FYI):
- `/Users/hendrixx./Desktop/untitled folder/scripts/lieb_oxford_closure.py`
  (lines 92, 124, 210)
