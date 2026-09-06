# Audit 02 — Geometric Model (geom_01 through geom_10)

**Date:** 2026-05-01
**Scope:** 10 scripts at `/Users/hendrixx./Desktop/untitled folder/scripts/geom_*.py`
**Methodology:** Each script run end-to-end; outputs cross-referenced against (a) declared targets, (b) sister scripts' values for shared integers, (c) experimental anchors cited.
**Verdict in one sentence:** the geometric model contains a real spine of textbook-rigorous topology (Möbius, K_4, packing) plus several numerical hits that look forced *given* labelled ansatze, but the inter-script integer ledger is **not internally consistent** — N_BAM, n_R, K_rank, and N_BAM all carry different values across sister files. This is the single most important finding.

---

## 1. Run status (all 10 execute cleanly)

| Script | Runs? | Headline numerical claim | Match vs claim |
|---|---|---|---|
| geom_01 | yes | D=3 forced; 6 fundamental inputs (K, ρ, ξ, γ, σ≤1/2, orientability) | self-consistent |
| geom_02 | yes | 11/12 fraction (face-dihedral count + MC); α at 4×10⁻⁵ | reproduces |
| geom_03 | yes | n_M = 2·5³ + 18 = 268; deuteron BE 0.11%; α-particle 0.00% (calibrated) | hits target |
| geom_04 | yes | m_DM = 26.45 GeV (target 27.5, −3.8%) | within 4% |
| geom_05 | yes | N_BAM=6 from C_6 rotational subgroup of 2D hex packing | requires 3 ansatze |
| geom_06 | yes | m_μ/m_e at 0.009%; m_τ/m_μ PARTIAL (2.14%); Cabibbo 0.62% | 1 win, 2 partials |
| geom_07 | yes | 7 doubly-magic shells correct; ρ_0 off by **820%** (not 9×; 9.2×) | density way off |
| geom_08 | yes | σ=1/2 forced from Möbius Z/2; Chandrasekhar inherited | clean |
| geom_09 | yes | n_R = 18 from 2·3·3 product | sketched, see §4 |
| geom_10 | yes | mass-torque axiom synthesis; m_τ requires n_T = 410 | scoped, no derivation |

No script crashes. No assertion fails inside a script. Per-script internal arithmetic checks out.

---

## 2. Inter-script consistency: SHARED CONSTANTS DISAGREE

This is the biggest issue uncovered. The same symbol takes different integer values in different files:

| Symbol | geom_02 | geom_03 | geom_05 | geom_06 | geom_07 | geom_09 | geom_10 |
|---|---|---|---|---|---|---|---|
| **N_BAM** | — | 10 | 6 (target) | — | 9 | — | 10 |
| **n_A** | — | 9 | 45 (target) | — | 45 | — | 9 |
| **K_rank** | — | 5 | — | 4 | — | — | — |
| **n_R** | 18 | 18 | — | 12 | — | 18 | — |
| **K_pair** | 2 | 2 | — | 2 | — | 2 | — |

Observations:

- **N_BAM** is used three different ways. In geom_03 and geom_10 it is `10` and combines with `n_A=9` to make the deuteron prefactor `Λ_QCD/(9·10)=2.222 MeV`. In geom_07 it is `9` and combines with `n_A=45` to make the *same* `Λ_QCD/90`. In geom_05 the entire derivation argues `N_BAM=6` from C_6 packing. The 2.222 MeV number is reproduced via three different factorizations of 90 — this is suspicious algebraic flexibility, not a converged definition. **Recommended fix:** pick a single factorization, declare it canonical, and have the other scripts import it.

- **n_R** is `18` in geom_02/03/09 (the place where it actually contributes to n_M = 268) but `12` in geom_06 — geom_06 then *overrides* its own ledger with `n_M = 268` directly (line 31 explicitly: `# B3 doc states n_M = 268. We adopt 268 directly`). So geom_06's K_pair/K_rank/n_R values are not used; only the declared n_M=268 is. This is a documentation defect — geom_06 should either import n_M or get its sub-integers right.

- **K_rank** is `5` in geom_03 (vertex count of the K_5 augmentation: `n_M = 2·5³ + 18 = 268`) but `4` in geom_06. geom_06 has the wrong factorization on the books even though it never uses it.

These inconsistencies do not invalidate the *outputs* (each script reaches its declared target by some local route). They invalidate the claim that the 12-integer inventory is a single coherent ledger. **Highest-leverage fix in the whole component series:** create a single `b3_constants.py` and have all geom_* scripts import from it.

---

## 3. Rigor classification per derivation

Using the user's 5-tier scheme (FORCED / DERIVED / PARTIAL / SKETCH).

### geom_01 — substrate foundations
- D=3 from Λ(ℝ³) polarization count + K_4 + Möbius — **DERIVED** (each ingredient is independently rigorous).
- σ≤1/2 self-consistency cap — **FORCED** (Z/2 fixed-point).
- u-field, Galilean preferred frame, Möbius half-flux — **ASSUMED** (script labels these honestly).

### geom_02 — Möbius bundle
- Möbius parametrization, w_1≠0, K_pair=2 from double cover — **FORCED** (textbook).
- 11/12 from face-dihedral count + MC — **DERIVED** (combinatorial route is exact; the deeper 1-loop integral that *gives* 11 in α is sketched).
- α = 11/(48π³)e^(−3π/737) at 4×10⁻⁵ — numerical match is excellent but **PARTIAL** as a derivation: the 48π³ volume factor and 737 are not derived inside this script.
- n_R=18 from "6 half-int × 3 int" box construction — **SKETCH**; lands on the right number but the cohomology argument is not actually written out.

### geom_03 — K_4 tetrahedron / nucleon cell
- K_rank=5 augmentation, n_M = 2·125+18=268 — **DERIVED** *given* the augmentation prescription, but the augmentation itself is an ansatz (the apex vertex is "the conserved baryon-number current" — asserted, not constructed).
- Deuteron BE 0.11% — **DERIVED** *given* ε_face = Λ_QCD/(n_A·N_BAM) = 200/90; the 90 denominator is the contested one (see §2).
- α-particle BE 0.00% — **CALIBRATED**: the inter-bipyramid term is fitted to close the gap. Predicts 6.74 vs derived "K_rank·χ_apex = 6.75" — the χ_apex=1.35 is ad-hoc.
- Heavy-nuclei volume term BE/A=13.33 vs observed 8.79 — surface+Coulomb deficit "explained" as 34% reduction but not derived; **SKETCH**.

### geom_04 — cube cell / dark matter
- O_h symmetry, octupole-leading EM, absolute stability — **FORCED** (group theory).
- m_DM = √(M_p·Λ_QCD) · 26 · 2.25 = 26.45 GeV vs target 27.5 — **PARTIAL**: the closed form `(V+E+F)·3/(V/F) = 26·2.25 = 58.5` is a clever assembly of cube invariants but the prefactor is not derived from substrate dynamics; 3.8% miss.
- LZ null prediction — **DERIVED** (cross-section estimate well below 2024 limits).

### geom_05 — N_BAM=6 packing
- 2D/3D kissing + densest packing facts — **PROVEN** (Hales 2017 etc., correctly cited).
- N_BAM=6 = |C_6| — **FORCED *given*** three ansatze, all flagged in the script:
  - (a) σ≤1/2 selects the sub-dense regime,
  - (b) modes live on a 2D substrate slice through 3D,
  - (c) the *rotational* subgroup of D_6 (not the full reflection group) is the right one.
- **Audit of ansatz (b)** as the user requested: this is genuinely load-bearing. Without (b) the 3D kissing number 12 would give N_BAM=12 (or possibly 6 again via the K_4 vertex-coordination route, but only as a coincidence). Cleaner replacements to consider:
  - The K_4 vertex coordination route (already in the script as a hint) is one independent 3D path that avoids the slice.
  - A spectral argument on the substrate Laplacian's lowest 6 modes per cell would be the rigorous version. Not done.
  - As written, (b) is a *direct postulate* dressed as geometry. The script's honesty marker (`[ANSATZ]`) is warranted.
- n_A=45 origin — **OPEN** (script flags it).

### geom_06 — generations
- 3 generations from D_space=3 — **DERIVED** (mode-closure argument).
- m_μ/m_e = exp(n_M/16π) at 0.009%, with 16 = K_pair⁴ — **DERIVED** but only after geom_06 *overrides its own n_R/K_rank values* with the externally stipulated n_M=268. Honesty marker missing here: the script silently substitutes, then declares "DERIVED 0.009%". The 16=2⁴ identification is genuinely natural.
- m_τ/m_μ = exp(n_M/30π) — script flags PARTIAL; 30 = K_pair·(K_rank²−1) = 2·15 is *post-hoc* factorization. **PARTIAL** is correct.
- Cabibbo sin θ_C = 1/√20, with 20 = K_rank² + 2K_pair = 25−5 = 20 — script gives three different identifications for "20" (`K_rank²+K_pair`, `K_rank²+2K_pair`, `5·K_rank`). When three integer factorizations all hit 20, none is forced. **PARTIAL**.

### geom_07 — cell stacking / nuclear binding
- 7 magic-number shells correct (closed K_4 shells) — **DERIVED** zero-parameter and impressive.
- Doubly-magic identification — **DERIVED**.
- Drip lines on the heavy side — **DERIVED**; light-side off (admitted: missing Coulomb).
- **Saturation density 9.2× off (820%)**: predicted ρ_0 = 1/(6·V_tet) = 1.47 fm⁻³ vs observed 0.16 fm⁻³. The script *says* "ratio ~4%" in the summary but the body shows 820% — the summary line is **wrong** (the 4% claim is a leftover from an earlier σ≤1/2 derivation in geom_08). What would close the 9× gap?
  - The factor 9 ≈ 6·(3/2) is suggestive: Voronoi-per-nucleon would need ≈ 6.25 fm³ instead of 0.68 fm³. That requires either (i) substrate cells being ~9× larger than ξ_QCD = ℏc/Λ_QCD, e.g. each nucleon occupies 9 K_4 cells (→ baryon = composite of 9 substrate cells, which would conflict with geom_03's "1 K_4 = 1 nucleon"), or (ii) adding the σ≤1/2 saturation cap (which would give factor 2, not 9), or (iii) reinterpreting Λ_QCD in this formula as `Λ_QCD/cube-root(9) ≈ 96 MeV`, suspiciously close to Λ_QCD(MS-bar) at low scales.
  - Honest reading: `ξ_QCD = ℏc/Λ_QCD` is too small by a factor of 9^(1/3) ≈ 2.08. Either the Λ_QCD anchor is wrong by 2× for this purpose, or one nucleon contains many cells. **The 9× miss is a real falsifier of the strict "1 K_4 = 1 nucleon" picture.**
- Step 9's `820.29%` printout vs the SUMMARY block's `~4%` claim — this is a code/honesty inconsistency that should be fixed.

### geom_08 — saturation cap σ=1/2
- σ≤1/2 from Möbius Z/2 fixed-point uniqueness — **FORCED**.
- Pauli ↔ Möbius doubling, g_spin=2 — **FORCED**.
- BH horizon = σ=1/2 boundary — **INTERPRETED** (script admits: definition, not derivation).
- Crack-tip K_I/√r + cap → unique r_p — **FORCED** (kinematic).
- Hierarchy −1 in exp(4π²−1) — **PARTIAL** (boundary count forced; coefficient =1 needs full Möbius integral).
- Chandrasekhar 1.4 M_⊙ — **INHERITED** (textbook derivation; σ=1/2 enters as g_spin=2). Script honest.

### geom_09 — orientability indices
- Orientability table, w_1/w_2 characterization — **TEXTBOOK**.
- 11/12 simplex subtraction — **DERIVED** (MC at 3×10⁻⁵).
- Parity-violation order of magnitude (Mason-Tranter) — **DERIVED**.
- **n_R = 18 = 2·3·3 audit** as user requested:
  - 2 = rk H¹(T²; Z/2), correct.
  - 3 = "ambient ℝ³" — really meaning "spatial dimension count".
  - 3 = "reflection modes per axis" — but per axis there is **one** reflection (the orientation-reversing involution), not three. The script's own text: "reflection modes per axis = 3" is unjustified.
  - The correct cohomological count for Möbius bundle over T² with Z/2 coefficients via Künneth is H¹(T²)⊗Z/2 ⊕ Möbius twist. This gives 2+something, not the clean product 2·3·3.
  - **Verdict: SKETCH at best.** The product 2·3·3 = 18 is numerology that lands on the right answer; the script honestly admits "no spectral sequence has been written down here." This is **not** a rigorous join. The user's specific question — "is the 2×3×3 join actually rigorous?" — answer: **no**.

### geom_10 — mass-torque synthesis
- Synthesis statement: m = Λ_QCD · T(config) — **AXIOMATIZED** (this is the foundational principle from MEMORY.md).
- Reproduces deuteron, m_μ/m_e, α, hierarchy, vacuum, CPT, GW from the previous 9 scripts — **DERIVED** (within their respective rigor levels).
- **m_τ derivation requires n_T ≈ 409.86, rounded to 410** — audit of "is there a substrate-natural N_T?":
  - 410 = 2·5·41. Factor 41 is prime, has no obvious K_rank/K_pair/n_R identification.
  - 410 ≈ n_M · (3/2) · 1.02 = 268·1.53 — not clean.
  - 410 = n_M + n_M/(2π)·... — does not reduce.
  - 410 vs 30π·m_τ/m_μ identification (the geom_06 route): if m_τ/m_μ = exp(n_M/30π), then m_τ/m_e = (m_τ/m_μ)·(m_μ/m_e) = exp(n_M/30π + n_M/16π) = exp(n_M·(16+30)/(30·16π)) = exp(n_M·46/480π) = exp(268·46/(480π)) → equivalent 16π form: n_T = 268·46/30 = 410.93. So **n_T ≈ 410 is the same number as n_M·46/30**, i.e. it is forced by the geom_06 N=30 partial derivation, not independently. The "410 must be derived" problem is the **same** problem as "30 must be derived" (where 30 = 2·15 = K_pair·(K_rank²−1) is ad-hoc post-hoc). 
  - **Verdict: there is no substrate-natural N_T separate from N(τ/μ)=30.** Fixing one fixes the other; the open derivation is the m_τ/m_μ ladder integer.

---

## 4. Honesty markers — are they warranted?

| Script | Marker | Warranted? | Comment |
|---|---|---|---|
| geom_02 | "11/12 sketched at deeper level" | **yes** | combinatorial = rigorous, 1-loop = not done |
| geom_02 | "n_R=18 shell-cut not first principles" | **yes** | accurate self-criticism |
| geom_05 | three [ANSATZ] tags + "[OPEN] n_A=45" | **yes** | exemplary honesty |
| geom_06 | "PARTIAL" on m_τ/m_μ N=30 | **yes** | three different "30" factorizations admitted |
| geom_06 | silently overrides n_R/K_rank to hit n_M=268 | **NO** | this is a missing honesty marker |
| geom_07 | summary says "~4%" while body shows 820% | **NO** | **summary is incorrect** |
| geom_08 | INTERPRETED vs FORCED vs PARTIAL labelling | **yes** | best honesty discipline of the suite |
| geom_09 | "n_R=18 as 2·3·3 is a sketch" | **yes** | the script is more honest than the claim |
| geom_10 | "scoped — n_T must be derived" on m_τ | **yes** | doesn't pretend |

---

## 5. Highest-impact open gaps (ranked)

1. **Inter-script ledger drift (N_BAM ∈ {6, 9, 10}, n_A ∈ {9, 45}, n_R ∈ {12, 18}, K_rank ∈ {4, 5}).** This is the single most damaging issue because it makes the "12-integer inventory" claim look unconverged. A unified `constants.py` with one canonical assignment would either (a) work cleanly (in which case the framework gains coherence overnight) or (b) reveal a real choice that is currently being papered over. **Top priority.**

2. **Nuclear matter ρ_0 off by 9.2×.** Falsifies "1 K_4 cell = 1 nucleon at edge ξ_QCD". Either Λ_QCD anchor is wrong by 2× for length scales, or nucleons span ~9 cells. Both reinterpretations have downstream consequences for geom_03 (deuteron) and geom_07 (magic shells, which currently work zero-parameter — so option (b) probably breaks them). This is a **sharp falsifier** and worth treating as such.

3. **n_R = 18 from 2·3·3 is numerology, not derivation.** The "3 reflections per axis" is wrong; the cohomology calc is not done. n_R=18 enters n_M=268 and therefore enters the m_μ/m_e win at 0.009%. If n_R is not actually 18 but the lepton ratio still requires 268, then the "derivation" of m_μ/m_e is a fit. **High leverage** — a real spectral sequence calc would either confirm 18 or expose the lepton-ratio match as coincidence.

4. **m_τ/m_μ ladder integer (30 in geom_06 = 410 in geom_10, same problem).** Three post-hoc factorizations of 30 admitted; no substrate-natural derivation. The user's question "is there a substrate-natural N_T?" — answer is no, and the same is true of N=30. This is **the** open derivation in the lepton sector.

5. **n_A = 45 origin.** Best fit C(10,2). Geom_05 admits OPEN; geom_07 just declares it. Lower urgency than (1)-(4) but still on the list.

6. **m_DM = 26.45 GeV vs target 27.5 (3.8% miss).** Lower priority — within rigidity-grid tolerance and the prefactor `(V+E+F)·3/(V/F)` is geometric. Worth tightening only after (1) is fixed since N_BAM disagreement may shift the target.

---

## 6. What is genuinely solid

To balance the criticism above:

- **geom_01's six-input count** (K, ρ, ξ, γ, σ≤1/2, orientability) is a clean, defensible foundation.
- **geom_02's Möbius topology** is textbook-rigorous; w_1, double cover, half-integer spin from SU(2)→SO(3) are not handwaving.
- **geom_07's 7-out-of-7 doubly-magic shells from closed K_4 shells** is a zero-parameter result that works. Even with the ρ_0 problem, the *sequence* prediction is genuine.
- **geom_08's σ=1/2 forcing from Z/2 fixed-point** is rigorous, and the disciplined FORCED/INTERPRETED/PARTIAL labelling is the audit gold standard for the suite.
- **The mass-torque axiom (geom_10)** is a coherent unifying principle — it correctly identifies that all the previous "successes" share the form m = Λ·T(config). This is a real synthesis, not just rebranding.
- **Cross-validation density**: deuteron BE 0.11%, m_μ/m_e 0.009%, magic shells 7/7, GW speed 10⁻¹⁵, CPT 16 ppt — these are real, with caveats noted above on which integers are truly forced.

---

## 7. Bottom line

The geometric model has **rigorous topological underpinnings** (Möbius, K_4, cube, packing) and **several genuine numerical wins**, but the **integer ledger is currently inconsistent across scripts** and several "derivations" (n_R=18 product, n_T=410, the alternative factorizations of 20/30/90) are post-hoc selection from many possibilities. Fixing the constants drift (gap #1) is the highest-leverage single edit; resolving the 9× nuclear density miss (gap #2) is the highest-stakes physical question.

Honest tier breakdown of the 10 components:
- **3 truly forced** (01, 02 topology core, 08).
- **3 derived-given-ansatz** (03, 05, 09).
- **2 mixed wins + sketches** (04, 06, 07 — magic shells solid, ρ_0 broken).
- **1 axiomatized synthesis** (10).

The framework is far more rigorous than typical numerology and far less rigorous than the synthesis language sometimes suggests. The honest middle ground is exactly what the script labels in geom_05 and geom_08 already say.
