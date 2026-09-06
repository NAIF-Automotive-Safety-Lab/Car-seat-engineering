# Falsifier Dependency Graph — Substrate Framework

**Author:** Analysis pass over RESULTS.md + 11 papers (01–11)
**Date:** 2026-05-01
**Question:** Which substrate predictions are correlated, and which are independent? If one falsifier fires, what propagates?

---

## 1. Catalog of explicit falsifiers across the corpus

Aggregated from RESULTS.md and papers 01–11. Each falsifier is tagged with the substrate axiom or derivation chain it directly tests.

| # | Falsifier (experiment + threshold) | Predicted value | Substrate axiom(s) tested | Year |
|---|---|---|---|---|
| F1 | DESI DR3 Σm_ν < 30 meV | 60.5 meV | de-saturation cosmology + 15/16 doubled-exterior + Möbius-Majorana ν chain | 2026 |
| F2 | DESI DR3+ w(z) ≠ −1 at >5σ | w = −1, w_a = 0 exactly | ρ_Λ = c₄·Λ_QCD⁴/M_Pl² as static geometric ratio (cap σ ≤ 1/2) | 2027–2030 |
| F3 | LiteBIRD primordial r > 10⁻³ | r = 0 (no inflaton) | substrate de-saturation cosmology replaces inflation | 2030 |
| F4 | LEGEND-1000 / nEXO 0νββ at m_ββ > 25 meV | m_ββ ∈ [0.06, 5.3] meV | Möbius Z/2 + NH spectrum + m_lightest = 2.26 meV | 2030 |
| F5 | LEGEND-1000 / nEXO Dirac confirmation (no 0νββ ever) | Majorana forced | Möbius Z/2 sheet-swap on neutral states | indefinite |
| F6 | Discovery of 4th-generation lepton at any collider | exactly 3 generations | spatial D = 3 axis-stack rule | open |
| F7 | TeV-scale superpartners / KK resonances at LHC/FCC | none | exp(4π² − 1) hierarchy without SUSY | 2030–2050 |
| F8 | μ-EDM detection (or n-EDM, e-EDM at strong CP-violating level) | θ_QCD = 0 from Möbius Z/2 | Möbius Z/2 forces vanishing CP-odd vacuum angle | ongoing |
| F9 | GW170817-class \|Δv_GW − c\|/c > 10⁻¹⁵ | v_GW = c structurally | one substrate Lagrangian, two transverse modes | future BNS+EM event |
| F10 | T2K/NOvA sin² θ_23 outside [0.527, 0.567] | 0.546 = ½ + 2πα | σ_max = 1/2 cap + Möbius angular measure | 2025–2027 |
| F11 | JUNO sin² θ_12 outside [0.296, 0.317] | 0.30649 = 42α | K_4 face-edge incidences (24) + n_R orbits (18) | 2025–2030 |
| F12 | Daya Bay sin² θ_13 outside [0.020, 0.025] | 0.02189 = 3α | (K_pair + 1) reflection-line orientations | ongoing |
| F13 | HL-LHC m_H drift outside [125.10, 125.45] GeV | 125.27 GeV | ε_EW · K_pair⁴ · Λ_QCD vacuum-strain mode | 2030 |
| F14 | LHC BR(H → invisible) > 0.005 at 5σ | 0 | no hidden-sector states couple to substrate vacuum-strain | 2030 |
| F15 | BMW HVP shifts → Δa_μ^B3 = 45×10⁻¹¹ tension > 3σ | +45.12×10⁻¹¹ | cone-bouncing drag-loop with ξ_μ = 9/125 | 2025–2030 |
| F16 | EHT-class horizon imaging deviates from r_s = 2GM/c² by >few % | σ = 1/2 cap = horizon | Möbius Z/2 fixed point identified with σ_max | 2025+ |
| F17 | LHC discovers extended Higgs sector (H′, A, H±, doubly-charged) | none | minimal substrate vacuum-strain mode | 2030+ |
| F18 | DUNE δ_CP outside 135° ± 20° | δ_CP ≈ 3π/4 (weak commitment) | Möbius holonomy vs generation-3 cycle (ambiguous) | 2030+ |
| F19 | KATRIN endpoint m_β > 0.45 eV (would falsify any sub-eV scheme) | m_β well below 0.45 eV | NH spectrum from m_lightest = 2.26 meV | 2025 final |
| F20 | JUNO confirms inverted hierarchy (IH) at >5σ | NH preferred (m_1 = 2.26 meV) | substrate NH ansatz | 2027 |

These 20 cluster into **5 high-leverage axiom buckets**:

- **A1 — Möbius Z/2 sheet-swap:** σ = 1/2 cap + Majorana ν + half-integer spin + θ_QCD = 0 + g_spin = 2 + Pauli + CPT m_p̄ = m_p
- **A2 — De-saturation cosmology (no inflation):** r = 0, w = −1 static, ρ_Λ derived
- **A3 — D_space = 3:** exactly 3 generations, 3 PMNS angles structure, K_pair⁴ = 16 (3+1)
- **A4 — One substrate Lagrangian, two transverse modes:** v_GW = c, no graviton mass, photon mass = 0
- **A5 — K_4 / K_5 simplex topology (n_M = 268, K_rank = 5, n_R = 18):** lepton ratios, mixing angles, Higgs mass, g − 2

---

## 2. Five-to-eight falsifier dependency chains

Each chain shows: **trigger** → **what propagates** → **structural reason**.

### Chain A — Möbius Z/2 chain (the largest connected cluster)

**Trigger candidates that all probe the same Z/2 axiom:**
- F4 (m_ββ > 25 meV detected)
- F5 (Dirac neutrino confirmed, eg by long-term 0νββ silence with m_lightest known to be > 5 meV)
- F8 (μ-EDM or n-EDM detection requiring θ_QCD ≠ 0)
- F16 (BH horizon at radius ≠ 2GM/c²)

**Propagates to:**
1. **σ_max = 1/2 cap fails** — ALL of paper 02's seven derivations fall together:
   - Schwarzschild radius coefficient (currently 10⁻¹² verified) becomes wrong
   - Chandrasekhar mass M_Ch = 1.4 M_sun loses derivation
   - g_spin = 2 (Pauli, paper 08 §6) loses substrate origin
   - Crack-tip Dugdale-Barenblatt cohesive zone agreement becomes coincidence
   - Half-integer fermion spin loses geometric origin (becomes a Dirac-equation postulate)
2. **Hierarchy −1 zero-mode subtraction fails** (paper 03 §3.2) → exp(4π² − 1) = 0.093% match degrades to exp(4π²) = 39.48 (1% off in log) — a noticeable but not catastrophic shift
3. **21 cm hyperfine identification** (paper 08 §8) loses Möbius-half-flux interpretation
4. **CPT m_p̄ = m_p at 16 ppt** (b3_passed_banner_tests) loses Z/2 sheet-swap origin
5. **Aharonov-Bohm phase = π** (RESULTS.md §quantum phenomena) becomes coincidence

**Falsifies:** A1 (the entire Möbius Z/2 axiom). This is the **largest single-axiom falsifier cluster** in the framework.

### Chain B — De-saturation cosmology chain

**Trigger candidates:**
- F3 (LiteBIRD detects clean primordial r > 10⁻³)
- F2 (DESI confirms w(z) ≠ −1 at >5σ)
- F1 (DESI DR3 Σm_ν < 30 meV)

**Propagates to:**
1. **Σm_ν = 60.5 meV chain breaks** (paper 04). The chain Σm_ν → ρ_Λ → H_0 becomes inconsistent.
2. **ρ_Λ = c₄·Λ_QCD⁴/M_Pl²** (paper 11) loses its calibration, but the underlying geometric ratio could survive if c₄ is reinterpreted; w = −1 prediction still fails for F2.
3. **H_0 = 71.92 km/s/Mpc** (paper 04 §2) loses its derivation chain — the SH0ES-side preference becomes circumstantial.
4. **σ_8 = 0.783** (paper 04) loses its derivation.
5. **ρ_Λ at 0.04% match** survives ρ_Λ falsification only if the geometric form (not the de-saturation interpretation) is preserved — but if w ≠ −1, the static-cap interpretation is dead.

**Falsifies:** A2. Note: F1 alone (Σm_ν < 30 meV) does NOT necessarily kill de-saturation cosmology — it could still kill only the m_lightest = 2.26 meV anchor while leaving the broader picture. F2 (dynamical w) is the **harder kill** because static geometric ratio is structurally rigid.

### Chain C — D_space = 3 chain (the spatial-dimension count)

**Trigger:** F6 (4th-generation fermion discovered)

**Propagates to:**
1. **K_pair⁴ = 16 doubling factor** (papers 01, 03, 10) loses the "spatial D + bundle direction = 4" derivation. Specifically:
   - Paper 01 m_μ/m_e = exp(268/16π): the "16" in the denominator is K_pair^(D+1) = 2^4. With D = 4 spatial, exponent becomes 2^5 = 32, ratio becomes exp(268/(32π)) = exp(2.666) = 14.39 — wildly wrong.
   - Paper 10 m_H = ε_EW · K_pair⁴ · Λ_QCD: the K_pair⁴ factor changes to K_pair^(D+1) = 32, doubling the predicted Higgs mass to ~250 GeV.
   - Paper 03 hierarchy 4π² Möbius-doubled S³ surface area: 4π² is the surface measure of unit S³ in 4D Euclidean (i.e., D_space + 1 = 4). With D = 4, surface measure shifts to S^4 area = 8π²/3 ≈ 26.32 vs 4π² ≈ 39.48 — log-space match degrades by ~30%.
2. **PMNS angles 42α / 3α / (½ + 2πα)** (paper 06) lose their generation-counting basis.
3. **Cabibbo angle 1/√20** survives only if K_rank = 5 still gives the correct K_5 closure — but the K_5 closure assumption depends on D = 3 (it's the smallest non-trivial closure of a Möbius bundle on K_4 in 3D space).

**Falsifies:** A3 (D = 3). The propagation is **catastrophic** because nearly every closed-form prediction in the framework uses K_pair⁴ = 16. Discovery of a 4th generation would essentially terminate the substrate framework.

### Chain D — Single-Lagrangian / two-transverse-modes chain

**Trigger:** F9 (any clean v_GW ≠ c detection at >10⁻¹⁵)

**Propagates to:**
1. **GW150914 chirp mass formula** (paper 05) loses its "v_GW = c structural identity" justification. The 0.2% chirp-mass agreement becomes a coincidence rather than a substrate prediction.
2. **Photon mass = 0** (RESULTS.md sharp predictions) loses structural origin.
3. **Stefan-Boltzmann σ_SB** (paper 07) survives because σ_SB depends only on (k_B, ℏ, c) — so if c itself is well-defined for transverse modes at thermal scales, σ_SB still works. But the **interpretation** "substrate IS the EM field" becomes shakier.
4. **Atomic Schrödinger equation derivation** (paper 08 Appendix A) survives at non-relativistic atomic scales but the unifying ontology weakens.

**Falsifies:** A4. Less catastrophic than Chain C because A4 is more about ontological unification than numerical predictions — the numerical agreements largely survive; the **explanatory power** is what dies.

### Chain E — K_4/K_5 simplex topology chain (the master integers)

**Trigger candidates:**
- F11 (JUNO sin² θ_12 outside [0.296, 0.317])
- F10 (T2K sin² θ_23 outside [0.527, 0.567])
- F12 (Daya Bay sin² θ_13 outside [0.020, 0.025])
- HL-LHC m_μ/m_e measured at <0.001% precision and substrate's 0.009% becomes a 9σ deviation

**Propagates to:**
1. **K_rank = 5 falls** if the Cabibbo 1/√20 fails alongside multiple PMNS angles. Then:
   - n_M = 2·5³ + 18 = 268 loses structural derivation; m_μ/m_e at 0.009% becomes coincidence.
   - α_s = (K_rank − 1)·α loses derivation.
   - g − 2 ξ_μ = (n_R/K_pair)·K_rank⁻³ = 9/125 loses derivation; Δa_μ^B3 = 45×10⁻¹¹ becomes a fitted number.
2. **n_R = 18 falls** if PMNS θ_12 = 42α (which uses 24+18) fails while θ_13 = 3α and θ_23 = ½ + 2πα survive. Then m_μ/m_e shifts (n_M = 268 ± n_R perturbation = 250 ± 9 → ratio shifts by exp(±9/16π) = exp(±0.179) = ±20%).

**Falsifies:** A5. This chain is the framework's most "soft" — single-angle failures don't necessarily cascade because the K_rank = 5 prediction has multiple cross-checks (Cabibbo, multiple PMNS, multiple mass formulas). Only **multi-angle simultaneous failure** would kill the topology.

### Chain F — Λ_QCD / M_Pl scale chain

**Trigger candidates:**
- F2 (w ≠ −1 at >5σ) breaks ρ_Λ = c₄·Λ_QCD⁴/M_Pl²
- F13 (m_H drift >3σ from 125.27 GeV) breaks ε_EW·K_pair⁴·Λ_QCD

**Propagates to:**
1. If Λ_QCD identification fails: the deuteron binding 2.222 MeV (paper 01 §2.3), high-T_c bound, baryon spectrum (b3_baryon_face_spin_v4 at 0.36%), σ_8 = 0.783, H_0 = 71.92 — **all fail simultaneously**. This is the **highest-leverage anchor** because Λ_QCD enters ≥ 10 derivations.
2. If M_Pl identification fails: hierarchy exp(4π² − 1) at 0.093%, ρ_Λ = c₄·Λ_QCD⁴/M_Pl² at 0.04%, m_H = ε_EW·K_pair⁴·Λ_QCD at 0.02%, H_0 chain — **all fail simultaneously**.

Paper 11 makes this explicit (§5.4 cross-check table): five derivations share Λ_QCD, M_Pl, v_EW; failure in one propagates to all five.

**Falsifies:** scale-anchor consistency (which is an emergent axiom from Λ_QCD and M_Pl being shared substrate constants). This is the single most powerful internal cross-check the framework offers.

### Chain G — Cosmology consistency triangle (Σm_ν, ρ_Λ, H_0)

**Trigger:** F1 confirmed AT THE SAME TIME as F2 holds (i.e., Σm_ν < 30 measured but w = −1 confirmed)

**Propagates to:**
1. Σm_ν derivation chain breaks (m_lightest 2.26 meV anchor invalid), but ρ_Λ chain still works → **partial substrate cosmology survives**.
2. H_0 = 71.92 km/s/Mpc derivation breaks (uses Σm_ν → ρ_Λ → H_0) → SH0ES preference becomes coincidence.
3. σ_8 = 0.783 derivation breaks.

**Verdict:** Σm_ν falsification alone is a **partial kill** of cosmology sector but does not necessarily propagate to particle-physics sector. It does NOT take down the Möbius Z/2 axiom because m_ββ ∈ [0.06, 5.3] meV could still be consistent with smaller m_lightest.

### Chain H — Higgs-sector closure (paper 10 cross-checks)

**Trigger:** F14 (BR(H → invisible) > 0.005) OR F17 (extended Higgs sector discovery)

**Propagates to:**
1. m_H = 125.27 GeV survives numerically (the calibrated factor ε_EW · K_pair⁴ doesn't depend on hidden-sector closure).
2. λ = 0.1293 survives numerically.
3. **The "minimal substrate vacuum-strain mode" interpretation dies** — substrate must be extended to accommodate hidden states. This propagates to:
   - Naturalness: the exp(4π² − 1) suppression argument requires no extra Higgs partners; discovery of any forces revisiting the topological closure argument.
   - SUSY exclusion vs substrate: substrate's null-prediction-passes-LHC argument weakens.

**Verdict:** numerical predictions survive but ontological claims weaken. Soft falsifier.

---

## 3. Load-bearing predictions (failure propagates broadly)

Ranked by **propagation breadth** (number of other predictions that fail with the trigger):

| Rank | Load-bearing prediction | If it fails, what else dies |
|---|---|---|
| 1 | **Möbius Z/2 axiom** (any of F4, F5, F8, F16) | σ = 1/2 cap, BH horizon r_s, M_Ch, half-integer spin, Pauli g = 2, CPT m_p̄ = m_p, 21cm interpretation, 7 distinct paper-02 derivations, hierarchy −1 zero-mode subtraction, Aharonov-Bohm π phase |
| 2 | **D_space = 3** (F6: 4th-generation discovery) | K_pair⁴ = 16 → ALL closed-form lepton ratios, Higgs mass, hierarchy 4π², PMNS structure |
| 3 | **Λ_QCD shared anchor** (multiple Λ_QCD-using predictions degrade simultaneously) | deuteron 2.222 MeV, baryon spectrum 0.36%, high-T_c bound, ρ_Λ at 0.04%, m_H at 0.02%, σ_8, H_0 chain, mixing angles |
| 4 | **M_Pl shared anchor** (multiple M_Pl-using predictions degrade simultaneously) | hierarchy 0.093%, ρ_Λ 0.04%, m_H 0.02%, H_0 chain |
| 5 | **K_rank = 5** (multi-angle failure across F10–F12) | n_M = 268 → m_μ/m_e, Cabibbo, all PMNS, α_s, g − 2 ξ_μ |
| 6 | **De-saturation cosmology** (F2 or F3) | w = −1, r = 0, Σm_ν chain, H_0 chain, σ_8 chain |
| 7 | **One Lagrangian / two transverse modes** (F9) | v_GW = c, photon mass = 0, GW chirp formula justification, ontological unification |

**Top observation:** the Möbius Z/2 axiom is the framework's single most load-bearing structural element. Falsifying it via 0νββ silence or μ-EDM detection would terminate the largest cluster of derivations — paper 02 lists seven, but in fact 10+ predictions across the corpus depend on it.

---

## 4. Independent predictions (failure does NOT propagate)

These are predictions whose falsification would damage at most one or two other claims:

| Independent prediction | Why it doesn't propagate |
|---|---|
| **Stefan-Boltzmann σ_SB** (paper 07) | Depends only on (k_B, ℏ, c) — no substrate-specific integers; failure means substrate ≠ EM-field-as-transverse-mode but doesn't kill K_4 topology |
| **Δa_μ^B3 = 45×10⁻¹¹** (paper 09) | Prediction is leading-order with HVP scenario dependence; failure would kill the cone-bouncing drag-loop diagram but not other K_4 predictions |
| **DUNE δ_CP ≈ 3π/4** (paper 06 §6.2) | Author flags this as a "weaker" prediction; failure leaves PMNS angles intact |
| **GW150914 chirp at 0.2%** (paper 05) | Inherits standard PN inspiral formula; substrate's contribution is structural v_GW = c, so chirp-mass match per se is not load-bearing |
| **He IP at 24.59 eV** (RESULTS.md Tier D) | Calibration sensitivity to Hylleraas multi-parameter fit; no substrate-specific integer at risk |
| **Lyman/Balmer/Paschen** (paper 08) | Standard Schrödinger spectrum; substrate derivation is downstream of long-wavelength limit; failure would mean broader physics is wrong |
| **Higgs invisible BR < 0.005** (F14) | Soft failure: ontological claim dies but m_H = 125.27 GeV survives |
| **Higgs extended sector absence** (F17) | Same as F14: numerical m_H survives, ontology weakens |

These are the framework's "outer crust": failures here do NOT cascade. They are useful **single-target** tests but not framework-terminators.

---

## 5. Smallest set of co-falsifiers that would terminate the framework

To definitively close out the substrate program, the minimum set of coherent simultaneous failures is:

### Minimum lethal set (3 falsifiers, must all fire)

1. **F1: DESI DR3 measures Σm_ν < 30 meV at high confidence**
   → kills de-saturation cosmology m_lightest anchor + Σm_ν chain
2. **F2: DESI DR3+ confirms w(z) ≠ −1 at >5σ**
   → kills ρ_Λ static geometric ratio, no rescue via reinterpretation
3. **F4 OR F5: 0νββ rules out Majorana neutrinos** (either positive m_ββ > 25 meV which is incompatible with substrate spectrum, OR Dirac confirmation)
   → kills Möbius Z/2 axiom → kills σ = 1/2 cap, BH horizon, Pauli, spin-½, paper 02's full cascade

If all three fire, both substrate cosmology (A2) AND the Möbius Z/2 axiom (A1) are dead simultaneously. Without A1 + A2, the closed-form lepton ratios become coincidence and the framework reduces to "fitted formulas with topological language" — i.e., effectively terminated as a unification claim.

### Even-smaller alternative lethal set (2 falsifiers)

1. **F6: 4th-generation fermion discovered at any energy**
   → kills D_space = 3 → kills K_pair⁴ = 16 → kills almost every closed-form lepton ratio, Higgs mass, hierarchy
2. **F2: w(z) ≠ −1 at >5σ**
   → kills ρ_Λ static interpretation → kills cosmology consistency triangle

Two falsifiers, but F6 alone is **almost a one-shot kill** because K_pair⁴ = 16 is structurally everywhere.

### Single-shot lethal candidates

- **F6 (4th-gen lepton)** — alone, this would force substrate to abandon D = 3 (a foundational axiom). Would reduce the framework to something so different it would no longer be the same theory.
- **F8 (μ-EDM detection at strong CP-violating level)** alone takes out the Möbius Z/2 axiom directly, with no possible workaround — if θ_QCD ≠ 0 at the substrate level, the Z/2 sheet-swap doesn't enforce its fixed point and σ_max = 1/2 loses derivation.
- **F16 (ngEHT-class horizon imaging shows r_s ≠ 2GM/c² at >few %)** would falsify the σ = 1/2 = horizon identification — this also kills A1 if BH horizon is the cleanest anchor for σ_max = 1/2.

---

## 6. Surprising structural observations

Three findings worth flagging as "predictions thought independent but actually linked":

### Observation 1 — Hierarchy and Higgs mass are correlated via M_Pl

Paper 03 (hierarchy at 0.093%) and Paper 10 (m_H at 0.02%) **both** use M_Pl as an anchor scale. They look independent (different observables, different methods), but if M_Pl-of-record were wrong by a few percent, **both** predictions would degrade together. This is the framework's strongest internal consistency check — the fact that M_Pl works at sub-percent for both observables is a 1-in-N² coincidence that is actually structural.

### Observation 2 — F1 (Σm_ν) and F4 (0νββ) test partially-overlapping axioms

Naively F1 and F4 look like two independent tests of the neutrino sector. They are not.

- F1 tests: m_lightest = 2.26 meV anchor + de-saturation cosmology
- F4 tests: Möbius Z/2 forces Majorana + NH spectrum + same m_lightest

The shared variable is m_lightest = 2.26 meV. **F1 and F4 cannot fail in the same direction without invalidating the m_lightest anchor in a way the substrate framework cannot recover from.**

If F1 measures Σm_ν < 30 meV AND F4 finds m_ββ in [5, 21] meV (the substrate's "win" band for 0νββ), the situation is **internally inconsistent for substrate** — m_lightest cannot simultaneously be ~0 (forced by F1) and 2.26 meV (required for the upper-band m_ββ).

This is a hidden cross-correlation: the framework predicts F1 PASS and F4 PASS together, OR F1 FAIL and F4 inconclusive together. It cannot survive F1 fail + F4 specific-band detection.

### Observation 3 — F3 (LiteBIRD r=0) and Σm_ν chain are bundled by de-saturation cosmology

Paper 04 §3.1 derives Σm_ν via the same de-saturation chain that predicts r = 0 (no inflation). If LiteBIRD detects r > 0 (F3), the de-saturation interpretation dies, which means:

- Σm_ν derivation loses its underlying mechanism (becomes a fitted curve)
- ρ_Λ at 0.04% loses the de-saturation context (geometric ratio could survive but the *reason* dies)
- H_0 = 71.92 km/s/Mpc derivation loses its anchor

So **F3 alone is a 3-prediction simultaneous failure cluster**, not a single-prediction failure. This is an unusually high-leverage "null result tested" — LiteBIRD measuring r = 0 is a genuinely strong substrate confirmation, not a routine null.

---

## 7. Summary table — falsifier dependency at a glance

| Falsifier | Direct kill | Cascade kill (1st order) | Cascade kill (2nd order) | Independence rank |
|---|---|---|---|---|
| F1 (Σm_ν < 30) | Σm_ν chain | H_0, σ_8, m_lightest | de-sat cosmology | Medium-coupled |
| F2 (w ≠ −1) | ρ_Λ static | de-sat cosmology, H_0, σ_8 | hierarchy if M_Pl chain breaks | High-coupled |
| F3 (r > 0) | de-sat cosmology | Σm_ν, ρ_Λ interpretation, H_0 | A1 indirectly | High-coupled |
| F4 (m_ββ > 25 OR Dirac) | Möbius Z/2 axiom | σ = 1/2 cap, BH horizon, Pauli, half-spin, 21cm | hierarchy −1, CPT, AB phase | Highest-coupled |
| F6 (4th-gen lepton) | D = 3 | K_pair⁴ = 16, ALL lepton ratios, m_H, hierarchy 4π² | mixing angles | Highest-coupled |
| F7 (SUSY/KK) | none (substrate predicts absence) | none | none | Independent (pass = no propagation) |
| F8 (μ-EDM) | Möbius Z/2 | same as F4 | same as F4 | Highest-coupled |
| F9 (v_GW ≠ c) | ontological unification | photon mass = 0, GW chirp interpretation | none numerical | Medium-coupled |
| F10–F12 (PMNS) | individually mild | K_rank = 5 if multi-angle | n_M = 268 → mass ratios | Medium-coupled if multi |
| F13 (m_H drift) | ε_EW · K_pair⁴ · Λ_QCD | scale-anchor consistency | nothing direct | Medium-coupled |
| F14, F17 (Higgs sector) | minimal-Higgs ontology | none numerical | none | Independent |
| F16 (BH horizon) | σ = 1/2 | same as F4 | same as F4 | Highest-coupled |
| F18 (DUNE δ_CP) | weak; topology not yet fixed | none | none | Independent |
| F19 (KATRIN m_β > 0.45) | NH spectrum | Σm_ν, m_lightest | de-sat cosmology | Medium-coupled |
| F20 (IH confirmed) | NH ansatz | Σm_ν shifts to ~105 meV | partial cosmology | Medium-coupled |

---

## 8. Bottom line

**Three structural takeaways for the framework:**

1. **The Möbius Z/2 axiom is the single most load-bearing element.** It governs σ = 1/2, BH horizon, half-spin, Pauli, Majorana, CPT, 21cm hyperfine, the −1 zero-mode in the hierarchy, and Aharonov-Bohm. F4, F5, F8, F16 all probe this axiom via different observables — and all must pass (or fail) together for the axiom to be intact.

2. **F6 (4th-generation discovery) is the single most catastrophic falsifier.** It kills D = 3, which kills K_pair⁴ = 16, which kills nearly every closed-form lepton/Higgs/hierarchy prediction in the corpus. Single-shot termination.

3. **The cosmology chain (Σm_ν, ρ_Λ, H_0, σ_8, w, r) is a tightly coupled triangle.** The smallest combination that cleanly terminates substrate is {F1 + F2 + F4} — Σm_ν below threshold AND w(z) dynamical AND Majorana ruled out. But F2 + F4 alone is enough to kill both A1 and A2 simultaneously, which would leave only the SM-inheritance derivations (atomic spectra, Stefan-Boltzmann, GW chirp formulas) as substrate's claim — and those are not unification claims.

**The framework's robustness signature:** if substrate is right, all the strong falsifiers should pass together in correlated bundles. If substrate is wrong in a structurally informative way, the failures should also propagate together (e.g., F1 + F2 fail simultaneously because de-saturation cosmology is wrong, OR F4 + F8 + F16 fail simultaneously because Möbius Z/2 is wrong). **Uncorrelated failures across these axiom buckets would be the most diagnostically valuable outcome** — they would falsify the framework AND tell the framework which axiom to revisit.

Conversely, the most diagnostically *uninformative* outcome would be everything passing at sub-percent for 5 more years — which is also the substrate framework's strong hypothesis. The 5–10-year experimental window (DESI DR3 2027, LiteBIRD 2030, LEGEND-1000 2030, JUNO 2027–2030, DUNE 2030+) will resolve most of these in clusters, not one at a time.
