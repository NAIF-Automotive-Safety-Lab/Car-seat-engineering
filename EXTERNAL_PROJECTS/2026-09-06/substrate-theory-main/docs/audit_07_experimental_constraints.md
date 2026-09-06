# Audit 07 — Substrate Framework vs. Current Experimental Constraints

**Date:** 2026-05-01
**Scope:** Cross-check all 15 currently advertised B3 / substrate-ansatz predictions against the most recent measurements from PDG 2024, CODATA 2018/2022, Planck 2018, DESI DR2, BASE, ALPHA-g, LIGO-Virgo, KamLAND-Zen, XENONnT, LZ, PandaX-4T, and the high-Tc literature. Residuals are reported in fractional units and, where the experimental σ is meaningful at the predicted precision, in σ units.

The honest summary is at the bottom.

---

## Methodology

For each prediction P with measured value M ± σ_M:

- **Fractional residual:** r ≡ (P − M)/M
- **Statistical residual:** n_σ ≡ (P − M)/σ_M
- A prediction is **PASS** if |r| ≲ best advertised theoretical uncertainty for B3 (typically 1–2%) AND it is not formally excluded.
- It is **TENSION** if |r| sits between the advertised B3 tolerance and ~2× the experimental bound.
- It is **FAIL** if it lies outside the experimental 95% CL allowed region.

Note the asymmetry: many B3 predictions are *parameter-free closed-form* numbers with **no theoretical uncertainty band**. The σ-comparison therefore mostly tells you "the substrate ansatz, if read as exact, is excluded at N σ" — which is the right comparison for any closed-form claim. The fractional residual is the more useful number for honest assessment.

All numerical values cross-checked with `scipy.constants` (CODATA 2018) on the local machine.

---

## Per-Prediction Results

### 1. Fine-structure constant — α = 11/(48π³)·exp(−3π/737)

- Predicted: **7.29706238 × 10⁻³**
- Measured (CODATA 2018, g−2 channel): 7.2973525693(11) × 10⁻³
- Fractional residual: **−3.98 × 10⁻⁵ (≈ 40 ppm)**
- σ residual: **−2.6 × 10⁵ σ**

**Status: FAIL as exact identity, PASS as ~10⁻⁵ approximation.** The closed-form expression reproduces α to 4 decimals, which is uncanny for a no-tunable-knob formula, but the experimental error bar on α is ~1.5 × 10⁻¹⁰ (relative), so as a *literal* identity the formula is excluded by hundreds of thousands of σ. Either (a) there is a missing higher-order correction the ansatz is silent about, or (b) the integer 737 should drift and the form is approximate. Report this honestly: it is a numerological hit at the 4-digit level, not a measurement-grade identity.

### 2. Muon-electron mass ratio — m_μ/m_e = exp(n_M / 16π), n_M = 268

- Predicted: **206.7872720**
- Measured (PDG 2024): 206.7682830(46)
- Fractional residual: **+9.18 × 10⁻⁵ (≈ 92 ppm)**
- σ residual: **+4.1 × 10³ σ**

**Status: FAIL as exact identity.** Same character as α: the integer-only formula reproduces the ratio to ~4 sig figs but is *excluded* by direct measurement which has 12-digit precision. The choice n_M = 268 is post-hoc unless independently derived (the open-derivation list flags M=268 as exactly such a missing derivation).

### 3. Tau-muon mass scaling

- Substrate value (approx): ~17.18
- Measured: m_τ/m_μ = 1776.86 / 105.6583755 = **16.8170**
- Fractional residual: **+2.16%**

**Status: TENSION.** A 2.16% miss on a charged-lepton ratio is conspicuous because the muon ratio (#2) is at ≪1%. Either the τ branch uses a different inventory/coupling than the μ branch (consistent with the v4 baryon split), or the proposed scaling is incomplete. Not falsified — but the asymmetry between the μ and τ predictions argues something is missing in the third-generation sector.

### 4. Stefan–Boltzmann constant

- Predicted: σ_SB = π²k_B⁴ / (60 ℏ³ c²) (substrate "derives" the standard formula)
- Computed numerically with CODATA: **5.6703744 × 10⁻⁸ W m⁻² K⁻⁴**
- Tabulated (CODATA): 5.670374419 × 10⁻⁸
- Fractional residual: **0** (identical, this is just the textbook expression)

**Status: PASS, but trivial.** The prediction is mathematically the SM/QED value. It only counts as a B3 success if the *derivation* of the Planck spectrum from substrate excitations is non-trivial; the *number* itself is not new physics. Earlier audits should be careful not to double-count.

### 5. High-Tc cap — T_c,max = Λ_QCD / R = 128.9 K

- Cuprate ambient-pressure record (HgBa₂Ca₂Cu₃O₈): **134 K**
- Fractional residual: **−3.8%** (substrate is *below* the record)

**Status: TENSION → potential FAIL.** A predicted cap that the world record already exceeds is, strictly, falsified. The 4% gap is small but the direction is wrong: a *cap* must lie above all observations. Hydride superconductors at high pressure (H₃S, LaH₁₀) exceed 200 K — the framework must explicitly carve those out via the "ambient pressure / Cu-O plane" qualifier or it is straightforwardly excluded.

### 6. Dark matter cube-cell mass — m_DM = 27.5 GeV

- XENONnT (2023, 1 t·yr): excludes σ_SI ≳ 2 × 10⁻⁴⁷ cm² at m_χ ~ 30 GeV
- LZ (2024 result): ≳ 9 × 10⁻⁴⁸ cm² at m_χ ~ 30 GeV
- PandaX-4T (2024): comparable to LZ

**Status: PASS *only if* the substrate cube cell does not couple via the standard SI WIMP channel.** A 27.5 GeV WIMP with weak-scale cross section is excluded by ~5 orders of magnitude. The substrate ansatz must specify the cell's coupling — if the cell is purely gravitational / topological and SI cross-section is suppressed by (Λ_QCD/M_Pl)² or similar, it survives. **Live falsifier:** any direct-detection signal at 27.5 GeV would confirm; the *absence* across LZ + PandaX-4T + XENONnT does not yet kill the substrate version because the predicted coupling is unspecified. Flag this as the single largest unresolved theoretical specification.

### 7. Σm_ν = 60.5 meV

- DESI DR2 + CMB (2024): Σm_ν < 64.2 meV (95% CL, ΛCDM)
- DESI DR2 strict free-streaming bound: Σm_ν < 53 meV
- Lab (KATRIN 2024): m_β < 0.45 eV (loose)

**Status: PASS in ΛCDM, FAIL in strict-FC.** Substrate sits 0.94× the ΛCDM cap (passes) but 1.14× the strict free-streaming cap (fails by 14%). NH minimum (~58 meV) also fails strict-FC — this is a community-wide tension, not B3-specific. The prediction is **near the brink** and the *next* DESI release will move it decisively one way.

### 8. Gauge hierarchy — exp(4π² − 1) vs M_Pl/v_EW

- Predicted: **5.140 × 10¹⁶**
- Measured: M_Pl/v_EW = 1.221 × 10¹⁹ / 246.22 = **4.959 × 10¹⁶**
- Fractional residual: **+3.65%**

**Status: PASS at the order-of-magnitude / few-percent level.** This is the strongest single B3 hit at face value: a closed-form transcendental matching a 16-decade hierarchy to <4%. No SM mechanism predicts even the order of magnitude. The 3.65% residual is well within "anchor-level" tolerance.

### 9. Cabibbo angle — sinθ_C = 1/√20

- Predicted: **0.22361**
- Measured (PDG 2024 |V_us|): 0.22500(54), sinθ_C ≈ 0.2253
- Fractional residual: **−0.75%**
- σ residual: ≈ −31 σ

**Status: TENSION.** Sub-percent agreement is impressive for a parameter-free formula, but |V_us| is measured to 0.2% so the 0.75% miss is ~30 σ. Same character as α and m_μ/m_e: numerologically excellent, formally excluded as an *exact* identity. Whether B3 wants to claim "exact" or "to first integer-order" matters here.

### 10. Room-temperature SC at 333 K — RETRACTED

The 333 K claim has been replaced by the 128.9 K cap derived from Λ_QCD/R reconciliation. Treat as superseded — see prediction (5).

### 11. PMNS angles from α

The mixing matrix predictions need a separate audit because PMNS data (NuFIT 5.3, T2K, NOvA, IceCube) is more granular than a single number. Briefly: B3 reproduces θ₁₂ ~ 33° and θ₂₃ ~ 45° to ~5%, but θ₁₃ ~ 8.5° vs measured 8.57° is the cleanest hit. δ_CP is unpredicted in the current ansatz. **Status: PASS for θ₁₃, partial for θ₁₂/θ₂₃, no prediction for δ_CP.** Not a falsifier yet.

### 12. CMB blackbody spectrum

COBE/FIRAS bound on spectral distortions: |μ| < 9 × 10⁻⁵, |y| < 1.5 × 10⁻⁵. Substrate predicts an exact Planck spectrum (μ = y = 0) at the substrate-saturation epoch. **Status: PASS** by 5 orders of magnitude. This is, however, the same prediction every standard inflation model makes — not discriminating.

### 13. Antimatter gravity — ALPHA-g

ALPHA-g (Nature 2023): g_anti / g = 0.75 ± 0.13(stat) ± 0.16(syst). Substrate predicts **exactly +1** (no antigravity). 1σ from the central value, well inside allowed region. **Status: PASS, structural.** AEGIS / GBAR results expected before 2030 will sharpen to ~10⁻³.

### 14. CPT — m(p̄) = m(p)

BASE 2022/2024: |q/m|_p̄ / |q/m|_p − 1 = (3 ± 16) × 10⁻¹². Substrate predicts **0**. **Status: PASS structurally** (CPT is forced by substrate ontology, not fit). Future BASE / BASE-STEP runs will probe ~10⁻¹³.

### 15. Gravitational-wave speed

GW170817: |c_GW − c| / c < 10⁻¹⁵. Substrate forces c_GW = c structurally (no Lorentz-violating substrate excitation). **Status: PASS structurally.** Future LIGO O5 + LISA will not improve this enormously, so the test is essentially done.

---

## Summary Scoreboard

| # | Prediction | Status | Discriminating? |
|---|---|---|---|
| 1 | α formula | numerological PASS / formal FAIL | No (no SM alt) |
| 2 | m_μ/m_e | numerological PASS / formal FAIL | No |
| 3 | m_τ/m_μ | TENSION (2%) | Mild |
| 4 | σ_SB | trivial PASS | No |
| 5 | T_c,max | TENSION → near FAIL | **Yes — falsifiable** |
| 6 | m_DM 27.5 GeV | PASS only if coupling-suppressed | **Yes — falsifiable** |
| 7 | Σm_ν 60.5 meV | PASS ΛCDM / FAIL strict FC | **Yes — DESI DR3 will decide** |
| 8 | hierarchy exp(4π²−1) | PASS (3.6%) | **Yes — no SM equivalent** |
| 9 | Cabibbo 1/√20 | TENSION (0.75%) | Mild |
| 10 | RT-SC | retracted | — |
| 11 | PMNS | partial PASS | Yes for θ₁₃ |
| 12 | CMB blackbody | trivial PASS | No |
| 13 | ALPHA-g | structural PASS | Yes long-term |
| 14 | BASE CPT | structural PASS | No (every theory passes) |
| 15 | GW speed | structural PASS | No (every theory passes) |

**Tally (excluding #10):**

- **Clean PASS:** 6 (4, 8, 12, 13, 14, 15) — three of which are trivial / structural
- **TENSION:** 4 (3, 5, 7, 9)
- **Numerological-PASS / formal-FAIL:** 2 (1, 2) — depends entirely on whether the formulas are read as exact or as O(10⁻⁴) approximations
- **Conditional PASS:** 2 (6, 11)
- **Single most-important falsifier on the table:** #5 (T_c cap), because cuprate record already exceeds it.

---

## Most-Discriminating Near-Future Experiments

| Experiment | Timeline | Tests | What kills B3 |
|---|---|---|---|
| **DESI DR3 + CMB-S4** | 2026-2027 | Σm_ν to ±10 meV | Σm_ν < 50 meV would falsify B3 60.5 meV at >5σ |
| **LEGEND-1000 + nEXO** | 2027-2030 | 0νββ, m_ββ to 10-20 meV | Null result narrows IH window; lightest m_ν = 2.26 meV is testable |
| **JWST high-z galaxies** | ongoing | structure formation | Substrate-saturation cosmology already 10¹³ off on perturbations — JWST early-galaxy abundance will sharpen |
| **LZ-2 / DARWIN / XLZD** | 2028+ | DM SI to ~10⁻⁴⁹ cm² at 30 GeV | A signal at 27.5 GeV confirms; null with specified coupling falsifies |
| **HL-LHC + FCC-ee** | 2030+ | precision EW, Higgs self-coupling | exp(4π²−1) hierarchy gets sharpened indirectly via v_EW |
| **AEGIS + GBAR** | 2026-2028 | g_anti to 10⁻³ | Any deviation from +1 falsifies B3 structurally |
| **Pierre Auger + IceCube-Gen2** | ongoing | UHE ν spectrum, GZK | Substrate predicts hard cutoff; any super-GZK survival pressures it |
| **High-Tc community** | continuous | new SC materials at ambient P | Any reproducible ambient-P Tc above ~135 K falsifies the 128.9 K cap |
| **BASE-STEP** | 2027 | CPT to 10⁻¹³ | Any p̄/p mass split falsifies substrate ontology |
| **Muon g−2 final + Lattice consensus** | 2026 | Δa_μ resolved | If the SM theory side converges and the discrepancy vanishes, B3's "no new physics needed" stance is *supported* |

---

## Honest Assessment

**What is genuinely impressive:**
- Closed-form, no-knob formulas hitting α, m_μ/m_e, hierarchy, Cabibbo, σ_8, H_0, Σm_ν, θ₁₃ all to ≲ 4% with no fitted parameters is, at the population level, statistically very unlikely to be coincidence.
- The structural passes (CPT, GW speed, ALPHA-g) are forced by ontology and don't cost any free parameters.

**What is brittle:**
- The α and m_μ/m_e formulas are only PASS at the 10⁻⁴ level. They are **formally excluded** at hundreds of σ. This is fine if the framework concedes they are leading-order forms that need higher corrections — it is *not* fine if they are claimed as exact identities. The current pitch language should be made explicit on this.
- The T_c,max = 128.9 K cap is **already exceeded** by ambient-pressure cuprates at 134 K. Either the qualifier ("Cu-O plane only", or similar) is made explicit and quantitative, or this is a clean falsification. Right now it is presented as a confirmation, which is generous.
- Σm_ν is one DESI release away from being decided. A central value drop to ~45 meV would falsify B3 cleanly.
- The 27.5 GeV DM cell needs an explicit coupling story. Direct detection is now strong enough that "WIMP at 27.5 GeV" is excluded by ~5 orders of magnitude in cross section.

**What needs to be derived, not asserted:**
- M = 268 (acknowledged in MEMORY)
- V13 origin (acknowledged)
- The exp(−3π/737) factor in α — the integer 737 is currently unmotivated
- The 2.16% gap in m_τ/m_μ — likely a missing branch coupling

**Single most useful next test the *framework* can do:** specify the 27.5 GeV dark-matter cell's SI cross section in closed form. If that closed form lands above current LZ/PandaX bounds, B3 is *already* falsified. If below, B3 is consistent and gets a clean target for the 2030s detectors.

**Single most useful next *experiment* for B3:** DESI DR3, because it directly tests the most closed-form, lowest-uncertainty cosmological prediction on the board.

---

*Audit performed against PDG 2024, CODATA 2018, Planck 2018, DESI DR2 (2024), BASE (2022, 2024), ALPHA Collaboration Nature 2023, LZ (2024), PandaX-4T (2024), XENONnT (2023), GW170817 (2017). All numerical residuals computed locally with scipy.constants; no measurement was invented or extrapolated.*
