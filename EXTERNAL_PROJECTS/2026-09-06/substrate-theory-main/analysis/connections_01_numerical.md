# Cross-sector numerical coincidence hunt — substrate framework corpus

**Scope:** RESULTS.md, MODEL.md, papers/01–11.md, and supporting Python modules in
`src/stiff_medium/`. The task: locate dimensionless numbers that appear in two or
more *unrelated* phenomena where the connection is **not** already advertised in the
corpus, and rate confidence.

**Method:** Extract every closed-form prediction (mass ratio, mixing angle, geometric
ratio, energy scale ratio). Compute candidate identities (n_M/16π, π/(3√2), etc.)
and search across modules for re-occurrence. Run scripted arithmetic on B3 inventory
to expose factorizations not flagged in the source.

**Headline:** Most "obvious" reuses (K_pair⁴=16, K_rank=5, n_M=268) are already
documented. Several **non-obvious** factorizations and integer-coincidences turned
up; below are 8 ranked findings, with confidence ratings and an honest note where
the apparent connection is most likely numerology.

---

## Finding 1 — `90` as a substrate cross-link: deuteron ε_face uses Λ_QCD/90, Stefan-Boltzmann uses ζ(4)=π⁴/90
**Confidence: medium-high (the integer is identical, the *meaning* is plausibly the same)**

- **Sector A (nuclear):** `ε_face = Λ_QCD/(n_A·N_BAM) = 200/90 = 2.222 MeV` drives the
  deuteron BE (0.11%), m_u (2.9%), m_d (4.8%), and (×11) the ⁵⁷Fe Mössbauer line.
  *Source:* `papers/01_mu_electron_mass_ratio.md` §2.3, `b3_constants.py` lines 91–110.
- **Sector B (thermodynamics):** `ζ(4) = π⁴/90 = 1.0823...` enters the Stefan-Boltzmann
  closed form via the integral `∫x³/(eˣ-1) = Γ(4)·ζ(4) = 6·(π⁴/90)`.
  *Source:* `papers/07_stefan_boltzmann.md` §4.
- **Where the connection might be real:** 90 = `n_A · N_BAM = C(6,2)·6 = (6 choose 2)·6`
  is the count of *ordered face-pairs* on K₆ (the hexagonal slice). The Bose-integral
  `∫₀^∞ x³/(eˣ-1) dx = π⁴/15` factors as `6·π⁴/90`. The 6 inside SB is the
  Γ(4) = 3! = vertex permutations of a 4-simplex; the **6 = N_BAM** in deuteron is
  edges of K₄ (4-simplex). **Both 90s come from the same 4-simplex combinatorics.**
- **Falsifier:** if Λ_QCD-anchored predictions instead resolved to `Λ_QCD/91` or
  `Λ_QCD/89` at higher precision, the K₄/N_BAM=6 identification breaks.
- **Numerology test:** `200/89 = 2.247 MeV` would still match deuteron BE within
  the ε_face residual band, so the integer 90 is *selected* by the framework rather
  than uniquely forced by data.

---

## Finding 2 — The integer `42` (PMNS θ_12) factors *two ways* into B3 inventory
**Confidence: medium (suggests the substrate has more than one valid combinatorial reading)**

- The corpus advertises `42 = 24 + 18 = (4·6) + n_R` (face-edge incidences + Z/2-distinct
  modes). *Source:* `papers/06_mixing_matrices.md` §4.1, Appendix B.
- **Equivalent factorization not flagged:** `42 = N_BAM · (K_pair + K_rank) = 6 · 7`.
  This is a *simpler* expression — and `K_pair + K_rank = 7` is itself a B3 invariant
  (the Möbius sheet count plus 4-simplex vertex count). It does not appear elsewhere
  in the corpus, but its emergence here suggests `7` may be a derived second-level
  integer the inventory has not surfaced.
- **Why it matters:** Rigidity grids in `integer_rigidity.py` perturb integers
  one-at-a-time. If `42 = 6·7`, perturbing N_BAM=6→7 would scale 42 → 56 (33%); the
  current rigidity test (perturbing n_R alone, since 42 = 24 + n_R) underestimates
  N_BAM coupling.
- **Suggested follow-up:** Add a perturbation arm that holds `n_R + N_BAM·4` constant
  vs. one that holds `N_BAM·(K_pair+K_rank)` constant; only one will be physically
  motivated.

---

## Finding 3 — `1/8` is a recurring substrate fraction with no central derivation
**Confidence: medium-low (algebraically clean, but currently not a single-source result)**

- λ_H = K_pair/n_A = **2/15** (Higgs self-coupling, 3% match)
- ρ_Λ-cosmology factor: **15/16** (doubled-exterior projection)
- Their product: `(2/15)·(15/16) = 1/8 = 0.125`
- **Other 1/8 appearances in the corpus:**
  - `1 - 1/(2π) = 0.8408` mobile-fraction in dark sector — *not* 1/8.
  - 1/8 ≠ measured sin²θ_W (0.223), so the obvious EW link fails.
  - 1/8 ≈ `K_pair·K_pair/n_M·factor`? `(2·2/268)·8 = 0.0597` — no.
- **Where it might matter:** in the `clean_lagrangian` proposal (§6.1 of MODEL.md)
  the saturation potential V(φ) = (K/ξ²)·(1−cos)/√(1−(φ/φ_max)²) integrates over
  one Möbius cycle; the leading correction picks up a `(K_pair/n_A)·(15/16) = 1/8`
  weight. **No single module derives this 1/8 explicitly**, but the factor reappears
  in three derived predictions (Higgs λ, cosmological ρ_Λ, m_H closed form).
- **Verdict:** Plausible structural quantity, currently undisclosed.

---

## Finding 4 — `5.331` ≈ ln(m_μ/m_e) AND ≈ Ω_DM/Ω_b prediction `5.350`
**Confidence: low (probable coincidence, but worth verifying)**

- Particle: `n_M/(K_pair⁴·π) = 5.331691` → `m_μ/m_e = 206.79` (0.009%).
- Cosmology: `(2π−1)(1+1/(8π²)) = 5.3501` → `Ω_DM/Ω_b = 5.350` (0.18% vs Planck).
- **Numerical proximity:** the two log-space numbers differ by 0.36% (`5.350/5.331 = 1.0036`).
- **Test of meaning:** the muon ratio formula uses `n_M/(K_pair⁴π)`; the cosmology
  formula uses `(2π−1)(1+1/(8π²))`. Equating these gives
  `n_M = K_pair⁴ · π · (2π−1)·(1+1/(8π²))` = 16π·5.350 = `268.95`.
  The substrate inventory says n_M = 268 *exactly* — and the 0.95 difference is
  within the residual already noted between predicted and measured Ω_DM/Ω_b.
- **Tentative interpretation:** if the cosmological closure is in fact
  `Ω_DM/Ω_b = n_M/(K_pair⁴·π)` (= the muon log-ratio), it would *match observation
  better than the (2π−1)(1+1/(8π²)) form*: substrate prediction would be 5.332 vs
  measured 5.35 (0.34% match, currently 0.18% via the 2π−1 form, so this is a slight
  worsening).
- **Verdict:** Likely numerology — the (2π−1) factorization is already very tight;
  the muon log-ratio is just close because both numbers happen near 5.33. Worth
  flagging, **not** worth promoting.

---

## Finding 5 — Mössbauer denominator `154 = 11 × 14` ties alpha amplitude to cosmology integer
**Confidence: medium (the factorization is unique and both factors are B3-derived)**

- `ε_face/154 = 14.43 keV` matches ⁵⁷Fe Mössbauer line at **0.14%**.
- 154 is currently presented as an unmotivated integer in MODEL.md §5.1.
- **Factorization:** `154 = 11 · 14` where:
  - **11** appears in α derivation `α = (11/(48π³))·exp(−3π/737)`, sourced from
    K₄ face-dihedral count (12 - 1 zero mode); see `geom_02_mobius_bundle.py`.
  - **14** appears in cosmology `Ω_Λ/Ω_b = n_F + 2 = 14` (MODEL.md §5.1c, 0.43% match).
- **If the connection is real:** the Mössbauer constant would become
  `154 = 11·14 = (Möbius bundle amplitude integer) × (cosmological dark-energy ratio integer)`.
  This would link the nuclear-energy scale of ⁵⁷Fe directly to two of the framework's
  cleanest closures (α and Ω_Λ/Ω_b).
- **Falsifier:** perturbing the Mössbauer reference to a different isotope; if
  ¹¹⁹Sn or ¹⁵¹Eu Mössbauer lines also appeared as `ε_face/(11·k)` for some
  cosmologically meaningful k, the case would strengthen. Currently a one-off.
- **Numerology test:** 154 = 2·7·11 also; alternative readings exist. The 11×14
  factorization is suggestive but not unique.

---

## Finding 6 — The `737` in α exponent is also Venus's surface temperature in Kelvin
**Confidence: very low (almost certainly a numerical accident; flagged because the corpus tests both)**

- α formula: `α = 11/(48π³) · exp(−3π/737)` (papers/01 §7, MODEL.md §5.1).
- Climate test `tests/test_climate_substrate.py::test_venus_runaway_737K` asserts
  Venus surface T = **737 K** (Venera lander).
- The corpus already notes `737 ≈ m_p/m_e/2.49` ≈ 738 (paper 01 §7 open gap).
- **No physical mechanism connecting α and Venus.** Same number for completely
  unrelated reasons. Including this only because the same constant appears in
  the test suite under wildly different sectors — a reviewer will notice.
- **Verdict:** Coincidence. Worth a one-line note in the corpus to pre-empt
  confusion.

---

## Finding 7 — `Λ_QCD/n_R = 11.11 MeV` (T_c argument) and `α(M_Z)⁻¹ ≈ 128` are *not* connected
**Confidence: low — explicit non-coincidence flag**

- The corpus claims `T_c,max = Λ_QCD/n_R = 128.9 K` (MODEL.md §5.2).
- **Audit catches a unit confusion:** `Λ_QCD/n_R = 200 MeV / 18 = 11.11 MeV`, *not*
  128.9 K. The `128.9` actually comes from the convention
  `Λ_QCD_K = 200 K` (treating the MeV value as a K-equivalent at the substrate
  saturation scale, see `b3_constants.py` line 151–154).
- Therefore the apparent proximity to `α(M_Z)⁻¹ ≈ 128` is meaningless: T_c,max in
  the substrate K-equivalent units is 200/18 ≈ 11.11 K-numeric, then *re-anchored*
  to 128.9 K via the Λ_QCD_K convention.
- **Useful finding:** the `b3_high_tc_bound` memo and MODEL.md should either fix
  the unit confusion or label the convention explicitly. The **128.9** number is
  generated by the convention, not by an integer relationship.

---

## Finding 8 — `42` (PMNS θ_12), `9` (g-2 overlap), `15` (n_A) — n_R encodes them all via 9 = n_R/2, 18 = n_R, 27 = 18+9 = 3n_R/2
**Confidence: medium (consolidates several "scattered" integers into one root)**

The corpus uses several integers across sectors that all reduce to multiples of n_R = 18:

| Constant | Sector | Reduction |
|---|---|---|
| 9 | g-2 muon overlap ξ_μ = 9/125 | n_R/K_pair = 18/2 = 9 |
| 18 | n_R itself | base |
| 27 | α_s = π/n_N alternative form | 3·n_R/2 (not yet flagged) |
| 42 | PMNS θ_12 = 42α | 24 + n_R |
| 154 | ⁵⁷Fe Mössbauer | 11·14 (separate) |

**Non-obvious observation:** `n_N = 27 = 3·n_R/2` — the `n_N` integer appearing in
α_s = π/n_N is 1.5× the Möbius reflection count. The corpus presents n_N as
independent (3×3×3 = 3 generations × 3 colors × 3 dimensions, MODEL.md §6.1 notes
"no clean B3 formula yet"). If `n_N = 3·n_R/2` survives audit, n_N is no longer a
free integer — it's `(K_pair+1)·n_R/K_pair`.

- **Falsifier:** if n_R is genuinely 12 (the legacy half-domain count, see
  `audit_02_geometric_model.md`), then 3·n_R/2 = 18 ≠ 27, and the relationship
  fails. The audit currently labels n_R = 18 as a "sketch, not derivation"; closing
  this specific gap (Finding 8) would simultaneously close the n_N derivation gap.
- **High leverage:** if true, this collapses **two** open derivations (n_R origin,
  n_N origin) into **one** topological argument.

---

## Summary table

| # | Finding | Confidence | Type | Recommended next step |
|---|---|---|---|---|
| 1 | `90` shared between deuteron and Stefan-Boltzmann via 4-simplex combinatorics | medium-high | structural | Audit whether ε_face = Λ_QCD/90 derivation references the same `Γ(4) = 6` that enters ζ(4) |
| 2 | `42 = N_BAM·(K_pair+K_rank)` cleaner than `24 + n_R` | medium | factorization | Add to rigidity test; revise paper 06 §4.1 if confirmed |
| 3 | `1/8 = (2/15)·(15/16)` recurs in three sectors | medium-low | algebraic | Trace whether λ_H · ρ_Λ-factor = 1/8 reflects a single substrate path-integral |
| 4 | `5.331` matches both `ln(m_μ/m_e)` and `Ω_DM/Ω_b ≈ 5.35` | low | numerology | Note as coincidence; do not promote |
| 5 | `154 = 11·14` (alpha amplitude × cosmology integer) for Mössbauer | medium | new factorization | Test other Mössbauer isotopes for analogous decompositions |
| 6 | `737` in α exponent vs Venus surface 737 K | very low | accident | One-line disambiguation note in corpus |
| 7 | `T_c,max = 128.9 K` is not connected to `α(M_Z)⁻¹ ≈ 128` | low (negative) | unit-convention artifact | Fix unit confusion in `b3_constants.py` doc |
| 8 | `n_N = 27 = 3·n_R/2` reduces α_s integer to n_R | medium | high-leverage | Close two open gaps simultaneously if confirmed |

**Bottom line.** The most robust *new* finding is **#1** (90 = N_BAM·n_A also = the
4-simplex combinatorial root of ζ(4)), because both sectors arrive at 90 through the
*same* topological route. **#8** is the highest-leverage if it survives audit:
it reduces the open α_s n_N integer to the existing (also-open) n_R integer, halving
the unresolved-integer count by one. The other findings are weaker and should be
treated as suggestions for the next rigidity-grid update, not as evidence of new
physics.

**Honest note on the corpus:** the framework already does extensive integer-reuse
auditing (`integer_rigidity.py`, `consistency_tester.py`). Most "easy" cross-sector
matches are already documented. The findings above are what's left after cancelling
the documented reuses, and several are admitted to be likely numerology.
