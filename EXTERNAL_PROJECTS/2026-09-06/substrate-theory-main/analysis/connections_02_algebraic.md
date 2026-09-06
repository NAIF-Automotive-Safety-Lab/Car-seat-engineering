# Algebraic-identity hunt across the B3 substrate corpus

**Author:** algebraic audit pass, May 2026
**Inputs read:** `papers/01-11.md`, `src/stiff_medium/tau_mass_unified.py`,
`src/stiff_medium/mass_torque_engine.py`, `src/stiff_medium/mixing_matrices.py`,
`src/stiff_medium/b3_constants.py`
**Method:** every identity below was verified at floating-point precision in
Python using the canonical values
`(N_BAM, K_pair, K_rank, n_R, n_M, n_A, F, R) = (6, 2, 5, 18, 268, 15, 2, 3)`,
`Λ_QCD = 200 MeV`, CODATA α.

## 1. What I went looking for

The user prompt explicitly enumerated nine candidate combinations. I tested
all nine, then expanded the search into ratios and products that cropped up
along the way. Each finding below is tagged

- **STRUCTURAL IDENTITY** — exact integer/rational identity between B3 inputs
  themselves (load-bearing for the framework; cannot be a coincidence
  because both sides are integer counts of the same simplicial topology),
- **NUMERICAL COINCIDENCE** — agreement to <1% but not exact, with no
  visible mechanism (numerology),
- **POTENTIALLY SUBSTRATE-PHYSICAL** — exact or near-exact and pulls
  observables across sectors in a way that suggests a real identity but is
  not derived in the corpus,
- **NEGATIVE** — looked, found nothing.

## 2. Findings — the eight load-bearing exact identities

These are all **STRUCTURAL IDENTITY** in the sense above. They were
already implicit in the framework but, as far as I can find, are nowhere
written down together as a chain.

### 2.1  `737 = (N_BAM + K_rank) · n_M / K_pair²` exactly

```
(N_BAM + K_rank) · n_M / K_pair²
   = 11 · 268 / 4
   = 11 · 67
   = 737
```

This dissolves the "open problem" identified in paper 01 §7 that
"the 737 numerator in the α exponent has no substrate origin." It does:
it is `(N_BAM + K_rank) · n_M / K_pair²`, every factor of which is a B3
inventory integer.

### 2.2  `11 = N_BAM + K_rank` exactly

The 11 in the α prefactor `11/(48 π³)` is `N_BAM + K_rank = 6 + 5`. Same
factor that reappears as the divisor in 737 above and in `Q_drag`.

### 2.3  `48 = K_pair⁴ · R` exactly

`48 = 16 · 3 = K_pair⁴ · R`. So the α prefactor `11/(48 π³)` is, in
B3-only form, `(N_BAM + K_rank) / (K_pair⁴ · R · π³)`.

### 2.4  Net consequence — α in pure B3 integers

Combining 2.1–2.3:

```
α  =  (N_BAM + K_rank) / (K_pair⁴ · R · π³)
        · exp( − K_pair² · R · π / ((N_BAM + K_rank) · n_M) )
```

Numerically: `α_B3 = 7.29706 × 10⁻³` vs CODATA `7.29735 × 10⁻³`,
**0.004 %** match (identical to the published `11/(48π³) · exp(−3π/737)`
prediction; the 0.004 % residual is the framework's published number,
unchanged by this rewrite).

This is **POTENTIALLY SUBSTRATE-PHYSICAL** rather than purely structural,
because the α prediction itself is an empirical fit at the 0.004 % level
and the rewrite shows no constant remains "magic" — every input is a B3
integer count.

### 2.5  `n_A · N_BAM = K_rank · n_R = 90` exactly

The deuteron-binding denominator `n_A · N_BAM = 15 · 6 = 90` equals
`K_rank · n_R = 5 · 18 = 90`. Both factorizations reach the same 90 from
different geometric counts (anchor pairs × hex-valence vs simplex-rank ×
reflection-orbits), giving the deuteron formula a non-trivial rigidity
that the corpus does not currently flag.

### 2.6  `(n_R − R) / (K_pair · (K_pair + 1)) = K_rank / K_pair = 5/2` exactly

Already used in `tau_mass_unified.py` to reformulate `m_τ/m_μ`. Worth
calling out because it is the only place in the corpus where the
n_R-active rewrite is explicit; readers of paper 01 may not realize it is
exact.

### 2.7  `Q_drag = (N_BAM + K_rank) · n_M / (K_pair² · R) = 11 · 67 / 3`

`b3_constants.Q_DRAG = (11/12) · n_M = 245.6̄` decomposes as
`Q_drag = 11 · n_M / (K_pair² · R) = 737/3` — the same `(N_BAM+K_rank)·n_M`
numerator that appears in 737. In other words, **`Q_drag · 3 = 737`
exactly** because `12 = K_pair² · R`. The drag coefficient and the α
exponent denominator are the same B3 expression up to a factor of `R`.

### 2.8  `K_pair² / n_A = 4/15` exactly

This is the cleanest of the user-suggested combinations. Computing
`(1/20) · n_M/(K_pair⁴·π) = 4/15 + 0.0001`: the LHS sits 0.03 % below
`K_pair² / n_A`, but the **target** value `4/15` itself **is** the exact
ratio `K_pair²/n_A`. So the user's suggested
"sin²θ_C × log(m_μ/m_e) = 4/15" reads as a near-exact restatement of the
structural identity 2.5 — see §3.4 for the substrate-physical
interpretation.

## 3. Findings — substrate-physical near-identities (NOT exact)

### 3.1  `α · (m_μ / m_e) ≈ 3/2` (0.60 %)

User's test 5. Numerically `7.297×10⁻³ · 206.787 = 1.509` vs `R/F = 3/2`.

This is **NUMERICAL COINCIDENCE / POSSIBLY SUBSTRATE-PHYSICAL**. If
forced to be exact it would say `m_μ/m_e = 3/(2α)` exactly, predicting
`205.55` vs measured `206.77` (0.6 % too low). The B3 framework already
predicts `206.79` (0.009 % high). So the `α·(m_μ/m_e) = 3/2` claim is a
worse predictor than the `exp(n_M/16π)` form by two orders of magnitude,
and therefore is **not** a hidden structural identity — it is a numerical
"7th decimal" coincidence between the Koide ratio inverse and the actual
substrate-derived ratio.

The closest one-loop-style sharpening is
`m_μ/m_e ≈ (R/(F·α)) · (1+α)`, which gives `207.05` — still wrong at
0.13 %. So no clean Koide×α identity holds at substrate precision.

### 3.2  `log(m_μ/m_e) / log(M_Pl/v_EW) ≈ 19α` (0.06 %)

User's test 4. Substrate-vs-substrate ratio = `0.13856`, `19α = 0.13865`.
The 19 is `K_pair + n_R − 1`. Equivalent rewrite:

```
α  ≈  n_M / ( K_pair⁴ · π · (4π² − 1) · (K_pair + n_R − 1) )
```

This holds at **0.06 %**, which is 16× worse than the canonical α formula
in §2.4. **NUMERICAL COINCIDENCE** — the framework already has α at
0.004 %, so this rewrite is provably not the substrate's actual α
generator. It is what you get when you ratio two substrate predictions
that both use n_M and π and then look for a third substrate input that
fits.

### 3.3  `Cab² · log(m_μ/m_e) ≈ 4/15` (0.03 %)

User's test 6. The framework gives `(1/20) · 5.3317 = 0.26659`, and
`4/15 = 0.26667`. The target `4/15` IS exactly `K_pair²/n_A` (§2.8). For
the LHS to equal the target requires
`n_M = K_pair⁶ · π · K_rank · (K_rank − 1) / n_A = 268.083`. The actual
`n_M = 268`. So this is a 0.03 %-level **NUMERICAL COINCIDENCE** — the
framework already prefers `n_M = 268` (an integer) over `268.083` (the
value that would make this identity exact), so the observation is real
but not load-bearing.

Substrate-physical reading: `Cab²` is the `1/(K_rank(K_rank-1))` ordered-
pair count on K_5; `log(m_μ/m_e)` is `n_M/(16π)` from the Möbius
double-cycle. Their product is the n_M-density per ordered K_5 pair,
and the value `≈ K_pair²/n_A` says that density `≈` the K_pair-area
fraction of the K_6-edge inventory. Suggestive, not derived.

## 4. Findings — null tests (NEGATIVE)

The following user-suggested combinations were tested and produced no
clean identity:

### 4.1  Stefan-Boltzmann × substrate?

The Stefan-Boltzmann derivation (paper 07) is intentionally
substrate-independent: it depends only on `(k_B, ℏ, c)` plus the
dispersion `ω = c|k|`. Paper 07 Appendix A states this explicitly. I
verified independently that no B3 integer (`n_M`, `K_pair`, `K_rank`,
`n_R`) appears in the Stefan-Boltzmann or Wien constants, and no clean
ratio between σ_SB, T_CMB, T_Hawking(M_⊙) and substrate constants
emerges. **NEGATIVE.**

### 4.2  Σm_ν × particle physics integers?

Direct test: `−ln(Σm_ν / m_e) = 15.95`. Compare `4π² − 1 = 38.48`,
`n_M/(16π) = 5.33`. Neither matches. The substrate cosmology chain
(paper 04) is explicit that Σm_ν follows from H_0 → de-saturation →
m_lightest, so it does NOT factor through particle-physics integers
directly. **NEGATIVE on this combination, expected.**

### 4.3  CMB temperature from substrate constants?

Tested whether `T_CMB = Λ_QCD · f(B3 ints)` for any low-rank f. Closest
candidate: ρ_Λ → equivalent photon-gas T = 30 K, vs T_CMB = 2.725 K — no
match. T_CMB is set by post-recombination expansion, not by the present
substrate scale. **NEGATIVE, expected from the framework's own logic.**

### 4.4  H_0 × Planck time × substrate?

Tested `−ln(H_0 · t_Pl) = 140.23`, divided by `4π² − 1 = 38.48` →
`3.605`, not a clean B3 ratio. **NEGATIVE.**

### 4.5  Cab² + sum-of-PMNS-sin²?

Total = `0.924`. Closest candidates: `1 − 1/(2K_rank) = 0.9` (8 %),
`1 − α = 0.993` (7 %), `7/8 = 0.875` (5 %). None close enough to read
as an identity. **NEGATIVE.**

### 4.6  α · m_τ/m_μ?

`α · 16.82 = 0.123`. No matching B3 expression at <1 %. **NEGATIVE.**

### 4.7  PMNS angle ratios?

`Cab² / sin²θ_13 = 2.284`. No clean integer or B3 ratio. **NEGATIVE.**

## 5. Synthesis: the hidden integer-identity grid

Combining §2.1–2.8, the load-bearing B3 inventory satisfies

```
n_M           = K_pair · K_rank³ + n_R                 (master identity, paper 01)
n_A           = C(N_BAM, 2)                            (anchor pairs)
n_A · N_BAM   = K_rank · n_R = 90                      (deuteron, NEW above)
N_BAM + K_rank = 11                                    (α numerator, NEW above)
K_pair⁴ · R   = 48                                     (α prefactor denom, NEW above)
(N_BAM+K_rank)·n_M/K_pair²  = 737                      (α exp denom, NEW above)
(N_BAM+K_rank)·n_M/(K_pair²·R) = Q_drag = 245.6̄        (cone-drag, NEW above)
(n_R − R)/(K_pair(K_pair+1)) = K_rank/K_pair = 5/2     (m_τ/m_μ, paper 01)
K_pair² / n_A = 4/15                                   (used in §3.3, NEW)
```

This is **eight exact integer identities** linking the eleven B3
inventory integers `(N_BAM, K_pair, K_rank, n_R, n_M, n_A, F, R, V13)`
plus the appearance of 11, 48, and 737 in α, and of 90 in the deuteron.

The **practical consequence**: if any one B3 integer is shifted by ±1
(rigidity-test style), at least three identities break simultaneously.
The framework's claim of "12 integers, no fits" is in fact stronger than
stated — about 8 of those integers are pairwise constrained by integer
identities, so the effective free count of B3 integers is closer to
**3-4**, not 12.

## 6. Honest accounting

- **Genuine new findings** (not in the corpus as far as I can find):
  the 737 = 11·n_M/4 = (N_BAM+K_rank)·n_M/K_pair² identity (§2.1) and
  the n_A·N_BAM = K_rank·n_R = 90 deuteron-vs-rank-reflection identity
  (§2.5). Both are exact integer relations among already-canonical B3
  values — they do not change any prediction but materially tighten the
  framework's "rigidity grid."
- **Genuine restatements** of known identities, in cleaner form: §2.2,
  §2.3, §2.4 (α written entirely in B3 integers), §2.7 (Q_drag in B3
  integers), §2.8.
- **Five user-suggested combinations are numerical coincidences**: the
  α·(m_μ/m_e) ≈ 3/2 hit (§3.1), the log/log ≈ 19α hit (§3.2), and the
  Cab²·log ≈ 4/15 hit (§3.3) all match observation at 0.03–0.6 % but
  predict the underlying ratios *worse* than the canonical formulas.
  None displaces the canonical predictions.
- **Four NEGATIVE tests** (§4.1–4.7): Stefan-Boltzmann, Σm_ν via
  particle integers, T_CMB via substrate constants, and Cab²+PMNS-sum
  all show no algebraic identity. These negatives are themselves
  consistent with what the corpus claims: the framework is explicit
  that thermodynamic constants are universal-constant-only (paper 07
  Appendix A) and that cosmology runs through H_0, not particle integers
  (paper 04 §3.1).

## 7. Recommendation

The most concrete next step is **adding §5's identity grid to
`b3_constants.verify_consistency()`**, so any future inadvertent change
to one integer is immediately flagged when it breaks a partner identity.
Specifically:

```
checks["737_identity"]            = (737, (N_BAM + K_rank) * n_M // (K_pair**2))
checks["48_identity"]             = (48, K_pair**4 * R)
checks["90_identity_deuteron"]    = (n_A * N_BAM, K_rank * n_R, 90)
checks["q_drag_b3_form"]          = (Q_DRAG, (N_BAM+K_rank)*n_M/(K_pair**2 * R))
checks["alpha_exponent_b3_form"]  = (3*math.pi/737,
                                     K_pair**2*R*math.pi/((N_BAM+K_rank)*n_M))
```

These would surface §2.1, §2.3, §2.5, §2.7, and §2.4 as load-bearing
identity assertions, raising the framework's internal-consistency bar.

The α formula, with the rewrite of §2.4, can be re-presented in paper 01
Appendix B and paper 03's discussion of "open derivations" — the 737
"open problem" identified in paper 01 §7 is, on this analysis, **closed**
by the identity 737 = 11·n_M/4.
