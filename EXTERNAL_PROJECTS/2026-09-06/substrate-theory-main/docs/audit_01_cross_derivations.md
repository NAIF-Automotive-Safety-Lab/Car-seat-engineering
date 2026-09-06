# Audit 01 — Cross-Derivation Internal Consistency

**Date:** 2026-05-01
**Scope:** Substrate framework's "fully end-to-end" derivations and their
consistency across the modules in `src/stiff_medium/`.
**Method:** Each named configuration was actually executed (not just read)
through the live code, observed values pulled from the engine's PDG/CODATA
constants, and the underlying integers/primitives traced across dependent
modules.

---

## 1. Per-derivation numerical results

The `MassTorque` engine's bundled `report()` was run verbatim. All nine
registered configurations pass within their declared tolerance:

| name           | predicted        | observed        | rel_err   | tol     | pass |
|----------------|------------------|------------------|-----------|---------|------|
| deuteron       | 2.22222 MeV      | 2.22457 MeV      | 1.06e-3   | 2.0e-3  | Y    |
| alpha (BE)     | 28.3333 MeV      | 28.295 MeV       | 1.36e-3   | 5.0e-2  | Y    |
| electron       | 0.510999 MeV     | 0.510999 MeV     | 0 (anchor)| 1.0e-9  | Y    |
| muon           | 105.668 MeV      | 105.658 MeV      | 9.18e-5   | 5.0e-3  | Y    |
| tau            | 1793.62 MeV      | 1776.86 MeV      | 9.43e-3   | 5.0e-2  | Y    |
| higgs          | 123 000 MeV      | 125 250 MeV      | 1.80e-2   | 5.0e-2  | Y    |
| hierarchy      | 5.140e16         | 4.963e16         | 3.57e-2   | 5.0e-2  | Y    |
| fine_structure | 7.29706e-3       | 7.29735e-3       | 3.98e-5   | 1.0e-2  | Y    |
| t_c_max        | 128.9 K          | 128.9 K          | 0         | 1.0e-2  | Y    |

### Independent recomputation of the headline closed forms

Recomputed from scratch, no shared module state:

| Closed form                                  | value          | reference         | residual |
|----------------------------------------------|----------------|-------------------|----------|
| α = 11/(48π³)·exp(−3π/737)                   | 1/137.04145    | CODATA 1/137.03600| 3.98e-5  |
| m_μ/m_e = exp(268/(16π))                     | 206.78727      | PDG 206.768       | 9.32e-5  |
| M_Pl/v_EW = exp(4π² − 1)                     | 5.140e16       | obs 4.963e16      | 3.57e-2  |
| σ_SB = π²k_B⁴/(60ℏ³c²)                       | 5.6703744e-8   | CODATA            | < 1e-15  |
| Tetrahedral arccos(−1/3)                     | 109.47122°     | exact             | 0        |
| FCC packing π/(3√2)                          | 0.74048        | exact             | 0        |
| Cabibbo sin θ_C = 1/√20                      | 0.22361        | PDG 0.2257        | 9.27e-3  |
| Deuteron BE = Λ_QCD/90 = 200/90              | 2.2222 MeV     | obs 2.22457 MeV   | 1.06e-3  |

All eight closed forms reproduce as advertised. The largest pure-closed-form
residual is the gravitational hierarchy at ~3.6 % (within its declared
tolerance band) and Cabibbo at ~0.9 % (slightly worse than the often-quoted
"1.3 %" from the master card, but still single-digit-percent). Stefan-
Boltzmann is exact because the formula *is* the textbook definition.

---

## 2. Möbius bundle 11/12 amplitude — cross-module check

`mobius_k4_numerical.run_all` returns three independent estimators of the
ratio that drives the α prefactor:

| estimator         | value            |
|-------------------|------------------|
| analytic          | 0.916666666 (=11/12) |
| combinatorial     | 0.916666666      |
| Riemann quadrature| 0.916666666 (machine eps) |
| Monte-Carlo (5e5) | 0.91574 ± 1.0e-3 |

Combinatorial / quadrature / analytic agree to machine precision; MC sample
agrees within 1 σ. The 11/12 factor *is* the 11 in the α numerator (since
the full prefactor is 11/(48π³) = (11/12)·(1/(4π³)) up to the trivial-bundle
normalisation). The α derivation in `mass_torque_engine` and the geometric
derivation in `mobius_k4_numerical` therefore use the **same** integer.

---

## 3. Drag-as-mass formula — cross-module consistency

The "drag → mass" identity m c² = ℏω_b should be the same across the six
listed modules. Comparing the actual code:

| module                   | ω_b formula                                  |
|--------------------------|----------------------------------------------|
| mass_torque_engine.py    | γ · ξ · ρ · √K · f_int                       |
| drag_mass_generator.py   | (c/ξ) · √(1 + α_drag · γ̃),  γ̃ = γξ/√(Kρ)    |
| lattice_substrate_2d.py  | γ · √(K/ρ) / ξ                               |
| primitive_anchoring.py   | γ = ω_e · ξ / c    (solved from m_e anchor)  |
| generation_map.py        | inherits ω_e baseline; no own ω_b form       |
| bound_state_spectrum.py  | uses confining-potential μ(K, ρ, ξ, γ); no explicit ω_b |

**Finding (mild inconsistency):** the three modules that *do* spell out
ω_b use three different functional forms:

* `mass_torque_engine`:  ω_b ∝ γ · √K           (in substrate-natural units, ρ=ξ=1)
* `drag_mass_generator`: ω_b ∝ √(1 + γ̃) → at small γ̃, ω_b ≈ const + ½γ̃
* `lattice_substrate_2d`:ω_b ∝ γ · √K / ρ^{1/2} · 1/ξ

In the unit-baseline (K=ρ=ξ=γ=1) all three reduce to ω_b = 1 (or O(1)) and
the electron anchor pin makes the difference invisible at the engine's
verification step. Off the unit baseline they would disagree. None of the
configurations ever evaluated in this audit *uses* an off-baseline primitive
set, so the drag-share column reported by `MassTorque.report()` is correct
but doesn't probe this difference. **The three formulas should be reconciled
to a single canonical ω_b; right now they are mutually consistent only at
the anchor point.**

The drag mass identity itself
    m_e (drag) = 0.510998951 MeV
recovered by `PrimitiveAnchoring(anchor='electron_drag')` matches PDG to
~1 ppb (machine round-off), and `cross_check_drag_alpha` returns
α_drag = 0.0072970624 — identical to the value used in
`mass_torque_engine._t_alpha_em` to 12 sig figs.

---

## 4. Integer-set consistency across modules

Every module that names the integer set uses the *same* values:

| integer | mass_torque_engine | generation_map | (others searched) |
|---------|--------------------|----------------|-------------------|
| n_M     | 268                | 268            | (no override seen)|
| N_BAM   | 6                  | —              | —                 |
| K_pair  | 2                  | 2              | —                 |
| K_rank  | 5                  | 5              | —                 |
| n_R     | 18                 | —              | —                 |
| n_A     | 45                 | —              | —                 |
| F       | 2                  | —              | —                 |
| R       | 3                  | —              | —                 |

A `grep` for `N_BAM\s*=`, `n_M\s*=`, `K_pair\s*=`, etc. across all of
`src/stiff_medium/` produced **zero** competing assignments. There is **no
case anywhere of N_BAM=6 in one place and N_BAM=7 in another**, no n_M
collision, no K_pair conflict. The 12-integer set is held in exactly one
place per module and is the same set everywhere.

---

## 5. Contradiction found: two formulas for τ/μ

`mass_torque_engine._t_tau` and `generation_map.lepton_ratio(3,2)` are both
documented as the canonical lepton-tower step from generation 2 to 3, but
they implement **different formulas**:

| source                    | formula                                | value     | rel_err vs 16.817 |
|---------------------------|----------------------------------------|-----------|--------------------|
| `generation_map`          | exp(n_M / (K_rank·(K_rank+1)·π))       | 17.17695  | 2.14 %             |
| `mass_torque_engine._t_tau` | exp(n_M/(K_pair⁴ π) − K_rank/K_pair) | 16.97413  | 0.93 %             |
| observed (PDG)            | m_τ / m_μ                              | 16.81703  | —                  |

Both pass their per-module tolerance, but they cannot both be the canonical
derivation. The `mass_torque_engine` form is a tighter fit and the comment
inside that file (lines 339–345) explicitly notes that the
`(K_rank·n_R − n_M)` form was rejected as negative and a "simpler stable
form" was substituted — i.e. the `_t_tau` formula is acknowledged in-source
to be a heuristic, not a derivation. **This is the cleanest contradiction
in the cross-derivation set: pick one formula and propagate it.**

The same ambiguity does **not** afflict μ/e: both modules use
exp(n_M/(K_pair⁴ π)) = exp(268/16π), which agrees with PDG at 9.3e-5.

---

## 6. Cabibbo consistency: cabibbo_substrate ↔ mixing_matrices

`MixingMatrices().lam` (Wolfenstein λ) returns 0.22360679… = 1/√20 exactly,
matching the closed form quoted in the audit prompt. PDG sin θ_C is
0.2257, residual 0.93 %. Both modules report the same B3 prediction; no
cross-module conflict.

The full PMNS comparison (executed live):

| observable          | predicted   | observed   | residual |
|---------------------|-------------|------------|----------|
| sin θ_C             | 0.22361     | 0.2253     | -0.75 %  |
| |V_us|              | 0.22361     | 0.2243     | -0.31 %  |
| |V_cb|              | 0.04055     | 0.041      | -1.10 %  |
| |V_ub|              | 0.003520    | 0.00382    | -7.86 %  |
| J_CKM               | 2.842e-5    | 3.18e-5    | -10.62 % |
| sin² θ_12 (PMNS)    | 0.30649     | 0.307      | -0.17 %  |
| sin² θ_13 (PMNS)    | 0.02189     | 0.022      | -0.49 %  |
| sin² θ_23 (PMNS)    | 0.54585     | 0.572      | -4.57 %  |
| δ_CP (PMNS)         | 2.356 (3π/4)| 3.84       | -38.6 %  |
| J_PMNS              | 0.02349     | 0.0329     | -28.6 %  |

The PMNS angles themselves (12, 13) are excellent; θ_23 is mediocre and
the CP phase is the well-known weak point of the substrate ansatz.

---

## 7. Σm_ν ↔ PMNS coupling

The Σm_ν = 60.5 meV prediction (B3 cosmology, hubble derivation) is **not**
encoded in `mixing_matrices.py`; that module only computes the angles and
Jarlskog invariants and never references the absolute neutrino mass scale.
Hence there is no consistency check to perform in the present codebase
between Σm_ν and PMNS — they are computed in disjoint modules and would
*only* couple via mass-squared splittings, which `mixing_matrices` does
not output. **No contradiction, but also no positive cross-check.** This
is a gap, not a failure.

---

## 8. Summary

### Pass cleanly
* All 9 `MassTorque` registry items reproduce observed values within
  declared tolerance (α at 4e-5, μ/e at 9e-5, deuteron at 1e-3, hierarchy
  at 3.6 %, etc.).
* α closed form 11/(48π³)·exp(−3π/737) and the geometric Möbius bundle
  ratio 11/12 agree to machine precision via three independent estimators.
* All twelve B3 integers (n_M=268, N_BAM=6, K_pair=2, K_rank=5, n_R=18,
  n_A=45, F=2, R=3, plus the four primitives K, ρ, ξ, γ defaulting to
  unity) are assigned in exactly one place per module and never disagree
  numerically across the codebase.
* Drag-mass anchor: `PrimitiveAnchoring(anchor='electron_drag')` solves
  γ = 7.76e20 Hz from m_e c² = ℏγ, recovers m_e at 1 ppb, and yields the
  same α value as `mass_torque_engine._t_alpha_em` to 12 sig figs.
* Cabibbo λ = 1/√20 returned by `MixingMatrices` matches the closed form
  quoted in the audit prompt; PMNS θ_12/θ_13 match observation at sub-1 %.
* Stefan-Boltzmann reproduces CODATA at 1e-15 (it is the textbook formula).
* Tetrahedral angle and FCC packing are mathematical identities and pass
  trivially.

### Contradictions / inconsistencies
1. **τ/μ ratio is computed two different ways.**
   `generation_map`: exp(n_M/(K_rank(K_rank+1)π)) = 17.18 (2.1 % off).
   `mass_torque_engine`: exp(n_M/(K_pair⁴ π) − K_rank/K_pair) = 16.97 (0.9 % off).
   The in-source comment of `mass_torque_engine` admits its form is a
   heuristic. These need to be reconciled.

2. **Three different functional forms of ω_b** appear in
   `mass_torque_engine`, `drag_mass_generator`, and `lattice_substrate_2d`.
   At the unit-baseline electron anchor they all collapse to O(1) and the
   verifier doesn't see the difference, but off-baseline they would disagree.
   They should be unified to a single canonical drag-bouncing law.

### Numerical residuals worth noting
* Hierarchy M_Pl/v_EW: 3.6 % (just inside the 5 % tolerance, but the
  largest residual among the headline derivations).
* Cabibbo sin θ_C: 0.93 % (the closed form gives 0.22361 vs PDG 0.2257).
* τ/m_μ heuristic: 0.93 % (mass_torque_engine) or 2.1 % (generation_map).
* Higgs mass: 1.8 % (uses an explicit 5/4 calibration coefficient inside
  the formula, so this is partially fitted, not pure derivation).
* δ_CP (PMNS) and J_PMNS: 28–39 % off — known weak point.

### Gaps (not failures)
* Σm_ν = 60.5 meV (cosmology) and the PMNS angles in `mixing_matrices.py`
  live in disjoint modules and are not cross-checked anywhere in the
  codebase. A future module computing mass-squared splittings from Σm_ν
  + ordering would close this loop.
* The drag-as-mass machinery is exercised only at the electron anchor; a
  multi-particle off-baseline test would discriminate the three ω_b forms.

### Bottom line
The framework's headline closed forms reproduce as advertised, the
12-integer set is held consistently across all modules, and the α/Möbius
cross-check is rigorous. The two real consistency issues are (a) the
duplicate τ/μ formula and (b) the three competing ω_b laws. Neither
invalidates a published result, but both are genuine technical-debt items
that should be resolved before the next round of "end-to-end" claims.
