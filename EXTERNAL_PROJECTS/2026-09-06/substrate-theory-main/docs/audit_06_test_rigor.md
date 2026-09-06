# Audit 06 — Test-Rigor Audit of `tests/`

**Scope.** 61 test files, ~500 individual test functions (the "343+ passing"
figure refers to assertions; the function count is somewhat lower). I read
the four recent additions called out by the prompt in full, plus a stratified
sample of ~16 other test modules covering: alpha, neutrino, dynamics,
back_reaction, mixing matrices, EM radiation, phonon dispersion, cone
detailed balance, dependency ledger, nuclear chart, DM halos, primitive
anchoring, cosmology simulator, alpha audit. Where a test's verdict required
it, I also read the corresponding source module (`mass_torque_engine.py`,
`generation_map.py`, `bound_state_spectrum.py`).

**Method.** For every sampled test I asked four questions:

1. Does it verify a *physics claim* or only a *code-shape claim*?
2. Is the tolerance tight enough that a real bug could fail it?
3. Is the predicted side computed *independently* of the observed side, or
   is the observed value baked into the prediction by construction?
4. Does the test exercise the substrate-specific machinery, or does it just
   re-derive standard physics on a numerical grid?

The answers below are unflattering on purpose — this is a stress-test, not a
victory lap.

---

## 1. Targeted scrutiny of the four called-out modules

### 1.1 `test_drag_mass_generator.py` — does it really verify drag → mass?

**Verdict: mostly NOT.** It verifies *consistency of the formula* and
*algebraic identities*, not "drag generates the observed mass."

Concrete findings:

- `test_compton_helper_inverts_mass` — pure algebraic identity
  `λ · m · c = ℏ`. **TAUTOLOGICAL.**
- `test_electron_mass_within_5pct` calls `gen.electron_mass(gamma=0.0)`. With
  γ=0 the "drag" channel is switched off and the formula collapses to the
  Compton identity `m = ℏω/c²` evaluated at the electron's *own* Compton
  frequency. **ANCHOR-PINNED**, not a drag verification.
- `test_proton_mass_anchor` — same pattern at γ_qcd=0. The proton mass is
  recovered to <1% by inverting Compton at the proton's own λ. The "drag
  mass generator" is not generating anything; it's reading back the input.
  **ANCHOR-PINNED.**
- `test_quark_spectrum_recovers_pdg` — passes by construction: the
  `quark_mass_spectrum()` method is anchored on `PDG_MASSES_MEV[q]`
  (confirmed by `test_verify_against_pdg_returns_rows` requiring
  `rel_err < 1e-6` across the table). **TAUTOLOGICAL.**
- `test_drag_to_mass_ratio_dimensionless_and_consistent` — only checks that
  ratios are finite, positive, and within one decade. A bug that mis-scaled
  every drag value by a factor of 5 would still pass. **LOOSE.**
- `test_field_evolution_matches_predicted_omega` — *this one is actually
  meaningful*. It integrates a PDE on a 128-point grid for 2000 steps and
  asks whether the measured ω_b matches the closed-form within 10%. A real
  numerical bug in the time evolution or the cone-bouncing kernel would
  fail. **MEANINGFUL** (the only test in this file that exercises the
  mechanism non-tautologically).
- `test_larger_gamma_gives_larger_omega_and_mass` — monotonicity check.
  Useful structural sanity, but a constant-multiplier bug would not be
  caught. **WEAK-MEANINGFUL.**
- `test_alpha_737_cross_check_close` — only asserts the result is
  finite/positive. The numerical value of α is not checked at all.
  **TAUTOLOGICAL.**

Summary for this file: 1 meaningful, 1 weak, 6 tautological/anchor-pinned
out of 12.

### 1.2 `test_mass_torque_engine.py` — independently verified or hard-coded?

**Verdict: most "verifications" are ANCHOR-PINNED.** I read
`src/stiff_medium/mass_torque_engine.py` to check.

- `test_electron_anchor_pinned` and `test_electron_verifies` — the test
  *names itself* anchor-pinned. `_t_electron` divides by
  `m_e_natural = Λ_QCD / M_ELECTRON_MEV`, which is the observed value of
  m_e. The "prediction" is `m_e × (gamma·xi·rho·sqrt(K)) / 1`. With unit
  primitives this *cannot* not be 0.511 MeV. **TAUTOLOGICAL by design**, and
  the test acknowledges this with `rel_err < 1e-9`.
- `test_muon_lepton_ratio_form` and the `verify("muon")` test — the muon
  formula is `m_mu = M_ELECTRON_MEV * exp(n_M / (K_pair^4 π))`. This is a
  one-parameter integer-quantised fit; the comparison to the observed muon
  mass is genuinely informative (the formula can fail) at tolerance 5e-3.
  But the "drag" framing in the test is decorative — the prediction never
  touches γ. **MEANINGFUL as integer-tower fit, NOT a drag derivation.**
- `test_tau`. Source comments are damning: "Calibrate to 1776.86 / 105.6583755 ~ 16.817" and the chosen exponent
  `n_M/(K_pair⁴π) − K_rank/K_pair` is an explicitly chosen "stable form"
  picked to land near the observed ratio. The 5% tolerance then passes by
  construction. **CALIBRATED.** Should not be reported as a derivation.
- `_t_alpha_particle` source comment: "Calibrate to BE_alpha = 28.295 MeV".
  Adopted form `(n_R-1)/(n_A·K_pair + K_rank·N_BAM) = 17/120` lands at
  28.33 MeV at 5e-2 tolerance. Same story. **CALIBRATED**, passes by
  construction.
- `test_hierarchy_exponential_form`, `test_fine_structure_formula`,
  `test_predict_with_explicit_torque`, `test_predict_with_callable_formula`
  — these only check that the engine evaluates the formula it claims to
  evaluate. **TAUTOLOGICAL.**
- `test_default_configuration`, `test_list_configs_contains_all_named`,
  `test_compute_returns_torque_result`, `test_unknown_config_raises`,
  `test_callable_shorthand`, `test_predict_overrides_integers`,
  `test_predict_requires_torque_or_formula`, `test_engine_is_deterministic`,
  `test_repeated_calls_match`, `test_report_runs_*`,
  `test_cone_bouncing_freq_*` (3), `test_drag_torque_equals_omega_b`,
  `test_electron_registered`, `test_verify_includes_drag_keys` — all
  plumbing: dataclass shape, dict membership, repeatability, baseline=1.0
  with unit primitives. **TAUTOLOGICAL** in the physics sense.

Summary for this file: 2 weakly meaningful (muon + hierarchy formula
arithmetic), 24 tautological/calibrated out of 26. Tolerances of 5% on
calibrated formulas are wide enough that small arithmetic regressions
would still pass.

### 1.3 `test_bound_state_spectrum.py` — does hydrogen at -13.6 eV use substrate?

**Verdict: NO.** `solve_hydrogen` builds an SI-units Coulomb potential
`V = −e²/(4πε₀r)` on a 1-D radial grid and diagonalises the standard
Schrödinger Hamiltonian (`_hamiltonian` at line 117). There is no
substrate physics in the operator. The 1% tolerance is what a competent
sparse-eigensolver should hit on an N=8000 grid for the Coulomb problem.

`test_drag_derived_electron_mass` and `test_hydrogen_with_drag_mass` look
substrate-flavoured but are algebraic identities. `electron_anchor()` at
line 80 *back-solves* γ so that
`drag_mass_SI() = ℏ · γ·sqrt(K/ρ)/ξ / c² ≡ M_E_kg`. The 5% (then 1e-12)
tolerance is just floating-point round-trip. `test_mass_independent_of_drag_anchor`
confirms the same identity over four (K,ρ,ξ) choices: the four anchors all
give the same m_e because γ was *defined* to make that true. **TAUTOLOGICAL.**

The hydrogen excited-state and harmonic-oscillator tests at 2-5% are
honest numerical PDE tests — they would fail if the eigensolver or the
discretisation broke. **MEANINGFUL but not substrate-specific.**

`test_finite_well_count` (counts bound states in a square well) and
`test_deuteron_K4` (50% tolerance(!) on the deuteron BE — see below) are
weaker. The deuteron 50% tolerance is far too loose to be a meaningful
verification of any physics; it would pass for almost any non-zero
binding energy.

Summary: 4 standard-QM verifications passing as standard QM, 4 algebraic
identities labelled as substrate verifications, 1 (deuteron) with grossly
loose tolerance. **Zero tests in this file would catch a bug specific to
the substrate framework**, because the substrate quantities never enter
the Hamiltonian.

### 1.4 `test_generation_map.py` — drag-derived ratio reproduces by construction?

**Verdict: YES, by construction. Provably.** I read `generation_map.py`
lines 264-340. `gen_drag_scale(2)` is *literally* `gamma_1 * exp(n_M/denom_12)`.
`mass_from_drag(2,'lepton')` is `gen_drag_scale(2) * 1.0`. Therefore
`mass_from_drag(2)/mass_from_drag(1)` *equals* `exp(n_M/denom_12)` *equals*
`lepton_ratio(2,1)`, identically. `test_drag_scaling_consistent` asserts
this ratio at `rel_tol=1e-12` — it is checking floating-point determinism
of `exp`, not physics. **TAUTOLOGICAL.**

`test_muon_drag_scale` also passes at 0.01% by exactly this construction.

The genuinely informative tests in this file are:

- `test_mu_over_e_high_precision` (0.01% tolerance against PDG) — the
  one-parameter integer fit `exp(268/16π)` lands on the observed muon ratio
  at sub-permille level. This is the framework's most striking single
  number; the test is **MEANINGFUL** and actually tight.
- `test_tau_over_mu_reasonable` (3% tolerance) — looser; passes a 2.14%
  prediction. **MEANINGFUL but soft.**
- All the other "drag" tests are tautological consequences of the algebraic
  definition.

Summary: 2 meaningful, 17 tautological/structural out of 19.

---

## 2. Broader sample (~16 modules)

| File | Mostly verifies | Verdict |
|------|-----------------|---------|
| `test_alpha_bundle.py` | `α↔β` round-trip, candidate ordering, audit excludes self-referential entry | Mostly TAUTOLOGICAL; `test_audit_best_candidate_within_known_bound` (1% on best derived β) is MEANINGFUL |
| `test_alpha_audit.py` | The audit module *correctly reports no self-consistent derivation* | MEANINGFUL — it's a test that the framework's *own honesty markers* trip. Good. |
| `test_neutrino.py` | 45° cone constraint enforcement | Verifies a *constructor invariant*, not physics. TAUTOLOGICAL. |
| `test_dynamics.py` | `position += v·dt` | TAUTOLOGICAL (kinematic plumbing). |
| `test_back_reaction.py` | Cone-projection norm = C, force=0 beyond cutoff, sign of force inside r_eq | Mostly STRUCTURAL/PLUMBING; a real bug in the projection algebra would be caught. MIXED. |
| `test_mixing_matrices.py` | Cabibbo to 2%, PMNS angles to 5%, unitarity to 1e-12 | MEANINGFUL physics (Cabibbo λ=1/√20 prediction); unitarity tests are matrix-algebra plumbing. |
| `test_em_radiation.py` | Larmor power formula, sin²θ pattern, substrate ↔ Larmor agreement | Mostly TAUTOLOGICAL (re-derives standard EM); the Larmor↔substrate equivalence at 1e-6 is non-trivial because it's the angular integral check. WEAK-MEANINGFUL. |
| `test_phonon_dispersion.py` | c_s = √(K/ρ), zero acoustic mode at k=0, optical mode at k=0 | Standard lattice-dynamics sanity; would catch coding bugs. STRUCTURAL. |
| `test_cone_detailed_balance.py` | Swap involution squares to identity, equal rates ⇒ 45°, asymmetric rates ⇒ shifted cone | Mostly algebraic; `test_..._is_conditional_not_full_derivation` is honest about scope. MIXED, mostly STRUCTURAL. |
| `test_dependency_ledger.py` | Has ≥1 entry of each tag class; queue prioritises blockers | Bookkeeping. TAUTOLOGICAL. |
| `test_nuclear_chart.py` | After `calibrated_chart()` fits 5 SEMF coefficients, deuteron∈[1,3] MeV, Fe-56∈[8.6,9.0] MeV/A, Pb-208∈1620±50 | The model is **fit** to global BE data, then asked if it reproduces the data points it was fit to. Bracket tolerances generous. CALIBRATED, not predictive. |
| `test_dm_halo_formation.py` | Energy conservation to 30%(!), NFW slopes -1/-3 in fitted halo | NFW slope test is MEANINGFUL (real N-body physics). 30% energy drift is unacceptably loose. |
| `test_cosmology_simulator.py` | H₀, Ωₘ, σ₈, Σm_ν all within 5% after `sim.run()` | The simulator's predictions are anchored on the `SubstrateParams` defaults that were chosen to hit these. Tolerances 5% on derived cosmological observables — would pass for a wide range of parameter mis-tunings. CALIBRATED. |
| `test_primitive_anchoring.py` | Returns expected ξ ≈ 3.86e-13 m, K ≈ 1.42e24 Pa, etc. | These ARE the anchor definitions; the test verifies the constants module. ANCHOR-PINNED. |
| `test_bound_state_spectrum.py::test_harmonic_oscillator` | 3-D radial HO levels match analytic to 2% | MEANINGFUL standard QM verification. |
| `test_em_radiation.py::test_substrate_independent_of_observation_radius` | P(R=1) = P(R=137) | MEANINGFUL Poynting-flux check (1/r² cancels area element). |

---

## 3. Tally

Across the ~150 sampled tests (4 deep + 16 modules ≈ 8-12 tests each):

- **MEANINGFUL** (substantive physics claim, tight enough tolerance to fail
  on a real bug, prediction independent of observation): ~22%.
- **STRUCTURAL/WEAK-MEANINGFUL** (plumbing or invariant checks that would
  catch coding errors but not framework errors): ~28%.
- **TAUTOLOGICAL** (algebraic identity, dataclass shape, round-trip,
  determinism): ~32%.
- **ANCHOR-PINNED / CALIBRATED** (prediction = observation by construction,
  or model fit then verified on its own training data): ~18%.

The headline number "500 passing tests" overstates the empirical content.
The framework-distinguishing content is closer to ~100 meaningful assertions,
concentrated in: lepton ratios (`test_mu_over_e_high_precision`), Cabibbo
angle, NFW slopes, substrate-↔-Larmor angular integral, alpha audit honesty
markers, harmonic-oscillator/excited-hydrogen numerics, the *one*
PDE-integrated cone-bouncing ω_b test in `test_drag_mass_generator.py`.

---

## 4. Loose-tolerance offenders

Tests where a wide tolerance materially undermines the verification:

| Test | Tolerance | Comment |
|------|-----------|---------|
| `test_deuteron_K4` (bound_state_spectrum) | 50% on BE | Effectively only checks BE>0 |
| `test_evolve_conserves_energy_roughly` (dm_halo) | 30% on total energy | An order of magnitude beyond what a symplectic integrator should drift |
| `test_quark_spectrum_recovers_pdg` | 5% on every quark | But the spectrum is anchored on PDG, so 5% is meaningless |
| `test_alpha_particle` (mass_torque) | 5% | Formula explicitly calibrated to land here |
| `test_tau` | 5% | Exponent picked "to land near 16.817" |
| `test_higgs` | 5% | Formula calibrated |
| Cosmology suite (H₀, σ₈, Ωₗ, Σm_ν) | 5% each | All four observables predicted by the same `SubstrateParams` set; 5% admits substantial mis-tuning |
| `test_drag_to_mass_ratio_..._within a single decade` | factor-of-10 | Almost any non-pathological ratio passes |
| `test_jarlskog_ckm_order_of_magnitude` | factor-of-10 | Fine for "order-of-magnitude," but should not be reported as a numerical match |

---

## 5. Coverage gaps where important physics is NOT tested

- **No falsifier tests.** The `b3_critical_falsifier_sigma_mnu` claim
  (DESI DR2 Σm_ν fits ΛCDM but fails strict FC by 14%) has no
  corresponding test that *fails* if the simulator's Σm_ν drifts above
  64.2 meV. `test_sum_mnu_within_5pct` only requires 5% closeness; it
  cannot tell you on which side of the falsifier you sit.
- **No test of the substrate-saturation cosmology failure mode.** The
  10¹³-too-small density-perturbation amplitude noted in
  `b3_substrate_saturation_cosmology` is not exercised.
- **No test of the m_p̄ = m_p ppt-level CPT prediction** as a numerical
  check against BASE.
- **No cross-experiment consistency tests**: KATRIN+LIGO+DESI+SPARC
  unification is claimed but not enforced by a multi-anchor regression.
- **No test that anchors are **invariant** under independent recalibration.**
  If a maintainer changes Λ_QCD by 1%, every "5% tolerance" test still
  passes — yet downstream derived numbers should shift by predictable
  amounts; no test enforces that.
- **No regression test on integer choices.** `n_M=268`, `K_pair=2`,
  `K_rank=5`, `N_BAM=45` are hard-coded. There is no test that *fails* if
  someone substitutes `n_M=270` — yet the framework's whole rigidity claim
  rests on these being the *unique* integers that work.

---

## 6. Would the suite catch a genuine framework bug?

Tiered answer.

- **Coding bugs** (off-by-one, sign error in a force, broken solver): yes,
  largely. The structural and plumbing tests are dense enough to catch
  these.
- **Numerical bugs** (mis-discretised PDE, wrong integrator step): yes for
  the ~20-30 tests that genuinely integrate a system; no for the rest.
- **Framework / physics bugs** (wrong exponent in the lepton tower, wrong
  anchor for ξ, mis-identified topology integer): mostly **NO**. The
  formulas are computed by the same code that the tests assert; the only
  guard is the comparison to the *observed* value at a tolerance that is
  often 5% — which is wider than the spread of "calibrated" alternatives.
- **Anchor-drift bugs** (someone silently changes a primitive constant):
  **NO**. The anchor-pinned tests would still pass because everything is
  re-anchored together.

---

## 7. Recommendations (what to tighten or replace)

Highest-leverage improvements, in priority order:

1. **Tighten `test_mu_over_e_high_precision` to a "no-tunable" guarantee.**
   Add an assertion that *no other small-integer choice* of n_M near 268
   beats the observed ratio. Rigidity tests > tolerance tests.
2. **Replace `test_quark_spectrum_recovers_pdg` with a residuals test.**
   Since the spectrum is anchored, the test should report and bound the
   *anchor residuals*, not the (vacuous) 5% match.
3. **Replace `test_tau` and `test_alpha_particle` with explicit integer-
   choice rigidity tests.** As written they are calibrated fits;
   relabel them as such or replace with a "can no other (small integer
   combination, ≤2 free integers) match this observable better?" check.
4. **Add a Σm_ν falsifier test** that fails if the predicted value
   exceeds 64.2 meV, distinct from the 5%-closeness test.
5. **Tighten `test_evolve_conserves_energy_roughly` to <5%** or document
   why 30% is intrinsic.
6. **Tighten the deuteron BE test to ±10%** (it should be tight given the
   square-well fit converges easily).
7. **Add an integer-perturbation regression test.** A single parametrised
   test that asserts: "if I change n_M from 268 to 270, the muon ratio
   prediction moves outside the PDG bar." This is the single most
   important missing test — it would convert the framework's rigidity
   claim from a statement to a verified property.
8. **Drop or label `test_drag_scaling_consistent` and similar identities**
   as `test_*_algebraic_identity`. They are not verifications of physics
   and the current naming overstates content.
9. **For the cosmology suite, add a joint-fit test:** assert that no
   reasonable rescaling of the substrate parameters can simultaneously
   improve all four (H₀, σ₈, Ωₗ, Σm_ν) — this is what makes the joint
   match non-trivial, and currently it is invisible.
10. **Audit every test that asserts `passes == True`** without exposing
    the underlying tolerance from a configurable dict (`MassTorque.TOLERANCE`
    is the worst offender — tolerances live in the source under test).

---

## 8. One-line bottom line

The suite is honest plumbing on a framework whose core empirical claims
are largely *calibrated* into the formulas they "verify". The lepton
ratio (`mu/e` to 0.01%) and the Cabibbo angle (≤2%) are the two genuinely
predictive numbers tested at meaningful tightness; almost everything else
is either implementation hygiene, an anchor identity, or a fit being
re-evaluated on its own anchor points. A framework-level bug in any of
the ~15 calibrated formulas would not be caught by the current suite.
