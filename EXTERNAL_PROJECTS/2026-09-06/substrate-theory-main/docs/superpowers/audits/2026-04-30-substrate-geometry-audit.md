# Substrate Geometry Honesty Audit
**Date:** 2026-04-30
**Author:** code audit pass
**Question:** For each "predictive" module, how much is the substrate actually
producing geometry vs. how much is QFT/QM with substrate-labeled inputs?

## Categorization key

- **A — REAL substrate geometry**: solves substrate field equations / lattice
  problems / kink dynamics directly. Output is a geometric structure (kink
  profile, bound state wavefunction, lattice eigenvalue, RG flow) emerging
  from substrate dynamics.
- **B — SUBSTRATE-DERIVED FORMULA**: the formula itself comes from substrate
  mechanics (Y-junction geometry, half-flux holonomy, σ × ξ scaling, kink
  spin-spin contact). Substrate is producing the formula, not just feeding
  inputs.
- **C — INHERITED FORMULA + substrate input**: takes a known QFT/QM formula
  (Rydberg, Goldberger-Treiman, Cornell potential, Sirlin V-A, etc.) and
  plugs in substrate-derived numbers for the inputs. The structure of the
  formula is imported.

---

## Module-by-module audit

### 1. `src/stiff_medium/substrate_rg_running.py` — K(ξ) running
- **Category:** **A**
- **What it actually computes:** Solves a one-parameter ODE
  `dK/d(ln ξ) = β_K(K, ξ)` for several β-function families (power-law,
  log-running, QCD-like) anchored at three physical scales (electron Compton,
  QCD scale, Planck) and asks whether a single 2-parameter β can fit all
  three anchors.
- **Substrate content:** The whole framing — "K is a substrate stiffness
  field that runs with scale" — is substrate mechanical. The anchors are
  derived from substrate identities (`K_e = ℏc/ξ_e^4`, `K_Planck = c^7/(ℏG^2)`,
  `K_QCD = σ_QCD/σ_lattice`).
- **Genuine substrate prediction:** The fitted exponent `a ≈ 5.69` for the
  power-law `K(ξ) = K_e (ξ_e/ξ)^a` interpolating electron→QCD scales. This
  is a real RG-flow output, not a plug-in.

### 2. `src/stiff_medium/regge_spectrum.py` — hadron Regge masses
- **Category:** **mixed B/C**
- **What it actually computes:** Builds Regge trajectories `M² = (J − α₀)/α'`
  with slope `α' = 1/(2πσ)` (mesons) or `α'/2` (closed strings) using
  substrate-derived σ from K(ξ). Anchors α₀ to one hadron per channel, then
  predicts the rest.
- **Substrate content:** σ comes from the running K, which is substrate-A.
  But the formula `α' = 1/(2π σ)` is the Nambu-Goto string slope —
  inherited from string theory, not derived from substrate mechanics within
  this module. Channel intercepts α₀ are ONE-PARAMETER FIT per channel.
- **Genuine substrate prediction:** σ is substrate-derived, and the
  Y-junction baryon slope (α'_open same as meson) is justified empirically
  at this level. The trajectory shape itself is inherited string theory.

### 3. `src/stiff_medium/chromomagnetic_substrate.py` — N-Δ splitting
- **Category:** **B**
- **What it actually computes:** `Δm_NΔ = 4 ξ² σ^{3/2}` — derived from
  composing (8π/3) × σξ² (Möbius coupling α_M) × C_F=4/3 × σ^{3/2} contact
  density, divided by `m_K_eff² = σ`.
- **Substrate content:** The CLOSED-FORM `Δm_NΔ = 4 ξ² σ^{3/2}` is a
  substrate-mechanical formula: every factor (α_M = σξ², ψ(0)² ~ σ^{3/2},
  m_K_eff² = σ) is substrate-derived. The only "inherited" piece is the
  spin-spin contact STRUCTURE (Breit-Pauli) and the C_F = 4/3 colour Casimir,
  which §18.49 declares to be inherited from SU(3). Output: 313.8 MeV vs
  293 MeV empirical, +6.8% with no fit.
- **Genuine substrate prediction:** Yes — the closed form `4 ξ² σ^{3/2}` is
  a real substrate-mechanical answer to a chromomagnetic question.

### 4. `src/stiff_medium/glueball_closed_string.py` — α'_closed = 3/4
- **Category:** **B** (with a hand-picked branch)
- **What it actually computes:** Six topology candidates for closed-loop
  Regge slope (NG, Möbius antiperiodic, breathing zero-point, half-flux J+1,
  hybrid arithmetic mean, hybrid geometric mean). Each is derived from a
  mode-counting argument. The "winner" — substrate hybrid arithmetic — gives
  `α'_closed = (3/4) α'_open` matching lattice QCD.
- **Substrate content:** Each candidate's derivation is substrate-mechanical
  (Möbius half-flux modes, 3D breathing, hybrid open+closed kink loop).
  However: the 3/4 ratio comes specifically from POSTULATING that the
  closed kink loop is "half-closed + half-open" with arithmetic mean of the
  two slopes. That hybrid argument is substrate-flavoured but is the
  selection criterion for matching lattice; the module honestly tabulates
  all six.
- **Genuine substrate prediction:** Mode-counting on each topology is real
  substrate work. The choice of arithmetic-mean hybrid as physically correct
  rather than geometric-mean (1/√2) is post-hoc.

### 5. `src/stiff_medium/proton_radius_v3.py` — proton charge radius
- **Category:** **A** (for the spread term) + **C** (for the rest)
- **What it actually computes:** Solves the radial Schrödinger equation
  `-u''/(2μ) + σr u = Eu` numerically (sparse eigsh) AND variationally
  (Gaussian) for one effective constituent in the linear potential V = σr,
  to get `⟨r²⟩_3body`. Adds Lüscher transverse string fluctuations and
  intrinsic kink width to assemble r_p².
- **Substrate content:** The Schrödinger problem in V = σr with μ = √σ/2 is
  a genuinely solved substrate problem (the Airy ground state in the
  substrate confining potential). The Lüscher and intrinsic-kink-width
  contributions are inherited from string theory / sech² profile.
- **Genuine substrate prediction:** ⟨r²⟩ from the Schrödinger eigenproblem
  in the substrate confining potential — that's category A. The composition
  with Lüscher logs is C.

### 6. `src/stiff_medium/hyperon_spectrum.py` — hyperons
- **Category:** **B/C mixed**, leans **C**
- **What it actually computes:** Sakharov-De Rújula-Georgi-Glashow octet/
  decuplet mass formula `m = Σ m_q + ΔE_chromo` with chromomag shift from
  eq. (2). Anchors m_q_struct from proton (eq. 7) and m_s_struct from Λ (one
  flavor parameter). Predicts Σ, Ξ, Ω.
- **Substrate content:** The chromomag pair coupling K_substrate = (8/3)
  σ^{3/2} ξ² is substrate-B. The Sakharov-DRGG SU(6) spin-flavor coefficients
  (-3/4 N, +1/4 - 1 Σ, etc.) are imported from the standard quark model and
  declared to be inherited via §18.49 SU(3) gauge inheritance — this is C.
- **Genuine substrate prediction:** Numerical chromomag couplings inserted
  into a fixed SU(6) Clebsch-Gordan structure. The SU(6) algebra is not
  derived inside this module.

### 7. `src/stiff_medium/nucleon_magnetic_moments.py` — μ_p, μ_n
- **Category:** **C**
- **What it actually computes:** `μ_p = (4/3) μ_u − (1/3) μ_d` and
  `μ_n = (4/3) μ_d − (1/3) μ_u` — Wigner SU(2) spin-flavour formulas — with
  Dirac magnetic moments `μ_q = e_q ℏ/(2 m_K c)` and substrate-derived
  m_K_eff candidates plugged in.
- **Substrate content:** The choice m_K_eff = √σ (chromomagnetic substrate)
  is substrate-derived. Quark charges (+2/3, −1/3) come from §18.49 SU(3)
  inheritance.
- **Genuine substrate prediction:** Just the input m_K_eff. The Wigner
  algebra (4/3, −1/3) and the Dirac μ = eℏ/(2mc) are entirely inherited.

### 8. `src/stiff_medium/pion_decay_constant.py` — f_π = ½σξ
- **Category:** **B**
- **What it actually computes:** Tables ~20 dimensional combinations of σ and
  ξ at the QCD scale (√σ/(2π), √(σ/(4π²)), σξ², etc.) and reports the
  fractional error of each vs f_π = 92.4 MeV. The cleanest match is
  `f_π = ½ σξ_QCD = 91.22 MeV` (−1.3% error, no fit).
- **Substrate content:** ½σξ is a substrate-mechanical quantity — energy per
  Möbius half-flux quantum stored in one coherence-domain of confining flux.
  This is purely σ × ξ × geometric ½, no QFT formula imported.
- **Genuine substrate prediction:** Yes — `f_π = ½σξ` is a clean
  substrate-derived formula. The CHOICE of ½σξ from a list of dimensional
  candidates is a "best-of-many" selection but each candidate is a
  substrate-natural quantity.

### 9. `src/stiff_medium/strange_quark_sector.py` — f_K
- **Category:** **B**
- **What it actually computes:** `f_K = ½ σ ξ × cosh(½ ln(m_K/m_π))` with
  the cosh factor representing the substrate's "geometric-mean tension" of a
  string with two unequal-stiffness endpoints. m_K, m_π are empirical inputs.
- **Substrate content:** The cosh form is derived from the symmetric
  substrate argument that a string with two different-stiffness endpoints
  has σ_eff/σ = ½(√(m_K/m_π) + √(m_π/m_K)) = cosh(½ ln(m_K/m_π)). This is
  substrate-mechanical. The light/strange mass ratio enters as an empirical
  numerical input.
- **Genuine substrate prediction:** The cosh formula itself, predicting
  f_K/f_π = 1.206 vs 1.193 empirical (+1.3%) without a fit parameter.

### 10. `src/stiff_medium/cabibbo_angle.py` — Cabibbo
- **Category:** **C** (best matches) / honest "no derivation" (rest)
- **What it actually computes:** Catalogues many candidate substrate-topological
  formulas for sin(θ_C). Best match: GIM identity sin(θ_C) = √(m_d/m_s).
  Module's honest verdict: the angle is NOT yet derived from substrate
  topology. The §18.49.7 spec admits this is empirical input.
- **Substrate content:** The TESTED candidate forms (π/N, Z₃ × half-flux
  combos, foot-phase residuals) are substrate-flavoured but none cleanly
  reproduce 13.04°. The GIM identity is the best numerical match but it
  reduces θ_C to the d/s mass ratio, which itself is an empirical input.
- **Genuine substrate prediction:** None — this is honestly flagged by the
  module as an open problem.

### 11. `src/stiff_medium/isospin_breaking.py` — m_d − m_u
- **Category:** **mixed B/C**
- **What it actually computes:** `(m_d − m_u) = α × (a₁/3) √σ ≈ 2.41 MeV`
  where (a₁/3) √σ ≈ 330 MeV is the Airy-kinetic constituent mass (from the
  σr Schrödinger problem). Empirical 2.5 MeV; +4% error without fit.
- **Substrate content:** The Airy mass from the σr Schrödinger problem is
  category-A substrate work (re-used here). The α × m factor comes from a
  Möbius orientation cost argument (§18.30) — substrate-mechanical, but
  qualitative. EM Coulomb self-energy (3/5)αℏc Q²/r is a known QED formula
  imported.
- **Genuine substrate prediction:** The Airy mass scale IS substrate
  geometry; combining it with α × Q² is a hybrid substrate + QED argument.

### 12. `src/stiff_medium/atomic_transitions.py` — hydrogen spectrum
- **Category:** **C** (most explicit case in the framework)
- **What it actually computes:** Rydberg `R_y = ½ m_e c² α²`, Lyman-α =
  3R_y/4, Balmer-α = 5R_y/36, 21-cm hyperfine `(8/3) g_p α² (m_e/m_p) R_y`,
  fine structure `α² R_y / 16`, Welton-form Lamb shift `(α/π) α⁴ R_y ln(1/α²)`.
- **Substrate content:** Just the inputs. α is substrate-derived (photon
  renormalization on Möbius bundle). m_e is substrate fundamental. m_p is
  substrate Y-junction. g_p is substrate SU(2) spin-flavour. The FORMULAS
  (Bohr 1/n², hyperfine Fermi contact, fine structure α²R_y/16, Welton heuristic
  for Lamb) are entirely classic QM/QED.
- **Genuine substrate prediction:** None within this module — it is a
  recombination of substrate-derived inputs into inherited Bethe-Salpeter
  formulas. The module's docstring is honest about this: "the hydrogen atom
  is therefore a *recombination* of substrate-derived quantities".

### 13. `src/stiff_medium/nuclear_binding.py` — light nuclei
- **Category:** **C**
- **What it actually computes:** Goldberger-Treiman `g_πNN = m_N g_A / f_π`
  with substrate f_π. One-pion-exchange potential (textbook Ericson-Weise
  form). Deuteron variational Gaussian. Bethe-Weizsäcker SEMF with
  dimensionally-derived coefficients.
- **Substrate content:** f_π is substrate-derived (½σξ). g_A is empirical.
  m_N empirical. The whole pipeline (GT relation, OPE potential, SEMF) is
  inherited from nuclear physics.
- **Genuine substrate prediction:** None — substrate inputs into a textbook
  nuclear-physics pipeline.

### 14. `src/stiff_medium/neutron_beta_decay.py` — neutron lifetime
- **Category:** **C**
- **What it actually computes:** Standard V-A formula `Γ_n = G_F²|V_ud|²
  (1+3g_A²) m_e^5 f(ε₀) / (2π³)` with phase-space integral f from the
  Sirlin closed form. Substrate inputs: m_e, Δm_np = m_n - m_p, ε₀.
- **Substrate content:** Δm_np uses the substrate prediction (1.293 MeV
  comes from §18.63.3 nucleon mass splitting, also substrate-B). m_e is
  substrate fundamental. SU(6) prediction g_A = 5/3 is substrate.
- **Genuine substrate prediction:** None — V-A formula and Sirlin
  phase-space integral are textbook QFT. The module honestly says: "the
  only piece coming from the collider/electroweak sector is the dimensionful
  Fermi coupling G_F itself."

### 15. `src/stiff_medium/bbn_predictions.py` — BBN abundances
- **Category:** **C**
- **What it actually computes:** Friedmann radiation-era H(T), weak rate
  Γ(n↔p) ~ G_F²(1+3g_A²)T^5, freeze-out condition H = Γ → T_F, n/p exp(-Δm/T_F),
  Y_p = 2(n/p)/(1+(n/p)).
- **Substrate content:** Δm is substrate (§18.63.3). m_e, α substrate.
  Everything else (FRW, weak rate, freeze-out condition) is standard
  cosmology/QFT.
- **Genuine substrate prediction:** None — this is BBN with substrate Δm
  plugged in. The module is honest: "G_F (electroweak, collider sector —
  outside §18.75.1 natural domain)".

### 16. `src/stiff_medium/heavy_quarkonium.py` — charmonium / bottomonium
- **Category:** **A** (eigenproblem) + **C** (hyperfine)
- **What it actually computes:** Solves the radial Schrödinger equation in
  the substrate Cornell potential `V(r) = σr − (4/3) α_M / r` numerically
  (sparse Lanczos, N_r = 400) for 1S, 2S, 1P levels. Anchors m_c, m_b from
  spin-averaged 1S. Predicts ψ(2S), χ_cJ centre-of-gravity. Hyperfine
  splittings from `(32π/9) α_M |ψ(0)|² / m_q²`.
- **Substrate content:** The Cornell potential `σr − (4/3) α_M /r` has
  substrate origin: σ from K(ξ), α_M = σξ² from Möbius coupling, C_F=4/3
  from §18.49 SU(3). The radial Schrödinger eigenproblem is genuinely solved
  on a substrate-determined potential.
- **Genuine substrate prediction:** The 2S, 1P-cog energies in the substrate
  Cornell potential — this is real category-A work. The hyperfine formula
  `(32π/9) α |ψ|² / m²` is the inherited one-gluon-exchange Breit-Pauli
  form (C).

### 17. `src/stiff_medium/foot_phase_derivation.py` — Foot Q=2/3
- **Category:** **A** (Foot ↔ Koide structure) + honest "no derivation" for δ
- **What it actually computes:** `m_n = M (1+√2 cos(2πn/3 + δ))²` and
  enumerates topological candidates for δ. Honestly verifies Q=2/3 is
  structural to the parametrization for ANY δ. Best topological match for
  δ ≈ 1.404π is `7π/5 = 1.4π` (within 0.3% numerically) but has no clean
  topological derivation.
- **Substrate content:** The Foot/Koide algebra (Q=2/3, √2 cone amplitude,
  Z₃ vertex closure) is laid out cleanly with §18.10/18.30 substrate
  motivation. The 45° cone connection (cos 45 = 1/√2) is substrate-A.
  However δ = 1.404π itself is NOT derivable from these primitives alone —
  the module is explicit about this.
- **Genuine substrate prediction:** Q = 2/3 is automatic. δ is honestly
  flagged as not yet derived.

### 18. `src/stiff_medium/mobius_dirac_vertex.py` — Möbius Dirac eigenvalues
- **Category:** **A**
- **What it actually computes:** Builds the squared-Dirac operator
  `H = -∂_r² + ℓ_eff(ℓ_eff+1)/r² + m(r)² + m'(r)` with Möbius half-flux
  ℓ_eff ∈ {-1/2, 1/2, 3/2} and kink mass profile m(r) = M tanh(r/ξ). Solves
  three radial eigenvalue problems via sparse shift-invert eigsh (N_r ≤ 500).
  Returns the three κ_n eigenvalues, normalises, asks if they match Foot.
- **Substrate content:** Pure substrate. The radial Dirac operator on a
  half-flux Möbius bundle with kink mass profile is the §18.30 program
  literally. The eigenproblem output is substrate geometry.
- **Genuine substrate prediction:** The κ_n eigenvalues themselves — three
  numbers from a Möbius-Dirac eigenproblem.

### 19. `src/stiff_medium/wilson_fermions.py` — Wilson lattice for Jackiw-Rebbi
- **Category:** **A**
- **What it actually computes:** Builds the 3D Wilson-Dirac operator D_W as
  a sparse complex CSR matrix on an N³ lattice (4·N³ DOF), with kink mass
  profile m(x) = M tanh(x/ξ), Wilson term r/(2a) ∇² to suppress doublers.
  Diagonalises H² = D_W† D_W via sparse Lanczos to extract the lowest
  eigenvalue (Jackiw-Rebbi zero-mode).
- **Substrate content:** Pure substrate. This is solving the Dirac equation
  in a kink background on a 3D lattice. PBC, APBC, and Dirichlet boundary
  conditions are systematically tested. Output: the JR zero-mode and its
  approach to zero as N increases.
- **Genuine substrate prediction:** The lattice eigenvalue λ_0 → 0 in the
  continuum limit — this is direct substrate-field-equation output.

### 20. `src/stiff_medium/multi_kink.py` — multi-kink configurations
- **Category:** **B** (kink-antikink potential is substrate; composite
  binding is a phenomenological pairwise sum)
- **What it actually computes:** Defines Kink and MultiKinkComposite
  dataclasses (positions, charges, chiralities). Uses the §18.48.3 pairwise
  potential `V(R) = -32(K/ξ) e^(-R/ξ)` (large R) / `+8K/ξ ln(R/ξ)` (small R)
  to compute pairwise binding energies. Constructs DM dimer, baryon 3-kink,
  meson 2-kink composites.
- **Substrate content:** The kink-antikink potential V(R) is substrate-
  derived (sech² kink overlap). Pairwise sum + zero-point correction is a
  semi-empirical composite model — not fully derived field theory.
- **Genuine substrate prediction:** The pairwise V(R) shape and its
  equilibrium length R_eq ≈ ξ. The composite assembly is mostly
  bookkeeping.

---

## Summary

### Tally by category
| Category | Count | Modules |
|---|---|---|
| **A** (real substrate geometry) | **5** | substrate_rg_running, mobius_dirac_vertex, wilson_fermions, foot_phase_derivation (Q=2/3 piece), heavy_quarkonium (radial Cornell eigenproblem) |
| **B** (substrate-derived formula) | **6** | chromomagnetic_substrate, glueball_closed_string, pion_decay_constant, strange_quark_sector, multi_kink, regge_spectrum (substrate σ side) |
| **C** (inherited formula + substrate input) | **6** | nucleon_magnetic_moments, atomic_transitions, nuclear_binding, neutron_beta_decay, bbn_predictions, hyperon_spectrum |
| Mixed / honest "open" | **3** | proton_radius_v3 (A spread + C composition), isospin_breaking (B/C), cabibbo_angle (no derivation) |

### Counts after assigning mixed modules to dominant category
- **Category A**: 6 (mobius_dirac_vertex, wilson_fermions, substrate_rg_running, foot_phase_derivation, heavy_quarkonium eigensolver, proton_radius_v3 spread)
- **Category B**: 6 (chromomagnetic_substrate, glueball_closed_string, pion_decay_constant, strange_quark_sector, multi_kink, regge_spectrum)
- **Category C**: 7 (nucleon_magnetic_moments, atomic_transitions, nuclear_binding, neutron_beta_decay, bbn_predictions, hyperon_spectrum, isospin_breaking)
- **Honest open**: 1 (cabibbo_angle)

### Honest reading

The framework's "substrate is producing the geometry" claim is strongest in
modules that explicitly solve substrate field equations:

- **wilson_fermions** and **mobius_dirac_vertex**: directly diagonalise
  Dirac operators on substrate backgrounds.
- **substrate_rg_running**: an actual RG flow ODE.
- **proton_radius_v3** and **heavy_quarkonium**: radial Schrödinger
  eigenproblems in substrate-derived potentials.

The **B-class** modules are intermediate — they produce closed-form
formulas (`f_π = ½σξ`, `Δm_NΔ = 4 ξ² σ^{3/2}`, `f_K/f_π = cosh(½ ln(m_K/m_π))`,
`α'_closed = 3/4 α'_open`) that are dimensionally and structurally
substrate-mechanical. These are the "best-case" predictions: the formula
itself encodes substrate geometry, not just the input numbers.

The **C-class** modules are honest about what they are: substrate-derived
inputs (α, m_e, m_p, μ_p, f_π, Δm_np) plugged into known QFT/QM pipelines
(Bohr formulas, Goldberger-Treiman, V-A neutron decay, Friedmann +
freeze-out for BBN). For these, the substrate is producing the inputs but
not the formula structure.

### What this means for the "28+ predictions" claim

- A meaningful fraction (~6 modules, ~30%) is real substrate geometry
  emerging from substrate dynamics or eigenproblems.
- A second similar fraction (~6 modules, ~30%) is substrate-derived
  formulas. These are non-trivial substrate predictions.
- A roughly equal fraction (~7 modules, ~35%) is inherited formulas with
  substrate inputs. These are LESS strong as evidence that substrate
  "produces the geometry," but they ARE meaningful tests that the
  substrate's input values (α, m_e, m_p, f_π, Δm_np) are consistent with
  observation.
- 1 module honestly flags itself as not yet derived (Cabibbo).

The cross-scale consistency claim is real — the same σ, ξ, α, m_e drive
predictions from atomic spectra to BBN — but the framework's "substrate
geometry" is best in the A modules and intermediate in the B modules.
The C modules are properly described as "substrate inputs into
inherited physics" rather than "substrate produces the result."
