# The Stiff Medium Model

A working substrate-mechanical framework built on a single 3D elastic medium. Stable matter, electromagnetism, gravity, particle physics, and cosmology are treated as patterns of one medium, with several sectors already quantitative and several still open.

**Read first:** [`MODEL.md`](MODEL.md) — the unified theory in a single document.

**Full spec:** [`docs/superpowers/specs/2026-04-29-stiff-medium-theory-design.md`](docs/superpowers/specs/2026-04-29-stiff-medium-theory-design.md) — derivations, tests, successes, boundary results, and the current weak-sector stress tests through §18.95.

---

## The encompassing Lagrangian

The current candidate unifying expression (§18.45):

```
ℒ = ½ρ(∂_t φ)² − ½K|∇φ|² − V(φ)              [substrate]
  + ψ̄(iℏγ^μ(∂_μ + ieA_μ) − g_Y φ)ψ           [fermion + EM + Yukawa]
  − ¼ F_μν F^μν                              [bundle]
  + λ(x)[(∂_zφ)² − |∇_⊥φ|²]                 [45° cone]
```

with `V(φ) = (K/ξ²)(1 − cos(φ/ξ))/√(1 − (φ/φ_max)²) − ε_0` and Möbius half-flux topology.

The current compact-geometric candidate (§18.84) rewrites the cone multiplier as substrate null geometry and folds the dark sector into a symmetric-traceless neutral-stress term:

```
ℒ_geo = ½ρ(D_tφ)² − ½K g_cone^ij D_iφ D_jφ − V(φ)
      + ψ̄(iℏγ^a e_a^μD_μ − g_Yφ)ψ − ¼F_A²
      − ½Kα² Tr(ST(strain)²)
```

with `D = d + ieA_EM + iA_Möbius`. This is a candidate simplification, not a completed derivation.

**10-20 working parameters**, depending on how many effective-sector terms are counted (vs ~30 in SM + ΛCDM + GR). The goal is to derive the observed sectors from this substrate structure; flavour mixing, Planck-scale closure, and precision cosmology remain active gaps.

See `scripts/encompassing_lagrangian.py` for structure visualization.

## Status

The strongest benchmark checks are in atomic physics, QED inheritance, gravity/EM hierarchy, strong-field GR limits, and QCD-scale mechanics. Later sections also document failures and boundaries, especially in lepton hierarchy, CKM/flavour mixing, Planck-scale UV closure, and one-shot Big-Bang baryogenesis.

| Domain | Representative checks | Best agreement |
|---|---|---|
| Atomic / chemistry | 8 core checks | Hydrogen E_1s exact |
| Universal physics | 5 core checks | Gravity/EM ratio 0.06% |
| Strong-field GR | 5 benchmark checks | Mercury precession <1% |
| Particle physics / QED | 8 benchmark checks | Michel parameters 0.01% |
| QCD-scale mechanics | growing set | f_π, f_K, proton radius, magnetic moments |
| Open/boundary sectors | active gaps | CKM, lepton hierarchy, Planck scale, CMB/Hubble |

The compression claim is now treated as provisional: the substrate core is compact, but open sectors must be closed without adding ad-hoc parameters.

---

## What's in this repo

```
MODEL.md          ← Unified theory document (start here)
docs/             ← Full spec with derivations, verifications, history
src/stiff_medium/ ← Core simulation modules
scripts/          ← Demonstration & verification scripts
tests/            ← Unit tests for core modules
```

### Core simulation modules (`src/stiff_medium/`)

- `neutrino.py`, `three_d.py`: 45° cone propagation primitive
- `dynamics.py`: time evolution
- `back_reaction.py`: medium response forces
- `mobius_dynamics.py`: half-flux topology / spin-½
- `atomic.py`: multi-electron N-body Coulomb + Pauli + EM damping
- `spinor.py`: Möbius spinor structure
- `em_field.py`, `em_field_3d.py`: EM wave propagation (1D, 3D)
- `detector.py`: bound-state tracking
- `visualize.py`: matplotlib animation

### Key demonstration scripts

**Atomic physics:**
- `helium_variational.py` — variational He, Li⁺, Be²⁺ ground states
- `hartree_radial.py` — Madelung 4s/3d ordering for K
- `madelung_rule.py` — Slater rules across periodic table
- `magnesium_screened.py` — heavy atom (Z=12) with orbital screening
- `h2_lcao.py` — H₂ molecular bonding via LCAO-MO
- `atomic_emission_spectroscopy.py` — Lyman-α and emission lines

**Universal physics:**
- `mass_energy_equivalence.py` — E = mc² as kinematic identity
- `gravity_static_deflection.py` — 1/r² law from substrate Poisson
- `g_from_substrate.py` — gravity/EM hierarchy at 0.06%
- `cone_bouncing_mass.py` — §18.35 cone-bouncing mass mechanism

**Strong-field GR:**
- `strong_field_gravity.py` — light bending, Mercury precession, GPS, Pound-Rebka, BH horizon

**Particle physics:**
- `lepton_dirac_solver.py` — single-kink Dirac spectrum
- `multi_kink_dirac.py` — multi-kink configurations
- `lepton_koide_in_model.py` — Koide relation in our model
- `muon_decay_spectrum.py` — Michel spectrum, V-A coupling
- `precision_qed_tests.py` — tau lifetime, electron g-2, Lamb shift, 21cm

**Cosmology:**
- `cosmology_numerical.py` — dark matter, dark energy, BH formation, inflation
- `cyclic_cosmology_timescales.py` — end-state cycle lengths
- `cmb_phase_transition.py` — CMB as substrate de-saturation

**Audit / gap tests:**
- `alpha_audit.py`, `tests/test_alpha_audit.py` — canonical fine-structure constant audit; current conclusion is no derivation yet, with the best fixed two-loop scan point low by `0.528795527` in inverse-alpha (`0.385880739%`) vs CODATA 2022
- `dependency_ledger_report.py` — dependency/status ledger for model claims
- `missing_piece_hypotheses_report.py` — weak-sector hypothesis tests for UV, leptons, CKM, cosmology transfer, dark matter, and matter orientation
- `mechanism_trials_report.py` — concrete mechanism trials for the current weak sectors
- `substrate_polarization_dm_test.py` — strict no-DM/pure-polarization boundary test
- `dark_stress_hybrid_test.py` — neutral-kink / substrate-polarization hybrid dark-stress test
- `dark_stress_parameter_closure_test.py` — dark-stress abundance, mobile fraction, memory, and halo-radius closure candidates
- `dark_stress_memory_clock_test.py` — polarization-memory clock trials
- `dark_stress_scale_closure_test.py` — derives the kpc coherence clock from `α³(c/H0)/√3` and `αc/√5`
- `dark_stress_factor_scan_test.py` — checks whether the dark-scale factors are unique or numerological
- `dark_stress_cluster_dynamics_test.py` — cluster lensing split, memory offset, and dark-stress causal-horizon audit
- `neutral_stress_tensor_modes_test.py` — symmetric-traceless stress projector and `√5` dark-speed mode count
- `dark_stress_transport_test.py` — finite-speed 1D mobile-kink / locked-polarization transport profile
- `neutral_coupling_suppression_test.py` — narrows `α` dark-speed suppression to `K_eff/K = α²`
- `dark_stress_em_darkness_test.py` — operational EM-darkness gate: no emission, absorption, or reflection channel
- `geometric_action_compaction_test.py` — compact candidate action with cone metric, Möbius connection, and neutral stress
- `cone_variational_origin_test.py` — equal-partition elastic penalty selecting the 45° cone
- `cone_lattice_microgeometry_test.py` — lattice-invariant audit showing the cone quartic requires self-dual exchange
- `cone_self_dual_exchange_test.py` — paired dual-branch mechanism for cancelling cone bias conditionally
- `cone_detailed_balance_test.py` — local detailed-balance route to exact 50/50 dual-branch weights
- `cone_swap_generator_origin_test.py` — elastic-cell automorphism route to the swap-degenerate cone generator
- `cone_diamond_cell_geometry_test.py` — saturated diamond spring-cell candidate for the branch-swap automorphism
- `cone_diamond_cell_selection_test.py` — minimal graph audit for when the diamond cell is uniquely selected
- `cone_anchor_induced_exchange_test.py` — Schur-complement route from shared finite anchors to effective L-T exchange
- `cone_two_anchor_origin_test.py` — phase-slip endpoint neutrality route to paired finite anchors
- `cone_phase_slip_lattice_test.py` — discrete lattice phase-slip segment and stiffness-ratio robustness audit
- `cone_saturated_bond_selection_test.py` — energetic check showing the barrier alone delocalizes phase slip

**EM:**
- `em_propagation_test.py`, `em_3d_test.py` — wave propagation (1D, 3D)
- `em_3d_spectroscopy.py` — 3D atomic absorption with geometric falloff

---

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Run any demonstration:
```bash
python scripts/g_from_substrate.py        # gravity from substrate (0.06% match)
python scripts/strong_field_gravity.py    # 5 GR tests
python scripts/precision_qed_tests.py     # QED precision tests
python scripts/cmb_phase_transition.py    # CMB reinterpretation
```

Run tests:
```bash
pytest
```

---

## What the model claims

The substrate is the only fundamental thing. Everything else is intended to emerge as a stable or transient substrate pattern:

```
Substrate
  │
  ├─ Charge-asymmetric channel ──→ Electromagnetism
  ├─ Charge-symmetric channel ───→ Gravity
  ├─ Cone-bouncing frequency ────→ Mass
  ├─ Möbius half-flux ───────────→ Spin-½, Pauli exclusion
  ├─ Local saturation σ=½ ───────→ Black holes (no singularity)
  ├─ Universe-wide saturation ───→ Long bleed-off era, not a singular beginning
  ├─ De-saturation phase shift ──→ CMB mark + clean radiation/matter split
  ├─ Multi-kink composites ──────→ Dark matter
  ├─ Baseline strain σ₀ ─────────→ Dark energy
  └─ Saturation/dissipation ─────→ Cosmic cycles
```

**One mechanism, many phenomena.** The model aims to unify what standard physics treats as separate: EM and gravity, wave and particle behavior, lepton excitations, black-hole saturation, CMB de-saturation, and cosmic cycles.

---

## Open problems

The current high-value gaps:

1. Numerical α from the substrate Lagrangian
2. Specific lepton excitation energies (Δ₁, Δ₂), now reframed as vertex eigenvalue ratios `κ_μ/κ_e ≈ 4.28e4` and `κ_τ/κ_e ≈ 1.21e7`
3. CKM/PMNS and the flavour-mixing operator, not more angle numerology
4. Current-quark masses and SU(3)-breaking renormalization
5. Matter-sector orientation selection / inheritance, not one-shot Big-Bang baryogenesis
6. Planck-scale UV completion and the missing dimensionless input `χ_UV ≈ 4.2e-23`
7. Saturated bleed-off law, critical de-saturation threshold, `f_vis <= 4e-4`, and full CMB/Hubble fit from derived transfer windows
8. Quantitative dark substrate-stress sector: derive the neutral phase-space measure, mobile split, `α³/√3` coherence filter, second-order neutral stiffness `K_eff/K = α²`, coupled `ρ_kink/ρ_pol` dynamics, and full lensing-map signatures; current finite-speed transport and EM-darkness audits require heavy neutral mobile stress plus ultra-low-frequency coherent polarization, not pure polarization
9. Strong-field GR full nonlinear regime

These are bounded problems, but several may require additional substrate physics rather than only more algebra.

Current audit: strongest in QCD-scale mechanics and standard-physics containment; weakest in lepton hierarchy, CKM/flavour mixing, Planck/UV closure, precision cosmology transfer functions, and dark matter quantification.

---

## Philosophy

> **The medium is the eternal entity; matter patterns are temporary.**

The substrate is more fundamental than any matter configuration. Particles, atoms, stars, galaxies — even our entire universe — are temporary patterns of strain in an eternal elastic medium. Cosmic cycles don't reset the substrate; they reset the patterns ON it.

The substrate parameters (K, ρ, ξ, ...) are properties of the eternal medium, constant across all cycles. Anthropic selection isn't needed — every cycle has the same physics.

---

*See `MODEL.md` for the current working framework.*
