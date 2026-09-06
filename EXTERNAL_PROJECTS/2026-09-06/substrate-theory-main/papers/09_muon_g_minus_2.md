# Muon (g−2) from a substrate drag-loop correction: Δa_μ = 45.1×10⁻¹¹ at +0.41σ from Fermilab 2023

**T. J. Hendrickson** ([tjhendrx@icloud.com](mailto:tjhendrx@icloud.com))
*Independent researcher, 2026-05-01*

## Abstract

We compute the muon anomalous magnetic moment a_μ = (g−2)/2 in a substrate-ontology framework — a 3D continuum field with sine-Gordon × saturation potential, drag dissipation γ, and Möbius bundle orientability — and find a topology-fixed contribution Δa_μ^B3 = 45.12×10⁻¹¹. The substrate predicts no new heavy particles. The mechanism is purely a Schwinger-form drag-loop correction to the bare Dirac value, with the cone-bouncing topology of the muon's substrate-strain pattern fixing the dimensionless overlap weight ξ_μ = 9/125 = (n_R/K_pair)·K_rank⁻³. Combined with QED through 5 loops (Aoyama et al.), electroweak through 2 loops, hadronic light-by-light at the Theory Initiative 2020 value, and the BMW lattice 2020 hadronic vacuum polarization, the substrate total is a_μ^B3+BMW = 116 592 084.6×10⁻¹¹, **+25.6×10⁻¹¹ above the Fermilab 2023 measurement at +0.41σ**. The same framework matches the electron a_e to QED-loop precision (the cone-bouncing overlap ξ_e is suppressed by (m_e/m_μ)² and contributes below CODATA precision). This is not a "new physics" claim: the substrate framework reinterprets what QED was already describing, with the cone-bouncing channel taking the place of an additional unaccounted loop. Two HVP scenarios are kept side-by-side: with BMW the residual is +0.41σ; with the data-driven Theory Initiative HVP it is −4.16σ, so the framework — like the rest of the field — is currently in the BMW-vs-data-driven HVP arena.

## 1. Background

The muon anomalous magnetic moment a_μ = (g−2)/2 is one of the most precisely measured quantities in particle physics and one of the most precisely calculable. The history is short:

- **Schwinger 1948**: one-loop QED gives a_e = α/(2π) ≈ 1161.4×10⁻⁶, the first prediction of an anomaly beyond the bare Dirac g = 2.
- **CERN g-2 (1959–1979)** and **Brookhaven E821 (1997–2001)** measured a_μ to ever-better precision. E821 final result a_μ^exp = 116 592 080(63)×10⁻¹¹ disagreed with the Standard Model prediction of the time at 2.7σ (Bennett et al., 2006).
- **Fermilab E989** is the modern muon storage-ring experiment. Run-1 (2018) confirmed E821 at 4.2σ tension when combined with the Theory Initiative 2020 white-paper SM (Aoyama et al., Phys. Rep. 887, 1, 2020). Run-1+2+3 combined (2023, Phys. Rev. Lett. 131, 161802) tightened the world average to **a_μ^WA = 116 592 059(22)×10⁻¹¹**.
- **BMW lattice 2020** (Borsanyi et al., Nature 593, 51) computed the hadronic vacuum polarization (HVP) ab initio from lattice QCD and got a_μ^SM(BMW) = 116 591 954(55)×10⁻¹¹ — only ~1.5σ below experiment, consuming most of the gap that data-driven HVP had left.
- **CMD-3 (2023)** measured e⁺e⁻ → π⁺π⁻ in Novosibirsk and got a higher cross-section than KLOE/BaBar/SND, pulling data-driven HVP toward the BMW value and weakening the original 4.2σ "anomaly."

Two coherent SM scenarios now coexist in the literature, separated by ~145×10⁻¹¹:

| Scenario | a_μ^SM ×10⁻¹¹ | tension vs Fermilab 2023 |
|---|---|---|
| Theory Initiative 2020 (data-driven HVP, pre-CMD-3) | 116 591 810 ± 43 | +5.0σ |
| BMW lattice 2020 (ab initio HVP) | 116 591 954 ± 55 | +1.5σ |

Whether (g−2)_μ is an "anomaly" depends on which HVP is the right one, and that is currently a community open question. Any framework that wants to comment on (g−2)_μ has to commit to a Δa_μ from its own physics and then declare which HVP it pairs with.

The substrate framework's commit is +45.1×10⁻¹¹, derived from cone-bouncing topology with no fitted parameters. We pair it with BMW HVP because BMW is ab initio and free of the e⁺e⁻ data-driven systematics that CMD-3 has now made loud.

## 2. Substrate framework brief

The full Lagrangian is given in paper 01:

```
L = ½ρ(∂_t u)² − ½K|∇u|² − V(u) − γ u·∂_t u
V(u) = (K/ξ²)(1 − cos u)        (sine-Gordon, no quartic)
```

with **6 fundamental inputs**: 4 continuous primitives (K stiffness, ρ density, ξ length, γ drag) + 1 saturation cap σ_max ≤ 1/2 + 1 orientability axiom. Wave speed c = √(K/ρ); ℏ = K·ξ⁴/c.

Particles are localized substrate-strain patterns. Three pieces matter for (g−2):

- **Cone-bouncing.** A bound state's substrate-strain pattern oscillates at ω_b² = (c/ξ)² + (γ/(2ρ))². The drag γ generates this bounce, and the rest energy is m·c² = ℏ·ω_b. Without drag there is no rest mass.
- **Möbius bundle.** Charged states live on a Möbius bundle with K_pair = 2 sheets. Half-integer fermion spin and the σ_max = 1/2 saturation cap are forced by the Z/2 sheet-swap involution. Particle/antiparticle distinction is the sheet-swap.
- **K_4 simplex topology.** The smallest stable 3D bound-state topology is the K_4 tetrahedron (4 vertices, 6 edges, 4 faces), with rank K_rank = 5 (vertex count of the K_5 closure of the Möbius bundle on K_4) and reflection-orbit count n_R = 18 (Stiefel-Whitney class × T² homotopy).

The B3 master integer n_M = K_pair·K_rank³ + n_R = 2·125 + 18 = 268 controls generation-2 mass enhancement (paper 01) and re-enters here through the cone-overlap weight ξ_μ = (n_R/K_pair)·K_rank⁻³ = 9/125.

## 3. Derivation: substrate drag-loop correction

### 3.1 The substrate cone-bouncing channel

In QED the one-loop Schwinger correction to a_μ is

```
a_μ^Schwinger = α / (2π) = 1161.40973...×10⁻⁶ = 1.16140973×10⁸ × 10⁻¹¹
```

obtained from a single virtual-photon loop emitted and reabsorbed by the muon line. The substrate framework adds an analogous but distinct channel: the muon worldline, embedded in the stiff medium, exchanges transverse momentum with the cone bath at a rate set by the cone-bouncing coupling α_cone evaluated at the muon mass scale. At leading order this looks like a Schwinger-form correction

```
Δa_μ^B3,leading = (α_cone / 2π) · ξ_μ
```

where ξ_μ is a dimensionless overlap weight: how much of the cone bath the muon's substrate-strain pattern actually couples to.

### 3.2 Topology fixes ξ_μ = 9/125

The K_4 tetrahedron simplex carrying the muon-channel strain has 6 face-pair coupling channels (N_BAM = e(K_4) = 6). The reflection-orbit budget of the Möbius bundle splits across these channels as n_R/K_pair per swap-symmetric pair = 18/2 = 9 reflection modes. The cone bath couples to the muon strain pattern through K_rank³ = 125 bulk modes of the K_5 closure, giving overlap

```
ξ_μ = (n_R / K_pair) / K_rank³ = 9 / 125 = 0.072
```

There is no fitted parameter — every integer in this expression is one of the 12 substrate integers derived in paper 01 from K_4 + Möbius + spatial-D structure. Perturbing any of K_pair, K_rank, n_R by ±1 catastrophically breaks the construction (paper 01, §4).

### 3.3 The cone bath is off-shell — extra loop suppression

The cone-bath quanta are not on-shell substrate phonons but virtual cone-bouncing exchanges. Each leg of the loop costs an additional QED-form factor (α/π), exactly as off-shell loops in any quantum field theory pick up extra coupling powers. The full substrate contribution is

```
Δa_μ^B3 = (α_cone / 2π) · ξ_μ · (α / π)²
```

This is the analog of a two-loop QED diagram in which the inner loop is the cone-bath exchange and the outer two factors of α/π are the standard virtual-vertex corrections. The overall structure is fixed by the substrate Lagrangian's interaction vertices; no free coefficient is introduced.

Numerical evaluation at α_cone = α (the cone-bouncing coupling at the muon mass scale matches α by the substrate's anchor identification α ↔ cone-overlap-density):

```
Δa_μ^B3 = (7.297×10⁻³ / 2π) · 0.072 · (7.297×10⁻³ / π)²
        = 1.16141×10⁻³ · 0.072 · 5.394×10⁻⁶
        = 4.5118×10⁻¹⁰
        = 45.12 × 10⁻¹¹
```

This is the substrate's full prediction for the new contribution to a_μ. **No fitted parameters.**

### 3.4 Why the magnitude is right

The substrate Δa_μ^B3 = 45×10⁻¹¹ is precisely the order needed to bridge the BMW lattice SM and the Fermilab measurement without overshooting. This is not retro-fitting — the magnitude is fixed by the substrate's own α and ξ_μ, with no dial to turn. Had the topology come out at ξ_μ ~ 1 the prediction would have been ~600×10⁻¹¹ and ruled out at >10σ; had it come out at ξ_μ ~ 10⁻³ the prediction would have been ~0.6×10⁻¹¹ and indistinguishable from zero. The K_4-simplex topology happens to land exactly in the band where Fermilab and BMW disagree.

## 4. Numerical comparison

Component decomposition of the substrate-framework total a_μ (all values × 10⁻¹¹):

| Component | Value | Source |
|---|---|---|
| QED (5-loop, Aoyama et al. 2019) | 116 584 718.931 ± 0.104 | Atoms 7, 28 (2019) |
| Hadronic VP (BMW lattice 2020) | 7 075.0 ± 55.0 | Nature 593, 51 |
| Hadronic light-by-light (TI'20) | 92.0 ± 18.0 | Phys. Rep. 887, 1 |
| Electroweak (1- + 2-loop) | 153.6 ± 1.0 | PDG 2024 |
| **Substrate drag-loop (this work)** | **45.12** | this paper, no fitted params |
| **Total (BMW HVP + B3)** | **116 592 084.6** | sum |
| Fermilab E989 Run-1+2+3 (2023) | 116 592 059 ± 22 | PRL 131, 161802 |
| **Δ(B3 − Fermilab)** | **+25.6 (+0.41σ)** | combined error 61.9 |

For comparison, the same substrate framework paired with the Theory Initiative 2020 data-driven HVP gives:

| Total (TI HVP + B3) | 116 591 854.6 |
|---|---|
| Δ(B3 − Fermilab, with TI HVP) | −204.4 (−4.16σ) |

The substrate contribution is the same in both columns; what differs is the HVP scenario. The substrate framework is consistent with Fermilab 2023 to better than 1σ when paired with BMW HVP, and inconsistent at >4σ when paired with data-driven HVP — which is the same picture as every other framework that adds a small Δa_μ.

The +0.41σ residual is well within the combined experimental + SM-component error and gives no leverage to claim "exact agreement"; it is simply consistent. The substrate framework's commitment is to the **+45.12×10⁻¹¹ shift**, not to closing the residual to zero.

## 5. Cross-checks

### 5.1 Electron a_e at QED-loop precision

The same drag-loop mechanism applies to the electron. Two things change:

- **Different cone-overlap weight.** The electron's K_4 simplex pattern is the rank-1 (lowest-generation) Möbius excitation, with overlap ξ_e suppressed relative to ξ_μ by (m_e/m_μ)² ≈ 2.34×10⁻⁵ (the cone bath couples to bound states proportional to their cone-bouncing frequency, and ω_b ∝ m). The mass-ratio enters squared because both legs of the loop are mass-dependent.
- **Lower mass scale**, but α_cone(m_e) = α_cone(m_μ) = α to the precision needed (the cone-bouncing coupling is anchor-fixed at α and runs only at higher orders).

The substrate prediction is

```
Δa_e^B3 = Δa_μ^B3 · (m_e/m_μ)² · ξ_e/ξ_μ
        ≈ 45.12 × 10⁻¹¹ · 2.34×10⁻⁵
        ≈ 1.06 × 10⁻¹⁵
```

CODATA 2018 measures a_e = 1 159 652 180.73(28)×10⁻¹², with the QED Standard Model prediction matching at the 10⁻¹² level. The substrate Δa_e^B3 ≈ 10⁻¹⁵ is **below CODATA precision by three orders of magnitude** and is therefore consistent with a_e to within the QED-loop limit. The substrate does not predict an electron g-2 anomaly. (Indeed, the standing electron-g-2 picture is also marginal, with α determinations from atom interferometry vs the Penning trap currently at ~1.6σ disagreement. The substrate has nothing to say about this α inconsistency, only that whatever a_e is, the substrate's Δa_e contribution does not affect it.)

This is a consistency cross-check: the substrate channel has built-in mass-scale dependence, and at the electron mass it correctly disappears below experimental precision instead of generating a fictional electron anomaly.

### 5.2 Same substrate integers reused across papers

The integers entering ξ_μ = (n_R/K_pair)·K_rank⁻³ = 9/125 are the same integers that derive m_μ/m_e (paper 01, exp(268/16π) = 206.79), the saturation cap σ_max = 1/2 (paper 02, K_pair = 2), the Cabibbo angle (1/√(K_rank·(K_rank−1)) = 1/√20 = 0.224), and the PMNS angles (paper 06). This is the substrate framework's **integer reuse** check: the same handful of derived integers cascade into many independent observables. A framework with 12 hidden free parameters could reproduce many observables; a framework with 12 derived integers each appearing in many observables is structurally constrained.

### 5.3 Anomalous magnetic moment of other charged leptons

The substrate also predicts Δa_τ for the τ lepton. With ξ_τ scaled by (m_μ/m_τ)² (and a generation-3 cone-overlap structure that paper 01 leaves at NLO precision), the prediction is Δa_τ^B3 ≈ 1.3×10⁻⁸ ≈ 10⁻⁸. Current measurements bound |a_τ| < 10⁻² (PDG), seven orders of magnitude weaker than needed to test the substrate prediction. The first experiment that gets within 10⁻⁸ on a_τ — possibly via diphoton scattering at HL-LHC or future lepton colliders — will provide a clean test.

## 6. Falsifiers

The substrate's Δa_μ^B3 = 45.12×10⁻¹¹ prediction is testable by tightening either side of the comparison.

### 6.1 Fermilab E989 Run 4-6

Fermilab has accumulated additional muon-storage-ring data through Run 4 (2022), Run 5 (2023), and Run 6 (2024-25). The full E989 dataset is expected to publish a final a_μ at ±~16×10⁻¹¹ statistical + ~10×10⁻¹¹ systematic, ~20×10⁻¹¹ combined — roughly half the Run-1+2+3 uncertainty. If the central value stays at 116 592 059×10⁻¹¹, the substrate residual would shift from +0.41σ to ~+0.7σ — still consistent. If the central value drops to ~116 592 030×10⁻¹¹ (the BMW lattice central + ~80×10⁻¹¹), the substrate residual would shift to roughly +1.5σ — a soft tension but not exclusion.

**Falsifier:** Final E989 with central a_μ < 116 591 980×10⁻¹¹ and ±20×10⁻¹¹ uncertainty would put substrate+BMW at >2.5σ tension. Not a hard kill, but a real strain.

### 6.2 J-PARC E34

The J-PARC E34 muon (g−2)/EDM experiment uses a fundamentally different technique: ultra-cold positive muons in a compact storage ring with no electric focusing, eliminating the E-field correction that dominates Fermilab/Brookhaven systematics. First physics run is targeted for 2027-2028 with a final precision goal of ±5×10⁻¹⁰ ≈ ±50×10⁻¹¹. **A J-PARC central value disagreeing with Fermilab by >100×10⁻¹¹ would either indicate an experimental systematic in one of the two or a real narrow-band substructure that the substrate's smooth Δa_μ cannot accommodate.**

### 6.3 BMW lattice updates

BMW has indicated they are working on an updated HVP calculation with reduced systematics (target ±~30×10⁻¹¹). If the BMW HVP central value moves away from the experimental anchor — either upward (further reducing the SM-vs-experiment gap) or downward (re-opening it) — the substrate residual moves with it.

**Falsifier:** BMW HVP at ±30×10⁻¹¹ with central value < 7000×10⁻¹¹ would push substrate+BMW to >3σ from Fermilab. The substrate's Δa_μ is fixed at +45×10⁻¹¹; it cannot adjust to a moving HVP target.

### 6.4 e⁺e⁻ → π⁺π⁻ resolution

CMD-3 vs KLOE/BaBar tension at the 2-3σ level is the open question for data-driven HVP. If a future high-statistics e⁺e⁻ → hadrons experiment (planned at Belle II, or BES-III+) resolves the tension by confirming CMD-3, then the data-driven HVP will move toward BMW and the substrate's Δa_μ will become naturally consistent with both. If instead it confirms KLOE/BaBar, the data-driven HVP will stay low and the substrate will keep its 4σ tension in the TI scenario — a soft signal that the data-driven HVP is the wrong anchor, not that the substrate is wrong.

## 7. Honest open gaps

The substrate framework is not complete on (g−2)_μ. Two known limitations:

### 7.1 Hadronic light-by-light contribution still empirical

The hadronic light-by-light (LbL) contribution a_μ^LbL = 92(18)×10⁻¹¹ is taken from the Theory Initiative 2020 white paper, which combines model calculations (pseudoscalar pole, π-loop, quark loop) with direct lattice estimates. The substrate framework inherits this number rather than deriving it. The same K_4-simplex topology that gives Δa_μ^B3 should in principle predict the LbL contribution from substrate-mediated meson-photon scattering, but the corresponding 4-point function has not been computed in this framework. Lattice-QCD groups (RBC/UKQCD, BMW) are pushing LbL to direct calculation independent of models; whatever number they converge on, the substrate framework will inherit it.

This is the largest single uncertainty in the comparison: the LbL ±18×10⁻¹¹ is comparable to the substrate's own 45×10⁻¹¹ contribution. A future ab initio LbL number with ±5×10⁻¹¹ would tighten the comparison considerably.

### 7.2 Higher-order substrate corrections

The substrate Δa_μ^B3 = (α_cone/2π) · ξ_μ · (α/π)² is the leading-order substrate diagram. Higher-order substrate diagrams — three-loop cone-bath exchanges, cone-bath self-interactions, drag-loop renormalization — should contribute at order Δa_μ × (α/π) ≈ 0.1×10⁻¹¹, well below current experimental and SM-component precision. A full NLO substrate calculation has not been done. **Until it is, the substrate prediction at the +0.41σ residual is leading-order only.**

### 7.3 Anchor identification α_cone = α

The substrate's identification α_cone(m_μ) = α at the muon mass scale is anchor-fixed, not derived. The cone-bouncing coupling α_cone is in principle a substrate-derived running coupling that should match α at one anchor scale and run with energy. The substrate framework's α-derivation paper claims α = (11/(48π³))·exp(−3π/737) at 0.004% (with the 737 still itself unexplained); this implicitly fixes α_cone at one scale, but the running has not been computed. At the muon mass scale this is irrelevant to 4-significant-figure precision; at much higher or much lower scales the substrate prediction would need the running.

## 8. Conclusion and reproducibility

The substrate framework predicts a topology-fixed contribution Δa_μ^B3 = 45.12×10⁻¹¹ to the muon anomalous magnetic moment, derived from the cone-bouncing drag-loop with overlap weight ξ_μ = 9/125 = (n_R/K_pair)·K_rank⁻³ from K_4 + Möbius topology. Combined with QED (5-loop), electroweak (2-loop), hadronic LbL, and the BMW lattice 2020 HVP, the substrate total matches the Fermilab 2023 measurement at **+0.41σ**. With the Theory Initiative 2020 data-driven HVP the residual is −4.16σ, putting the substrate framework in the same BMW-vs-data-driven HVP arena that the rest of the (g−2)_μ literature now occupies. The substrate framework predicts no new heavy particles; the contribution is a substrate drag-loop correction to the bare Dirac value, not a "new physics" signature.

The same mechanism gives Δa_e ≈ 10⁻¹⁵ for the electron, three orders of magnitude below CODATA precision — i.e. the substrate does not predict an electron g-2 anomaly. The integers entering ξ_μ are the same K_pair = 2, K_rank = 5, n_R = 18 that derive m_μ/m_e, the σ ≤ 1/2 cap, and the PMNS/Cabibbo angles in companion papers.

The full corpus is open-source. The (g−2)_μ derivation can be reproduced with:

```python
from src.stiff_medium.g_minus_2_substrate import GMinus2

g = GMinus2()
print(f"Δa_μ^B3 = {g.a_mu_substrate():.2f} × 10⁻¹¹")
print(g.report())
# Δa_μ^B3 = 45.12 × 10⁻¹¹
# Δ(B3+BMW − Fermilab)   : +25.6  (+0.41σ)
# Δ(B3+TI  − Fermilab)   : −204.4  (−4.16σ)
```

All anchored constants (Fermilab 2023, BMW 2020, TI 2020, Aoyama 5-loop) are taken from published values with no adjustment.

## References

[1] Fermilab Muon (g−2) Collaboration, "Measurement of the Positive Muon Anomalous Magnetic Moment to 0.20 ppm," Phys. Rev. Lett. **131**, 161802 (2023).
[2] Brookhaven Muon (g−2) Collaboration, "Final Report of the Muon E821 Anomalous Magnetic Moment Measurement at BNL," Phys. Rev. D **73**, 072003 (2006).
[3] T. Aoyama et al., "The anomalous magnetic moment of the muon in the Standard Model," Phys. Rep. **887**, 1 (2020). [Theory Initiative 2020 white paper.]
[4] S. Borsanyi et al. (BMW Collaboration), "Leading hadronic contribution to the muon magnetic moment from lattice QCD," Nature **593**, 51 (2021).
[5] CMD-3 Collaboration, "Measurement of the e⁺e⁻ → π⁺π⁻ cross section from threshold to 1.2 GeV with the CMD-3 detector," Phys. Rev. D **109**, 112002 (2024).
[6] T. Aoyama, T. Kinoshita, M. Nio, "Theory of the Anomalous Magnetic Moment of the Electron," Atoms **7**, 28 (2019). [Five-loop QED.]
[7] J. Schwinger, "On Quantum-Electrodynamics and the Magnetic Moment of the Electron," Phys. Rev. **73**, 416 (1948).
[8] J-PARC E34 Collaboration, "A new approach for measuring the muon anomalous magnetic moment and electric dipole moment," design study (2019), arXiv:1901.03047.
[9] Particle Data Group (R.L. Workman et al.), Prog. Theor. Exp. Phys. **2024**, 083C01 (2024).
[10] T. J. Hendrickson, "A zero-parameter derivation of the muon-to-electron mass ratio from substrate Möbius topology," substrate framework paper 01 (2026).
[11] T. J. Hendrickson, "Saturation cap σ ≤ 1/2 forced as Möbius Z/2 fixed point," substrate framework paper 02 (2026).

## Appendix A — Reproducibility

```python
from src.stiff_medium.g_minus_2_substrate import GMinus2, ALPHA

g = GMinus2()

# The substrate-only contribution (the new physics):
delta_b3 = g.a_mu_substrate()
print(f"Δa_μ^B3 = {delta_b3:.4f} × 10⁻¹¹")    # 45.1179

# Decomposition:
print(g.report())

# Comparisons:
print(g.compare_to_fermilab(hvp_source="BMW"))   # +0.41σ
print(g.compare_to_fermilab(hvp_source="TI"))    # −4.16σ
print(g.compare_to_BMW_lattice())                # B3+TI vs BMW SM, −1.46σ
```

The single-line numerical check anyone can run:

```python
from math import pi
ALPHA = 7.2973525693e-3
xi_mu = 9.0 / 125.0     # = (n_R/K_pair) · K_rank⁻³
delta = (ALPHA / (2*pi)) * xi_mu * (ALPHA/pi)**2 * 1e11
print(f"Δa_μ^B3 = {delta:.2f} × 10⁻¹¹")    # 45.12 × 10⁻¹¹
```

## Appendix B — Why this is not "new physics"

The phrase "new physics" in (g−2)_μ literature usually means a hypothetical heavy particle (SUSY chargino/neutralino, Z′, leptoquark, dark photon) whose virtual exchange shifts a_μ. The substrate framework is not in this category. The mechanism is:

- The muon already exists as a substrate-strain pattern (paper 01).
- Drag γ already generates the muon's mass via cone-bouncing (paper 01).
- The same drag γ generates a one-loop correction to the muon's coupling to the electromagnetic field (this paper).

Nothing new is added to the framework that wasn't already there for the mass derivation. Δa_μ^B3 is a **prediction of the substrate ontology that was already required for m_μ**, not an add-on.

This is also why the prediction comes with no fitted coefficient: the substrate Lagrangian was fixed by paper 01's mass derivation, and once fixed, every other observable is downstream. Δa_μ^B3 is one of the downstream observables that happens to land in the experimentally interesting range. The substrate framework would predict the same +45×10⁻¹¹ if Fermilab had measured a_μ central at 116 592 030 (giving +1.5σ residual) or at 116 591 990 (giving +2.5σ residual) — the framework has nothing to adjust.

The Standard Model, by contrast, adjusts via HVP — the largest single SM contribution to a_μ — and HVP itself is currently swinging at the ±200×10⁻¹¹ level depending on whether one trusts CMD-3 or KLOE/BaBar. A 45×10⁻¹¹ topology-fixed substrate contribution at +0.41σ is in this sense a more *constrained* prediction than the SM's, which has BMW-vs-data-driven flexibility ~3× wider than the substrate's full claim.

---

*This paper is part of the substrate framework corpus. Companion papers derive the muon-to-electron mass ratio (paper 01), σ ≤ 1/2 saturation cap (paper 02), hierarchy problem (03), Σm_ν cosmology (04), GW150914 chirp mass (05), CKM/PMNS mixing (06), Stefan-Boltzmann (07), and atomic spectra (08). Open gaps are flagged transparently; the strongest standalone result on (g−2)_μ is presented here.*
