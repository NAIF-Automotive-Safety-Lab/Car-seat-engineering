# 10 Concrete Client Scenarios

Each scenario follows the same structure: **the pain**, **the substrate fit**, **the deliverable**, **the indicative price**.

These are scoped at pilot size ($5k–$25k) unless noted. Larger integrations ($50k–$200k) and enterprise licenses are quoted separately.

---

## 1. Battery manufacturer screening 100k cathodes

**The pain.** Industrial DFT screens of cathode candidates burn $200k–$1M per year of cluster time, and >95% of candidates fail downstream stability or capacity gates.

**The substrate fit.** The substrate cohesive-energy and lattice-stiffness modules (`substrate_cohesive_energy.py`, `substrate_elasticity.py`, `crystal_substrate.py`) score crystal candidates in microseconds against the same six substrate primitives that already match Cornell string tension and BCS gap ratio. Use it as the pre-DFT filter.

**The deliverable.** A scoring function `f(composition, lattice) → (stability, capacity_proxy)` benchmarked against your last 12 months of DFT outputs. We define a Pareto-front filter that passes the top 1–5% of candidates to your existing DFT pipeline.

**Indicative price.** $15k–$25k pilot. Expected outcome: 20–50× cluster cost reduction per round of screening.

---

## 2. Pharma pre-screening 1B compounds

**The pain.** Modern combinatorial chemistry generates 10⁹+ candidate molecules. Even ultra-fast docking pipelines (Schrödinger Glide, AutoDock Vina) cost cents per pose, so a billion-compound screen is prohibitive.

**The substrate fit.** The substrate molecular-bond and solvation modules (`molecular_bond_substrate.py`, `solvation_substrate.py`, `pka_substrate.py`) give closed-form bond energies and solvation free energies. Use them to filter the 10⁹ space down to ~10⁵ chemically-realistic leads before docking.

**The deliverable.** A Python pre-filter module that takes SMILES strings in, emits a substrate-derived "physicality score" out, calibrated against a reference set of your historical actives and inactives. Achieves >100× compression with <5% loss of true actives in the calibration validation.

**Indicative price.** $25k pilot or $75k–$150k production integration. Expected outcome: $0.001 per compound vs $0.10+ for docking — 100× pre-screen savings.

---

## 3. Aerospace real-time fatigue monitoring

**The pain.** Structural-health monitoring on rotating or vibrating airframe components (rotor blades, wing roots, jet-engine fan blades) can't run a full FEM mesh refresh in flight. Operators rely on coarse vibration signatures and conservative inspection intervals.

**The substrate fit.** The substrate fracture model (`fracture_substrate_test.py`, `stress_loading.py`, `topological_defect_perturbations.py`) evaluates crack-tip stress intensity in microseconds — fast enough for embedded monitoring at kilohertz sampling rates without a GPU.

**The deliverable.** A C-callable substrate fracture inference function (port from Python via ctypes or RustPython), benchmarked on your last 3 years of inspection data. Produces a real-time "remaining safe-cycle estimate" with confidence interval.

**Indicative price.** $25k pilot, $100k–$200k embedded port. Expected outcome: extend inspection intervals by 20–50% with no safety regression.

---

## 4. Semiconductor design for novel materials

**The pain.** A typical semiconductor design loop tries hundreds of compositions. Full DFT band structure for each is hours-to-days. The wait kills iteration speed.

**The substrate fit.** The substrate semiconductor module (`semiconductor_substrate.py`, `crystal_substrate.py`, `phonon_dispersion.py`) gives closed-form band-edge energies, electron/hole effective masses, and phonon dispersions in seconds.

**The deliverable.** A design-of-experiments tool that produces a substrate-derived band-structure summary for any candidate composition in <1 second. Calibrated against your reference materials, it captures band gap to ±0.1 eV — sufficient for "should we even DFT this?" go/no-go.

**Indicative price.** $15k pilot. Expected outcome: 10× faster early-funnel iteration, 50%+ reduction in wasted DFT runs.

---

## 5. Cosmology forecasting for grant proposals

**The pain.** A competitive observational-cosmology grant needs a defensible forecast: what's your expected H_0, σ_8, Σm_ν, dark-energy equation-of-state under various survey designs? Standard MCMC chains run for days and need ~14 nuisance parameters.

**The substrate fit.** The substrate cosmology stack (`cosmology_simulator.py`, `hubble_paired.py`, `sigma_mnu_falsifier.py`, `bao_sound_horizon.py`, `cmb_paired.py`) gives substrate-derived predictions for H_0, σ_8, and Σm_ν directly from substrate primitives. Predictions hit Planck precision on ρ_Λ (0.04%) and live within both Planck and SH0ES on H_0.

**The deliverable.** A grant-ready forecasting script that accepts your survey design (effective area, redshift coverage, sensitivity) and emits substrate-derived parameter forecasts plus a tension diagnostic vs ΛCDM. Includes companion paper draft section ready for the proposal.

**Indicative price.** $5k–$15k. Expected outcome: shipped grant proposal in 1–2 weeks.

---

## 6. Superconductor research lab — high-T_c candidate filter

**The pain.** Synthesizing candidate high-temperature superconductors is expensive and slow. Wishful-thinking candidates (e.g. predicted T_c far above the ambient-pressure record) waste furnace time.

**The substrate fit.** The substrate superconductivity stack (`superconductivity_substrate.py`, `bcs_gap_ratio_test.py`, `substrate_eliashberg.py`, `substrate_multiband.py`) reproduces the BCS gap ratio 2π/e^γ exactly and gives a hard upper bound T_c,max ≈ Λ_QCD/R ≈ 130 K at ambient pressure — within ~4% of the 31-year HgBaCa₂Cu₃O₈ record (134 K).

**The deliverable.** A T_c-feasibility filter that scores any proposed material against the substrate upper bound and the substrate-Eliashberg gap-ratio prediction. Identifies which candidates have any chance of beating the ambient-pressure record vs which are physically ruled out.

**Indicative price.** $5k–$15k pilot. Expected outcome: 20–50% reduction in synthesis attempts on infeasible candidates.

---

## 7. University physics teaching platform

**The pain.** Undergraduate physics courses jump from one isolated formula to the next — hydrogen energy levels, Stefan–Boltzmann, BCS gap, Cornell string tension — with no unifying picture.

**The substrate fit.** The substrate engine *is* the unifying picture: one Lagrangian, six inputs, ~30 sectors. The visuals folder ships 125 figures and animations spanning atomic, hadron, BCS, cosmology, biology.

**The deliverable.** A teaching-platform license + curriculum-integration consulting: notebooks that walk students through one substrate primitive at a time, showing how the same Lagrangian produces hydrogen energy, the cosmic microwave background, and the BCS gap. White-label option available.

**Indicative price.** Open-source self-serve free; $20k–$50k integration with onboarding and in-class support.

---

## 8. HPC cost reduction (1000× cheaper than DFT)

**The pain.** Your group's HPC bill is $500k–$5M/year, dominated by DFT. Most runs are early-funnel screens, not paper-quality calculations.

**The substrate fit.** Replace 90%+ of early-funnel DFT screens with substrate scoring. The substrate engine gives you a closed-form proxy at <1 ms per candidate vs ~10 minutes for DFT.

**The deliverable.** A pre-DFT filter that ingests your candidate queue, scores everything in substrate, and routes only the top fraction to DFT. Audit against your last 6 months of DFT outputs to set the cutoff threshold conservatively.

**Indicative price.** $25k pilot, $100k–$200k production integration. Expected outcome: 5–20× reduction in cluster spend without measurably worse hit rate.

---

## 9. Open-source community — GitHub stars / pull requests

**The pain.** Independent researchers and small academic groups want to benchmark new ideas against substrate-derived baselines without the friction of a full DFT install.

**The substrate fit.** Free, MIT-licensed, runs on a laptop. The 2,697-test suite is the trust signal: any user can audit the framework end-to-end.

**The deliverable.** Public GitHub repo with full source, papers, visualizations. Issues triaged best-effort. Pull requests welcome for new sector modules under the same A/B/C derivability discipline.

**Indicative price.** $0. Value: community amplification, citation pipeline, contributor pool that I can later draw from for paid engagements.

---

## 10. AI methodology consulting for R&D departments

**The pain.** R&D engineering managers see what one independent researcher built with AI-assisted multi-agent dispatch in a year, and want their team of 8 to deliver at the same breadth — but their team is using ChatGPT one tab at a time.

**The substrate fit.** The substrate engine *is* the case study. The methodology — Claude Code + parallel sub-agent worktrees + strict A/B/C derivability tagging + reproducibility hygiene — is independently sellable, regardless of whether your domain is physics.

**The deliverable.** A 6–8 week training engagement: workshop on dispatch design, hands-on rebuild of a small substrate-style scoreboard in your domain, written playbook for your team. Plus optional self-paced course seats for new hires going forward.

**Indicative price.** $20k–$50k bespoke training; $2k self-paced course; $5k/seat live cohort. Expected outcome: 2–5× output per engineer on AI-amenable R&D work.

---

## How to start

For each scenario above, the entry point is the same:

→ **Email:** `tjhendrx@icloud.com`
→ **Subject:** `Substrate engine — [scenario name]`
→ **Body:** one paragraph describing your specific situation

I'll respond within 48 hours with either a scoped pilot quote or an honest "this isn't a substrate problem, here's who else to ask."
