# Cold-Outreach Email Templates

Five templates, each <200 words, each focused on **one pain point**, **one substrate solution**, and **one next step**. Customize the bracketed fields, send, and follow up after 5 business days.

---

## Template 1 — R&D Director at a deeptech startup

**Subject:** Cutting your DFT screen from $200k/yr to $20k/yr

Hi [first name],

I saw [Company]'s recent [press release / paper / job posting] on [specific material/molecule/component] and noticed you're running a heavy [DFT / lattice / FEM] pre-screen on candidates.

I built a substrate physics engine that replaces that pre-screen step with closed-form Python expressions — same six inputs as the standard model + ΛCDM + GR, but computing in microseconds where DFT takes core-hours. It already reproduces the Cornell string tension exactly, the BCS gap ratio exactly, and the dark-energy density at 0.04%. The full code, 2,697 passing tests, and 11 papers are open-source and run on a laptop.

For a deeptech R&D budget, the typical pilot lands at $15k–$25k for a 4-week engagement: I take one of your candidate-screening problems, benchmark a substrate scoring function against your last 12 months of DFT output, and write the pilot-to-production transition plan.

Worth a 30-minute call this week? I'll tell you within those 30 minutes whether the engine fits your problem. If it doesn't, I'll say so.

Best,
[Your name]
tjhendrx@icloud.com

---

## Template 2 — University Materials Science PI

**Subject:** Closed-form pre-screen for your [cathode / catalyst / oxide] candidates

Dear Professor [last name],

I read your group's [recent paper / preprint] on [specific material system] and was struck by the candidate-space size you're working through.

I've built an open-source substrate physics engine that gives closed-form scoring for crystal cohesive energy, lattice stiffness, and band edges — same six substrate primitives that already reproduce the Cornell string tension, BCS gap ratio, and Stefan–Boltzmann constant. The point is not to replace your DFT runs, but to *pre-filter* your candidate list before the DFT queue.

The engine is MIT-licensed (241k LOC, 2,697 tests). Your students can clone and run any of the ~30 sector demos in ~30 minutes. I'm happy to do a one-hour Zoom seminar with your group, walking through how to map your candidate-screening pipeline onto the substrate scoring function — no charge, no commitment.

If it's useful afterward, we can discuss either co-authorship on a methods paper or a small-scale consulting engagement for your group's compute budget.

Best,
[Your name]
tjhendrx@icloud.com
GitHub: [your-handle]/substrate-theory

---

## Template 3 — Superconductor Research Lab Head

**Subject:** Hard upper bound on T_c for ambient-pressure candidates

Dr. [last name],

I noticed your lab's [recent paper / talk] proposing [specific high-T_c candidate]. I wanted to share one number that may be worth checking before the next furnace run.

The substrate physics framework I've been building yields a hard upper bound at ambient pressure:

> **T_c,max ≈ Λ_QCD / R ≈ 130 K**

— matched against the 31-year HgBaCa₂Cu₃O₈ record of 134 K at about 4%. The same framework reproduces the BCS gap ratio 2Δ/k_B T_c = 2π/e^γ = 3.528 exactly and matches 6/9 elemental superconductors within 10% of measured 2Δ/k_B T_c with zero per-material parameters.

The full implementation is open-source: see `bcs_gap_ratio_test.py`, `substrate_eliashberg.py`, `superconductivity_substrate.py`. Your students can run it on a laptop in under a minute.

Would you be open to a 30-minute conversation? I'd like to learn whether this bound is useful as a candidate-pruning filter for your group, and if it is, we can discuss either co-authorship or a short consulting engagement.

Best,
[Your name]
tjhendrx@icloud.com

---

## Template 4 — Pharma Chemistry Director

**Subject:** Pre-screen 1B compounds before docking — closed-form, in Python

Hi [first name],

If [Company]'s discovery pipeline is anything like the industry default, you're staring at a 10⁹+ chemical space and a docking pipeline that costs cents per compound. The math doesn't support running everything.

I've built a substrate physics engine with closed-form modules for molecular bond energies, solvation free energies, and pK_a — `molecular_bond_substrate.py`, `solvation_substrate.py`, `pka_substrate.py`. They evaluate in microseconds and run on a laptop. The substrate framework already reproduces dozens of physical-chemistry quantities (Stefan–Boltzmann to 0.0005%, BCS gap ratio exactly, dark energy to 0.04%) from six fundamental inputs.

The pilot offer is concrete: 4 weeks, $15k–$25k. I take a calibration set of your historical actives and inactives, build a substrate "physicality score" pre-filter, and benchmark its true-positive retention vs how aggressively it compresses the candidate list. Typical target: 100× compression with <5% true-active loss.

Worth a 30-minute call?

Best,
[Your name]
tjhendrx@icloud.com

---

## Template 5 — Structural Engineering Firm

**Subject:** Microsecond fatigue inference for embedded health monitors

Hi [first name],

I noticed [Firm] does structural health work on [bridges / wind-turbine blades / jet-engine fan blades / specific component]. The standard FEM-mesh refresh cycle puts a hard floor on how often you can re-evaluate fatigue state — usually scheduled inspections rather than live monitoring.

I've built a substrate physics engine with a fracture/fatigue module (`fracture_substrate_test.py`, `stress_loading.py`, `topological_defect_perturbations.py`) that evaluates crack-tip stress intensity in microseconds. The whole thing runs on a laptop today — and it ports cleanly to embedded C via ctypes or RustPython for in-flight or in-structure deployment.

A typical 4-week pilot ($15k–$25k) takes one of your existing inspection-data sets, benchmarks the substrate fatigue prediction against your incumbent FEM workflow, and produces a transition plan for embedded deployment. The longer integration ($100k–$200k) ships the embedded port itself.

Worth a 30-minute conversation? I'll be straightforward: if the substrate engine doesn't fit your specific application, I'll tell you in the first 10 minutes.

Best,
[Your name]
tjhendrx@icloud.com

---

## Sending discipline

- **Subject lines** are specific dollar amounts or specific predictions. No "quick question" or "interested in a chat." 
- **Open with their work**, not yours. One sentence proving you read their paper / press release / job posting.
- **One pain, one solution, one ask.** Resist the urge to mention everything the engine does. Pick the slice that matches their domain.
- **Include the floor offer** ($5k–$25k pilot) so the price is never a surprise on the call.
- **Make it easy to say no.** "I'll tell you in 30 minutes whether it fits" is a stronger CTA than "let me know if you're interested."
- **Follow up once after 5 business days.** Don't follow up more than twice. If they don't bite, they're not the buyer.
