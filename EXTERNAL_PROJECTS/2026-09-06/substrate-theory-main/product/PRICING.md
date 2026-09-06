# Pricing & Business Model

Five engagement tiers, ordered from lowest commitment to highest. Pick the one that matches your problem; you can always graduate up.

---

## 1. Open Source — Free

**MIT-licensed code on GitHub.** The full Substrate Physics Engine: 241,800 lines of Python, 2,697 passing tests, 125 visualizations, 11 companion papers.

| What you get | What you don't get |
|---|---|
| Clone the repo, run `pytest`, audit every Category-A claim | Personalized derivation help |
| Use any of the ~30 sector demos (atomic, hadron, BCS, cosmology, BBN, ...) in your own work | Custom extensions for your domain |
| Cite freely in publications and grant proposals | SLA, support contract, or response guarantee |
| File GitHub issues — best-effort response | Priority bug fixes |

**Who this is for:** academic researchers, students, hobbyists, anyone evaluating whether the engine fits their problem before buying time.

**Price:** $0. **Time to first result:** ~30 minutes from clone to first plot.

---

## 2. Consulting — $250–$500/hr

**Hourly engagement, typically $5k–$25k per pilot.**

This is the most common entry point. You bring a hard prediction problem; I scope it, run the substrate calculation, benchmark against your incumbent tool, and write up a closed report with a transition plan.

| Tier | Rate | Typical scope |
|------|-----:|---------------|
| **Discovery** ($250/hr) | $250/hr | Mapping your problem to substrate modules; "can this help?" answer in ~10 hours |
| **Build-out** ($400/hr) | $400/hr | Custom calculation, benchmark vs DFT/lattice/FEM/MCMC, written report |
| **Architect** ($500/hr) | $500/hr | Pilot-to-production handoff; codebase walkthroughs, integration plan |

**Typical 4-week pilot:** $5k–$25k depending on depth.
- $5k pilot: one well-scoped calculation, single benchmark, one-page report
- $15k pilot: full substrate scoreboard for 3–5 of your candidate quantities, written report with derivability tags
- $25k pilot: end-to-end integration prototype + transition plan + on-call hours during your team's evaluation

**Deliverables (every pilot):**
- Substrate scoreboard with A/B/C derivability tags for every quantity in scope
- Benchmark numbers vs your incumbent tool (speed, accuracy, cost)
- Reproducible Python scripts in your repo
- Written final report (PDF, ~5–15 pages)
- 1-hour wrap call

**Who this is for:** R&D teams with a specific bottleneck — "we burn $200k/year on this DFT screen; can we cut it?" — who want a defensible answer before committing further.

---

## 3. Custom Integration — $50k–$200k

**Fixed-scope engagement, 8–12 weeks, deliverable is a Python package you own.**

For teams that want a domain-specific `substrate-XX` extension shipped as production code: substrate-battery, substrate-polymer, substrate-fluiddrag, substrate-fracture, substrate-pharma, etc.

| Tier | Price | Scope |
|------|------:|-------|
| **Standard integration** | $50k–$100k | One domain extension; ~3 closed-form models; full test suite; documentation; 4 weeks of post-ship support |
| **Production integration** | $100k–$200k | Multi-module domain extension; integration with your existing pipeline (DFT/FEM/proprietary); CI/CD setup; runbook; 12 weeks of post-ship support |

**What "you own" means:**
- The extension repo lives in your private GitHub org
- You hold the source code under your standard contributor agreement
- I retain the right to upstream non-proprietary derivations into the open-source core
- No royalty, no per-seat fee, no SaaS dependency

**Who this is for:** R&D directors who've completed a successful pilot and want a hardened production tool.

---

## 4. Enterprise License — Custom

**Annual license for use in production products you ship.**

When the substrate engine becomes part of *your* product (a screening SaaS, a design tool, an embedded health-monitor), use shifts from internal R&D to commercial deployment. That's a different conversation.

- **Annual license:** custom, typically $50k–$500k/year depending on deployment scale
- **Includes:** quarterly check-in, priority bug fixes, written assurance of derivation provenance for any quantity you ship to your customers
- **Compliance:** I'll provide written attestations of which quantities are Category-A (derived) vs Category-C (anchor) for your regulatory or scientific-defensibility documentation

**Who this is for:** companies embedding substrate calculations in a commercial product.

---

## 5. Course & Training — $2,000 (self-paced)

**Self-paced course on the AI-assisted multi-agent methodology.**

The substrate engine was built in about one year by a single independent researcher using Claude Code, parallel sub-agent worktrees, and a strict A/B/C derivability discipline. That methodology is reproducible.

| Tier | Price | Scope |
|------|------:|-------|
| **Self-paced course** | $2,000 | Video + written modules: dispatch design, A/B/C tagging, derivability scoreboards, reproducibility hygiene; lifetime access |
| **Live cohort training** | $5,000/seat | Same content + 4 weekly live sessions + cohort Slack |
| **Bespoke team training** | $20k–$50k | On-site or remote; tailored to your domain; up to 15 engineers |

**Who this is for:** R&D engineering managers who want their internal teams to ship substrate-style breadth from a small headcount.

---

## 6. White-label SaaS — $1k–$10k/month

**Hosted pre-screening platform under your brand.**

For teams who want substrate pre-screening as a service rather than as code. I host the inference endpoint; you hit the API from your existing pipeline.

| Tier | Price | Scope |
|------|------:|-------|
| **Starter** | $1,000/mo | 100k API calls/mo; one domain (e.g. battery cathode scoring); shared infra |
| **Growth** | $3,000/mo | 1M API calls/mo; up to 3 domains; dedicated worker |
| **Scale** | $10,000/mo | Unlimited calls; custom domain; SLA; private deployment in your VPC available |

**Who this is for:** product teams that want substrate scoring as a managed service without taking on the integration work.

---

## Bundles & graduation paths

The pricing tiers compose:

- **Pilot → Integration:** finish a $15k pilot, then commit to a $75k integration. The pilot fee is credited against the integration price.
- **Integration → Enterprise license:** if your custom integration ships in a commercial product, the enterprise license fee replaces additional one-time integration billing.
- **Course bundled with engagements:** every $50k+ integration includes one self-paced course seat for your internal lead.

---

## What I won't do

- **Fixed-price open-ended scopes.** I'll quote a fixed price for a *specified* deliverable. Open-ended "build us X eventually" is hourly only.
- **Equity-only engagements.** I will accept partial equity in exchange for a *discount* on cash, but cash baseline is non-negotiable: this engine took a year of full-time work to build.
- **Sole-practitioner long-term retainers.** I'm one person. I cap commitments at ~60% capacity to protect the open-source core. After ~3 active engagements I refer the next one to a wait list.

---

## Next step

→ **Email:** `tjhendrx@icloud.com`
→ **Subject:** `Substrate engine — pricing for [your problem]`

Tell me one paragraph about the prediction problem. I'll come back within 48 hours with a quoted scope or a clean "this isn't a substrate problem, here's who you should call instead."
