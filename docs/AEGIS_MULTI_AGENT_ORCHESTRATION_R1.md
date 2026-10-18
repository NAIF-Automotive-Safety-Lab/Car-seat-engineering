# V7-AEGIS Multi-Agent Orchestration R1

## Purpose
Unify the engineering control plane around CAR-SEAT-ENGINEERING without collapsing authority boundaries or converting tool output into safety evidence.

## Authority matrix
| Actor | Authority |
|---|---|
| Peter | Architecture, commands, independent verification, final gate decision |
| JON_CHATGPT | Design development, engineering proposals, controlled modifications |
| Manus | Execution, test development, simulation execution, evidence generation |
| Claude | Independent adversarial review; no merge or baseline authority |
| Floot | Web control plane, orchestration surface, runtime gateway |
| GitHub | Canonical versioned source of truth |

## Hard prohibitions
1. R4.1 is immutable unless an explicitly authorized engineering-change record creates a separate controlled derivative.
2. No silent mutation of V7 baseline geometry or historical evidence packages.
3. Candidate supplier data never becomes confirmed evidence without external provenance.
4. Software execution never becomes physical-test evidence.
5. Simulation output never becomes validation evidence without the required input and verification chain.
6. No automatic merge of artifact-changing work.
7. Every executable result must retain input identity, source commit, environment, command, output/hash, and limitation state.
8. Missing evidence produces BLOCKED / NOT_PROVEN, never an invented PASS.

## Orchestration flow
`Peter command -> John design/proposal -> Manus execution -> GitHub evidence/versioning -> Claude adversarial review -> Peter gate`

Floot is the observable control surface for the chain and may dispatch/record actions, but it cannot grant evidence authority.

## Canonical anchors
- Repository: `NAIF-Automotive-Safety-Lab/Car-seat-engineering`
- Branch: `main`
- Current canonical commit at contract creation: `8c0c00a7dc2428d73d96b81f91b1e8826c592c26`
- Immutable R4.1 SHA-256: `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

## Integration state at R1 creation
- GitHub connector: authenticated and operational for controlled repository actions.
- Floot project: `V7-AEGIS — Engineering Test, Simulation & Evidence Intelligence Platform` exists and contains existing John/GitHub bridge surfaces.
- Anthropic/Claude: credential connection requested; until completed, Claude lane remains BLOCKED.
- Direct JON_CHATGPT API: not available through current Floot resources; the John bridge therefore remains a controlled message/ledger surface, not a hidden live API.
- Manus: treated as an external execution lane; no fabricated direct API connection is declared.

## Release gate
This contract does not authorize manufacturer release, CAE safety validation, physical validation, regulatory compliance, or geometry change. It only establishes the control-plane architecture and evidence boundaries.