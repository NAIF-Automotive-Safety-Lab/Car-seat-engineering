# substrate_physics — quickstart

A single-file CLI + Python library that wraps the most useful predictions of
the **B3 / Stiff-Medium substrate framework** into one tool.

Install: just drop `substrate_physics.py` somewhere on your `$PATH` (or this
folder) and run it. The only hard dependency is the Python standard library;
`numpy` is optional and `colorama` is used for prettier terminal output if it
happens to be installed.

```
python substrate_physics.py info
python substrate_physics.py list
python substrate_physics.py predict alpha_em
python substrate_physics.py predict ie carbon
python substrate_physics.py predict bandgap silicon
python substrate_physics.py predict tc_max
```

## CLI commands

| Command | Description |
|---|---|
| `predict TOPIC [args...]` | Run one prediction. Emits coloured text by default; `--json` for machine-readable output. |
| `list` | List all available predictions. |
| `info` | Show framework info: integers, anchors, derivability classes. |
| `batch IN.csv OUT.csv` | Run many predictions. Each input row: `command,arg1,arg2,...`. |

## Library use

```python
from substrate_physics import SubstratePhysics

sp = SubstratePhysics()

p = sp.predict_lepton_mass_ratio()
print(p.value, p.unit, p.precision_estimate)

print(sp.predict_atomic_ie("C").value)
print(sp.predict_bandgap("silicon").to_dict())
print(sp.predict_alpha_em().value)
print(sp.predict_hadron_mass("p"))
```

Each predictor returns a `Prediction` dataclass with fields:

- `name`, `value`, `unit`
- `precision_estimate`  (e.g. `"0.004%"`)
- `category`            (`A` substrate-derivable / `B` derivable in principle / `C` empirical anchor)
- `source_module`       (which substrate file holds the derivation)
- `measured_value`      (PDG / CODATA / NIST value if known)
- `source_reference`    (where the measured value came from)
- `notes`, `extras`

## What is exposed

| Method | Returns |
|---|---|
| `predict_lepton_mass_ratio()` | m_μ / m_e from inventory integers |
| `predict_tau_over_electron_ratio()` | m_τ / m_e from lepton tower |
| `predict_atomic_ie(element)` | First ionisation energy (H..Ar) |
| `predict_bandgap(material)` | Semiconductor bandgap |
| `predict_madelung(crystal)` | Madelung constant |
| `predict_fracture_zone(K_I, sigma_y)` | Irwin plastic-zone radius (mm) |
| `predict_hadron_mass(name)` | Hadron mass from face-spin v4 |
| `predict_lifetime(particle)` | Particle lifetime (s) |
| `predict_debye_temperature(material)` | Θ_D in K |
| `predict_bcs_gap_ratio(material)` | 2 Δ / k_B T_c |
| `predict_tc_max()` | Substrate T_c ceiling (128.9 K) |
| `predict_dark_matter_mass()` | Cube-cell DM (27.5 GeV) |
| `predict_neutrino_sum()` | Σ m_ν (60.5 meV) |
| `predict_alpha_em()` | Fine-structure constant (closed form) |
| `predict_hierarchy()` | exp(4 π² − 1) gauge-hierarchy ratio |
| `predict_string_tension()` | Cornell σ from inventory |
| `predict_grueneisen(material)` | γ_G |

## Performance

Each prediction runs in well under 10 ms — the substrate's analytic
closed-form character means there is no iterative solver in the hot path.

## Provenance

All formulae trace to one of the modules in `src/stiff_medium/`. The
canonical 12 inventory integers (N_BAM=6, K_pair=2, K_rank=5, n_R=18,
n_M=268, n_A=15, F=2, R=3, V13=13) and Λ_QCD=200 MeV anchor are exposed
through `SubstratePhysics().constants`.

See `examples.py` for 10 worked examples and `QUICKSTART.md` for a 1-page
intro.
