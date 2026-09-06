# substrate_physics — 1-page quickstart

A single Python file that lets you compute B3 substrate-framework
predictions from your shell or your own scripts.

## 30-second tour

```
$ cd "/Users/hendrixx./Desktop/Substrate Theory/cli"

$ python substrate_physics.py info                # show framework anchors
$ python substrate_physics.py list                # show every prediction

$ python substrate_physics.py predict alpha_em            # fine-structure
$ python substrate_physics.py predict m_mu_over_m_e       # lepton ratio
$ python substrate_physics.py predict bandgap silicon
$ python substrate_physics.py predict ie carbon
$ python substrate_physics.py predict tc_max              # 128.9 K ceiling
$ python substrate_physics.py predict hadron p
$ python substrate_physics.py predict lifetime muon
$ python substrate_physics.py predict fracture 50 600     # K_I, sigma_y
```

Add `--json` to any `predict` call to get machine-readable output.

## Python API in one block

```python
from substrate_physics import SubstratePhysics

sp = SubstratePhysics()

# any predictor returns a Prediction dataclass
p = sp.predict_alpha_em()
print(p.value, p.precision_estimate)         # -> 0.0072970624 0.0040%
print(p.to_dict())                           # full provenance dict
```

## Batch processing

Write `input.csv` with one prediction per line:

```csv
# command,arg1,arg2
alpha_em
m_mu_over_m_e
ie,carbon
bandgap,silicon
hadron,p
fracture,50,600
```

Run:

```
python substrate_physics.py batch input.csv output.csv
```

`output.csv` will contain one row per prediction with value, unit,
precision, category and source-module provenance.

## Three things to know

1. **Zero parameters.** Every prediction uses the framework's 12 audited
   inventory integers (n_M=268, K_pair=2, K_rank=5, ...) and the
   Λ_QCD = 200 MeV anchor. No per-prediction tuning.
2. **Derivability tagged.** Each `Prediction.category` is `A` (substrate-
   derivable today), `B` (derivable in principle, calculation pending), or
   `C` (irreducible empirical anchor).
3. **Speed.** Each prediction is closed-form; expect well under 10 ms per
   call.

## Where the formulae live

The `source_module` field on each `Prediction` points to the file in
`src/stiff_medium/` that derives the result.
