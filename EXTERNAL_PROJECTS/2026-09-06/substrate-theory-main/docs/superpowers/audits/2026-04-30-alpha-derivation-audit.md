# Fine-Structure Constant Audit
**Date:** 2026-04-30
**Question:** How close is the current `alpha_derivation` / `alpha_two_loop`
work to the low-energy fine-structure constant?

## Target

Use the current NIST/CODATA 2022 recommended inverse fine-structure constant:

```text
alpha^-1 = 137.035999177(21)
```

Source: NIST CODATA value for the inverse fine-structure constant,
https://physics.nist.gov/cgi-bin/cuu/Value?eqalphinv=

The executable code uses `stiff_medium.physical_constants` as the single source
of truth for this target. `scripts/alpha_audit.py` fails fast if the derivation
and two-loop audit sections drift from that shared constant.

## Reproduction Commands

Run from the repository root:

```bash
PYTHONPATH=src python3 scripts/alpha_audit.py
PYTHONPATH=src pytest tests/test_alpha_audit.py -q
```

The old section-level script names are retained as wrappers:

```bash
PYTHONPATH=src python3 scripts/alpha_derivation_test.py
PYTHONPATH=src python3 scripts/alpha_two_loop_test.py
```

## Results

`alpha_audit.py derivation`:

- RG running now round-trips reproducibly. Running from `m_e` up to
  `Q_substrate = 27 GeV` gives `alpha^-1(Q_substrate) = 129.157263733`;
  running that value back to `m_e` returns `137.035999177`.
- The beta value that gives `alpha_bare = 1/137.035999177` directly is
  `beta^2 = 3.910354*pi`.
- The Higgs/W breather constraint gives `beta^2 = 4.547836*pi`.
- These differ by `0.637482*pi`.
- The Higgs/W beta is in the attractive/free Coleman regime
  (`g_Thirring <= 0`) and gives no positive `alpha_bare` in the current map.

`alpha_audit.py two-loop`:

- Main one-loop-scan beta used by the script:
  `beta^2 = 3.908957*pi`, `alpha^-1_2loop = 134.9082`.
  This is low by `2.1278` in inverse-alpha, or `1.553%`.
- Best fixed grid point printed by the two-loop scan:
  `beta^2 = 3.910000*pi`, `alpha^-1_2loop = 136.507203650`.
  This is low by `0.528795527` in inverse-alpha, or `0.385880739%`.
  Equivalently, `alpha = 0.007325620724`, which is high by
  `0.000028268159` relative to `1/137.035999177`.
- A continuous beta fit finds
  `beta^2 = 3.910339619191*pi`, `alpha^-1_2loop = 137.035999177000`.
  This is an exact calibration to the target, not a prediction.
- The independently motivated Higgs/W beta remains incompatible with positive
  alpha in the naive Coleman-bosonization route.

## Conclusion

The current branch does not derive `alpha = 1/137.035999177` from first
principles.

The best finite number from the fixed two-loop scan is close but still off by
`0.528795527` in inverse-alpha (`0.385880739%`), or high by
`0.000028268159` in alpha itself. The exact match appears only when beta is
continuously fitted to the target. The independent Higgs/W beta does not
produce a finite positive alpha in the current Coleman map.

The honest remaining task is the missing Möbius-bundle effective-action
calculation, not another report that treats beta fitting as a derivation.
