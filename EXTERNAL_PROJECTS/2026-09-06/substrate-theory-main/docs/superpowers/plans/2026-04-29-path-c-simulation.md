# Path C: Lattice-Gas Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a 2D simulation that tests whether the Stiff-Medium Confinement Theory's displacement-only neutrino dynamics produce stable bound states (electron formation) under direct execution of the rules — no fitting, no fudge factors. Per the spec at `docs/superpowers/specs/2026-04-29-stiff-medium-theory-design.md` §2 (Methodology) and §15 (Roadmap).

**Architecture:** Python 3.11+ with numpy for arrays, matplotlib for animation, pytest for tests. Particles (neutrinos) are dataclasses with position and velocity vectors constrained to 45° in the xy plane. The simulation loop propagates particles, detects coordinate-overlap conflicts, and resolves them by displacing positions while leaving velocities unchanged (the load-bearing rule from spec §5). A bound-state detector flags persistence near a shared center as a candidate "electron." The experiment runs two neutrinos on collision course and reports what happens — bound state forms, or doesn't.

**Tech Stack:** Python 3.11+, numpy ≥ 1.24, matplotlib ≥ 3.7, pytest ≥ 7. No other dependencies.

**Out of scope for v1 (Path C):** 3D, bi-pyramid / nucleon formation, atomic structure, lepton stress-loading, gravitational coupling. These are deferred until v1 either (a) produces stable electrons or (b) tells us the rules need revision.

---

## File structure

```
src/stiff_medium/
├── __init__.py
├── neutrino.py        # Neutrino dataclass with 45° + speed-c validation
├── dynamics.py        # propagate(), detect_overlap(), displace(), step()
├── detector.py        # BoundStateTracker — flags persistent proximity
└── visualize.py       # animate(history) — matplotlib FuncAnimation

tests/
├── __init__.py
├── test_neutrino.py
├── test_dynamics.py
├── test_detector.py
└── test_visualize.py

scripts/
└── electron_formation.py   # Main experiment: 2 neutrinos on collision course

pyproject.toml
README.md
```

Each file has one responsibility. `dynamics.py` is the load-bearing rule logic; everything else is plumbing or analysis. Tests parallel source files exactly.

---

## Tasks

### Task 1: Project scaffolding

**Files:**
- Create: `pyproject.toml`
- Create: `README.md`
- Create: `src/stiff_medium/__init__.py`
- Create: `tests/__init__.py`
- Create: `.gitignore`

- [ ] **Step 1: Create `pyproject.toml`**

```toml
[project]
name = "stiff_medium"
version = "0.1.0"
description = "Path C lattice-gas simulation for Stiff-Medium Confinement Theory"
requires-python = ">=3.11"
dependencies = [
    "numpy>=1.24",
    "matplotlib>=3.7",
]

[project.optional-dependencies]
dev = ["pytest>=7"]

[build-system]
requires = ["setuptools>=64", "wheel"]
build-backend = "setuptools.build_meta"

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
```

- [ ] **Step 2: Create `README.md`**

```markdown
# Stiff-Medium Path C Simulation

2D lattice-gas simulation testing whether the displacement-only neutrino dynamics described in `docs/superpowers/specs/2026-04-29-stiff-medium-theory-design.md` produce stable electron formation.

## Setup
```
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Run experiment
```
python scripts/electron_formation.py
```

## Test
```
pytest
```

## What this tests

Two neutrinos on collision course at 45° angles. If the rules produce a bound orbital state, that's a positive result for the theory. If they don't, the rules need revision.
```

- [ ] **Step 3: Create `.gitignore`**

```
.venv/
__pycache__/
*.egg-info/
.pytest_cache/
*.pyc
.DS_Store
```

- [ ] **Step 4: Create `src/stiff_medium/__init__.py`**

```python
"""Stiff-Medium Confinement Theory — Path C lattice-gas simulation."""
```

- [ ] **Step 5: Create `tests/__init__.py`**

Empty file:

```python
```

- [ ] **Step 6: Install and smoke test**

Run:
```
python -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"
```
Expected: install completes without error.

Run:
```
pytest --collect-only
```
Expected: "no tests ran" (no tests yet) without import errors.

- [ ] **Step 7: Commit**

If git is initialized:
```
git add pyproject.toml README.md .gitignore src/ tests/
git commit -m "chore: scaffold Path C simulation project"
```

If git is NOT initialized, run `git init` first then the above.

---

### Task 2: Neutrino dataclass with validation

**Files:**
- Create: `src/stiff_medium/neutrino.py`
- Create: `tests/test_neutrino.py`

The neutrino is a 2D particle: position (x, y), velocity vector (vx, vy). Velocity must be at 45° to one axis (i.e., |vx| = |vy|) and have magnitude c (we use c = 1 in natural units). Per spec §5.

- [ ] **Step 1: Write the failing test for valid construction**

`tests/test_neutrino.py`:

```python
import numpy as np
import pytest
from stiff_medium.neutrino import Neutrino, C


def test_construct_valid_neutrino():
    n = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([C / np.sqrt(2), C / np.sqrt(2)]),
    )
    assert np.allclose(n.position, [0.0, 0.0])
    assert np.allclose(n.velocity, [C / np.sqrt(2), C / np.sqrt(2)])
```

- [ ] **Step 2: Run test to confirm it fails**

Run:
```
pytest tests/test_neutrino.py::test_construct_valid_neutrino -v
```
Expected: FAIL with "ModuleNotFoundError: No module named 'stiff_medium.neutrino'".

- [ ] **Step 3: Implement minimal Neutrino**

`src/stiff_medium/neutrino.py`:

```python
"""Neutrino: 1D-vector particle in the 2D simulation. See spec §5."""

from dataclasses import dataclass
import numpy as np

C: float = 1.0  # natural units; c is the medium's wave speed


@dataclass
class Neutrino:
    """A neutrino is a position + velocity vector at 45° to an axis, speed c.

    The velocity vector is never reoriented during simulation (spec §5).
    Only the position is displaced when overlap conflicts arise.
    """

    position: np.ndarray  # shape (2,), float
    velocity: np.ndarray  # shape (2,), float, |vx|=|vy|, magnitude C

    def __post_init__(self) -> None:
        if self.position.shape != (2,):
            raise ValueError(f"position must be shape (2,), got {self.position.shape}")
        if self.velocity.shape != (2,):
            raise ValueError(f"velocity must be shape (2,), got {self.velocity.shape}")
        if not np.isclose(abs(self.velocity[0]), abs(self.velocity[1])):
            raise ValueError(
                f"velocity must be 45° to an axis (|vx|=|vy|), got {self.velocity}"
            )
        speed = float(np.linalg.norm(self.velocity))
        if not np.isclose(speed, C):
            raise ValueError(f"velocity magnitude must be C={C}, got {speed}")
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```
pytest tests/test_neutrino.py::test_construct_valid_neutrino -v
```
Expected: PASS.

- [ ] **Step 5: Add validation tests**

Append to `tests/test_neutrino.py`:

```python
def test_reject_non_45_velocity():
    with pytest.raises(ValueError, match="45"):
        Neutrino(
            position=np.array([0.0, 0.0]),
            velocity=np.array([1.0, 0.0]),  # along x-axis, not 45°
        )


def test_reject_wrong_speed():
    with pytest.raises(ValueError, match="magnitude"):
        Neutrino(
            position=np.array([0.0, 0.0]),
            velocity=np.array([1.0, 1.0]),  # 45° but magnitude √2 ≠ C
        )


def test_accept_all_four_45_directions():
    s = C / np.sqrt(2)
    for vx, vy in [(s, s), (s, -s), (-s, s), (-s, -s)]:
        n = Neutrino(
            position=np.array([0.0, 0.0]),
            velocity=np.array([vx, vy]),
        )
        assert np.isclose(np.linalg.norm(n.velocity), C)
```

- [ ] **Step 6: Run all neutrino tests**

Run:
```
pytest tests/test_neutrino.py -v
```
Expected: 4 tests PASS.

- [ ] **Step 7: Commit**

```
git add src/stiff_medium/neutrino.py tests/test_neutrino.py
git commit -m "feat(neutrino): add Neutrino dataclass with 45°+speed-c validation"
```

---

### Task 3: Free propagation

**Files:**
- Modify: `src/stiff_medium/dynamics.py` (create new)
- Create: `tests/test_dynamics.py`

A neutrino in free motion (no overlap) advances its position by `velocity * dt` each step. Velocity is unchanged. Per spec §5: "Vectors never reorient."

- [ ] **Step 1: Write the failing test**

`tests/test_dynamics.py`:

```python
import numpy as np
from stiff_medium.neutrino import Neutrino, C
from stiff_medium.dynamics import propagate


def test_propagate_advances_position_by_velocity_times_dt():
    s = C / np.sqrt(2)
    n = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    dt = 0.1
    moved = propagate(n, dt)
    assert np.allclose(moved.position, [s * dt, s * dt])


def test_propagate_preserves_velocity():
    s = C / np.sqrt(2)
    n = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, -s]),
    )
    moved = propagate(n, 1.0)
    assert np.allclose(moved.velocity, [s, -s])
```

- [ ] **Step 2: Run tests to confirm they fail**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: ImportError on `propagate`.

- [ ] **Step 3: Implement `propagate`**

`src/stiff_medium/dynamics.py`:

```python
"""Dynamics rules. See spec §5 for the load-bearing displacement rule."""

import numpy as np
from stiff_medium.neutrino import Neutrino


def propagate(n: Neutrino, dt: float) -> Neutrino:
    """Advance position by velocity*dt. Velocity unchanged (spec §5)."""
    return Neutrino(
        position=n.position + n.velocity * dt,
        velocity=n.velocity.copy(),
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: 2 tests PASS.

- [ ] **Step 5: Commit**

```
git add src/stiff_medium/dynamics.py tests/test_dynamics.py
git commit -m "feat(dynamics): add free propagation"
```

---

### Task 4: Overlap detection

**Files:**
- Modify: `src/stiff_medium/dynamics.py`
- Modify: `tests/test_dynamics.py`

Two neutrinos "would overlap a coordinate" when their positions are within a threshold radius `r_overlap`. This is the trigger for the displacement rule.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_dynamics.py`:

```python
from stiff_medium.dynamics import detect_overlap


def test_detect_overlap_when_close():
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([C / np.sqrt(2), C / np.sqrt(2)]),
    )
    n2 = Neutrino(
        position=np.array([0.05, 0.05]),
        velocity=np.array([-C / np.sqrt(2), -C / np.sqrt(2)]),
    )
    assert detect_overlap(n1, n2, r_overlap=0.1) is True


def test_no_overlap_when_far():
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([C / np.sqrt(2), C / np.sqrt(2)]),
    )
    n2 = Neutrino(
        position=np.array([10.0, 10.0]),
        velocity=np.array([-C / np.sqrt(2), -C / np.sqrt(2)]),
    )
    assert detect_overlap(n1, n2, r_overlap=0.1) is False
```

- [ ] **Step 2: Run tests to confirm they fail**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: ImportError on `detect_overlap`.

- [ ] **Step 3: Implement `detect_overlap`**

Append to `src/stiff_medium/dynamics.py`:

```python
def detect_overlap(a: Neutrino, b: Neutrino, r_overlap: float) -> bool:
    """Return True if two neutrinos are within r_overlap of each other.

    This is the trigger for the displacement rule (spec §5).
    """
    distance = float(np.linalg.norm(a.position - b.position))
    return distance < r_overlap
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```
git add src/stiff_medium/dynamics.py tests/test_dynamics.py
git commit -m "feat(dynamics): add overlap detection"
```

---

### Task 5: Displacement rule (the load-bearing one)

**Files:**
- Modify: `src/stiff_medium/dynamics.py`
- Modify: `tests/test_dynamics.py`

When two neutrinos overlap, displace them apart along the line connecting their positions. **Velocities are NOT changed** (spec §5: "Vectors never reorient"). This is the central rule the simulation tests.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_dynamics.py`:

```python
from stiff_medium.dynamics import displace


def test_displace_pushes_neutrinos_apart():
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([0.05, 0.0]),
        velocity=np.array([-s, s]),
    )
    moved1, moved2 = displace(n1, n2, push=0.1)
    new_dist = float(np.linalg.norm(moved1.position - moved2.position))
    old_dist = float(np.linalg.norm(n1.position - n2.position))
    assert new_dist > old_dist


def test_displace_preserves_both_velocities():
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([0.05, 0.0]),
        velocity=np.array([-s, s]),
    )
    moved1, moved2 = displace(n1, n2, push=0.1)
    assert np.allclose(moved1.velocity, n1.velocity)
    assert np.allclose(moved2.velocity, n2.velocity)


def test_displace_handles_coincident_positions_safely():
    """If positions are identical, displacement direction is degenerate;
    must not crash, must not produce NaN."""
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([1.0, 1.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([1.0, 1.0]),
        velocity=np.array([-s, -s]),
    )
    moved1, moved2 = displace(n1, n2, push=0.1)
    assert np.all(np.isfinite(moved1.position))
    assert np.all(np.isfinite(moved2.position))
    assert not np.allclose(moved1.position, moved2.position)
```

- [ ] **Step 2: Run tests to confirm they fail**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: ImportError on `displace`.

- [ ] **Step 3: Implement `displace`**

Append to `src/stiff_medium/dynamics.py`:

```python
def displace(
    a: Neutrino, b: Neutrino, push: float
) -> tuple[Neutrino, Neutrino]:
    """Push two overlapping neutrinos apart along the line connecting them.

    Velocities are NOT changed (spec §5). If positions coincide exactly,
    use the velocity-difference vector as a fallback; if that's also zero,
    fall back to (1, 0).
    """
    diff = b.position - a.position
    norm = float(np.linalg.norm(diff))

    if norm < 1e-12:
        # Coincident: use velocity difference, then a fixed fallback.
        diff = b.velocity - a.velocity
        norm = float(np.linalg.norm(diff))
    if norm < 1e-12:
        diff = np.array([1.0, 0.0])
        norm = 1.0

    unit = diff / norm
    shift = unit * (push / 2.0)

    moved_a = Neutrino(position=a.position - shift, velocity=a.velocity.copy())
    moved_b = Neutrino(position=b.position + shift, velocity=b.velocity.copy())
    return moved_a, moved_b
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```
git add src/stiff_medium/dynamics.py tests/test_dynamics.py
git commit -m "feat(dynamics): add displacement rule (vectors preserved)"
```

---

### Task 6: Simulation step (compose propagate + detect + displace)

**Files:**
- Modify: `src/stiff_medium/dynamics.py`
- Modify: `tests/test_dynamics.py`

A `step` function takes a list of neutrinos and a dt, propagates each, then resolves any pairwise overlaps. Returns the new state.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_dynamics.py`:

```python
from stiff_medium.dynamics import step


def test_step_propagates_isolated_neutrinos():
    s = C / np.sqrt(2)
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([10.0, 10.0]),
        velocity=np.array([-s, -s]),
    )
    new_state = step([n1, n2], dt=1.0, r_overlap=0.1, push=0.1)
    # Far apart, so they just propagate freely.
    assert np.allclose(new_state[0].position, [s, s])
    assert np.allclose(new_state[1].position, [10.0 - s, 10.0 - s])


def test_step_displaces_overlapping_neutrinos():
    s = C / np.sqrt(2)
    # Place them so propagation will bring them within r_overlap.
    n1 = Neutrino(
        position=np.array([0.0, 0.0]),
        velocity=np.array([s, s]),
    )
    n2 = Neutrino(
        position=np.array([s * 1.0, s * 1.0]),  # right where n1 will be after dt=1
        velocity=np.array([-s, -s]),
    )
    new_state = step([n1, n2], dt=1.0, r_overlap=0.5, push=0.2)
    # After step, they should be > r_overlap apart and have unchanged velocities.
    dist = float(np.linalg.norm(new_state[0].position - new_state[1].position))
    assert np.allclose(new_state[0].velocity, n1.velocity)
    assert np.allclose(new_state[1].velocity, n2.velocity)
    # Displacement happened at least.
    assert dist > 0.0
```

- [ ] **Step 2: Run tests to confirm they fail**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: ImportError on `step`.

- [ ] **Step 3: Implement `step`**

Append to `src/stiff_medium/dynamics.py`:

```python
def step(
    neutrinos: list[Neutrino],
    dt: float,
    r_overlap: float,
    push: float,
) -> list[Neutrino]:
    """One simulation step: propagate all, then resolve pairwise overlaps.

    Overlap resolution iterates pairwise (O(n²) — fine for small n; spec v1
    only requires 2-particle experiments).
    """
    moved = [propagate(n, dt) for n in neutrinos]

    # Resolve pairwise overlaps. Single pass is enough for n=2; for larger
    # n, repeat until no overlaps remain (left for v2).
    for i in range(len(moved)):
        for j in range(i + 1, len(moved)):
            if detect_overlap(moved[i], moved[j], r_overlap):
                moved[i], moved[j] = displace(moved[i], moved[j], push)

    return moved
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```
pytest tests/test_dynamics.py -v
```
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```
git add src/stiff_medium/dynamics.py tests/test_dynamics.py
git commit -m "feat(dynamics): add step function composing propagate+detect+displace"
```

---

### Task 7: Bound-state tracker

**Files:**
- Create: `src/stiff_medium/detector.py`
- Create: `tests/test_detector.py`

A bound state is when two neutrinos remain within a threshold distance for at least N consecutive steps. This is our operational definition of "electron formation" for v1 (spec §6: stable orbital pattern).

- [ ] **Step 1: Write the failing tests**

`tests/test_detector.py`:

```python
import numpy as np
from stiff_medium.neutrino import Neutrino, C
from stiff_medium.detector import BoundStateTracker


def _n(x, y, vx, vy):
    return Neutrino(
        position=np.array([x, y], dtype=float),
        velocity=np.array([vx, vy], dtype=float),
    )


def test_tracker_starts_unbound():
    s = C / np.sqrt(2)
    tracker = BoundStateTracker(r_bound=0.5, persistence=10)
    a, b = _n(0, 0, s, s), _n(10, 10, -s, -s)
    assert tracker.update([a, b]) is False


def test_tracker_flags_bound_after_persistence_steps():
    s = C / np.sqrt(2)
    tracker = BoundStateTracker(r_bound=0.5, persistence=3)
    a, b = _n(0, 0, s, s), _n(0.1, 0, -s, -s)
    assert tracker.update([a, b]) is False  # step 1 close
    assert tracker.update([a, b]) is False  # step 2 close
    assert tracker.update([a, b]) is True   # step 3 close — flagged


def test_tracker_resets_on_separation():
    s = C / np.sqrt(2)
    tracker = BoundStateTracker(r_bound=0.5, persistence=3)
    a, b = _n(0, 0, s, s), _n(0.1, 0, -s, -s)
    far_a, far_b = _n(0, 0, s, s), _n(10, 10, -s, -s)
    tracker.update([a, b])
    tracker.update([a, b])
    tracker.update([far_a, far_b])  # separation resets counter
    # Need 3 more close steps to flag.
    assert tracker.update([a, b]) is False
    assert tracker.update([a, b]) is False
    assert tracker.update([a, b]) is True
```

- [ ] **Step 2: Run tests to confirm they fail**

Run:
```
pytest tests/test_detector.py -v
```
Expected: ImportError on `BoundStateTracker`.

- [ ] **Step 3: Implement `BoundStateTracker`**

`src/stiff_medium/detector.py`:

```python
"""Bound-state detection. Spec §6: a stable orbital pattern is two
neutrinos persisting near each other long enough that they're not just
passing through. Operationally: within r_bound for `persistence` steps."""

import numpy as np
from stiff_medium.neutrino import Neutrino


class BoundStateTracker:
    """Flags when two neutrinos have been within r_bound for N consecutive steps."""

    def __init__(self, r_bound: float, persistence: int) -> None:
        if r_bound <= 0:
            raise ValueError("r_bound must be positive")
        if persistence < 1:
            raise ValueError("persistence must be >= 1")
        self.r_bound = r_bound
        self.persistence = persistence
        self._consecutive = 0

    def update(self, neutrinos: list[Neutrino]) -> bool:
        """Update tracker with current state. Returns True if a bound
        state has been detected (persistence threshold reached this step).

        v1: only handles the 2-neutrino case; ignores n != 2.
        """
        if len(neutrinos) != 2:
            self._consecutive = 0
            return False

        dist = float(np.linalg.norm(neutrinos[0].position - neutrinos[1].position))
        if dist < self.r_bound:
            self._consecutive += 1
        else:
            self._consecutive = 0

        return self._consecutive >= self.persistence
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```
pytest tests/test_detector.py -v
```
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```
git add src/stiff_medium/detector.py tests/test_detector.py
git commit -m "feat(detector): add BoundStateTracker"
```

---

### Task 8: Visualization (matplotlib animation)

**Files:**
- Create: `src/stiff_medium/visualize.py`
- Create: `tests/test_visualize.py`

Take a list of states (each state = list of neutrinos at one timestep) and produce a matplotlib animation. Plot positions as dots, velocity vectors as arrows. Color changes when bound state is detected.

- [ ] **Step 1: Write a smoke test for the animator**

`tests/test_visualize.py`:

```python
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for tests

from stiff_medium.neutrino import Neutrino, C
from stiff_medium.visualize import animate


def test_animate_returns_figure_without_error():
    s = C / np.sqrt(2)
    history = [
        [
            Neutrino(np.array([float(i) * 0.1, 0.0]), np.array([s, s])),
            Neutrino(np.array([1.0 - float(i) * 0.1, 0.0]), np.array([-s, s])),
        ]
        for i in range(20)
    ]
    bound_flags = [False] * 15 + [True] * 5

    fig, anim = animate(history, bound_flags, xlim=(-2, 2), ylim=(-2, 2))
    assert fig is not None
    assert anim is not None
```

- [ ] **Step 2: Run test to confirm it fails**

Run:
```
pytest tests/test_visualize.py -v
```
Expected: ImportError on `animate`.

- [ ] **Step 3: Implement `animate`**

`src/stiff_medium/visualize.py`:

```python
"""Matplotlib animation for the simulation. v1: 2D scatter of positions,
arrows for velocities, color change when bound state is flagged."""

from typing import Sequence
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

from stiff_medium.neutrino import Neutrino


def animate(
    history: Sequence[Sequence[Neutrino]],
    bound_flags: Sequence[bool],
    xlim: tuple[float, float] = (-5, 5),
    ylim: tuple[float, float] = (-5, 5),
    interval_ms: int = 50,
) -> tuple[plt.Figure, FuncAnimation]:
    """Create animation; caller can plt.show() or save."""
    if len(history) != len(bound_flags):
        raise ValueError("history and bound_flags must be same length")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_title("Stiff-Medium Path C: neutrino dynamics")
    ax.grid(True, alpha=0.3)

    scatter = ax.scatter([], [], s=80, c="steelblue")
    quivers: list = []
    label = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top")

    def init():
        scatter.set_offsets(np.zeros((0, 2)))
        return [scatter, label]

    def update(frame_idx: int):
        for q in quivers:
            q.remove()
        quivers.clear()

        state = history[frame_idx]
        positions = np.array([n.position for n in state])
        velocities = np.array([n.velocity for n in state])

        scatter.set_offsets(positions)
        scatter.set_color("crimson" if bound_flags[frame_idx] else "steelblue")

        q = ax.quiver(
            positions[:, 0], positions[:, 1],
            velocities[:, 0], velocities[:, 1],
            angles="xy", scale_units="xy", scale=2, color="gray", alpha=0.6,
        )
        quivers.append(q)

        label.set_text(
            f"step {frame_idx}/{len(history) - 1}  "
            f"{'BOUND' if bound_flags[frame_idx] else 'free'}"
        )
        return [scatter, q, label]

    anim = FuncAnimation(
        fig, update, frames=len(history),
        init_func=init, interval=interval_ms, blit=False, repeat=False,
    )
    return fig, anim
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```
pytest tests/test_visualize.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```
git add src/stiff_medium/visualize.py tests/test_visualize.py
git commit -m "feat(visualize): add 2D animation of simulation history"
```

---

### Task 9: Electron formation experiment

**Files:**
- Create: `scripts/electron_formation.py`

Compose everything into the actual experiment. Two neutrinos on collision course at 45° angles (one heading NE, one heading NW, paths intersecting at origin). Run for 1000 steps. Plot result. Print whether bound state formed.

- [ ] **Step 1: Create `scripts/electron_formation.py`**

```python
"""Path C experiment: do two colliding neutrinos form a bound state?

Setup:
- Neutrino A starts at (-2, 0) heading NE (vx=+s, vy=+s)
- Neutrino B starts at (+2, 0) heading NW (vx=-s, vy=+s)
- Their paths cross near (0, 2).
- Run 1000 steps, dt=0.01.

Per spec §2 (no correction loops): no parameter tuning beyond the basic
overlap radius and push amplitude — those represent the medium's
discretization and stiffness, NOT free fitting parameters.

Outcomes:
- BOUND: tracker flags persistent proximity → simulation produces an
  electron-like state, positive result.
- UNBOUND: particles pass through and separate → the displacement-only
  rule alone is insufficient; theory needs revision (per §2 methodology,
  this is a real falsification signal, not a tuning opportunity).
"""

import numpy as np
import matplotlib.pyplot as plt

from stiff_medium.neutrino import Neutrino, C
from stiff_medium.dynamics import step
from stiff_medium.detector import BoundStateTracker
from stiff_medium.visualize import animate


def main() -> None:
    s = C / np.sqrt(2)

    a = Neutrino(
        position=np.array([-2.0, 0.0]),
        velocity=np.array([s, s]),
    )
    b = Neutrino(
        position=np.array([2.0, 0.0]),
        velocity=np.array([-s, s]),
    )

    # These are the two simulation parameters: medium discretization
    # (r_overlap) and stiffness response (push). Both should be small
    # relative to the dynamics scale (c*dt).
    DT = 0.01
    R_OVERLAP = 0.05
    PUSH = 0.05
    R_BOUND = 0.5
    PERSISTENCE = 50  # consecutive steps within R_BOUND
    N_STEPS = 1000

    state = [a, b]
    tracker = BoundStateTracker(r_bound=R_BOUND, persistence=PERSISTENCE)
    history: list[list[Neutrino]] = [state]
    bound_flags: list[bool] = [False]

    bound_first_seen = -1
    for k in range(N_STEPS):
        state = step(state, dt=DT, r_overlap=R_OVERLAP, push=PUSH)
        flagged = tracker.update(state)
        history.append(state)
        bound_flags.append(flagged)
        if flagged and bound_first_seen < 0:
            bound_first_seen = k + 1

    if bound_first_seen >= 0:
        print(f"RESULT: BOUND state first detected at step {bound_first_seen}.")
    else:
        print("RESULT: NO bound state detected over 1000 steps.")
        print("       The displacement-only rule did not produce a stable orbit.")
        print("       Per spec §2 methodology, this is a falsification signal,")
        print("       not an invitation to tune parameters.")

    fig, _anim = animate(
        history, bound_flags,
        xlim=(-3, 3), ylim=(-3, 3),
        interval_ms=20,
    )
    plt.show()


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the experiment**

Run:
```
python scripts/electron_formation.py
```
Expected: opens an animation window, prints either "BOUND state first detected at step N" or "NO bound state detected" with the explicit falsification framing.

- [ ] **Step 3: Record the result in README**

Append a "## Result of v1 experiment" section to `README.md` with whatever the simulation reported. Be honest — record both the outcome and the parameter values used (DT, R_OVERLAP, PUSH, R_BOUND, PERSISTENCE).

Example template (fill in actual numbers):

```markdown
## Result of v1 experiment

Run on: 2026-04-29
Parameters: DT=0.01, R_OVERLAP=0.05, PUSH=0.05, R_BOUND=0.5, PERSISTENCE=50, N_STEPS=1000.

Outcome: [BOUND at step N | NO bound state observed]

Interpretation: [If BOUND: the rules generate a stable bound state from collision-only displacement, supporting Path A architecture. If UNBOUND: the rules as specified are insufficient; spec §5 needs revision, candidate causes include ... (medium back-reaction, wave emission, restoring force).]
```

- [ ] **Step 4: Run all tests one more time**

Run:
```
pytest -v
```
Expected: all tests pass.

- [ ] **Step 5: Commit**

```
git add scripts/electron_formation.py README.md
git commit -m "feat: add electron formation experiment + record v1 result"
```

---

## What "done" looks like for Path C v1

- All tests pass.
- The experiment script runs end-to-end and produces an animation + a printed verdict.
- The README contains the experiment's actual outcome (BOUND or UNBOUND) with parameter values, written honestly.

**Either outcome is a valid finding under the §2 methodology:**

- **BOUND** = positive result: the displacement-only rule alone produces a stable orbital state. Path A spec is supported. Move on to expanding the simulation (3D, bi-pyramids, lepton stress-loading).
- **UNBOUND** = falsification signal: the rule as currently specified does NOT produce stable orbits. Path A spec §5 needs direct revision (candidates: medium back-reaction with restoring force; wave emission and re-absorption; topological constraint not yet specified). No correction loop allowed.

In either case, the spec gets updated based on what we observed, and we move forward with reality grounding the theory.

---

## Self-review checklist

After implementing all tasks:

- [ ] All tests pass: `pytest -v`
- [ ] No TODOs, TBDs, or placeholders left in code or plan
- [ ] Method names consistent across files (`propagate`, `detect_overlap`, `displace`, `step`, `update`, `animate`)
- [ ] README reflects the actual experiment outcome (not pre-filled)
- [ ] No parameter tuning was done to force a desired outcome (per §2 methodology)
